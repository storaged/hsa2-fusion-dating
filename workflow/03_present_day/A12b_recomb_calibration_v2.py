"""A12b (v2 of A12): recombination calibration with log-scale predictors and an explicit end-calibration score.

A12 v1 finding: linear predictors fit the genome bulk but under-predict the excess at real ends (0.18 vs 0.46
in the first 2 Mb), so model choice by overall deviance is dominated by the interior. Here:
  * predictors on log scale (log(rate + 1% of median)), single and paternal+maternal pairs;
  * score 1 = deviance gain (all end + interior windows);
  * score 2 = end calibration: pooled observed vs predicted log-odds excess at real ends in distance bins
    0-2, 2-5, 5-10, 10-20 Mb (RMS error, the criterion that matters for the flank expectation);
  * 95% CI by Poisson-weighted 5-Mb block bootstrap (fast, weights on windows).
Input: results/A12/A12_windows.tsv.gz (written by A12). usage: python A12b_... results/A12
"""
import sys
import numpy as np
import pandas as pd
from scipy import optimize

OUT = sys.argv[1]
rng = np.random.default_rng(122)
d = pd.read_csv(f"{OUT}/A12_windows.tsv.gz", sep="\t")
RAW = [f"{v}_{s}" for s in ("pat", "mat") for v in ("nco", "co", "dsb")]
d = d.dropna(subset=RAW).copy()
for c in RAW:
    L = np.log(d[c].clip(lower=0) + 0.01 * d[c].median())
    d[c + "_L"] = (L - L.mean()) / L.std()
train = d[d.cls.isin(["end", "interior"])].copy()
ends = train[train.cls == "end"]
inter = train[train.cls == "interior"]
BINS = [(0, 2), (2, 5), (5, 10), (10, 20)]


def fit(df, cols, w=None):
    k = df.WS.values.astype(float); n = (df.WS + df.SW).values.astype(float)
    w = np.ones(len(df)) if w is None else w
    off = np.log(df.parent_W.values / df.parent_S.values)
    X = np.column_stack([np.ones(len(df))] + [df[c].values for c in cols])
    f = lambda b: -(w * (k * (off + X @ b) - n * np.logaddexp(0, off + X @ b))).sum()
    g = lambda b: -X.T @ (w * (k - n / (1 + np.exp(-(off + X @ b)))))
    r = optimize.minimize(f, np.zeros(X.shape[1]), jac=g, method="BFGS")
    return r.x, r.fun


def pooled(df):
    return np.log(df.WS.sum() / df.SW.sum()) - np.log(df.parent_W.sum() / df.parent_S.sum())


def pred_pooled(df, cols, b):
    """model-implied pooled log-odds excess of a set of windows (expected WS and SW summed)"""
    off = np.log(df.parent_W.values / df.parent_S.values)
    X = np.column_stack([np.ones(len(df))] + [df[c].values for c in cols])
    p = 1 / (1 + np.exp(-(off + X @ b))); n = (df.WS + df.SW).values
    return np.log((n * p).sum() / (n * (1 - p)).sum()) - np.log(df.parent_W.sum() / df.parent_S.sum())


base_i = pooled(inter)
_, nll0 = fit(train, [])
models = [[c + "_L"] for c in RAW] + [["nco_pat_L", "nco_mat_L"], ["co_pat_L", "co_mat_L"], ["dsb_pat_L", "dsb_mat_L"]]
rows = []
for cols in models:
    b, nll = fit(train, cols)
    bi = pred_pooled(inter, cols, b)
    cal = []
    for lo, hi in BINS:
        e = ends[(ends.dist >= lo * 1e6) & (ends.dist < hi * 1e6)]
        cal.append((pooled(e) - base_i, pred_pooled(e, cols, b) - bi))
    rms = np.sqrt(np.mean([(o - p) ** 2 for o, p in cal]))
    rows.append(dict(model="+".join(cols), coefs=" ".join(f"{x:.3f}" for x in b[1:]), dev_gain=2 * (nll0 - nll),
                     end_cal_rms=rms, **{f"obs_{lo}-{hi}": o for (lo, hi), (o, _) in zip(BINS, cal)},
                     **{f"pred_{lo}-{hi}": p for (lo, hi), (_, p) in zip(BINS, cal)}))
M = pd.DataFrame(rows).sort_values("end_cal_rms")
M.to_csv(f"{OUT}/A12b_models.tsv", sep="\t", index=False)
pd.set_option("display.width", 250)
print("== models ranked by end calibration (observed vs predicted excess at real ends by distance)")
print(M.round(3).to_string(index=False))

best = M.iloc[0].model.split("+")
b_hat, _ = fit(train, best)
bi = pred_pooled(inter, best, b_hat)
# telomere kernel of the best predictors: median per 1-Mb distance bin over real-end windows
ends_b = ends.assign(dbin=(ends.dist // 1e6).astype(int))
Kmed = ends_b.groupby("dbin")[best].median()


def flank_rows(b, bi_):
    out = []
    for D in (2, 5, 10):
        for side, cl in (("2a", ["fusion_2a"]), ("2b", ["fusion_2b"]), ("both", ["fusion_2a", "fusion_2b"])):
            f = d[d.cls.isin(cl) & (d.dist < D * 1e6)]
            ff = f.copy()
            ff[best] = Kmed.reindex((f.dist // 1e6).astype(int)).values   # as if still a typical end
            e = ends[ends.dist < D * 1e6]
            out.append(dict(first_Mb=D, side=side, obs_flank=pooled(f) - base_i,
                            pred_flank_now=pred_pooled(f, best, b) - bi_,
                            pred_flank_if_end=pred_pooled(ff, best, b) - bi_,
                            obs_ends=pooled(e) - base_i, pred_ends=pred_pooled(e, best, b) - bi_))
    o = pd.DataFrame(out)
    o["rho_rec"] = o.obs_flank / o.pred_flank_if_end
    return o


F = flank_rows(b_hat, bi)
blocks = train.block.values
ub = np.unique(blocks)
boot = []
for _ in range(200):
    wb = dict(zip(ub, rng.poisson(1.0, len(ub))))
    w = np.array([wb[x] for x in blocks], float)
    bb, _ = fit(train, best, w)
    boot.append(flank_rows(bb, pred_pooled(inter, best, bb)).rho_rec.values)
boot = np.array(boot)
F["rho_lo"], F["rho_hi"] = np.percentile(boot, 2.5, axis=0), np.percentile(boot, 97.5, axis=0)
F.to_csv(f"{OUT}/A12b_flank.tsv", sep="\t", index=False)
print(f"\n== best model: {'+'.join(best)}; flank expectation (rho_rec = observed / expected-if-still-an-end; "
      f"95% block-bootstrap CI on the calibration only)")
print(F.round(3).to_string(index=False))

pe = []
for eid, g in ends[ends.dist < 5e6].groupby("end_id"):
    if (g.WS + g.SW).sum() >= 200:
        pe.append(dict(end_id=eid, obs=pooled(g) - base_i, pred=pred_pooled(g, best, b_hat) - bi, n=(g.WS + g.SW).sum()))
pe = pd.DataFrame(pe); pe["resid"] = pe.obs - pe.pred
vs = (4 / pe.n).mean()
print(f"\n== per-end (first 5 Mb, {len(pe)} ends): corr(obs,pred) {np.corrcoef(pe.obs, pe.pred)[0,1]:.2f}; "
      f"mean resid {pe.resid.mean():.3f}; between-end var explained "
      f"{1 - max(pe.resid.var() - vs, 0) / (pe.obs.var() - vs):.2f}; residual SD beyond sampling "
      f"{np.sqrt(max(pe.resid.var() - vs, 0)):.3f}")
pe.to_csv(f"{OUT}/A12b_per_end.tsv", sep="\t", index=False)
print("A12b_DONE")
