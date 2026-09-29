"""A12: calibrate the human-branch gBGC signal against measured present-day recombination (deCODE).

Question 1: which present-day recombination measure predicts the W->S excess on the human branch
            across real chromosome ends and the interior? (paternal vs maternal; NCO vs CO vs DSB)
Question 2: given that calibration, what excess would the fusion flanks carry had they stayed ends
            for the whole human lineage, and what fraction of it do they carry? (rho_rec)
Question 3: how much of the heterogeneity between real ends does present-day recombination explain?
            (this sets the prior width sigma_A of the end strengths in the model)

Model (non-CpG, human branch, 100-kb hs1 windows; deCODE 1-Mb bins joined via hs1->hg38 liftover):
  WS_w | WS_w + SW_w ~ Binomial(pi_w),  logit pi_w = log(pW_w / pS_w) + a + b * x_w
  quasi-binomial dispersion from Pearson chi2; 95% CI by 5-Mb block bootstrap.
Fusion flanks (hs1 chr2 within 30 Mb of the fusion) and acrocentric p-arm ends are excluded from the fit.
usage: python A12_recomb_calibration.py results/A3/windows.tsv.gz data results/A12
"""
import os
import subprocess
import sys
import tempfile

import numpy as np
import pandas as pd
from scipy import optimize

win_f, DATA, OUT = sys.argv[1:4]
os.makedirs(OUT, exist_ok=True)
rng = np.random.default_rng(12)
FUSION = (113_940_058 + 114_049_496) / 2
LO = os.path.join(os.path.dirname(sys.executable), "liftOver")

# ---------------- human-branch windows, classes as in A3 ----------------
d = pd.read_csv(win_f, sep="\t", dtype={"chrom": str, "sp_name": str})
d = d[(d.branch == "human") & d.chrom.str.match(r"^chr\d+$") & (d.n_callable > 20_000) & d.sp_dist_end.notna()].copy()
d["mid"] = d.win_start + 50_000
d["dist"] = d.sp_dist_end
d["cls"] = np.where(d.dist < 30e6, "end", "interior")
d["end_id"] = d.sp_name.astype(str) + "_" + d.sp_arm_end
fl = (d.chrom == "chr2") & ((d.mid - FUSION).abs() < 30e6)
d.loc[fl, "cls"] = np.where(d.loc[fl, "mid"] < FUSION, "fusion_2a", "fusion_2b")
d.loc[fl, "dist"] = (d.loc[fl, "mid"] - FUSION).abs()
acro = d.chrom.isin(["chr13", "chr14", "chr15", "chr21", "chr22"]) & (d.sp_arm_end == "p") & (d.cls == "end")
d.loc[acro, "cls"] = "end_acro_p"
d = d[(d.WS + d.SW) > 0]

# ---------------- join deCODE via hs1 -> hg38 liftover of window midpoints ----------------
with tempfile.TemporaryDirectory() as t:
    bed = d[["chrom", "mid"]].assign(e=d.mid + 1, name=np.arange(len(d)))
    bed.to_csv(f"{t}/in.bed", sep="\t", header=False, index=False)
    subprocess.run([LO, f"{t}/in.bed", f"{DATA}/chains/hs1ToHg38.over.chain.gz", f"{t}/out.bed", f"{t}/un"], check=True)
    lo = pd.read_csv(f"{t}/out.bed", sep="\t", header=None, names=["c38", "p38", "e38", "name"])
d = d.reset_index(drop=True)
d = d.join(lo.set_index("name"))
d = d[d.c38 == d.chrom].copy()
d["bin"] = (d.p38 // 1_000_000) * 1_000_000 + 500_000

VARS = {"map": "nco", "cMperMb": "co", "DSB": "dsb"}
for sex in ("pat", "mat"):
    m = pd.read_csv(f"{DATA}/decode/maps.{sex}.tsv", sep="\t", comment="#").rename(columns={"Chr": "chrom", "pos": "bin"})
    m = m[["chrom", "bin"] + list(VARS)].rename(columns={k: f"{v}_{sex}" for k, v in VARS.items()})
    d = d.merge(m, on=["chrom", "bin"], how="left")
PRED = [f"{v}_{s}" for s in ("pat", "mat") for v in VARS.values()]
d["block"] = d.chrom + ":" + (d.win_start // 5_000_000).astype(str)
d.to_csv(f"{OUT}/A12_windows.tsv.gz", sep="\t", index=False)


# ---------------- binomial GLM with offset ----------------
def fit(df, xcols):
    k = df.WS.values.astype(float); n = (df.WS + df.SW).values.astype(float)
    off = np.log(df.parent_W.values / df.parent_S.values)
    X = np.column_stack([np.ones(len(df))] + [df[c].values for c in xcols])

    def nll(beta):
        eta = off + X @ beta
        return -(k * eta - n * np.logaddexp(0, eta)).sum()

    def grad(beta):
        p = 1 / (1 + np.exp(-(off + X @ beta)))
        return -X.T @ (k - n * p)
    r = optimize.minimize(nll, np.zeros(X.shape[1]), jac=grad, method="BFGS")
    p = 1 / (1 + np.exp(-(off + X @ r.x)))
    phi = (((k - n * p) ** 2) / (n * p * (1 - p))).sum() / (len(df) - X.shape[1])
    H = (X * (n * p * (1 - p))[:, None]).T @ X
    se = np.sqrt(np.diag(np.linalg.inv(H)) * phi)
    return r.x, se, r.fun, phi


def pooled_excess(df):
    """pooled log-odds of W->S vs S->W, corrected for W/S opportunity"""
    return np.log(df.WS.sum() / df.SW.sum()) - np.log(df.parent_W.sum() / df.parent_S.sum())


train = d[d.cls.isin(["end", "interior"])].dropna(subset=PRED).copy()
for c in PRED:  # standardise within the training set so slopes are comparable
    mu, sd = train[c].mean(), train[c].std()
    for df in (train, d):
        df[c + "_z"] = (df[c] - mu) / sd

rows = []
_, _, nll0, _ = fit(train, [])
for c in PRED:
    b, se, nll, phi = fit(train, [c + "_z"])
    rows.append(dict(model=c, slope_per_sd=b[1], se=se[1], dev_gain=2 * (nll0 - nll), dispersion=phi))
for combo in (["nco_pat_z", "nco_mat_z"], ["co_pat_z", "co_mat_z"], ["nco_pat_z", "co_pat_z"]):
    b, se, nll, phi = fit(train, combo)
    rows.append(dict(model="+".join(combo), slope_per_sd=np.nan, se=np.nan, dev_gain=2 * (nll0 - nll), dispersion=phi,
                     detail="; ".join(f"{c}={bb:.4f}({s:.4f})" for c, bb, s in zip(combo, b[1:], se[1:]))))
mods = pd.DataFrame(rows).sort_values("dev_gain", ascending=False)
mods.to_csv(f"{OUT}/A12_model_comparison.tsv", sep="\t", index=False)
print("== Q1: which recombination measure predicts the human-branch W->S excess (ends + interior)?")
print(mods.round(4).to_string(index=False))

# ---------------- Q2: calibrated expectation at the fusion flanks ----------------
best = mods[~mods.model.str.contains(r"\+")].iloc[0].model
xb = best + "_z"
print(f"\nbest single predictor: {best}")
ends = d[d.cls == "end"].dropna(subset=[xb])
inter = d[d.cls == "interior"].dropna(subset=[xb])
# telomere kernel of the best predictor at real ends (median over ends per 1-Mb distance bin)
ends_b = ends.assign(dbin=(ends.dist // 1e6).astype(int))
K = ends_b.groupby("dbin")[xb].median()
x_int = inter[xb].mean()


def rho_rec(dd, D, slope):
    f = d[d.cls.isin(dd) & (d.dist < D * 1e6)]
    e = ends[ends.dist < D * 1e6]
    eps_f = pooled_excess(f) - pooled_excess(inter)
    eps_e = pooled_excess(e) - pooled_excess(inter)
    pred_f = slope * (K.reindex((f.dist // 1e6).astype(int)).values - x_int).mean()   # as if still an end
    pred_e = slope * (e[xb].mean() - x_int)
    return eps_f, eps_e, pred_f, pred_e


b_hat = fit(train, [xb])[0][1]
blocks = train.block.unique()
boot = []
for _ in range(200):
    bs = rng.choice(blocks, len(blocks))
    tb = pd.concat([train[train.block == x] for x in bs])
    boot.append(fit(tb, [xb])[0][1])
b_lo, b_hi = np.percentile(boot, [2.5, 97.5])
out = []
for D in (2, 5, 10):
    for side, dd in (("2a", ["fusion_2a"]), ("2b", ["fusion_2b"]), ("both", ["fusion_2a", "fusion_2b"])):
        eps_f, eps_e, pred_f, pred_e = rho_rec(dd, D, b_hat)
        out.append(dict(first_Mb=D, side=side, excess_flank=eps_f, excess_ends_obs=eps_e, excess_ends_pred=pred_e,
                        pred_flank_if_end=pred_f, rho_rec=eps_f / pred_f,
                        rho_rec_lo=eps_f / (pred_f * b_hi / b_hat), rho_rec_hi=eps_f / (pred_f * b_lo / b_hat)))
out = pd.DataFrame(out)
out.to_csv(f"{OUT}/A12_flank_expectation.tsv", sep="\t", index=False)
print(f"\n== Q2: slope {b_hat:.4f} per SD (block-bootstrap 95% CI {b_lo:.4f}-{b_hi:.4f})")
print("calibration check: observed vs predicted excess at real ends; retained fraction at the flanks (slope CI only)")
print(out.round(3).to_string(index=False))

# ---------------- Q3: between-end heterogeneity explained ----------------
pe = []
for eid, g in ends[ends.dist < 5e6].groupby("end_id"):
    if (g.WS + g.SW).sum() < 200:
        continue
    pe.append(dict(end_id=eid, obs=pooled_excess(g) - pooled_excess(inter), pred=b_hat * (g[xb].mean() - x_int),
                   n=(g.WS + g.SW).sum()))
pe = pd.DataFrame(pe)
pe["resid"] = pe.obs - pe.pred
se_bin = np.sqrt(1 / pe.n * 4)  # rough sampling SD of a log-odds from n events (p~0.5)
v_obs, v_res, v_samp = pe.obs.var(), pe.resid.var(), (se_bin ** 2).mean()
pe.to_csv(f"{OUT}/A12_per_end.tsv", sep="\t", index=False)
print(f"\n== Q3: {len(pe)} real ends (first 5 Mb). Var(obs excess) {v_obs:.4f}; Var(residual) {v_res:.4f}; "
      f"sampling var ~{v_samp:.4f}")
print(f"   share of between-end variance (beyond sampling) explained by {best}: "
      f"{1 - max(v_res - v_samp, 0) / max(v_obs - v_samp, 1e-9):.2f}; corr(obs,pred) = {np.corrcoef(pe.obs, pe.pred)[0, 1]:.2f}")
print("A12_DONE")
