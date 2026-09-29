"""A4 v1: pooled model for T_off: a linear Gaussian mixed model over all ancestral ends and lineages.

Linear regime (manuscript Prop. 1; justified by A12c: present-day B <= ~0.25 at ends). Observations:
  y[u, l, j] = composition-corrected W->S vs S->W log-odds of unit u on lineage l in distance bin j,
               minus the pooled log-odds of lineage l's interior;  sampling variance v = 1/WS + 1/SW.
  units u: real human ends (hs1 chrom arm; acrocentric p-arms and chr2 near the fusion excluded) and the two
           fusion flanks F2a, F2b (hs1 windows within D_FLANK of the junction; human distance = distance to the
           junction, i.e. to the former telomere);
  lineages l: human (Anc4 -> hs1), Pan (panstem + chimp, or + bonobo), gorilla (Anc3 -> gorilla);
  bins j: lineage-specific distance to that lineage's physical telomere: 0-1, 1-2, 2-5, 5-10, 10-15, 15-20, 20-25 Mb.
Model:
  y[u,l,j] = a_u * K[l,j] * S[u,l] + e,   e ~ N(0, v + tau_l^2)
  a_u ~ N(mu, sigma^2) (end strength, shared by the lineages: same locus), integrated out analytically;
  K[human, 0-1 Mb] = 1 (scale); all other K free;
  S = 1 except human-lineage observations of F2a and F2b: S = rho = 1 - T_off / T_s (one fusion event).
Inference: maximum likelihood; profile likelihood for T_off on a grid (95% CI by LR, chi2_1 = 3.84).
Validation ('--pseudo'): each real end in turn is given its own free rho (true T = 0); report the coverage of its
95% profile CI and the bias; also random pairs of ends sharing one rho (as the two flanks do).
usage: python A4_pooled_model.py windows.tsv.gz out_dir [pan=chimp|bonobo] [--pseudo] [--no-pan] [--no-gorilla]
"""
import os
import sys
from multiprocessing import Pool

import numpy as np
import pandas as pd
from scipy import optimize

args = [a for a in sys.argv[1:] if not a.startswith("--")]
flags = {a for a in sys.argv[1:] if a.startswith("--")}
win_f, OUT = args[:2]
PAN = args[2] if len(args) > 2 else "chimp"
os.makedirs(OUT, exist_ok=True)
TS = 6.0
F = (113_940_058 + 114_049_496) / 2
D_FLANK = 10e6
ACRO = {"chr13", "chr14", "chr15", "chr21", "chr22"}
EDGES = np.array([0, 1, 2, 5, 10, 15, 20, 25]) * 1e6
NB = len(EDGES) - 1
LINS = ["human", "pan"] + ([] if "--no-gorilla" in flags else ["gorilla"])
if "--no-pan" in flags:
    LINS = [l for l in LINS if l != "pan"]

# ---------------- data ----------------
cols = ["chrom", "win_start", "branch", "n_callable", "parent_W", "parent_S", "WS", "SW", "sp_dist_end", "sp_arm_end"]
d = pd.read_csv(win_f, sep="\t", usecols=cols)
d = d[d.chrom.str.match(r"^chr\d+$") & (d.n_callable > 20_000)]
H = d[d.branch == "human"].set_index(["chrom", "win_start"])
tab = {"human": H[["WS", "SW", "parent_W", "parent_S", "sp_dist_end"]].copy()}
ps = d[d.branch == "panstem"].set_index(["chrom", "win_start"])
tip = d[d.branch == PAN].set_index(["chrom", "win_start"])
idx = ps.index.intersection(tip.index)
pan = ps.loc[idx, ["parent_W", "parent_S"]].copy()
pan["WS"] = ps.loc[idx, "WS"] + tip.loc[idx, "WS"]
pan["SW"] = ps.loc[idx, "SW"] + tip.loc[idx, "SW"]
pan["sp_dist_end"] = tip.loc[idx, "sp_dist_end"]           # present Pan karyotype distance
tab["pan"] = pan
tab["gorilla"] = d[d.branch == "gorilla"].set_index(["chrom", "win_start"])[["WS", "SW", "parent_W", "parent_S", "sp_dist_end"]]

mid = H.index.get_level_values(1) + 50_000
chrom = H.index.get_level_values(0)
near2 = (chrom == "chr2") & (np.abs(mid - F) < 30e6)
unit = pd.Series("", index=H.index)
is_end = (H.sp_dist_end < 25e6) & ~near2 & ~((H.sp_arm_end == "p") & chrom.isin(ACRO))
unit[is_end] = (chrom[is_end] + H.sp_arm_end[is_end]).values
fl = (chrom == "chr2") & (np.abs(mid - F) < D_FLANK)
unit[fl] = np.where(mid[fl] < F, "F2a", "F2b")
hdist = H.sp_dist_end.copy()
hdist[fl] = np.abs(mid[fl] - F)                              # human flank: distance to the former telomere
interior = (H.sp_dist_end >= 30e6) & ~near2


def lo(x):
    return np.log(x.WS.sum() / x.SW.sum()) - np.log(x.parent_W.sum() / x.parent_S.sum())


obs = []   # (unit, lineage, bin, y, v)
for l in LINS:
    t = tab[l]
    ii = t.index.intersection(H.index[interior])
    base = lo(t.loc[ii])
    dist = hdist if l == "human" else t.sp_dist_end
    for u in sorted(set(unit) - {""}):
        w = unit.index[unit == u].intersection(t.index)
        if len(w) == 0:
            continue
        dd = dist.reindex(w)
        b = np.digitize(dd.values, EDGES) - 1
        for j in range(NB):
            x = t.loc[w[b == j]]
            if len(x) == 0 or x.WS.sum() < 20 or x.SW.sum() < 20:
                continue
            obs.append((u, l, j, lo(x) - base, 1 / x.WS.sum() + 1 / x.SW.sum()))
O = pd.DataFrame(obs, columns=["unit", "lin", "bin", "y", "v"])
O.to_csv(f"{OUT}/A4_observations_{PAN}.tsv", sep="\t", index=False)
UNITS = sorted(O.unit.unique())
LI = {l: i for i, l in enumerate(LINS)}
groups = {u: (O[O.unit == u].lin.map(LI).values, O[O.unit == u].bin.values, O[O.unit == u].y.values,
              O[O.unit == u].v.values) for u in UNITS}
print(f"{len(UNITS)} units ({sum(u.startswith('F') for u in UNITS)} flanks), {len(O)} observations, lineages {LINS}",
      flush=True)

# ---------------- likelihood ----------------
# theta = [mu, log_sigma, K (len(LINS)*NB - 1, K[human,0] fixed at 1), log_tau (len(LINS))]
NK = len(LINS) * NB - 1


def unpack(th):
    mu, ls = th[0], th[1]
    K = np.concatenate([[1.0], th[2:2 + NK]]).reshape(len(LINS), NB)
    tau2 = np.exp(2 * th[2 + NK:2 + NK + len(LINS)])
    return mu, np.exp(ls), K, tau2


def unit_ll(u, mu, s2, K, tau2, rho_of_unit):
    li, bj, y, v = groups[u]
    c = K[li, bj].copy()
    r = rho_of_unit.get(u)
    if r is not None:
        c = np.where(li == LI["human"], c * r, c)
    Dg = v + tau2[li]
    # y ~ N(mu c, s2 c c' + diag(Dg)); Sherman-Morrison
    Di = 1 / Dg
    resid = y - mu * c
    cDc = np.sum(c * c * Di)
    q = np.sum(resid * resid * Di) - s2 * np.sum(c * resid * Di) ** 2 / (1 + s2 * cDc)
    logdet = np.sum(np.log(Dg)) + np.log(1 + s2 * cDc)
    return -0.5 * (q + logdet + len(y) * np.log(2 * np.pi))


def nll(th, rho_of_unit):
    mu, s, K, tau2 = unpack(th)
    return -sum(unit_ll(u, mu, s * s, K, tau2, rho_of_unit) for u in UNITS)


def th0():
    K0 = np.tile(np.array([1, .9, .7, .45, .3, .2, .15]), (len(LINS), 1)).ravel()[1:]
    return np.concatenate([[0.4, np.log(0.15)], K0, np.full(len(LINS), np.log(0.05))])


def fit(rho_of_unit, start=None):
    r = optimize.minimize(nll, th0() if start is None else start, args=(rho_of_unit,), method="L-BFGS-B")
    return r.fun, r.x


GRID = np.round(np.arange(-3.0, TS + 1e-9, 0.1), 2)      # T_off grid (negative allowed for pseudo-fusions)


def profile(units_with_rho):
    """profile NLL over T for a set of units sharing one rho = 1 - T/TS"""
    best, start, prof = None, None, []
    for T in GRID:
        f, x = fit({u: 1 - T / TS for u in units_with_rho}, start)
        start = x
        prof.append(f)
    prof = np.array(prof)
    k = prof.argmin()
    ok = GRID[prof - prof[k] <= 1.92]
    return GRID[k], ok.min(), ok.max(), prof


def main():
  if "--pseudo" not in flags:
      That, lo_, hi_, prof = profile(["F2a", "F2b"])
      f, x = fit({"F2a": 1 - That / TS, "F2b": 1 - That / TS})
      mu, s, K, tau2 = unpack(x)
      res = dict(pan=PAN, lineages="+".join(LINS), T_hat=That, T_lo=lo_, T_hi=hi_, mu=mu, sigma=s,
                 cv_between_ends=s / mu, **{f"tau_{l}": np.sqrt(tau2[i]) for i, l in enumerate(LINS)})
      for side in ("F2a", "F2b"):
          t1, l1, h1, _ = profile([side])
          res.update({f"T_{side}": t1, f"T_{side}_lo": l1, f"T_{side}_hi": h1})
      pd.DataFrame([res]).to_csv(f"{OUT}/A4_fit_{PAN}_{'+'.join(LINS)}.tsv", sep="\t", index=False)
      pd.DataFrame({"T": GRID, "nll": prof}).to_csv(f"{OUT}/A4_profile_{PAN}_{'+'.join(LINS)}.tsv", sep="\t", index=False)
      print(pd.Series(res).to_string())
      print("K (rows = lineages", LINS, "; cols = 0-1,1-2,2-5,5-10,10-15,15-20,20-25 Mb):")
      print(np.round(K * mu, 3))
      print("A4_DONE")
  else:
      real = [u for u in UNITS if not u.startswith("F")]
      rng = np.random.default_rng(4)
      pairs = [tuple(rng.choice(real, 2, replace=False)) for _ in range(int(os.environ.get("A4_NPAIRS", 40)))]
      tasks = [[u] for u in real] + [list(p) for p in pairs]
      with Pool(int(os.environ.get("A4_PROCS", 12))) as pool:
          out = pool.map(profile, tasks)
      rows = [dict(units="+".join(t), n=len(t), T_hat=o[0], T_lo=o[1], T_hi=o[2], covers0=o[1] <= 0 <= o[2])
              for t, o in zip(tasks, out)]
      R = pd.DataFrame(rows)
      R.to_csv(f"{OUT}/A4_pseudo_{PAN}_{'+'.join(LINS)}.tsv", sep="\t", index=False)
      for n in (1, 2):
          x = R[R.n == n]
          print(f"pseudo-fusions with {n} end(s): {len(x)}; coverage of T=0 by 95% CI = {x.covers0.mean():.2f}; "
                f"median T_hat {x.T_hat.median():.2f}; mean CI width {(x.T_hi - x.T_lo).mean():.2f}")
      print("A4_PSEUDO_DONE")


if __name__ == "__main__":
    main()
