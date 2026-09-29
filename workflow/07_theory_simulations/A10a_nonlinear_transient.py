"""A10a: exact Wright-Fisher computation of the gBGC switch-off estimator's bias (manuscript Prop. 2 remark).

Expected single-genome branch counts of W->S and S->W derived alleles that arose on the lineage after the split
(time 0 = split, forward to the present at T_s), under a time-varying gBGC strength B(t) at an ex-end:
  full end:       B(t) = B_end throughout                   -> reference log-ratio L_full
  switch-off:     B(t) = B_end * S(t),  S = 1 before, 0 after, with
                  'step'  : instantaneous at T_off
                  'ramp'  : linear decline over 4N generations centred on T_off (neutral-sojourn-like; midpoint = T_off)
Counts are derived-on-branch alleles in one present genome relative to one ancestral genome sampled at the split,
including polymorphisms segregating at the split (stationary burn-in).
Statistic: L = log(E[WS] / E[SW]) (neutral: 0).  rho = L / L_full,  T_hat = T_s (1 - rho)  [primary; the neutral
divergence depth T_eff, which includes ancestral coalescence, over-scales the estimate: reported as alternative],
bias = T_hat - T_off (T_off counted back from the present), for all counts and for fixed-only counts.
Genic selection p' = p(1+s)/(1+ps), s = B/(2N_genes)... here s = B / (4N) with N diploids (2N genes).
Diffusion scaling: results depend on B and on times in units of N; checked at two population sizes.
usage: python A10a_nonlinear_transient.py results/A10 [N=200] [Ts_in_N=8]
"""
import os
import sys
from functools import lru_cache

import numpy as np
from scipy.stats import binom

OUT = sys.argv[1]
N = int(sys.argv[2]) if len(sys.argv) > 2 else 200
TSN = float(sys.argv[3]) if len(sys.argv) > 3 else 8.0
os.makedirs(OUT, exist_ok=True)
G = 2 * N
i = np.arange(G + 1)
xg = i / G
TS = int(TSN * N)


@lru_cache(maxsize=None)
def P(Bq):
    s = Bq / 20.0 / (4 * N)
    ps = np.clip(xg * (1 + s) / (1 + xg * s), 0, 1)
    return binom.pmf(i[None, :], G, ps[:, None])


def Bgrid(B):
    return int(round(B * 20))   # 0.05 grid (memory: <= ~100 cached matrices)


_stat = {}


def stationary(Bq):
    """expected frequency-mass vector of one mutation class at mutation-selection-drift stationarity"""
    if Bq not in _stat:
        v = np.zeros(G + 1)
        for _ in range(12 * N):
            v[1] += 1.0
            v = v @ P(Bq)
        v[0] = 0.0; v[G] = 0.0                    # lost / fixed before the split contribute nothing to the branch
        _stat[Bq] = v
    return _stat[Bq]


def run(Bt, B_pre):
    """Expected derived-on-branch counts in one present genome relative to one ancestral genome sampled at the split.
    Includes mutations segregating at the split (weight 1-y if the ancestor lacks them; y if it carries them, which
    then appear as the reverse substitution class when lost from the present genome). Returns per class
    (all, fixed-only) counts."""
    res = {}
    for sign, key in ((1, "WS"), (-1, "SW")):
        v0 = stationary(sign * Bgrid(B_pre))
        a, b, new = v0 * (1 - xg), v0 * xg, np.zeros(G + 1)
        for t in range(TS):
            new[1] += 1.0
            Pt = P(sign * Bgrid(Bt[t]))
            new, a, b = new @ Pt, a @ Pt, b @ Pt
        res[key] = dict(same=(new @ xg + a @ xg, new[G] + a[G]), rev=(b @ (1 - xg), b[0]))
    return {"WS": tuple(res["WS"]["same"][k] + res["SW"]["rev"][k] for k in (0, 1)),
            "SW": tuple(res["SW"]["same"][k] + res["WS"]["rev"][k] for k in (0, 1))}


def L(o, idx):
    return np.log(o["WS"][idx] / o["SW"][idx])


neu = run(np.zeros(TS), 0.0)
TEFF = (neu["WS"][0] * G, neu["WS"][1] * G)    # effective branch depth in generations (all, fixed-only); cohorts enter at 1/G
print(f"T_s = {TS} gen; effective depth all = {TEFF[0]:.0f} gen, fixed-only = {TEFF[1]:.0f} gen", flush=True)
rows = []
for B in (0.5, 1.0, 2.0, 3.0, 5.0):
    full = run(np.full(TS, B), B)
    for kind in ("step", "ramp"):
        for toff_frac in (0.15, 0.3, 0.5, 0.7):
            toff = toff_frac * TS                    # generations before present
            t_fwd = np.arange(TS)                    # generation index from split
            age = TS - t_fwd                         # generations before present
            if kind == "step":
                S = (age > toff).astype(float)
            else:
                S = np.clip((age - (toff - 2 * N)) / (4 * N), 0, 1)
            o = run(B * S, B)
            for idx, lab in ((0, "all"), (1, "fixed_only")):
                rho = L(o, idx) / L(full, idx)
                That = TS * (1 - rho)                    # primary: scale by the split time T_s
                That_eff = TEFF[idx] * (1 - rho)         # alternative: scale by the neutral divergence depth
                rows.append(dict(N=N, Ts_in_N=TSN, Teff_in_N=TEFF[idx] / N, B=B, kind=kind, Toff_in_N=toff / N, rho=rho,
                                 That_in_N=That / N, bias_in_N=(That - toff) / N,
                                 bias_Teff_in_N=(That_eff - toff) / N, counts=lab))
            print(f"B={B} {kind} Toff/N={toff / N:.2f}  bias(all)={rows[-2]['bias_in_N']:+.3f} N  "
                  f"bias(fixed)={rows[-1]['bias_in_N']:+.3f} N", flush=True)
import pandas as pd
R = pd.DataFrame(rows)
R.to_csv(f"{OUT}/A10a_bias_N{N}_Ts{TSN:g}.tsv", sep="\t", index=False)
print("bias (T_hat - T_off) in units of N; T_hat = T_s (1 - rho)")
print(R.pivot_table(index=["B", "kind", "counts"], columns="Toff_in_N", values="bias_in_N").round(3).to_string())
print("alternative scaling by neutral divergence depth T_eff")
print(R.pivot_table(index=["B", "kind", "counts"], columns="Toff_in_N", values="bias_Teff_in_N").round(3).to_string())
print("A10a_DONE")
