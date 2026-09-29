"""A10d: does a population bottleneck bias the switch-off estimator? Exact Wright-Fisher with a changing size.

As A10a (single-genome derived-on-branch counts relative to an ancestral genome sampled at the split; stationary
start; T_hat = T_s (1 - rho)), but the population size follows a schedule N(t). gBGC acts through the per-generation
conversion bias b (constant), so the population-scaled strength B(t) = 4 N(t) b drops during a bottleneck.
Scenario units: baseline N0; T_s = 12 N0 (~6 Mya for N0 = 2e4, g = 25 y); switch-off T_off = 5.6 N0 (~2.8 Mya).
Bottleneck: size N0 * f for a duration d * N0 generations, centred at t_b before the present, with
t_b in {1.8 N0 (~0.9 Mya, the proposed Early-Middle Pleistocene bottleneck), 5.6 N0 (during the switch-off),
8 N0 (before)}; f in {0.2, 0.1, 0.05}; d = 0.2. B0 = 4 N0 b in {0.5, 1}. Step and gradual (4 N0 ramp) switch-off.
Reference ('full end') and switch-off runs share the same demography, as real ends and flanks share one history.
usage: python A10d_bottleneck.py results/A10 [N0=300]
"""
import os
import sys
from functools import lru_cache

import numpy as np
import pandas as pd
from scipy.stats import binom

OUT = sys.argv[1]
N0 = int(sys.argv[2]) if len(sys.argv) > 2 else 300
os.makedirs(OUT, exist_ok=True)
TS = 12 * N0


@lru_cache(maxsize=None)
def P(G1, G2, sq):
    """transition from G1 to G2 gene copies with genic selection s = sq / 1e6"""
    s = sq / 1e6
    x = np.arange(G1 + 1) / G1
    ps = np.clip(x * (1 + s) / (1 + x * s), 0, 1)
    return binom.pmf(np.arange(G2 + 1)[None, :], G2, ps[:, None])


def sgrid(b):
    return int(round(b * 1e6))


def stationary(G, sq):
    v = np.zeros(G + 1)
    for _ in range(12 * G // 2):
        v[1] += 1.0
        v = v @ P(G, G, sq)
    v[0] = v[G] = 0.0
    return v


_st = {}


def run(bt, sizes):
    """bt: conversion bias per generation (forward from split); sizes: diploid N per generation. Returns E[WS], E[SW]."""
    res = {}
    for sign, key in ((1, "WS"), (-1, "SW")):
        G0 = 2 * sizes[0]
        k = (G0, sign * sgrid(bt[0]))
        if k not in _st:
            _st[k] = stationary(G0, sign * sgrid(bt[0]))
        v0 = _st[k]
        xg = np.arange(G0 + 1) / G0
        a, b_, new = v0 * (1 - xg), v0 * xg, np.zeros(G0 + 1)
        G = G0
        for t in range(TS):
            G2 = 2 * sizes[t]
            new[1] += G / G0                       # 2N_t mu new mutations, each at frequency 1/(2N_t)
            Pt = P(G, G2, sign * sgrid(bt[t]))
            new, a, b_ = new @ Pt, a @ Pt, b_ @ Pt
            G = G2
        xg = np.arange(G + 1) / G
        res[key] = (new @ xg + a @ xg, b_ @ (1 - xg))
    return res["WS"][0] + res["SW"][1], res["SW"][0] + res["WS"][1]


rows = []
age = TS - np.arange(TS)                           # generations before present
for B0 in (0.5, 1.0):
    b = B0 / (4 * N0)
    for bn in ("none", 1.8, 5.6, 8.0):
        for f in ((1.0,) if bn == "none" else (0.2, 0.1, 0.05)):
            sizes = np.full(TS, N0)
            if bn != "none":
                sizes[np.abs(age - bn * N0) <= 0.1 * N0] = max(int(N0 * f), 5)
            full = run(np.full(TS, b), sizes)
            Lf = np.log(full[0] / full[1])
            for kind in ("step", "ramp"):
                toff = 5.6 * N0
                S = (age > toff).astype(float) if kind == "step" else np.clip((age - (toff - 2 * N0)) / (4 * N0), 0, 1)
                o = run(b * S, sizes)
                rho = np.log(o[0] / o[1]) / Lf
                That = TS * (1 - rho)
                rows.append(dict(N0=N0, B0=B0, bottleneck_at_N0=bn, size_frac=f, kind=kind, T_off_N0=5.6,
                                 T_hat_N0=That / N0, bias_N0=(That - toff) / N0))
                print(rows[-1], flush=True)
R = pd.DataFrame(rows)
R.to_csv(f"{OUT}/A10d_bottleneck_N{N0}.tsv", sep="\t", index=False)
pd.set_option("display.width", 200)
print(R.pivot_table(index=["B0", "kind"], columns=["bottleneck_at_N0", "size_frac"], values="bias_N0").round(3).to_string())
print("A10d_DONE")
