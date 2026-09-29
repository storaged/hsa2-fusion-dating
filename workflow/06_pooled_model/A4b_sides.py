"""A4b: do the two fusion flanks record one switch-off event? Likelihood-ratio test in the pooled model (A4).

H0: one shared rho for F2a and F2b (one event);  H1: separate rho_a, rho_b.
LR = 2 (NLL_H0 - NLL_H1), df = 1. The H1 optimum is found by alternating profile maximisation over T_a and T_b.
Calibration: the same LR for random pairs of real ends (true T = 0 for both, so H0 holds exactly) gives an empirical
null distribution; the flank LR is reported with both the chi2_1 and the empirical p-value.
usage: python A4b_sides.py windows.tsv.gz out_dir [pan] [n_pairs=60] [procs=16]
"""
import sys
from multiprocessing import Pool

import numpy as np
import pandas as pd
from scipy.stats import chi2

win_f, OUT = sys.argv[1:3]
PAN = sys.argv[3] if len(sys.argv) > 3 else "chimp"
NPAIRS = int(sys.argv[4]) if len(sys.argv) > 4 else 60
PROCS = int(sys.argv[5]) if len(sys.argv) > 5 else 16
sys.argv = ["A4_pooled_model.py", win_f, OUT, PAN]
import A4_pooled_model as M  # noqa: E402  (loads the data)


def prof1(fixed, unit):
    """profile NLL over T for `unit`, other rho fixed as in `fixed`; returns (T_min, nll_min, start)"""
    best = (None, np.inf, None); start = None
    for T in M.GRID:
        f, x = M.fit({**fixed, unit: 1 - T / M.TS}, start)
        start = x
        if f < best[1]:
            best = (T, f, x)
    return best


def lr_test(pair):
    a, b = pair
    # H0
    T0, _, _, prof = M.profile([a, b])
    nll0 = prof.min()
    # H1: alternate
    Ta, Tb = T0, T0
    for _ in range(3):
        Ta, _, _ = prof1({b: 1 - Tb / M.TS}, a)
        Tb, nll1, _ = prof1({a: 1 - Ta / M.TS}, b)
    LR = max(2 * (nll0 - nll1), 0.0)
    return dict(pair=f"{a}+{b}", T_shared=T0, T_a=Ta, T_b=Tb, LR=LR, p_chi2=chi2.sf(LR, 1))


if __name__ == "__main__":
    real = [u for u in M.UNITS if not u.startswith("F")]
    rng = np.random.default_rng(44)
    pairs = [("F2a", "F2b")] + [tuple(rng.choice(real, 2, replace=False)) for _ in range(NPAIRS)]
    with Pool(PROCS) as pool:
        res = pool.map(lr_test, pairs)
    R = pd.DataFrame(res)
    fl, null = R.iloc[0], R.iloc[1:]
    R["is_flank"] = [True] + [False] * len(null)
    R.to_csv(f"{OUT}/A4b_sides_{PAN}.tsv", sep="\t", index=False)
    p_emp = (1 + (null.LR >= fl.LR).sum()) / (1 + len(null))
    print(f"flanks: T_shared {fl.T_shared:.1f}; T_2a {fl.T_a:.1f}; T_2b {fl.T_b:.1f}; LR {fl.LR:.2f}; "
          f"p(chi2_1) {fl.p_chi2:.3g}; p(empirical, {len(null)} real-end pairs) {p_emp:.3g}")
    print(f"null (real-end pairs): LR median {null.LR.median():.2f}, 95th pct {null.LR.quantile(.95):.2f}; "
          f"fraction with p_chi2 < 0.05: {(null.p_chi2 < 0.05).mean():.2f}")
    print("A4b_DONE")
