"""A6b: present-day gBGC strength B = 4 Ne b from derived allele frequency spectra (DAF).

Model (Glemin et al. 2015 Genome Res; Eyre-Walker et al. 2006 style nuisance parameters):
  unfolded SFS over i = 1..n-1 derived copies (n = sample size of haploid genomes, projected),
  neutral (GC-conservative WW/SS) :  E[k_i] = theta_N * r_i / i
  W->S                                 E[k_i] = theta_WS * r_i * f_i(+B)
  S->W                                 E[k_i] = theta_SW * r_i * f_i(-B)
  f_i(B) = int_0^1 C(n,i) x^i (1-x)^(n-i) * (1 - exp(-B(1-x))) / (x (1-x) (1 - exp(-B))) dx
  r_i: free demographic/linked-selection distortion shared across classes (r_1 = 1)
  polarisation error eps: fraction of sites with swapped ancestral/derived state (class also swaps WS<->SW)
Poisson likelihood; B, eps, thetas, r_i by maximum likelihood; 95% CI by profile likelihood on B.
Input: TSV with columns region_class, cls (WS/SW/WW/SS), dac (derived allele count), n (haploid sample size)
       (projected to a common n by hypergeometric down-sampling).
usage: python A6_gbgc_sfs.py snps_polarised.tsv.gz out.tsv [n_proj=40]
"""
import sys

import numpy as np
import pandas as pd
from scipy import integrate, optimize
from scipy.special import comb, gammaln
from scipy.stats import hypergeom


def project(dac, n, m):
    """hypergeometric projection of derived counts (dac of n) to sample size m -> expected SFS over 0..m"""
    sfs = np.zeros(m + 1)
    for (a, nn), c in pd.Series(1, index=pd.MultiIndex.from_arrays([dac, n])).groupby(level=[0, 1]).size().items():
        sfs += c * hypergeom(nn, a, m).pmf(np.arange(m + 1))
    return sfs[1:m]  # drop monomorphic classes


_fcache = {}
BGRID = np.round(np.arange(-20, 20.0001, 0.05), 4)


def _f_exact(B, m):
    out = np.empty(m - 1)
    for i in range(1, m):
        if abs(B) < 1e-6:
            out[i - 1] = 1.0 / i
            continue
        g = lambda x: comb(m, i) * x ** i * (1 - x) ** (m - i) * (-np.expm1(-B * (1 - x))) / (x * (1 - x) * (-np.expm1(-B)))
        out[i - 1] = integrate.quad(g, 0, 1, limit=200)[0]
    return out


def f_i(B, m):
    """smooth in B: cubic interpolation of log f_i over a precomputed grid (so finite-difference gradients work)"""
    if m not in _fcache:
        from scipy.interpolate import CubicSpline
        tab = np.log(np.array([_f_exact(b, m) for b in BGRID]))
        _fcache[m] = CubicSpline(BGRID, tab, axis=0)
    return np.exp(_fcache[m](np.clip(B, BGRID[0], BGRID[-1])))


def negll(params, kN, kWS, kSW, m, fixB=None):
    j = 0
    B = fixB if fixB is not None else params[j]; j += fixB is None
    eps = 1 / (1 + np.exp(-params[j])); j += 1
    tN, tWS, tSW = np.exp(params[j:j + 3]); j += 3
    r = np.concatenate(([1.0], np.exp(params[j:])))
    i = np.arange(1, m)
    eN = tN * r / i
    fp, fm = f_i(B, m), f_i(-B, m)
    eWS = tWS * r * fp
    eSW = tSW * r * fm
    # polarisation error: observed class i <- true class m-i of the opposite type
    oN = (1 - eps) * eN + eps * eN[::-1]
    oWS = (1 - eps) * eWS + eps * eSW[::-1]
    oSW = (1 - eps) * eSW + eps * eWS[::-1]
    ll = 0.0
    for k, e in ((kN, oN), (kWS, oWS), (kSW, oSW)):
        e = np.maximum(e, 1e-300)
        ll += np.sum(k * np.log(e) - e - gammaln(k + 1))
    return -ll


def fit(kN, kWS, kSW, m):
    x0 = np.concatenate(([0.5, -4.0], np.log([kN.sum(), kWS.sum(), kSW.sum()]) - np.log(np.sum(1 / np.arange(1, m))), np.zeros(m - 2)))
    res = optimize.minimize(negll, x0, args=(kN, kWS, kSW, m), method="L-BFGS-B",
                            bounds=[(-20, 20), (-12, 0)] + [(None, None)] * 3 + [(-5, 5)] * (m - 2))
    Bhat, nll = res.x[0], res.fun
    # profile CI
    def prof(B):
        r = optimize.minimize(negll, res.x[1:], args=(kN, kWS, kSW, m, B), method="L-BFGS-B",
                              bounds=[(-12, 0)] + [(None, None)] * 3 + [(-5, 5)] * (m - 2))
        return r.fun - nll - 1.92
    lo = optimize.brentq(prof, Bhat - 15, Bhat) if prof(Bhat - 15) > 0 else np.nan
    hi = optimize.brentq(prof, Bhat, Bhat + 15) if prof(Bhat + 15) > 0 else np.nan
    return Bhat, lo, hi, 1 / (1 + np.exp(-res.x[1]))


if __name__ == "__main__":
    inp, out = sys.argv[1], sys.argv[2]
    m = int(sys.argv[3]) if len(sys.argv) > 3 else 40
    d = pd.read_csv(inp, sep="\t")
    d = d[(d.dac > 0) & (d.dac < d.n) & (d.n >= m)]
    rows = []
    for rc, g in d.groupby("region_class"):
        k = {c: project(g[g.cls.isin(cs)].dac.values, g[g.cls.isin(cs)].n.values, m)
             for c, cs in (("N", ["WW", "SS"]), ("WS", ["WS"]), ("SW", ["SW"]))}
        B, lo, hi, eps = fit(k["N"], k["WS"], k["SW"], m)
        rows.append(dict(region_class=rc, n_snps=len(g), B=B, B_lo=lo, B_hi=hi, eps=eps))
        print(rows[-1], flush=True)
    pd.DataFrame(rows).to_csv(out, sep="\t", index=False)
