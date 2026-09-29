import sys, numpy as np
sys.path.insert(0, "scripts")
from A6_gbgc_sfs import f_i, fit
rng = np.random.default_rng(3)
m = 30
i = np.arange(1, m)
r = np.exp(rng.normal(0, 0.1, m - 1)); r[0] = 1
for Btrue in (0.0, 0.4, 1.0, 2.0):
    est = []
    for rep in range(3):
        kN = rng.poisson(20000 * r / i); kWS = rng.poisson(15000 * r * f_i(Btrue, m)); kSW = rng.poisson(15000 * r * f_i(-Btrue, m))
        B, lo, hi, eps = fit(kN.astype(float), kWS.astype(float), kSW.astype(float), m)
        est.append((round(B, 2), round(lo, 2), round(hi, 2), round(eps, 4)))
    print("B_true", Btrue, est)
