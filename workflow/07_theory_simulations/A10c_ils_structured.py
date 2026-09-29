"""A10c: can lineage sorting at the fusion flanks date the fusion? Structured-coalescent predictions (msprime).

Three haploid samples: H (human, fused), C (chimp, unfused), G (gorilla). Discordance D = P(gene tree != ((H,C),G)).
Times (generations, g = 25 y): T_HC = 6.0 Mya, T_HG = 8.6 Mya (calibrations as in Yoo, Munson & Eichler 2026).
Ancestral HC effective size N_A calibrated so that the genome-wide D matches the observed 0.323 (A16d).

H_poly (ancestral fused/unfused polymorphism): between T_HC and T_HC + Delta the HC ancestor is structured into
  F (fused, size f N_A) and U (unfused, size (1-f) N_A); H descends from F, C from U.
  A lineage at genetic distance y (Morgans) from the junction switches arrangement background by recombination in
  heterozygotes: backward migration rates F->U = c*y*(1-f), U->F = c*y*f (c = recombination in heterozygotes
  relative to normal, <= 1: trivalent pairing may suppress it). Before T_HC + Delta: panmictic (size N_A).
H_sub (ancestral substrate): no structure; the local ancestral effective size is k * N_A along the whole region
  (high-recombination subtelomeres with weak linked selection); gives a flat, one-sided-compatible elevation.
Outputs: D(y) for a grid of Delta, c and y; D for H_sub as a function of k; the Delta and k needed to reproduce the
observed 2b elevation (fold ~1.3-1.8 over 1-10 Mb; ~2-2.6 in the first 0.2 Mb) and whether H_poly can produce a
flat profile over ~10 cM.
usage: python A10c_ils_structured.py results/A10 [reps=20000] [procs=8]
"""
import os
import sys
from multiprocessing import Pool

import msprime
import numpy as np
import pandas as pd

OUT = sys.argv[1]
REPS = int(sys.argv[2]) if len(sys.argv) > 2 else 20000
PROCS = int(sys.argv[3]) if len(sys.argv) > 3 else 8
os.makedirs(OUT, exist_ok=True)
GEN = 25
T_HC, T_HG = 6.0e6 / GEN, 8.6e6 / GEN
D_OBS = 0.323


def topo_discordant(ts):
    t = ts.first()
    # samples 0=H, 1=C, 2=G; concordant if MRCA(H,C) is younger than MRCA(H,G)
    return t.tmrca(0, 1) >= t.tmrca(0, 2) - 1e-9


def D_panmictic(NA, reps, seed):
    dem = msprime.Demography()
    for p in ("H", "C", "G", "HC", "HCG"):
        dem.add_population(name=p, initial_size=NA)
    dem.add_population_split(time=T_HC, derived=["H", "C"], ancestral="HC")
    dem.add_population_split(time=T_HG, derived=["HC", "G"], ancestral="HCG")
    it = msprime.sim_ancestry(samples={"H": 1, "C": 1, "G": 1}, ploidy=1, demography=dem, num_replicates=reps,
                              random_seed=seed)
    return np.mean([topo_discordant(ts) for ts in it])


def D_poly(args):
    NA, delta_gen, f, c, y, reps, seed = args
    dem = msprime.Demography()
    for p, n in (("H", NA), ("C", NA), ("G", NA), ("F", f * NA), ("U", (1 - f) * NA), ("HC", NA), ("HCG", NA)):
        dem.add_population(name=p, initial_size=n)
    dem.add_population_split(time=T_HC, derived=["H"], ancestral="F")
    dem.add_population_split(time=T_HC, derived=["C"], ancestral="U")
    m = c * y
    if m > 0:
        dem.add_migration_rate_change(time=T_HC, rate=m * (1 - f), source="F", dest="U")
        dem.add_migration_rate_change(time=T_HC, rate=m * f, source="U", dest="F")
    dem.add_population_split(time=T_HC + delta_gen, derived=["F", "U"], ancestral="HC")
    dem.add_migration_rate_change(time=T_HC + delta_gen, rate=0, source="F", dest="U")
    dem.add_migration_rate_change(time=T_HC + delta_gen, rate=0, source="U", dest="F")
    if T_HC + delta_gen < T_HG:
        dem.add_population_split(time=T_HG, derived=["HC", "G"], ancestral="HCG")
    else:  # polymorphism reaching back beyond the gorilla split: G joins the unstructured ancestor later
        dem.add_population_split(time=T_HC + delta_gen + 1, derived=["HC", "G"], ancestral="HCG")
    dem.sort_events()
    it = msprime.sim_ancestry(samples={"H": 1, "C": 1, "G": 1}, ploidy=1, demography=dem, num_replicates=reps,
                              random_seed=seed)
    return dict(Delta_Myr=delta_gen * GEN / 1e6, f=f, c=c, y_cM=y * 100, D=np.mean([topo_discordant(ts) for ts in it]))


if __name__ == "__main__":
    # 1. calibrate N_A to the genome-wide discordance (haploid samples: pair coalescence rate 1/N,
    #    so D = 2/3 exp(-dT / N_A) analytically; checked by simulation on a bracket)
    dT = T_HG - T_HC
    NA0 = dT / np.log((2 / 3) / D_OBS)
    cal = []
    for k, NA in enumerate(NA0 * np.array([0.8, 0.9, 1.0, 1.1, 1.25])):
        cal.append((NA, D_panmictic(NA, REPS, 100 + k)))
    NA = float(NA0)   # analytic; the simulated bracket above is reported as a check (D at NA0 should be ~D_OBS)
    print(f"N_A calibration: {[(round(a), round(b, 3)) for a, b in cal]} -> N_A = {NA:.0f}", flush=True)
    # 2. H_sub: local ancestral size k * N_A
    sub = [dict(k=k, D=D_panmictic(k * NA, REPS, 200 + i)) for i, k in enumerate((1.0, 1.2, 1.5, 2.0, 3.0))]
    S = pd.DataFrame(sub); S["fold"] = S.D / D_OBS
    S.to_csv(f"{OUT}/A10c_Hsub.tsv", sep="\t", index=False)
    print("H_sub (local N_A multiplier k):"); print(S.round(3).to_string(index=False), flush=True)
    # 3. H_poly grid
    grid = []
    seed = 1000
    for delta in (0.25, 0.5, 1.0, 2.0, 3.0):
        for f in (0.5,):
            for c in (1.0, 0.1, 0.01, 0.001):
                for y_cM in (0.0, 0.001, 0.01, 0.05, 0.1, 0.3, 1.0, 3.0, 10.0):
                    seed += 1
                    grid.append((NA, delta * 1e6 / GEN, f, c, y_cM / 100, REPS, seed))
    with Pool(PROCS) as pool:
        res = pool.map(D_poly, grid)
    Pz = pd.DataFrame(res); Pz["fold"] = Pz.D / D_OBS
    Pz.to_csv(f"{OUT}/A10c_Hpoly.tsv", sep="\t", index=False)
    print("H_poly: fold over genome-wide D, by distance from the junction (cM)")
    print(Pz.pivot_table(index=["Delta_Myr", "c"], columns="y_cM", values="fold").round(2).to_string())
    print("A10c_DONE")
