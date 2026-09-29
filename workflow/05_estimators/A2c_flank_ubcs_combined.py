"""A2c: fusion-flank UBCS, human lineage (Anc4->hs1) vs full Pan lineages (Anc4->chimp = panstem+chimp,
Anc4->bonobo = panstem+bonobo), same hs1 1-Mb windows (orthologous).  T2T data, B = 2 bins, all SNDs.

Estimators (per flank side and flank width D):
  rho_excess = (U_H,flank - U_H,interior) / (U_P,flank - U_P,interior)
  rho_2022   = sum U_H,flank / sum U_P,flank, divided by the mean of the same ratio over real ends (legacy form)
  T = T_eff * (1 - rho), T_eff in {6.0, 6.5, 7.0}
Also reports the Pan-internal temporal change: panstem-only vs recent (chimp/bonobo) UBCS per Myr.
Paired 1-Mb block bootstrap (1000).
"""
import glob

import numpy as np
import pandas as pd

import sys
rng = np.random.default_rng(3)
# optional args: <ubcs_dir> <chrom.sizes> <fusion_mid> <out_tsv>  (defaults: CHM13 run)
A2DIR = sys.argv[1] if len(sys.argv) > 1 else "results/A2"
SIZES = sys.argv[2] if len(sys.argv) > 2 else "data/ref/hs1.chrom.sizes"
F = float(sys.argv[3]) if len(sys.argv) > 3 else (113_940_058 + 114_049_496) / 2
OUT = sys.argv[4] if len(sys.argv) > 4 else "results/A2/A2c_flank_ubcs_combined.tsv"
cols = ["chrom", "begin", "end", "branch", "B", "cpg", "p", "exp", "obs", "UBCS", "n"]
d = pd.concat([pd.read_csv(f, sep="\t", header=None, names=cols) for f in glob.glob(f"{A2DIR}/*.B2.tsv")])
sizes = pd.read_csv(SIZES, sep="\t", header=None, names=["chrom", "len"]).set_index("chrom")["len"]
w = d.pivot_table(index=["chrom", "begin"], columns="branch", values="UBCS").reset_index()
need = ["human", "panstem+chimp", "panstem+bonobo"]
w = w.dropna(subset=need)
w["mid"] = w.begin + 5e5
w["dist_end"] = np.minimum(w.mid, w.chrom.map(sizes) - w.mid)
near_f = (w.chrom == "chr2") & ((w.mid - F).abs() < 30e6)
interior = w[(w.dist_end > 30e6) & ~near_f]
base = interior[need].mean()
ACRO = {"chr13", "chr14", "chr15", "chr21", "chr22"}
w["arm"] = np.where(w.mid < w.chrom.map(sizes) / 2, "p", "q")

print("interior mean UBCS/Mb:", base.round(1).to_dict())
rows = []
for D in (5e6, 10e6, 15e6):
    ends = w[(w.dist_end < D) & ~near_f & ~(w.chrom.isin(ACRO) & (w.arm == "p"))].copy()
    ends["end_id"] = ends.chrom + ends.arm
    fl = w[(w.chrom == "chr2") & ((w.mid - F).abs() < D)].copy()
    fl["side"] = np.where(fl.mid < F, "2a", "2b")
    for ref in ["panstem+chimp", "panstem+bonobo"]:
        ctl = ends.groupby("end_id").apply(lambda x: x.human.sum() / x[ref].sum())
        for side in ["2a", "2b", "both"]:
            f = fl if side == "both" else fl[fl.side == side]

            def stat(x, c):
                ex = (x.human.mean() - base.human) / (x[ref].mean() - base[ref])
                r22 = (x.human.sum() / x[ref].sum()) / c.mean()
                return ex, r22
            pt = stat(f, ctl)
            bs = []
            for _ in range(1000):
                xb = f.iloc[rng.integers(0, len(f), len(f))]
                cb = ctl.iloc[rng.integers(0, len(ctl), len(ctl))]
                bs.append(stat(xb, cb))
            bs = np.array(bs)
            for j, nm in enumerate(["rho_excess", "rho_2022"]):
                lo, hi = np.percentile(bs[:, j], [2.5, 97.5])
                rows.append(dict(D_Mb=D / 1e6, ref=ref, side=side, stat=nm, rho=pt[j], lo=lo, hi=hi,
                                 T65=6.5 * (1 - pt[j]), T65_lo=6.5 * (1 - hi), T65_hi=6.5 * (1 - lo),
                                 U_human=f.human.mean(), U_ref=f[ref].mean(), ends_ratio_mean=ctl.mean()))
r = pd.DataFrame(rows)
r.to_csv(OUT, sep="\t", index=False)
pd.set_option("display.width", 230)
print(r.round(2).to_string(index=False))

# temporal change inside Pan at the flanks: UBCS per Mb on stem vs recent branches, relative to interior
if "panstem" in w.columns:
    print("\nPan-internal (flank 0-10 Mb, excess over interior):")
    fl = w[(w.chrom == "chr2") & ((w.mid - F).abs() < 10e6)].dropna(subset=["panstem", "chimp", "bonobo"])
    ib = interior.dropna(subset=["panstem", "chimp", "bonobo"])[["panstem", "chimp", "bonobo"]].mean()
    for side, x in fl.groupby(np.where(fl.mid < F, "2a", "2b")):
        print(side, {b: round(x[b].mean() - ib[b], 1) for b in ["panstem", "chimp", "bonobo"]})
