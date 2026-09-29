"""A2b: clustered UBCS (T2T, B=2 bins, all SNDs incl. CpG as in 2022) at the fusion flanks.
Per branch: excess = mean UBCS/Mb in flank minus interior mean; rho = human excess / Pan excess
(same hs1 windows = orthologous). Also 2022-style human/chimp ratio. Paired 1-Mb block bootstrap."""
import glob, numpy as np, pandas as pd
rng = np.random.default_rng(2)
cols = ["chrom","begin","end","branch","B","cpg","p","exp","obs","UBCS","n"]
d = pd.concat([pd.read_csv(f, sep="\t", header=None, names=cols) for f in glob.glob("results/A2/*.B2.tsv")])
F = (113_940_058 + 114_049_496) / 2
sizes = pd.read_csv("data/ref/hs1.chrom.sizes", sep="\t", header=None, names=["chrom","len"]).set_index("chrom")["len"]
d["len"] = d.chrom.map(sizes); mid = d.begin + 5e5
d["dist_end"] = np.minimum(mid, d.len - mid)
w = d.pivot_table(index=["chrom","begin"], columns="branch", values="UBCS").dropna().reset_index()
w["mid"] = w.begin + 5e5
w["dist_end"] = np.minimum(w.mid, w.chrom.map(sizes) - w.mid)
interior = w[(w.dist_end > 30e6) & ~((w.chrom == "chr2") & ((w.mid - F).abs() < 30e6))]
base = interior[["human","chimp","bonobo"]].mean()
ends = w[w.dist_end < 5e6]
print("interior mean UBCS/Mb:", base.round(1).to_dict())
print("real ends (first 5 Mb) mean UBCS/Mb:", ends[["human","chimp","bonobo"]].mean().round(1).to_dict())
rows = []
for D in (5e6, 10e6, 15e6):
    fl = w[(w.chrom == "chr2") & ((w.mid - F).abs() < D)].copy()
    fl["side"] = np.where(fl.mid < F, "2a", "2b")
    for side in ["2a","2b","both"]:
        f = fl if side == "both" else fl[fl.side == side]
        def stat(x):
            eh = x.human.mean() - base.human
            out = {}
            for ref in ["chimp","bonobo"]:
                out[ref] = eh / (x[ref].mean() - base[ref])
            out["pan_mean"] = eh / (np.mean([x.chimp.mean() - base.chimp, x.bonobo.mean() - base.bonobo]))
            out["ratio2022"] = x.human.sum() / x.chimp.sum()
            return out
        pt = stat(f)
        bs = [stat(f.iloc[rng.integers(0, len(f), len(f))]) for _ in range(1000)]
        for k in pt:
            v = np.array([b[k] for b in bs]); lo, hi = np.percentile(v, [2.5, 97.5])
            rows.append(dict(D_Mb=D/1e6, side=side, stat=k, value=pt[k], lo=lo, hi=hi,
                             T_6_5=6.5*(1-pt[k]) if k!="ratio2022" else np.nan,
                             T_lo=6.5*(1-hi) if k!="ratio2022" else np.nan, T_hi=6.5*(1-lo) if k!="ratio2022" else np.nan,
                             human_flank=f.human.mean(), chimp_flank=f.chimp.mean(), bonobo_flank=f.bonobo.mean()))
r = pd.DataFrame(rows)
r.to_csv("results/A2/A2b_flank_ubcs.tsv", sep="\t", index=False)
pd.set_option("display.width", 220)
print(r.round(2).to_string(index=False))
