"""A14b: primary statistic (non-CpG W->S vs S->W log-odds, W/S-composition corrected) on the human branch,
fusion flanks vs real ends vs interior, computed identically for several A1 runs (assembly/haplotype cells).
Classes on hs1: flank 2a/2b (<= D Mb from the fusion), real ends (<= D Mb from a telomere; acrocentric p-arms and
chr2 within 30 Mb of the fusion excluded), interior (>= 30 Mb from ends and fusion).
frac = flank excess / real-end excess (both over interior); T = T_s (1 - frac), T_s = 6.0; 1-Mb block bootstrap.
usage: python A14_logodds_cells.py data/ref/hs1.chrom.sizes out.tsv name1=chunkdir1 name2=chunkdir2 ...
"""
import glob, sys
import numpy as np, pandas as pd
sizes_f, out = sys.argv[1:3]
cells = dict(a.split("=") for a in sys.argv[3:])
rng = np.random.default_rng(14)
F = (113_940_058 + 114_049_496) / 2
sizes = pd.read_csv(sizes_f, sep="\t", header=None, names=["chrom", "len"]).set_index("chrom")["len"]
ACRO = {"chr13", "chr14", "chr15", "chr21", "chr22"}


def load(cd):
    rows = []
    for f in glob.glob(f"{cd}/*.subs.tsv.gz"):
        s = pd.read_csv(f, sep="\t", usecols=["chrom", "pos", "branch", "cls", "cpg"])
        s = s[(s.branch == "human") & (s.cpg == 0) & s.cls.isin(["WS", "SW"])]
        s["mb"] = s.pos // 1_000_000
        rows.append(s.groupby(["chrom", "mb", "cls"]).size().unstack(fill_value=0))
    sub = pd.concat(rows).groupby(level=[0, 1]).sum()
    cal = pd.concat([pd.read_csv(f, sep="\t") for f in glob.glob(f"{cd}/*.callable.tsv")])
    cal = cal[cal.branch == "human"].assign(mb=lambda x: x.win_start // 1_000_000)
    cal = cal.groupby(["chrom", "mb"])[["parent_W", "parent_S"]].sum()
    d = sub.join(cal, how="inner").reset_index()
    d = d[d.chrom.str.match(r"^chr\d+$")]
    mid = d.mb * 1e6 + 5e5
    L = d.chrom.map(sizes)
    d["dp"], d["dq"] = mid, L - mid
    d["dfus"] = np.where(d.chrom == "chr2", (mid - F).abs(), np.inf)
    return d


def classify(d, D):
    near2 = (d.chrom == "chr2") & (d.dfus < 30e6)
    fl = np.where((d.chrom == "chr2") & (d.dfus < D * 1e6), np.where(d.mb * 1e6 + 5e5 < F, "2a", "2b"), "")
    endp = (d.dp < D * 1e6) & ~d.chrom.isin(ACRO)
    endq = d.dq < D * 1e6
    end = (endp | endq) & ~near2
    inter = (d.dp >= 30e6) & (d.dq >= 30e6) & ~near2
    return np.select([fl != "", end, inter], [fl, "end", "interior"], "x")


def ex(x):
    return np.log(x.WS.sum() / x.SW.sum()) - np.log(x.parent_W.sum() / x.parent_S.sum())


res = []
for name, cd in cells.items():
    d = load(cd)
    for D in (2, 5, 10):
        c = classify(d, D)
        g = {k: d[c == k] for k in ("2a", "2b", "end", "interior")}
        g["both"] = d[np.isin(c, ["2a", "2b"])]

        def stat(g):
            e_end = ex(g["end"]) - ex(g["interior"])
            return {s: (ex(g[s]) - ex(g["interior"])) / e_end for s in ("2a", "2b", "both")}, e_end
        est, e_end = stat(g)
        boot = []
        for _ in range(500):
            gb = {k: v.sample(len(v), replace=True, random_state=rng.integers(1e9)) for k, v in g.items()}
            boot.append(stat(gb)[0])
        for s in ("2a", "2b", "both"):
            bs = np.array([b[s] for b in boot])
            lo, hi = np.percentile(bs, [2.5, 97.5])
            res.append(dict(cell=name, first_Mb=D, side=s, end_excess=e_end, frac=est[s], lo=lo, hi=hi,
                            T=6 * (1 - est[s]), T_lo=6 * (1 - hi), T_hi=6 * (1 - lo)))
    print(name, "done", flush=True)
r = pd.DataFrame(res)
r.to_csv(out, sep="\t", index=False)
pd.set_option("display.width", 200)
print(r.round(3).to_string(index=False))
