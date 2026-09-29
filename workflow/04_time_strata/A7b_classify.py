"""A7b: time-stratified human-branch substitutions (Anc4 -> hs1) using archaic hominins and 1KGP-CHM13.

Strata (derived = CHM13 base on the human branch):
  poly         : modern derived allele frequency (all 2504) < 0.95        -> recent / still segregating
  fixed_shared : DAF >= 0.95 and every called archaic genome carries the derived allele (>= 1 Neandertal
                 and Denisova called)                                     -> fixed mostly before the modern-archaic split
  fixed_human  : DAF >= 0.95 and every called archaic genome is homozygous ancestral (same calling
                 requirement)                                             -> fixed after the split (human-specific)
  other        : everything else (mixed archaic states, too few calls)
Sites absent from the 1KGP VCF are monomorphic in 1KGP, i.e. fixed derived (DAF = 1).
Region classes (hs1): fusion_2a / fusion_2b (<= 10 Mb from the fusion site), end (<= 10 Mb from a real end,
acrocentric p-arms excluded), interior (>= 30 Mb from ends and fusion).
Metric: W->S share = WS / (WS + SW) among non-CpG substitutions; 1-Mb block bootstrap CI.
usage: python A7b_classify.py results/A7 data/ref/hs1.chrom.sizes
"""
import glob
import sys

import numpy as np
import pandas as pd

O, sizes_f = sys.argv[1:3]
rng = np.random.default_rng(7)
F = (113_940_058 + 114_049_496) / 2
ACRO = {"chr13", "chr14", "chr15", "chr21", "chr22"}
COMP = str.maketrans("ACGT", "TGCA")

# sites and their hg19 placement
b = pd.read_csv(f"{O}/human_subs.hg19.bed", sep="\t", header=None, names=["c19", "s19", "e19", "name", "sc", "strand"])
f = b.name.str.split("|", expand=True)
b["id"], b["anc"], b["der"], b["cls"], b["cpg"] = f[0], f[1], f[2], f[3], f[4].astype(int)
b["chrom"] = b.id.str.split(":").str[0]
b["pos"] = b.id.str.split(":").str[1].astype(int) + 1          # hs1 1-based
b["der19"] = np.where(b.strand == "-", b.der.str.translate(COMP), b.der)
b["key19"] = b.c19.str.replace("chr", "", regex=False) + ":" + b.e19.astype(str)

# archaic genotypes
arch = {}
for fn in glob.glob(f"{O}/archaic/*.tsv"):
    g = fn.split("/")[-1].split(".")[0]
    a = pd.read_csv(fn, sep="\t", header=None, names=["c", "p", "ref", "alt", "gt", "gq"], dtype=str)
    a["key19"] = a.c + ":" + a.p
    arch.setdefault(g, []).append(a)
state = {}
for g, parts in arch.items():
    a = pd.concat(parts).drop_duplicates("key19").set_index("key19")
    m = b[["key19", "der19"]].join(a, on="key19", how="left")
    gt = m["gt"].fillna("./.")
    called = ~gt.str.contains(r"\.", regex=True) & (pd.to_numeric(m["gq"], errors="coerce").fillna(0) >= 20)
    alleles = [m.ref.fillna("N"), m.alt.fillna(".")]
    g1 = gt.str[0].map({"0": 0, "1": 1}).fillna(-1).astype(int)
    g2 = gt.str[-1].map({"0": 0, "1": 1}).fillna(-1).astype(int)
    base1 = np.where(g1 == 0, alleles[0], np.where(g1 == 1, alleles[1], "N"))
    base2 = np.where(g2 == 0, alleles[0], np.where(g2 == 1, alleles[1], "N"))
    nder = (base1 == m.der19.values).astype(int) + (base2 == m.der19.values).astype(int)
    state[g] = pd.DataFrame({"called": called.values, "nder": np.where(called, nder, -1)}, index=b.index)
nean = [g for g in state if g != "Denisova"]
called_n = sum(state[g].called for g in nean)
called_d = state["Denisova"].called if "Denisova" in state else pd.Series(False, index=b.index)
all_called = [state[g] for g in state]
n_called = sum(s.called.astype(int) for s in all_called)
n_carry = sum(((s.nder > 0) & s.called).astype(int) for s in all_called)
n_homanc = sum(((s.nder == 0) & s.called).astype(int) for s in all_called)
ok_calls = (called_n > 0) & called_d

# modern DAF
k = pd.concat([pd.read_csv(fn, sep="\t", header=None, names=["chrom", "pos", "ref", "alt", "an", "ac", "an_afr", "ac_afr"])
               for fn in glob.glob(f"{O}/kgp/*.tsv")]).drop_duplicates(["chrom", "pos"])
m = b[["chrom", "pos", "der"]].merge(k, on=["chrom", "pos"], how="left")
daf = np.where(m.alt == m.der, m.ac / m.an, np.where(m.ref == m.der, (m.an - m.ac) / m.an, np.nan))
daf = np.where(m.an.isna(), 1.0, daf)          # not in VCF -> monomorphic -> fixed derived
b["daf"] = daf
b["stratum"] = "other"
b.loc[b.daf < 0.95, "stratum"] = "poly"
fixed = b.daf >= 0.95
b.loc[fixed & ok_calls & (n_carry == n_called), "stratum"] = "fixed_shared"
b.loc[fixed & ok_calls & (n_homanc == n_called), "stratum"] = "fixed_human"

# region classes on hs1
sizes = pd.read_csv(sizes_f, sep="\t", header=None, names=["chrom", "len"]).set_index("chrom")["len"]
b["len"] = b.chrom.map(sizes)
dist = np.minimum(b.pos, b.len - b.pos)
arm_p = b.pos < b.len / 2
b["region"] = np.where(dist <= 10e6, "end", np.where(dist >= 30e6, "interior", "mid"))
b.loc[(b.region == "end") & b.chrom.isin(ACRO) & arm_p, "region"] = "acro_p"
nf = (b.chrom == "chr2") & ((b.pos - F).abs() <= 10e6)
b.loc[nf, "region"] = np.where(b.pos[nf] < F, "fusion_2a", "fusion_2b")
b.loc[(b.chrom == "chr2") & ((b.pos - F).abs() < 30e6) & ~nf & (b.region == "interior"), "region"] = "mid"
b["blk"] = b.chrom + ":" + (b.pos // 1_000_000).astype(str)
b.to_csv(f"{O}/A7_sites_classified.tsv.gz", sep="\t", index=False,
         columns=["chrom", "pos", "anc", "der", "cls", "cpg", "daf", "stratum", "region"])

x = b[(b.cpg == 0) & b.cls.isin(["WS", "SW"])]
rows = []
for st in ["fixed_shared", "fixed_human", "poly", "other", "ALL"]:
    xs = x if st == "ALL" else x[x.stratum == st]
    for rg in ["fusion_2a", "fusion_2b", "end", "interior"]:
        xr = xs[xs.region == rg]
        n = len(xr)
        if n == 0:
            continue
        share = (xr.cls == "WS").mean()
        blks = xr.blk.unique()
        grp = xr.groupby("blk").cls.agg(lambda v: ((v == "WS").sum(), len(v)))
        ws = np.array([g[0] for g in grp]); tot = np.array([g[1] for g in grp])
        bs = []
        for _ in range(500):
            i = rng.integers(0, len(ws), len(ws))
            bs.append(ws[i].sum() / tot[i].sum())
        lo, hi = np.percentile(bs, [2.5, 97.5])
        rows.append(dict(stratum=st, region=rg, n=n, WS_share=share, lo=lo, hi=hi))
r = pd.DataFrame(rows)
r.to_csv(f"{O}/A7_summary.tsv", sep="\t", index=False)
pd.set_option("display.width", 200)
print(b.stratum.value_counts().to_string())
print(r.round(4).to_string(index=False))
