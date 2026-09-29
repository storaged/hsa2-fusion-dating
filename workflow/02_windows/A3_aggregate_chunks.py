"""A3a: aggregate A1 chunk outputs into per-(branch, 100-kb hs1 window) tables with species placement.

For every branch and hs1 100-kb window:
  counts of WS/SW/WW/SS substitutions, split by parent-CpG context, parsimony support and repeat flag
  callable sites and parent W / S counts (from A1 .callable.tsv)
  species placement: modal species chromosome and median species coordinate of the substitutions,
  the species chromosome length, and the distance (bp) to the nearest end of that chromosome.
Human and panstem branches are placed on hs1 (panstem: chimp placement is added from the chimp branch
of the same window when available, since the Pan ancestor karyotype = chimp for these purposes).

usage: python A3_aggregate_chunks.py results/A1/chunks data/asm_reports results/A3/windows.tsv.gz
"""
import glob
import os
import sys

import numpy as np
import pandas as pd

chunk_dir, rep_dir, out = sys.argv[1:4]
WIN = 100_000

# accession -> species chromosome name
names = {}
for f in glob.glob(f"{rep_dir}/*.report.txt"):
    for line in open(f):
        if line.startswith("#"):
            continue
        x = line.rstrip("\n").split("\t")
        if len(x) > 4 and x[1] == "assembled-molecule" and x[4] != "na":
            names[x[4]] = x[0]

frames = []
for sf in sorted(glob.glob(f"{chunk_dir}/*.subs.tsv.gz")):
    cf = sf.replace(".subs.tsv.gz", ".callable.tsv")
    if not os.path.exists(cf):
        continue
    s = pd.read_csv(sf, sep="\t", dtype={"chrom": str, "sp_chrom": str})
    c = pd.read_csv(cf, sep="\t", dtype={"chrom": str})
    if c.empty:
        continue
    s["win_start"] = s.pos // WIN * WIN
    s["cls_cpg"] = s.cls + np.where(s.cpg == 1, "_cpg", "")
    cnt = (s.pivot_table(index=["chrom", "win_start", "branch"], columns="cls_cpg", values="pos",
                         aggfunc="size", fill_value=0))
    pars = s[s.pars == 1].pivot_table(index=["chrom", "win_start", "branch"], columns="cls", values="pos",
                                      aggfunc="size", fill_value=0).add_prefix("pars_")
    rep = s[s.repeat == 1].pivot_table(index=["chrom", "win_start", "branch"], columns="cls", values="pos",
                                       aggfunc="size", fill_value=0).add_prefix("rep_")
    place = (s.groupby(["chrom", "win_start", "branch"])
              .agg(sp_chrom=("sp_chrom", lambda v: v.mode().iat[0]), sp_pos=("sp_pos", "median"),
                   sp_len=("sp_len", "max"), sp_chrom_frac=("sp_chrom", lambda v: (v == v.mode().iat[0]).mean())))
    t = c.set_index(["chrom", "win_start", "branch"]).join([cnt, pars, rep, place], how="left")
    frames.append(t.reset_index())

d = pd.concat(frames, ignore_index=True).fillna({k: 0 for k in []})
num = [col for col in d.columns if col.split("_")[0] in ("WS", "SW", "WW", "SS", "pars", "rep")]
d[num] = d[num].fillna(0).astype(int)
d["sp_name"] = d.sp_chrom.map(names).fillna(d.sp_chrom)
d["sp_dist_end"] = np.minimum(d.sp_pos, d.sp_len - d.sp_pos)
d["sp_arm_end"] = np.where(d.sp_pos < d.sp_len / 2, "p", "q")
os.makedirs(os.path.dirname(out), exist_ok=True)
d.to_csv(out, sep="\t", index=False, compression="gzip")
print(d.shape)
print(d.groupby("branch")[[c for c in ["n_callable", "WS", "SW", "WS_cpg", "SW_cpg"] if c in d.columns]].sum())
