"""A7a: human-branch substitutions (8-way T2T A1, Anc4 -> hs1) -> target files for the archaic and 1KGP queries.

Outputs (results/A7/):
  human_subs.hs1.bed      hs1 0-based BED, name = id|anc|der|cls|cpg|pars, strand '+'
  hs1 targets per chromosome for 1KGP:   targets_hs1/<chrom>.tsv  (chrom, pos 1-based)
After liftOver (A7_run.sh): human_subs.hg19.bed -> targets_hg19/<N>.tsv (no 'chr', archaic naming)
usage: python A7a_sites.py results/A1/chunks results/A7
"""
import glob
import os
import sys

import pandas as pd

chunks, out = sys.argv[1:3]
os.makedirs(f"{out}/targets_hs1", exist_ok=True)
fr = []
for f in glob.glob(f"{chunks}/*.subs.tsv.gz"):
    s = pd.read_csv(f, sep="\t", usecols=["chrom", "pos", "branch", "anc", "der", "cls", "cpg", "pars"], dtype={"chrom": str})
    fr.append(s[s.branch == "human"])
s = pd.concat(fr).drop_duplicates(["chrom", "pos"]).sort_values(["chrom", "pos"])
s = s[s.chrom.str.match(r"^chr(\d+|X)$")]
s["name"] = (s.chrom + ":" + s.pos.astype(str) + "|" + s.anc + "|" + s.der + "|" + s.cls + "|" + s.cpg.astype(str)
             + "|" + s.pars.astype(str))
bed = pd.DataFrame({"c": s.chrom, "s": s.pos, "e": s.pos + 1, "n": s.name, "sc": 0, "st": "+"})
bed.to_csv(f"{out}/human_subs.hs1.bed", sep="\t", header=False, index=False)
for c, g in s.groupby("chrom"):
    pd.DataFrame({"c": g.chrom, "p": g.pos + 1}).to_csv(f"{out}/targets_hs1/{c}.tsv", sep="\t", header=False, index=False)
print(len(s), "human-branch substitutions;", s.cls.value_counts().to_dict())
