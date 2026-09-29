"""A8: is ILS elevated at ALL ancestral chromosome ends, or specifically at the fusion flanks?

Input: A1 *.ils.tsv (per 100-kb hs1 window: n_valid, HP_G, HG_P, PG_H site patterns; orangutan outgroup)
ILS proxy: discordant fraction D = (HG_P + PG_H) / (HP_G + HG_P + PG_H).
Classes (hs1 coordinates; human chromosome ends are ancestral ends except at the fusion):
  end      : window within 30 Mb of a hs1 chromosome end (acrocentric p-arms flagged separately)
  fusion   : hs1 chr2 within 30 Mb of the fusion site (2a proximal / 2b distal)
  interior : otherwise
Outputs: profile by 1-Mb distance bin (ends vs fusion sides), per-end values for the first 1/5 Mb,
percentile of the fusion flanks among ends; windows are weighted by informative sites.
Also joins deCODE recombination (hg38 windows lifted approximately by nearest hs1 position is NOT done here;
see A8b) — this script is purely positional.
usage: python A8_ils_profiles.py results/A1/chunks data/ref/hs1.chrom.sizes results/A8
"""
import glob
import os
import sys

import numpy as np
import pandas as pd

chunk_dir, sizes_f, outdir = sys.argv[1:4]
os.makedirs(outdir, exist_ok=True)
F = (113_940_058 + 114_049_496) / 2
ACRO = {"chr13", "chr14", "chr15", "chr21", "chr22"}
sizes = pd.read_csv(sizes_f, sep="\t", header=None, names=["chrom", "len"]).set_index("chrom")["len"]

d = pd.concat([pd.read_csv(f, sep="\t", dtype={"chrom": str}) for f in glob.glob(f"{chunk_dir}/*.ils.tsv")])
d = d[d.chrom.str.match(r"^chr\d+$")].groupby(["chrom", "win_start"], as_index=False).sum()
d["inf"] = d.HP_G + d.HG_P + d.PG_H
d = d[d.inf >= 50]
mid = d.win_start + 50_000
d["len"] = d.chrom.map(sizes)
d["dist_p"], d["dist_q"] = mid, d.len - mid
d["dist"] = np.minimum(d.dist_p, d.dist_q)
d["arm"] = np.where(d.dist_p < d.dist_q, "p", "q")
d["cls"] = np.where(d.dist < 30e6, "end", "interior")
d.loc[(d.cls == "end") & d.chrom.isin(ACRO) & (d.arm == "p"), "cls"] = "end_acro_p"
on_f = (d.chrom == "chr2") & ((mid - F).abs() < 30e6)
d.loc[on_f, "cls"] = np.where(mid[on_f] < F, "fusion_2a", "fusion_2b")
d.loc[on_f, "dist"] = (mid[on_f] - F).abs()
d["end_id"] = d.chrom + d.arm
d["dbin"] = (d.dist // 1e6).astype(int)


def disc(g):
    return (g.HG_P.sum() + g.PG_H.sum()) / g.inf.sum()


gm = disc(d)
prof = (d[d.cls != "interior"].groupby(["cls", "dbin"])
        .apply(lambda g: pd.Series(dict(D=disc(g), inf=g.inf.sum(), nwin=len(g)))).reset_index())
prof["fold_vs_genome"] = prof.D / gm
prof.to_csv(f"{outdir}/A8_profile.tsv", sep="\t", index=False)

rows = []
for K in (0.5, 1, 2, 5):
    sub = d[d.dist < K * 1e6]
    pe = sub[sub.cls == "end"].groupby("end_id").apply(lambda g: pd.Series(dict(D=disc(g), inf=g.inf.sum())))
    pe = pe[pe.inf >= 200]
    for side in ("fusion_2a", "fusion_2b"):
        f = sub[sub.cls == side]
        if f.empty:
            continue
        fv = disc(f)
        rows.append(dict(first_Mb=K, side=side, D_fusion=fv, fusion_inf=f.inf.sum(), D_ends_median=pe.D.median(),
                         D_ends_q10=pe.D.quantile(.1), D_ends_q90=pe.D.quantile(.9), n_ends=len(pe),
                         pct_among_ends=(pe.D < fv).mean(), D_genome=gm, D_interior=disc(d[d.cls == "interior"])))
    pe.to_csv(f"{outdir}/A8_per_end_first{K}Mb.tsv", sep="\t")
summ = pd.DataFrame(rows)
summ.to_csv(f"{outdir}/A8_summary.tsv", sep="\t", index=False)
pd.set_option("display.width", 200)
print(summ.round(3).to_string(index=False))
print(prof[prof.dbin < 6].round(3).to_string(index=False))
