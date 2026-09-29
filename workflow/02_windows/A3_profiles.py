"""A3b: telomere-distance profiles of W->S / S->W substitution rates per branch, and the
retained telomere-type signal R at the human fusion flanks.

Metrics (non-CpG only):
  rWS = WS / parent_W     rSW = SW / parent_S     gBGC index g = rWS / rSW
Distances: distance to the nearest end of the chromosome in the genome where the branch ends
(species coordinate). For the human branch, windows on hs1 chr2 within 30 Mb of the fusion site
are classified separately as ex-ends ('fusion_2a' proximal, 'fusion_2b' distal).
R = (value at fusion flank - interior) / (value at real ends - interior), distance-matched.
usage: python A3_profiles.py results/A3/windows.tsv.gz results/A3
"""
import sys

import numpy as np
import pandas as pd

inp, outdir = sys.argv[1:3]
FUSION_HS1 = (113_940_058 + 114_049_496) / 2
MAXD = 30e6
d = pd.read_csv(inp, sep="\t", dtype={"chrom": str, "sp_chrom": str, "sp_name": str})
d = d[d.chrom.str.match(r"^chr\d+$")]                       # autosomes (hs1)
d = d[(d.n_callable > 20_000) & d.sp_dist_end.notna()]
d["WS"] = d["WS"]; d["SW"] = d["SW"]

# classification
d["dist"] = d.sp_dist_end
d["cls"] = np.where(d.dist < MAXD, "end", "interior")
d["end_id"] = d.sp_name.astype(str) + "_" + d.sp_arm_end
hum = d.branch.isin(["human"])
on2 = hum & (d.chrom == "chr2") & ((d.win_start + 50_000 - FUSION_HS1).abs() < MAXD)
d.loc[on2, "cls"] = np.where(d.loc[on2, "win_start"] + 50_000 < FUSION_HS1, "fusion_2a", "fusion_2b")
d.loc[on2, "dist"] = (d.loc[on2, "win_start"] + 50_000 - FUSION_HS1).abs()
# acrocentric p-arm ends (hs1 13,14,15,21,22 p) are rDNA/satellite-bounded: flag
acro = d.sp_name.str.contains(r"hsa(13|14|15|21|22)$|^chr(13|14|15|21|22)$", regex=True) & (d.sp_arm_end == "p")
d.loc[acro & (d.cls == "end"), "cls"] = "end_acro_p"
d["dbin"] = (d.dist // 1e6).astype(int)

def rates(g):
    return pd.Series(dict(WS=g.WS.sum(), SW=g.SW.sum(), pW=g.parent_W.sum(), pS=g.parent_S.sum(),
                          n=g.n_callable.sum(), nwin=len(g)))

prof = d[d.cls.isin(["end", "fusion_2a", "fusion_2b"])].groupby(["branch", "cls", "dbin"]).apply(rates).reset_index()
inter = d[d.cls == "interior"].groupby("branch").apply(rates).reset_index()
for t in (prof, inter):
    t["rWS"] = t.WS / t.pW; t["rSW"] = t.SW / t.pS; t["g"] = t.rWS / t.rSW
prof.to_csv(f"{outdir}/A3_profile_by_distance.tsv", sep="\t", index=False)
inter.to_csv(f"{outdir}/A3_interior.tsv", sep="\t", index=False)

# per-end aggregates for the first K Mb (for between-end variability)
rows = []
for K in (2, 5, 10, 15):
    sub = d[d.dist < K * 1e6]
    per_end = sub[sub.cls == "end"].groupby(["branch", "end_id"]).apply(rates).reset_index()
    per_end["rWS"] = per_end.WS / per_end.pW; per_end["g"] = per_end.rWS / (per_end.SW / per_end.pS)
    per_end = per_end[per_end.n > 0.3 * K * 1e6]  # need reasonable callable coverage
    for br in per_end.branch.unique():
        pe = per_end[per_end.branch == br]
        it = inter[inter.branch == br].iloc[0]
        base = dict(rWS=it.rWS, g=it.g)
        for side in ("fusion_2a", "fusion_2b"):
            f = sub[(sub.branch == br) & (sub.cls == side)]
            if f.empty:
                continue
            fr = rates(f)
            fv = dict(rWS=fr.WS / fr.pW, g=(fr.WS / fr.pW) / (fr.SW / fr.pS))
            for m in ("rWS", "g"):
                rows.append(dict(branch=br, first_Mb=K, side=side, metric=m, fusion=fv[m],
                                 ends_mean=pe[m].mean(), ends_sd=pe[m].std(), n_ends=len(pe),
                                 interior=base[m], R=(fv[m] - base[m]) / (pe[m].mean() - base[m]),
                                 pct_among_ends=(pe[m] < fv[m]).mean()))
    per_end.to_csv(f"{outdir}/A3_per_end_first{K}Mb.tsv", sep="\t", index=False)
summ = pd.DataFrame(rows)
summ.to_csv(f"{outdir}/A3_fusion_summary.tsv", sep="\t", index=False)
pd.set_option("display.width", 220)
print(inter.round(5).to_string(index=False))
print(summ.round(4).to_string(index=False))
