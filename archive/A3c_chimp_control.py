"""A3c: chimp/bonobo/panstem branch values at the loci orthologous to the human fusion flanks,
compared with the same branch's real chromosome ends at the same species-distance to the end."""
import numpy as np, pandas as pd
d = pd.read_csv("results/A3/windows.tsv.gz", sep="\t", dtype={"chrom": str, "sp_chrom": str, "sp_name": str})
d = d[(d.n_callable > 20000) & d.sp_dist_end.notna() & d.chrom.str.match(r"^chr\d+$")]
F = (113_940_058 + 114_049_496) / 2
rows = []
for br in ["chimp", "bonobo", "gorilla", "panstem", "human"]:
    x = d[d.branch == br]
    if br == "panstem":   # place on chimp coordinates
        pl = d[d.branch == "chimp"].set_index(["chrom", "win_start"])[["sp_name", "sp_dist_end"]]
        x = x.drop(columns=["sp_name", "sp_dist_end"]).join(pl, on=["chrom", "win_start"]).dropna(subset=["sp_dist_end"])
    x = x.copy(); x["dist"] = x.sp_dist_end
    fl = x[(x.chrom == "chr2") & ((x.win_start + 5e4 - F).abs() < 15e6)].copy()
    fl["side"] = np.where(fl.win_start + 5e4 < F, "2a", "2b")
    inter = x[x.dist > 30e6]
    base = (inter.WS.sum() / inter.parent_W.sum()) / (inter.SW.sum() / inter.parent_S.sum())
    for side, g in fl.groupby("side"):
        med_d = g.dist.median()
        # real ends of this branch (excluding chr2 flanks) at matched distance band
        lo, hi = g.dist.quantile(.1), g.dist.quantile(.9)
        e = x[(x.dist >= lo) & (x.dist <= hi) & ~((x.chrom == "chr2") & ((x.win_start + 5e4 - F).abs() < 15e6))]
        ge = (e.WS.sum() / e.parent_W.sum()) / (e.SW.sum() / e.parent_S.sum())
        gf = (g.WS.sum() / g.parent_W.sum()) / (g.SW.sum() / g.parent_S.sum())
        rows.append(dict(branch=br, side=side, sp_chrom=g.sp_name.mode().iat[0], median_dist_to_end_Mb=med_d / 1e6,
                         g_flank=gf, g_matched_ends=ge, g_interior=base, R=(gf - base) / (ge - base)))
print(pd.DataFrame(rows).round(3).to_string(index=False))
