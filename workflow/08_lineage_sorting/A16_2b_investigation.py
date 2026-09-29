"""A16: why does the 2b flank retain less ex-telomere signal than 2a? Four computational checks.

(a) Retained-fraction profile along each flank, 1-Mb steps to 30 Mb, vs real ends at matched distance
    (human branch, composition-corrected W/S log-odds; excess over interior).
(b) Data quality along the flanks: callable fraction and repeat share of substitutions near the junction,
    vs real ends at matched distance (is 2b sequence missing / unalignable near the junction?).
(c) Was the ancestral 2B end weak? Excess on OTHER lineages (gorilla, panstem, chimp, bonobo) at the hs1 windows
    of each flank, vs that lineage's own real ends at matched species distance (sp_dist_end), as a percentile.
    If the 2B orthologue is a weak end in gorilla (independent lineage) too, a low human 2b fraction need not
    mean an older fusion.
(d) ILS around the junction (A1 *.ils.tsv, orangutan outgroup): discordant fraction D in 0.1-Mb bins to 2 Mb
    and 1-Mb bins to 10 Mb, per side, and in genetic distance (deCODE sex-averaged CO, from A12 windows);
    symmetric/junction-centred (ancestral polymorphism) or one-sided/flat (substrate)?
usage: python A16_2b_investigation.py results/A3/windows.tsv.gz results/A1/chunks results/A12/A12_windows.tsv.gz results/A16
"""
import glob
import os
import sys

import numpy as np
import pandas as pd

win_f, chunk_dir, a12_f, OUT = sys.argv[1:5]
os.makedirs(OUT, exist_ok=True)
F = (113_940_058 + 114_049_496) / 2
ACRO = {"chr13", "chr14", "chr15", "chr21", "chr22"}
pd.set_option("display.width", 220)
d = pd.read_csv(win_f, sep="\t", dtype={"sp_name": str})
d = d[d.chrom.str.match(r"^chr\d+$")]
d["mid"] = d.win_start + 50_000
d["dfus"] = np.where(d.chrom == "chr2", (d.mid - F).abs(), np.inf)
d["side"] = np.where(d.dfus < 30e6, np.where(d.mid < F, "2a", "2b"), "")
d["near2"] = (d.chrom == "chr2") & (d.dfus < 30e6)
d["end_id"] = d.sp_name.astype(str) + "_" + d.sp_arm_end.astype(str)


def ex(x):
    return np.log(x.WS.sum() / x.SW.sum()) - np.log(x.parent_W.sum() / x.parent_S.sum())


# ---------- (a) profile along flanks (human) ----------
h = d[(d.branch == "human") & (d.n_callable > 20_000)]
inter = h[(h.sp_dist_end >= 30e6) & ~h.near2]
ends = h[(h.sp_dist_end < 30e6) & ~h.near2 & ~(h.chrom.isin(ACRO) & (h.sp_arm_end == "p"))]
base = ex(inter)
rows = []
for k in range(30):
    e = ends[(ends.sp_dist_end >= k * 1e6) & (ends.sp_dist_end < (k + 1) * 1e6)]
    ee = ex(e) - base
    for s in ("2a", "2b"):
        f = h[(h.side == s) & (h.dfus >= k * 1e6) & (h.dfus < (k + 1) * 1e6)]
        if len(f) >= 3 and (f.WS.sum() + f.SW.sum()) > 50:
            rows.append(dict(Mb=k, side=s, n_win=len(f), flank_excess=ex(f) - base, ends_excess=ee,
                             frac=(ex(f) - base) / ee, n_subs=int(f.WS.sum() + f.SW.sum())))
A = pd.DataFrame(rows)
A.to_csv(f"{OUT}/A16a_profile.tsv", sep="\t", index=False)
print("== (a) retained fraction along the flanks (1-Mb bins from the junction)")
print(A.pivot(index="Mb", columns="side", values="frac").round(2).T.to_string())

# ---------- (b) data quality near the junction ----------
hq = d[d.branch == "human"]
qrows = []
for k in range(10):
    e = hq[(hq.sp_dist_end >= k * 1e6) & (hq.sp_dist_end < (k + 1) * 1e6) & ~hq.near2]
    for s in ("2a", "2b"):
        f = hq[(hq.side == s) & (hq.dfus >= k * 1e6) & (hq.dfus < (k + 1) * 1e6)]
        qrows.append(dict(Mb=k, side=s, callable_frac=f.n_callable.sum() / (len(f) * 1e5) if len(f) else np.nan,
                          ends_callable_frac=e.n_callable.sum() / (len(e) * 1e5),
                          rep_share=(f.rep_WS.sum() + f.rep_SW.sum()) / max(f.WS.sum() + f.SW.sum(), 1),
                          ends_rep_share=(e.rep_WS.sum() + e.rep_SW.sum()) / (e.WS.sum() + e.SW.sum())))
Q = pd.DataFrame(qrows)
Q.to_csv(f"{OUT}/A16b_quality.tsv", sep="\t", index=False)
print("\n== (b) callable fraction / repeat share near the junction (flank vs real ends at same distance)")
print(Q.round(3).to_string(index=False))

# ---------- (c) orthologous end strength on other lineages ----------
crow = []
for br in ("gorilla", "panstem", "chimp", "bonobo"):
    b = d[(d.branch == br) & (d.n_callable > 20_000) & d.sp_dist_end.notna()]
    b_int = b[(b.sp_dist_end >= 30e6) & ~b.near2]
    bb = ex(b_int)
    b_ends = b[(b.sp_dist_end < 30e6) & ~b.near2]
    for s in ("2a", "2b"):
        for D in (5, 10):
            f = b[(b.side == s) & (b.dfus < D * 1e6)]
            if len(f) < 5:
                continue
            sd = f.sp_dist_end.median()
            # this lineage's real ends at matched species distance (+-2.5 Mb), per end
            per = []
            for e, g in b_ends[(b_ends.sp_dist_end - sd).abs() < 2.5e6].groupby("end_id"):
                if (g.WS.sum() + g.SW.sum()) >= 100:
                    per.append(ex(g) - bb)
            per = np.array(per)
            fx = ex(f) - bb
            crow.append(dict(lineage=br, side=s, first_Mb=D, sp_dist_Mb=sd / 1e6, sp_end=f.end_id.mode().iat[0],
                             flank_excess=fx, ends_matched_median=np.median(per) if len(per) else np.nan,
                             n_ends=len(per), pct_among_ends=(per < fx).mean() if len(per) else np.nan,
                             frac_of_matched=fx / np.median(per) if len(per) else np.nan))
C = pd.DataFrame(crow)
C.to_csv(f"{OUT}/A16c_orthologue_strength.tsv", sep="\t", index=False)
print("\n== (c) end strength of the flank orthologues on other lineages (vs that lineage's ends at matched distance)")
print(C.round(3).to_string(index=False))

# ---------- (d) ILS around the junction ----------
ils = pd.concat([pd.read_csv(f, sep="\t", dtype={"chrom": str}) for f in glob.glob(f"{chunk_dir}/chr2_*.ils.tsv")])
ils = ils[ils.chrom == "chr2"].copy()
ils["mid"] = ils.win_start + 50_000
ils["dfus"] = (ils.mid - F).abs()
ils["side"] = np.where(ils.mid < F, "2a", "2b")
allils = pd.concat([pd.read_csv(f, sep="\t", dtype={"chrom": str}) for f in glob.glob(f"{chunk_dir}/*.ils.tsv")])
Dg = (allils.HG_P.sum() + allils.PG_H.sum()) / (allils.HP_G.sum() + allils.HG_P.sum() + allils.PG_H.sum())
a12 = pd.read_csv(a12_f, sep="\t", usecols=["chrom", "win_start", "co_pat", "co_mat"])
a12 = a12[a12.chrom == "chr2"]
ils = ils.merge(a12, on=["chrom", "win_start"], how="left")
ils["cM_win"] = ((ils.co_pat + ils.co_mat) / 2).fillna(((ils.co_pat + ils.co_mat) / 2).median()) * 0.1
ils = ils.sort_values("mid")
for s, asc in (("2a", False), ("2b", True)):
    m = ils.side == s
    ils.loc[m, "cM_from_junction"] = ils[m].sort_values("mid", ascending=asc).cM_win.cumsum().reindex(ils[m].index)
drows = []
bins = [(i / 10, (i + 1) / 10) for i in range(20)] + [(k, k + 1) for k in range(2, 10)]
for s in ("2a", "2b"):
    for lo, hi in bins:
        x = ils[(ils.side == s) & (ils.dfus >= lo * 1e6) & (ils.dfus < hi * 1e6)]
        n = x.HP_G.sum() + x.HG_P.sum() + x.PG_H.sum()
        if n > 0:
            drows.append(dict(side=s, Mb_lo=lo, Mb_hi=hi, n_inf=int(n), D=(x.HG_P.sum() + x.PG_H.sum()) / n,
                              fold=((x.HG_P.sum() + x.PG_H.sum()) / n) / Dg, cM_mid=x.cM_from_junction.median()))
Dt = pd.DataFrame(drows)
Dt.to_csv(f"{OUT}/A16d_ils_junction.tsv", sep="\t", index=False)
print(f"\n== (d) ILS discordant fraction around the junction (genome D = {Dg:.3f})")
print(Dt.pivot_table(index=["Mb_lo", "Mb_hi"], columns="side", values=["fold", "cM_mid"]).round(2).to_string())
print("A16_DONE")
