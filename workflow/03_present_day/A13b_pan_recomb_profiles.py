"""A13b: present-day (LD-based) recombination in Pan/gorilla at the orthologues of the human fusion flanks,
versus that species' own real chromosome ends and interior.

Inputs: results/A13/pan_rho.<sp>.hs1_100kb.tsv (A13a; rho/kb on hs1 100-kb windows)
        results/A3/windows.tsv.gz (branch = <sp>: distance of each hs1 window to the physical end in the species'
        own T2T karyotype, sp_dist_end, and which end, sp_name + sp_arm_end)
Metrics: rho normalised by the species' interior median (windows >= 30 Mb from any end).
  * profile over physical distance to the telomere (1-Mb bins) at real ends;
  * orthologues of the human fusion flanks (hs1 chr2 within 10 Mb of the fusion), per side;
  * cap proxy: for each species end, the distance from the physical end to the first window orthologous to hs1
    (the unaligned terminal block, i.e. cap + unalignable subtelomere). Ends with a block > 3 Mb are 'capped'.
    Profiles are also computed over the distance from that boundary, to ask whether end-type recombination
    follows the physical telomere (pushed away by caps) or the cap boundary.
usage: python A13b_pan_recomb_profiles.py <species> results/A13 results/A3/windows.tsv.gz
"""
import sys
import numpy as np
import pandas as pd

sp, OUT, win_f = sys.argv[1:4]
FUSION = (113_940_058 + 114_049_496) / 2
rho = pd.read_csv(f"{OUT}/pan_rho.{sp}.hs1_100kb.tsv", sep="\t")
w = pd.read_csv(win_f, sep="\t", usecols=["chrom", "win_start", "branch", "n_callable", "sp_name", "sp_dist_end", "sp_arm_end"],
                dtype={"sp_name": str})
w = w[(w.branch == sp) & w.chrom.str.match(r"^chr\d+$") & w.sp_dist_end.notna()]
d = w.merge(rho, on=["chrom", "win_start"], how="inner")
d = d[d.cov_bp > 20_000].copy()
d["end_id"] = d.sp_name + "_" + d.sp_arm_end
mid = d.win_start + 50_000
d["flank"] = np.where((d.chrom == "chr2") & ((mid - FUSION).abs() < 10e6), np.where(mid < FUSION, "fl_2a", "fl_2b"), "")
d["d_fus"] = np.where(d.flank != "", (mid - FUSION).abs(), np.nan)
inter = d[(d.sp_dist_end >= 30e6) & (d.flank == "")]
norm = inter.rho_kb.median()
d["r"] = d.rho_kb / norm
# cap proxy: unaligned terminal block per species end
first = d.groupby("end_id").sp_dist_end.min().rename("term_block")
d = d.join(first, on="end_id")
d["d_boundary"] = d.sp_dist_end - d.term_block
ends = d[(d.sp_dist_end < 30e6) & (d.flank == "")]
capped_ids = first[first > 3e6].index

print(f"== {sp}: {len(d)} windows with rho; interior median rho/kb {norm:.4f}; "
      f"{ends.end_id.nunique()} ends ({len(set(capped_ids) & set(ends.end_id))} with terminal unaligned block > 3 Mb)")
prof = []
for lab, sub, col in (("all_ends_phys", ends, "sp_dist_end"),
                      ("uncapped_phys", ends[~ends.end_id.isin(capped_ids)], "sp_dist_end"),
                      ("capped_phys", ends[ends.end_id.isin(capped_ids)], "sp_dist_end"),
                      ("capped_from_boundary", ends[ends.end_id.isin(capped_ids)], "d_boundary")):
    for lo, hi in ((0, 2), (2, 5), (5, 10), (10, 15), (15, 20), (20, 25), (25, 30)):
        x = sub[(sub[col] >= lo * 1e6) & (sub[col] < hi * 1e6)]
        if len(x) >= 5:
            prof.append(dict(set=lab, bin_Mb=f"{lo}-{hi}", n_win=len(x), n_ends=x.end_id.nunique(),
                             median_r=x.r.median(), mean_r=x.r.mean()))
prof = pd.DataFrame(prof)
fl = []
for side in ("fl_2a", "fl_2b"):
    x = d[d.flank == side]
    for lo, hi in ((0, 2), (0, 5), (0, 10)):
        y = x[(x.d_fus >= lo * 1e6) & (x.d_fus < hi * 1e6)]
        if len(y):
            # matched real-end reference at the same physical distance in this species
            dd = y.sp_dist_end.median()
            ref = ends[(ends.sp_dist_end - dd).abs() < 2.5e6]
            near = ends[ends.sp_dist_end < 2e6]
            fl.append(dict(side=side, human_dist_Mb=f"{lo}-{hi}", n_win=len(y), median_r=y.r.median(), mean_r=y.r.mean(),
                           sp_dist_end_Mb=dd / 1e6, sp_end=y.end_id.mode().iat[0], term_block_Mb=y.term_block.median() / 1e6,
                           ends_same_dist_median_r=ref.r.median(), ends_0_2Mb_median_r=near.r.median(),
                           pct_among_same_dist_windows=(ref.r < y.r.median()).mean()))
fl = pd.DataFrame(fl)
pd.set_option("display.width", 220)
print(prof.round(3).to_string(index=False))
print(fl.round(3).to_string(index=False))
prof.to_csv(f"{OUT}/A13b_{sp}_profile.tsv", sep="\t", index=False)
fl.to_csv(f"{OUT}/A13b_{sp}_flanks.tsv", sep="\t", index=False)
