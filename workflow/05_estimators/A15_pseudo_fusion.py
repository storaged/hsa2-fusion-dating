"""A15: calibrate the T_off estimators on real chromosome ends ("pseudo-fusions", true T = 0).

Primary statistic: non-CpG W->S vs S->W log-odds, W/S-composition corrected, excess over the lineage's
interior (>= 30 Mb from ends, excluding chr2 within 30 Mb of the fusion). hs1 100-kb windows from A3.
Two references:
  within-human : frac_e = human excess(end e) / human excess(all other real ends, pooled)
  Pan-referenced: frac_e = human excess(e) / Pan-lineage excess at the same hs1 windows
                  (Pan lineage = panstem + chimp, or panstem + bonobo; counts summed)
T_e = T_s (1 - frac_e). Real ends have T = 0, so the distribution of T_e is the estimator's error
distribution including end heterogeneity. Pseudo-fusion pairs (two random real ends, pooled) mimic the
two-sided fusion estimate. The fusion flank estimate is then placed against these distributions:
  empirical CI for the flank: T_flank - quantiles(T_e)  (bias- and heterogeneity-aware).
usage: python A15_pseudo_fusion.py results/A3/windows.tsv.gz results/A15 [n_pairs=2000] [T_s=6.0]
"""
import os
import sys

import numpy as np
import pandas as pd

win_f, OUT = sys.argv[1:3]
NP = int(sys.argv[3]) if len(sys.argv) > 3 else 2000
TS = float(sys.argv[4]) if len(sys.argv) > 4 else 6.0
os.makedirs(OUT, exist_ok=True)
rng = np.random.default_rng(15)
F = (113_940_058 + 114_049_496) / 2
ACRO = {"chr13", "chr14", "chr15", "chr21", "chr22"}
cols = ["chrom", "win_start", "branch", "n_callable", "parent_W", "parent_S", "WS", "SW", "sp_dist_end", "sp_arm_end"]
d = pd.read_csv(win_f, sep="\t", usecols=cols)
d = d[d.chrom.str.match(r"^chr\d+$") & (d.n_callable > 20_000)]
# lineages on hs1 windows: human; Pan = panstem + chimp / panstem + bonobo
piv = d.pivot_table(index=["chrom", "win_start"], columns="branch", values=["WS", "SW", "parent_W", "parent_S"], aggfunc="sum")
lin = {}
lin["human"] = piv.xs("human", axis=1, level=1)
for tip in ("chimp", "bonobo"):
    x = piv.xs("panstem", axis=1, level=1)[["WS", "SW"]] + piv.xs(tip, axis=1, level=1)[["WS", "SW"]]
    x[["parent_W", "parent_S"]] = piv.xs("panstem", axis=1, level=1)[["parent_W", "parent_S"]]
    lin[f"panstem+{tip}"] = x
hum = d[d.branch == "human"].set_index(["chrom", "win_start"])
L = pd.DataFrame(index=lin["human"].index)
L["mid"] = L.index.get_level_values(1) + 50_000
L["chrom"] = L.index.get_level_values(0)
L["dist"] = hum.sp_dist_end.reindex(L.index)
L["arm"] = hum.sp_arm_end.reindex(L.index)
L["end_id"] = L.chrom + L.arm.fillna("")
near2 = (L.chrom == "chr2") & ((L.mid - F).abs() < 30e6)
L["dfus"] = np.where(L.chrom == "chr2", (L.mid - F).abs(), np.inf)


def ex(x):
    return np.log(x.WS.sum() / x.SW.sum()) - np.log(x.parent_W.sum() / x.parent_S.sum())


rows, flank_rows = [], []
for D in (2, 5, 10):
    is_end = (L.dist < D * 1e6) & ~near2 & ~((L.arm == "p") & L.chrom.isin(ACRO))
    is_int = (L.dist >= 30e6) & ~near2
    fl = {"2a": (L.dfus < D * 1e6) & (L.mid < F), "2b": (L.dfus < D * 1e6) & (L.mid >= F)}
    fl["both"] = fl["2a"] | fl["2b"]
    ends = L[is_end].end_id.unique()
    for refname in ("within_human", "panstem+chimp", "panstem+bonobo"):
        H = lin["human"].dropna()
        R = H if refname == "within_human" else lin[refname].dropna()
        idx = H.index.intersection(R.index)
        H, R = H.loc[idx], R.loc[idx]
        m_end, m_int = is_end.reindex(idx).fillna(False), is_int.reindex(idx).fillna(False)
        eid = L.end_id.reindex(idx)
        hi, ri = ex(H[m_int]), ex(R[m_int])
        per_end = {}
        for e in ends:
            me = m_end & (eid == e)
            if (H[me].WS.sum() + H[me].SW.sum()) < 300:
                continue
            h_e = ex(H[me]) - hi
            if refname == "within_human":
                others = m_end & (eid != e)
                r_e = ex(H[others]) - hi
            else:
                r_e = ex(R[me]) - ri
            per_end[e] = (h_e, r_e, me)
        T_e = pd.Series({e: TS * (1 - h / r) for e, (h, r, _) in per_end.items()})
        # pseudo-fusion pairs: pool two random real ends (as the two flanks are pooled)
        keys = list(per_end)
        Tp = []
        for _ in range(NP):
            a, b = rng.choice(len(keys), 2, replace=False)
            mm = per_end[keys[a]][2] | per_end[keys[b]][2]
            h = ex(H[mm]) - hi
            if refname == "within_human":
                r = ex(H[m_end & ~mm]) - hi
            else:
                r = ex(R[mm]) - ri
            Tp.append(TS * (1 - h / r))
        Tp = np.array(Tp)
        # fusion flanks with the same reference
        for side, msk in fl.items():
            mf = msk.reindex(idx).fillna(False)
            h = ex(H[mf]) - hi
            r = (ex(H[m_end]) - hi) if refname == "within_human" else (ex(R[mf]) - ri)
            Tf = TS * (1 - h / r)
            null = Tp if side == "both" else T_e.values
            q = np.percentile(null, [2.5, 50, 97.5])
            flank_rows.append(dict(first_Mb=D, reference=refname, side=side, T_raw=Tf,
                                   null_median=q[1], null_sd=np.std(null),
                                   T_debiased=Tf - q[1], T_emp_lo=Tf - q[2], T_emp_hi=Tf - q[0],
                                   pct_of_flank_in_null=(null < Tf).mean()))
        rows.append(dict(first_Mb=D, reference=refname, n_ends=len(T_e), Te_median=T_e.median(), Te_mean=T_e.mean(),
                         Te_sd=T_e.std(), Te_q025=T_e.quantile(.025), Te_q975=T_e.quantile(.975),
                         pair_median=np.median(Tp), pair_sd=Tp.std(), pair_q025=np.percentile(Tp, 2.5),
                         pair_q975=np.percentile(Tp, 97.5)))
        T_e.rename("T_e").to_csv(f"{OUT}/A15_Te_{refname}_{D}Mb.tsv", sep="\t")
    print(f"D={D} done", flush=True)
R1, R2 = pd.DataFrame(rows), pd.DataFrame(flank_rows)
R1.to_csv(f"{OUT}/A15_null_distributions.tsv", sep="\t", index=False)
R2.to_csv(f"{OUT}/A15_flank_calibrated.tsv", sep="\t", index=False)
pd.set_option("display.width", 220)
print("== estimator on real ends (true T = 0): single ends and pseudo-fusion pairs")
print(R1.round(3).to_string(index=False))
print("== fusion flanks, calibrated against the real-end null (T_debiased = T_raw - null median; empirical 95% CI)")
print(R2.round(3).to_string(index=False))
