"""A5: Is subtelomere-type recombination still active at the HSA2 fusion flanks today?

deCODE 1-Mb maps (Palsson et al. 2025; CO from Halldorsson et al. 2019), GRCh38:
  NCO map, CO cM/Mb, DSB map, oNCO map; paternal and maternal.
Window classes:
  end_k      : k-th Mb from a real telomere (non-acrocentric arms; skip unassembled edges)
  fusion_k   : k-th Mb from the fusion site (hg38 chr2:113,515,527-113,624,768), 2a (proximal) or 2b (distal)
  interior   : >= 20 Mb from any telomere, centromere and the fusion site
Output: per-class summaries and a per-distance profile; empirical percentile of the fusion flank
among real chromosome ends at the same distance.
"""
import sys
import numpy as np
import pandas as pd

DATA = sys.argv[1] if len(sys.argv) > 1 else "data"
OUT = sys.argv[2] if len(sys.argv) > 2 else "results/A5"
import os
os.makedirs(OUT, exist_ok=True)

FUSION = (113_515_527 + 113_624_768) / 2
# approximate hg38 centromere midpoints (UCSC cytoBand acen), Mb
CEN = {1: 123.4, 2: 93.9, 3: 90.9, 4: 50.0, 5: 48.8, 6: 59.8, 7: 60.1, 8: 45.2, 9: 43.0, 10: 39.8,
       11: 53.4, 12: 35.5, 13: 17.7, 14: 17.2, 15: 19.0, 16: 36.8, 17: 25.1, 18: 18.5, 19: 26.2,
       20: 28.1, 21: 12.0, 22: 15.0}
ACRO = {13, 14, 15, 21, 22}

sizes = pd.read_csv(f"{DATA}/ref/hg38.chrom.sizes", sep="\t", header=None, names=["chr", "len"]).set_index("chr")["len"]


def load(sex):
    m = pd.read_csv(f"{DATA}/decode/maps.{sex}.tsv", sep="\t", comment="#")
    m = m.rename(columns={"Chr": "chr"})
    m["sex"] = sex
    return m


m = pd.concat([load("pat"), load("mat")])
m = m[m.chr.isin([f"chr{i}" for i in range(1, 23)])].copy()
m["cn"] = m.chr.str[3:].astype(int)
m["len"] = m.chr.map(sizes)
m["d_p"] = m.pos / 1e6            # Mb from p-telomere (window centre)
m["d_q"] = (m.len - m.pos) / 1e6  # Mb from q-telomere
m["d_fus"] = np.where(m.cn == 2, (m.pos - FUSION) / 1e6, np.nan)
m["d_cen"] = np.abs(m.pos / 1e6 - m.cn.map(CEN))

rows = []
for _, r in m.iterrows():
    cls, k = "other", np.nan
    if r.cn == 2 and abs(r.d_fus) < 20:
        cls, k = ("fusion_2b" if r.d_fus > 0 else "fusion_2a"), int(abs(r.d_fus))
    elif r.d_p < 20 and r.cn not in ACRO:
        cls, k = "end_p", int(r.d_p)
    elif r.d_q < 20:
        cls, k = "end_q", int(r.d_q)
    elif r.d_p >= 20 and r.d_q >= 20 and r.d_cen >= 10:
        cls, k = "interior", np.nan
    rows.append((cls, k))
m["cls"], m["k"] = zip(*rows)
m["end_id"] = np.where(m.cls.str.startswith("end"), m.chr + m.cls.str[-1], None)

vars_ = ["cMperMb", "map", "DSB", "oNCO"]
m.to_csv(f"{OUT}/A5_windows_classified.tsv", sep="\t", index=False)

# 1. profile by distance (0..14 Mb) for ends vs fusion flanks
prof = []
for sex in ["pat", "mat"]:
    for k in range(15):
        ends = m[(m.sex == sex) & m.cls.str.startswith("end") & (m.k == k)]
        inter = m[(m.sex == sex) & (m.cls == "interior")]
        for side in ["fusion_2a", "fusion_2b"]:
            f = m[(m.sex == sex) & (m.cls == side) & (m.k == k)]
            if f.empty:
                continue
            for v in vars_:
                ev = ends[v].dropna()
                fv = f[v].values[0]
                prof.append(dict(sex=sex, k_Mb=k, side=side, var=v, fusion=fv, ends_median=ev.median(),
                                 ends_q25=ev.quantile(.25), ends_q75=ev.quantile(.75),
                                 interior_median=inter[v].median(),
                                 pct_among_ends=(ev < fv).mean()))
prof = pd.DataFrame(prof)
prof.to_csv(f"{OUT}/A5_profile_by_distance.tsv", sep="\t", index=False)

# 2. summary over first 5 and 10 Mb
summ = []
for sex in ["pat", "mat"]:
    for K in [2, 5, 10]:
        sub = m[(m.sex == sex) & (m.k < K)]
        inter = m[(m.sex == sex) & (m.cls == "interior")]
        per_end = sub[sub.cls.str.startswith("end")].groupby("end_id")[vars_].mean()
        for side in ["fusion_2a", "fusion_2b"]:
            f = sub[sub.cls == side][vars_].mean()
            for v in vars_:
                summ.append(dict(sex=sex, first_Mb=K, side=side, var=v, fusion=f[v],
                                 ends_mean=per_end[v].mean(), ends_median=per_end[v].median(),
                                 interior_mean=inter[v].mean(),
                                 pct_among_ends=(per_end[v] < f[v]).mean(),
                                 rel_position=(f[v] - inter[v].mean()) / (per_end[v].mean() - inter[v].mean())))
summ = pd.DataFrame(summ)
summ.to_csv(f"{OUT}/A5_summary.tsv", sep="\t", index=False)
pd.set_option("display.width", 200)
print(summ.round(3).to_string(index=False))
