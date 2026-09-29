#!/bin/bash
# A12c: present-day gBGC strength B (SFS fit, A6 fitter) per bin of present-day paternal DSB rate
# (A12 calibration's best end predictor), for real-end and interior regions; AFR and ALL samples.
# SNPs from A6 regions (hs1), joined to deCODE via the A12 hs1 100-kb window table.
set -euo pipefail
P=${HSA2_ROOT:-$HOME/projects/hsa2_followup}; cd $P
PY=${HSA2_PY:-python}
NPROJ=${1:-40}
mkdir -p results/A12c
$PY - <<'PYEOF'
import glob, numpy as np, pandas as pd
w = pd.read_csv("results/A12/A12_windows.tsv.gz", sep="\t", usecols=["chrom", "win_start", "dsb_pat"]).dropna()
d = pd.concat([pd.read_csv(f, sep="\t") for f in glob.glob("results/A6/regions/*.tsv")])
d = d[(d.region_class.isin(["end", "interior"])) & (d.cpg == 0) & (d.pars == 1)]
d["win_start"] = (d.pos // 100_000) * 100_000
d = d.merge(w, on=["chrom", "win_start"], how="inner")
q = np.quantile(w.dsb_pat, [0, .5, .75, .9, .97, 1])
d["rbin"] = pd.cut(d.dsb_pat, q, include_lowest=True, labels=["q0-50", "q50-75", "q75-90", "q90-97", "q97-100"]).astype(str)
d["region_class"] = d.rbin
for pop in ("all", "afr"):
    x = d.rename(columns={f"dac_{pop}": "dac", f"n_{pop}": "n"})[["region_class", "cls", "dac", "n"]]
    x.to_csv(f"results/A12c/fit_input_{pop}.tsv.gz", sep="\t", index=False)
print("dsb_pat quantile edges:", np.round(q, 3))
print(d.groupby(["rbin", "cls"]).size().unstack())
PYEOF
for pop in all afr; do
  $PY scripts/A6_gbgc_sfs.py results/A12c/fit_input_$pop.tsv.gz results/A12c/B_by_dsbbin_$pop.tsv $NPROJ > results/A12c/fit_$pop.log 2>&1 &
done
wait
cat results/A12c/B_by_dsbbin_afr.tsv results/A12c/B_by_dsbbin_all.tsv
echo A12c_DONE
