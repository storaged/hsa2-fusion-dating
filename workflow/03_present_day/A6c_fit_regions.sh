#!/bin/bash
# A6c: pool polarised SNPs per region class (non-CpG, parsimony-supported) and fit B = 4Ne*b for ALL and AFR.
set -euo pipefail
P=${HSA2_ROOT:-$HOME/projects/hsa2_followup}
cd $P
PY=${HSA2_PY:-python}
$PY - <<'EOF'
import glob, pandas as pd
d = pd.concat([pd.read_csv(f, sep="\t") for f in glob.glob("results/A6/regions/*.tsv")])
d = d[(d.region_class != "test") & (d.cpg == 0) & (d.pars == 1)]
both = d[d.region_class.str.startswith("fusion")].assign(region_class="fusion_both")
d = pd.concat([d, both])
for pop in ("all", "afr"):
    x = d.rename(columns={f"dac_{pop}": "dac", f"n_{pop}": "n"})[["region_class", "cls", "dac", "n"]]
    x.to_csv(f"results/A6/fit_input_{pop}.tsv.gz", sep="\t", index=False)
print(d.groupby(["region_class", "cls"]).size().unstack())
EOF
for pop in all afr; do
  nice -n 15 $PY scripts/A6_gbgc_sfs.py results/A6/fit_input_$pop.tsv.gz results/A6/B_by_region_$pop.tsv 40 > results/A6/fit_$pop.log 2>&1 &
done
wait
echo A6_FIT_DONE
