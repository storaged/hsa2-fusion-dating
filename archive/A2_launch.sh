#!/bin/bash
cd ${HSA2_ROOT:-$HOME/projects/hsa2_followup}
mkdir -p results/A2
PY=${HSA2_PY:-python}
for b in human chimp bonobo; do for c in $(seq 1 22); do echo "chr$c $b"; done; done | \
  xargs -P 10 -L 1 sh -c "nice -n 15 $PY scripts/A2_ubcs_run.py results/A1/chunks \$0 \$1 2 0 results/A2/\$0.\$1.B2.tsv" > logs/A2_run.log 2>&1
echo A2_DONE >> logs/A2_run.log
