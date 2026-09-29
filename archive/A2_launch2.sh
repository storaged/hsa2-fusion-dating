#!/bin/bash
cd ${HSA2_ROOT:-$HOME/projects/hsa2_followup}
PY=${HSA2_PY:-python}
for b in panstem+chimp panstem+bonobo panstem gorilla; do for c in $(seq 1 22); do echo "chr$c $b"; done; done | \
  xargs -P 6 -L 1 sh -c "nice -n 15 $PY scripts/A2_ubcs_run.py results/A1/chunks \$0 \$1 2 0 results/A2/\$0.\$1.B2.tsv" > logs/A2_run2.log 2>&1
echo A2_2_DONE >> logs/A2_run2.log
