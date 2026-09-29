#!/bin/bash
# Assembly comparison, cell "hg38 reference, same T2T-apes 8-way HAL":
# clustered UBCS (B=2) on the hg38-referenced A1 output for the human and full Pan lineages,
# then the fusion-flank comparison (hg38 fusion site chr2:113,515,527-113,624,768, midpoint 113,570,147).
cd ${HSA2_ROOT:-$HOME/projects/hsa2_followup}
mkdir -p results/A2_hg38
PY=${HSA2_PY:-python}
for b in human panstem+chimp panstem+bonobo; do for c in $(seq 1 22); do echo "chr$c $b"; done; done | \
  xargs -P 8 -L 1 sh -c "nice -n 15 $PY scripts/A2_ubcs_run.py results/A1_hg38/chunks \$0 \$1 2 0 results/A2_hg38/\$0.\$1.B2.tsv" > logs/A2_hg38.log 2>&1
$PY -W ignore scripts/A2c_flank_ubcs_combined.py results/A2_hg38 data/ref/hg38.chrom.sizes 113570147 \
  results/A2_hg38/A2c_flank_ubcs_combined_hg38.tsv > results/A2_hg38/A2c_stdout.txt 2>&1
echo A2_HG38_DONE >> logs/A2_hg38.log
