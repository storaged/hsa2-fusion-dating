#!/bin/bash
# A14: haplotype noise floor. For the two 16-way role mappings (h16a: CHM13 + primary/mat apes;
# h16b: HG002-mat + other ape haplotypes): clustered UBCS cell (as A2_hg38.sh) and the primary
# W/S log-odds statistic (A14_logodds_cells.py) for 8-way CHM13, h16a, h16b with one method.
cd ${HSA2_ROOT:-$HOME/projects/hsa2_followup}
PY=${HSA2_PY:-python}
for cell in h16a h16b; do
  mkdir -p results/A2_$cell
  for b in human panstem+chimp panstem+bonobo; do for c in $(seq 1 22); do echo "chr$c $b"; done; done | \
    xargs -P 8 -L 1 sh -c "nice -n 15 $PY scripts/A2_ubcs_run.py results/A1_$cell/chunks \$0 \$1 2 0 results/A2_$cell/\$0.\$1.B2.tsv" > logs/A2_$cell.log 2>&1
  $PY -W ignore scripts/A2c_flank_ubcs_combined.py results/A2_$cell data/ref/hs1.chrom.sizes 113994777 \
    results/A2_$cell/A2c_flank_ubcs_combined_$cell.tsv > results/A2_$cell/A2c_stdout.txt 2>&1
  echo A2_${cell}_DONE >> logs/A2_$cell.log
done
echo A14_UBCS_DONE
