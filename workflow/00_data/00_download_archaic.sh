#!/bin/bash
# D9: archaic hominin genomes (hg19, uniformly processed, Pruefer et al. 2017; Mafessoni et al. 2020)
cd ${HSA2_ROOT:-$HOME/projects/hsa2_followup}/data/archaic
B=https://cdna.eva.mpg.de/neandertal
for c in $(seq 1 22) X; do
  for g in Altai Vindija33.19 Denisova; do
    mkdir -p $g; wget -q -c -P $g $B/Vindija/VCF/$g/chr${c}_mq25_mapab100.vcf.gz $B/Vindija/VCF/$g/chr${c}_mq25_mapab100.vcf.gz.tbi
  done
  mkdir -p Chagyrskaya; wget -q -c -P Chagyrskaya $B/Chagyrskaya/VCF/chr${c}.noRB.vcf.gz $B/Chagyrskaya/VCF/chr${c}.noRB.vcf.gz.tbi
done
echo ARCHAIC_DONE
