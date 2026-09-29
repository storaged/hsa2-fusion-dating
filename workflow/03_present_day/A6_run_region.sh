#!/bin/bash
# usage: A6_run_region.sh chrom start length class
# 1KGP-CHM13 biallelic SNPs (AN/AC for all 2504 and for AFR) in the region, polarised by Anc4 (8-way HAL).
set -euo pipefail
P=${HSA2_ROOT:-$HOME/projects/hsa2_followup}
export PATH=$P/envs/hsa2/bin:$P/data/tools/cactus-bin-v2.9.3/bin:$PATH   # conda bcftools (with plugins) first
export LD_LIBRARY_PATH=$P/data/tools/cactus-bin-v2.9.3/lib:${LD_LIBRARY_PATH:-}
c=$1; s=$2; l=$3; cls=$4
out=$P/results/A6/regions/${c}_${s}.tsv
[ -s $out ] && exit 0
mkdir -p $P/results/A6/regions
tmp=$(mktemp -d /tmp/a6.XXXXXX)
BCF=$P/data/kgp/1KGP.CHM13v2.0.whole_genome.recalibrated.snp_indel.pass.phased.native_maps.biallelic.2504.bcf.gz
bcftools view -r $c:$((s+1))-$((s+l)) -v snps -m2 -M2 -Ou $BCF | \
  bcftools +fill-tags -Ou -- -t AN,AC | \
  bcftools +fill-tags -Ou -- -S $P/data/kgp/groups_afr.txt -t AN,AC | \
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AN\t%AC\t%AN_AFR\t%AC_AFR\n' > $tmp/snps.tsv
hal2maf $P/data/cactus8/8-t2t-apes-2023v2.hal $tmp/x.maf --refGenome hs1 --refSequence $c --start $s --length $l --noDupes \
  --targetGenomes Anc4,GCA_028858775.2,GCA_029289425.2,GCA_029281585.2
python $P/scripts/A6a_polarize.py $tmp/x.maf $tmp/snps.tsv $tmp/pol.tsv
awk -v C=$cls 'BEGIN{OFS="\t"} NR==1{print $0,"region_class"; next}{print $0,C}' $tmp/pol.tsv > $out
rm -rf $tmp
