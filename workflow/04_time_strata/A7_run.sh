#!/bin/bash
# A7 driver: human-branch substitutions -> hg19 liftover -> archaic genotypes (4 genomes) + 1KGP-CHM13 frequencies.
set -euo pipefail
P=${HSA2_ROOT:-$HOME/projects/hsa2_followup}
export PATH=$P/envs/hsa2/bin:$PATH
cd $P
O=results/A7
mkdir -p $O/targets_hg19 $O/archaic $O/kgp
python scripts/A7a_sites.py results/A1/chunks $O
liftOver -bedPlus=6 $O/human_subs.hs1.bed data/chains/hs1ToHg19.over.chain.gz $O/human_subs.hg19.bed $O/human_subs.hg19.unmapped
awk 'BEGIN{OFS="\t"}{c=$1; sub(/^chr/,"",c); print c, $3, $4, $6 > "'$O'/targets_hg19/" c ".tsv"}' $O/human_subs.hg19.bed
for f in $O/targets_hg19/*.tsv; do sort -k2,2n $f -o $f; done
# archaic genotypes: 4 genomes x chromosomes, 8 parallel
jobs=""
for c in $(seq 1 22) X; do
  [ -s $O/targets_hg19/$c.tsv ] || continue
  cut -f1,2 $O/targets_hg19/$c.tsv > $O/targets_hg19/$c.pos
  for g in Altai Vindija33.19 Denisova; do echo "$g $c data/archaic/$g/chr${c}_mq25_mapab100.vcf.gz"; done
  echo "Chagyrskaya $c data/archaic/Chagyrskaya/chr${c}.noRB.vcf.gz"
done | xargs -P 6 -L 1 sh -c 'nice -n 15 bcftools query -T '$O'/targets_hg19/$1.pos -f "%CHROM\t%POS\t%REF\t%ALT[\t%GT\t%GQ]\n" $2 > '$O'/archaic/$0.$1.tsv'
# modern 1KGP frequencies at hs1 positions (all + AFR), per chromosome
BCF=data/kgp/1KGP.CHM13v2.0.whole_genome.recalibrated.snp_indel.pass.phased.native_maps.biallelic.2504.bcf.gz
ls $O/targets_hs1/*.tsv | xargs -n1 basename | sed 's/.tsv//' | xargs -P 6 -I{} sh -c \
  "nice -n 15 bcftools view -r {} -T $O/targets_hs1/{}.tsv -v snps -Ou $BCF | bcftools +fill-tags -Ou -- -t AN,AC | bcftools +fill-tags -Ou -- -S data/kgp/groups_afr.txt -t AN,AC | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AN\t%AC\t%AN_AFR\t%AC_AFR\n' > $O/kgp/{}.tsv"
echo A7_EXTRACT_DONE
