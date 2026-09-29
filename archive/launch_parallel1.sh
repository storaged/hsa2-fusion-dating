#!/bin/bash
cd ${HSA2_ROOT:-$HOME/projects/hsa2_followup}
# 1. A1 with hg38 as reference (same 8-way HAL)
awk -v C=5000000 '$1 ~ /^chr([0-9]+|X)$/ {for(s=0;s<$2;s+=C){l=(s+C>$2)?$2-s:C; print $1, s, l}}' data/ref/hg38.chrom.sizes > results/A1_hg38_chunks.txt
T=Anc3,Anc4,Anc5,GCA_028858775.2,GCA_029289425.2,GCA_029281585.2,GCA_028885655.2,GCA_028885625.2
nohup bash -c "cat results/A1_hg38_chunks.txt | xargs -P 8 -L 1 nice -n 15 scripts/A1_run_chunk_ref.sh hg38 data/cactus8/8-t2t-apes-2023v2.hal results/A1_hg38/chunks $T; echo A1_HG38_DONE" > logs/A1_hg38.log 2>&1 &
# 2. 16-way diploid HAL for haplotype comparison
mkdir -p data/cactus16
B=https://cgl.gi.ucsc.edu/data/cactus/t2t-apes/16-t2t-apes-2023v2
nohup nice -n 10 bash -c "wget -q -c -P data/cactus16 $B/16-t2t-apes-2023v2.hal $B/16-t2t-apes-2023v2.hal.md5 $B/16-t2t-apes-2023v2.README.md; echo HAL16_DONE" > logs/D_cactus16.log 2>&1 &
# 3. A8 ILS profiles
envs/hsa2/bin/python -W ignore scripts/A8_ils_profiles.py results/A1/chunks data/ref/hs1.chrom.sizes results/A8 > results/A8/A8_stdout.txt 2>&1 || mkdir -p results/A8
[ -s results/A8/A8_stdout.txt ] || envs/hsa2/bin/python -W ignore scripts/A8_ils_profiles.py results/A1/chunks data/ref/hs1.chrom.sizes results/A8 > results/A8/A8_stdout.txt 2>&1
cat results/A8/A8_stdout.txt
wc -l results/A1_hg38_chunks.txt
