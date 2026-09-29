#!/bin/bash
set -euo pipefail
cd ${HSA2_ROOT:-$HOME/projects/hsa2_followup}/data
B=https://cgl.gi.ucsc.edu/data/cactus/t2t-apes/8-t2t-apes-2023v2
for f in 8-t2t-apes-2023v2.hal 8-t2t-apes-2023v2.hal.md5 8-t2t-apes-2023v2.nh 8-t2t-apes-2023v2.README.md 8-t2t-apes-2023v2.hs1.maf.gz 8-t2t-apes-2023v2.hs1.maf.gz.md5 8-t2t-apes-2023v2.hs1.maf.coverage.tsv; do
  wget -q -c -P cactus8 $B/$f
done
wget -q -c -P ref https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.fa.gz
wget -q -c -P ref https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.chrom.sizes
# Cactus static binaries (hal tools)
cd tools && wget -q -c https://github.com/ComparativeGenomicsToolkit/cactus/releases/download/v2.9.3/cactus-bin-v2.9.3.tar.gz && tar xzf cactus-bin-v2.9.3.tar.gz
echo DONE_CORE
