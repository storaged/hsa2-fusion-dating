#!/bin/bash
# usage: A1_run_chunk.sh chrom start length
set -euo pipefail
P=${HSA2_ROOT:-$HOME/projects/hsa2_followup}
export PATH=$P/data/tools/cactus-bin-v2.9.3/bin:$PATH
export LD_LIBRARY_PATH=$P/data/tools/cactus-bin-v2.9.3/lib:${LD_LIBRARY_PATH:-}
source ~/miniconda3/etc/profile.d/conda.sh; conda activate $P/envs/hsa2
c=$1; s=$2; l=$3; id=${c}_${s}
out=$P/results/A1/chunks/$id
[ -s $out.ils.tsv ] && exit 0
tmp=$(mktemp -d /tmp/a1.XXXXXX)
hal2maf $P/data/cactus8/8-t2t-apes-2023v2.hal $tmp/x.maf --refGenome hs1 --refSequence $c --start $s --length $l --noDupes   --targetGenomes Anc3,Anc4,Anc5,GCA_028858775.2,GCA_029289425.2,GCA_029281585.2,GCA_028885655.2,GCA_028885625.2
python $P/scripts/A1_branch_substitutions.py $tmp/x.maf $out
rm -rf $tmp
