#!/bin/bash
# usage: A1_run_chunk_ref.sh REF HAL OUTDIR TARGETS chrom start length [ROLE_CONFIG.json]
# Generic A1 chunk runner for the assembly/haplotype comparison (reference genome, HAL and role mapping configurable).
set -euo pipefail
P=${HSA2_ROOT:-$HOME/projects/hsa2_followup}
export PATH=$P/data/tools/cactus-bin-v2.9.3/bin:$PATH
export LD_LIBRARY_PATH=$P/data/tools/cactus-bin-v2.9.3/lib:${LD_LIBRARY_PATH:-}
ref=$1; hal=$2; outdir=$3; targets=$4; c=$5; s=$6; l=$7; cfg=${8:-$ref}
out=$outdir/${c}_${s}
[ -s $out.ils.tsv ] && exit 0
mkdir -p $outdir
tmp=$(mktemp -d /tmp/a1r.XXXXXX)
hal2maf $hal $tmp/x.maf --refGenome $ref --refSequence $c --start $s --length $l --noDupes --targetGenomes $targets
${HSA2_PY:-python} $P/scripts/A1_branch_substitutions.py $tmp/x.maf $out $cfg
rm -rf $tmp
