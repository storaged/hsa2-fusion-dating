#!/bin/bash
# A13a: great-ape LD recombination maps (Stevison et al. 2016 MBE; GARMaps, zenodo 13975), hg18 coordinates
#       -> hg38 -> hs1, averaged (length-weighted rho/kb) onto the A3 100-kb hs1 windows.
# Output: results/A13/pan_rho.<species>.hs1_100kb.tsv  (chrom, win_start, rho_kb, cov_bp)
set -euo pipefail
P=${HSA2_ROOT:-$HOME/projects/hsa2_followup}
cd $P
PY=${HSA2_PY:-python}
LO=${LIFTOVER:-liftOver}
G=data/pan_recomb/great-ape-recombination
O=results/A13
mkdir -p $O/tmp
for sp in Bonobo Chimp Gorilla; do
  s=$(echo $sp | tr A-Z a-z)
  # hg18 columns are in Mb; keep filter==0, positive-length intervals
  for f in "$G/Final $sp Map"/*.gz; do zcat "$f" | awk -F'\t' 'NR>1 && $13==0 {s=int($10*1e6); e=int($11*1e6); if (e>s) printf "%s\t%d\t%d\t%s\n", $9, s, e, $6}'; done \
    | sort -k1,1 -k2,2n > $O/tmp/$s.hg18.bed
  $LO -bedPlus=3 $O/tmp/$s.hg18.bed data/chains/hg18ToHg38.over.chain.gz $O/tmp/$s.hg38.bed $O/tmp/$s.hg18.unmapped
  $LO -bedPlus=3 $O/tmp/$s.hg38.bed data/chains/hg38ToHs1.over.chain.gz $O/tmp/$s.hs1.bed $O/tmp/$s.hg38.unmapped
  $PY - "$O/tmp/$s.hs1.bed" "$O/pan_rho.$s.hs1_100kb.tsv" <<'EOF'
import sys, pandas as pd, numpy as np
b = pd.read_csv(sys.argv[1], sep="\t", header=None, names=["chrom", "s", "e", "rho"])
b = b[(b.e - b.s) < 50_000]            # drop intervals stretched by liftover
b["w"] = (b.s // 100_000) * 100_000
b["len"] = b.e - b.s
g = b.groupby(["chrom", "w"]).apply(lambda x: pd.Series(dict(rho_kb=np.average(x.rho, weights=x.len), cov_bp=x.len.sum())))
g.reset_index().rename(columns={"w": "win_start"}).to_csv(sys.argv[2], sep="\t", index=False)
print(sys.argv[2], len(g))
EOF
  echo "$s: hg18 $(wc -l < $O/tmp/$s.hg18.bed)  hg38 $(wc -l < $O/tmp/$s.hg38.bed)  hs1 $(wc -l < $O/tmp/$s.hs1.bed)"
done
echo A13a_DONE
