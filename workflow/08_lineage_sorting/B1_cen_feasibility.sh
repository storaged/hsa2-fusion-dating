#!/bin/bash
# B1 feasibility: is the 2q21 inactive alpha-satellite remnant (CHM13) usable as a divergence clock?
# Locate alpha satellite (RepeatMasker ALR/Alpha) on chr2 q; extract the remnant and a same-size piece of the
# active cen2 HOR array; tile into 2-kb windows; all-vs-all minimap2 (asm20, --eqx) within each; report the
# identity distribution of off-diagonal alignments (remnant vs active). Remnant copies should be more diverged
# (homogenisation stopped at inactivation) if the clock is readable.
set -euo pipefail
cd ${HSA2_ROOT:-$HOME/projects/hsa2_followup}
O=results/B1; mkdir -p $O
RM=data/censat/chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.bed
awk -F'\t' '$1=="chr2" && $2>100000000 && $0 ~ /ALR|Alpha/' $RM | sort -k2,2n > $O/alpha_2q_all.bed
# the remnant = the cluster within 500 kb of the largest 2q alpha hit (isolated single monomers elsewhere excluded)
c=$(awk -F'\t' '{print $3-$2, int(($2+$3)/2)}' $O/alpha_2q_all.bed | sort -k1,1nr | head -1 | cut -d' ' -f2)
awk -F'\t' -v c=$c '($2 > c-500000) && ($3 < c+500000)' $O/alpha_2q_all.bed > $O/remnant_alpha.bed
awk -F'\t' '{s+=$3-$2} END {print "remnant alpha bp:", s}' $O/remnant_alpha.bed
rs=$(head -1 $O/remnant_alpha.bed | cut -f2); re=$(tail -1 $O/remnant_alpha.bed | cut -f3)
echo "remnant span chr2:$rs-$re ($(( (re-rs)/1000 )) kb)"
len=$((re-rs))
# active array piece: middle of the largest hor_2 block in CenSat
read as ae < <(awk -F'\t' '$1=="chr2" && $4 ~ /^hor_2_/ {print $2, $3, $3-$2}' data/censat/chm13v2.0_censat_v2.1.bed | sort -k3,3nr | head -1 | awk -v L=$len '{m=int(($1+$2)/2); print m-int(L/2), m+int(L/2)}')
echo "active piece chr2:$as-$ae"
for nm in remnant active; do
  if [ $nm = remnant ]; then s=$rs; e=$re; else s=$as; e=$ae; fi
  samtools faidx data/ref/hs1.fa chr2:$s-$e > $O/$nm.fa
  ${HSA2_PY:-python} - $O/$nm.fa $O/$nm.tiles.fa <<'PY'
import sys
seq = "".join(l.strip() for l in open(sys.argv[1]) if not l.startswith(">"))
with open(sys.argv[2], "w") as f:
    for i in range(0, len(seq) - 2000 + 1, 1000):
        f.write(f">t{i}\n{seq[i:i + 2000]}\n")
PY
  minimap2 -x asm20 -X -c --eqx -t 4 $O/$nm.tiles.fa $O/$nm.tiles.fa 2>/dev/null > $O/$nm.ava.paf
  ${HSA2_PY:-python} - $O/$nm.ava.paf $nm <<'PY'
import sys, numpy as np
ids = []
for l in open(sys.argv[1]):
    x = l.split("\t")
    q0, t0 = int(x[0][1:]), int(x[5][1:])
    if abs(q0 - t0) < 2000 or int(x[10]) < 1000:
        continue                      # skip self/overlapping tiles and short alignments
    ids.append(int(x[9]) / int(x[10]))
ids = np.array(ids)
if len(ids):
    print(f"{sys.argv[2]}: {len(ids)} off-diagonal alignments; identity median {np.median(ids):.4f}, "
          f"IQR {np.percentile(ids, 25):.4f}-{np.percentile(ids, 75):.4f}, divergence median {1 - np.median(ids):.4f}")
else:
    print(f"{sys.argv[2]}: no off-diagonal alignments")
PY
done
echo B1_DONE
