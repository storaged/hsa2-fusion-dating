#!/bin/bash
# Launch the three parallel runs: haplotype noise floor (16-way, 2 role mappings), A6 SNP polarisation, A7 archaic strata.
cd ${HSA2_ROOT:-$HOME/projects/hsa2_followup}
chmod +x scripts/*.sh
PY=${HSA2_PY:-python}
HAL16=data/cactus16/16-t2t-apes-2023v2.hal
for cfg in h16a h16b; do
  T=$($PY -c "import json;print(','.join(v for k,v in json.load(open('scripts/$cfg.json')).items() if k!='hs1'))")
  nohup bash -c "cat results/A1/chunks.txt | xargs -P 4 -L 1 sh -c 'nice -n 15 scripts/A1_run_chunk_ref.sh hs1 $HAL16 results/A1_$cfg/chunks $T \$0 \$1 \$2 scripts/$cfg.json'; echo A1_${cfg}_DONE" > logs/A1_$cfg.log 2>&1 &
done
nohup bash -c "cat results/A6/regions.txt | xargs -P 5 -L 1 nice -n 15 scripts/A6_run_region.sh; echo A6_EXTRACT_DONE" > logs/A6_run.log 2>&1 &
nohup bash scripts/A7_run.sh > logs/A7_run.log 2>&1 &
sleep 5
echo "hal2maf: $(pgrep -u $USER -fc '[h]al2maf')  bcftools: $(pgrep -u $USER -fc '[b]cftools')  A7a: $(pgrep -u $USER -fc '[A]7a_sites')"
