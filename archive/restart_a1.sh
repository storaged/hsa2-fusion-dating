#!/bin/bash
cd ${HSA2_ROOT:-$HOME/projects/hsa2_followup}
pkill -u krzysiek -f '[x]args -P 12 -L 1'
pkill -u krzysiek -f '[A]1_run_chunk.sh'
pkill -u krzysiek -f '[h]al2maf /home/krzysiek/projects/hsa2_followup'
pkill -u krzysiek -f '[A]1_branch_substitutions'
sleep 3
echo "remaining:"; ps -u krzysiek -o pid,cmd | grep -E '[h]al2maf|[A]1_' | head -3
sed -i 's#\[ -s \$out.callable.tsv \] && exit 0#[ -s $out.ils.tsv ] \&\& exit 0#' scripts/A1_run_chunk.sh
grep -n 'exit 0' scripts/A1_run_chunk.sh
rm -f results/A1/chunks/*; rm -rf /tmp/a1.*
nohup bash -c "cat results/A1/chunks.txt | xargs -P 12 -L 1 nice -n 15 scripts/A1_run_chunk.sh" > logs/A1_run.log 2>&1 &
sleep 2; echo restarted; ps -u krzysiek -o pid,cmd | grep -c '[h]al2maf'
