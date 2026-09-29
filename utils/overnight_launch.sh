#!/bin/bash
# Overnight batch 2026-09-28: all jobs niced, BLAS threads capped; bookkeeping via run_job.sh; monitor every 10 min.
cd ${HSA2_ROOT:-$HOME/projects/hsa2_followup}
export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 MKL_NUM_THREADS=4
PY=${HSA2_PY:-python}
R=scripts/run_job.sh
rm -rf results/A10_test results/A15_test results/A14/test_chunks results/A14/test_windows.tsv.gz
mkdir -p results/A10 results/A15
nohup bash $R A15_pseudo_fusion   $PY -W ignore scripts/A15_pseudo_fusion.py results/A3/windows.tsv.gz results/A15 2000 6.0 </dev/null &
nohup bash $R A10a_N200_Ts8        $PY -W ignore scripts/A10a_nonlinear_transient.py results/A10 200 8 </dev/null &
nohup bash $R A10a_N200_Ts16       $PY -W ignore scripts/A10a_nonlinear_transient.py results/A10 200 16 </dev/null &
nohup bash $R A10a_N400_Ts8        $PY -W ignore scripts/A10a_nonlinear_transient.py results/A10 400 8 </dev/null &
nohup bash $R A10c_ils_structured  $PY -W ignore scripts/A10c_ils_structured.py results/A10 50000 16 </dev/null &
nohup bash $R A12c_B_by_recomb     bash scripts/A12c_B_by_recomb.sh 40 </dev/null &
nohup bash $R A7c_composition      $PY scripts/A7c_composition_correct.py results/A7/A7_summary.tsv results/A12/A12_windows.tsv.gz results/A7/A7c_corrected.tsv </dev/null &
nohup bash $R HAP_chain bash -c '
  set -e
  for c in h16a h16b; do
    mkdir -p results/A3_$c results/A15_$c results/A16_$c
    '"$PY"' -W ignore scripts/A3_aggregate_chunks.py results/A1_$c/chunks data/asm_reports results/A3_$c/windows.tsv.gz
    '"$PY"' -W ignore scripts/A15_pseudo_fusion.py results/A3_$c/windows.tsv.gz results/A15_$c 2000 6.0
    '"$PY"' -W ignore scripts/A12_recomb_calibration.py results/A3_$c/windows.tsv.gz data results/A12_$c > results/A12_$c.log 2>&1 || true
    '"$PY"' -W ignore scripts/A16_2b_investigation.py results/A3_$c/windows.tsv.gz results/A1_$c/chunks results/A12/A12_windows.tsv.gz results/A16_$c
    echo "$c done"
  done' </dev/null &
sleep 5
nohup bash scripts/overnight_monitor.sh > /dev/null 2>&1 </dev/null &
echo launched
