#!/bin/bash
# Smoke tests (~1 min): numerical checks of the model identities and small versions of the theory/simulation steps.
# Each step must finish and print its completion marker.
set -euo pipefail
PY=${PYTHON:-python}
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
step() {  # step <n> <title> <marker> <command...>
  local n=$1 title=$2 marker=$3; shift 3
  echo "[$n/3] $title"
  "$@" > "$T/out.txt" 2>&1 || { echo "  FAILED"; tail -20 "$T/out.txt"; exit 1; }
  grep -q "$marker" "$T/out.txt" || { echo "  FAILED: no '$marker' in output"; tail -20 "$T/out.txt"; exit 1; }
  echo "  ok"
}
step 1 "exact Wright-Fisher identities (first-order, no-lag identity)" "ratio" \
  $PY -W ignore workflow/07_theory_simulations/M0_model_checks_exact.py
step 2 "estimator bias, exact Wright-Fisher, small population" "A10a_DONE" \
  $PY -W ignore workflow/07_theory_simulations/A10a_nonlinear_transient.py "$T" 40 4
awk '/bias \(T_hat/{f=1} f && /step +all/{print "  step switch-off, single-genome counts, bias (N units):", $0; exit}' "$T/out.txt"
step 3 "structured coalescent (msprime), few replicates" "A10c_DONE" \
  $PY -W ignore workflow/07_theory_simulations/A10c_ils_structured.py "$T" 200 2
echo "SMOKE TESTS PASSED"
