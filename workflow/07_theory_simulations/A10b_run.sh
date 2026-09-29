#!/bin/bash
# A10b driver: SLiM switch-off grid (bias x T_off x replicates), then the estimator on pooled replicates.
set -euo pipefail
cd ${HSA2_ROOT:-$HOME/projects/hsa2_followup}
O=results/A10b; mkdir -p $O
PROCS=${1:-32}; REPS=${2:-40}
seed=1000
for bias in 0.25 0.5 1.0; do for toff in 0 1000 2000 3000 4000; do for r in $(seq 1 $REPS); do
  seed=$((seed+1)); echo "$bias $toff $r $seed"; done; done; done | \
  xargs -P $PROCS -L 1 sh -c '${SLIM:-slim} -s $3 -d N=1000 -d SEG=50000 -d MU=2.5e-6 -d RHI=1.85e-6 -d RLO=1e-8 -d NCOF=0.9 -d TRACT=300 -d BIAS=$0 -d BURN=10000 -d TS=4000 -d TOFF=$1 -d REP=$2 scripts/A10b_slim_switchoff.slim 2>/dev/null | grep ^RES' > $O/raw.tsv
wc -l $O/raw.tsv
${HSA2_PY:-python} - <<'PY'
import numpy as np, pandas as pd
r = pd.read_csv("results/A10b/raw.tsv", sep="\t", header=None, names=["tag","rep","bias","toff","seg","ws","sw","nW","nS"])
TS = 4000
rng = np.random.default_rng(10)
def lo(x):
    return np.log(x.ws.sum() / x.sw.sum()) - np.log(x.nW.sum() / x.nS.sum())
def est(g):
    e, f, i = (g[g.seg == k] for k in (0, 1, 2))
    rho = (lo(f) - lo(i)) / (lo(e) - lo(i))
    return TS * (1 - rho), lo(e) - lo(i)
rows = []
for (b, t), g in r.groupby(["bias", "toff"]):
    T, ex = est(g)
    reps = g.rep.unique()
    bs = [est(pd.concat([g[g.rep == x] for x in rng.choice(reps, len(reps))]))[0] for _ in range(500)]
    rows.append(dict(bias=b, T_true=t, n_reps=len(reps), end_excess=ex, T_hat=T, T_lo=np.percentile(bs, 2.5),
                     T_hi=np.percentile(bs, 97.5), covered=np.percentile(bs, 2.5) <= t <= np.percentile(bs, 97.5)))
R = pd.DataFrame(rows)
R["bias_gen"] = R.T_hat - R.T_true
R.to_csv("results/A10b/A10b_estimates.tsv", sep="\t", index=False)
pd.set_option("display.width", 200)
print(R.round(3).to_string(index=False))
print(f"coverage of true T_off by 95% bootstrap CI: {R.covered.mean():.2f}; mean |bias| {R.bias_gen.abs().mean():.0f} gen "
      f"({R.bias_gen.abs().mean() / TS:.3f} of T_s)")
PY
echo A10b_DONE
