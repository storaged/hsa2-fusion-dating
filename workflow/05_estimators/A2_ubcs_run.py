"""A2: UBCS per 1-Mb hs1 region for one branch and one chromosome, from A1 chunk outputs.

Follows tytus getUBCSFastStats.py / SNDWindow.py:
  SNDs  = substitutions on the branch (optionally excluding parent-CpG sites)
  biased = W->S;  p-hat = biased / all SNDs per 1-Mb region
  window 300 bp, number_of_bins B (bin size 300/B); for B == 300 the window is compressed to
  representative windows (legacy `compress`); dense configurations use the exact enumeration fallback
  only when affordable, otherwise the legacy fallback to 20 bins is applied and logged.
Output columns: chrom begin end branch B cpg_excluded p expected_BCS actual_BCS UBCS n_snd
usage: python A2_ubcs_run.py <chunk_dir> <chrom> <branch> <B> <exclude_cpg 0/1> <out.tsv>
"""
import glob
import math
import sys
from collections import defaultdict

import numpy as np
import pandas as pd

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from ubcs import MIN_SIZE, WTS_MIN_FRAC, _need, compress, prob_bc_dp

chunk_dir, chrom, branch, B, excl_cpg, out = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5]), sys.argv[6]
WINDOW, REGION = 300, 1_000_000

files = sorted(glob.glob(f"{chunk_dir}/{chrom}_*.subs.tsv.gz"), key=lambda f: int(f.rsplit("_", 1)[1].split(".")[0]))
s = pd.concat([pd.read_csv(f, sep="\t", usecols=["pos", "branch", "cls", "cpg"]) for f in files])
s = s[s.branch.isin(branch.split("+"))]   # e.g. "panstem+chimp" = Anc4 -> chimp lineage
if excl_cpg:
    s = s[s.cpg == 0]
s = s.drop_duplicates("pos").sort_values("pos")
pos = s.pos.to_numpy(np.int64)
biased = (s.cls == "WS").to_numpy()
region = pos // REGION
p_reg = pd.Series(biased).groupby(region).mean().to_dict()
n_reg = pd.Series(biased).groupby(region).size().to_dict()


def window_counts(i, nb):
    bin_size = WINDOW // nb
    mid_bin = pos[i] // bin_size
    lo = max(0, (mid_bin - nb + 1) * bin_size)
    hi = (mid_bin + nb) * bin_size  # exclusive
    a, b = np.searchsorted(pos, lo, "left"), np.searchsorted(pos, hi, "left")
    bins = (pos[a:b] // bin_size - mid_bin + nb - 1).astype(int)
    wc = np.bincount(bins, minlength=2 * nb - 1)[: 2 * nb - 1]
    wb = np.bincount(bins, weights=biased[a:b], minlength=2 * nb - 1)[: 2 * nb - 1]
    return wc, wb, b - a


exp_ = defaultdict(float)
obs_ = defaultdict(int)
fallback = 0
memo = {}
for i in range(len(pos)):
    r = region[i]
    nb = B
    wc, wb, ntot = window_counts(i, nb)
    if ntot < MIN_SIZE:
        continue
    # observed BC
    for k in range(nb):
        size = wc[k:k + nb].sum()
        if size >= MIN_SIZE and wb[k:k + nb].sum() >= _need(size, WTS_MIN_FRAC):
            obs_[r] += 1
            break
    bins = compress(list(wc), nb) if nb == WINDOW else list(wc)
    if nb == WINDOW and bins and max(sum(bins[k:k + len(bins) // 2 + 1]) for k in range(len(bins) // 2 + 1)) >= 23:
        fallback += 1                                   # legacy fallback to 20 bins
        wc20, _, _ = window_counts(i, 20)
        bins = list(wc20)
    key = (tuple(bins), round(p_reg[r], 6))
    if key not in memo:
        memo[key] = prob_bc_dp(bins, p_reg[r])
    exp_[r] += memo[key]

with open(out, "w") as f:
    for r in sorted(n_reg):
        f.write(f"{chrom}\t{r * REGION}\t{(r + 1) * REGION - 1}\t{branch}\t{B}\t{excl_cpg}\t{p_reg[r]:.6f}"
                f"\t{exp_[r]:.4f}\t{obs_[r]}\t{obs_[r] - exp_[r]:.4f}\t{n_reg[r]}\n")
print(f"{chrom} {branch} B={B}: {len(pos)} SNDs, fallback={fallback}, memo={len(memo)}", file=sys.stderr)
