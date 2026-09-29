"""UBCS v2 engine: probability that a substitution is Biased-Clustered (BC) under the null of no
association between W->S bias and clustering (Dreszer et al. 2007; Poszewiecka et al. 2022).

`prob_bc_dp` is a faithful port of tytus `SNDWindow.get_prob_of_bcs` (the 2022 dynamic programme,
Algorithm 3). `prob_bc_exact` enumerates bias counts per bin (exact, exponential; used for testing and
for dense clusters where the legacy code fell back to coarser bins).

Bins: the 2n-1 consecutive bin counts around the focal substitution; cluster k (k = 0..n-1) is
bins[k : k + n] and always contains the middle bin (index n - 1), which holds the focal substitution.
A cluster is BC if it has >= MIN_SIZE substitutions and >= WTS_MIN_FRAC of them are biased.
"""
import math
from functools import lru_cache
from itertools import product

MIN_SIZE = 5
WTS_MIN_FRAC = 0.8


def _need(size, frac):
    # smallest number of biased substitutions making a cluster of `size` biased (legacy: n_b/size >= frac)
    return math.ceil(frac * size - 1e-12)


@lru_cache(maxsize=None)
def _binom(n, k, p):
    if k < 0 or k > n:
        return 0.0
    return math.comb(n, k) * p ** k * (1 - p) ** (n - k)


def _binom_from(n, k0, p):
    return sum(_binom(n, k, p) for k in range(max(k0, 0), n + 1))


def prob_bc_dp(bins, p, min_size=MIN_SIZE, frac=WTS_MIN_FRAC):
    """Port of the 2022 DP (SNDWindow.get_prob_of_bcs)."""
    bins = tuple(bins)
    if not bins:
        return 0.0
    n = len(bins) // 2 + 1
    prob = 0.0
    first = sum(bins[:n])
    if first >= min_size:
        prob = _binom_from(first, _need(first, frac), p)
    mem = None  # None == default 1
    for k in range(1, n):
        cond_counts = bins[k:k + n - 1]
        cluster_size = sum(bins[k:k + n])
        prev_size = sum(bins[k - 1:k + n - 1])
        prev_mem, mem = mem, {}
        for cond in product(*(range(c + 1) for c in cond_counts)):
            cs = sum(cond)
            cprob = 1.0
            for c, x in zip(cond_counts, cond):
                cprob *= _binom(c, x, p)
            a = _binom_from(bins[k + n - 1], max(_need(cluster_size, frac) - cs, 0), p)
            if prev_size < min_size:
                ub = bins[k - 1] + 1
            else:
                ub = min(_need(prev_size, frac) - cs, bins[k - 1] + 1)
            na = 0.0
            for x in range(0, max(ub, 0)):
                key = (x,) + cond[:-1]
                na += _binom(bins[k - 1], x, p) * (1.0 if prev_mem is None else prev_mem.get(key, 0.0))
            if cluster_size >= min_size:
                prob += a * na * cprob
            mem[cond] = na
    return prob


def prob_bc_exact(bins, p, min_size=MIN_SIZE, frac=WTS_MIN_FRAC):
    """Exact by enumeration of biased counts per bin."""
    bins = tuple(bins)
    if not bins:
        return 0.0
    n = len(bins) // 2 + 1
    sizes = [sum(bins[k:k + n]) for k in range(n)]
    need = [_need(s, frac) for s in sizes]
    tot = 0.0
    for xs in product(*(range(b + 1) for b in bins)):
        w = 1.0
        for b, x in zip(bins, xs):
            w *= _binom(b, x, p)
        for k in range(n):
            if sizes[k] >= min_size and sum(xs[k:k + n]) >= need[k]:
                tot += w
                break
    return tot


def compress(window_counts, n_bins, min_size=MIN_SIZE):
    """Legacy compression of a (2*n_bins-1)-vector of bin counts to representative windows."""
    starts, seen = [], set()
    for i in range(n_bins):
        seg = window_counts[i:i + n_bins]
        if sum(seg) >= min_size:
            ind = tuple(i + j for j, v in enumerate(seg) if v > 0)
            if ind not in seen:
                seen.add(ind)
                starts.append(i)
    starts = starts + [s + n_bins for s in starts]
    return [sum(window_counts[s:e]) for s, e in zip(starts, starts[1:])]


if __name__ == "__main__":
    import random
    random.seed(1)
    worst = 0.0
    for trial in range(400):
        n = random.randint(1, 4)
        bins = [random.randint(0, 4) for _ in range(2 * n - 1)]
        bins[n - 1] = max(bins[n - 1], 1)
        p = random.choice([0.2, 0.35, 0.5, 0.7])
        a, b = prob_bc_dp(bins, p), prob_bc_exact(bins, p)
        worst = max(worst, abs(a - b))
        if abs(a - b) > 1e-9:
            print("MISMATCH", bins, p, a, b)
    print("max |DP - exact| over 400 random configurations:", worst)
