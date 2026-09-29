"""A3d: retained GC-bias excess on the human branch relative to Pan branches at the SAME orthologous
windows (hs1 chr2 fusion flanks, 0-D Mb). excess_b = g_flank_b / g_interior_b - 1.
Ratio rho = excess_human / excess_ref; T_off(+lag) ~ T_eff * (1 - rho). Block bootstrap over 1-Mb blocks
(paired: same blocks for both branches)."""
import numpy as np, pandas as pd
rng = np.random.default_rng(1)
d = pd.read_csv("results/A3/windows.tsv.gz", sep="\t", dtype={"chrom": str})
d = d[(d.n_callable > 20000) & d.chrom.str.match(r"^chr\d+$")]
F = (113_940_058 + 114_049_496) / 2
def g(x): return (x.WS.sum() / x.parent_W.sum()) / (x.SW.sum() / x.parent_S.sum())
out = []
for D in (5e6, 10e6, 15e6):
    fl = d[(d.chrom == "chr2") & ((d.win_start + 5e4 - F).abs() < D)].copy()
    fl["side"] = np.where(fl.win_start + 5e4 < F, "2a", "2b")
    fl["blk"] = fl.win_start // 1_000_000
    gi = {b: g(d[(d.branch == b) & ~((d.chrom == "chr2") & ((d.win_start - F).abs() < 30e6))]) for b in ["human", "chimp", "bonobo", "panstem"]}
    for side in ["2a", "2b", "both"]:
        f = fl if side == "both" else fl[fl.side == side]
        blocks = f.blk.unique()
        for ref in ["chimp", "bonobo", "pan_mean"]:
            def rho(bl):
                s = f[f.blk.isin(bl)] if bl is not None else f
                if bl is not None:
                    s = pd.concat([f[f.blk == k] for k in bl])
                eh = g(s[s.branch == "human"]) / gi["human"] - 1
                if ref == "pan_mean":
                    er = np.mean([g(s[s.branch == b]) / gi[b] - 1 for b in ["chimp", "bonobo"]])
                else:
                    er = g(s[s.branch == ref]) / gi[ref] - 1
                return eh / er
            pt = rho(None)
            bs = [rho(rng.choice(blocks, len(blocks), replace=True)) for _ in range(500)]
            lo, hi = np.percentile(bs, [2.5, 97.5])
            out.append(dict(D_Mb=D / 1e6, side=side, ref=ref, rho=pt, lo=lo, hi=hi,
                            T_if_Teff6_5=6.5 * (1 - pt), T_lo=6.5 * (1 - hi), T_hi=6.5 * (1 - lo)))
print(pd.DataFrame(out).round(2).to_string(index=False))
