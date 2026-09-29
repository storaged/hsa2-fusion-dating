"""A7c: W/S-composition correction of the A7 time-strata summary (the raw W->S shares understate the end excess
because ends are GC-rich). Composition (parent W/S) per region class from the A12 human-branch window table.
Output: log-odds excess over interior per stratum and region; fraction of the real-end excess carried by the flanks.
usage: python A7c_composition_correct.py results/A7/A7_summary.tsv results/A12/A12_windows.tsv.gz results/A7/A7c_corrected.tsv
"""
import sys
import numpy as np
import pandas as pd
summ_f, win_f, out = sys.argv[1:4]
w = pd.read_csv(win_f, sep="\t")
w["reg"] = np.select([w.cls.isin(["fusion_2a", "fusion_2b"]) & (w.dist < 10e6), (w.cls == "end") & (w.dist < 10e6),
                      w.cls == "interior"], [w.cls, "end", "interior"], "x")
comp = w[w.reg != "x"].groupby("reg")[["parent_W", "parent_S"]].sum()
comp["logWS"] = np.log(comp.parent_W / comp.parent_S)
s = pd.read_csv(summ_f, sep="\t").rename(columns={"region": "reg"}).merge(comp["logWS"], left_on="reg", right_index=True)
lg = lambda p: np.log(p / (1 - p))
for c in ("WS_share", "lo", "hi"):
    s["x_" + c] = lg(s[c]) - s.logWS
base = s[s.reg == "interior"].set_index("stratum").x_WS_share
for a, b in (("exc", "x_WS_share"), ("exc_lo", "x_lo"), ("exc_hi", "x_hi")):
    s[a] = s[b] - s.stratum.map(base)
endx = s[s.reg == "end"].set_index("stratum").exc
for a, b in (("frac_of_end", "exc"), ("frac_lo", "exc_lo"), ("frac_hi", "exc_hi")):
    s[a] = s[b] / s.stratum.map(endx)
s.to_csv(out, sep="\t", index=False)
print(s[["stratum", "reg", "n", "WS_share", "exc", "frac_of_end", "frac_lo", "frac_hi"]].round(3).to_string(index=False))
