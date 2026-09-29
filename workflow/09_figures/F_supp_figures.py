"""Supplementary Figures S2-S5 (same style as F_make_figures.py). usage: python F_supp_figures.py"""
import os

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import font_manager

R, OUT = "results", "results/figures"
for f in ("Regular", "Bold", "Italic"):
    p = f"/usr/share/fonts/truetype/lato/Lato-{f}.ttf"
    if os.path.exists(p):
        font_manager.fontManager.addfont(p)
MM = 1 / 25.4
C = dict(human="#0090a3", pan="#c0770f", gorilla="#4854b3", f2b="#a8204f", f2a="#e07a98", ends="#8a8f98",
         band="#dfe2e6", interior="#b9bec5", ink="#1f2328", ink2="#5b6068", grid="#eceef1")
TEAL = ["#b8e3e8", "#6cc3cf", "#0090a3", "#006675", "#00414b"]
mpl.rcParams.update({"font.family": ["Lato", "DejaVu Sans"], "font.size": 7, "axes.labelsize": 7, "xtick.labelsize": 6.5,
                     "ytick.labelsize": 6.5, "legend.fontsize": 6, "axes.linewidth": 0.6, "axes.edgecolor": C["ink2"],
                     "xtick.color": C["ink2"], "ytick.color": C["ink2"], "axes.spines.top": False,
                     "axes.spines.right": False, "axes.grid": True, "grid.color": C["grid"], "grid.linewidth": 0.5,
                     "axes.axisbelow": True, "legend.frameon": False, "pdf.fonttype": 42, "lines.linewidth": 1.2,
                     "mathtext.fontset": "custom", "mathtext.rm": "Lato", "mathtext.it": "Lato:italic"})


def panel(ax, letter, title):
    ax.text(-0.02, 1.04, letter, transform=ax.transAxes, fontsize=9, fontweight="bold", ha="right", va="bottom")
    ax.set_title(title, loc="left", fontsize=7, fontweight="bold", pad=4, x=0.02)


def save(fig, name):
    fig.savefig(f"{OUT}/{name}.pdf", bbox_inches="tight", pad_inches=0.02)
    fig.savefig(f"{OUT}/{name}.png", bbox_inches="tight", pad_inches=0.02, dpi=300)
    plt.close(fig)


# ---- S2: single-reference estimators on real ends (true T = 0) and the flanks
fig, axs = plt.subplots(1, 3, figsize=(180 * MM, 55 * MM), sharey=True, gridspec_kw=dict(wspace=0.12))
fl = pd.read_csv(f"{R}/A15/A15_flank_calibrated.tsv", sep="\t")
rng = np.random.default_rng(1)
for ax, D in zip(axs, (2, 5, 10)):
    for j, (ref, lab, col) in enumerate((("within_human", "human ends", C["human"]),
                                          ("panstem+chimp", "Pan stem + chimp.", C["pan"]),
                                          ("panstem+bonobo", "Pan stem + bonobo", C["pan"]))):
        te = pd.read_csv(f"{R}/A15/A15_Te_{ref}_{D}Mb.tsv", sep="\t", index_col=0).iloc[:, 0].clip(-8, 8)
        ax.scatter(te, j + rng.uniform(-0.18, 0.18, len(te)), s=6, color=C["ends"], alpha=0.8, lw=0)
        f = fl[(fl.reference == ref) & (fl.first_Mb == D) & (fl.side == "both")].iloc[0]
        ax.plot(f.T_raw, j, "s", color=C["f2b"], ms=4.5, mec="white", mew=0.5)
        ax.plot(np.median(te), j, "|", color=col, ms=12, mew=1.6)
    ax.axvline(0, color=C["ink"], lw=0.7)
    ax.set_yticks(range(3), ["human ends", "$\\it{Pan}$ (+chimp.)", "$\\it{Pan}$ (+bonobo)"])
    ax.set_xlim(-8.3, 8.3); ax.grid(axis="y", visible=False)
    ax.set_xlabel("estimated $T$ (Mya); truth = 0 for real ends")
    panel(ax, "abc"[[2, 5, 10].index(D)], f"first {D} Mb")
save(fig, "FigS2_single_reference")

# ---- S3: recombination calibration (A12b)
fig, axs = plt.subplots(1, 2, figsize=(150 * MM, 55 * MM), gridspec_kw=dict(wspace=0.35))
m = pd.read_csv(f"{R}/A12/A12b_models.tsv", sep="\t").set_index("model")
bins = ["0-2", "2-5", "5-10", "10-20"]
ax = axs[0]
x = np.arange(4)
ax.plot(x, [m.iloc[0][f"obs_{b}"] for b in bins], "o-", color=C["ink"], label="observed (real ends)")
for mod, lab, col, ls in (("dsb_pat_L", "predicted: paternal DSB", C["human"], "-"), ("dsb_mat_L", "predicted: maternal DSB", C["ends"], (0, (3, 1.5))),
                          ("co_pat_L", "predicted: paternal CO", TEAL[1], (0, (1, 1)))):
    ax.plot(x, [m.loc[mod, f"pred_{b}"] for b in bins], color=col, ls=ls, marker="s", ms=2.6, label=lab)
ax.set_xticks(x, [b + " Mb" for b in bins]); ax.set_ylabel("W→S log-odds excess over interior")
ax.legend(fontsize=5.8); ax.set_ylim(0, 0.55)
panel(ax, "a", "End signal predicted by present-day recombination")
ax = axs[1]
pe = pd.read_csv(f"{R}/A12/A12b_per_end.tsv", sep="\t")
ax.scatter(pe.pred, pe.obs, s=10, color=C["ends"], lw=0)
lim = [min(pe.pred.min(), pe.obs.min()) - 0.03, max(pe.pred.max(), pe.obs.max()) + 0.03]
ax.plot(lim, lim, color=C["ink2"], lw=0.6, ls=(0, (2, 2)))
r = np.corrcoef(pe.pred, pe.obs)[0, 1]
ax.text(0.97, 0.03, f"39 real ends, first 5 Mb\nr = {r:.2f}; residual SD\nbeyond sampling ≈ 0.095",
        transform=ax.transAxes, fontsize=5.8, va="bottom", ha="right", color=C["ink2"])
ax.set_xlabel("predicted excess (paternal DSB)"); ax.set_ylabel("observed excess")
panel(ax, "b", "Heterogeneity between ends")
save(fig, "FigS3_recomb_calibration")

# ---- S4: ape LD recombination vs own telomere distance
fig, axs = plt.subplots(1, 3, figsize=(180 * MM, 52 * MM), sharey=True, gridspec_kw=dict(wspace=0.12))
for ax, sp, col in zip(axs, ("chimp", "bonobo", "gorilla"), (C["pan"], C["pan"], C["gorilla"])):
    pr = pd.read_csv(f"{R}/A13/A13b_{sp}_profile.tsv", sep="\t")
    mids = {"0-2": 1, "2-5": 3.5, "5-10": 7.5, "10-15": 12.5, "15-20": 17.5, "20-25": 22.5, "25-30": 27.5}
    for st, lab, ls in (("uncapped_phys", "uncapped ends", "-"), ("capped_phys", "capped ends (physical distance)", (0, (3, 1.5))),
                        ("capped_from_boundary", "capped ends (from cap boundary)", (0, (1, 1)))):
        x = pr[pr.set == st]
        ax.plot([mids[b] for b in x.bin_Mb], x.median_r, color=col if st == "uncapped_phys" else C["ink2"], ls=ls, marker="o",
                ms=2.2, lw=1.0, label=lab)
    fl = pd.read_csv(f"{R}/A13/A13b_{sp}_flanks.tsv", sep="\t")
    for _, rr in fl[fl.human_dist_Mb == "0-5"].iterrows():
        if rr.sp_dist_end_Mb < 30:
            mk, c2 = ("o", C["f2a"]) if rr.side == "fl_2a" else ("s", C["f2b"])
            ax.plot(rr.sp_dist_end_Mb, rr.median_r, mk, color=c2, ms=4.5, mec="white", mew=0.5)
    ax.axhline(1, color=C["interior"], lw=0.8, ls=(0, (3, 2)))
    ax.set_xlabel(f"distance to the {sp} telomere (Mb)")
    panel(ax, "abc"[["chimp", "bonobo", "gorilla"].index(sp)], {"chimp": "chimpanzee", "bonobo": "bonobo", "gorilla": "gorilla"}[sp])
axs[0].set_ylabel("LD recombination (ρ / interior median)")
axs[0].legend(fontsize=5.6, loc="upper right")
save(fig, "FigS4_ape_recombination")

# ---- S5: the two flanks in detail (A16a/b)
fig, axs = plt.subplots(1, 2, figsize=(150 * MM, 50 * MM), gridspec_kw=dict(wspace=0.35))
a = pd.read_csv(f"{R}/A16/A16a_profile.tsv", sep="\t")
ax = axs[0]
for side, col, mk in (("2a", C["f2a"], "o"), ("2b", C["f2b"], "s")):
    x = a[(a.side == side) & (a.Mb < 20)]
    ax.plot(x.Mb + 0.5, x.frac.clip(-1.5, 2.5), color=col, marker=mk, ms=2.6, lw=0.9, label=f"flank {side}")
ax.axhline(1, color=C["ends"], lw=0.7, ls=(0, (3, 2))); ax.axhline(0, color=C["ink2"], lw=0.6)
ax.set_xlabel("distance from the junction (Mb)"); ax.set_ylabel("fraction of real-end excess")
ax.legend(fontsize=5.8, loc="lower left"); ax.set_ylim(-1.6, 2.6)
panel(ax, "a", "Retained signal along the flanks (1-Mb bins)")
q = pd.read_csv(f"{R}/A16/A16b_quality.tsv", sep="\t")
ax = axs[1]
for side, col, mk in (("2a", C["f2a"], "o"), ("2b", C["f2b"], "s")):
    x = q[q.side == side]
    ax.plot(x.Mb + 0.5, x.callable_frac, color=col, marker=mk, ms=2.6, lw=0.9, label=f"flank {side}")
x = q[q.side == "2a"]
ax.plot(x.Mb + 0.5, x.ends_callable_frac, color=C["ends"], lw=1.2, label="real ends (same distance)")
ax.set_ylim(0.6, 1.0); ax.set_xlabel("distance from the junction / telomere (Mb)"); ax.set_ylabel("callable fraction")
ax.legend(fontsize=5.8, loc="lower right")
panel(ax, "b", "Data quality near the junction")
save(fig, "FigS5_flanks_detail")
print("ok")
