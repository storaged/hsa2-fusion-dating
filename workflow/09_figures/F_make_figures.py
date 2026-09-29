"""Publication figures and table data for the HSA2 switch-off manuscript.

Reads results/* on amor, computes the statistics shown in the figures, writes results/figures/{Fig*.pdf,png} and
results/figures/stats.json + table TSVs. Style: compact panels, Inter 7 pt, recessive axes, direct labels,
validated palette (scripts: dataviz validator; all-pairs pass for teal/amber/indigo/crimson).
usage: python F_make_figures.py [--only 1,2,3,4,S1]
"""
import json
import os
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import font_manager
from scipy import stats

R = "results"
OUT = f"{R}/figures"
os.makedirs(OUT, exist_ok=True)
ONLY = set(sys.argv[sys.argv.index("--only") + 1].split(",")) if "--only" in sys.argv else {"1", "2", "3", "4", "S1"}
TS = 6.0
MM = 1 / 25.4

# ---------------- palette & style ----------------
C = dict(human="#0090a3", pan="#c0770f", gorilla="#4854b3", f2b="#a8204f", f2a="#e07a98",
         ends="#8a8f98", ends_band="#dfe2e6", interior="#b9bec5", ink="#1f2328", ink2="#5b6068", grid="#eceef1")
TEAL_RAMP = ["#b8e3e8", "#6cc3cf", "#0090a3", "#006675", "#00414b"]      # sequential (magnitude), one hue
SLATE_RAMP = ["#b9c0cc", "#7e889a", "#4b5566"]                          # neutral sequential for model scenarios
for _f in ("Regular", "Bold", "Semibold", "Italic", "Medium"):
    for _d in ("/usr/share/fonts/truetype/lato",):
        if os.path.exists(f"{_d}/Lato-{_f}.ttf"):
            font_manager.fontManager.addfont(f"{_d}/Lato-{_f}.ttf")
fams = {f.name for f in font_manager.fontManager.ttflist}
FONT = ["Lato", "DejaVu Sans"] if "Lato" in fams else ["DejaVu Sans"]
mpl.rcParams.update({
    "font.family": FONT, "font.size": 7, "axes.titlesize": 7.5, "axes.labelsize": 7, "xtick.labelsize": 6.5,
    "ytick.labelsize": 6.5, "legend.fontsize": 6.5, "axes.linewidth": 0.6, "axes.edgecolor": C["ink2"],
    "axes.labelcolor": C["ink"], "xtick.color": C["ink2"], "ytick.color": C["ink2"], "text.color": C["ink"],
    "xtick.major.width": 0.6, "ytick.major.width": 0.6, "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.color": C["grid"],
    "grid.linewidth": 0.5, "axes.axisbelow": True, "legend.frameon": False, "lines.linewidth": 1.2,
    "lines.markersize": 3.5, "pdf.fonttype": 42, "savefig.dpi": 300, "figure.dpi": 150,
    "axes.titleweight": "bold", "axes.titlelocation": "left", "axes.titlepad": 4,
    "mathtext.fontset": "custom", "mathtext.rm": "Lato", "mathtext.it": "Lato:italic", "mathtext.bf": "Lato:bold",
})
ST = {}   # statistics for the manuscript


def panel(ax, letter, title=None):
    ax.text(-0.02, 1.03, letter, transform=ax.transAxes, fontsize=9, fontweight="bold", ha="right", va="bottom")
    if title:
        ax.set_title(title, loc="left", fontsize=7, fontweight="bold", pad=4, x=0.02)


def note(ax, s, x=0.98, y=0.96, ha="right", va="top", **kw):
    ax.text(x, y, s, transform=ax.transAxes, ha=ha, va=va, fontsize=6, color=C["ink2"], **kw)


def save(fig, name):
    fig.savefig(f"{OUT}/{name}.pdf", bbox_inches="tight", pad_inches=0.02)
    fig.savefig(f"{OUT}/{name}.png", bbox_inches="tight", pad_inches=0.02, dpi=300)
    plt.close(fig)


def se_from_ci(lo, hi):
    return (hi - lo) / 3.92


# ======================= Figure 1: present-day status =======================
if "1" in ONLY:
    fig, axs = plt.subplots(1, 3, figsize=(180 * MM, 54 * MM), gridspec_kw=dict(width_ratios=[1.1, 1, 1.08], wspace=0.5))
    # a) paternal CO rate vs distance (deCODE, 1-Mb bins)
    ax = axs[0]
    w = pd.read_csv(f"{R}/A5/A5_windows_classified.tsv", sep="\t")
    w = w[w.sex == "pat"]
    e = w[w.cls.str.startswith("end")].groupby("k").cMperMb
    k = np.arange(15)
    med, q1, q3 = e.median().reindex(k), e.quantile(.25).reindex(k), e.quantile(.75).reindex(k)
    ax.fill_between(k + .5, q1, q3, color=C["ends_band"], lw=0, label="_")
    ax.plot(k + .5, med, color=C["ends"], lw=1.4)
    ax.text(10.8, med.iloc[10] + 0.35, "real ends\n(median, IQR)", color=C["ink2"], fontsize=6, ha="left", va="bottom")
    im = w[w.cls == "interior"].cMperMb.median()
    ax.axhline(im, color=C["interior"], lw=0.9, ls=(0, (3, 2)))
    ax.text(14.9, im - 0.1, "interior", color=C["ink2"], fontsize=6, ha="right", va="top")
    for side, col, mk in (("fusion_2a", C["f2a"], "o"), ("fusion_2b", C["f2b"], "s")):
        f = w[w.cls == side].groupby("k").cMperMb.mean().reindex(k)
        ax.plot(k + .5, f, color=col, marker=mk, ms=2.8, lw=1.1)
    ax.text(1.2, 0.05, "flank 2a", color=C["f2a"], fontsize=6.2, fontweight="bold")
    ax.text(5.2, 1.12, "flank 2b", color=C["f2b"], fontsize=6.2, fontweight="bold")
    ax.set_xlim(0, 15); ax.set_ylim(0, 5.6)
    ax.set_xlabel("distance from telomere / fusion junction (Mb)")
    ax.set_ylabel("paternal crossover rate (cM/Mb)")
    s5 = pd.read_csv(f"{R}/A5/A5_summary.tsv", sep="\t")
    s5 = s5[(s5.sex == "pat") & (s5.first_Mb == 5)]
    ST["fig1a"] = {r.side + "_" + r["var"]: dict(fusion=r.fusion, ends_mean=r.ends_mean, interior=r.interior_mean,
                                                 rel=r.rel_position, pct=r.pct_among_ends) for _, r in s5.iterrows()}
    n_ends = int(w[w.cls.str.startswith("end")].end_id.nunique())
    note(ax, f"flanks < all {n_ends} real ends\n(first 5 Mb; CO, NCO)\nrank $p$ < {1 / (n_ends + 1):.3f}", y=0.97)
    panel(ax, "a", "Recombination today")
    ST["n_real_ends_decode"] = n_ends

    # b) present-day B (AFR)
    ax = axs[1]
    b = pd.read_csv(f"{R}/A6/B_by_region_afr.tsv", sep="\t").set_index("region_class")
    q = pd.read_csv(f"{R}/A12c/B_by_dsbbin_afr.tsv", sep="\t").set_index("region_class")
    rows = [("real ends", b.loc["end"], C["ends"], "o"), ("interior", b.loc["interior"], C["interior"], "o"),
            ("flank 2a", b.loc["fusion_2a"], C["f2a"], "o"), ("flank 2b", b.loc["fusion_2b"], C["f2b"], "s")]
    qlab = {"q0-50": "0–50", "q50-75": "50–75", "q75-90": "75–90", "q90-97": "90–97", "q97-100": "97–100"}
    y = 0
    yt, yl = [], []
    for lab, r, col, mk in rows:
        ax.errorbar(r.B, y, xerr=[[r.B - r.B_lo], [r.B_hi - r.B]], fmt=mk, color=col, ms=3.5, elinewidth=1, capsize=0)
        yt.append(y); yl.append(lab); y -= 1
    y -= 0.6
    ax.text(-0.215, y + 0.5, "paternal DSB rate (percentile)", fontsize=6, color=C["ink2"], va="center")
    for i, (key, lab) in enumerate(qlab.items()):
        r = q.loc[key]
        ax.errorbar(r.B, y, xerr=[[r.B - r.B_lo], [r.B_hi - r.B]], fmt="o", color=TEAL_RAMP[i], ms=3.5,
                    elinewidth=1, mec=C["human"] if i == 0 else TEAL_RAMP[i], mew=0.5)
        yt.append(y); yl.append(lab); y -= 1
    ax.axvline(0, color=C["ink2"], lw=0.6)
    ax.set_yticks(yt, yl); ax.set_ylim(y + 0.4, 0.6); ax.grid(axis="y", visible=False)
    ax.set_xlim(-0.22, 0.32)
    ax.set_xlabel("present-day gBGC strength $B = 4N_e b$")
    fb = b.loc[["fusion_2a", "fusion_2b"]]
    zs = {}
    for nm, ref in (("ends", "end"), ("interior", "interior")):
        for s_ in ("fusion_2a", "fusion_2b"):
            d_ = b.loc[s_, "B"] - b.loc[ref, "B"]
            se = np.hypot(se_from_ci(b.loc[s_, "B_lo"], b.loc[s_, "B_hi"]), se_from_ci(b.loc[ref, "B_lo"], b.loc[ref, "B_hi"]))
            zs[f"{s_}_vs_{nm}"] = dict(diff=d_, z=d_ / se, p=2 * stats.norm.sf(abs(d_ / se)))
    ST["fig1b"] = dict(B={k_: b.loc[k_, ["B", "B_lo", "B_hi"]].to_dict() for k_ in b.index}, tests=zs,
                       dsb={k_: q.loc[k_, ["B", "B_lo", "B_hi"]].to_dict() for k_ in q.index})
    pmax_end = max(zs["fusion_2a_vs_ends"]["p"], zs["fusion_2b_vs_ends"]["p"])
    note(ax, f"flanks vs ends\n$p$ ≤ {pmax_end:.0e}\nvs interior: n.s.", x=0.99, y=0.8)
    panel(ax, "b", "gBGC strength today")

    # c) retained fraction by time stratum (A7c)
    ax = axs[2]
    s7 = pd.read_csv(f"{R}/A7/A7c_corrected.tsv", sep="\t")
    strata = [("poly", "poly-\nmorphic"), ("fixed_shared", "fixed &\nshared with\narchaics"), ("ALL", "all")]
    for j, (st_, lab) in enumerate(strata):
        for dx, side, col, mk in ((-0.13, "fusion_2a", C["f2a"], "o"), (0.13, "fusion_2b", C["f2b"], "s")):
            r = s7[(s7.stratum == st_) & (s7.reg == side)].iloc[0]
            ax.errorbar(j + dx, r.frac_of_end, yerr=[[r.frac_of_end - r.frac_lo], [r.frac_hi - r.frac_of_end]],
                        fmt=mk, color=col, ms=3.5, elinewidth=1)
    ax.axhline(0, color=C["ink2"], lw=0.6); ax.axhline(1, color=C["ends"], lw=0.8, ls=(0, (3, 2)))
    ax.text(2.45, 1.02, "real ends", color=C["ink2"], fontsize=6, ha="right", va="bottom")
    ax.text(2.45, 0.02, "interior", color=C["ink2"], fontsize=6, ha="right", va="bottom")
    ax.set_xticks(range(3), [l for _, l in strata], fontsize=6.2)
    ax.set_xlim(-0.5, 2.5); ax.set_ylim(-0.65, 1.2); ax.grid(axis="x", visible=False)
    ax.set_ylabel("fraction of real-end W→S excess")
    ax.plot([], [], "o", color=C["f2a"], label="flank 2a"); ax.plot([], [], "s", color=C["f2b"], label="flank 2b")
    ax.legend(loc="lower right", fontsize=6, handletextpad=0.2, borderaxespad=0.1)
    ST["fig1c"] = s7[s7.reg.str.startswith("fusion")][["stratum", "reg", "n", "frac_of_end", "frac_lo", "frac_hi"]].to_dict("records")
    panel(ax, "c", "Signal by time stratum")
    save(fig, "Fig2_present_day")

# ======================= Figure 2: pooled model =======================
if "2" in ONLY:
    sys.argv = ["A4_pooled_model.py", f"{R}/A3/windows.tsv.gz", f"{R}/A4", "chimp"]
    sys.path[:0] = ["scripts", "workflow/06_pooled_model"]   # flat install or repository layout
    import A4_pooled_model as M  # noqa: E402
    fitrow = pd.read_csv(f"{R}/A4/A4_fit_chimp_human+pan+gorilla.tsv", sep="\t").iloc[0]
    That = fitrow.T_hat
    _, x = M.fit({"F2a": 1 - That / TS, "F2b": 1 - That / TS})
    mu, sig, K, tau2 = M.unpack(x)
    O = M.O
    ctr = np.array([0.5, 1.5, 3.5, 7.5, 12.5, 17.5, 22.5])
    fig = plt.figure(figsize=(180 * MM, 118 * MM))
    gs = fig.add_gridspec(2, 2, wspace=0.34, hspace=0.95, width_ratios=[1, 1])
    # a) kernels and flank observations
    ax = fig.add_subplot(gs[0, 0])
    for l, lab in (("human", "human"), ("pan", "$\\it{Pan}$ (stem + chimp.)"), ("gorilla", "gorilla")):
        li = M.LI[l]
        ob = O[(O.lin == l) & ~O.unit.str.startswith("F")]
        g = ob.groupby("bin").y
        ax.errorbar(ctr[g.median().index] + {"human": -0.25, "pan": 0, "gorilla": 0.25}[l], g.median(),
                    yerr=[g.median() - g.quantile(.25), g.quantile(.75) - g.median()], fmt="o", ms=2.6,
                    color=C[l], elinewidth=0.7, alpha=0.9)
        ax.plot(ctr, mu * K[li], color=C[l], lw=1.2, label=lab + " real ends")
    rho = 1 - That / TS
    ax.plot(ctr[:4], mu * K[M.LI["human"]][:4] * rho, color=C["f2b"], lw=1.1, ls=(0, (3, 1.5)),
            label=f"human flanks, fit: ρ = {1 - That / TS:.2f}")
    for u, col, mk in (("F2a", C["f2a"], "o"), ("F2b", C["f2b"], "s")):
        ob = O[(O.lin == "human") & (O.unit == u)]
        ax.errorbar(ctr[ob.bin] + (0.35 if u == "F2b" else -0.35), ob.y, yerr=1.96 * np.sqrt(ob.v + tau2[0]), fmt=mk,
                    ms=3.2, color=col, elinewidth=0.8, mec="white", mew=0.4)
    ax.plot([], [], "o", color=C["f2a"], ms=3, label="flank 2a (human)"); ax.plot([], [], "s", color=C["f2b"], ms=3, label="flank 2b (human)")
    ax.legend(loc="upper right", fontsize=5.8, handlelength=1.8, borderaxespad=0.1, labelspacing=0.3)
    ax.set_xlim(0, 25); ax.set_ylim(-0.2, 0.8)
    ax.set_xlabel("distance to the lineage's own telomere (Mb)")
    ax.set_ylabel("W→S log-odds excess")
    note(ax, f"39 real ends; between-end CV {sig / mu:.2f}", x=0.98, y=0.03, va="bottom")
    panel(ax, "a", "Telomere kernels by lineage (pooled model)")
    ST["fig2a"] = dict(mu=mu, sigma=sig, cv=sig / mu, K=K.tolist(), rho=rho, tau=np.sqrt(tau2).tolist())
    # b) profile likelihoods
    ax = fig.add_subplot(gs[0, 1])
    confs = [("chimp_human+pan+gorilla", "A4", "main (CHM13; Pan = chimpanzee)", C["ink"], 1.6, "-"),
             ("bonobo_human+pan+gorilla", "A4", "Pan = bonobo", SLATE_RAMP[2], 0.9, "-"),
             ("chimp_human+gorilla", "A4", "no Pan", SLATE_RAMP[1], 0.9, (0, (3, 1.5))),
             ("chimp_human", "A4", "human only", SLATE_RAMP[0], 0.9, (0, (1, 1))),
             ("chimp_human+pan+gorilla", "A4_h16a", "haplotype set h16a", SLATE_RAMP[2], 0.9, (0, (4, 1, 1, 1))),
             ("chimp_human+pan+gorilla", "A4_h16b", "haplotype set h16b", SLATE_RAMP[1], 0.9, (0, (4, 1, 1, 1, 1, 1)))]
    for tag, d_, lab, col, lw, ls in confs:
        p = pd.read_csv(f"{R}/{d_}/A4_profile_{tag}.tsv", sep="\t")
        p = p[p["T"] >= 0]
        ax.plot(p["T"], p.nll - p.nll.min(), color=col, lw=lw, ls=ls, label=lab)
    ax.axhline(1.92, color=C["f2b"], lw=0.7, ls=(0, (2, 2)))
    ax.text(2.7, 2.2, "95% (Δ = 1.92)", color=C["f2b"], fontsize=6, ha="center")
    ax.axvspan(5, 6, color="#f3e3e8", lw=0, zorder=0)
    ax.text(5.5, 12.3, ">5 Mya (Yang et al.)", ha="center", fontsize=6, color=C["f2b"], clip_on=False)
    ax.set_xlim(0, 6); ax.set_ylim(0, 12)
    ax.set_xlabel("switch-off time $T_{off}$ (Mya)"); ax.set_ylabel("profile −log L (relative)")
    ax.legend(loc="upper center", fontsize=5.8, handlelength=2.4, borderaxespad=0.0, ncol=3, columnspacing=0.9,
              bbox_to_anchor=(0.45, -0.25))
    panel(ax, "b", "Profile likelihood of $T_{off}$")
    # c) forest plot
    ax = fig.add_subplot(gs[1, 0])
    fr = []
    for tag, d_, lab in (("chimp_human+pan+gorilla", "A4", "pooled model, main"), ("bonobo_human+pan+gorilla", "A4", "Pan = bonobo"),
                         ("chimp_human+gorilla", "A4", "no Pan lineage"), ("chimp_human", "A4", "human lineage only"),
                         ("chimp_human+pan+gorilla", "A4_h16a", "haplotype set h16a"), ("chimp_human+pan+gorilla", "A4_h16b", "haplotype set h16b")):
        r = pd.read_csv(f"{R}/{d_}/A4_fit_{tag}.tsv", sep="\t").iloc[0]
        fr.append((lab, r.T_hat, r.T_lo, r.T_hi, C["ink"] if lab.startswith("pooled") else C["ink2"], "o"))
    r = pd.read_csv(f"{R}/A4/A4_fit_chimp_human+pan+gorilla.tsv", sep="\t").iloc[0]
    fr.append(("flank 2a alone", r.T_F2a, r.T_F2a_lo, r.T_F2a_hi, C["f2a"], "o"))
    fr.append(("flank 2b alone", r.T_F2b, r.T_F2b_lo, r.T_F2b_hi, C["f2b"], "s"))
    a15 = pd.read_csv(f"{R}/A15/A15_flank_calibrated.tsv", sep="\t")
    for ref, lab in (("within_human", "single reference: human ends"), ("panstem+chimp", "single reference: Pan orthologue")):
        x_ = a15[(a15.reference == ref) & (a15.side == "both") & (a15.first_Mb == 10)].iloc[0]
        fr.append((lab, x_.T_debiased, x_.T_emp_lo, x_.T_emp_hi, C["ends"], "D"))
    for i, (lab, t, lo_, hi_, col, mk) in enumerate(fr):
        yy = -i - (0.5 if i >= 6 else 0) - (0.5 if i >= 8 else 0)
        ax.plot([lo_, hi_], [yy, yy], color=col, lw=1.1, solid_capstyle="round")
        ax.plot(t, yy, mk, color=col, ms=3.8 if i else 4.5, mec="white", mew=0.4)
        ax.text(-0.1, yy, lab, ha="right", va="center", fontsize=6.2, color=C["ink"], transform=ax.get_yaxis_transform())
    ax.axvspan(5, 7.3, color="#f3e3e8", lw=0, zorder=0)
    ax.axvline(0.9, color=C["ink2"], lw=0.6, ls=(0, (1, 1.5)))
    ax.text(0.9, 1.35, "2022 estimate", fontsize=5.8, color=C["ink2"], va="bottom", ha="center")
    ax.text(6.15, 1.35, ">5 Mya (Yang et al.)", fontsize=5.8, color=C["f2b"], va="bottom", ha="center")
    ax.set_yticks([]); ax.grid(axis="y", visible=False); ax.spines["left"].set_visible(False)
    ax.set_xlim(-1.5, 7.3); ax.set_ylim(-len(fr) - 0.6, 2.2)
    ax.set_xlabel("$T_{off}$ (Mya), estimate and 95% interval")
    panel(ax, "c", "Estimates across configurations")
    ax.title.set_position((-0.62, 1))
    ST["fig2c"] = [dict(label=a, T=b_, lo=c_, hi=d_) for a, b_, c_, d_, _, _ in fr]
    # d) calibration caterpillar
    ax = fig.add_subplot(gs[1, 1])
    ps = pd.read_csv(f"{R}/A4/A4_pseudo_chimp_human+pan+gorilla.tsv", sep="\t")
    s1 = ps[ps.n == 1].sort_values("T_hat").reset_index(drop=True)
    cov = s1.covers0.mean(); nc = int(s1.covers0.sum())
    ci = stats.binomtest(nc, len(s1)).proportion_ci(method="exact")
    s2 = ps[ps.n == 2]
    for i, r in s1.iterrows():
        col = C["ends"] if r.covers0 else C["ink"]
        ax.plot([r.T_lo, r.T_hi], [i, i], color=col, lw=0.9)
        ax.plot(r.T_hat, i, "o", color=col, ms=2.4)
    ax.axvline(0, color=C["ink"], lw=0.8)
    yF = len(s1) + 1.5
    ax.plot([fitrow.T_lo, fitrow.T_hi], [yF, yF], color=C["f2b"], lw=1.4)
    ax.plot(fitrow.T_hat, yF, "s", color=C["f2b"], ms=4)
    ax.text(fitrow.T_hi + 0.15, yF, "fusion flanks", color=C["f2b"], fontsize=6, va="center")
    ax.set_yticks([]); ax.grid(axis="y", visible=False)
    ax.set_xlim(-3.2, 6); ax.set_ylim(-1, yF + 1.5)
    ax.set_xlabel("estimated $T$ (Mya); truth = 0 for real ends")
    ax.set_ylabel("real ends as mock fusions")
    note(ax, f"covers $T$ = 0\n{nc}/{len(s1)} = {cov:.2f}\n(95% CI\n{ci.low:.2f}–{ci.high:.2f})\n\n"
             f"pairs\n{int(s2.covers0.sum())}/{len(s2)} = {s2.covers0.mean():.2f}\n\nmedian $\\hat{{T}}$\n{s1.T_hat.median():.1f}",
         x=0.995, y=0.03, va="bottom")
    ST["fig2d"] = dict(cov_single=cov, n_single=len(s1), n_cov=nc, ci=[ci.low, ci.high], cov_pair=s2.covers0.mean(),
                       n_pair=len(s2), median_T=s1.T_hat.median(), width=(s1.T_hi - s1.T_lo).mean())
    panel(ax, "d", "Calibration on real ends (true T = 0)")
    save(fig, "Fig3_pooled_model")

# ======================= Figure 3: lineage sorting =======================
if "3" in ONLY:
    fig, axs = plt.subplots(1, 2, figsize=(180 * MM, 56 * MM), gridspec_kw=dict(wspace=0.3))
    ax = axs[0]
    d16 = pd.read_csv(f"{R}/A16/A16d_ils_junction.tsv", sep="\t")
    d16["mid"] = (d16.Mb_lo + d16.Mb_hi) / 2
    a8 = pd.read_csv(f"{R}/A8/A8_profile.tsv", sep="\t")
    ae = a8[(a8.cls == "end") & (a8.dbin >= 0) & (a8.dbin < 10)]
    ax.plot(ae.dbin + 0.5, ae.fold_vs_genome, color=C["ends"], lw=1.3)
    ax.plot([], [], color=C["ends"], lw=1.3, label="real ends (mean)")
    ax.axhline(1, color=C["interior"], lw=0.8, ls=(0, (3, 2)))
    for side, col, mk in (("2a", C["f2a"], "o"), ("2b", C["f2b"], "s")):
        x_ = d16[d16.side == side]
        ax.plot(x_.mid, x_.fold, color=col, marker=mk, ms=2.4, lw=0.9, label=f"flank {side}")
    ax.legend(loc="upper left", fontsize=6, handlelength=1.8, borderaxespad=0.1)
    ax.set_xscale("log"); ax.set_xlim(0.04, 11)
    ax.set_xticks([0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10], ["0.05", "0.1", "0.2", "0.5", "1", "2", "5", "10"])
    ax.set_ylim(0.0, 3.7)
    ax.set_xlabel("distance from the fusion junction (Mb)")
    ax.set_ylabel("ILS discordance (fold over genome)")
    one = d16[(d16.Mb_hi <= 10) & (d16.Mb_lo >= 1)]
    a_, b_ = one[one.side == "2a"].fold, one[one.side == "2b"].fold
    mw = stats.mannwhitneyu(b_, a_, alternative="greater")
    pe = pd.read_csv(f"{R}/A8/A8_per_end_first5Mb.tsv", sep="\t")
    s5b = d16[(d16.side == "2b") & (d16.Mb_hi <= 5)]
    D2b = np.average(s5b.D, weights=s5b.n_inf)
    s5a = d16[(d16.side == "2a") & (d16.Mb_hi <= 5)]
    D2a = np.average(s5a.D, weights=s5a.n_inf)
    pct_b, pct_a = (pe.D < D2b).mean(), (pe.D < D2a).mean()
    note(ax, f"2b > 2a (1–10 Mb bins): Mann–Whitney $p$ = {mw.pvalue:.1e}\n"
             f"first 5 Mb vs {len(pe)} real ends: 2b {100 * pct_b:.0f}th, 2a {100 * pct_a:.0f}th percentile", y=0.02, va="bottom")
    ST["fig3a"] = dict(mw_p=mw.pvalue, pct_2b=pct_b, pct_2a=pct_a, D2b=D2b, D2a=D2a, n_ends=len(pe))
    panel(ax, "a", "Observed: lineage sorting around the junction")
    ax = axs[1]
    hp = pd.read_csv(f"{R}/A10/A10c_Hpoly.tsv", sep="\t")
    hs = pd.read_csv(f"{R}/A10/A10c_Hsub.tsv", sep="\t")
    ax.axhspan(1.3, 1.6, color="#f3dbe2", lw=0, zorder=0)
    ax.text(0.0012, 1.55, "observed 2b plateau", color=C["f2b"], fontsize=6, va="top")
    ax.axhline(1.0, color=C["f2a"], lw=1.0, ls=(0, (1, 1)))
    ax.text(0.0012, 0.96, "observed 2a", color=C["f2a"], fontsize=6, va="top")
    for i, dl in enumerate((1.0, 2.0, 3.0)):
        for c_, ls in ((0.001, "-"), (1.0, (0, (3, 1.5)))):
            x_ = hp[(hp.Delta_Myr == dl) & (hp.c == c_) & (hp.y_cM > 0)]
            ax.plot(x_.y_cM, x_.fold, color=SLATE_RAMP[i], lw=1.1, ls=ls)
        x0 = hp[(hp.Delta_Myr == dl) & (hp.c == 0.001) & (hp.y_cM == 0.001)].fold.iloc[0]
        ax.text(0.0011, x0 + 0.02, f"Δ = {dl:g} Myr", color=SLATE_RAMP[i], fontsize=6, ha="left", va="bottom")
    for _, r in hs[hs.k.isin([1.5, 2.0])].iterrows():
        ax.plot([0.02, 12], [r.fold, r.fold], color=C["gorilla"], lw=1.0, ls=(0, (5, 1.5)))
        ax.text(12.5, r.fold, f"k = {r.k:g}", color=C["gorilla"], fontsize=6, va="center")
    ax.set_xscale("log"); ax.set_xlim(0.0003, 30)
    ax.set_ylim(0.8, 2.35)
    ax.set_xlabel("genetic distance from the junction (cM)")
    ax.set_ylabel("predicted ILS fold")
    ax.plot([], [], color=SLATE_RAMP[2], lw=1.1, label="polymorphism, recomb. 0.1% of normal")
    ax.plot([], [], color=SLATE_RAMP[2], lw=1.1, ls=(0, (3, 1.5)), label="polymorphism, normal recomb.")
    ax.plot([], [], color=C["gorilla"], lw=1.0, ls=(0, (5, 1.5)), label="larger ancestral $N_e$ (k-fold)")
    ax.legend(loc="upper right", fontsize=5.8, handlelength=2.2, borderaxespad=0.1)
    panel(ax, "b", "Simulated: what could produce it")
    save(fig, "Fig4_ils")

# ======================= Figure 4: validation =======================
if "4" in ONLY:
    fig, axs = plt.subplots(1, 3, figsize=(180 * MM, 52 * MM), gridspec_kw=dict(wspace=0.42))
    ax = axs[0]
    a10 = pd.read_csv(f"{R}/A10/A10a_bias_N400_Ts8.tsv", sep="\t")
    for (cnt, kind), col, ls, lab in ((("all", "step"), C["human"], "-", "genome, step"),
                                      (("all", "ramp"), C["human"], (0, (3, 1.5)), "genome, gradual"),
                                      (("fixed_only", "step"), C["ends"], "-", "fixed, step"),
                                      (("fixed_only", "ramp"), C["ends"], (0, (3, 1.5)), "fixed, gradual")):
        for B_, mk in ((1.0, "o"), (5.0, "^")):
            x_ = a10[(a10.counts == cnt) & (a10.kind == kind) & (a10.B == B_)].sort_values("Toff_in_N")
            ax.plot(x_.Toff_in_N, x_.bias_in_N, color=col, ls=ls, marker=mk, ms=2.6, lw=1,
                    label=lab if B_ == 1.0 else None)
    ax.axhline(0, color=C["ink"], lw=0.6)
    ax.set_ylim(-1.3, 0.72)
    ax.set_xlabel("true $T_{off}$ ($N$ generations)"); ax.set_ylabel("bias of $\\hat{T}$ ($N$ generations)")
    ax.legend(loc="center", bbox_to_anchor=(0.5, 0.5), fontsize=5.6, handlelength=1.8, ncol=2, columnspacing=0.8, handletextpad=0.4)
    note(ax, "exact Wright–Fisher, $T_s = 8N$\ncircles $B$ = 1, triangles $B$ = 5", y=0.99)
    panel(ax, "a", "Estimator bias (exact)")
    ST["fig4a"] = dict(max_abs_bias_all_step=a10[(a10.counts == "all") & (a10.kind == "step")].bias_in_N.abs().max(),
                       max_abs_bias_all_ramp=a10[(a10.counts == "all") & (a10.kind == "ramp")].bias_in_N.abs().max(),
                       fixed_range=[a10[a10.counts == "fixed_only"].bias_in_N.min(), a10[a10.counts == "fixed_only"].bias_in_N.max()])
    ax = axs[1]
    if os.path.exists(f"{R}/A10b/A10b_estimates.tsv"):
        sb = pd.read_csv(f"{R}/A10b/A10b_estimates.tsv", sep="\t")
        TSg = 4000
        for i, (bias, col) in enumerate(zip(sorted(sb.bias.unique()), [TEAL_RAMP[1], TEAL_RAMP[2], TEAL_RAMP[4]])):
            x_ = sb[sb.bias == bias].sort_values("T_true")
            dx = (i - 1) * 0.02
            ax.errorbar(x_.T_true / TSg + dx, x_.T_hat / TSg, yerr=[(x_.T_hat - x_.T_lo) / TSg, (x_.T_hi - x_.T_hat) / TSg],
                        fmt="o", color=col, ms=3, elinewidth=0.9, label=f"conversion bias {bias:g}")
        ax.plot([-0.1, 1.1], [-0.1, 1.1], color=C["ink2"], lw=0.7, ls=(0, (2, 2)))
        ax.legend(loc="upper left", fontsize=5.8)
        note(ax, f"covers truth: {sb.covered.sum()}/{len(sb)}\nmean |bias| {sb.bias_gen.abs().mean() / TSg:.3f} $T_s$",
             x=0.98, y=0.2)
        ST["fig4b"] = dict(coverage=sb.covered.mean(), n=len(sb), n_cov=int(sb.covered.sum()),
                           mean_abs_bias_Ts=sb.bias_gen.abs().mean() / TSg,
                           end_excess={str(k_): v for k_, v in sb.groupby("bias").end_excess.mean().items()})
    ax.set_xlim(-0.1, 1.1); ax.set_ylim(-0.3, 1.3)
    ax.set_xlabel("true $T_{off}/T_s$"); ax.set_ylabel("estimated $T_{off}/T_s$")
    panel(ax, "b", "Forward simulation (SLiM)")
    ax = axs[2]
    f4b = f"{R}/A4/A4b_sides_chimp.tsv"
    if os.path.exists(f4b):
        sd = pd.read_csv(f4b, sep="\t")
        fl, null = sd[sd.is_flank].iloc[0], sd[~sd.is_flank]
        bins = np.linspace(0, max(null.LR.max(), fl.LR) * 1.05, 18)
        ax.hist(null.LR, bins=bins, color=C["ends_band"], edgecolor=C["ends"], lw=0.5)
        ax.axvline(fl.LR, color=C["f2b"], lw=1.4)
        pe_ = (1 + (null.LR >= fl.LR).sum()) / (1 + len(null))
        ax.text(fl.LR - 0.15, ax.get_ylim()[1] * 0.95, f"flanks 2a vs 2b\nLR = {fl.LR:.1f}\nempirical $p$ = {pe_:.2g}\n($\\chi^2_1$ $p$ = {fl.p_chi2:.2g})",
                color=C["f2b"], fontsize=6, va="top", ha="right")
        ax.set_xlabel("LR, separate vs shared $T$"); ax.set_ylabel(f"real-end pairs (n = {len(null)})")
        ST["fig4c"] = dict(LR=fl.LR, p_chi2=fl.p_chi2, p_emp=pe_, n_null=len(null), T_2a=fl.T_a, T_2b=fl.T_b,
                           T_shared=fl.T_shared, null_median=null.LR.median(), null_q95=null.LR.quantile(.95),
                           null_frac_chi2_sig=(null.p_chi2 < 0.05).mean())
    panel(ax, "c", "One event on both sides?")
    save(fig, "Fig5_validation")

# ======================= Figure S1: orthologue end strength =======================
if "S1" in ONLY:
    fig, axs = plt.subplots(1, 2, figsize=(150 * MM, 48 * MM), gridspec_kw=dict(wspace=0.55))
    ax = axs[0]
    rows = []
    for cell, path in (("CHM13 set", f"{R}/A16/A16c_orthologue_strength.tsv"), ("h16a", f"{R}/A16_h16a/A16c_orthologue_strength.tsv"),
                       ("h16b", f"{R}/A16_h16b/A16c_orthologue_strength.tsv")):
        c16 = pd.read_csv(path, sep="\t")
        c16 = c16[(c16.lineage != "panstem") & (c16.first_Mb == 10) & c16.pct_among_ends.notna()]
        for _, r in c16.iterrows():
            rows.append((cell, r.lineage, r.side, r.pct_among_ends))
    Rr = pd.DataFrame(rows, columns=["cell", "lineage", "side", "pct"])
    labs = []
    for i, ((lin, side), g) in enumerate(Rr.groupby(["lineage", "side"], sort=False)):
        col = C["pan"] if lin in ("chimp", "bonobo") else C["gorilla"]
        ax.plot(g.pct * 100, [i] * len(g), "o" if side == "2a" else "s", color=col, ms=3.3, alpha=0.9, mec="white", mew=0.4)
        labs.append(f"{lin} {side.upper()}")
    ax.axvline(50, color=C["ends"], lw=0.7, ls=(0, (3, 2)))
    ax.set_yticks(range(len(labs)), labs); ax.grid(axis="y", visible=False)
    ax.set_xlim(0, 100); ax.set_xlabel("percentile among that lineage's real ends")
    note(ax, "orthologues of the flanks, first 10 Mb\n(three assembly/haplotype sets)", x=0.02, y=0.99, ha="left", va="top")
    panel(ax, "a", "Substitution end-strength")
    ax = axs[1]
    rr = []
    for sp in ("chimp", "bonobo", "gorilla"):
        pr = pd.read_csv(f"{R}/A13/A13b_{sp}_profile.tsv", sep="\t")
        pk = pr[pr.set == "all_ends_phys"].median_r.max()
        fl = pd.read_csv(f"{R}/A13/A13b_{sp}_flanks.tsv", sep="\t")
        for _, r in fl[fl.human_dist_Mb == "0-5"].iterrows():
            if r.sp_dist_end_Mb < 30:
                rr.append((sp, r.side.replace("fl_", "").upper(), (r.median_r - 1) / (pk - 1)))
    Z = pd.DataFrame(rr, columns=["sp", "side", "frac"])
    for i, r in Z.iterrows():
        col = C["gorilla"] if r.sp == "gorilla" else C["pan"]
        ax.barh(i, r.frac, color=col, height=0.62)
        ax.text(r.frac + 0.03, i, f"{r.frac:.2f}", va="center", fontsize=6, color=C["ink"])
    ax.axvline(1, color=C["ends"], lw=0.7, ls=(0, (3, 2)))
    ax.set_yticks(range(len(Z)), [f"{a} {b}" for a, b in zip(Z.sp, Z.side)]); ax.grid(axis="y", visible=False)
    ax.set_xlim(0, 1.8); ax.set_xlabel("LD recombination excess, fraction of peak end")
    panel(ax, "b", "Recombination today at the orthologues")
    save(fig, "FigS1_orthologues")
    ST["figS1"] = dict(pct=Rr.to_dict("records"), recomb=Z.to_dict("records"))

prev = {}
if os.path.exists(f"{OUT}/stats.json"):
    prev = json.load(open(f"{OUT}/stats.json"))
prev.update(json.loads(json.dumps(ST, default=float)))
json.dump(prev, open(f"{OUT}/stats.json", "w"), indent=1, default=float)
print("figures:", sorted(f for f in os.listdir(OUT) if f.endswith(".pdf")))
