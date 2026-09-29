"""Figure 1 (conceptual overview), Nature style: panel letters only, glyphs instead of labels; the caption explains.
Glyphs: capped rods = ancestral acrocentric chromosomes (black caps = telomeres); rod with a red band = fused
chromosome 2 (junction); rod with an open waist = relic centromere; stacked blocks = segmental duplications;
crossed tree = lineage sorting; chip = W->S substitutions; data/model/test icons in d.
usage: python F_schematic.py results/figures
"""
import os
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch, Rectangle

OUT = sys.argv[1] if len(sys.argv) > 1 else "results/figures"
os.makedirs(OUT, exist_ok=True)
for f in ("Regular", "Bold", "Italic"):
    p = f"/usr/share/fonts/truetype/lato/Lato-{f}.ttf"
    if os.path.exists(p):
        font_manager.fontManager.addfont(p)
MM = 1 / 25.4
C = dict(human="#0090a3", pan="#c0770f", gorilla="#4854b3", f2b="#a8204f", ends="#8a8f98", band="#e6e8eb",
         rod="#cfd3d9", interior="#b9bec5", ink="#1f2328", ink2="#5b6068", orang="#9aa0a8", pink="#f3dbe2")
mpl.rcParams.update({"font.family": ["Lato", "DejaVu Sans"], "font.size": 7, "axes.linewidth": 0.6,
                     "axes.edgecolor": C["ink2"], "text.color": C["ink"], "xtick.color": C["ink2"],
                     "ytick.color": C["ink2"], "xtick.labelsize": 6.5, "ytick.labelsize": 6.5,
                     "axes.spines.top": False, "axes.spines.right": False, "pdf.fonttype": 42,
                     "mathtext.fontset": "custom", "mathtext.rm": "Lato", "mathtext.it": "Lato:italic"})
FIG_W, FIG_H = 180 * MM, 118 * MM


# ---------------------------------------------------------------- glyphs (unit box, equal aspect)
def icon(fig, rect, equal=True):
    a = fig.add_axes(rect)
    a.set_xlim(0, 1); a.set_ylim(0, 1)
    if equal:
        a.set_aspect("equal", adjustable="datalim")
    a.axis("off")
    return a


def chrom_icon(fig, x, y, w, kind, **kw):
    """chromosome glyph with a fixed physical thickness (figure fraction box of height RH)"""
    a = icon(fig, [x, y, w, RH], equal=False)
    if kind == "ancestral":
        rod(a, 0.0, 0.5, 0.47, caps=[0.39]); rod(a, 0.53, 0.5, 0.47, caps=[0.61])
    elif kind == "fused":
        rod(a, 0.0, 0.5, 1.0, junction=kw.get("junction", 0.44))
    elif kind == "relic":
        rod(a, 0.0, 0.5, 1.0, junction=0.28, waist=0.66)
    elif kind == "end":
        rod(a, 0.0, 0.5, 1.0, caps=[0.09])
    elif kind == "plain":
        rod(a, 0.0, 0.5, 1.0)
    return a


RH = 0.022   # chromosome glyph thickness (figure fraction)


def rod(a, x0, y, w, h=None, caps=(), junction=None, waist=None, col=None, lw=5.5):
    """chromosome drawn as a round-capped thick line (thickness in points, independent of box aspect)"""
    e = 0.06
    a.plot([x0 + e, x0 + w - e], [y, y], color=col or C["rod"], lw=lw, solid_capstyle="round", clip_on=False)
    for cx in caps:
        a.plot([cx - 0.022, cx + 0.022], [y, y], color=C["ink"], lw=lw, solid_capstyle="butt", clip_on=False)
    if junction is not None:
        a.plot([junction - 0.014, junction + 0.014], [y, y], color=C["f2b"], lw=lw * 1.6, solid_capstyle="butt", clip_on=False)
    if waist is not None:
        a.plot(waist, y, "o", ms=lw * 0.75, mfc="white", mec=C["ink2"], mew=0.6)


def g_ancestral(a):
    rod(a, 0.02, 0.5, 0.44, caps=[0.43]); rod(a, 0.54, 0.5, 0.44, caps=[0.57])


def g_fused(a):
    rod(a, 0.05, 0.5, 0.9, junction=0.44)


def g_relic(a):
    rod(a, 0.05, 0.5, 0.9, junction=0.3, waist=0.66)


def g_sd(a):
    for x, y, c in ((0.08, 0.62, C["ends"]), (0.36, 0.62, C["f2b"]), (0.22, 0.38, C["ends"]),
                    (0.5, 0.38, C["ends"]), (0.64, 0.62, C["ends"])):
        a.add_patch(Rectangle((x, y - 0.09), 0.26, 0.18, fc=c, ec="white", lw=0.6))


def g_ils(a):
    a.plot([0.12, 0.5, 0.88], [0.18, 0.82, 0.18], color=C["ends"], lw=1.1)
    a.plot([0.5, 0.31], [0.18, 0.5], color=C["ends"], lw=1.1)
    for xx in (0.12, 0.5, 0.88):
        a.plot(xx, 0.18, "o", ms=2.3, color=C["ink2"])


def g_ws(a):
    a.add_patch(FancyBboxPatch((0.02, 0.26), 0.96, 0.48, boxstyle="round,pad=0,rounding_size=0.14", fc=C["pink"], ec="none"))
    a.text(0.5, 0.5, "W→S", ha="center", va="center", fontsize=6.3, color=C["f2b"], fontweight="bold")


def g_genomes(a):
    for y, c in ((0.8, C["human"]), (0.6, C["pan"]), (0.4, C["gorilla"]), (0.2, C["orang"])):
        rod(a, 0.06, y, 0.88, col=c, lw=2.6)


def g_map(a):
    xx = np.linspace(0.05, 0.95, 80)
    yy = 0.22 + 0.55 * np.exp(-((xx - 0.08) / 0.1) ** 2) + 0.5 * np.exp(-((xx - 0.93) / 0.09) ** 2) + 0.04 * np.sin(45 * xx)
    a.plot(xx, yy, color=C["human"], lw=1.0)
    a.plot([0.05, 0.95], [0.14, 0.14], color=C["ink2"], lw=0.6)


def g_population(a):
    for i in range(5):
        for j in range(4):
            a.add_patch(Circle((0.12 + 0.19 * i, 0.2 + 0.2 * j), 0.06, fc=C["human"] if (i * 3 + j) % 4 else C["ends"], ec="none"))


def g_archaic(a):
    a.plot([0.08, 0.92], [0.5, 0.5], color=C["ink2"], lw=5, solid_capstyle="round", ls=(0, (1.2, 0.8)), alpha=0.6)
    rod(a, 0.14, 0.5, 0.26, col=C["ends"], lw=3.2); rod(a, 0.52, 0.5, 0.32, col=C["ends"], lw=3.2)


def g_model(a):
    xx = np.linspace(0.08, 0.92, 50)
    a.plot(xx, 0.18 + 0.64 * np.exp(-(xx - 0.08) / 0.28), color=C["ink"], lw=1.2)
    rng = np.random.default_rng(3)
    for x0 in np.linspace(0.12, 0.86, 9):
        a.plot(x0, 0.18 + 0.64 * np.exp(-(x0 - 0.08) / 0.28) + rng.normal(0, 0.07), "o", ms=2.0, color=C["ends"])


def g_target(a):
    for r, c in ((0.42, C["band"]), (0.28, "white"), (0.14, C["band"])):
        a.add_patch(Circle((0.5, 0.5), r, fc=c, ec=C["ink2"], lw=0.5))
    a.plot([0.5, 0.5], [0.04, 0.96], color=C["ink"], lw=0.9)


def g_dice(a):
    a.add_patch(FancyBboxPatch((0.14, 0.14), 0.72, 0.72, boxstyle="round,pad=0,rounding_size=0.12", fc="white",
                               ec=C["ink2"], lw=0.8))
    for x0, y0 in ((0.33, 0.67), (0.5, 0.5), (0.67, 0.33)):
        a.add_patch(Circle((x0, y0), 0.06, fc=C["ink"], ec="none"))


def arrow(fig, p0, p1, head=True, col=None, lw=0.8):
    fig.patches.append(FancyArrowPatch(p0, p1, transform=fig.transFigure, arrowstyle="-|>" if head else "-",
                                       mutation_scale=7, color=col or C["ink2"], lw=lw))


def letter(fig, x, y, s):
    fig.text(x, y, s, fontsize=9, fontweight="bold", ha="left", va="top")


fig = plt.figure(figsize=(FIG_W, FIG_H))

# ---------------------------------------------------------------- a: three events on one timeline
ax = fig.add_axes([0.04, 0.60, 0.42, 0.32])
t = np.linspace(0, 7, 700)
T_orig, T_fix, T_off = 3.3, 2.3, 2.8
x = np.clip((T_orig - t) / (T_orig - T_fix), 0, 1)
x = 0.5 - 0.5 * np.cos(np.pi * x)
for lo, hi in ((5.5, 6.3), (0.55, 0.75)):
    ax.axvspan(lo, hi, color=C["band"], lw=0, zorder=0)
ax.plot(t, 1 - x, color=C["human"], lw=1.8)
ax.plot(t, x, color=C["f2b"], lw=1.8)
for tt in (T_orig, T_off, T_fix):
    ax.plot([tt, tt], [-0.04, 1.04], color=C["ink2"], lw=0.6, ls=(0, (2, 2)))
ax.text(T_orig + 0.08, -0.1, "$T_{orig}$", ha="right", va="top", fontsize=6.6)
ax.text(T_off, -0.1, "$T_{off}$", ha="center", va="top", fontsize=6.6)
ax.text(T_fix - 0.08, -0.1, "$T_{fix}$", ha="left", va="top", fontsize=6.6)
ax.plot([3.8, 1.5], [1.22, 1.22], color=C["f2b"], lw=3, solid_capstyle="round")
ax.text(2.65, 1.29, "this study (95%)", ha="center", va="bottom", fontsize=6, color=C["f2b"])
ax.annotate("", xy=(7.0, 1.22), xytext=(5.0, 1.22), arrowprops=dict(arrowstyle="-|>", color=C["gorilla"], lw=1.2))
ax.text(6.0, 1.29, "proposed >5 Mya", ha="center", va="bottom", fontsize=6, color=C["gorilla"])
ax.text(5.9, -0.05, "human–chimp\nsplit", ha="center", va="top", fontsize=5.6, color=C["ink2"], linespacing=1.0)
ax.text(0.65, 0.5, "modern–\narchaic\nsplit", ha="center", va="center", fontsize=5.6, color=C["ink2"], linespacing=1.0,
        bbox=dict(boxstyle="square,pad=0.15", fc="white", ec="none"))
ax.text(4.6, 0.93, "end-type recombination", ha="center", va="top", fontsize=6, color=C["human"])
ax.text(1.4, 0.93, "fused chromosome", ha="center", va="top", fontsize=6, color=C["f2b"])
ax.set_xlim(7, 0); ax.set_ylim(-0.32, 1.42)
ax.set_yticks([]); ax.spines["left"].set_visible(False)
ax.set_xticks(range(0, 8)); ax.set_xlabel("Mya", labelpad=1)
ax.spines["bottom"].set_position(("data", -0.32))
bb = ax.get_position()
fx = lambda v: bb.x0 + bb.width * (7 - v) / 7
fy = lambda v: bb.y0 + bb.height * (v + 0.32) / 1.74
chrom_icon(fig, fx(5.05), fy(0.66), 0.085, "ancestral")         # under the teal plateau (older side)
chrom_icon(fig, fx(2.2), fy(0.66), 0.07, "fused")             # under the crimson plateau (recent side)
letter(fig, 0.0, 0.985, "a")

# ---------------------------------------------------------------- b: what each line of evidence dates
ax = fig.add_axes([0.62, 0.60, 0.36, 0.32])
rows = [("sd", 5.0, None, C["ends"], "arrow"), ("ils", 8.6, 6.0, C["ends"], "bar"),
        ("ws", 3.8, 1.5, C["f2b"], "bar"), ("relic", 0.65, None, C["human"], "arrow")]
ys = [3, 2, 1, 0]
blab = {"sd": "age of junction duplications", "ils": "lineage sorting at the flanks",
        "ws": "W→S footprint: switch-off (this study)", "relic": "relic centromere in archaic genomes"}
for (kind, a0, a1, col, st), y in zip(rows, ys):
    ax.text(9.2, y + 0.22, blab[kind], ha="left", va="bottom", fontsize=6, color=col if kind in ("ws", "relic") else C["ink2"])
    if st == "arrow":
        ax.annotate("", xy=(9.2, y), xytext=(a0, y), arrowprops=dict(arrowstyle="-|>", color=col, lw=2.2))
    else:
        ax.plot([a1, a0], [y, y], color=col, lw=4.5, solid_capstyle="round")
ax.axvspan(5.5, 6.3, color=C["band"], lw=0, zorder=0)
ax.set_xlim(9.3, 0); ax.set_ylim(-0.6, 3.75)
ax.set_yticks([]); ax.spines["left"].set_visible(False)
ax.set_xticks([0, 2, 4, 6, 8]); ax.set_xlabel("Mya", labelpad=1)
ax.text(5.9, -0.55, "split", ha="center", va="bottom", fontsize=5.6, color=C["ink2"])
bb = ax.get_position()
for (kind, *_), y in zip(rows, ys):
    yc = bb.y0 + bb.height * (y + 0.6) / 4.35
    if kind == "relic":
        chrom_icon(fig, bb.x0 - 0.085, yc - RH / 2, 0.07, "relic")
    else:
        {"sd": g_sd, "ils": g_ils, "ws": g_ws}[kind](icon(fig, [bb.x0 - 0.085, yc - 0.03, 0.07, 0.06]))
letter(fig, 0.515, 0.985, "b")

# ---------------------------------------------------------------- c: the clock
chrom_icon(fig, 0.06, 0.32, 0.16, "ancestral")
fig.text(0.14, 0.355, "ancestral 2A + 2B", ha="center", va="bottom", fontsize=6, color=C["ink2"])
fig.text(0.14, 0.195, "human chromosome 2", ha="center", va="top", fontsize=6, color=C["ink2"])
arrow(fig, (0.14, 0.305), (0.14, 0.255))
chrom_icon(fig, 0.06, 0.22, 0.16, "fused")
ax = fig.add_axes([0.29, 0.08, 0.18, 0.34])
ts = np.linspace(0, 6, 200)
ax.plot(ts, ts * 0.1, color=C["ends"], lw=1.5)
ax.plot(ts, np.minimum(ts, 6 - T_off) * 0.1, color=C["f2b"], lw=1.8)
ax.plot(ts, ts * 0, color=C["interior"], lw=1.3, ls=(0, (3, 2)))
ax.axvline(6 - T_off, color=C["ink2"], lw=0.6, ls=(0, (2, 2)))
ax.annotate("", xy=(6.3, 0.6), xytext=(6.3, 0.32), arrowprops=dict(arrowstyle="<->", color=C["ink"], lw=0.8))
ax.text(6.5, 0.46, "$T_{off}/T_s$", fontsize=6.4, va="center")
ax.set_xticks([0, 6 - T_off, 6], ["$T_s$", "$T_{off}$", "0"]); ax.set_yticks([])
ax.set_xlim(0, 8.2); ax.set_ylim(-0.08, 0.72)
ax.spines["bottom"].set_bounds(0, 6)
ax.set_ylabel("W→S excess", fontsize=6.6)
bb = ax.get_position()
xr = lambda v: bb.x0 + bb.width * v / 8.2
yr = lambda v: bb.y0 + bb.height * (v + 0.08) / 0.8
# line-end glyphs: real end (capped), fusion flank (junction), interior (plain)
for yv, kind, dy, lab, col in ((0.6, "end", 0.008, "real ends", C["ends"]), (0.32, "fused", -0.03, "fusion flanks", C["f2b"]),
                               (0.0, "plain", 0.006, "interior", C["ink2"])):
    chrom_icon(fig, xr(6.1), yr(yv) + dy, 0.045, kind, junction=0.08)
    fig.text(xr(6.1) + 0.05, yr(yv) + dy + RH / 2, lab, ha="left", va="center", fontsize=6, color=col)
letter(fig, 0.0, 0.47, "c")

# ---------------------------------------------------------------- d: data and approach
ax = fig.add_axes([0.54, 0.06, 0.13, 0.37]); ax.axis("off"); ax.set_xlim(0, 3.2); ax.set_ylim(1.9, 5.3)
tips = {"H": (5.0, C["human"]), "C": (4.35, C["pan"]), "B": (3.75, C["pan"]), "G": (3.05, C["gorilla"]),
        "O": (2.4, C["orang"])}
xt = 2.5
anc4, anc5, anc3, root = (1.45, 4.4), (1.95, 4.05), (0.95, 3.7), (0.45, 3.05)
seg = [(anc4, (xt, tips["H"][0]), C["human"], 2.0), (anc4, anc5, C["pan"], 2.0), (anc5, (xt, tips["C"][0]), C["pan"], 2.0),
       (anc5, (xt, tips["B"][0]), C["pan"], 2.0), (anc3, anc4, C["ends"], 1.1), (anc3, (xt, tips["G"][0]), C["gorilla"], 2.0),
       (root, anc3, C["ends"], 1.1), (root, (xt, tips["O"][0]), C["orang"], 1.1)]
for (x0, y0), (x1, y1), col, lw in seg:
    ax.plot([x0, x0, x1], [y0, y1, y1], color=col, lw=lw, solid_capstyle="round")
for k, (y, col) in tips.items():
    ax.text(xt + 0.1, y, k, va="center", fontsize=6.8, color=col, fontweight="bold")
for p in (anc4, anc5, anc3):
    ax.plot(*p, "o", ms=3.2, color="white", mec=C["ink"], mew=0.8, zorder=3)
ix, iw = 0.795, 0.06
ys_d = (0.35, 0.275, 0.2, 0.125)
for fn, y, lab in zip((g_genomes, g_map, g_population, g_archaic), ys_d,
                     ("T2T ape genomes", "recombination maps", "1000 Genomes", "archaic genomes")):
    fn(icon(fig, [ix, y, iw, 0.055], equal=fn not in (g_genomes, g_map, g_archaic)))
    fig.text(ix - 0.006, y + 0.0275, lab, ha="right", va="center", fontsize=5.8, color=C["ink2"])
    arrow(fig, (ix + iw + 0.004, y + 0.0275), (0.874, 0.265), head=False, col=C["band"], lw=0.9)
g_model(icon(fig, [0.872, 0.23, 0.075, 0.07]))
arrow(fig, (0.9, 0.225), (0.878, 0.165))
arrow(fig, (0.92, 0.225), (0.952, 0.165))
g_target(icon(fig, [0.845, 0.105, 0.055, 0.055]))
g_dice(icon(fig, [0.93, 0.105, 0.055, 0.055]))
fig.text(0.91, 0.305, "pooled model", ha="center", va="bottom", fontsize=6, color=C["ink"])
fig.text(0.8725, 0.098, "calibration\non real ends", ha="center", va="top", fontsize=5.8, color=C["ink2"], linespacing=1.0)
fig.text(0.9575, 0.098, "simulations", ha="center", va="top", fontsize=5.8, color=C["ink2"])
letter(fig, 0.515, 0.47, "d")

fig.savefig(f"{OUT}/Fig1_overview.pdf", bbox_inches="tight", pad_inches=0.02)
fig.savefig(f"{OUT}/Fig1_overview.png", bbox_inches="tight", pad_inches=0.02, dpi=300)
print("ok")
