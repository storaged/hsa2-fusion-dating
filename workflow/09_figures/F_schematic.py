"""Figure 1 (conceptual overview), Nature style: panel letters only, glyphs instead of labels; the caption explains.
Glyphs: capped rods = ancestral acrocentric chromosomes (black caps = telomeres); rod with a red band = fused
chromosome 2 (junction); rod with an open waist = relic centromere; stacked blocks = segmental duplications;
crossed tree = lineage sorting; chip = W->S substitutions; data/model/test icons in d.
Panel a: G-banded ideograms (GRCh38 cytobands, UCSC) of the ancestral 2A/2B and human chromosome 2, with public-domain
PhyloPic silhouettes (see assets/SOURCES.txt).
usage: python F_schematic.py results/figures [assets_dir=data/silhouettes]
"""
import os
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch, Rectangle

OUT = sys.argv[1] if len(sys.argv) > 1 else "results/figures"
ASSETS = sys.argv[2] if len(sys.argv) > 2 else "data/silhouettes"
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
FIG_W, FIG_H = 180 * MM, 88 * MM


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
    a.plot([x0 + e, x0 + w - e], [y, y], color=col or C["rod"], lw=lw, solid_capstyle="butt", clip_on=False)
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
    """three lineage kernels sharing one end strength; fusion-flank points below (retained fraction)"""
    xx = np.linspace(0.06, 0.94, 60)
    for col, amp in ((C["pan"], 0.72), (C["gorilla"], 0.66), (C["human"], 0.6)):
        a.plot(xx, 0.12 + amp * np.exp(-(xx - 0.06) / 0.3), color=col, lw=1.1)
    for x0 in (0.12, 0.22, 0.34, 0.48):
        a.plot(x0, 0.12 + 0.32 * np.exp(-(x0 - 0.06) / 0.3), "s", ms=2.3, color=C["f2b"], mec="none")
    a.plot([0.04, 0.96], [0.08, 0.08], color=C["ink2"], lw=0.5)


def g_target(a):
    """calibration: intervals of real ends treated as fusions, crossing the true value 0"""
    rng = np.random.default_rng(5)
    for k, y in enumerate(np.linspace(0.12, 0.88, 8)):
        c = rng.normal(0.5, 0.12)
        col = C["ends"] if k != 6 else C["ink"]
        a.plot([c - 0.22, c + 0.22] if k != 6 else [0.62, 0.92], [y, y], color=col, lw=1.0, solid_capstyle="butt")
    a.plot([0.5, 0.5], [0.04, 0.96], color=C["ink"], lw=0.9)


def g_dice(a):
    """simulations: estimates on the identity line"""
    a.plot([0.08, 0.92], [0.08, 0.92], color=C["ink2"], lw=0.6, ls=(0, (2, 1.5)))
    for v in (0.15, 0.33, 0.5, 0.67, 0.85):
        a.plot([v, v], [v - 0.08, v + 0.08], color=C["human"], lw=0.9)
        a.plot(v, v + 0.01, "o", ms=2.4, color=C["human"], mec="none")
    a.plot([0.08, 0.08, 0.92], [0.92, 0.08, 0.08], color=C["ink2"], lw=0.6)


def arrow(fig, p0, p1, head=True, col=None, lw=0.8):
    fig.patches.append(FancyArrowPatch(p0, p1, transform=fig.transFigure, arrowstyle="-|>" if head else "-",
                                       mutation_scale=7, color=col or C["ink2"], lw=lw))


def letter(fig, x, y, s):
    fig.text(x, y, s, fontsize=9, fontweight="bold", ha="left", va="top")


# ---------------------------------------------------------------- ideograms and silhouettes
STAIN = {"gneg": "#f4f5f7", "gpos25": "#cfd2d7", "gpos50": "#a9adb4", "gpos75": "#7d828a", "gpos100": "#4a4f57",
         "gvar": "#b9bdc4", "stalk": "#dfe2e6", "acen": "#d9dce1"}
CB = [l.rstrip("\n").split("\t") for l in open(f"{ASSETS}/chr2_cytoband.tsv")]
CB = [(int(a) / 1e6, int(b) / 1e6, st) for _, a, b, _, st in CB]
FUS, CEN2A, CEN2B, CHRLEN = 113.57, 93.9, 132.0, 242.19


def ideogram(ax, x, y0, y1, start, end, width, telomeres=(True, True), junction=None, relic=None, cen=None):
    """vertical G-banded chromosome for hg38 chr2 interval [start, end] Mb drawn from y0 (top) to y1 (bottom)"""
    sc = (y1 - y0) / (end - start)
    Y = lambda mb: y0 + (mb - start) * sc
    for a, b, st in CB:
        a, b = max(a, start), min(b, end)
        if b <= a:
            continue
        ax.add_patch(Rectangle((x - width / 2, Y(a)), width, Y(b) - Y(a), fc=STAIN.get(st, "#f4f5f7"), ec="none"))
    ax.plot([x - width / 2, x - width / 2], [Y(start), Y(end)], color=C["ink2"], lw=0.5)
    ax.plot([x + width / 2, x + width / 2], [Y(start), Y(end)], color=C["ink2"], lw=0.5)
    if cen is not None:                                       # centromeric constriction: two white notches
        h = 3.5 * abs(sc)
        for sgn in (-1, 1):
            xe = x + sgn * (width / 2 + 0.002)
            ax.add_patch(mpl.patches.Polygon([(xe, Y(cen) - h), (x + sgn * width * 0.12, Y(cen)), (xe, Y(cen) + h)],
                                             closed=True, fc="white", ec=C["ink2"], lw=0.5, zorder=3))
    cap = 0.012 * abs(y1 - y0) / 0.6
    if telomeres[0]:
        ax.add_patch(Rectangle((x - width / 2, Y(start) - np.sign(sc) * cap), width, np.sign(sc) * cap, fc=C["ink"], ec="none"))
    if telomeres[1]:
        ax.add_patch(Rectangle((x - width / 2, Y(end)), width, np.sign(sc) * cap, fc=C["ink"], ec="none"))
    if junction is not None:
        ax.add_patch(Rectangle((x - width * 0.75, Y(junction) - 0.6 * cap), width * 1.5, 1.2 * cap, fc=C["f2b"], ec="none", zorder=4))
    if relic is not None:
        ax.plot(x + width * 0.95, Y(relic), marker="<", ms=3.2, color=C["ink2"], mec="none")


def ideogram_h(ax, y, x0, x1, start, end, height, telomeres=(True, True), junction=None, relic=None, cen=None):
    """horizontal G-banded chromosome for hg38 chr2 [start, end] Mb, drawn from x0 (left) to x1 (right)"""
    sc = (x1 - x0) / (end - start)
    X = lambda mb: x0 + (mb - start) * sc
    for a, b, st in CB:
        a, b = max(a, start), min(b, end)
        if b > a:
            ax.add_patch(Rectangle((X(a), y - height / 2), X(b) - X(a), height, fc=STAIN.get(st, "#f4f5f7"), ec="none"))
    ax.add_patch(Rectangle((x0, y - height / 2), x1 - x0, height, fc="none", ec=C["ink2"], lw=0.5))
    cap = 0.012
    if telomeres[0]:
        ax.add_patch(Rectangle((x0 - cap, y - height / 2), cap, height, fc=C["ink"], ec="none"))
    if telomeres[1]:
        ax.add_patch(Rectangle((x1, y - height / 2), cap, height, fc=C["ink"], ec="none"))
    if cen is not None:
        w = 3.5 * sc
        for sgn in (-1, 1):
            ye = y + sgn * (height / 2 + 0.004)
            ax.add_patch(mpl.patches.Polygon([(X(cen) - w, ye), (X(cen), y + sgn * height * 0.12), (X(cen) + w, ye)],
                                             closed=True, fc="white", ec=C["ink2"], lw=0.5, zorder=3))
    if junction is not None:
        ax.add_patch(Rectangle((X(junction) - 0.006, y - height * 0.8), 0.012, height * 1.6, fc=C["f2b"], ec="none", zorder=4))
    if relic is not None:
        ax.plot(X(relic), y - height * 0.95, marker="^", ms=3.2, color=C["ink2"], mec="none", clip_on=False)
    return X


def silhouette(fig, name, rect, color="#8a8f98"):
    img = plt.imread(f"{ASSETS}/{name}.png")
    rgba = np.zeros_like(img)
    rgb = mpl.colors.to_rgb(color)
    rgba[..., 0], rgba[..., 1], rgba[..., 2] = rgb
    rgba[..., 3] = img[..., 3] if img.shape[-1] == 4 else 1 - img[..., 0]
    a = fig.add_axes(rect); a.imshow(rgba); a.set_aspect("equal"); a.axis("off")
    return a


fig = plt.figure(figsize=(FIG_W, FIG_H))
TOP, MID = 0.975, 0.5

# ---------------------------------------------------------------- a: the fusion (top to bottom) and the clock
ax = fig.add_axes([0.02, 0.50, 0.27, 0.46]); ax.axis("off"); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
H = 0.075                                     # ideogram thickness (axes fraction)
L0, L1 = 0.03, 0.97                           # chr2 extent; 2A + gap + 2B spans the same
gap = 0.08
sc = (L1 - L0 - gap) / CHRLEN
a1 = L0 + FUS * sc
y_anc, y_hum = 0.66, 0.30
X2a = ideogram_h(ax, y_anc, L0, a1, 0, FUS, H, cen=CEN2A)
X2b = ideogram_h(ax, y_anc, a1 + gap, L1, FUS, CHRLEN, H, cen=CEN2B)
ax.text(L0, y_anc - H / 2 - 0.035, "2A", ha="left", va="top", fontsize=6.3)
ax.text(L1, y_anc - H / 2 - 0.035, "2B", ha="right", va="top", fontsize=6.3)
# funnel: the two chromosomes narrow onto the junction of chromosome 2
xj = L0 + FUS * (L1 - L0) / CHRLEN
ax.add_patch(mpl.patches.Polygon([(L0, y_anc - H / 2 - 0.1), (L1, y_anc - H / 2 - 0.1), (xj + 0.035, y_hum + H / 2 + 0.085),
                                  (xj - 0.035, y_hum + H / 2 + 0.085)], closed=True, fc=C["band"], ec="none", zorder=0))
ax.text(xj, (y_anc - H / 2 - 0.1 + y_hum + H / 2 + 0.085) / 2 + 0.01, "fusion", ha="center", va="center", fontsize=6, color=C["ink2"])
X2 = ideogram_h(ax, y_hum, L0, L1, 0, CHRLEN, H, junction=FUS, relic=CEN2B, cen=CEN2A)
ax.text(X2(FUS), y_hum + H / 2 + 0.015, "2q13", ha="center", va="bottom", fontsize=5.8, color=C["f2b"])
ax.text(X2(CEN2B) + 0.02, y_hum - H / 2 - 0.05, "relic centromere", ha="left", va="top", fontsize=5.8, color=C["ink2"])
ax.text(L0, y_hum - H / 2 - 0.05, "human chr2", ha="left", va="top", fontsize=6.3)
bb = ax.get_position()
fxa = lambda v: bb.x0 + bb.width * v
fya = lambda v: bb.y0 + bb.height * v
silhouette(fig, "chimp", [fxa(0.2), fya(0.77), 0.036, 0.075]); silhouette(fig, "bonobo", [fxa(0.39), fya(0.775), 0.052, 0.07])
silhouette(fig, "gorilla", [fxa(0.62), fya(0.77), 0.046, 0.075])
silhouette(fig, "human", [fxa(0.37), fya(-0.02), 0.018, 0.09]); silhouette(fig, "neanderthal", [fxa(0.47), fya(-0.02), 0.018, 0.09])
# clock: square panel on the right
axc = fig.add_axes([0.345, 0.585, 0.145, 0.30])
ts = np.linspace(0, 6, 200)
T_off = 2.8
axc.plot(ts, ts * 0.1, color=C["ends"], lw=1.5)
axc.plot(ts, np.minimum(ts, 6 - T_off) * 0.1, color=C["f2b"], lw=1.8)
axc.plot(ts, ts * 0, color=C["interior"], lw=1.3, ls=(0, (3, 2)))
axc.axvline(6 - T_off, color=C["ink2"], lw=0.6, ls=(0, (2, 2)))
axc.annotate("", xy=(6.25, 0.6), xytext=(6.25, 0.32), arrowprops=dict(arrowstyle="<->", color=C["ink"], lw=0.8))
axc.text(6.4, 0.46, "lost:\n$T_{off}/T_s$", fontsize=6, va="center", linespacing=1.1)
axc.text(3.0, 0.45, "real ends", fontsize=6, color=C["ends"], ha="right")
axc.text(4.6, 0.26, "fusion\nflanks", fontsize=6, color=C["f2b"], ha="center", va="top", linespacing=1.0)
axc.text(4.6, 0.02, "interior", fontsize=6, color=C["ink2"], ha="center", va="bottom")
axc.set_xticks([0, 6 - T_off, 6], ["split", "$T_{off}$", "now"]); axc.set_yticks([])
axc.set_xlim(0, 8.3); axc.set_ylim(-0.05, 0.72); axc.spines["bottom"].set_bounds(0, 6)
axc.set_ylabel("W→S excess (human)", fontsize=6.3, labelpad=2)
letter(fig, 0.0, TOP, "a")

# ---------------------------------------------------------------- b: three events on one timeline
ax = fig.add_axes([0.56, 0.615, 0.43, 0.33])
t = np.linspace(0, 7, 700)
T_orig, T_fix = 3.3, 2.3
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
letter(fig, 0.515, TOP, "b")

# ---------------------------------------------------------------- c: what each line of evidence dates
ax = fig.add_axes([0.12, 0.08, 0.36, 0.34])
rows = [("sd", 5.0, None, C["ends"], "arrow"), ("ils", 8.6, 6.0, C["ends"], "bar"),
        ("ws", 3.8, 1.5, C["f2b"], "bar"), ("relic", 0.65, None, C["human"], "arrow")]
blab = {"sd": "age of junction duplications", "ils": "lineage sorting at the flanks",
        "ws": "W→S footprint: switch-off (this study)", "relic": "relic centromere in archaic genomes"}
ys = [3, 2, 1, 0]
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
        {"sd": g_sd, "ils": g_ils, "ws": g_ws}[kind](icon(fig, [bb.x0 - 0.09, yc - 0.04, 0.08, 0.08]))
letter(fig, 0.0, 0.47, "c")

# ---------------------------------------------------------------- d: data and approach, top-to-bottom flow
ax = fig.add_axes([0.53, 0.05, 0.13, 0.38]); ax.axis("off"); ax.set_xlim(0, 3.2); ax.set_ylim(1.9, 5.3)
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
# flow with orthogonal connectors
def seg(p0, p1, head=False, col=C["ink2"]):
    fig.patches.append(FancyArrowPatch(p0, p1, transform=fig.transFigure, arrowstyle="-|>" if head else "-",
                                       mutation_scale=7, color=col, lw=0.8, shrinkA=0, shrinkB=0))


xs = [0.715, 0.785, 0.855, 0.925]
iw, ih, yi = 0.052, 0.045, 0.355
xc = 0.82
for fn, xx, lab in zip((g_genomes, g_map, g_population, g_archaic), xs,
                       ("T2T apes", "recombination\nmaps", "1000\nGenomes", "archaic\ngenomes")):
    fn(icon(fig, [xx - iw / 2, yi, iw, ih], equal=fn not in (g_genomes, g_map, g_archaic)))
    fig.text(xx, yi + ih + 0.008, lab, ha="center", va="bottom", fontsize=5.6, color=C["ink2"], linespacing=1.0)
    seg((xx, yi - 0.006), (xx, 0.325))
bus_in = 0.325
seg((xs[0], bus_in), (xs[-1], bus_in))
ym, mh = 0.225, 0.07
seg((xc, bus_in), (xc, ym + mh + 0.004), head=True)
g_model(icon(fig, [xc - 0.04, ym, 0.08, mh], equal=False))
fig.text(xc + 0.048, ym + mh / 2, "pooled model\n39 ends × 3 lineages", ha="left", va="center", fontsize=5.8, color=C["ink"], linespacing=1.1)
bus_out = 0.19
xcal, xsim = 0.765, 0.875
seg((xc, ym - 0.004), (xc, bus_out)); seg((xcal, bus_out), (xsim, bus_out))
yo, oh = 0.1, 0.055
seg((xcal, bus_out), (xcal, yo + oh + 0.004), head=True); seg((xsim, bus_out), (xsim, yo + oh + 0.004), head=True)
g_target(icon(fig, [xcal - 0.03, yo, 0.06, oh], equal=False)); g_dice(icon(fig, [xsim - 0.03, yo, 0.06, oh], equal=False))
fig.text(xcal, yo - 0.008, "calibration\n95% coverage", ha="center", va="top", fontsize=5.6, color=C["ink2"], linespacing=1.05)
fig.text(xsim, yo - 0.008, "simulations\n15/15 recovered", ha="center", va="top", fontsize=5.6, color=C["ink2"], linespacing=1.05)
letter(fig, 0.515, 0.47, "d")

fig.savefig(f"{OUT}/Fig1_overview.pdf", bbox_inches="tight", pad_inches=0.02)
fig.savefig(f"{OUT}/Fig1_overview.png", bbox_inches="tight", pad_inches=0.02, dpi=300)
print("ok")
