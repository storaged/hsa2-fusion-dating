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
FIG_W, FIG_H = 180 * MM, 100 * MM


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

# ---------------------------------------------------------------- a: the problem and the clock
ax = fig.add_axes([0.0, 0.535, 0.30, 0.325]); ax.axis("off"); ax.set_xlim(0, 1); ax.set_ylim(1.06, -0.02)
yt, yb = 0.0, 0.97
# ancestral 2A and 2B (apes); telomeres at all ends, own centromeres
ideogram(ax, 0.13, yt, yt + (yb - yt) * FUS / CHRLEN, 0, FUS, 0.07, cen=CEN2A)
ideogram(ax, 0.27, yt + 0.06, yt + 0.06 + (yb - yt) * (CHRLEN - FUS) / CHRLEN, FUS, CHRLEN, 0.07, cen=CEN2B)
ax.text(0.13, yt + (yb - yt) * FUS / CHRLEN + 0.03, "2A", ha="center", va="top", fontsize=6.3)
ax.text(0.27, yt + 0.06 + (yb - yt) * (CHRLEN - FUS) / CHRLEN + 0.03, "2B", ha="center", va="top", fontsize=6.3)
ax.annotate("", xy=(0.56, 0.5), xytext=(0.38, 0.5), arrowprops=dict(arrowstyle="-|>", color=C["ink2"], lw=0.9))
ax.text(0.47, 0.46, "fusion", ha="center", va="bottom", fontsize=6, color=C["ink2"])
# human chromosome 2
ideogram(ax, 0.72, yt, yb, 0, CHRLEN, 0.07, junction=FUS, relic=CEN2B, cen=CEN2A)
ax.text(0.72, yb + 0.03, "human chr2", ha="center", va="top", fontsize=6.3)
ax.text(0.80, yt + (yb - yt) * FUS / CHRLEN - 0.005, "2q13 junction", ha="left", va="bottom", fontsize=5.8, color=C["f2b"])
ax.text(0.80, yt + (yb - yt) * CEN2B / CHRLEN + 0.01, "relic\ncentromere", ha="left", va="top", fontsize=5.8, color=C["ink2"], linespacing=1.0)
# silhouettes: apes above the ancestral pair, humans above chr2
silhouette(fig, "chimp", [0.022, 0.875, 0.04, 0.075]); silhouette(fig, "gorilla", [0.064, 0.875, 0.05, 0.075])
silhouette(fig, "human", [0.196, 0.87, 0.02, 0.085]); silhouette(fig, "neanderthal", [0.218, 0.87, 0.02, 0.085])
# clock: W->S excess along the human branch
axc = fig.add_axes([0.33, 0.57, 0.13, 0.33])
ts = np.linspace(0, 6, 200)
T_off = 2.8
axc.plot(ts, ts * 0.1, color=C["ends"], lw=1.5)
axc.plot(ts, np.minimum(ts, 6 - T_off) * 0.1, color=C["f2b"], lw=1.8)
axc.plot(ts, ts * 0, color=C["interior"], lw=1.3, ls=(0, (3, 2)))
axc.axvline(6 - T_off, color=C["ink2"], lw=0.6, ls=(0, (2, 2)))
axc.annotate("", xy=(6.25, 0.6), xytext=(6.25, 0.32), arrowprops=dict(arrowstyle="<->", color=C["ink"], lw=0.8))
axc.text(6.4, 0.46, "lost:\n$T_{off}/T_s$", fontsize=6, va="center", linespacing=1.1)
axc.text(3.1, 0.47, "real ends", fontsize=6, color=C["ends"], ha="right")
axc.text(4.9, 0.27, "fusion\nflanks", fontsize=6, color=C["f2b"], ha="center", va="top", linespacing=1.0)
axc.text(4.9, 0.02, "interior", fontsize=6, color=C["ink2"], ha="center", va="bottom")
axc.set_xticks([0, 6 - T_off, 6], ["split", "$T_{off}$", "now"]); axc.set_yticks([])
axc.set_xlim(0, 8.3); axc.set_ylim(-0.05, 0.72); axc.spines["bottom"].set_bounds(0, 6)
axc.set_ylabel("W→S excess (human)", fontsize=6.3, labelpad=2)
letter(fig, 0.0, TOP, "a")

# ---------------------------------------------------------------- b: three events on one timeline
ax = fig.add_axes([0.56, 0.60, 0.43, 0.33])
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
        {"sd": g_sd, "ils": g_ils, "ws": g_ws}[kind](icon(fig, [bb.x0 - 0.085, yc - 0.03, 0.07, 0.06]))
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
# flow
xs = [0.715, 0.785, 0.855, 0.925]
iw, ih, yi = 0.052, 0.05, 0.345
for fn, xx, lab in zip((g_genomes, g_map, g_population, g_archaic), xs,
                       ("T2T apes", "recombination\nmaps", "1000\nGenomes", "archaic\ngenomes")):
    fn(icon(fig, [xx - iw / 2, yi, iw, ih], equal=fn not in (g_genomes, g_map, g_archaic)))
    fig.text(xx, yi + ih + 0.008, lab, ha="center", va="bottom", fontsize=5.6, color=C["ink2"], linespacing=1.0)
    arrow(fig, (xx, yi - 0.005), (0.82, 0.29), head=False, col=C["interior"], lw=0.8)
ym = 0.215
g_model(icon(fig, [0.79, ym, 0.06, 0.07]))
fig.text(0.855, ym + 0.035, "pooled model\n(39 real ends,\n3 lineages)", ha="left", va="center", fontsize=5.8, color=C["ink"], linespacing=1.05)
arrow(fig, (0.81, ym - 0.005), (0.775, 0.145)); arrow(fig, (0.83, ym - 0.005), (0.865, 0.145))
g_target(icon(fig, [0.75, 0.09, 0.05, 0.05])); g_dice(icon(fig, [0.84, 0.09, 0.05, 0.05]))
fig.text(0.775, 0.083, "calibration\non real ends", ha="center", va="top", fontsize=5.6, color=C["ink2"], linespacing=1.0)
fig.text(0.865, 0.083, "simulations", ha="center", va="top", fontsize=5.6, color=C["ink2"])
letter(fig, 0.515, 0.47, "d")

fig.savefig(f"{OUT}/Fig1_overview.pdf", bbox_inches="tight", pad_inches=0.02)
fig.savefig(f"{OUT}/Fig1_overview.png", bbox_inches="tight", pad_inches=0.02, dpi=300)
print("ok")
