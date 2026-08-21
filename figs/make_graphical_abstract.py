"""Overview figure and table-of-contents graphic (drawn from primitives, no AI image
generation, so every mark is exact and reproducible).

The two are deliberately different objects. The overview figure is the article's Figure 1 and
carries the pipeline; the table-of-contents graphic carries the single comparison the article
turns on, at the size a contents entry is printed.

OVERVIEW, left to right: every alchemical edge of the benchmark arrives with a per-edge error bar
that the estimator computes rather than learns, composed from the two legs' uncertainties; that
bar is validated against three independent protocol repeats; and it is then the null of a
cycle-closure quality-control fit. On real data (1143 edges in the 48 systems with a cycle, of
1145 over 49 systems) the test flags a small, chemically sensible minority of networks (6, 6 and
3 of 48 on the three repeats); a stand-in bar built from a learned head's measured
overconfidence flags nearly everything (42 of 48). The contrast is carried by grid area, not by a
sentence.

What the overview deliberately does NOT draw: forward and reverse work samples entering the BAR
bottleneck to give B/I^2. That path is the closed form of the theory section and is exercised on
the alchemtest windows, where the works are released; the benchmark ships per-leg MBAR
uncertainties and no works, so drawing it as the benchmark pipeline would depict a computation
that was never run there.

The overview is authored at the article's full text width (paperstyle.FULL = 6.5 in) and in the
article's semantic palette, so LaTeX reproduces it at scale 1.0 and a colour means here what it
means in every other figure: OURS for this article's calibrated quantity and everything computed
from it (the per-edge bar, the replicate check that validates that same bar, and the flags it
raises), FOIL for the overconfident stand-in only, REF/MUTED for structure that is present but
not the point. The contents graphic keeps its own pastel ground, which is sized and printed
elsewhere.

Deterministic: no randomness, no network access. Run:
    python figs/make_graphical_abstract.py
"""
from __future__ import annotations

import pathlib
import sys

import matplotlib as mpl
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch, Rectangle

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from paperstyle import (  # noqa: E402
    FOIL,
    FULL,
    INK,
    MUTED,
    OURS,
    REF,
    finish,
    panel,
    tint,
    use_paper_style,
)

# ---- pastel palette of the table-of-contents graphic (that artwork only) ----
BG = "#FAF7F2"
BLUE = "#A8C5DA"
SAGE = "#B5CDB7"
CORAL = "#E0A9A0"
INK_TOC = "#3A4454"

ANNOT = 8.5


def _darken(hexcolor: str, factor: float = 0.7):
    r, g, b = mcolors.to_rgb(hexcolor)
    return (r * factor, g * factor, b * factor)


def _tint(hexcolor: str, amount: float = 0.78):
    r, g, b = mcolors.to_rgb(hexcolor)
    br, bgc, bb = mcolors.to_rgb(BG)
    return (r + (br - r) * amount, g + (bgc - g) * amount, b + (bb - b) * amount)


def _style():
    """rcParams of the contents graphic: its own pastel ground, and no tight bbox."""
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans", "Arial"],
        "mathtext.fontset": "dejavusans",
        "font.size": ANNOT,
        "savefig.dpi": 300,
        "savefig.facecolor": BG,
        "savefig.bbox": "standard",
        "savefig.pad_inches": 0.1,
        "figure.facecolor": BG,
        "axes.unicode_minus": True,
    })


def _canvas(w, h):
    fig, ax = plt.subplots(figsize=(w, h))
    ax.set_xlim(0, w)
    ax.set_ylim(0, h)
    ax.set_position([0, 0, 1, 1])
    ax.axis("off")
    return fig, ax


def _writer(ax, words):
    def txt(x, y, s, size, color=INK_TOC, weight="normal", style="normal",
            ha="center", va="center", zorder=6, count=True):
        ax.text(x, y, s, ha=ha, va=va, fontsize=size, color=color,
                 fontweight=weight, fontstyle=style, zorder=zorder)
        if count:
            words.extend(s.replace("/", " ").split())
    return txt


# ======================================================================================
# Figure 1 (fig:overview). Geometry is in INCHES on a full-text-width canvas: every axes
# is placed by _stage()/_ax_in() from inches, so one unit of data is one inch on the page
# and the panels share one baseline, one height and one card geometry by construction.
# ======================================================================================
OW = FULL

# The outcome row is laid out first, because everything above it is stacked on top of it.
# 48 systems as 12 x 4 rather than 8 x 6: the two waffles deliver two integers, and a wide
# flat block spends a third less page height on them while keeping the comparison instant
# (6/48 is half of the first row of four; 42/48 is three and a half of the four).
COLS, ROWS = 12, 4                       # 48 systems with a cycle
PITCH, CELL = 0.215, 0.190
GW = COLS * PITCH - (PITCH - CELL)
GH = ROWS * PITCH - (PITCH - CELL)
WAFFLE_GAP = 0.90
WAFFLE_X0 = (OW - (2 * GW + WAFFLE_GAP)) / 2
WAFFLE_Y0 = 0.02
WAFFLE_H = GH + 0.26                     # the grid, plus the line the fraction is set on

ROW1_BOT = WAFFLE_Y0 + WAFFLE_H + 0.68   # the three pipeline cards, above the outcome row
ROW1_TOP = ROW1_BOT + 1.68
CARD_R = 0.10                            # one corner radius for every card
FAN_Y = ROW1_BOT - 0.24                  # the connector that carries stage c into the row below
TITLE_Y = ROW1_TOP + 0.58
SUBTITLE_Y = ROW1_TOP + 0.32
OH = TITLE_Y + 0.18

CARD_FACE = tint(REF, 0.90)
CARD_EDGE = tint(REF, 0.55)
NODE_FACE = tint(REF, 0.72)
EMPTY = tint(REF, 0.82)

TITLE_PT, BODY_PT, SMALL_PT, FRAC_PT = 10.5, 8.5, 7.5, 10.0


def _ax_in(fig, x0, y0, w, h):
    """An axes placed in inches, whose data units are inches (so circles are circles)."""
    ax = fig.add_axes([x0 / OW, y0 / OH, w / OW, h / OH])
    ax.set_xlim(0, w)
    ax.set_ylim(0, h)
    ax.axis("off")
    return ax


def _card(ax, w, h, inset=0.03):
    ax.add_patch(FancyBboxPatch(
        (inset, inset), w - 2 * inset, h - 2 * inset,
        boxstyle=f"round,pad=0,rounding_size={CARD_R}",
        linewidth=0.8, edgecolor=CARD_EDGE, facecolor=CARD_FACE, zorder=1))


def _fig_arrow(fig, p0, p1, scale=11.0):
    fig.add_artist(FancyArrowPatch(
        (p0[0] / OW, p0[1] / OH), (p1[0] / OW, p1[1] / OH), transform=fig.transFigure,
        arrowstyle="-|>", mutation_scale=scale, linewidth=1.2, color=INK,
        shrinkA=0, shrinkB=0, zorder=2))


def _fig_line(fig, xs, ys):
    fig.add_artist(Line2D(
        [x / OW for x in xs], [y / OH for y in ys], transform=fig.transFigure,
        color=INK, linewidth=1.1, solid_capstyle="round", solid_joinstyle="round",
        zorder=2))


def _edge_glyph(ax, cx, cy, half=0.52, node_r=0.095, bar=0.24):
    """One alchemical edge: two ligand nodes, the estimate, and its error bar."""
    ax.plot([cx - half, cx + half], [cy, cy], color=MUTED, linewidth=1.6,
            zorder=3, solid_capstyle="round")
    for sgn in (-1, +1):
        ax.add_patch(Circle((cx + sgn * half, cy), node_r, facecolor=NODE_FACE,
                            edgecolor=REF, linewidth=0.9, zorder=5))
    ax.errorbar([cx], [cy], yerr=[bar], fmt="o", ms=4.2, color=OURS, ecolor=OURS,
                elinewidth=1.4, capsize=4.0, capthick=1.4, zorder=6)


def _repeat_glyph(ax, cx, cy, spread=0.32, bar=0.22):
    """The three independent protocol repeats the per-edge bar is checked against.

    Drawn in OURS: these are the calibrated per-edge bar itself, measured against the
    repeats, which is the quantity Figure 4 plots on the real data.
    """
    for dx in (-spread, 0.0, +spread):
        ax.errorbar([cx + dx], [cy], yerr=[bar], fmt="o", ms=3.6, color=OURS,
                    ecolor=OURS, elinewidth=1.3, capsize=3.4, capthick=1.3, zorder=5)


def _network(ax):
    """Small perturbation network: one cycle highlighted, one edge flagged."""
    a = (0.64, 1.22)
    b = (1.22, 1.25)
    c = (1.26, 0.73)
    d = (0.68, 0.70)
    e = (1.70, 0.99)
    f = (0.20, 0.95)

    for p, q in ((b, e), (d, f)):
        ax.plot([p[0], q[0]], [p[1], q[1]], color=MUTED, linewidth=0.9,
                zorder=3, solid_capstyle="round")
    for p, q in ((a, b), (b, c), (d, a)):
        ax.plot([p[0], q[0]], [p[1], q[1]], color=INK, linewidth=1.8,
                zorder=4, solid_capstyle="round")

    ax.plot([d[0], c[0]], [d[1], c[1]], color=OURS, linewidth=2.4,
            zorder=4, solid_capstyle="round")

    for node in (a, b, c, d, e, f):
        ax.add_patch(Circle(node, 0.085, facecolor=NODE_FACE, edgecolor=REF,
                            linewidth=0.9, zorder=6))

    mx, my = (c[0] + d[0]) / 2, (c[1] + d[1]) / 2
    ax.plot([mx, mx], [my - 0.05, 0.42], color=OURS, linewidth=1.0, zorder=5)
    return mx, 0.39


def _waffle(ax, n_flag, total, flag_color):
    """COLS x ROWS unit squares, filled in reading order, top-left first.

    Both waffles fill the same way, so the flagged block of one is directly comparable
    with the flagged block of the other: the contrast is area, not a sentence.
    """
    top = GH + 0.24
    for r in range(ROWS):
        for col in range(COLS):
            idx = r * COLS + col
            ax.add_patch(Rectangle(
                (col * PITCH, top - r * PITCH - CELL), CELL, CELL,
                facecolor=flag_color if idx < n_flag else EMPTY,
                edgecolor="none", zorder=5))
    ax.text(GW / 2, 0.14, f"{n_flag}/{total}", ha="center", va="top",
            fontsize=FRAC_PT, fontweight="bold", color=INK, zorder=6)


def overview():
    """Figure 1: the pipeline the article actually runs on the benchmark."""
    fig = plt.figure(figsize=(OW, OH))
    words: list[str] = []

    def head(ax, letter, title):
        panel(ax, letter, title)
        words.extend(title.split())

    def say(x, y, s, size, *, ax=None, count=True, **kw):
        kw.setdefault("ha", "center")
        kw.setdefault("va", "center")
        kw.setdefault("color", INK)
        if ax is None:
            fig.text(x / OW, y / OH, s, fontsize=size, **kw)
        else:
            ax.text(x, y, s, fontsize=size, zorder=8, **kw)
        if count:
            words.extend(s.replace("/", " ").split())

    # ---- title ----
    say(OW / 2, TITLE_Y, "A closed-form error bar turns cycle closure into a quality-control test",
        TITLE_PT, fontweight="bold")
    say(OW / 2, SUBTITLE_Y, "1143 real binding edges, 48 systems with a cycle",
        SMALL_PT, fontstyle="italic", color=REF)

    # ==== the three stages: one baseline, one height, one card geometry ====
    h = ROW1_TOP - ROW1_BOT
    xa, wa = 0.05, 1.90
    xb, wb = 2.35, 1.80
    xc, wc = 4.55, 1.90

    # stage a: the per-edge bar the benchmark ships, computed not learned
    ax_a = _ax_in(fig, xa, ROW1_BOT, wa, h)
    _card(ax_a, wa, h)
    head(ax_a, "A", "per-edge error bar")
    _edge_glyph(ax_a, wa / 2, 1.18)
    say(wa / 2, 0.68,
        r"$\mathrm{se}_e=\sqrt{\mathrm{se}_{\mathrm{complex}}^{2}"
        r"+\mathrm{se}_{\mathrm{solvent}}^{2}}$", BODY_PT, ax=ax_a, count=False)
    say(wa / 2, 0.34, "computed, not learned", SMALL_PT, ax=ax_a, color=OURS)

    _fig_arrow(fig, (xa + wa + 0.05, (ROW1_TOP + ROW1_BOT) / 2),
               (xb - 0.05, (ROW1_TOP + ROW1_BOT) / 2))

    # stage b: what makes that bar trustworthy
    ax_b = _ax_in(fig, xb, ROW1_BOT, wb, h)
    _card(ax_b, wb, h)
    head(ax_b, "B", "validated on replicates")
    _repeat_glyph(ax_b, wb / 2, 1.18)
    say(wb / 2, 0.68, r"$\mathrm{reported/replicate}=1.41$", BODY_PT, ax=ax_b, count=False)

    _fig_arrow(fig, (xb + wb + 0.05, (ROW1_TOP + ROW1_BOT) / 2),
               (xc - 0.05, (ROW1_TOP + ROW1_BOT) / 2))

    # stage c: the outcome -- a perturbation network under the cycle-closure fit
    ax_c = _ax_in(fig, xc, ROW1_BOT, wc, h)
    _card(ax_c, wc, h)
    head(ax_c, "C", "cycle-closure fit")
    lx, ly = _network(ax_c)
    say(lx, ly, "flagged", SMALL_PT, ax=ax_c, color=OURS, va="top")

    # ==== the outcome row, wired to stage c so it reads as the same pipeline ====
    cx_left = WAFFLE_X0 + GW / 2
    cx_right = WAFFLE_X0 + GW + WAFFLE_GAP + GW / 2
    _fig_line(fig, [cx_left, xc + wc / 2, xc + wc / 2], [FAN_Y, FAN_Y, ROW1_BOT])
    for cx in (cx_left, cx_right):
        _fig_arrow(fig, (cx, FAN_Y), (cx, FAN_Y - 0.20), scale=10.0)

    ax_d = _ax_in(fig, WAFFLE_X0, WAFFLE_Y0, GW, WAFFLE_H)
    head(ax_d, "D", "calibrated")
    _waffle(ax_d, 6, 48, OURS)

    ax_e = _ax_in(fig, WAFFLE_X0 + GW + WAFFLE_GAP, WAFFLE_Y0, GW, WAFFLE_H)
    head(ax_e, "E", "overconfident stand-in")
    _waffle(ax_e, 42, 48, FOIL)
    words.extend(["6", "48", "42", "48"])

    n_words = len(words)
    print(f"word count on overview figure: {n_words}")
    assert n_words <= 40, f"too many words on figure: {n_words}"
    # the axes are placed in inches by _ax_in; a layout engine would move them
    finish(fig, "graphical_abstract", layout=None)
    plt.close(fig)


def toc():
    """Table-of-contents graphic: one object carrying the article's actual claim.

    A contents entry travels without its caption, so it must not carry a comparison the body
    scopes. It previously showed the flag counts under a calibrated bar against an overconfident
    stand-in, which Section 4 now scopes to heads five to eleven times too small. What it shows
    instead is the split the paper is about and how little of a real network's error falls on the
    side a cycle can see.
    """
    w, h = 3.3, 1.85
    fig, ax = _canvas(w, h)
    words: list[str] = []
    txt = _writer(ax, words)

    txt(w / 2, h - 0.17, "What can a thermodynamic cycle see?", 10.0, weight="bold")

    sage, coral = _darken(SAGE, 0.72), _darken(CORAL, 0.62)
    x0, x1 = 0.30, w - 0.82
    for y, frac, top_label, right_label in (
        (h - 0.72, 0.325, "the error space", "a third"),
        (h - 1.26, 0.006, "its error against experiment", "under 1%"),
    ):
        ax.add_patch(plt.Rectangle((x0, y), x1 - x0, 0.20, facecolor=_tint(BLUE, 0.86),
                                   edgecolor="none", zorder=1))
        ax.add_patch(plt.Rectangle((x0, y), max((x1 - x0) * frac, 0.012), 0.20,
                                   facecolor=sage, edgecolor="none", zorder=2))
        txt(x0, y + 0.29, top_label, 8.0, ha="left", color=_darken(BLUE, 0.45))
        txt(x1 + 0.06, y + 0.10, right_label, 8.5, ha="left", va="center",
            color=sage, weight="bold")

    txt(x0, 0.15, "auditable", 8.0, ha="left", color=sage)
    txt(x1, 0.15, "invisible to any cycle", 8.0, ha="right", color=coral)

    n_words = len(words)
    print(f"word count on table-of-contents graphic: {n_words}")
    assert n_words <= 24, f"too many words on the contents graphic: {n_words}"
    for ext in ("pdf", "png"):
        fig.savefig(ROOT / "figs" / f"toc_graphic.{ext}")
    plt.close(fig)
    print("wrote figs/toc_graphic.(pdf|png)")


def main():
    use_paper_style()
    overview()
    _style()
    toc()


if __name__ == "__main__":
    main()
