"""The single source of visual truth for this article's figures.

Import this module from every ``figs/make_*.py`` that produces a main-text or SI figure,
call :func:`use_paper_style` once at the top of ``main()``, size the figure with
:func:`figsize`, head every panel with :func:`panel`, and write it with :func:`finish`.
Nothing else in the repo should set ``matplotlib.rcParams``.

Physical size — the measured fact this module is built on
---------------------------------------------------------
``\\textwidth`` was measured (``\\the\\textwidth`` in both preambles) as **469.755pt =
exactly 6.5 in** in BOTH builds, ``paper/paper_jctc.tex`` and ``paper/paper_draft.tex``.
A figure authored 6.5 in wide and included at ``width=\\textwidth`` is therefore
reproduced by LaTeX at **scale 1.0**: its nominal point sizes are its rendered point
sizes, so 9 pt here is 9 pt on the page, next to 10 pt body text. Figures authored at
7.4–8.0 in are shrunk 12–19 % on inclusion, which is exactly why their type reads small.
Author at :data:`FULL`; never compensate by inflating font sizes.

The same arithmetic fixes the one other inclusion width this article uses. A figure
included at ``width=0.7\\textwidth`` is reproduced at scale ``0.7 * 6.5 / W``, so it comes
back at scale 1.0 only when it was authored ``0.7 * 6.5 = 4.55`` in wide: that is
:data:`NARROW`, for the genuinely narrow figure that would look absurd stretched across
the text block. There are exactly two canonical authored widths, :data:`FULL` and
:data:`NARROW`, and each is paired with one inclusion width in the ``.tex``.

The decidable colour rule
-------------------------
A colour means the SAME THING in every figure of this article, and the rule is decidable:
two people applying it to the same series must reach the same answer. "OURS vs FOIL" alone
did not, because two of our own methods sometimes compete and one series is sometimes split
by sign. The rule is therefore stated by quantity family, not by rhetorical role:

  OURS  #0072B2 blue    the calibrated per-edge bar and ANYTHING computed from it: the sandwich,
                        the calibrated detector, its flags, its |z|-guided removal rule, the
                        measured error against experiment. One quantity family, one colour.
  FOIL  #D55E00 orange  ONLY the overconfident / miscalibrated stand-in and real learned-variance
                        heads. Never a random baseline. Never our own estimator on systems where
                        it happens to under-report.
  ALT   #E69F00 amber   a second NAMED method that is neither ours-by-construction nor the foil,
                        e.g. the theorem-derived influence-ranked repair rule.
  MUTED #8C8C8C grey    random / null baselines, and de-emphasised data.
  REF   #6E6E6E grey    reference LINES only: ideal, chance, y=x, nominal level. Never a series.
  GOOD  #009E73 green   currently used in exactly one place and nowhere else; it is being retired
                        from the main-text figures.

  Two refinements the first pass needed and did not have:
  * a family with several members takes ONE hue and is told apart by dash pattern, marker and
    tint. Three learned heads are three FOIL tints, not three hues; otherwise the family reads
    as three unrelated things and a hue is spent per member.
  * a NAMED estimator this article corrects, as opposed to a null, is not grey. Grey means
    random, null, or de-emphasised data, and nothing else. The corrected rival takes THIRD, so
    that grey carries one meaning on facing pages.

The consequences, stated so two agents cannot disagree:

* the |z|-guided / residual-ranked removal rule is OURS blue in EVERY figure it appears in;
* the random-removal baseline is MUTED grey in EVERY figure it appears in;
* the influence-ranked (theorem-derived) rule is ALT amber;
* a single series must not be split into OURS and FOIL by the sign of its values; if a sign or a
  threshold crossing must be visible, use ``tint(OURS, 0.45)`` for the other side, same hue;
* one name for the overconfident construct in figure text: "stand-in" (never "stress test").
  A genuinely trained head is a "learned head" and keeps that name.

:data:`THIRD` remains available for a third named method, after ALT.

The hexes are the Okabe–Ito colour-blind-safe set and are the ones the repo's figures
already use, so adopting this module changes no figure's hue.

The minimum in-axes type size
-----------------------------
Because a figure authored at :data:`FULL` is reproduced at scale 1.0, the nominal point
size of a label IS its point size on the page, sitting beside 10 pt body text. An audit
found in-panel annotations at 6.5 pt in four SI figures — the smallest type anywhere in
the article, 35 % below the body text, and below what a printed page or a referee's
screen renders comfortably. :data:`MIN_ANNOTATION_PT` is therefore the floor: **no text
drawn inside an axes goes below it**. That covers in-panel annotations, tick labels,
axis labels and legend entries — everything :func:`check_min_type` walks — so a script
or a test can assert ``check_min_type(fig) == []`` before the figure is delivered. The
floor is the size of the panel subtitle and the house legend, which are already the
smallest type the style permits; nothing needs to be smaller than the smallest thing.

The annotation-figure exemption
-------------------------------
The semantic palette above governs figures whose colours encode METHODS. One figure in
this article does not: ``figs/analyze_qc_structures.py`` draws two protein structures,
where the colours are annotations of chemical OBJECTS — a catalytic aspartate dyad,
crystallographic waters, a recognition asparagine, the ligand. Forcing the semantic
palette onto it made it strictly less legible: the dyad and the waters both became blue,
and the asparagine became a pale-blue square barely separable from the blue water
circles, where the original red / green / blue / grey was instantly readable and is the
structural-biology field's own convention.

The exemption, stated so it cannot be applied by anyone who feels like it: a figure
**whose colours annotate objects rather than methods** keeps everything else in this
module — the house typography, the spines, the :func:`panel` headings, the canonical
widths and :func:`finish` — but chooses a legible categorical scheme instead of the
semantic palette. The palette does not apply there, because forcing it would make one
hue mean "a residue" in this figure and "the calibrated estimator" in every other one,
which is precisely the collision the palette exists to prevent. An exempt figure must
therefore contain no method series at all; the moment it compares two estimators, it is
back under the semantic rule.

:data:`CATEGORICAL` is the scheme such a figure should draw from — **for annotation
figures only**, never for a series that a reader could mistake for a method. Its hues
are Paul Tol's colour-blind-safe "bright" qualitative set, deliberately chosen so that
none of them is a hex from the semantic palette: an exempt figure's colours cannot then
be confused, by eye or by grep, with OURS / FOIL / ALT / THIRD.

Run ``python -m pytest tests/test_paperstyle.py`` after touching anything here.
"""
from __future__ import annotations

import pathlib
from collections.abc import Sequence

import matplotlib as mpl
from cycler import cycler
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.text import Text

__all__ = [
    "ALT",
    "CATEGORICAL",
    "CYCLE",
    "FIGDIR",
    "FOIL",
    "FULL",
    "GOOD",
    "HALF",
    "INK",
    "MIN_ANNOTATION_PT",
    "MUTED",
    "NARROW",
    "OURS",
    "REF",
    "THIRD",
    "check_min_type",
    "despine",
    "figsize",
    "finish",
    "legend",
    "panel",
    "reference_line",
    "tint",
    "use_paper_style",
]

# --------------------------------------------------------------------------------------
# 1. Canonical sizes (inches).
# --------------------------------------------------------------------------------------
FULL = 6.5   # = \textwidth (469.755 pt) in both builds; include at width=\textwidth
HALF = 3.2   # two figures side by side, each at width=0.48\textwidth
NARROW = 4.55  # = 0.7 x 6.5; include at width=0.7\textwidth, also reproduced at scale 1.0

FIGDIR = pathlib.Path(__file__).resolve().parent


def figsize(n_panels: int, height: float = 2.6) -> tuple[float, float]:
    """Return the canonical ``(width, height)`` for a full-width row of ``n_panels``.

    The width is always :data:`FULL`: a row of panels spans the text block, whether it
    holds one panel or four. ``n_panels`` is taken only to make the caller's intent
    explicit (and to reject a nonsense row); pick ``height`` so each panel lands near a
    2:1.6 aspect — roughly 2.2–2.6 in for a 3-panel row, 3.0–3.4 in for a 2-panel row.
    """
    if n_panels < 1:
        raise ValueError(f"n_panels must be >= 1, got {n_panels}")
    if height <= 0:
        raise ValueError(f"height must be > 0, got {height}")
    return (FULL, float(height))


# --------------------------------------------------------------------------------------
# 2. The semantic palette (Okabe–Ito). Meaning, not hue — see the module docstring.
# --------------------------------------------------------------------------------------
OURS = "#0072B2"   # the calibrated / sandwich / physics quantity: this article's own
FOIL = "#D55E00"   # the overconfident / learned / stand-in quantity it is contrasted with
REF = "#6E6E6E"    # reference lines: ideal, chance, y = x, nominal level, null
ALT = "#E69F00"    # a second baseline
THIRD = "#CC79A7"  # a third baseline
GOOD = "#009E73"   # an explicitly good / positive outcome
INK = "#222222"    # text and spines
MUTED = "#8C8C8C"  # de-emphasised data

#: The prop_cycle order: our quantity first, its foil second, then baselines.
CYCLE: tuple[str, ...] = (OURS, FOIL, REF, ALT, THIRD, GOOD)

#: For ANNOTATION FIGURES ONLY — see "The annotation-figure exemption" in the module
#: docstring. A figure whose colours annotate chemical or structural OBJECTS (residues,
#: waters, a ligand) rather than methods draws its hues from here, in this order, and the
#: semantic palette above does not apply to it. Paul Tol's colour-blind-safe "bright"
#: qualitative set: blue, red, green, yellow, purple, grey — none of them a hex from the
#: semantic palette, so an object colour can never be read as a method colour.
CATEGORICAL: tuple[str, ...] = (
    "#4477AA",  # blue
    "#EE6677",  # red
    "#228833",  # green
    "#CCBB44",  # yellow
    "#AA3377",  # purple
    "#BBBBBB",  # grey
)


def tint(color: str, amount: float = 0.5) -> str:
    """Lighten ``color`` toward white by ``amount`` in [0, 1], returning a hex string.

    ``amount=0`` returns the colour unchanged, ``amount=1`` returns white. Use it for
    bar fills, confidence bands and histogram faces so a fill and its line carry the
    same meaning at different weights — never to introduce a new meaning.
    """
    if not 0.0 <= amount <= 1.0:
        raise ValueError(f"amount must be in [0, 1], got {amount}")
    r, g, b = (int(color[i:i + 2], 16) for i in (1, 3, 5))
    mix = tuple(round(c + (255 - c) * amount) for c in (r, g, b))
    return "#{:02X}{:02X}{:02X}".format(*mix)


# --------------------------------------------------------------------------------------
# 3. The rcParams.
# --------------------------------------------------------------------------------------
def use_paper_style() -> None:
    """Install the house rcParams. Idempotent; call once per figure script."""
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 9,
        "axes.labelsize": 9.5,
        "axes.titlesize": 9.5,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 7.5,
        # mathtext must match the sans body text, not fall back to a serif face
        "mathtext.fontset": "dejavusans",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.edgecolor": INK,
        "axes.linewidth": 0.8,
        "axes.labelcolor": INK,
        "axes.axisbelow": True,
        "axes.prop_cycle": cycler(color=list(CYCLE)),
        "text.color": INK,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "xtick.major.size": 3,
        "ytick.major.size": 3,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "xtick.color": INK,
        "ytick.color": INK,
        "xtick.labelcolor": INK,
        "ytick.labelcolor": INK,
        "legend.frameon": False,
        "figure.dpi": 120,
        "savefig.dpi": 300,
        # the canvas is delivered uncropped at exactly FULL wide; see finish()
        "savefig.bbox": "standard",
        "savefig.pad_inches": 0.02,
    })


# --------------------------------------------------------------------------------------
# 4. The one panel-heading convention.
# --------------------------------------------------------------------------------------
_HEAD_PT = 4.0        # points above the axes top for the heading baseline
_HEAD_PT_SUB = 13.0   # ... when a subtitle occupies the line below it
_LETTER_GAP_PT = 11.0  # points from the letter's anchor to the title's anchor
_LETTER_CHAR_PT = 5.5  # ... widened by this much per EXTRA character in the letter


def panel(ax: Axes, letter: str, title: str | None = None, subtitle: str | None = None) -> None:
    """Draw the house panel heading: a bold ``letter``, optionally a short ``title``.

    The heading sits above the axes at the top-left, in axes coordinates with a
    points offset, so it lands in the same place at the same size in every figure of
    the article regardless of axes size or data limits. Use this instead of
    ``ax.set_title`` everywhere.

    ``title`` must be ONE short line (a handful of words). It is drawn as a single
    unwrapped line by design: nothing here shortens or wraps it, so an over-long title
    will simply run past the panel. Shortening it is the caller's job — move the detail
    into the axis labels or the caption.

    ``subtitle`` renders one further short line below the heading, in :data:`REF` at a
    smaller size; supplying it shifts the heading line up to make room.

    ``letter`` is a PANEL LETTER — one character, or two for a sub-panel such as ``A1``
    — and it is drawn bold. Passing a sentence there is rejected rather than drawn: a
    figure script once called ``panel(ax, "dose-response to $\\sigma$ miscalibration")``,
    which set the whole heading bold and broke the mathtext sigma's weight mid-line,
    because bold is a face the mathtext font set does not have. A heading's words belong
    in ``title``, which is the second parameter and is drawn at normal weight.
    """
    if len(letter) > 2 or any(ch.isspace() for ch in letter):
        raise ValueError(
            f"letter must be a panel letter of at most 2 characters with no whitespace, "
            f"got {letter!r}; pass the words as the `title` parameter instead: "
            f"panel(ax, 'A', {letter!r})"
        )
    dy = _HEAD_PT_SUB if subtitle else _HEAD_PT
    # The gap is measured from the letter's left edge, so a two-character sub-panel letter
    # needs one character's width more of it; at 11.0 flat, "B1" and its title collided.
    gap = _LETTER_GAP_PT + (len(letter) - 1) * _LETTER_CHAR_PT
    ax.annotate(
        letter,
        xy=(0.0, 1.0), xycoords="axes fraction",
        xytext=(0.0, dy), textcoords="offset points",
        ha="left", va="baseline", fontsize=9.5, fontweight="bold",
        color=INK, annotation_clip=False, clip_on=False,
    )
    if title:
        ax.annotate(
            title,
            xy=(0.0, 1.0), xycoords="axes fraction",
            xytext=(gap, dy), textcoords="offset points",
            ha="left", va="baseline", fontsize=9.5, fontweight="normal",
            color=INK, annotation_clip=False, clip_on=False,
        )
    if subtitle:
        ax.annotate(
            subtitle,
            xy=(0.0, 1.0), xycoords="axes fraction",
            xytext=(0.0, _HEAD_PT), textcoords="offset points",
            ha="left", va="baseline", fontsize=7.5,
            color=REF, annotation_clip=False, clip_on=False,
        )


# --------------------------------------------------------------------------------------
# 5. Saving.
# --------------------------------------------------------------------------------------
def finish(
    fig: Figure,
    path_stem: str,
    formats: Sequence[str] = ("pdf", "png"),
    layout: str | None = "constrained",
) -> None:
    """Save ``fig`` as ``figs/<path_stem>.<fmt>`` for each format, at the AUTHORED width.

    Takes a stem, not a path: every figure of this article lives next to this module.

    The canvas is never cropped. ``bbox_inches='tight'`` used to trim the canvas down to
    its content, which took a different amount off each figure depending on how far its
    y-label sat from the edge: the delivered PDFs measured 469.0, 467.4, 467.4, 464.9,
    463.7 and 451.2 pt across the six main-text figures. Each then had to be stretched
    back to ``\textwidth`` on inclusion, at a different scale factor, so a nominal 9.5 pt
    heading rendered at 9.52 pt in one figure and 9.89 pt in another — a 4 % spread, which
    silently defeats this module's whole premise. Instead the content is fitted INSIDE a
    fixed canvas: ``layout`` (default ``'constrained'``) is set on the figure as its layout
    engine, and the save passes ``bbox_inches=None`` under ``savefig.bbox='standard'`` so
    no cropping happens and ``savefig.pad_inches`` is irrelevant. Every regenerated PDF now
    measures **468.0 pt** wide (6.5 in × 72 PDF pt/in, exactly; the 469.755 in the module
    docstring is the same length in TeX pt of 1/72.27 in), verified by reading the
    ``/MediaBox`` of each delivered file: 468.0 +/- 1.0 pt.

    Pass ``layout=None`` to opt out — for a figure that places its axes manually (the
    graphical abstract places them in inches) and would be re-laid-out by a layout engine.
    A caller that opts out owns its own margins and must still deliver a 6.5 in canvas;
    nothing is cropped for it either.
    """
    if "/" in path_stem or path_stem.endswith((".pdf", ".png")):
        raise ValueError(f"path_stem must be a bare stem, got {path_stem!r}")
    if layout is not None:
        fig.set_layout_engine(layout)
    written = []
    # 'standard' is matplotlib's name for "do not crop"; bbox_inches=None otherwise falls
    # back to rcParams['savefig.bbox'], which the house style used to set to 'tight'.
    with mpl.rc_context({"savefig.bbox": "standard"}):
        for fmt in formats:
            out = FIGDIR / f"{path_stem}.{fmt}"
            fig.savefig(out, dpi=mpl.rcParams["savefig.dpi"], bbox_inches=None)
            written.append(out.name)
    print(f"wrote {', '.join(written)} -> {FIGDIR}")


# --------------------------------------------------------------------------------------
# 6. Small helpers.
# --------------------------------------------------------------------------------------
def reference_line(ax: Axes, kind: str, value: float = 0.0, **kw: object) -> Line2D:
    """Draw a house-styled reference line (ideal / chance / y = x / nominal / null).

    ``kind`` is ``'hline'`` or ``'vline'`` (drawn at ``value``), or ``'diagonal'``
    (y = x, drawn with ``axline`` so it tracks the limits). Keywords pass through, so
    ``label=`` and a different ``zorder`` still work; the colour is :data:`REF` and the
    line sits below the data.
    """
    style: dict[str, object] = {"color": REF, "linestyle": ":", "linewidth": 1.0, "zorder": 0.5}
    style.update(kw)
    if kind == "hline":
        return ax.axhline(value, **style)  # type: ignore[arg-type]
    if kind == "vline":
        return ax.axvline(value, **style)  # type: ignore[arg-type]
    if kind == "diagonal":
        return ax.axline((0.0, 0.0), slope=1.0, **style)  # type: ignore[arg-type]
    raise ValueError(f"kind must be 'hline', 'vline' or 'diagonal', got {kind!r}")


def legend(ax: Axes, **kw: object) -> object:
    """House legend: frameless and tight. ``loc=`` and friends pass through."""
    style: dict[str, object] = {
        "frameon": False,
        "handlelength": 1.6,
        "borderaxespad": 0.4,
        "labelspacing": 0.35,
    }
    style.update(kw)
    return ax.legend(**style)  # type: ignore[arg-type]


#: The floor for any text drawn INSIDE an axes, in points — see the module docstring.
#: At scale 1.0 this is a real 7.5 pt beside 10 pt body text; 6.5 pt annotations, which an
#: audit found in four SI figures, are below what the page and the referee's screen carry.
MIN_ANNOTATION_PT = 7.5


def _axes_text_artists(ax: Axes) -> list[Text]:
    """Every text artist that belongs to ``ax``: annotations, ticks, labels, legend."""
    artists: list[Text] = [*ax.texts, ax.title, ax.xaxis.label, ax.yaxis.label]
    artists.extend(ax.get_xticklabels())
    artists.extend(ax.get_yticklabels())
    leg = ax.get_legend()
    if leg is not None:
        artists.extend(leg.get_texts())
        leg_title = leg.get_title()
        if leg_title is not None:
            artists.append(leg_title)
    return artists


def check_min_type(fig: Figure, floor: float = MIN_ANNOTATION_PT) -> list[tuple[str, float]]:
    """Return the ``(text, size)`` pairs in ``fig`` drawn inside an axes below ``floor``.

    A figure is compliant when this returns an empty list, so a figure script or a test
    asserts ``assert check_min_type(fig) == []`` just before :func:`finish`. Empty and
    whitespace-only artists (matplotlib makes a tick label object per tick whether or not
    it carries text) are skipped, and the comparison is strict, so type exactly at the
    floor — the panel subtitle and the house legend, both 7.5 pt — is not flagged.
    """
    offenders: list[tuple[str, float]] = []
    for ax in fig.axes:
        for artist in _axes_text_artists(ax):
            label = artist.get_text()
            if not label.strip():
                continue
            size = float(artist.get_fontsize())
            if size < floor:
                offenders.append((label, size))
    return offenders


def despine(ax: Axes) -> None:
    """Remove the top and right spines. Idempotent.

    rcParams already do this for axes created normally; use it for axes that escape
    them, e.g. the second axes made by ``ax.twinx()``.
    """
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
