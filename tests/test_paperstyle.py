# tests/test_paperstyle.py
"""The house figure style is a contract: these are the clauses the six figure scripts rely on."""
import re

import matplotlib as mpl
import matplotlib.pyplot as plt
import pytest
from figs import paperstyle
from figs.paperstyle import (
    ALT,
    CATEGORICAL,
    CYCLE,
    FOIL,
    FULL,
    GOOD,
    HALF,
    INK,
    MIN_ANNOTATION_PT,
    MUTED,
    NARROW,
    OURS,
    REF,
    THIRD,
    check_min_type,
    despine,
    figsize,
    finish,
    legend,
    panel,
    reference_line,
    tint,
    use_paper_style,
)

mpl.use("Agg")


@pytest.fixture(autouse=True)
def _isolate_rcparams():
    """These tests mutate global rcParams; hand the suite back exactly what it had."""
    saved = mpl.rcParams.copy()
    yield
    mpl.rcParams.update(saved)


HEX = re.compile(r"^#[0-9A-Fa-f]{6}$")
ALL_COLORS = {
    "OURS": OURS, "FOIL": FOIL, "REF": REF, "ALT": ALT, "THIRD": THIRD,
    "GOOD": GOOD, "INK": INK, "MUTED": MUTED,
}


# ---------------------------------------------------------------- palette
def test_every_exported_colour_is_a_hex_string():
    for name, value in ALL_COLORS.items():
        assert isinstance(value, str), name
        assert len(value) == 7, f"{name}={value!r} is not 7 chars"
        assert HEX.match(value), f"{name}={value!r} is not #RRGGBB"


def test_the_semantic_colours_are_mutually_distinct():
    """A colour carries a meaning; two meanings must never share a hue."""
    semantic = [OURS, FOIL, REF, ALT, THIRD, GOOD]
    assert len({c.upper() for c in semantic}) == 6
    # the cycle is exactly those six, in the documented order (ours first, its foil second)
    assert CYCLE == (OURS, FOIL, REF, ALT, THIRD, GOOD)


# ---------------------------------------------------------------- tint
def _rgb(hexstr):
    return tuple(int(hexstr[i:i + 2], 16) for i in (1, 3, 5))


def test_tint_moves_strictly_toward_white():
    for base in (OURS, FOIL, GOOD, INK):
        light = _rgb(tint(base, 0.4))
        dark = _rgb(base)
        for c_light, c_dark in zip(light, dark, strict=True):
            assert c_light >= c_dark
        assert light != dark, f"tint({base}, 0.4) did not lighten"
        # and further at a larger amount
        lighter = _rgb(tint(base, 0.8))
        assert sum(lighter) > sum(light)


def test_tint_endpoints():
    assert tint(OURS, 0.0) == OURS                      # identity at 0
    assert tint(FOIL, 1.0).lower() == "#ffffff"         # white at 1
    assert tint(tint(OURS, 0.0), 0.0) == OURS           # idempotent at 0


def test_tint_rejects_out_of_range():
    with pytest.raises(ValueError):
        tint(OURS, 1.5)


# ---------------------------------------------------------------- sizes
def test_full_width_is_the_measured_textwidth():
    # \the\textwidth = 469.755pt in BOTH paper/paper_jctc.tex and paper/paper_draft.tex,
    # and 469.755pt / 72.27 pt-per-inch = 6.500 in. Authoring at FULL => scale 1.0.
    assert FULL == 6.5
    assert HALF < FULL
    assert figsize(3, 2.4) == (6.5, 2.4)
    assert figsize(1)[0] == FULL          # one panel still spans the text block
    with pytest.raises(ValueError):
        figsize(0)


def test_narrow_width_is_seven_tenths_of_the_textwidth():
    r"""The only other inclusion width in the article: 0.7 x 6.5 = 4.55 in => scale 1.0."""
    assert NARROW == 4.55
    assert NARROW == pytest.approx(0.7 * FULL)
    assert HALF < NARROW < FULL


# ---------------------------------------------------------------- rcParams
def test_use_paper_style_sets_named_rcparams():
    mpl.rcdefaults()
    use_paper_style()
    assert mpl.rcParams["font.size"] == 9
    assert mpl.rcParams["axes.spines.right"] is False
    assert mpl.rcParams["savefig.dpi"] == 300
    # ACS asks for Helvetica or Arial in artwork, and mathtext is artwork: the custom set is
    # what keeps everything inside $...$ in the same face as the labels around it. Naming the
    # set is not enough on its own, so `make acscheck` reads the rendered PDFs as well.
    assert mpl.rcParams["mathtext.fontset"] == "custom"
    assert mpl.rcParams["mathtext.rm"] == "Arial"
    assert mpl.rcParams["mathtext.it"] == "Arial:italic"
    assert mpl.rcParams["font.sans-serif"][0] == "Arial"
    assert [d["color"] for d in mpl.rcParams["axes.prop_cycle"]] == list(CYCLE)


def test_use_paper_style_is_idempotent():
    use_paper_style()
    keys = ["font.size", "axes.labelsize", "legend.fontsize", "xtick.direction",
            "axes.linewidth", "savefig.pad_inches"]
    before = {k: mpl.rcParams[k] for k in keys}
    use_paper_style()
    after = {k: mpl.rcParams[k] for k in keys}
    assert before == after


# ---------------------------------------------------------------- panel headings
def test_panel_adds_exactly_the_expected_text_artists():
    fig, axes = plt.subplots(1, 3)
    panel(axes[0], "A")
    assert len(axes[0].texts) == 1

    panel(axes[1], "B", "sandwich vs foil")
    assert len(axes[1].texts) == 2

    panel(axes[2], "C", "stopping", subtitle="30 calls")
    assert len(axes[2].texts) == 3
    plt.close(fig)


def test_panel_does_not_use_set_title_and_is_bold_at_the_left():
    fig, ax = plt.subplots()
    panel(ax, "A", "ours")
    assert ax.get_title() == "", "panel() must not go through ax.set_title"
    letter = ax.texts[0]
    assert letter.get_text() == "A"
    assert letter.get_fontweight() == "bold"
    assert letter.get_ha() == "left"
    assert letter.get_position()[0] == 0.0   # flush with the axes' left edge
    assert ax.texts[1].get_fontweight() == "normal"
    plt.close(fig)


def test_panel_subtitle_is_reference_grey_and_smaller():
    fig, ax = plt.subplots()
    panel(ax, "A", "ours", subtitle="n = 30")
    head, sub = ax.texts[0], ax.texts[2]
    assert sub.get_color() == REF
    assert sub.get_fontsize() < head.get_fontsize()
    plt.close(fig)


# ---------------------------------------------------------------- panel() guard
def test_panel_rejects_a_sentence_passed_as_the_letter():
    """make_figP4b called panel(ax, "dose-response to $\\sigma$ miscalibration").

    The whole heading then rendered bold and the mathtext sigma broke weight mid-line.
    A panel letter is a letter; the words belong in the title parameter.
    """
    fig, ax = plt.subplots()
    with pytest.raises(ValueError) as excinfo:
        panel(ax, r"dose-response to $\sigma$ miscalibration")
    # the message must say where a sentence does belong
    assert "title" in str(excinfo.value)

    # two ways of being wrong, each caught on its own: too long, and whitespace
    with pytest.raises(ValueError):
        panel(ax, "abc")            # 3 characters, no whitespace
    with pytest.raises(ValueError):
        panel(ax, "A ")             # 2 characters, but one is a space
    assert len(ax.texts) == 0, "a rejected heading must draw nothing"
    plt.close(fig)


def test_panel_still_accepts_a_one_or_two_character_letter():
    fig, axes = plt.subplots(1, 2)
    panel(axes[0], "A", "dose-response to $\\sigma$ miscalibration")
    assert axes[0].texts[0].get_text() == "A"
    assert axes[0].texts[0].get_fontweight() == "bold"
    assert axes[0].texts[1].get_fontweight() == "normal"   # the sentence is NOT bold

    panel(axes[1], "A1", "sub-panel")                       # two characters is legal
    assert axes[1].texts[0].get_text() == "A1"
    plt.close(fig)


# ---------------------------------------------------------------- minimum type size
def test_min_annotation_pt_is_the_floor_the_house_style_already_respects():
    """7.5 pt: the panel subtitle and the house legend, the smallest type the style sets."""
    assert MIN_ANNOTATION_PT == 7.5
    use_paper_style()
    assert mpl.rcParams["legend.fontsize"] >= MIN_ANNOTATION_PT
    assert mpl.rcParams["xtick.labelsize"] >= MIN_ANNOTATION_PT
    assert mpl.rcParams["font.size"] >= MIN_ANNOTATION_PT


def test_check_min_type_finds_the_deliberately_small_label_and_nothing_else():
    """The 6.5 pt in-panel annotation the audit found by eye, now found by a test."""
    use_paper_style()
    fig, ax = plt.subplots(figsize=figsize(1, 2.0))
    ax.plot([0.0, 1.0], [0.0, 1.0], color=OURS, label="ours")
    ax.set_ylabel("y")
    legend(ax, loc="lower right")
    panel(ax, "A", "throwaway", subtitle="n = 30")          # subtitle sits AT the floor
    ax.text(0.5, 0.5, "6.5 pt annotation", fontsize=6.5)    # the one offender
    fig.canvas.draw()                                       # so tick labels carry text

    offenders = check_min_type(fig)
    assert offenders == [("6.5 pt annotation", 6.5)]
    plt.close(fig)


def test_check_min_type_passes_a_figure_that_respects_the_floor():
    use_paper_style()
    fig, ax = plt.subplots(figsize=figsize(1, 2.0))
    ax.plot([0.0, 1.0], [0.0, 1.0], color=OURS, label="ours")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    legend(ax, loc="lower right")
    panel(ax, "A", "throwaway", subtitle="n = 30")
    ax.text(0.5, 0.5, "at the floor", fontsize=MIN_ANNOTATION_PT)   # strict: not flagged
    fig.canvas.draw()

    assert check_min_type(fig) == []
    plt.close(fig)


def test_check_min_type_also_walks_ticks_labels_and_the_legend():
    """"Inside an axes" is not only ax.text: a 6 pt tick label is 6 pt type on the page."""
    use_paper_style()
    fig, ax = plt.subplots(figsize=figsize(1, 2.0))
    ax.plot([0.0, 1.0], [0.0, 1.0], color=OURS, label="ours")
    ax.set_xticks([0.0, 1.0])
    ax.set_xticklabels(["lo", "hi"], fontsize=6.0)
    ax.set_ylabel("y", fontsize=7.0)
    legend(ax, loc="lower right", fontsize=6.5)
    fig.canvas.draw()

    found = dict(check_min_type(fig))
    assert found["lo"] == 6.0 and found["hi"] == 6.0
    assert found["y"] == 7.0
    assert found["ours"] == 6.5
    plt.close(fig)


# ---------------------------------------------------------------- annotation exemption
def test_categorical_is_a_short_ordered_set_of_valid_distinct_hexes():
    """For annotation figures only -- objects, not methods; see the module docstring."""
    assert isinstance(CATEGORICAL, tuple)
    assert 4 <= len(CATEGORICAL) <= 6
    for value in CATEGORICAL:
        assert isinstance(value, str)
        assert HEX.match(value), f"{value!r} is not #RRGGBB"
    assert len({c.upper() for c in CATEGORICAL}) == len(CATEGORICAL), "hues must be distinct"


def test_categorical_shares_no_hue_with_the_semantic_palette():
    """An object colour must not be readable -- by eye or by grep -- as a method colour."""
    semantic = {c.upper() for c in ALL_COLORS.values()}
    assert semantic.isdisjoint({c.upper() for c in CATEGORICAL})


# ---------------------------------------------------------------- helpers
def test_reference_line_is_house_styled_and_below_the_data():
    fig, ax = plt.subplots()
    line = reference_line(ax, "hline", 1.0)
    assert line.get_color() == REF
    assert line.get_linestyle() == ":"
    assert line.get_zorder() < 1.0
    assert reference_line(ax, "vline", 0.5).get_color() == REF
    assert reference_line(ax, "diagonal").get_color() == REF
    with pytest.raises(ValueError):
        reference_line(ax, "wiggle")
    plt.close(fig)


def test_legend_is_frameless():
    fig, ax = plt.subplots()
    ax.plot([0, 1], [0, 1], label="ours")
    leg = legend(ax, loc="lower right")
    assert leg.get_frame_on() is False
    plt.close(fig)


def test_despine_is_idempotent():
    fig, ax = plt.subplots()
    despine(ax)
    despine(ax)
    assert ax.spines["top"].get_visible() is False
    assert ax.spines["right"].get_visible() is False
    assert ax.spines["left"].get_visible() is True
    plt.close(fig)


# ---------------------------------------------------------------- the delivered canvas
#: 6.5 in × 72 PDF pt/in. PDF points are 1/72 in; the 469.755 quoted for \textwidth is the
#: same length measured in TeX pt of 1/72.27 in.
FULL_PT = 468.0

#: 4.55 in x 72 = 327.6 pt, i.e. 0.7 x 468.0 -- the canvas a 0.7\textwidth figure is authored on.
NARROW_PT = 327.6


def _media_box_width_pt(path):
    """Width of page 1's /MediaBox, in PDF points. pypdf if it is installed, else bytes."""
    try:
        import pypdf
    except ImportError:
        box = re.search(rb"/MediaBox\s*\[([^\]]*)\]", path.read_bytes())
        assert box, f"no /MediaBox in {path}"
        x0, _y0, x1, _y1 = (float(v) for v in box.group(1).split())
        return x1 - x0
    page = pypdf.PdfReader(str(path)).pages[0]
    return float(page.mediabox.width)


def _probe(tmp_path, stem, ylabel, **kw):
    """Render a throwaway figure through finish() and return its delivered width in pt."""
    fig, ax = plt.subplots(figsize=figsize(1, 2.0))
    ax.plot([0.0, 1.0], [0.0, 1.0], color=OURS)
    ax.set_ylabel(ylabel)
    panel(ax, "A", "throwaway")
    finish(fig, stem, formats=("pdf",), **kw)
    plt.close(fig)
    return _media_box_width_pt(tmp_path / f"{stem}.pdf")


def test_finish_delivers_exactly_the_authored_width(tmp_path, monkeypatch):
    r"""The premise of the module: authored at FULL, included at \textwidth, scale 1.0.

    A cropped canvas has to be stretched back on inclusion by a figure-dependent factor,
    so a nominal 9.5 pt heading is no longer a rendered 9.5 pt heading.
    """
    monkeypatch.setattr(paperstyle, "FIGDIR", tmp_path)
    use_paper_style()
    assert abs(_probe(tmp_path, "probe", "y") - FULL_PT) <= 1.0


def test_finish_delivers_the_narrow_width_for_a_narrow_figure(tmp_path, monkeypatch):
    r"""Authored at NARROW, included at width=0.7\textwidth, also reproduced at scale 1.0."""
    monkeypatch.setattr(paperstyle, "FIGDIR", tmp_path)
    use_paper_style()

    fig, ax = plt.subplots(figsize=(NARROW, 3.0))
    ax.plot([0.0, 1.0], [0.0, 1.0], color=OURS)
    ax.set_ylabel("y")
    panel(ax, "A", "throwaway")
    finish(fig, "probe_narrow_canvas", formats=("pdf",))
    plt.close(fig)

    delivered = _media_box_width_pt(tmp_path / "probe_narrow_canvas.pdf")
    assert abs(delivered - NARROW_PT) <= 1.0


def test_finish_width_does_not_depend_on_how_far_the_ylabel_sits_from_the_edge(
        tmp_path, monkeypatch):
    """This is exactly what bbox_inches='tight' got wrong: 469.0 ... 451.2 pt across six figures."""
    monkeypatch.setattr(paperstyle, "FIGDIR", tmp_path)
    use_paper_style()
    narrow = _probe(tmp_path, "probe_narrow", "y")
    wide = _probe(tmp_path, "probe_wide", "an emphatically long y-axis label (kcal/mol)")
    assert abs(narrow - wide) < 0.01
    assert abs(wide - FULL_PT) <= 1.0


def test_finish_installs_the_layout_engine_and_a_caller_can_opt_out(tmp_path, monkeypatch):
    """Content is fitted inside the fixed canvas; a caller placing axes by hand may opt out."""
    monkeypatch.setattr(paperstyle, "FIGDIR", tmp_path)
    use_paper_style()

    fig, ax = plt.subplots(figsize=figsize(1, 2.0))
    ax.plot([0.0, 1.0], [0.0, 1.0])
    finish(fig, "engine", formats=("pdf",))
    assert fig.get_layout_engine() is not None
    plt.close(fig)

    # layout=None leaves the caller's own placement alone -- and still does not crop
    manual = plt.figure(figsize=(FULL, 2.0))
    manual.add_axes([0.12, 0.18, 0.84, 0.72]).plot([0.0, 1.0], [0.0, 1.0])
    finish(manual, "manual", formats=("pdf",), layout=None)
    assert manual.get_layout_engine() is None
    plt.close(manual)
    assert abs(_media_box_width_pt(tmp_path / "manual.pdf") - FULL_PT) <= 1.0


def test_finish_still_rejects_a_path_that_is_not_a_bare_stem(tmp_path, monkeypatch):
    monkeypatch.setattr(paperstyle, "FIGDIR", tmp_path)
    fig = plt.figure(figsize=(FULL, 2.0))
    with pytest.raises(ValueError):
        finish(fig, "sub/dir/figX")
    with pytest.raises(ValueError):
        finish(fig, "figX.pdf")
    plt.close(fig)
