"""Fig P4b -- dose-response of the cycle-closure QC to sigma miscalibration.

P4 showed the learned head's measured profile (all ratios in 0.09-0.20x) flags 42/48 systems
versus 6/48 calibrated, and that this is driven by MAGNITUDE rather than heterogeneity: any ratio
below ~0.20 inflates every system's chi-square by >=25x, which saturates the Benjamini-Hochberg
test. That confirms the referees' point that flagging is arithmetic in that band, and it leaves
open the question this script answers: at what shrink magnitude does the test actually stop being
selective?

Sweeps a frozen grid of UNIFORM se scale factors through the identical GLS + BH-FDR path and
reports the flagged fraction at each. Two frozen readouts:
  s_onset = the largest scale whose flagged count strictly exceeds the calibrated (x1) count
  s50     = the largest scale at which at least half the analysed systems are flagged
The learned head's measured band (min/max of the rank-transferred profile) is marked on the curve,
so the reader can see where a real head sits relative to the transition.

DESCRIPTIVE, not a hypothesis test: the whole curve is reported whatever its shape.

Run:  PYTHONPATH=src python figs/make_figP4b.py   (or `make figP4b`).  Deterministic.
"""
from __future__ import annotations

import csv
import math
import pathlib
import sys
from collections import defaultdict

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
from matplotlib.ticker import NullLocator  # noqa: E402

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from bar.qc import benjamini_hochberg, chi2_sf, gls_network  # noqa: E402
from bar.sigma_profile import PROFILE_POINTS, rank_transfer  # noqa: E402
from paperstyle import (  # noqa: E402
    FOIL, MUTED, NARROW, OURS, finish, legend, panel, reference_line, tint, use_paper_style,
)

DATA = ROOT / "data" / "openfe_replicates" / "combined_pymbar4_edge_data.csv"
FIGDIR = pathlib.Path(__file__).resolve().parent
# FROZEN pre-registration
SCALES = [2.0, 1.5, 1.3, 1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.35, 0.3, 0.25, 0.20, 0.15, 0.10]
# Semantic colours (paperstyle). The swept curve is NOT one method: it runs from the
# overconfident stand-in end (x0.15 and below) to the conservative end (x1.3, x2), and this
# article gives those ends opposite meanings, so painting the whole sweep in OURS blue would
# claim the calibrated bar's hue for the stand-in's own regime. The continuous curve is
# therefore de-emphasised data (MUTED), and the two sigma models the article NAMES are marked
# on it in the hues Fig L panel B gives them: FOIL for the x0.15 stand-in (and for the learned
# head's measured band, the same overconfident family), OURS for x1 as shipped. The two
# vertical markers are REF reference lines, told apart by dash pattern (see
# paperstyle.reference_line).
NAMED_SCALES = {0.15: "stand-in", 1.0: "as shipped"}   # the same two names Fig L panel B uses


def _f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return math.nan


def edge_val(r, k):
    cD = _f(r[f"complex_repeat_{k}_DG (kcal/mol)"]); cd = _f(r[f"complex_repeat_{k}_dDG (kcal/mol)"])
    sD = _f(r[f"solvent_repeat_{k}_DG (kcal/mol)"]);  sd = _f(r[f"solvent_repeat_{k}_dDG (kcal/mol)"])
    if any(math.isnan(v) for v in (cD, cd, sD, sd)):
        return None
    return cD - sD, math.sqrt(cd * cd + sd * sd)


def edge_overlap(r, k):
    c = _f(r[f"complex_repeat_{k}_smallest_overlap"])
    s = _f(r[f"solvent_repeat_{k}_smallest_overlap"])
    if math.isnan(c) or math.isnan(s):
        return math.nan
    return min(c, s)


def load():
    by = defaultdict(list)
    for r in csv.DictReader(open(DATA)):
        by[r["system name"]].append(r)
    out = {}
    for name, rows in sorted(by.items()):
        e, ov = [], []
        for r in rows:
            v = edge_val(r, 0)
            o = edge_overlap(r, 0)
            if v and not math.isnan(o):
                e.append((r["ligand_A"], r["ligand_B"], v[0], v[1]))
                ov.append(o)
        if len(e) >= 3 and gls_network(e).dof >= 1:
            out[name] = (e, np.array(ov))
    return out


def flagged_at(systems, scale):
    """Flagged count under a UNIFORM se scale factor, through the identical GLS + BH-FDR path."""
    names, ps = [], []
    for name, (edges, _ov) in sorted(systems.items()):
        scaled = [(a, b, y, se * scale) for (a, b, y, se) in edges]
        fit = gls_network(scaled)
        if fit.dof < 1:
            continue
        names.append(name)
        ps.append(chi2_sf(fit.chi2, fit.dof))
    fl = benjamini_hochberg(ps)
    return int(fl.sum()), len(names)


def main():
    use_paper_style()
    systems = load()
    curve = [(s, *flagged_at(systems, s)) for s in SCALES]
    tot = curve[0][2]
    cal = next(n for s, n, _t in curve if s == 1.0)

    above = [s for s, n, _t in curve if n > cal]
    s_onset = max(above) if above else math.nan
    half = [s for s, n, _t in curve if n >= math.ceil(tot / 2)]
    s50 = max(half) if half else math.nan

    # where the learned head's measured band sits
    all_ov = np.concatenate([systems[n][1] for n in sorted(systems)])
    band = rank_transfer(all_ov)
    b_lo, b_hi = float(band.min()), float(band.max())

    # A single curve on one narrow panel: authored at NARROW = 0.7 x 6.5 in and included at
    # width=0.7\textwidth, so it is reproduced at scale 1.0 like every other figure here.
    fig, ax = plt.subplots(figsize=(NARROW, 3.1))
    xs = [c[0] for c in curve]
    ys = [100.0 * c[1] / c[2] for c in curve]
    # The band is the learned head's own measured range, so it is the FOIL hue; it is drawn
    # OPAQUE at a pale tint and behind everything (zorder 0) rather than as a translucent wash,
    # so it never dims the markers that sit inside it.
    ax.axvspan(b_lo, b_hi, color=tint(FOIL, 0.86), lw=0, zorder=0)
    ax.plot(xs, ys, "o-", color=MUTED, lw=1.6, ms=4.0, zorder=3)
    # the two named sigma models, in their own hues, on top of the de-emphasised sweep
    for sc, note in NAMED_SCALES.items():
        col = FOIL if sc < 1.0 else OURS
        yv = 100.0 * next(n for s_, n, _t in curve if s_ == sc) / tot
        ax.plot([sc], [yv], "o", color=col, ms=6.0, zorder=5)
        # x1 is labelled to the RIGHT of its marker: set above it the words ran into the
        # curve, which is still descending through that height just to the left.
        off, ha, va = ((0, 10), "center", "bottom") if sc < 1.0 else ((9, 0), "left", "center")
        ax.annotate(note, xy=(sc, yv), xytext=off, textcoords="offset points",
                    ha=ha, va=va, fontsize=7.5, color=col, zorder=6)
    # the guides stop below the legend strip: they span the full plot height otherwise, and
    # then no in-axes corner is free of them -- which is what drove the legend onto the band.
    keys = [Patch(facecolor=tint(FOIL, 0.86), edgecolor=FOIL, lw=0.8,
                  label=f"learned-head band {b_lo:.2f}–{b_hi:.2f}×"),
            reference_line(ax, "vline", 1.0, ymax=0.76, ls=(0, (1.0, 2.0)), lw=1.0,
                           label="calibrated ×1")]
    if math.isfinite(s50):
        keys.append(reference_line(ax, "vline", s50, ymax=0.76, ls=(0, (4.5, 2.0)), lw=1.0,
                                   label=f"$s_{{50}}$ = {s50:g}"))
    ax.set_xscale("log")
    # the default log locator labels only the two decade ends (10^-1, 10^0) and leaves the
    # sweep's own scale factors unnamed; name the grid points the reader is asked to read off.
    ax.set_xticks([0.1, 0.2, 0.3, 0.5, 1.0, 2.0])
    ax.set_xticklabels(["0.1", "0.2", "0.3", "0.5", "1", "2"])
    ax.xaxis.set_minor_locator(NullLocator())
    ax.set_xlim(0.085, 2.4)
    ax.set_xlabel("uniform se scale factor (lower = more overconfident)")
    ax.set_ylabel("% of systems flagged")
    ax.set_ylim(-4, 100)
    ax.set_yticks([0, 25, 50, 75, 100])
    ax.spines["left"].set_bounds(0, 100)
    # upper right: the curve is low wherever the scale is generous, so nothing but the two
    # vertical guides reached that corner, and those now stop below it. At lower left the
    # legend sat ON the learned-head band and its first entry straddled the band's own edge.
    legend(ax, handles=keys, loc="upper right", handlelength=2.6)
    # one panel, so the heading is the title alone, set at the axes' top-left corner where
    # every other heading in the article starts; the panel-letter slot is left empty, as in
    # the article's other single-panel figure. The words go in `title`, at normal weight:
    # passed as the LETTER they were set bold, which broke the mathtext sigma's weight.
    panel(ax, "", r"dose-response to $\sigma$ miscalibration")

    finish(fig, "figP4b_dose_response")

    print(f"\n[P4b] systems={tot}; calibrated (x1) flags {cal}/{tot}")
    print(f"[P4b] profile points: {PROFILE_POINTS}; measured band {b_lo:.4f}-{b_hi:.4f}")
    print(f"{'scale':>7} {'flagged':>9} {'pct':>6}")
    for s, n, t in curve:
        print(f"{s:>7g} {f'{n}/{t}':>9} {100.0 * n / t:>5.0f}%")
    print(f"[P4b] s_onset (largest scale flagging more than calibrated) = {s_onset:g}")
    print(f"[P4b] s50 (largest scale flagging >= half of {tot}) = {s50:g}")
    print(f"[P4b] the measured head band ({b_lo:.2f}-{b_hi:.2f}) lies "
          f"{'BELOW' if b_hi < s50 else 'AT OR ABOVE'} s50")


if __name__ == "__main__":
    main()
