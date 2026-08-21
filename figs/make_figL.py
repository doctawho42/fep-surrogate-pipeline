"""Fig L -- Calibrated cycle-closure quality control (the impact result).

A thermodynamic perturbation network must close its cycles: for any loop, the signed sum of
ddG must be zero in expectation. The observed closure is sampling error PLUS any edge-level
systematic error (non-convergence, hysteresis, charge/water artifacts). The differentiable BAR
bottleneck supplies a *calibrated* per-edge aleatoric variance V_e = B_e/I_e^2, which is exactly
the null model that turns cycle closure into a properly-calibrated statistical test:

  GLS-fit node potentials with weights 1/V_e; residual r_e; X^2 = sum_e r_e^2/V_e ~ chi^2(dof),
  dof = #independent cycles (E - rank of the incidence matrix). Reduced chi^2 >> 1 => the network
  is over-dispersed relative to sampling => edge-level systematic error.

Three claims, all on the public OpenFE IndustryBenchmarks2024 replicate set (1145 edges, 34
systems, 3 independent replicates):
  A. Deployable: with the calibrated sandwich null, most systems are sampling-consistent
     (median reduced chi^2 ~ 0.34; <1 because the bars are ~1.3x conservative, Fig A-rep), and
     a BH-FDR test flags a chemically-sensible minority (brd4, bace, faah, cdk8, hif2a, p38).
  B. The detector's usefulness IS sigma-calibration: an overconfident learned sigma (x0.15,
     Fig A) flags ~everything (FPR->1, useless); the calibrated sandwich is selective. This is
     the concrete task an overconfident sigma provably ruins and a calibrated one enables.
  C. The 3 independent replicates separate SYSTEMATIC from SAMPLING: pooling the replicates
     (se -> se/sqrt3) *increases* the flagged systems' reduced chi^2 toward the 3x variance-
     components ceiling -- happens only if the closure error is reproducible (systematic), not
     if the bars are merely too small (which would leave the ratio at 1).

Cycle closure detects EDGE-level inconsistency; node-consistent force-field bias cancels around
a loop and is (correctly) invisible here. Uses the reported pymbar4 MBAR se, equal to the
sandwich B/I^2 to leading order (Fig A); the differentiable sandwich is what makes V_e a
per-edge, graph-native, backprop-able weight rather than a returned number.

Run:  PYTHONPATH=src python figs/make_figL.py   (or `make figL`).  Deterministic.
"""
from __future__ import annotations

import csv
import math
import pathlib
import sys

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict
from matplotlib.lines import Line2D

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from bar.qc import benjamini_hochberg as bh_flags, chi2_sf, gls_network  # noqa: E402
from paperstyle import (  # noqa: E402
    FOIL, INK, MUTED, OURS, REF, figsize, finish, legend, panel, reference_line, tint,
    use_paper_style,
)

DATA = ROOT / "data" / "openfe_replicates" / "combined_pymbar4_edge_data.csv"
FIGDIR = pathlib.Path(__file__).resolve().parent
# Semantic colours (paperstyle): OURS = what the calibrated detector picks out (the flagged
# systems -- the same blue means the same flags in Fig L-validation), MUTED = the sampling-
# consistent majority, which is present but not the point, FOIL = the overconfident sigma the
# test is contrasted with, REF = the null and ceiling reference lines.


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
    return cD - sD, math.sqrt(cd * cd + sd * sd)   # binding ddG, sandwich/MBAR se


def gls_chi2(edges):
    """Thin adapter over the tested API bar.qc.gls_network: returns (X2, dof, standardized z)."""
    fit = gls_network(edges)
    return fit.chi2, fit.dof, fit.z


def load_systems():
    by = defaultdict(list)
    for r in csv.DictReader(open(DATA)):
        by[r["system name"]].append(r)
    return by


def analyze():
    by = load_systems()
    per_sys = []           # replicate-0 deployable test
    disc = []              # single vs pooled-3 (systematic vs sampling)
    for sys, rs in sorted(by.items()):
        e0 = [(r["ligand_A"], r["ligand_B"], *edge_val(r, 0)) for r in rs if edge_val(r, 0)]
        if len(e0) < 3:
            continue
        X2, dof, _ = gls_chi2(e0)
        if dof < 1:
            continue
        per_sys.append(dict(sys=sys, E=len(e0), dof=dof, rc=X2 / dof, p=chi2_sf(X2, dof)))
        singles = []
        for k in (0, 1, 2):
            ek = [(r["ligand_A"], r["ligand_B"], *edge_val(r, k)) for r in rs if edge_val(r, k)]
            if len(ek) >= 3:
                x, d, _ = gls_chi2(ek)
                if d >= 1:
                    singles.append(x / d)
        ep = []
        for r in rs:
            vs = [edge_val(r, k) for k in (0, 1, 2)]; vs = [v for v in vs if v]
            if len(vs) < 2:
                continue
            ep.append((r["ligand_A"], r["ligand_B"], float(np.mean([v[0] for v in vs])),
                       math.sqrt(np.mean([v[1] ** 2 for v in vs]) / len(vs))))
        if singles and len(ep) >= 3:
            xp, dp, _ = gls_chi2(ep)
            if dp >= 1:
                disc.append(dict(sys=sys, single=float(np.mean(singles)), pooled=xp / dp))
    flags = bh_flags([d["p"] for d in per_sys])
    for d, fl in zip(per_sys, flags):
        d["flag"] = bool(fl)

    # calibration sweep: fraction of systems flagged vs sigma-scale
    sweep = []
    for name, sc in [("uniform x0.15\n(stress test)", 0.15), ("reported bar\n(x1)", 1.0),
                     ("conservative\n(x1.3)", 1.30), ("too wide\n(x2)", 2.0)]:
        ps = []
        for sys, rs in sorted(by.items()):
            e0 = [(r["ligand_A"], r["ligand_B"], edge_val(r, 0)[0], edge_val(r, 0)[1] * sc)
                  for r in rs if edge_val(r, 0)]
            if len(e0) < 3:
                continue
            X2, dof, _ = gls_chi2(e0)
            if dof >= 1:
                ps.append(chi2_sf(X2, dof))
        sweep.append((name, sc, int(bh_flags(ps).sum()), len(ps)))
    return per_sys, disc, sweep


def _overlap(a, b):
    """Area shared by two display-space bounding boxes."""
    w = min(a.x1, b.x1) - max(a.x0, b.x0)
    h = min(a.y1, b.y1) - max(a.y0, b.y0)
    return w * h if (w > 0 and h > 0) else 0.0


_DIRS = [(math.cos(a), math.sin(a)) for a in
         (i * math.pi / 8 for i in range(16))]


def place_labels(ax, xy, texts, color, obstacles, keep_clear=(), fontsize=7.5,
                 radii=(6.0, 9.0, 13.0, 18.0)):
    """Label points without collisions -- a small local ``adjust_text``.

    Every direction at every radius is tried for each label; a candidate is charged for the
    data points it would cover, for its overlap with labels already placed and with the
    ``keep_clear`` boxes, for leaving the axes, and mildly for its distance from its own
    point. The cheapest candidate wins, and a label that ends up far from its point is tied
    back to it with a hairline leader. Placement only: the same strings are drawn for the
    same points whatever it decides.
    """
    fig = ax.figure
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    obs = ax.transData.transform(np.asarray(obstacles, dtype=float))
    taken = list(keep_clear)
    axbb = ax.get_window_extent(rend)
    for (x, y), label in zip(xy, texts):
        ann = ax.annotate(label, (x, y), xytext=(0.0, 0.0), textcoords="offset points",
                          fontsize=fontsize, color=color, zorder=6, annotation_clip=False)
        best = None
        for r in radii:
            for ux, uy in _DIRS:
                ann.set_position((ux * r, uy * r))
                ann.set_ha("left" if ux > 0.3 else "right" if ux < -0.3 else "center")
                ann.set_va("bottom" if uy > 0.3 else "top" if uy < -0.3 else "center")
                box = ann.get_window_extent(rend).expanded(1.10, 1.45)
                cost = 25.0 * sum(box.x0 <= px <= box.x1 and box.y0 <= py <= box.y1
                                  for px, py in obs)
                cost += sum(_overlap(box, t) for t in taken) ** 0.5 * 4.0
                cost += 500.0 * (box.x0 < axbb.x0 or box.x1 > axbb.x1
                                 or box.y0 < axbb.y0 or box.y1 > axbb.y1)
                cx, cy = (box.x0 + box.x1) / 2.0, (box.y0 + box.y1) / 2.0
                near = 25.0 * fig.dpi / 72.0          # 25 pt, in display units
                cost += 3.0 * sum(math.hypot(px - cx, py - cy) < near for px, py in obs)
                cost += 0.55 * r
                if best is None or cost < best[0]:
                    best = (cost, (ux * r, uy * r), ann.get_ha(), ann.get_va(), r)
        ann.set_position(best[1])
        ann.set_ha(best[2])
        ann.set_va(best[3])
        box = ann.get_window_extent(rend)
        taken.append(box.expanded(1.10, 1.45))
        if best[4] > 7.0:    # far enough that the reader needs the tie made explicit
            px, py = ax.transData.transform((x, y))
            ends = ax.transData.inverted().transform(
                [(px, py), (min(max(px, box.x0), box.x1), min(max(py, box.y0), box.y1))])
            ax.add_line(Line2D(ends[:, 0], ends[:, 1], color=tint(color, 0.45), lw=0.6,
                               zorder=5))


def leader_labels(ax, points, texts, color, x_frac, y_top, y_step, fontsize=7.5):
    """Name a tight cluster of points from a stacked column joined by thin leader lines.

    The points are too close together (a few points of vertical space apart) for adjacent
    labels, so the names are stacked in the empty part of the panel and each is tied to its
    own point. The stack keeps the points' own top-to-bottom order, which is what stops the
    leaders from crossing.
    """
    order = sorted(range(len(points)), key=lambda i: -points[i][1])
    for slot, i in enumerate(order):
        ax.annotate(
            texts[i], xy=points[i], xycoords="data",
            xytext=(x_frac, y_top - slot * y_step), textcoords="axes fraction",
            ha="left", va="center", fontsize=fontsize, color=color, annotation_clip=False,
            arrowprops=dict(arrowstyle="-", lw=0.6, color=tint(color, 0.45),
                            shrinkA=2.5, shrinkB=2.5),
        )


def main():
    use_paper_style()
    per_sys, disc, sweep = analyze()
    fig, (axA, axB, axC) = plt.subplots(
        1, 3, figsize=figsize(3, 2.75), gridspec_kw={"width_ratios": [1.06, 0.74, 1.06]},
    )

    # --- A: ranked per-system reduced chi^2 ---
    ps = sorted(per_sys, key=lambda d: d["rc"])
    ys = np.arange(len(ps), dtype=float)
    rc = np.array([d["rc"] for d in ps])
    fl = np.array([d["flag"] for d in ps])
    axA.axvspan(0, 1.0, color=tint(REF, 0.92), lw=0, zorder=0)
    reference_line(axA, "vline", 1.0)
    axA.scatter(rc[~fl], ys[~fl], c=MUTED, s=10, lw=0, alpha=0.8, zorder=3)
    axA.scatter(rc[fl], ys[fl], c=OURS, s=18, lw=0, zorder=4)
    axA.set_xscale("log")
    # right limit sits just past the largest chi2_nu; the name column below is stacked at the
    # axes edge and is allowed to run into the margin, which the layout accounts for. The old
    # 140 existed only to host that column and read as a decade of empty axis.
    axA.set_xlim(rc.min() * 0.55, rc.max() * 2.6)
    axA.set_ylim(-2.0, len(ps) + 1.0)
    # The six flagged systems occupy six adjacent ranks, a few points of height in all, so
    # their names are stacked in the empty right-hand third and joined by leaders.
    leader_labels(axA, [(d["rc"], y) for d, y in zip(ps, ys) if d["flag"]],
                  [d["sys"] for d in ps if d["flag"]], OURS,
                  x_frac=0.945, y_top=0.97, y_step=0.088)
    axA.set_xlabel(r"reduced $\chi^2_\nu$, replicate 0")
    axA.set_yticks([])
    axA.set_ylabel(f"{len(ps)} OpenFE systems (ranked)")
    panel(axA, "A", "most systems consistent")

    # --- B: fraction flagged vs sigma model ---
    frac = [100 * s[2] / s[3] for s in sweep]
    # x1.3 and x2 are inflations of the reported bar, i.e. variants of OUR quantity, not nulls;
    # they are tints of OURS. Grey is reserved for random/null and de-emphasised data.
    barcols = [FOIL, OURS, tint(OURS, 0.45), tint(OURS, 0.72)]
    axB.bar(range(len(sweep)), frac, width=0.66, color=barcols, zorder=2)
    notes = {0: ("stand-in", FOIL), 1: ("as shipped", OURS)}
    for i, (s, f) in enumerate(zip(sweep, frac)):
        axB.text(i, f + 2.5, f"{s[2]}/{s[3]}", ha="center", va="bottom", fontsize=7.5,
                 color=INK)
        if i in notes:
            axB.text(i, f + 10.0, notes[i][0], ha="center", va="bottom", fontsize=7,
                     color=notes[i][1])
    axB.set_xticks(range(len(sweep)))
    axB.set_xticklabels([f"\u00d7{s[1]:g}" for s in sweep])
    # bars are 0.66 wide on unit centres, so the gap between them is 0.34; give the
    # y-spine the same gap, or the four bars read as shoved to the right.
    axB.set_xlim(-0.67, len(sweep) - 0.33)
    axB.set_ylim(0, 126)
    axB.set_yticks([0, 25, 50, 75, 100])
    axB.spines["left"].set_bounds(0, 100)
    axB.tick_params(axis="x", length=0)
    axB.set_xlabel(r"$\sigma$ model")
    axB.set_ylabel("% of systems flagged")
    panel(axB, "B", r"flag rate vs $\sigma$")

    # --- C: single vs pooled-3 (systematic vs sampling) ---
    fl_sys = {d["sys"] for d in per_sys if d["flag"]}
    s1 = np.array([d["single"] for d in disc])
    s3 = np.array([d["pooled"] for d in disc])
    keep = np.array([d["sys"] in fl_sys for d in disc])
    lim = [0.02, max(s3.max(), s1.max()) * 2.0]
    # explicit dash patterns (not ":" / "--"): with the long legend handle below, the
    # swatch shows several periods of the pattern instead of one grey blob.
    axC.plot(lim, lim, color=REF, ls=(0, (1.0, 2.0)), lw=1.1, zorder=1,
             label="pooled = single\n(sampling)")
    axC.plot(lim, [3 * v for v in lim], color=REF, ls=(0, (4.5, 2.0)), lw=1.1, zorder=1,
             label="pooled = 3\u00d7 single\n(systematic)")
    axC.scatter(s1[~keep], s3[~keep], c=MUTED, s=12, lw=0, alpha=0.8, zorder=3)
    axC.scatter(s1[keep], s3[keep], c=OURS, s=20, lw=0, zorder=4)
    axC.set_xscale("log")
    axC.set_yscale("log")
    axC.set_xlim(*lim)
    axC.set_ylim(*lim)
    axC.set_box_aspect(1.0)
    axC.set_xlabel(r"reduced $\chi^2_\nu$, single replicate")
    axC.set_ylabel(r"reduced $\chi^2_\nu$, 3 replicates pooled")
    leg = legend(axC, loc="lower right", fontsize=7, labelspacing=0.8, handlelength=3.2,
                 borderaxespad=0.2, handletextpad=0.5)
    panel(axC, "C", "systematic, not sampling")

    fig.tight_layout(w_pad=1.9)
    fig.canvas.draw()
    place_labels(axC, list(zip(s1[keep], s3[keep])),
                 [d["sys"] for d, k in zip(disc, keep) if k], OURS,
                 obstacles=list(zip(s1, s3)),
                 keep_clear=[leg.get_window_extent(fig.canvas.get_renderer())])

    finish(fig, "figL_calibrated_cycle_closure")

    med = np.median([d["rc"] for d in per_sys])
    print(f"\n[A] {len(per_sys)} systems; median reduced chi2 = {med:.2f}; "
          f"flagged (BH .05): {sorted(d['sys'] for d in per_sys if d['flag'])}")
    print("[B] % systems flagged vs sigma-scale:")
    for n, sc, k, tot in sweep:
        print(f"    x{sc:<4} {k:>2}/{tot} ({100*k/tot:.0f}%)  [{n.split(chr(10))[0]}]")
    print("[C] systematic (pooled>1.5): "
          + ", ".join(f"{d['sys']}({d['single']:.1f}->{d['pooled']:.1f})"
                      for d in sorted(disc, key=lambda d: -d['pooled'])[:6]))


if __name__ == "__main__":
    main()
