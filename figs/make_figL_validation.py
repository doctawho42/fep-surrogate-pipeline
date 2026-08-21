"""Fig L (validation) -- removal counts, and cross-replicate reproduction of the residuals.

Panel A is an illustration of the repair procedure and is NOT evidence of causal localization:
the statistic driven down is sum_e z_e^2 and the removal order is descending z_e^2, so a world
with no localized systematic error reproduces the same contrast. The evidential panel is B.

Two data-internal validations of the calibrated cycle-closure detector (Fig L), short of new MD:

  A. REPAIR TEST (causal localization). For each flagged system, greedily remove the highest-|z|
     edge until reduced chi^2 <= 1. The detector's guided removals reach consistency in a handful
     of edges, whereas removing edges at random needs many more -- i.e. the flags point at the
     largest contributors to the statistic; this is arithmetic, not evidence (see the docstring).

  B. OUT-OF-SAMPLE REPRODUCTION (held-out replicates). Fit each of the 3 independent replicates'
     networks separately and compare per-edge standardized residuals. A systematic error
     reproduces across independent runs (positive residual correlation); sampling noise does not
     (correlation ~0). Measured r ~ +0.30..0.42 over 1143 edges -> the flags predict which edges
     misbehave in UNSEEN runs.

These are retrospective, data-internal predictive checks; the decisive prospective test (re-run /
repair a flagged edge with fresh MD and confirm the cycle closes) is left as the natural next
experiment. Run: PYTHONPATH=src python figs/make_figL_validation.py  (or `make figLval`).
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

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from bar.qc import gls_network, repair_order  # noqa: E402
from paperstyle import (  # noqa: E402
    FULL, INK, MUTED, OURS, finish, legend, panel, reference_line, use_paper_style,
)

DATA = ROOT / "data" / "openfe_replicates" / "combined_pymbar4_edge_data.csv"
FIGDIR = pathlib.Path(__file__).resolve().parent
FLAGGED = ["brd4", "bace", "faah", "cdk8", "hif2a", "p38"]
# repair test is only meaningful with several cycles; brd4 has just one (removing any cycle-edge
# trivially breaks it), so the guided-vs-random contrast is shown on the multi-cycle systems.
REPAIR = ["bace", "faah", "cdk8", "hif2a", "p38"]
# Semantic colours (paperstyle): OURS = the calibrated detector and anything computed from it --
# its flags and its |z|-guided removal rule, in BOTH panels. MUTED = the random-removal null
# baseline, and edges that are present but not the point. No FOIL here: nothing in this figure
# is the overconfident stand-in.


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


def load():
    by = defaultdict(list)
    for r in csv.DictReader(open(DATA)):
        by[r["system name"]].append(r)
    return by


def reduced_chi2(edges):
    if len({a for a, *_ in edges} | {b for _, b, *_ in edges}) < 2:
        return 0.0
    fit = gls_network(edges)
    return fit.reduced_chi2 if fit.dof > 0 else 0.0


def repair_test(by):
    rng = np.random.default_rng(0)
    out = []
    for sysname in REPAIR:
        e0 = [(r["ligand_A"], r["ligand_B"], *edge_val(r, 0)) for r in by[sysname] if edge_val(r, 0)]
        removed, _ = repair_order(e0, target_reduced_chi2=1.0)
        guided = len(removed)
        rand = []
        for _ in range(40):
            ee = list(e0); c = 0
            while reduced_chi2(ee) > 1.0 and len(ee) > 1:
                ee.pop(int(rng.integers(len(ee)))); c += 1
            rand.append(c)
        out.append((sysname, len(e0), guided, int(np.median(rand)), int(np.max(rand))))
    return out


def out_of_sample(by):
    Z = {0: {}, 1: {}, 2: {}}
    flagged_edges = set()
    for sysname, rs in by.items():
        for k in (0, 1, 2):
            ek = [(r["ligand_A"], r["ligand_B"], *edge_val(r, k)) for r in rs if edge_val(r, k)]
            if len(ek) < 3:
                continue
            fit = gls_network(ek)
            if fit.dof < 1:
                continue
            for (a, b, _, _), zz in zip(ek, fit.z):
                Z[k][(sysname, a, b)] = zz
    keys = [k for k in Z[0] if k in Z[1] and k in Z[2]]
    z0 = np.array([Z[0][k] for k in keys]); z1 = np.array([Z[1][k] for k in keys])
    z2 = np.array([Z[2][k] for k in keys])
    corrs = (np.corrcoef(z0, z1)[0, 1], np.corrcoef(z0, z2)[0, 1], np.corrcoef(z1, z2)[0, 1])
    for k in keys:
        if k[0] in FLAGGED:
            flagged_edges.add(k)
    fl = np.array([k[0] in FLAGGED for k in keys])
    return z0, z1, z2, corrs, fl


def main():
    use_paper_style()
    by = load()
    rep = repair_test(by)
    z0, z1, z2, corrs, fl = out_of_sample(by)

    fig, (axA, axB) = plt.subplots(
        1, 2, figsize=(FULL, 3.0), constrained_layout=True,
        gridspec_kw={"width_ratios": [1.42, 1.0]},
    )

    # A: repair test
    names = [r[0] for r in rep]
    ypos = np.arange(len(rep), dtype=float)
    guided = [r[2] for r in rep]
    randmed = [r[3] for r in rep]
    axA.barh(ypos + 0.19, guided, height=0.33, color=OURS, label="guided by |z| (the flags)")
    axA.barh(ypos - 0.19, randmed, height=0.33, color=MUTED, label="random removal (median)")
    for y, r in zip(ypos, rep):
        axA.text(r[2] + 0.7, y + 0.19, str(r[2]), va="center", fontsize=7.5, color=OURS)
        axA.text(r[3] + 0.7, y - 0.19, str(r[3]), va="center", fontsize=7.5, color=MUTED)
    axA.set_yticks(ypos)
    axA.set_yticklabels(names)
    axA.tick_params(axis="y", length=0)
    axA.set_ylim(-0.62, len(rep) - 1 + 1.18)
    axA.set_xlim(0, max(max(guided), max(randmed)) * 1.10)
    axA.set_xlabel(r"edges removed to reach reduced $\chi^2_\nu \leq 1$")
    legend(axA, loc="upper right", handlelength=1.1, handleheight=0.9, borderaxespad=0.2)
    panel(axA, "A", "guided removal versus random")

    # B: out-of-sample residual reproduction
    axB.scatter(z0[~fl], z1[~fl], s=7, c=MUTED, alpha=0.55, lw=0, zorder=2,
                label="edge, other systems")
    axB.scatter(z0[fl], z1[fl], s=10, c=OURS, alpha=0.7, lw=0, zorder=3,
                label="edge, flagged systems")
    lim = 6
    reference_line(axB, "diagonal")
    axB.set_xlim(-lim, lim)
    axB.set_ylim(-lim, lim)
    axB.set_xticks(np.arange(-lim, lim + 1, 2))
    axB.set_yticks(np.arange(-lim, lim + 1, 2))
    axB.set_box_aspect(1)
    axB.set_xlabel(r"standardized residual $z_e$, replicate 0")
    axB.set_ylabel(r"$z_e$, held-out replicate 1")
    # The annotation carries NO patch behind it: it is placed in a corner that the data leave
    # empty. |z_e| < 4.3 everywhere, and no edge has z_0 < -0.9 together with z_1 > 3.2, so the
    # block's four lines clear every marker; the nearest are (-0.54, +3.15) and (-3.21, +2.86).
    # A white box over the scatter is what this replaces -- never put a patch over data.
    axB.text(0.03, 0.995,
             f"r(0,1) = {corrs[0]:+.2f}\n"
             f"r(0,2) = {corrs[1]:+.2f}\n"
             f"r(1,2) = {corrs[2]:+.2f}\n"
             f"sampling noise \u2192 0",
             transform=axB.transAxes, va="top", ha="left", fontsize=7.5, linespacing=1.2,
             color=INK, zorder=5)
    legend(axB, loc="lower right", markerscale=1.6, handletextpad=0.4, borderaxespad=0.2)
    panel(axB, "B", "residuals persist on held-out runs")

    finish(fig, "figL_validation")
    print("\n[A] repair test (guided vs random edges to reach reduced chi2<=1):")
    for s, E, g, rm, rx in rep:
        print(f"    {s:8s}: guided {g:>2d} vs random median {rm:>2d} (max {rx:>2d}) of {E} edges")
    print(f"[B] out-of-sample per-edge residual correlation across replicates: "
          f"r(0,1)={corrs[0]:+.2f}, r(0,2)={corrs[1]:+.2f}, r(1,2)={corrs[2]:+.2f}  (n={z0.size})")


if __name__ == "__main__":
    main()
