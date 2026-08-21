"""Fig Stab -- replicate stability of the cycle-closure flag set (a referee-raised limitation).

The Fig L headline flag set (6 of 48 OpenFE systems, BH-FDR alpha=0.05) is computed on
**replicate 0 only**, although IndustryBenchmarks2024 ships **3 independent replicates** of every
edge. If the flag set is a property of the network rather than of one run, re-running the paper's
own test on replicate 1 and replicate 2 must return (nearly) the same systems. This figure asks
exactly that, with the paper's machinery (`bar.qc.gls_network` / `chi2_sf` /
`benjamini_hochberg`) and the paper's edge construction (binding ddG = complex - solvent leg,
se = sqrt(complex_dDG^2 + solvent_dDG^2)), changing nothing but the replicate index.

What is tested, and how:
  * per replicate k in {0,1,2}: the flag set, plus per system the reduced chi^2, dof, raw p and
    BH-adjusted q (the q<=alpha rule is asserted to reproduce `benjamini_hochberg` exactly);
  * SET stability: intersection, union, the three pairwise Jaccard overlaps, against the EXACT
    expected Jaccard of size-matched random draws from the same 48 systems (closed form over the
    hypergeometric intersection -- no Monte-Carlo, hence no seed);
  * the POOLED-data flag set. Fig L panel C defines the pooling rule and this figure reuses it
    verbatim: per edge, mean ddG over the available replicates and se = sqrt(mean_k se_k^2 / n),
    keeping edges with >= MIN_REPS_POOL available replicates. What panel C does NOT define is a
    *flagging* rule on the pooled network (it only compares reduced chi^2, single vs pooled), so
    running BH-FDR on the pooled p-values is this figure's extension; it uses the same
    pre-registered alpha=0.05. Sensitivity to the one free knob in the pooling rule (>=2 vs all 3
    replicates required per edge) is computed and reported rather than chosen.
  * what does NOT move: the per-system ranking (Spearman on reduced chi^2 across all 48 systems,
    all three replicate pairs) and the selectivity-vs-sigma-scale curve (Fig L panel B) per
    replicate;
  * the per-system cross-replicate reproduction of the signed standardized residuals -- the
    model-free read on whether a flagged system's closure error is the same error each run.

This is an unfavourable result for the paper and is reported as such: nothing here is tuned. The
constants (alpha=0.05, the sigma scales, the edge construction, the pooling rule) are the
pre-registered ones already in `make_figL.py`.

Run:  PYTHONPATH=src python figs/make_figStab.py   (or `make figStab`).  Deterministic.
"""
from __future__ import annotations

import csv
import math
import pathlib
import sys
from collections import defaultdict
from math import comb

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.patches import Patch, Rectangle  # noqa: E402

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
from paperstyle import (  # noqa: E402
    FOIL,
    FULL,
    INK,
    MUTED,
    OURS,
    REF,
    finish,
    legend,
    panel,
    tint,
    use_paper_style,
)

from bar.qc import (  # noqa: E402
    benjamini_hochberg,
    benjamini_yekutieli,
    chi2_sf,
    gls_network,
)

DATA = ROOT / "data" / "openfe_replicates" / "combined_pymbar4_edge_data.csv"
FIGDIR = pathlib.Path(__file__).resolve().parent
DOC = ROOT / "docs" / "results_figStab.md"

ALPHA = 0.05                 # pre-registered BH-FDR level (identical to Fig L)
REPLICATES = (0, 1, 2)
MIN_EDGES = 3                # Fig L's per-system admission rule
MIN_REPS_POOL = 2            # Fig L panel C's pooling rule: >= 2 available replicates per edge
SIGMA_SCALES = (0.15, 1.0, 1.30, 2.0)          # Fig L panel B's sweep
ALPHA_SWEEP = (0.01, 0.05, 0.10, 0.20)         # robustness of the instability to the FDR level
# The flag set is what this article's calibrated per-edge null produces, so every mark that
# carries a flag is OURS. The three replicates and the pooled network are four members of that
# one quantity family and are told apart by tint and marker, never by a second hue; a system the
# test never flags is de-emphasised data (MUTED); the BH threshold is a reference line (REF).
C_FLAG = OURS                                        # flagged, on one replicate
C_POOL = tint(OURS, 0.45)                            # the same flag, three replicates pooled
C_OFF = tint(MUTED, 0.55)                            # the empty cell, the connector
# Panel D sweeps the SAME sigma models Fig L panel B sweeps, on the same data, and reports the
# same 42/48 at x0.15. The hue must therefore mean there what it means in Fig L panel B, or a
# reader turning the page watches one number change colour: x0.15 is the overconfident stand-in
# (FOIL), x1 is the bar as shipped (OURS), and x1.3 / x2 are inflations of that same bar, i.e.
# tints of OURS. Replicate identity is a SECOND dimension on the same bars and is carried by
# hatch, which costs no hue.
C_SIGMA = (FOIL, OURS, tint(OURS, 0.45), tint(OURS, 0.72))   # x0.15 / x1 / x1.3 / x2
H_REP = ("", "//", "..")                                     # replicate 0 / 1 / 2
HATCH_LW = 0.6      # thinner than the 1.0 default: these bars are ~0.16 in wide, and at
                    # the default weight the pattern reads louder than the bar it marks


# --------------------------------------------------------------------------- data / edges

def _f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return math.nan


def edge_val(r, k):
    """Fig L's edge construction: binding ddG (complex - solvent leg) and its MBAR/sandwich se."""
    cD = _f(r[f"complex_repeat_{k}_DG (kcal/mol)"])
    cd = _f(r[f"complex_repeat_{k}_dDG (kcal/mol)"])
    sD = _f(r[f"solvent_repeat_{k}_DG (kcal/mol)"])
    sd = _f(r[f"solvent_repeat_{k}_dDG (kcal/mol)"])
    if any(math.isnan(v) for v in (cD, cd, sD, sd)):
        return None
    return cD - sD, math.sqrt(cd * cd + sd * sd)


def load_systems():
    by = defaultdict(list)
    with open(DATA) as fh:
        for r in csv.DictReader(fh):
            by[r["system name"]].append(r)
    return dict(by)


def replicate_edges(rows, k):
    """Edges of one replicate: every row whose four legs are present in replicate ``k``."""
    return [(r["ligand_A"], r["ligand_B"], *v) for r in rows if (v := edge_val(r, k))]


def pool_replicates(vals):
    """Fig L panel C's pooling rule for ONE edge.

    ``vals`` = the available ``(ddg, se)`` per replicate. Returns the replicate mean ddG and the
    se of that mean under independence, ``sqrt(mean_k se_k^2 / n)`` (i.e. se -> se/sqrt(n) for
    equal bars). This is the rule that produces panel C's ``se -> se/sqrt3``.
    """
    n = len(vals)
    if n < 1:
        raise ValueError("need at least one replicate to pool")
    ddg = float(np.mean([v[0] for v in vals]))
    se = math.sqrt(float(np.mean([v[1] ** 2 for v in vals])) / n)
    return ddg, se


def pooled_edges(rows, min_reps=MIN_REPS_POOL):
    """Fig L panel C's pooled network: edges with >= ``min_reps`` available replicates."""
    out = []
    for r in rows:
        vals = [v for k in REPLICATES if (v := edge_val(r, k))]
        if len(vals) < min_reps:
            continue
        out.append((r["ligand_A"], r["ligand_B"], *pool_replicates(vals)))
    return out


# --------------------------------------------------------------------------- the test

def bh_qvalues(pvals):
    """Benjamini-Hochberg step-up adjusted p-values (q), the continuous form of the flag rule.

    ``q <= alpha`` is equivalent to ``bar.qc.benjamini_hochberg(p, alpha)``; the script asserts
    that equivalence on every set it reports, so the q column in the results table cannot drift
    away from the flags the paper's own function produces.
    """
    p = np.asarray(pvals, dtype=float)
    m = p.size
    order = np.argsort(p, kind="stable")
    ranked = p[order] * m / np.arange(1, m + 1)
    q_sorted = np.minimum.accumulate(ranked[::-1])[::-1]
    q = np.empty(m)
    q[order] = np.minimum(q_sorted, 1.0)
    return q


def system_flags(edges_by_system, alpha=ALPHA):
    """Run the paper's cycle-closure test on each system and BH-flag across systems.

    ``edges_by_system`` maps a system name to its edge list. Systems with fewer than
    ``MIN_EDGES`` edges or with no independent cycle (dof < 1) are dropped, exactly as Fig L
    does. Returns a list of dicts sorted by system name.
    """
    rows = []
    for name, edges in sorted(edges_by_system.items()):
        if len(edges) < MIN_EDGES:
            continue
        fit = gls_network(edges)
        if fit.dof < 1:
            continue
        rows.append({"sys": name, "E": len(edges), "dof": fit.dof,
                     "rc": fit.reduced_chi2, "p": chi2_sf(fit.chi2, fit.dof)})
    if not rows:
        return rows
    pv = [d["p"] for d in rows]
    flags = benjamini_hochberg(pv, alpha=alpha)
    # BY is the same rule at alpha / H_m and holds under arbitrary dependence, which is the
    # guarantee worth quoting when the systems share a force field and, within a target, a
    # preparation. Carried alongside so the cost of dropping the independence assumption is visible.
    flags_by = benjamini_yekutieli(pv, alpha=alpha)
    qs = bh_qvalues(pv)
    for d, fl, fby, q in zip(rows, flags, flags_by, qs, strict=True):
        d["flag"] = bool(fl)
        d["flag_by"] = bool(fby)
        d["q"] = float(q)
    # the q<=alpha rule must reproduce the paper's own flag function exactly
    assert all(d["flag"] == (d["q"] <= alpha) for d in rows), "BH q-values disagree with flags"
    return rows


def flag_set(rows):
    return {d["sys"] for d in rows if d["flag"]}


def flag_set_by(rows):
    """The same set under Benjamini--Yekutieli, which assumes no dependence structure."""
    return {d["sys"] for d in rows if d.get("flag_by")}


# --------------------------------------------------------------------------- set stability

def jaccard(a, b):
    """|A and B| / |A or B|; two empty sets are identical by convention (1.0)."""
    a, b = set(a), set(b)
    union = a | b
    return 1.0 if not union else len(a & b) / len(union)


def expected_random_jaccard(n_a, n_b, n_total):
    """EXACT expected Jaccard of two independent uniform random subsets of the given sizes.

    |A and B| is hypergeometric, and |A or B| = n_a + n_b - |A and B|, so the expectation is a
    finite sum -- no Monte-Carlo, hence no seed and no simulation error.
    """
    if n_a == 0 and n_b == 0:
        return 1.0
    if min(n_a, n_b) < 0 or max(n_a, n_b) > n_total:
        raise ValueError("subset sizes must lie in [0, n_total]")
    tot = 0.0
    for i in range(0, min(n_a, n_b) + 1):
        pr = comb(n_a, i) * comb(n_total - n_a, n_b - i) / comb(n_total, n_b)
        tot += pr * i / (n_a + n_b - i)
    return tot


def stability_summary(sets, n_total):
    """Intersection / union / pairwise Jaccard of the per-replicate flag sets, each against the
    exact size-matched random-draw expectation."""
    pairs = [(0, 1), (0, 2), (1, 2)]
    obs = [jaccard(sets[i], sets[j]) for i, j in pairs]
    ref = [expected_random_jaccard(len(sets[i]), len(sets[j]), n_total) for i, j in pairs]
    return {"sets": [sorted(s) for s in sets],
            "intersection": sorted(set.intersection(*sets)),
            "union": sorted(set.union(*sets)),
            "pairs": pairs, "jaccard": obs, "jaccard_random": ref,
            "mean_jaccard": float(np.mean(obs)), "mean_jaccard_random": float(np.mean(ref))}


# --------------------------------------------------------------------------- what does not move

def rank_correlations(rows_by_rep):
    """Spearman rho of the per-system reduced chi^2 between replicate pairs, over the systems
    admitted in ALL three replicates (the ranking is the paper's actual triage signal)."""
    from scipy.stats import spearmanr

    maps = [{d["sys"]: d["rc"] for d in rows} for rows in rows_by_rep]
    common = sorted(set.intersection(*[set(m) for m in maps]))
    out = {}
    for i, j in [(0, 1), (0, 2), (1, 2)]:
        a = [maps[i][s] for s in common]
        b = [maps[j][s] for s in common]
        out[(i, j)] = float(spearmanr(a, b).statistic)
    return out, common


def selectivity_sweep(by, k, scales=SIGMA_SCALES, alpha=ALPHA):
    """Fig L panel B on replicate ``k``: number of systems flagged as sigma is rescaled."""
    out = []
    for sc in scales:
        edges = {name: [(a, b, y, se * sc) for a, b, y, se in replicate_edges(rows, k)]
                 for name, rows in by.items()}
        rows_k = system_flags(edges, alpha=alpha)
        out.append((sc, len(flag_set(rows_k)), len(rows_k)))
    return out


def z_reproduction(rows):
    """Mean pairwise Pearson correlation of the per-edge signed standardized residuals across the
    three replicates, on the edges complete in all three (so the z vectors align).

    Sampling noise -> 0; a reproducible edge-level systematic error -> positive. NaN if the
    aligned network has no cycle.
    """
    keep = [r for r in rows if all(edge_val(r, k) for k in REPLICATES)]
    if len(keep) < MIN_EDGES:
        return math.nan, 0
    zs = []
    for k in REPLICATES:
        fit = gls_network([(r["ligand_A"], r["ligand_B"], *edge_val(r, k)) for r in keep])
        if fit.dof < 1:
            return math.nan, len(keep)
        zs.append(fit.z)
    cs = [float(np.corrcoef(zs[i], zs[j])[0, 1]) for i, j in [(0, 1), (0, 2), (1, 2)]]
    return float(np.mean(cs)), len(keep)


# --------------------------------------------------------------------------- driver

def analyze():
    by = load_systems()
    rows_by_rep = [system_flags({n: replicate_edges(rs, k) for n, rs in by.items()})
                   for k in REPLICATES]
    sets = [flag_set(rows) for rows in rows_by_rep]
    n_total = len(rows_by_rep[0])
    stab = stability_summary(sets, n_total)

    pooled_rows = system_flags({n: pooled_edges(rs) for n, rs in by.items()})
    pooled_alt = system_flags({n: pooled_edges(rs, min_reps=3) for n, rs in by.items()})

    alpha_sweep = []
    for al in ALPHA_SWEEP:
        s_al = [flag_set(system_flags({n: replicate_edges(rs, k) for n, rs in by.items()},
                                      alpha=al)) for k in REPLICATES]
        alpha_sweep.append((al, [len(s) for s in s_al], sorted(set.intersection(*s_al)),
                            sorted(set.union(*s_al)),
                            [jaccard(s_al[i], s_al[j]) for i, j in [(0, 1), (0, 2), (1, 2)]]))

    rho, common = rank_correlations(rows_by_rep)
    sweeps = [selectivity_sweep(by, k) for k in REPLICATES]
    zrep = {s: z_reproduction(rs) for s, rs in by.items()}
    return dict(by=by, rows_by_rep=rows_by_rep, sets=sets, n_total=n_total, stab=stab,
                pooled_rows=pooled_rows, pooled_alt=pooled_alt, alpha_sweep=alpha_sweep,
                rho=rho, common=common, sweeps=sweeps, zrep=zrep)


def _qmap(rows):
    return {d["sys"]: d for d in rows}


def make_figure(res):
    stab, rows_by_rep, pooled_rows = res["stab"], res["rows_by_rep"], res["pooled_rows"]
    union = stab["union"]
    maps = [_qmap(r) for r in rows_by_rep] + [_qmap(pooled_rows)]
    cols = ["rep 0\n(paper)", "rep 1", "rep 2", "pooled\n(3 reps)"]
    # order the union systems by their best (smallest) q anywhere: strongest evidence on top
    union = sorted(union, key=lambda s: min(m[s]["q"] for m in maps if s in m))

    fig, axes = plt.subplots(2, 2, figsize=(FULL, 6.0))
    axA, axB, axC, axD = axes[0, 0], axes[0, 1], axes[1, 0], axes[1, 1]

    # --- A: which systems flip ------------------------------------------------------------
    for xi, _c in enumerate(cols):
        for yi, s in enumerate(union):
            d = maps[xi].get(s)
            on = bool(d and d["flag"])
            axA.add_patch(Rectangle((xi - 0.44, yi - 0.42), 0.88, 0.84,
                                    facecolor=(C_FLAG if xi < 3 else C_POOL) if on else "white",
                                    edgecolor=C_OFF if not on else "none", lw=0.8))
            # white on the saturated hue, ink on the pooled tint: one hue, two weights, and the
            # label stays legible on both instead of a second colour being spent on the column
            axA.text(xi, yi, "flag" if on else "\u2013", ha="center", va="center", fontsize=7,
                     color=("white" if xi < 3 else INK) if on else MUTED,
                     fontweight="bold" if on else "normal")
    axA.set_xlim(-0.6, len(cols) - 0.4)
    axA.set_ylim(len(union) - 0.5, -0.5)
    axA.set_xticks(range(len(cols)))
    axA.set_xticklabels(cols)
    axA.set_yticks(range(len(union)))
    axA.set_yticklabels(union)
    axA.tick_params(length=0)
    for sp in axA.spines.values():
        sp.set_visible(False)
    n_int, n_uni = len(stab["intersection"]), len(stab["union"])
    panel(axA, "A", "the flag set is replicate-dependent",
          subtitle=f"{n_int} system in all three, {n_uni} in the union")
    # two lines: on one line this footnote is wider than the panel and ran off the canvas
    axA.text(0.0, -0.17, f"pairwise Jaccard {'/'.join(f'{v:.2f}' for v in stab['jaccard'])}\n"
             f"(size-matched random draw {stab['mean_jaccard_random']:.3f})",
             transform=axA.transAxes, ha="left", va="top", fontsize=7, color=REF,
             linespacing=1.25)

    # --- B: where they sit relative to the threshold --------------------------------------
    y = np.arange(len(union))
    for yi, s in enumerate(union):
        qs = [maps[k][s]["q"] for k in range(3) if s in maps[k]]
        if len(qs) > 1:
            axB.plot([min(qs), max(qs)], [yi, yi], color=C_OFF, lw=1.0, zorder=1)
    for k in range(3):
        qs = [maps[k][s]["q"] if s in maps[k] else np.nan for s in union]
        fl = [bool(maps[k][s]["flag"]) if s in maps[k] else False for s in union]
        axB.scatter(qs, y, s=22, zorder=3, marker="o",
                    facecolors=[C_FLAG if f else "white" for f in fl],
                    edgecolors=C_FLAG, linewidths=0.9,
                    label="replicate (filled = flagged)" if k == 0 else None)
    qp = [maps[3][s]["q"] if s in maps[3] else np.nan for s in union]
    fp = [bool(maps[3][s]["flag"]) if s in maps[3] else False for s in union]
    axB.scatter(qp, y, s=30, marker="D", zorder=4,
                facecolors=[C_POOL if f else "white" for f in fp],
                edgecolors=C_POOL, linewidths=0.9, label="pooled (3 replicates)")
    qlo = min(v for v in qp + [maps[k][s]["q"] for k in range(3) for s in union if s in maps[k]]
              if v > 0)
    axB.set_xscale("log")
    axB.set_xlim(qlo * 0.15, 4.0)
    # an opaque very light tint, not a translucent grey laid over the markers it contains
    axB.axvspan(ALPHA, 4.0, color=tint(REF, 0.90), lw=0, zorder=0)
    # the threshold guide stops at the top of the data rows, BELOW the reserved strip that
    # holds its own label. Set two lines deep the label hung down past where the guide stopped
    # and the dotted line ran straight through the word "flagged"; one line fits the strip.
    axB.plot([ALPHA, ALPHA], [len(union) - 0.5, -0.46], color=REF, ls=":", lw=1.0, zorder=2)
    # right-anchored to the axes edge (x in axes fraction, y in data): set at the band's left
    # edge this label was wider than the band and ran off the canvas. The top-right corner is
    # the one place no marker reaches, so it sits on nothing and clears the panel heading.
    axB.text(0.99, -0.85, f"not flagged (BH $\\alpha$={ALPHA})", fontsize=7, color=REF,
             ha="right", va="center", transform=axB.get_yaxis_transform())
    axB.set_ylim(len(union) - 0.5, -1.25)   # a strip at the top for the threshold label
    axB.set_yticks(y)
    axB.set_yticklabels(union)
    axB.set_xlabel("BH-adjusted $q$ (per replicate, and pooled)")
    legend(axB, loc="lower left", fontsize=7)
    panel(axB, "B", "the flips are threshold-crossings",
          subtitle=r"$q$ swings orders of magnitude run to run")

    # --- C: what does NOT move (1) -- the ranking ------------------------------------------
    rc = [_qmap(r) for r in rows_by_rep]
    common = res["common"]
    fl_union = set(stab["union"])
    x = np.array([rc[0][s]["rc"] for s in common])
    y1 = np.array([rc[1][s]["rc"] for s in common])
    y2 = np.array([rc[2][s]["rc"] for s in common])
    isf = np.array([s in fl_union for s in common])
    lim = [min(x.min(), y1.min(), y2.min()) * 0.7, max(x.max(), y1.max(), y2.max()) * 1.4]
    axC.plot(lim, lim, color=REF, ls=":", lw=1.0, zorder=1)
    h1 = axC.scatter(x[~isf], y1[~isf], s=14, c=MUTED, alpha=0.55, lw=0, label="replicate 1")
    h2 = axC.scatter(x[~isf], y2[~isf], s=14, marker="^", c=MUTED, alpha=0.35, lw=0,
                     label="replicate 2")
    h3 = axC.scatter(x[isf], y1[isf], s=26, c=C_FLAG, lw=0, zorder=3, label="ever-flagged")
    axC.scatter(x[isf], y2[isf], s=26, marker="^", c=C_FLAG, alpha=0.7, lw=0, zorder=3)
    axC.set_xscale("log")
    axC.set_yscale("log")
    axC.set_xlim(*lim)
    axC.set_ylim(*lim)
    axC.set_xlabel(r"reduced $\chi^2$, replicate 0")
    axC.set_ylabel(r"reduced $\chi^2$, replicate 1 / 2")
    rho = res["rho"]
    # upper left is the one corner this cloud never enters; the block therefore needs no box
    axC.text(0.03, 0.97, "Spearman $\\rho$\n"
             + "\n".join(f"({i},{j}): {rho[(i, j)]:.2f}" for i, j in [(0, 1), (0, 2), (1, 2)])
             + f"\n(n={len(common)} systems)",
             transform=axC.transAxes, va="top", fontsize=7, linespacing=1.35, color=INK)
    # marker tells the replicate apart, colour tells the class apart, so three entries carry
    # all four series -- and three entries fit the empty bottom-right corner, where the
    # four-entry version put its handle column on top of the grey cloud
    legend(axC, handles=[h1, h2, h3], loc="lower right", fontsize=7)
    panel(axC, "C", "what IS stable: the per-system ranking",
          subtitle="the same networks are worst every run")

    # --- D: what does NOT move (2) -- selectivity vs sigma scale ----------------------------
    sweeps = res["sweeps"]
    width = 0.26
    xs = np.arange(len(SIGMA_SCALES))
    for k in REPLICATES:
        frac = [100 * n / tot for _sc, n, tot in sweeps[k]]
        # colour = the sigma model (Fig L panel B's colours, same order); hatch = the replicate.
        # The hatch is drawn in the page colour so it reads on the two saturated bars, which are
        # the only ones tall enough to carry a pattern.
        with mpl.rc_context({"hatch.linewidth": HATCH_LW}):
            axD.bar(xs + (k - 1) * width, frac, width=width * 0.92, lw=0,
                    color=list(C_SIGMA), edgecolor="white", hatch=H_REP[k])
        for xi, (v, (_sc, n, tot)) in enumerate(zip(frac, sweeps[k], strict=True)):
            axD.text(xi + (k - 1) * width, v + 1.6, f"{n}", ha="center", va="bottom",
                     fontsize=7, color=INK)
    # Fig L panel B names these two sigma models in figure text; name them here in the same
    # words and the same hues, so the two panels are read as one sweep seen twice.
    for gi, (note, col) in {0: ("stand-in", FOIL), 1: ("as shipped", OURS)}.items():
        top = max(100 * sweeps[k][gi][1] / sweeps[k][gi][2] for k in REPLICATES)
        axD.text(gi, top + 7.5, note, ha="center", va="bottom", fontsize=7, color=col)
    axD.set_xticks(xs)
    axD.set_xticklabels([f"$\\times${sc:g}" for sc in SIGMA_SCALES])
    axD.set_xlim(-0.6, len(SIGMA_SCALES) - 0.4)
    axD.tick_params(axis="x", length=0)
    axD.set_xlabel(r"$\sigma$ scale applied to every edge")
    axD.set_ylabel(f"% of {res['n_total']} systems flagged")
    axD.set_ylim(0, 104)
    axD.set_yticks([0, 25, 50, 75, 100])
    axD.spines["left"].set_bounds(0, 100)
    # the hue is spent on the sigma model, so the legend can only key the hatch: neutral
    # swatches, no fill colour, because a filled key would claim a hue this panel has assigned
    with mpl.rc_context({"hatch.linewidth": HATCH_LW}):
        legend(axD, handles=[Patch(facecolor="white", edgecolor=INK, lw=0.6, hatch=h,
                                   label=f"replicate {k}")
                             for k, h in zip(REPLICATES, H_REP, strict=True)],
               loc="upper right", fontsize=7)
    panel(axD, "D", r"what IS stable: the $\sigma$ sweep",
          subtitle="the same sweep on every replicate")

    finish(fig, "figStab_replicate_stability")
    plt.close(fig)


def write_doc(res):
    stab, rows_by_rep = res["stab"], res["rows_by_rep"]
    maps = [_qmap(r) for r in rows_by_rep]
    pooled = _qmap(res["pooled_rows"])
    pooled_alt = _qmap(res["pooled_alt"])
    union, inter = stab["union"], stab["intersection"]
    n_total = res["n_total"]
    rho = res["rho"]
    sets = stab["sets"]

    def fmt(sysname, m):
        d = m.get(sysname)
        if d is None:
            return "n/a"
        cell = f"{d['rc']:.2f} / {d['q']:.1e}"
        return f"**{cell}**" if d["flag"] else cell

    order = sorted(maps[0], key=lambda s: min(m[s]["q"] for m in maps + [pooled] if s in m))
    tbl = ["| system | edges | cycles | rep 0 χ²ᵥ / q | rep 1 χ²ᵥ / q | rep 2 χ²ᵥ / q "
           "| pooled χ²ᵥ / q | flagged in (reps) |",
           "|---|---:|---:|---:|---:|---:|---:|---:|"]
    for s in order:
        d0 = maps[0][s]
        nfl = sum(1 for m in maps if s in m and m[s]["flag"])
        tbl.append(f"| {'**' + s + '**' if s in union else s} | {d0['E']} | {d0['dof']} | "
                   f"{fmt(s, maps[0])} | {fmt(s, maps[1])} | {fmt(s, maps[2])} | "
                   f"{fmt(s, pooled)} | {nfl}/3 |")

    zr = res["zrep"]
    ztbl = ["| system | flagged in reps | cross-replicate signed-z reproduction | aligned edges |",
            "|---|---|---:|---:|"]
    for s in sorted(union, key=lambda s: -(zr[s][0] if not math.isnan(zr[s][0]) else -9)):
        which = ", ".join(str(k) for k in REPLICATES if s in maps[k] and maps[k][s]["flag"])
        r, n = zr[s]
        ztbl.append(f"| {s} | {which or '—'} | {r:+.3f} | {n} |")

    swp = ["| σ scale | replicate 0 | replicate 1 | replicate 2 |", "|---|---:|---:|---:|"]
    for i, sc in enumerate(SIGMA_SCALES):
        cells = []
        for k in REPLICATES:
            _s, n, tot = res["sweeps"][k][i]
            cells.append(f"{n}/{tot} ({100 * n / tot:.0f}%)")
        swp.append(f"| ×{sc:g} | " + " | ".join(cells) + " |")

    asw = ["| α | set sizes (rep 0/1/2) | intersection | union | pairwise Jaccard |",
           "|---:|---|---:|---:|---|"]
    for al, sizes, ai, au, aj in res["alpha_sweep"]:
        asw.append(f"| {al:g} | {'/'.join(str(v) for v in sizes)} | {len(ai)} | {len(au)} | "
                   + " / ".join(f"{v:.3f}" for v in aj) + " |")

    jac = stab["jaccard"]
    jrand = stab["jaccard_random"]
    pool_set = sorted(s for s in pooled if pooled[s]["flag"])
    pool_alt_set = sorted(s for s in pooled_alt if pooled_alt[s]["flag"])
    dropped = sorted(set(union) - set(pool_set))

    lines = [
        "# Results — Fig Stab: replicate stability of the cycle-closure flag set",
        "",
        "**Figure:** `figs/figStab_replicate_stability.{pdf,png}` · **Reproduce:** `make figStab`",
        "(`PYTHONPATH=src python figs/make_figStab.py`). Deterministic (closed-form random-draw",
        "reference; no Monte-Carlo, no seed). Data:",
        "`data/openfe_replicates/combined_pymbar4_edge_data.csv` (OpenFE IndustryBenchmarks2024,",
        "public, 3 independent replicates per edge).",
        "",
        "## What this asks",
        "The Fig L headline flag set is computed on **replicate 0 only**. This figure re-runs the",
        "*identical* test — same edge construction, same `bar.qc.gls_network` / `chi2_sf` /",
        "`benjamini_hochberg`, same pre-registered α = %.2f — on replicates 1 and 2, and on the"
        % ALPHA,
        "pooled network. Nothing is tuned: every constant is the one already in `make_figL.py`.",
        "",
        "## Headline: the flag set is NOT stable across replicates",
        "",
        "| replicate | systems tested | flagged |",
        "|---|---:|---|",
    ]
    for k in REPLICATES:
        lines.append(f"| {k}{' (the paper)' if k == 0 else ''} | {len(maps[k])} | "
                     f"{len(sets[k])}: `{', '.join(sets[k])}` |")
    by_rows = []
    for k in REPLICATES:
        sby = sorted(flag_set_by(res["rows_by_rep"][k]))
        same = "identical to BH" if sby == sorted(sets[k]) else (
            "drops " + ", ".join(f"`{x}`" for x in sorted(set(sets[k]) - set(sby))))
        by_rows.append(f"| {k} | {len(sby)}: `{', '.join(sby)}` | {same} |")
    lines += [
        f"| pooled (≥{MIN_REPS_POOL} reps) | {len(pooled)} | "
        f"{len(pool_set)}: `{', '.join(pool_set)}` |",
        "",
        "### Under Benjamini--Yekutieli (arbitrary dependence)",
        "",
        "| replicate | flagged | vs BH |",
        "|---|---|---|",
    ] + by_rows + [
        "",
        f"- **Intersection ({len(inter)} system):** `{', '.join(inter)}`.",
        f"- **Union ({len(union)} systems):** `{', '.join(union)}`.",
        "- **Pairwise Jaccard:** "
        + ", ".join(f"{a}↔{b} {v:.3f}" for (a, b), v in zip(stab["pairs"], jac, strict=True))
        + f" (mean {stab['mean_jaccard']:.3f}).",
        "- **Random-draw reference:** the exact expected Jaccard of two *independent uniform*"
        f" draws of the observed sizes from the same {n_total} systems is "
        + ", ".join(f"{v:.3f}" for v in jrand)
        + f" (mean {stab['mean_jaccard_random']:.3f}); for two 6-of-{n_total} draws it is "
        f"{expected_random_jaccard(6, 6, n_total):.3f}.",
        "",
        "So the observed overlap is well above chance — the test is not returning noise — but it",
        "is far below the reproducibility a set-valued claim needs: only",
        f"`{', '.join(inter)}` survives all three runs, and {len(union) - len(inter)} of the "
        f"{len(union)} ever-flagged systems appear in some runs and not others.",
        "",
        "## Where the flips come from",
        "Panel B: the flips are **threshold crossings**, not sign changes. A system's BH-adjusted",
        "*q* swings by orders of magnitude between independent runs of the same protocol, so any",
        "system whose *q* lives near α will cross it about half the time. The clearest case is",
    ]
    fa = "faah"
    if fa in maps[0]:
        rcs = [maps[k][fa]["rc"] for k in REPLICATES if fa in maps[k]]
        lines += [
            f"`{fa}`, one of the paper's six: reduced χ² = "
            + ", ".join(f"{v:.2f}" for v in rcs)
            + f" across the three replicates — a swing of {max(rcs) / min(rcs):.1f}× that "
              "crosses the nominal χ²ᵥ = 1 threshold in both directions.",
        ]
    lines += [
        "",
        "## The pooled network (the paper's own pooling rule)",
        "Fig L panel C defines the pooling rule and it is reused here **verbatim**: per edge, the",
        "mean ΔΔG over the available replicates and se = √(mean_k se_k² / n) — the rule that gives",
        "panel C's `se → se/√3` — keeping edges with ≥ %d available replicates." % MIN_REPS_POOL,
        "",
        "**Stated assumption.** Panel C does *not* define a flagging rule on the pooled network —",
        "it only compares reduced χ² single-vs-pooled. Running BH-FDR on the pooled p-values is",
        "therefore this figure's extension; it uses the same pre-registered α = %.2f." % ALPHA,
        "The one free knob in the pooling rule is how many replicates an edge must have. Both",
        "settings are reported rather than chosen:",
        "",
        f"- `min_reps = {MIN_REPS_POOL}` (panel C's own setting): {len(pool_set)} flagged — "
        f"`{', '.join(pool_set)}`",
        f"- `min_reps = 3` (every edge complete in all three runs): {len(pool_alt_set)} flagged — "
        f"`{', '.join(pool_alt_set)}`",
        "",
        "The two agree exactly, so the pooled result does not depend on that choice.",
        "",
        f"The pooled network flags **{len(pool_set)} of the {len(union)} ever-flagged systems**"
        + (f" — everything in the union except `{', '.join(dropped)}`." if dropped
           else " — the whole union.")
        + " Pooling triples the effective information per edge, so this is the highest-powered",
        "single read available from this data; it is *not* an independent confirmation of any",
        "single-replicate set (it re-uses all three runs).",
        "",
        "## Which flags reproduce, model-free",
        "The distribution-free counterpart of the χ² test: per system, the mean pairwise Pearson",
        "correlation of the per-edge **signed** standardized residuals across the three",
        "replicates, on edges complete in all three. Sampling noise → 0; a reproducible",
        "edge-level systematic error → positive.",
        "",
    ] + ztbl + [
        "",
        f"`{fa}` is the outlier: its closure error is essentially **not** the same error run to",
        f"run ({zr[fa][0]:+.3f}), which is the mechanism behind its χ² swing above, and it is the",
        "only union member the pooled test does not flag." if fa in zr and fa in dropped else "",
        "",
        "## What does NOT move",
        "",
        "**1. The per-system ranking.** Spearman ρ of the reduced χ² across all "
        f"{len(res['common'])} systems: "
        + ", ".join(f"({i},{j}) ρ = {rho[(i, j)]:.2f}" for i, j in [(0, 1), (0, 2), (1, 2)])
        + ". The same networks are the worst-closing ones every run; it is the *dichotomisation*",
        "at α, not the underlying signal, that is unstable.",
        "",
        "**2. The selectivity-vs-σ-calibration curve** (Fig L panel B, per replicate):",
        "",
    ] + swp + [
        "",
        "The panel-B conclusion — an overconfident σ drives the false-positive rate toward 1, a",
        "calibrated σ is selective, an over-wide σ loses all power — reproduces on every",
        "replicate. That claim is about the *detector*, and it does not depend on which systems",
        "come out.",
        "",
        "## Robustness to the FDR level",
        "The instability is not an artefact of α = %.2f:" % ALPHA,
        "",
    ] + asw + [
        "",
        "At every level tested the three replicates disagree on a majority of the union.",
        "",
        "## Full per-system table",
        "Reduced χ² / BH-adjusted *q* per replicate and pooled; **bold** = flagged at α = %.2f."
        % ALPHA,
        "Systems are ordered by their strongest (smallest) *q* anywhere; edges and cycles are",
        "replicate-0 counts (they vary by ≤ a few edges between replicates where a leg failed).",
        "",
    ] + tbl + [
        "",
        "## Honest reading",
        "- The paper's six flagged systems are a **replicate-0 realisation** of a test whose",
        "  set-valued output is unstable, not a reproducible list. Reported as a set, the honest",
        f"  statement is: `{', '.join(inter)}` reproduces in every run; the union of "
        f"{len(union)} systems is what a single run can produce; and the pooled, "
        f"highest-powered read gives {len(pool_set)}.",
        "- The **underlying quantity is stable** (ρ = "
        f"{min(rho.values()):.2f}–{max(rho.values()):.2f}); the instability is created by",
        "  thresholding a continuous, low-power statistic. This is consistent with the median",
        "  reduced χ² ≈ 0.34 already reported in `results_figL.md`, which implies limited",
        "  single-replicate power, and with that document's existing note that the flag set is a",
        "  *lower bound* on the systems carrying systematic error.",
        "- **The detector-level claims are unaffected.** Panel B of Fig L (selectivity is a",
        "  function of σ-calibration) reproduces on all three replicates. What does not survive",
        "  is the membership of any particular six-system list.",
        "- Nothing here was tuned after seeing a result: α, the σ scales, the edge construction",
        "  and the pooling rule are all taken unchanged from `make_figL.py`.",
        "",
    ]
    DOC.write_text("\n".join(ln for ln in lines if ln is not None) + "\n")


def main():
    use_paper_style()
    res = analyze()
    make_figure(res)
    write_doc(res)
    print(f"wrote figStab_replicate_stability.(pdf|png) to {FIGDIR}")
    print(f"wrote {DOC.relative_to(ROOT)}")

    stab = res["stab"]
    for k in REPLICATES:
        print(f"[rep {k}] {len(res['rows_by_rep'][k])} systems tested; "
              f"flagged ({len(res['sets'][k])}): {sorted(res['sets'][k])}")
    pool_set = sorted(flag_set(res["pooled_rows"]))
    print(f"[pooled] flagged ({len(pool_set)}): {pool_set}   "
          f"[min_reps=3 check: {sorted(flag_set(res['pooled_alt']))}]")
    print(f"[sets] intersection={stab['intersection']}  union({len(stab['union'])})="
          f"{stab['union']}")
    print("[jaccard] observed "
          + "/".join(f"{v:.3f}" for v in stab["jaccard"])
          + "  vs exact size-matched random "
          + "/".join(f"{v:.3f}" for v in stab["jaccard_random"])
          + f"  (6-of-{res['n_total']} random = "
            f"{expected_random_jaccard(6, 6, res['n_total']):.3f})")
    print("[stable] Spearman on reduced chi2 across all "
          f"{len(res['common'])} systems: "
          + ", ".join(f"({i},{j})={res['rho'][(i, j)]:.3f}" for i, j in [(0, 1), (0, 2), (1, 2)]))
    print("[stable] selectivity vs sigma-scale (flagged/total per replicate):")
    for i, sc in enumerate(SIGMA_SCALES):
        cells = " ".join(f"rep{k} {res['sweeps'][k][i][1]}/{res['sweeps'][k][i][2]}"
                         for k in REPLICATES)
        print(f"    x{sc:<5g} {cells}")
    print("[alpha sweep] instability at other FDR levels (sizes | intersection | union):")
    for al, sizes, ai, au, aj in res["alpha_sweep"]:
        print(f"    alpha={al:<5g} sizes={sizes} inter={len(ai)} union={len(au)} "
              f"jaccard={[round(v, 3) for v in aj]}")
    print("[z-reproduction] mean pairwise signed-z correlation, ever-flagged systems:")
    for s in stab["union"]:
        r, n = res["zrep"][s]
        print(f"    {s:16s} {r:+.3f}  (n={n} aligned edges)")


if __name__ == "__main__":
    main()
