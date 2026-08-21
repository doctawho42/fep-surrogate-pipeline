"""Fig Cut -- P1 fixed-cutoff head-to-head + P2 chi^2 reconciliation (peer-review items).

P1: three detectors on the SAME systems (calibrated per-edge null / fixed 1.0 kcal/mol hysteresis
cutoff / one pooled se), scored by Mann-Whitney AUC against the cross-replicate reproducibility
anchor. PRE-REGISTERED: WIN iff the calibrated null beats BOTH foils and the paired bootstrap CI
excludes zero; otherwise TIE.

P2: is the gap between the closure-implied conservatism (median reduced chi^2 0.34 -> ~1.71x in se)
and the replicate-validated 1.25-1.41x explained by correlated residuals among edges sharing a
ligand endpoint? Compare the empirical pair correlation against the correlation the residual-maker
induces under the null, then recompute the implied factor against the effective dof.
PRE-REGISTERED: CONFIRMED iff the shared-node EXCESS correlation is positive with a CI excluding
zero AND the corrected factor closes at least half the 1.71 -> 1.41 gap (i.e. <= 1.56).
MANDATORY either way: re-run BH-FDR under the effective dof and report whether the flag set moved.

Run:  PYTHONPATH=src python figs/make_figCut.py   (or `make figCut`).  Deterministic.
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

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
from bar.detectors import (  # noqa: E402
    HYSTERESIS_CUTOFF,
    anchor_score,
    flag_calibrated,
    flag_fixed_cutoff,
    flag_fixed_se,
    paired_auc_bootstrap,
)
from bar.qc import benjamini_hochberg, chi2_sf, gls_network  # noqa: E402
from bar.residcorr import (  # noqa: E402
    effective_dof,
    empirical_pair_correlation,
    null_pair_correlation,
    pair_masks,
    residual_maker,
)
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from paperstyle import (  # noqa: E402
    ALT, INK, MUTED, OURS, REF, THIRD, check_min_type, figsize, finish, panel,
    reference_line, use_paper_style,
)

DATA = ROOT / "data" / "openfe_replicates" / "combined_pymbar4_edge_data.csv"
FIGDIR = pathlib.Path(__file__).resolve().parent
N_BOOT = 2000
FACTOR_CLOSURE = 1.71     # pre-registered: sqrt(1/0.34), the closure-implied se factor
FACTOR_REPLICATE = 1.41   # pre-registered: the replicate-validated se factor (raw)
FACTOR_HALFWAY = 1.56     # pre-registered CONFIRMED threshold (closes >= half the gap)
# Semantic colours (paperstyle): OURS = the calibrated per-edge null and its AUC bar; the two
# rival flagging rules are named methods that are neither ours nor an overconfident stand-in, so
# they take ALT and THIRD; REF grey is for the reference lines.
# Panel B's cloud is NOT a detector output. It is every system's reduced chi^2 under two dof
# conventions, and the panel singles out no subset of them, so the whole cloud is MUTED -- the
# same grey that carries "all systems" in Fig L(C), Fig Lev(A) and Fig Stab(C), where blue means
# "the systems the detector flagged". Drawing all 46 in OURS blue here invited exactly the
# reading -- these are the flags -- that this panel is not making.
C_OURS, C_CUTOFF, C_FIXEDSE = OURS, ALT, THIRD
C_CLOUD = MUTED


def _f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return math.nan


def edge_val(r, k):
    cD = _f(r[f"complex_repeat_{k}_DG (kcal/mol)"])
    cd = _f(r[f"complex_repeat_{k}_dDG (kcal/mol)"])
    sD = _f(r[f"solvent_repeat_{k}_DG (kcal/mol)"])
    sd = _f(r[f"solvent_repeat_{k}_dDG (kcal/mol)"])
    if any(math.isnan(v) for v in (cD, cd, sD, sd)):
        return None
    return cD - sD, math.sqrt(cd * cd + sd * sd)


def load_systems():
    by = defaultdict(list)
    for r in csv.DictReader(open(DATA)):
        by[r["system name"]].append(r)
    return by


def build(by):
    """Per system: replicate-0 edges (for the detectors) and the 3 aligned replicate edge lists
    (for the anchor). Only rows complete in ALL three replicates enter the aligned set."""
    rep0, aligned = {}, {}
    for name, rows in sorted(by.items()):
        e0 = [(r["ligand_A"], r["ligand_B"], *edge_val(r, 0)) for r in rows if edge_val(r, 0)]
        if len(e0) < 3 or gls_network(e0).dof < 1:
            continue
        rep0[name] = e0
        keep = [r for r in rows if all(edge_val(r, k) for k in (0, 1, 2))]
        if len(keep) >= 3:
            reps = [[(r["ligand_A"], r["ligand_B"], *edge_val(r, k)) for r in keep]
                    for k in (0, 1, 2)]
            if all(gls_network(e).dof >= 1 for e in reps):
                aligned[name] = reps
    return rep0, aligned


def run_p1(rep0, aligned):
    fa, fb, fc = flag_calibrated(rep0), flag_fixed_cutoff(rep0), flag_fixed_se(rep0)
    names = [n for n in sorted(rep0) if n in fa and n in fb and n in fc and n in aligned]
    anchor = np.array([anchor_score(aligned[n]) for n in names])
    A = np.array([fa[n] for n in names])
    B = np.array([fb[n] for n in names])
    C = np.array([fc[n] for n in names])
    out = paired_auc_bootstrap(A, B, C, anchor, n_boot=N_BOOT, seed=0)
    out.update(names=names, flags_a=A, flags_b=B, flags_c=C, anchor=anchor)
    return out


def _exclusion_reason(shared, disjoint) -> str:
    """Which half of the P2 eligibility guard (``not shared.any() or not disjoint.any()``)
    actually fired for an excluded system -- reflects today's data (all excluded systems are
    all-shared-node triangles) and stays correct if the other branch ever triggers."""
    if not shared.any() and not disjoint.any():
        return "no shared or disjoint edge pairs"
    if not shared.any():
        return "no shared-node edge pairs"
    return "no disjoint edge pairs"


def run_p2(aligned):
    """Pooled excess correlation + per-system effective dof, with the flag-set stability check."""
    emp_s, emp_d, null_s, null_d, per_sys = [], [], [], [], []
    excluded = []
    for name, reps in sorted(aligned.items()):
        edges = reps[0]
        M = residual_maker(edges)
        shared, disjoint = pair_masks(edges)
        if not shared.any() or not disjoint.any():
            excluded.append((name, _exclusion_reason(shared, disjoint)))
            continue
        z = np.vstack([gls_network(e).z for e in reps])
        es, ed = empirical_pair_correlation(z, shared), empirical_pair_correlation(z, disjoint)
        ns, nd = null_pair_correlation(M, shared), null_pair_correlation(M, disjoint)
        emp_s.append(es)
        emp_d.append(ed)
        null_s.append(ns)
        null_d.append(nd)
        per_sys.append((name, edges, M, shared, disjoint, es - ns, ed - nd))
    exc_s = np.array(emp_s) - np.array(null_s)
    exc_d = np.array(emp_d) - np.array(null_d)
    rng = np.random.default_rng(0)
    boot = np.array([np.mean(exc_s[rng.integers(0, exc_s.size, exc_s.size)])
                      for _ in range(N_BOOT)])
    lo, hi = (float(v) for v in np.percentile(boot, [2.5, 97.5]))

    rc_nom, rc_eff, names2, p_nom, p_eff = [], [], [], [], []
    rs, rd = float(np.mean(exc_s)), float(np.mean(exc_d))
    for name, edges, M, shared, disjoint, _es, _ed in per_sys:
        fit = gls_network(edges)
        de = effective_dof(M, shared, disjoint, rs, rd)
        if de <= 0:
            continue
        names2.append(name)
        rc_nom.append(fit.reduced_chi2)
        rc_eff.append(fit.chi2 / de)
        p_nom.append(chi2_sf(fit.chi2, fit.dof))
        p_eff.append(chi2_sf(fit.chi2, de))
    med_nom, med_eff = float(np.median(rc_nom)), float(np.median(rc_eff))
    fac_nom, fac_eff = math.sqrt(1.0 / med_nom), math.sqrt(1.0 / med_eff)
    flag_nom = {n for n, f in zip(names2, benjamini_hochberg(p_nom), strict=True) if f}
    flag_eff = {n for n, f in zip(names2, benjamini_hochberg(p_eff), strict=True) if f}
    confirmed = (rs > 0 and lo > 0 and fac_eff <= FACTOR_HALFWAY)
    return {"exc_shared": rs, "exc_disjoint": rd, "ci_lo": lo, "ci_hi": hi,
            "median_rc_nominal": med_nom, "median_rc_eff": med_eff,
            "factor_nominal": fac_nom, "factor_eff": fac_eff,
            "flag_nominal": sorted(flag_nom), "flag_eff": sorted(flag_eff),
            "flag_set_changed": flag_nom != flag_eff,
            "verdict": "CONFIRMED" if confirmed else "NOT-CONFIRMED",
            "rc_nom": rc_nom, "rc_eff": rc_eff,
            "n_systems": len(names2), "excluded": excluded}


def _signed(x, nd=3):
    """Format a signed number with the real Unicode minus sign, for figure text."""
    return f"{x:+.{nd}f}".replace("-", "\u2212")


def main():
    use_paper_style()
    rep0, aligned = build(load_systems())
    p1, p2 = run_p1(rep0, aligned), run_p2(aligned)

    fig, (axA, axB) = plt.subplots(1, 2, figsize=figsize(2, 3.15))
    labels = ["calibrated\n(ours)", f"fixed cutoff\n({HYSTERESIS_CUTOFF} kcal/mol)",
              "fixed se\n(pooled)"]
    vals = [p1["auc_a"], p1["auc_b"], p1["auc_c"]]
    axA.bar(range(3), vals, width=0.66, color=[C_OURS, C_CUTOFF, C_FIXEDSE], zorder=2)
    reference_line(axA, "hline", 0.5)
    # the chance line crosses all three bars, so its name goes in the clear strip to their
    # right rather than on top of one of them.
    axA.text(2.86, 0.512, "chance", ha="right", va="bottom", fontsize=7.5, color=REF)
    for i, v in enumerate(vals):
        axA.text(i, v + 0.015, f"{v:.2f}", ha="center", va="bottom", fontsize=8, color=INK)
    axA.set_xticks(range(3))
    axA.set_xticklabels(labels, fontsize=8)
    axA.set_xlim(-0.67, 2.9)
    axA.set_ylabel("AUC vs cross-replicate reproducibility")
    axA.set_ylim(0, 1.02)
    axA.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    axA.spines["left"].set_bounds(0, 1.0)
    axA.tick_params(axis="x", length=0)
    # The heading zone is a letter and one short line, like every other panel in this article.
    # The caveat that used to sit here as a third, centred, hard-wrapped block -- that the scored
    # label is standardized by the same reported se the calibrated rule weights by, so no
    # ordering among these bars is a discrimination gain -- is a sentence, not a heading; it is
    # carried by the caption of Figure~\ref{fig:Cut} in docs/paper_si.tex.
    panel(axA, "A", f"head-to-head detection AUC: {p1['verdict'].lower()}",
          subtitle=(f"\u0394 = {_signed(p1['diff'])} "
                    f"[{_signed(p1['ci_lo'])}, {_signed(p1['ci_hi'])}]"))

    lim = [min(p2["rc_nom"] + p2["rc_eff"]) * 0.8, max(p2["rc_nom"] + p2["rc_eff"]) * 1.2]
    # two reference lines: same REF grey, told apart by dash pattern.
    axB.plot(lim, lim, color=REF, ls=(0, (1.0, 2.0)), lw=1.0, zorder=1)
    axB.axhline(1.0, color=REF, ls=(0, (4.5, 2.0)), lw=1.0, zorder=1)
    axB.scatter(p2["rc_nom"], p2["rc_eff"], s=16, c=C_CLOUD, lw=0, alpha=0.8, zorder=3)
    axB.set_xscale("log")
    axB.set_yscale("log")
    axB.set_xlim(*lim)
    axB.set_ylim(*lim)
    axB.set_box_aspect(1.0)
    axB.text(lim[1] * 0.92, 1.0, r"$\chi^2_\nu=1$", ha="right", va="bottom", fontsize=7.5,
             color=REF)
    axB.set_xlabel(r"reduced $\chi^2$ (nominal dof)")
    axB.set_ylabel(r"reduced $\chi^2$ (effective dof)")
    panel(axB, "B", "shared-node residual correlation: "
                    f"{'supported' if p2['verdict'] == 'CONFIRMED' else 'not supported'}",
          subtitle=(f"se factor {p2['factor_nominal']:.2f}\u00d7 \u2192 "
                    f"{p2['factor_eff']:.2f}\u00d7, {len(p2['rc_nom'])} systems"))
    # the points lie on the identity, so the only clear ground in this panel is the wedge
    # above it: the pre-registration note is set there, opaque, and wrapped narrow enough
    # that its lower-right corner stays clear of the diagonal.
    axB.text(0.02, 0.985,
             f"pre-registered success needed the\ncorrected factor at or below "
             f"{FACTOR_HALFWAY}\u00d7, half-way\nfrom the closure-implied {FACTOR_CLOSURE}\u00d7 to\n"
             f"the replicate-validated {FACTOR_REPLICATE}\u00d7",
             transform=axB.transAxes, ha="left", va="top", fontsize=7.5, color=REF,
             linespacing=1.3)

    offenders = check_min_type(fig)
    assert not offenders, f"type below the house floor: {offenders}"
    finish(fig, "figCut_cutoff_benchmark")

    print(f"\n[P1] n_systems={len(p1['names'])}  "
          f"flagged: calibrated={int(p1['flags_a'].sum())}, "
          f"fixed-cutoff={int(p1['flags_b'].sum())}, fixed-se={int(p1['flags_c'].sum())}")
    print(f"[P1] AUC calibrated={p1['auc_a']:.3f}  fixed-cutoff={p1['auc_b']:.3f}  "
          f"fixed-se={p1['auc_c']:.3f}")
    print(f"[P1] paired diff A-max(B,C) = {p1['diff']:+.3f} "
          f"[{p1['ci_lo']:+.3f}, {p1['ci_hi']:+.3f}]")
    print(f"P1 VERDICT: {p1['verdict']}")
    disc = [(n, bool(a), bool(b), bool(c)) for n, a, b, c
            in zip(p1["names"], p1["flags_a"], p1["flags_b"], p1["flags_c"], strict=True)
            if not (a == b == c)]
    print(f"[P1] discordant systems ({len(disc)}): "
          + ", ".join(f"{n}(A={a},B={b},C={c})" for n, a, b, c in disc))

    exc_reasons = {r for _, r in p2["excluded"]}
    if len(exc_reasons) <= 1:
        reason = exc_reasons.pop() if exc_reasons else "n/a"
        print(f"\n[P2] n_systems={p2['n_systems']}  "
              f"excluded ({len(p2['excluded'])}, {reason}): "
              + ", ".join(n for n, _ in p2["excluded"]))
    else:
        print(f"\n[P2] n_systems={p2['n_systems']}  excluded ({len(p2['excluded'])}): "
              + ", ".join(f"{n} ({r})" for n, r in p2["excluded"]))
    print(f"[P2] excess corr shared-node={p2['exc_shared']:+.4f} "
          f"[{p2['ci_lo']:+.4f}, {p2['ci_hi']:+.4f}]  disjoint={p2['exc_disjoint']:+.4f}")
    print(f"[P2] median reduced chi2 {p2['median_rc_nominal']:.3f} (nominal) -> "
          f"{p2['median_rc_eff']:.3f} (effective dof)")
    print(f"[P2] implied se factor {p2['factor_nominal']:.2f}x -> {p2['factor_eff']:.2f}x "
          f"(replicate-validated {FACTOR_REPLICATE}x; CONFIRMED needs <= {FACTOR_HALFWAY})")
    print(f"[P2] BH-FDR flag set changed: {p2['flag_set_changed']}")
    print(f"[P2]   nominal : {p2['flag_nominal']}")
    print(f"[P2]   effective: {p2['flag_eff']}")
    print(f"P2 VERDICT: {p2['verdict']}")


if __name__ == "__main__":
    main()
