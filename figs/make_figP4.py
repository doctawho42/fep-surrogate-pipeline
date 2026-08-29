"""Fig P4 -- the QC calibration sweep under the REAL heterogeneous learned-sigma profile.

Fig L panel B shrank every edge's se by a uniform x0.15 to represent an overconfident learned
sigma, and reported that ~88% of systems get flagged. Referees objected that a uniform shrink is
near-mechanically forced. This script pushes the learned head's MEASURED, overlap-dependent profile
through the identical GLS + BH-FDR test, per edge, plus a shuffled control.

Arms: calibrated (x1); real heterogeneous profile (PRIMARY); shuffled profile (control, same ratio
multiset, association with overlap destroyed); uniform x0.15 (retained for continuity, but a STRESS
TEST, not the representative learned-head comparator).

PRE-REGISTERED: DEGRADES iff the real-profile arm flags at least 2x the calibrated arm's flagged
count; otherwise COMPARABLE, which means the "an overconfident sigma destroys the QC" claim is not
supported by a realistic profile and the manuscript must soften it.

Overlap caveat: Fig A's normalized BAR overlap and OpenFE's pymbar smallest_overlap are different
measurements, so the profile is transferred by RANK, not raw value (see bar/sigma_profile.py).

Run:  PYTHONPATH=src python figs/make_figP4.py   (or `make figP4`).  Deterministic.
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
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from bar.qc import benjamini_hochberg, chi2_sf, gls_network  # noqa: E402
from bar.sigma_profile import PROFILE_POINTS, rank_transfer, shuffled  # noqa: E402
from paperstyle import (  # noqa: E402
    FOIL, INK, MUTED, OURS, figsize, finish, panel, tint, use_paper_style,
)

DATA = ROOT / "data" / "openfe_replicates" / "combined_pymbar4_edge_data.csv"
FIGDIR = pathlib.Path(__file__).resolve().parent
SHUFFLE_SEED = 20260808          # frozen
UNIFORM_STRESS = 0.15            # the retained Fig L stress-test factor
DEGRADE_FACTOR = 2.0             # frozen pre-registered verdict threshold
# Semantic colours (paperstyle): OURS = the calibrated sandwich; FOIL = the real learned-variance
# head's profile and the uniform stand-in, which are two members of ONE family and so take one hue
# at two tints; MUTED = the shuffled control, which is a null.
ARM_COLOURS = (OURS, FOIL, MUTED, tint(FOIL, 0.55))
# Axis labels for the four arms. Kept separate from the printed arm names above so that the
# console output this script's results record quotes stays byte-identical.
ARM_TICKS = (
    "calibrated\nsandwich ×1",
    "real learned\nprofile",
    "shuffled\ncontrol",
    "uniform ×0.15\nstand-in",
)


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
    """The edge's overlap = the MIN of its two legs (the worse leg limits the edge). Frozen."""
    c = _f(r[f"complex_repeat_{k}_smallest_overlap"])
    s = _f(r[f"solvent_repeat_{k}_smallest_overlap"])
    if math.isnan(c) or math.isnan(s):
        return math.nan
    return min(c, s)


def load():
    """Per system: the replicate-0 edges plus each edge's overlap, keeping only systems with a cycle."""
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


def flagged(systems, scale_fn):
    """Flagged count under a per-system per-edge se scaling. `scale_fn(name, overlaps) -> ratios`."""
    names, ps = [], []
    for name, (edges, ov) in sorted(systems.items()):
        rat = scale_fn(name, ov)
        scaled = [(a, b, y, se * float(f)) for (a, b, y, se), f in zip(edges, rat, strict=True)]
        fit = gls_network(scaled)
        if fit.dof < 1:
            continue
        names.append(name)
        ps.append(chi2_sf(fit.chi2, fit.dof))
    fl = benjamini_hochberg(ps)
    return int(fl.sum()), len(names), {n for n, f in zip(names, fl, strict=True) if f}


def main():
    use_paper_style()
    systems = load()
    # Global percentile ranking (frozen): pool every edge, transfer, then split back per system.
    order = sorted(systems)
    all_ov = np.concatenate([systems[n][1] for n in order])
    all_rat = rank_transfer(all_ov)
    all_shuf = shuffled(all_rat, seed=SHUFFLE_SEED)
    idx, per_sys_rat, per_sys_shuf = 0, {}, {}
    for n in order:
        k = systems[n][1].size
        per_sys_rat[n] = all_rat[idx:idx + k]
        per_sys_shuf[n] = all_shuf[idx:idx + k]
        idx += k

    arms = [
        ("calibrated\nsandwich (x1)", lambda n, ov: np.ones(ov.size)),
        ("real learned\nprofile (primary)", lambda n, ov: per_sys_rat[n]),
        ("shuffled profile\n(control)", lambda n, ov: per_sys_shuf[n]),
        (f"uniform x{UNIFORM_STRESS}\n(stress test)", lambda n, ov: np.full(ov.size, UNIFORM_STRESS)),
    ]
    res = [(label, *flagged(systems, fn)) for label, fn in arms]
    cal_n = res[0][1]
    real_n = res[1][1]
    verdict = "DEGRADES" if real_n >= DEGRADE_FACTOR * cal_n else "COMPARABLE"

    fig, (ax, ax2) = plt.subplots(
        1, 2, figsize=figsize(2, 3.0), gridspec_kw={"width_ratios": [1.05, 1.0]},
    )
    frac = [100.0 * r[1] / r[2] for r in res]
    ax.bar(range(len(res)), frac, width=0.66, color=list(ARM_COLOURS), zorder=2)
    for i, r in enumerate(res):
        ax.text(i, frac[i] + 2.5, f"{r[1]}/{r[2]}", ha="center", va="bottom", fontsize=7.5,
                color=INK)
    ax.set_xticks(range(len(res)))
    ax.set_xticklabels(list(ARM_TICKS), fontsize=7.5)
    # bars are 0.66 wide on unit centres; give the y-spine the same 0.34 gap on each side.
    ax.set_xlim(-0.67, len(res) - 0.33)
    ax.tick_params(axis="x", length=0)
    ax.set_ylabel("% of systems flagged")
    ax.set_ylim(0, 100)
    ax.set_yticks([0, 25, 50, 75, 100])
    panel(ax, "A", r"flag rate vs per-edge $\sigma$")

    # Second panel: PROFILE_POINTS has only 4 knots, but np.interp is piecewise-linear between
    # them and each of the 1143 edges sits at its own distinct percentile, so the 1143 transferred
    # ratios are themselves (near-)continuous, not four repeated bar-chart-style levels -- this
    # histogram makes the "heterogeneous assignment, narrow deep-shrink band" argument visible
    # instead of prose-only, and marks Fig L panel B's uniform x0.15 stand-in for comparison.
    # The profile and the stand-in are one family: one hue, the histogram at a lighter tint so
    # the stand-in's marker line reads on top of it without any transparency.
    ax2.hist(all_rat, bins=40, color=tint(FOIL, 0.45), edgecolor="white", linewidth=0.5,
             zorder=2)
    ax2.axvline(UNIFORM_STRESS, color=FOIL, linestyle="--", linewidth=1.3, zorder=3)
    ax2.set_xlabel(r"transferred se ratio (real learned-$\sigma$ profile)")
    ax2.set_ylabel("edge count")
    ax2.set_ylim(0, 80)
    # the marker line is named where it stands, in the empty upper right, rather than in a
    # legend whose handle would sit a few points from the line itself and read as part of it.
    ax2.text(UNIFORM_STRESS * 1.02, 77.0, f"uniform ×{UNIFORM_STRESS}\nstand-in", ha="left",
             va="top", fontsize=7.5, color=FOIL, zorder=4)
    panel(ax2, "B", f"{all_rat.size} per-edge transferred ratios")

    finish(fig, "figP4_sigma_profile")

    print(f"\n[P4] profile points (overlap -> ratio): {PROFILE_POINTS}")
    print(f"[P4] edges={all_ov.size}  systems={res[0][2]}  "
          f"ratio range {all_rat.min():.3f}-{all_rat.max():.3f}")
    for label, n, tot, _ in res:
        print(f"[P4] {label.replace(chr(10), ' '):38s} flagged {n:>2}/{tot} ({100*n/tot:.0f}%)")
    print(f"[P4] real-vs-calibrated ratio: {real_n / cal_n:.2f}x "
          f"(DEGRADES needs >= {DEGRADE_FACTOR:g}x)")
    print(f"P4 VERDICT: {verdict}")
    only_real = res[1][3] - res[0][3]
    only_cal = res[0][3] - res[1][3]
    print(f"[P4] flagged by real profile but not calibrated: {sorted(only_real)}")
    print(f"[P4] flagged by calibrated but not real profile: {sorted(only_cal)}")

    # Every arm's flagged set, and the pairwise discordance between the three non-calibrated
    # arms (the SI/results doc cite these numbers; they must come from this script, not a
    # separate ad hoc check).
    for label, _n, _tot, flagset in res:
        print(f"[P4] {label.replace(chr(10), ' '):38s} flagged set: {sorted(flagset)}")
    real_s, shuf_s, unif_s = res[1][3], res[2][3], res[3][3]
    for a_name, a_s, b_name, b_s in (
        ("real", real_s, "shuffled", shuf_s),
        ("real", real_s, "uniform", unif_s),
        ("shuffled", shuf_s, "uniform", unif_s),
    ):
        print(f"[P4] {a_name} - {b_name} = {sorted(a_s - b_s)}   "
              f"{b_name} - {a_name} = {sorted(b_s - a_s)}")
    print(f"[P4] distinct transferred ratios (real arm): {len(set(all_rat.tolist()))}/{all_rat.size}")


if __name__ == "__main__":
    main()
