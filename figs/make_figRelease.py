"""Fig Release: how far between-preparation variance exceeds sampling variance.

Three passages of the article ask the same question and answer it the same way -- that public
data cannot supply the number. They are the disagreement with Wade et al. over whether the
reported bar is conservative, the observation that nothing flags once the bars are widened by
half, and the caveat that the three replicates share one starting structure, so what they
establish is reproducibility rather than physical origin. All three want one quantity: the size
of the variance between independent PREPARATIONS of the same target, relative to the sampling
variance the reported bar carries.

The benchmark supplies it in exactly one place. A system is named by target, and eight system
names union edges from more than one separately released data set. Where two releases share
ligands, the union carries cycles that exist in neither release alone, and those cycles close
only if the two preparations agree. Splitting the closure statistic into its within-release and
across-release blocks therefore measures the excess directly.

This script runs that split over EVERY multi-release system and all THREE replicates, and
reports how many systems actually carry the measurement. ``release_decomposition`` in
``figs/make_figHodge.py`` performs the split on replicate 0; this script imports it, reuses its
per-edge reader, and cross-checks its own replicate-0 output against it edge for edge.

Nothing here is tuned: the split is fixed by the incidence structure, there is no threshold to
choose, and the systems that carry no across-release cycle are reported as carrying none.

Run: PYTHONPATH=src python figs/make_figRelease.py   (or `make figRelease`)
"""
from __future__ import annotations

import csv
import math
import pathlib
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "figs"))
from paperstyle import (  # noqa: E402
    MUTED,
    OURS,
    REF,
    check_min_type,
    figsize,
    finish,
    legend,
    panel,
    tint,
    use_paper_style,
)

import make_figHodge as hodge  # noqa: E402  (the split this script generalises)
from bar.qc import chi2_sf, gls_network  # noqa: E402

DOC = ROOT / "docs" / "results_figRelease.md"
COMBINED = ROOT / "data" / "openfe_replicates" / "combined_pymbar4_edge_data.csv"
REPLICATES = (0, 1, 2)


# ------------------------------------------------------------------------------------------
# the split, generalised over replicates
# ------------------------------------------------------------------------------------------
def _rows_by_system() -> dict[str, list[dict]]:
    by: dict[str, list[dict]] = {}
    for row in csv.DictReader(COMBINED.open()):
        by.setdefault(row["system name"], []).append(row)
    return by


def _split(per_release: dict[str, list]) -> dict | None:
    """Within-release / across-release decomposition of one system's closure statistic.

    Identical in form to ``make_figHodge.release_decomposition``; the only change is that the
    edges arrive already grouped, so any replicate can be passed in. The releases' cycle spaces
    have disjoint edge supports, so they sit orthogonally inside the union's cycle space and
    both the chi-squared and the degrees of freedom add exactly.
    """
    edges = [e for block in per_release.values() for e in block]
    if len(per_release) < 2 or len(edges) < 3:
        return None
    union = gls_network(edges)
    if union.dof < 1:
        return None
    within_chi2 = within_dof = 0.0
    for block in per_release.values():
        if len(block) < 2:
            continue
        fit = gls_network(block)
        if fit.dof < 1:
            continue
        within_chi2 += fit.chi2
        within_dof += fit.dof
    across_dof = union.dof - int(within_dof)
    across_chi2 = union.chi2 - within_chi2
    return {"releases": len(per_release), "E": len(edges),
            "union_chi2": union.chi2, "union_dof": union.dof,
            "within_chi2": within_chi2, "within_dof": int(within_dof),
            "within_reduced": within_chi2 / within_dof if within_dof > 0 else math.nan,
            "across_chi2": across_chi2, "across_dof": across_dof,
            "across_reduced": across_chi2 / across_dof if across_dof > 0 else math.nan,
            "across_p": chi2_sf(across_chi2, across_dof) if across_dof > 0 else math.nan}


def _grouped(rows: list[dict], replicate: int) -> dict[str, list]:
    per_release: dict[str, list] = {}
    for r in rows:
        val = hodge._edge(r, replicate)
        if not val:
            continue
        per_release.setdefault(r["dataset_name"], []).append((r["ligand_A"], r["ligand_B"], *val))
    return per_release


def _shared(per_release: dict[str, list]) -> dict:
    """Ligands each pair of releases has in common, and the total appearing in more than one."""
    ligs = {name: {e[0] for e in block} | {e[1] for e in block}
            for name, block in per_release.items()}
    names = sorted(ligs)
    pairs = [(a, b, len(ligs[a] & ligs[b])) for i, a in enumerate(names) for b in names[i + 1:]]
    seen: dict[object, int] = {}
    for s in ligs.values():
        for lig in s:
            seen[lig] = seen.get(lig, 0) + 1
    return {"pairs": pairs, "max_pair": max((n for _, _, n in pairs), default=0),
            "n_multi": sum(1 for v in seen.values() if v > 1),
            "sizes": {k: len(v) for k, v in ligs.items()}}


def _loo_across(per_release: dict[str, list]) -> list[tuple[float, int]]:
    """Across-release reduced chi-squared under single-edge deletion, over every edge.

    Deleting an edge can remove an across-release cycle as well, so the degrees of freedom are
    carried alongside the value rather than assumed fixed.
    """
    out = []
    keys = list(per_release)
    for name in keys:
        for k in range(len(per_release[name])):
            trial = {n: (list(b) if n != name else [e for i, e in enumerate(b) if i != k])
                     for n, b in per_release.items()}
            trial = {n: b for n, b in trial.items() if b}
            res = _split(trial)
            if res and res["across_dof"] > 0:
                out.append((res["across_reduced"], res["across_dof"]))
    return out


def _inflate(per_release: dict[str, list], sigma: float) -> dict[str, list]:
    """The same edges with an extra independent variance ``sigma**2`` on every reported bar."""
    return {name: [(a, b, ddg, math.sqrt(se * se + sigma * sigma)) for a, b, ddg, se in block]
            for name, block in per_release.items()}


def between_sigma(per_release: dict[str, list], hi: float = 8.0) -> float:
    """Per-edge between-preparation term, in kcal/mol, that would close the across-release block.

    Solves ``across-release reduced chi-squared == 1`` for an extra independent variance added to
    every reported bar. Unlike the multiplicative widening factor this is not closed-form -- the
    added term changes the GLS weights and so the fit -- so it is bisected. Returns ``0.0`` when
    the block is already at or below one, which is the honest answer there: these data ask for no
    extra term at all.
    """
    def reduced(sigma: float) -> float:
        res = _split(_inflate(per_release, sigma))
        return res["across_reduced"] if res and res["across_dof"] > 0 else math.nan

    if not (reduced(0.0) > 1.0):
        return 0.0
    lo = 0.0
    if reduced(hi) > 1.0:
        return math.nan
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if reduced(mid) > 1.0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def null_calibration(per_release: dict[str, list], n_draw: int = 2000,
                     seed: int = 0) -> dict:
    """Does the across-release block have the null it is being read against?

    Synthetic data on the system's real topology and real reported bars: one node potential per
    ligand, SHARED by both releases, plus independent Gaussian noise of the reported size. Under
    that null the two preparations agree exactly and the across-release block should be a
    chi-squared on its stated degrees of freedom: mean reduced chi-squared 1, and 5 per cent of
    draws below p = 0.05. A departure here would mean the split's bookkeeping, not the data, is
    producing the excess. Nothing about this check is tuned; it is run once at a fixed seed.
    """
    rng = np.random.default_rng(seed)
    nodes = sorted({e[0] for b in per_release.values() for e in b}
                   | {e[1] for b in per_release.values() for e in b})
    index = {n: i for i, n in enumerate(nodes)}
    reduced, hits = [], 0
    for _ in range(n_draw):
        phi = rng.normal(size=len(nodes))
        draw = {name: [(a, b, phi[index[b]] - phi[index[a]] + rng.normal(0.0, se), se)
                       for a, b, _, se in block] for name, block in per_release.items()}
        res = _split(draw)
        if res is None or res["across_dof"] < 1:
            continue
        reduced.append(res["across_reduced"])
        hits += res["across_p"] < 0.05
    return {"n": len(reduced), "mean_reduced": float(np.mean(reduced)),
            "rate_p05": hits / len(reduced) if reduced else math.nan}


def measure() -> list[dict]:
    """The split for every multi-release system on every replicate."""
    by = _rows_by_system()
    out = []
    for name, rows in sorted(by.items()):
        if len({r["dataset_name"] for r in rows}) < 2:
            continue
        for replicate in REPLICATES:
            per_release = _grouped(rows, replicate)
            res = _split(per_release)
            if res is None:
                out.append({"system": name, "replicate": replicate, "measurable": False,
                            "reason": "fewer than two releases survive the replicate"})
                continue
            rec = {"system": name, "replicate": replicate,
                   "measurable": res["across_dof"] > 0, **res, **_shared(per_release)}
            if rec["measurable"]:
                rec["ratio"] = res["across_reduced"] / res["within_reduced"]
                rec["inflation"] = math.sqrt(res["across_reduced"])
                widened = res["across_chi2"] / 1.5 ** 2
                rec["reduced_at_1p5"] = widened / res["across_dof"]
                rec["p_at_1p5"] = chi2_sf(widened, res["across_dof"])
                loo = _loo_across(per_release)
                rec["loo_lo"] = min(v for v, _ in loo)
                rec["loo_hi"] = max(v for v, _ in loo)
                rec["loo_n"] = len(loo)
                rec["sigma_between"] = between_sigma(per_release)
                rec["null"] = null_calibration(per_release)
                rec["median_se"] = float(np.median([e[3] for b in per_release.values() for e in b]))
            out.append(rec)
    return out


def crosscheck() -> list[str]:
    """Replicate-0 output of this script against ``make_figHodge.release_decomposition``."""
    shipped = {r["system"]: r for r in hodge.release_decomposition()}
    notes = []
    for rec in measure():
        if rec["replicate"] != 0 or "union_chi2" not in rec:
            continue
        ref = shipped.get(rec["system"])
        if ref is None:
            notes.append(f"{rec['system']}: absent from the shipped decomposition")
            continue
        for key in ("union_chi2", "union_dof", "within_chi2", "within_dof",
                    "across_chi2", "across_dof"):
            a, b = rec[key], ref[key]
            if abs(float(a) - float(b)) > 1e-9:
                notes.append(f"{rec['system']}.{key}: {a} vs shipped {b}")
    return notes


def pooled(records: list[dict]) -> dict:
    """Across-release block pooled over the measurable systems, per replicate and overall."""
    out = {}
    for replicate in REPLICATES:
        rows = [r for r in records if r["replicate"] == replicate and r.get("measurable")]
        chi2 = sum(r["across_chi2"] for r in rows)
        dof = sum(r["across_dof"] for r in rows)
        within_chi2 = sum(r["within_chi2"] for r in rows)
        within_dof = sum(r["within_dof"] for r in rows)
        out[replicate] = {"n_systems": len(rows), "chi2": chi2, "dof": dof,
                          "reduced": chi2 / dof if dof else math.nan,
                          "p": chi2_sf(chi2, dof) if dof else math.nan,
                          "within_reduced": within_chi2 / within_dof if within_dof else math.nan,
                          "inflation": math.sqrt(chi2 / dof) if dof else math.nan}
    return out


# ------------------------------------------------------------------------------------------
# figure
# ------------------------------------------------------------------------------------------
def draw(records: list[dict], pool: dict) -> None:
    use_paper_style()
    fig, (axA, axB, axC) = plt.subplots(
        1, 3, figsize=figsize(3, 2.9), gridspec_kw={"width_ratios": [1.12, 1.0, 1.06]})

    # A: what carries the measurement at all -- shared ligands, and the cycles they create
    rep0 = sorted([r for r in records if r["replicate"] == 0 and "across_dof" in r],
                  key=lambda r: (r["across_dof"], r["max_pair"]))
    ys = np.arange(len(rep0))
    axA.barh(ys, [r["max_pair"] for r in rep0], 0.66, lw=0,
             color=[OURS if r["across_dof"] > 0 else MUTED for r in rep0])
    for y, r in zip(ys, rep0, strict=True):
        live = r["across_dof"] > 0
        axA.annotate(f"{r['across_dof']} cycle" + ("" if r["across_dof"] == 1 else "s"),
                     (r["max_pair"], y), xytext=(4, 0), textcoords="offset points",
                     va="center", ha="left", fontsize=7.5, color=OURS if live else REF)
    axA.set_yticks(ys)
    axA.set_yticklabels([r["system"] for r in rep0])
    axA.set_xlim(0, 9.8)
    axA.set_xticks([0, 2, 4, 6, 8])
    axA.set_xlabel("ligands shared between releases")
    panel(axA, "A", "where a joint cycle exists",
          subtitle="two shared ligands make one cycle")

    # B: the two blocks of the same closure statistic, three replicates each.
    # Both blocks are computed from the calibrated bar, so they are one colour family told
    # apart by weight, not two hues.
    live = sorted({r["system"] for r in records if r.get("measurable")})
    xs = np.arange(len(live), dtype=float)
    slot, w = 0.115, 0.10
    for k, rep in enumerate(REPLICATES):
        vals = {"within_reduced": [], "across_reduced": []}
        for name in live:
            rec = next((r for r in records if r["system"] == name and r["replicate"] == rep
                        and r.get("measurable")), None)
            for key in vals:
                vals[key].append(rec[key] if rec else np.nan)
        axB.bar(xs + (k - 2.6) * slot, vals["within_reduced"], w, color=tint(OURS, 0.62), lw=0,
                label="within release" if k == 0 else None)
        axB.bar(xs + (k + 0.6) * slot, vals["across_reduced"], w, color=OURS, lw=0,
                label="across releases" if k == 0 else None)
    axB.axhline(1.0, color=REF, lw=0.9, ls="--", zorder=3)
    # the only quadrant no bar enters: right of the last replicate, above the line
    axB.annotate("sampling only", (0.995, 1.0), xycoords=("axes fraction", "data"),
                 xytext=(0, 4), textcoords="offset points", ha="right", va="bottom",
                 fontsize=7.5, color=REF)
    axB.set_yscale("log")
    axB.set_ylim(0.02, 40.0)
    axB.set_yticks([0.1, 1.0, 10.0])
    axB.set_xticks(xs)
    axB.set_xticklabels(live)
    axB.set_ylabel(r"reduced $\chi^2$")
    axB.set_xlim(-0.55, len(live) - 0.45)
    panel(axB, "B", "the two blocks compared", subtitle="three replicates per system")
    legend(axB, loc="upper left")

    # C: how far the reported bar would have to be widened for the joint cycles to close
    grid = np.linspace(1.0, 4.0, 240)
    for name, style in zip(live, ("-", "--"), strict=False):
        for rep in REPLICATES:
            rec = next((r for r in records if r["system"] == name and r["replicate"] == rep
                        and r.get("measurable")), None)
            if rec is None:
                continue
            axC.plot(grid, rec["across_reduced"] / grid ** 2, style, color=OURS, lw=1.3,
                     alpha=0.8, label=name if rep == 0 else None)
    axC.axhline(1.0, color=REF, lw=0.9, ls="--", zorder=3)
    axC.axvline(1.5, color=REF, lw=0.8, ls=":", zorder=3)
    axC.annotate("bars widened\nby half", (1.55, 0.022), ha="left", va="bottom",
                 fontsize=7.5, color=REF, linespacing=1.3)
    axC.set_yscale("log")
    axC.set_xlim(1.0, 4.0)
    axC.set_ylim(0.02, 40.0)
    axC.set_yticks([0.1, 1.0, 10.0])
    axC.set_xticks([1, 2, 3, 4])
    # two lines: a single-line version of this label runs past the canvas at FULL width
    axC.set_xlabel("factor the reported bar\nis widened by")
    axC.set_ylabel(r"across-release reduced $\chi^2$")
    panel(axC, "C", "what would close them",
          subtitle=f"pooled, replicate 0: {pool[0]['inflation']:.1f}$\\times$")
    legend(axC, loc="upper right")

    small = check_min_type(fig)
    if small:
        raise AssertionError(f"type below the house floor: {small}")
    finish(fig, "figRelease_cross_release_variance")


# ------------------------------------------------------------------------------------------
# record
# ------------------------------------------------------------------------------------------
def _series(rows: list[dict], key: str, fmt: str) -> str:
    return ", ".join(format(r[key], fmt) for r in rows)


def _span(rows: list[dict], key: str, fmt: str) -> str:
    vals = [r[key] for r in rows]
    return f"{min(vals):{fmt}} to {max(vals):{fmt}}"


def _reading_signal(records: list[dict], live: list[str]) -> str:
    """The measured magnitude, per system, in one paragraph built from the records themselves."""
    parts = []
    for name in live:
        rows = _live_rows(records, name)
        parts.append(
            f"On `{name}` the across-release block runs {_series(rows, 'across_reduced', '.2f')} "
            f"over the three replicates, against a null of exactly 1, and its ratio to the "
            f"within-release block runs {_series(rows, 'ratio', '.1f')}.")
    return " ".join(parts) + (
        " Read the ratio with its denominator in mind: the within-release block is itself a "
        "finite-sample estimate and lands below 1 on some replicates, which inflates the ratio "
        "there. The across-release reduced chi-squared is the steadier read, since its null is "
        "exactly 1 by construction.")


def _reading_widening(records: list[dict], live: list[str]) -> str:
    """What the measurement says about the article's widened-bar stress test."""
    rows = _live_rows(records, live[0])
    median_se = float(np.median([r["median_se"] for r in rows]))
    return (
        f"For the widened-bar passage: multiplying every reported bar by 1.5 leaves `{live[0]}`'s "
        f"joint cycles at reduced chi-squared {_series(rows, 'reduced_at_1p5', '.2f')} "
        f"(p {_series(rows, 'p_at_1p5', '.0e')}), so a half-again wider bar does not absorb what "
        f"this measurement sees. Closing that block takes "
        f"{_span(rows, 'inflation', '.1f')}x, or an added per-edge term of "
        f"{_span(rows, 'sigma_between', '.2f')} kcal/mol against a median reported bar of "
        f"{median_se:.2f}.")


def _live_rows(records: list[dict], system: str) -> list[dict]:
    return [r for r in records if r["system"] == system and r.get("measurable")]


def _cycles_phrase(dof: int, system: str) -> str:
    return f"{dof} across-release cycle" + ("" if dof == 1 else "s") + f" on `{system}`"


def write_doc(records: list[dict], pool: dict, notes: list[str]) -> None:
    live = sorted({r["system"] for r in records if r.get("measurable")})
    dead = sorted({r["system"] for r in records
                   if "across_dof" in r and not r.get("measurable")} - set(live))
    lines = [
        "# Fig Release: the cross-release variance split, over every multi-release system",
        "",
        "Generated by `figs/make_figRelease.py` (`make figRelease`). The split is the one",
        "`release_decomposition` in `figs/make_figHodge.py` performs on replicate 0; this script",
        "imports that module, reuses its per-edge reader, and runs the split on all three",
        "replicates. Cross-check against the shipped function at replicate 0: "
        + ("no disagreement." if not notes else "DISAGREES -- " + "; ".join(notes)),
        "",
        "## What the benchmark can and cannot measure",
        "",
        "Eight of the benchmark's system names union edges from more than one separately",
        "released data set. Two shared ligands are needed before the union carries a cycle that",
        f"neither release carries alone, and only **{len(live)}** of the eight reach that:",
        "`" + "`, `".join(live) + "`. In the other "
        f"{len(dead)} (`" + "`, `".join(dead) + "`) the releases share at most one ligand and the",
        "union is disconnected or tree-joined, so no across-release cycle exists and no",
        "measurement is possible at any magnitude. **The sample size of this measurement is",
        f"n = {len(live)} systems**, carrying "
        + " and ".join(
            _cycles_phrase(next(r["across_dof"] for r in records
                                if r["system"] == s and r["replicate"] == 0), s)
            for s in live)
        + ".",
        "",
        "## The split, every multi-release system, every replicate",
        "",
        "| system | rep | releases | shared ligands | within chi2/dof | across chi2/dof | "
        "across reduced | p | ratio across/within |",
        "|---|---:|---:|---:|---|---|---:|---:|---:|",
    ]
    for r in records:
        if "across_dof" not in r:
            lines.append(f"| `{r['system']}` | {r['replicate']} | - | - | - | - | - | - | - |")
            continue
        ratio = f"{r['ratio']:.1f}x" if r.get("measurable") else "n/a"
        red = f"{r['across_reduced']:.2f}" if r.get("measurable") else "n/a"
        p = f"{r['across_p']:.1e}" if r.get("measurable") else "n/a"
        lines.append(
            f"| `{r['system']}` | {r['replicate']} | {r['releases']} | {r['max_pair']} "
            f"| {r['within_chi2']:.1f}/{r['within_dof']} "
            f"| {r['across_chi2']:.1f}/{r['across_dof']} | {red} | {p} | {ratio} |")
    lines += [
        "",
        "`shared ligands` is the largest number of ligands any two of that system's releases have",
        "in common; `ratio` is the across-release reduced chi-squared over the within-release one,",
        "which is the quantity the three passages need. A ratio of one would mean a second",
        "preparation of the same target adds nothing beyond sampling.",
        "",
        "## The measured magnitude",
        "",
        "| system | rep | across reduced | ratio | widening factor that closes it | "
        "reduced at 1.5x | p at 1.5x | added per-edge sigma that closes it | median reported se | "
        "leave-one-edge-out range |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---|",
    ]
    for r in records:
        if not r.get("measurable"):
            continue
        lines.append(
            f"| `{r['system']}` | {r['replicate']} | {r['across_reduced']:.2f} | {r['ratio']:.1f}x "
            f"| {r['inflation']:.2f}x | {r['reduced_at_1p5']:.2f} | {r['p_at_1p5']:.2e} "
            f"| {r['sigma_between']:.2f} | {r['median_se']:.2f} "
            f"| [{r['loo_lo']:.2f}, {r['loo_hi']:.2f}] over {r['loo_n']} deletions |")
    lines += [
        "",
        "The widening factor is exact rather than fitted: scaling every reported standard error by",
        "`lambda` divides the closure statistic by `lambda**2` and leaves the GLS fit unchanged,",
        "so",
        "the factor that brings the across-release block to a reduced chi-squared of one is the",
        "square root of that block's reduced chi-squared. The `1.5x` column is the same arithmetic",
        "at the widening the article's stress test applies.",
        "",
        "`added per-edge sigma` is the same statement in kcal/mol rather than as a factor: the",
        "independent per-edge term that, added in quadrature to every reported bar, brings the",
        "across-release block to a reduced chi-squared of one. It is bisected rather than solved,",
        "since adding it changes the GLS weights, and it is reported as `0.00` where the block is",
        "already at or below one, which is these data asking for no extra term at all. It is",
        "listed beside that system's median reported bar, which is what it has to be read against.",
        "",
        "## Is the null the block is read against the right one?",
        "",
        "Synthetic draws on each system's real topology and real reported bars, with one node",
        "potential per ligand shared by both releases plus independent noise of the reported size.",
        "Under that null the two preparations agree by construction, so the across-release block",
        "should sit at a reduced chi-squared of 1 and reject at 5 per cent. Fixed seed, run once.",
        "",
        "| system | rep | draws | mean reduced chi2 under the null | rate of p < 0.05 |",
        "|---|---:|---:|---:|---:|",
    ] + [
        f"| `{r['system']}` | {r['replicate']} | {r['null']['n']} "
        f"| {r['null']['mean_reduced']:.3f} | {r['null']['rate_p05']:.3f} |"
        for r in records if r.get("measurable")
    ] + [
        "",
        "## Pooled over the measurable systems",
        "",
        "| replicate | systems | across chi2/dof | reduced | p | within reduced "
        "| widening factor |",
        "|---:|---:|---|---:|---:|---:|---:|",
    ]
    for rep in REPLICATES:
        p = pool[rep]
        lines.append(
            f"| {rep} | {p['n_systems']} | {p['chi2']:.1f}/{p['dof']} | {p['reduced']:.2f} "
            f"| {p['p']:.1e} | {p['within_reduced']:.2f} | {p['inflation']:.2f}x |")
    lines += [
        "",
        "## Reading",
        "",
        "Every number above is reported; nothing was excluded and no criterion was chosen after",
        "seeing a result. The honest summary is that the benchmark supplies this measurement on",
        f"n = {len(live)} systems, that the two disagree with each other, and that one of them",
        "carries essentially all of the signal. The magnitude should be quoted with that n.",
        "",
        _reading_signal(records, live),
        "",
        _reading_widening(records, live),
        "",
        "Four caveats that limit what the number means, none of them removable from these data:",
        "",
        "1. The across-release block confounds between-preparation variance with anything else",
        "   that differs between two separately released data sets: force-field or protonation",
        "   choices, a different crystal structure, a different mapping. It is an upper bound on",
        "   the preparation term, not an isolate of it.",
        "2. The releases are not replicates of one protocol. Two data sets prepared by different",
        "   groups differ by more than a starting structure, so this bounds the replicate caveat",
        "   from above rather than resolving it.",
        "3. Both counts are small. The larger system's block rests on a handful of joint cycles",
        "   through seven shared ligands, and the leave-one-edge-out range above says how much of",
        "   it any single edge carries.",
        "4. This is not independent of anything the article already reports. `cdk8` is one of the",
        "   systems the closure test flags, and Section 4 already attributes that flag to the",
        "   cross-release join rather than to either release's simulations. The split here",
        "   quantifies that same event; it does not add a second one.",
    ]
    DOC.write_text("\n".join(lines) + "\n")
    print(f"wrote {DOC}")


def main() -> None:
    records = measure()
    notes = crosscheck()
    pool = pooled(records)
    draw(records, pool)
    write_doc(records, pool, notes)
    for r in records:
        if r.get("measurable"):
            print(f"{r['system']} rep{r['replicate']}: across {r['across_chi2']:.1f}/"
                  f"{r['across_dof']} reduced {r['across_reduced']:.2f} p {r['across_p']:.1e} "
                  f"ratio {r['ratio']:.2f} widen {r['inflation']:.2f}")
    print("cross-check:", "clean" if not notes else notes)


if __name__ == "__main__":
    main()
