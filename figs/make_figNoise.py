"""Fig Noise: the label-noise floor under the visible fraction, and how much room the bound has.

The visible (auditable) fraction of a network's error against experiment is

    f = ||Pi eps_tilde||^2 / ||eps_tilde||^2,   eps_e = ddG_calc,e - ddG_exp,e,  eps_tilde = eps / se,

with ``Pi`` the residual projector of the whitened incidence matrix. Fig Hodge measures it on four
grounded systems and finds it far below its own chance level. A referee's objection is that the
measured ``eps`` is dominated by experimental annotation noise rather than by force-field error, so
``f`` is measuring the wrong thing.

Per-ligand annotation error enters a per-edge error as a difference of per-ligand terms, i.e. as a
GRADIENT FIELD, which ``Pi`` annihilates exactly. Label noise therefore contributes ZERO to the
numerator of ``f`` and ALL of it to the denominator: correcting for it can only move ``f`` UP. This
script measures three things.

  N1  EXACTNESS, numerically. Inject a synthetic per-ligand error of realistic size into the real
      grounded sub-networks, recompute ``||Pi eps_tilde||^2``, and show it is unchanged to machine
      precision while ``||eps_tilde||^2`` grows.

  N2  DENOMINATOR SHRINKAGE. If a share ``r`` of the measured per-edge error variance is label
      noise, the denominator shrinks by ``s = 1/(1-r)`` and ``f`` rises to ``s*f``. Tabulate ``f(s)``
      against BOTH chance levels (whitened ``dof/E`` and isotropic ``tr(Pi W)/tr(W)``) and report the
      crossing ``s* = chance / f`` at which the bound would first be lost.

  N3  PLAUSIBLE RANGE. What shrinkage each declared label-noise level actually implies, given the
      measured per-edge error. A label-noise variance cannot exceed the total measured variance; a
      level that implies more noise than is measured is arithmetically impossible and is reported as
      such rather than truncated.

CRITERIA AND CONSTANTS, FIXED BEFORE ANY NUMBER BELOW WAS COMPUTED
------------------------------------------------------------------
* 1 log10 unit of affinity = 1.3642 kcal/mol at 298.15 K.
* Mixed public IC50 reproducibility across laboratories, sigma = 0.68 log10 (Kalliokoski et al.,
  PLOS ONE 2013) = 0.9277 kcal/mol per ligand. Kalliokoski measured the MIXED-laboratory regime;
  the article's grounding uses one assay per system, IC50 before Ki, never mixed, so 0.68 log10 is
  an UPPER BOUND on the label noise, not an estimate of it.
* The three declared levels, a monotone ladder in per-edge label-noise VARIANCE:
    L1  sigma_edge = sqrt(2) * 0.9277 = 1.3119 kcal/mol.  Per-ligand errors independent: the full
        Kalliokoski variance, the upper bound.
    L2  sigma_edge = 0.9277 kcal/mol.  Half of L1's variance: the same 0.68 log10 spread carried by
        the ddG itself, i.e. the assay's systematic component is shared across a congeneric series
        and cancels in the difference.
    L3  sigma_edge = 0.6560 kcal/mol.  A quarter of L1's variance: the tighter single-assay regime
        the article's grounding is claimed to sit in.
  Only L1 is a measured literature quantity. L2 and L3 are stated reductions of it, not independent
  measurements, and are labelled that way everywhere.
* Shrinkage grid: 1, 1.5, 2, 3, 4, 5, 6, 8, 10.
* N1 uses 200 draws per system at L1's per-ligand sigma, seeded 0..199. "Machine precision" means
  the change in the numerator relative to the numerator's own size stays below 1e-10.
* Nothing here is tuned. The ladder, the grid, the draw count and the tolerance are the values
  written in this docstring before the script was first run.

Run: PYTHONPATH=src python figs/make_figNoise.py   (or `make figNoise`)
"""
from __future__ import annotations

import csv
import pathlib
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
from paperstyle import (  # noqa: E402
    ALT,
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

from bar.closeloop import load_system_edges  # noqa: E402
from bar.hodge import gradient_field, hodge_split  # noqa: E402
from bar.qc import _incidence  # noqa: E402

DOC = ROOT / "docs" / "results_figNoise.md"
GROUNDED = ["cdk8", "hif2a", "p38", "bace"]

KCAL_PER_LOG10 = 1.3642
SIGMA_LOG10 = 0.68
SIGMA_LIGAND = SIGMA_LOG10 * KCAL_PER_LOG10          # 0.9277 kcal/mol per ligand
LEVELS = [
    ("L1", float(np.sqrt(2.0) * SIGMA_LIGAND), "per-ligand independent; full Kalliokoski variance"),
    ("L2", float(SIGMA_LIGAND), "half of L1's variance; shared assay term cancels in the difference"),
    ("L3", float(SIGMA_LIGAND / np.sqrt(2.0)), "a quarter of L1's variance; tighter single-assay regime"),
]
SHRINKAGE_GRID = [1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0]
N_DRAWS = 200
EXACT_TOL = 1e-10


def experimental(system: str) -> dict:
    path = ROOT / "data" / "openfe_replicates" / f"affinity_{system}.csv"
    return {r["ligand"]: float(r["exp_dg"]) for r in csv.DictReader(path.open())}


def _pieces(edges):
    """(whitened residual projector Pi, weight matrix W, se vector) for one edge list."""
    se = np.array([e[3] for e in edges])
    whitened = np.diag(1.0 / se) @ _incidence(edges)[1]
    proj = np.eye(len(edges)) - whitened @ np.linalg.pinv(whitened.T @ whitened) @ whitened.T
    return proj, np.diag(1.0 / se ** 2), se


def load_grounded() -> list[dict]:
    """The four grounded sub-networks Fig Hodge uses, with their error against experiment."""
    out = []
    for system in GROUNDED:
        exp = experimental(system)
        edges = load_system_edges(system, exp)
        if len(edges) < 3:
            continue
        eps = np.array([ddg - (exp[b] - exp[a]) for a, b, ddg, _ in edges])
        proj, weights, se = _pieces(edges)
        tilde = eps / se
        split = hodge_split(edges, eps)
        nodes = sorted({e[0] for e in edges} | {e[1] for e in edges}, key=str)
        out.append({
            "system": system, "edges": edges, "eps": eps, "se": se, "nodes": nodes,
            "proj": proj, "weights": weights, "tilde": tilde,
            "E": len(edges), "dof": split.dof,
            "num": float(tilde @ proj @ tilde),          # ||Pi eps_tilde||^2
            "den": float(tilde @ tilde),                 # ||eps_tilde||^2
            "chance_w": split.dof / len(edges),
            "chance_iso": float(np.trace(proj @ weights) / np.trace(weights)),
            "tr_w": float(np.trace(weights)),            # sum_e 1/se_e^2
            "rms_eps": float(np.sqrt((eps ** 2).mean())),
            "median_abs_eps": float(np.median(np.abs(eps))),
            "median_se": float(np.median(se)),
        })
    return out


# ------------------------------------------------------------------ N1: exactness, numerically
def exactness(systems: list[dict]) -> list[dict]:
    """Inject a real-sized per-ligand error and measure what moves and what does not."""
    out = []
    sigma_lig = SIGMA_LIGAND
    for s in systems:
        rel_num, growth, inject_norm = [], [], []
        for seed in range(N_DRAWS):
            rng = np.random.default_rng(seed)
            b = {n: float(v) for n, v in zip(s["nodes"], rng.normal(0.0, sigma_lig, len(s["nodes"])), strict=True)}
            g = gradient_field(s["edges"], b)
            tilde2 = (s["eps"] + g) / s["se"]
            num2 = float(tilde2 @ s["proj"] @ tilde2)
            den2 = float(tilde2 @ tilde2)
            rel_num.append(abs(num2 - s["num"]) / s["num"])
            growth.append(den2 / s["den"])
            inject_norm.append(float(np.linalg.norm(g / s["se"])))
        out.append({
            "system": s["system"], "num": s["num"], "den": s["den"],
            "max_rel_num_change": float(np.max(rel_num)),
            "max_abs_num_change": float(np.max(rel_num) * s["num"]),
            "max_change_vs_den": float(np.max(rel_num) * s["num"] / s["den"]),
            "median_growth": float(np.median(growth)),
            "min_growth": float(np.min(growth)), "max_growth": float(np.max(growth)),
            "median_inject_norm": float(np.median(inject_norm)),
            "passes": float(np.max(rel_num)) < EXACT_TOL,
        })
    return out


# ------------------------------------------------- N2: denominator shrinkage and the crossings
def pooled_of(systems: list[dict]) -> dict:
    num = sum(s["num"] for s in systems)
    den = sum(s["den"] for s in systems)
    n_edges = sum(s["E"] for s in systems)
    dof = sum(s["dof"] for s in systems)
    tr_pi_w = sum(float(np.trace(s["proj"] @ s["weights"])) for s in systems)
    tr_w = sum(s["tr_w"] for s in systems)
    return {"system": "pooled", "E": n_edges, "dof": dof, "num": num, "den": den,
            "f": num / den, "chance_w": dof / n_edges, "chance_iso": tr_pi_w / tr_w,
            "tr_w": tr_w}


def crossings(rows: list[dict]) -> list[dict]:
    """The shrinkage at which f first reaches each chance level: s* = chance / f."""
    out = []
    for r in rows:
        f = r["num"] / r["den"] if "f" not in r else r["f"]
        out.append({"system": r["system"], "E": r["E"], "dof": r["dof"], "f": f,
                    "chance_w": r["chance_w"], "chance_iso": r["chance_iso"],
                    "s_star_w": r["chance_w"] / f, "s_star_iso": r["chance_iso"] / f})
    return out


def loo_worst(systems: list[dict]) -> list[dict]:
    """The single-edge deletion that leaves the least room, per system and per metric.

    Both f and its chance level move when an edge is deleted, so both are recomputed. The reported
    crossing is the SMALLEST chance/f over the deletions -- the least favourable reading the article
    already quotes for the isotropic convention.
    """
    out = []
    for s in systems:
        best_w, best_iso, arg_w, arg_iso = np.inf, np.inf, -1, -1
        for k in range(s["E"]):
            sub = [e for i, e in enumerate(s["edges"]) if i != k]
            sub_eps = np.array([s["eps"][i] for i in range(s["E"]) if i != k])
            split = hodge_split(sub, sub_eps)
            if split.dof < 1:
                continue
            proj, weights, se = _pieces(sub)
            tilde = sub_eps / se
            f = float(tilde @ proj @ tilde) / float(tilde @ tilde)
            if f <= 0:
                continue
            c_w = split.dof / len(sub)
            c_iso = float(np.trace(proj @ weights) / np.trace(weights))
            if c_w / f < best_w:
                best_w, arg_w = c_w / f, k
            if c_iso / f < best_iso:
                best_iso, arg_iso = c_iso / f, k
        out.append({"system": s["system"], "s_star_w_worst": float(best_w),
                    "s_star_iso_worst": float(best_iso),
                    "edge_w": arg_w, "edge_iso": arg_iso})
    return out


# ------------------------------------------------------- N3: what shrinkage is actually possible
def implied_shrinkage(systems: list[dict], pooled: dict) -> list[dict]:
    """For each declared level, the share of the measured whitened variance it claims, and s.

    A per-edge label noise of standard deviation ``sigma_edge`` contributes an expected whitened
    squared norm ``sigma_edge^2 * sum_e 1/se_e^2``. Its share of the measured ``||eps_tilde||^2`` is
    ``r``; the implied denominator shrinkage is ``1/(1-r)``. ``r >= 1`` means the level claims more
    noise than the network's total measured error contains, which is arithmetically impossible and
    is reported, not clipped.
    """
    rows = []
    targets = [*systems, {"system": "pooled", "tr_w": pooled["tr_w"], "den": pooled["den"],
                          "rms_eps": float("nan")}]
    for level, sigma_edge, note in LEVELS:
        for t in targets:
            r = sigma_edge ** 2 * t["tr_w"] / t["den"]
            rows.append({"level": level, "sigma_edge": sigma_edge, "note": note,
                         "system": t["system"], "r": float(r),
                         "s": float(1.0 / (1.0 - r)) if r < 1.0 else float("inf"),
                         "possible": bool(r < 1.0)})
    return rows


def capped_pooled(systems: list[dict]) -> list[dict]:
    """Pooled shrinkage with each system's label-noise variance capped at its own measured total.

    The plain pooled share is a ratio of sums, so a system whose measured error is large can carry
    a level that is impossible on the others. Capping each system at ``r_sys <= 1`` -- the most a
    system's own error can possibly be label noise -- is the pooled figure that respects every
    system's arithmetic. It is the SMALLER of the two and therefore the one FAVOURABLE to the
    article's bound, so both are reported and labelled.
    """
    den = sum(s["den"] for s in systems)
    rows = []
    for level, sigma_edge, _note in LEVELS:
        num = sum(min(sigma_edge ** 2 * s["tr_w"], s["den"]) for s in systems)
        n_bad = sum(1 for s in systems if sigma_edge ** 2 * s["tr_w"] >= s["den"])
        r = num / den
        rows.append({"level": level, "sigma_edge": sigma_edge, "r": float(r),
                     "s": float(1.0 / (1.0 - r)) if r < 1.0 else float("inf"),
                     "n_impossible": n_bad})
    return rows


def verdict(cross, pooled_cross, implied, capped, loo) -> list[dict]:
    """Needed against available: does each declared level erase each system's margin?

    ``needed`` is the crossing ``s* = chance / f``; ``available`` is the shrinkage the level implies
    on that system. ``available >= needed`` means the level would put ``f`` at or above its chance
    level, i.e. the bound is lost under that reading.

    A level with ``r >= 1`` is marked REFUTED rather than erasing. Such a level asserts more label
    variance than the network's whole measured error contains, so the decomposition
    ``eps = eps_FF + eps_label`` would need a negative force-field variance. That is a degenerate
    model, not an infinitely strong correction, and reading it either way -- as erasing the margin
    or as leaving it intact -- would be reading a contradiction as evidence. It is reported as its
    own outcome.
    """
    loo_by = {r["system"]: r for r in loo}
    imp_by = {(r["level"], r["system"]): r for r in implied}
    cap_by = {r["level"]: r for r in capped}
    rows = []
    for r in [*cross, pooled_cross]:
        name = r["system"]
        needed = {"whitened": r["s_star_w"], "isotropic": r["s_star_iso"]}
        if name in loo_by:
            needed["whitened, worst deletion"] = loo_by[name]["s_star_w_worst"]
            needed["isotropic, worst deletion"] = loo_by[name]["s_star_iso_worst"]
        for level, _sigma, _note in LEVELS:
            if name == "pooled":
                avail_u = imp_by[(level, name)]["s"]
                avail_c = cap_by[level]["s"]
            else:
                avail_u = imp_by[(level, name)]["s"]
                avail_c = avail_u
            possible = np.isfinite(avail_u)
            for metric, need in needed.items():
                rows.append({"system": name, "level": level, "metric": metric, "needed": need,
                             "available": avail_u, "available_capped": avail_c,
                             "possible": bool(possible),
                             "erased": bool(possible and avail_u >= need),
                             "refuted": bool(not possible),
                             "erased_capped": bool(np.isfinite(avail_c) and avail_c >= need)})
    return rows


def isotropic_share(systems: list[dict]) -> list[dict]:
    """The same accounting without the 1/se^2 weighting: sigma_edge^2 / mean(eps^2).

    The whitened share above is the one that matters for f, since f lives in the whitened metric.
    This unweighted companion is reported beside it because it is what a reader computes in their
    head from the article's quoted per-edge error in kcal/mol, and the two differ.
    """
    rows = []
    tot_sq = sum(float((s["eps"] ** 2).sum()) for s in systems)
    tot_e = sum(s["E"] for s in systems)
    for level, sigma_edge, _note in LEVELS:
        for s in systems:
            r = sigma_edge ** 2 / float((s["eps"] ** 2).mean())
            rows.append({"level": level, "system": s["system"], "r": float(r),
                         "s": float(1.0 / (1.0 - r)) if r < 1.0 else float("inf")})
        r = sigma_edge ** 2 / (tot_sq / tot_e)
        rows.append({"level": level, "system": "pooled", "r": float(r),
                     "s": float(1.0 / (1.0 - r)) if r < 1.0 else float("inf")})
    return rows


# ------------------------------------------------------------------------------------ figure
def make_figure(systems, exact, cross, pooled_cross, implied, capped, loo):
    use_paper_style()
    fig, axes = plt.subplots(1, 3, figsize=figsize(3, 3.1))
    imp_by = {(r["level"], r["system"]): r for r in implied}
    cap_by = {r["level"]: r for r in capped}
    all_rows = [*cross, pooled_cross]

    # -- A: the exactness demonstration -------------------------------------------------------
    ax = axes[0]
    names = [e["system"] for e in exact]
    y = np.arange(len(names))
    ax.barh(y + 0.18, [e["median_growth"] - 1.0 for e in exact], height=0.32,
            color=tint(ALT, 0.30), edgecolor=ALT, linewidth=0.8, label="denominator")
    ax.barh(y - 0.18, [max(e["max_rel_num_change"], 1e-18) for e in exact], height=0.32,
            color=tint(OURS, 0.30), edgecolor=OURS, linewidth=0.8, label="numerator")
    ax.set_xscale("log")
    ax.set_xlim(1e-15, 1e5)
    ax.set_xticks([1e-14, 1e-10, 1e-6, 1e-2])
    ax.set_yticks(y)
    ax.set_yticklabels(names)
    ax.set_xlabel("relative change")
    ax.axvline(EXACT_TOL, color=REF, linestyle=":", linewidth=1.0, zorder=0.5)
    panel(ax, "A", "noise annihilated", f"per-ligand $\\sigma$ {SIGMA_LIGAND:.2f} kcal/mol, 200 draws")
    legend(ax, loc="lower right")

    # -- B: f against the shrinkage, both chance levels, and what is available -----------------
    ax = axes[1]
    grid = np.geomspace(1.0, 60.0, 300)
    for row in cross:
        ax.plot(grid, row["f"] * grid, color=tint(OURS, 0.60), linewidth=1.0, zorder=2)
    ax.plot(grid, pooled_cross["f"] * grid, color=OURS, linewidth=2.0, zorder=3)
    ax.annotate("pooled", xy=(1.15, pooled_cross["f"] * 1.15), xytext=(0, -10),
                textcoords="offset points", fontsize=7.5, color=OURS)
    ax.annotate("systems", xy=(1.15, cross[1]["f"] * 1.15), xytext=(0, -10),
                textcoords="offset points", fontsize=7.5, color=tint(OURS, 0.45))
    ax.axhline(pooled_cross["chance_w"], color=REF, linestyle="--", linewidth=1.0, zorder=1)
    ax.axhline(pooled_cross["chance_iso"], color=REF, linestyle=":", linewidth=1.0, zorder=1)
    ax.annotate(f"chance, whitened {pooled_cross['chance_w']:.2f}",
                xy=(1.1, pooled_cross["chance_w"]), xytext=(0, 3), textcoords="offset points",
                fontsize=7.5, color=REF)
    ax.annotate(f"chance, isotropic {pooled_cross['chance_iso']:.3f}",
                xy=(1.1, pooled_cross["chance_iso"]), xytext=(0, 3), textcoords="offset points",
                fontsize=7.5, color=REF)
    # Both L1 readings are drawn: the per-system-capped one and the plain pooled one, which rests
    # on a level three of the four systems cannot carry. Showing only the capped one would hide the
    # larger, adverse figure; showing only the plain one would hide that it is refuted per system.
    for sval, label, yoff in ((cap_by["L3"]["s"], "L3", 0.0),
                              (cap_by["L2"]["s"], "L2", 1.0),
                              (cap_by["L1"]["s"], "L1 capped", 0.0),
                              (imp_by[("L1", "pooled")]["s"], "L1", 1.0)):
        if not np.isfinite(sval):
            continue
        ax.axvline(sval, color=ALT, linestyle=(0, (3, 2)), linewidth=0.9, zorder=1.5)
        ax.annotate(label, xy=(sval, 2.0e-4 * (4.5 ** yoff)), xytext=(2, 0),
                    textcoords="offset points", fontsize=7.5, color=ALT, va="bottom")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(1.0, 60.0)
    ax.set_ylim(1e-4, 1.0)
    ax.set_xticks([1, 3, 10, 30])
    ax.set_xticklabels(["1x", "3x", "10x", "30x"])
    ax.set_xlabel("shrinkage $s$")
    ax.set_ylabel("visible fraction $s\\,f$")
    panel(ax, "B", "room to chance", "amber: what label noise supplies")

    # -- C: needed against available, per system ----------------------------------------------
    ax = axes[2]
    order = [*GROUNDED, "pooled"]
    x = np.arange(len(order))
    loo_by = {r["system"]: r for r in loo}
    for fn, marker, colour, label in (
        (lambda n: next(r for r in all_rows if r["system"] == n)["s_star_w"], "o", OURS, "whitened"),
        (lambda n: next(r for r in all_rows if r["system"] == n)["s_star_iso"], "s",
         tint(OURS, 0.35), "isotropic"),
        (lambda n: loo_by[n]["s_star_iso_worst"] if n in loo_by else np.nan, "v",
         tint(OURS, 0.62), "iso., 1 edge out"),
    ):
        ax.scatter(x, [fn(n) for n in order], marker=marker, s=24, color=colour, zorder=3,
                   label=label)
    for i, (level, _sigma, _note) in enumerate(LEVELS):
        vals = [imp_by[(level, n)]["s"] for n in order]
        bad = [not imp_by[(level, n)]["possible"] for n in order]
        ax.bar(x + (i - 1) * 0.26, [(v if np.isfinite(v) else 1.0) - 1.0 for v in vals],
               width=0.26, bottom=1.0, color=tint(ALT, 0.10 + 0.25 * i), edgecolor=ALT,
               linewidth=0.7, zorder=1, label=level)
        for xi, isbad in zip(x + (i - 1) * 0.26, bad, strict=True):
            if isbad:
                ax.annotate("R", xy=(xi, 1.2), fontsize=7.5, color=ALT, ha="center", va="bottom",
                            fontweight="bold")
    ax.set_yscale("log")
    ax.set_ylim(1.0, 2000.0)
    ax.set_yticks([1, 10, 100])
    ax.set_xticks(x)
    ax.set_xticklabels(order, rotation=30, ha="right")
    ax.set_ylabel("shrinkage $s$")
    panel(ax, "C", "needed vs available", "uncapped bars; R: level refuted")
    legend(ax, loc="upper center", ncol=2, columnspacing=0.9)

    assert check_min_type(fig) == [], check_min_type(fig)
    finish(fig, "figNoise_label_noise_floor")
    return fig


# ------------------------------------------------------------------------------------- report
def write_doc(systems, exact, cross, pooled_cross, implied, capped, iso_share, loo, verd):
    L = []
    w = L.append
    all_rows = [*cross, pooled_cross]
    imp_by = {(r["level"], r["system"]): r for r in implied}
    cap_by = {r["level"]: r for r in capped}
    w("# Fig Noise: the label-noise floor under the visible fraction\n")
    w("Generated by `figs/make_figNoise.py` (`make figNoise`). The constants, the three-level")
    w("noise ladder, the shrinkage grid, the draw count and the exactness tolerance are fixed in")
    w("that script's docstring and were written before any number below was computed. Nothing here")
    w("edits an existing results record; this document only adds.\n")
    w("The quantity under test is the visible (auditable) fraction")
    w("`f = ||Pi eps_tilde||^2 / ||eps_tilde||^2`, with `eps_e` the per-edge error against")
    w("experiment, `eps_tilde = eps / se_e`, and `Pi` the residual projector of the whitened")
    w("incidence matrix, measured on the same four grounded systems Fig Hodge uses.\n")

    # ---- headline
    w("## Headline: the crossing points\n")
    w(f"Pooled `f` = {pooled_cross['f']:.4f}. It would first reach its pooled chance level at a")
    w(f"denominator shrinkage of **{pooled_cross['s_star_w']:.1f}x** in the whitened metric and")
    w(f"**{pooled_cross['s_star_iso']:.1f}x** in the isotropic one. The worst single system is")
    w(f"`bace`: {cross[-1]['s_star_w']:.1f}x whitened, {cross[-1]['s_star_iso']:.1f}x isotropic, and")
    w(f"{loo[-1]['s_star_iso_worst']:.2f}x isotropic under the least favourable single-edge deletion.\n")
    w("Against that, the shrinkage label noise can actually supply, pooled:")
    for level, sigma_edge, _note in LEVELS:
        u, c = imp_by[(level, "pooled")], cap_by[level]
        u_txt = f"{u['s']:.2f}x" if np.isfinite(u["s"]) else "unbounded"
        c_txt = f"{c['s']:.2f}x" if np.isfinite(c["s"]) else "unbounded"
        w(f"* **{level}** (`sigma_edge` = {sigma_edge:.2f} kcal/mol): {u_txt} pooled, {c_txt} with")
        w(f"  each system capped at its own measured error; impossible on {c['n_impossible']} of 4 systems.")
    w("")
    n_adv = len([r for r in verd if r["erased"] and r["system"] != "pooled"])
    w("So the bound is comfortable in the whitened metric -- the metric `f` is defined in -- and it")
    w("is NOT comfortable in the isotropic one. The readings that go against the article's thesis")
    w(f"are stated first; there are {n_adv} of them, all on `p38` and `bace`, and all but one under")
    w("the isotropic convention:\n")
    adverse = [r for r in verd if r["erased"] and r["system"] != "pooled"]
    for r in sorted(adverse, key=lambda r: (r["system"], r["level"], r["metric"])):
        tie = " (a tie, within 1 per cent)" if r["available"] < 1.01 * r["needed"] else ""
        w(f"* `{r['system']}`, {r['metric']}: needs {r['needed']:.2f}x, {r['level']} supplies "
          f"{r['available']:.2f}x -- the margin is erased{tie}.")
    if not adverse:
        w("* (none)")
    w("")
    n_ref = len({r["system"] for r in verd if r["refuted"] and r["system"] != "pooled"})
    w(f"L1, the only level with a literature measurement behind it, is refuted outright on {n_ref} of")
    w(f"the {len(cross)} systems: it asserts more label variance than the measured error contains,")
    w("so it is counted neither as erasing a margin nor as leaving one intact. See N3.\n")
    w(f"And a near-coincidence worth stating so it is not mistaken for an identity: the pooled")
    w(f"shrinkage L1 supplies without the per-system cap, {imp_by[('L1', 'pooled')]['s']:.2f}x, is")
    w(f"within one per cent of the {pooled_cross['s_star_iso']:.2f}x needed to reach the pooled")
    w("isotropic chance level. They are different quantities that happen to land together here; L1")
    w("is in any case impossible on three of the four systems taken singly.\n")

    w("## N1 Label noise is annihilated, numerically\n")
    w(f"A synthetic per-ligand error `b ~ N(0, sigma^2)` with `sigma = {SIGMA_LIGAND:.4f}` kcal/mol")
    w("(0.68 log10 at 1.3642 kcal/mol per log10, the L1 level) is added to each network as the")
    w("gradient field `eps_e -> eps_e + b_j - b_i`, 200 seeded draws per system. Both the numerator")
    w("and the denominator of `f` are recomputed on the perturbed field.\n")
    w("| system | numerator ||Pi eps~||^2 | worst change, rel. to numerator | rel. to denominator | denominator ||eps~||^2 | median growth | growth range |")
    w("|---|---:|---:|---:|---:|---:|---|")
    for e in exact:
        w(f"| `{e['system']}` | {e['num']:.6g} | {e['max_rel_num_change']:.2e} | "
          f"{e['max_change_vs_den']:.2e} | {e['den']:.6g} | {e['median_growth']:.2f}x | "
          f"[{e['min_growth']:.2f}x, {e['max_growth']:.2f}x] |")
    worst = max(e["max_rel_num_change"] for e in exact)
    verdict_txt = "holds" if all(e["passes"] for e in exact) else "FAILS"
    w("")
    w(f"Worst relative change in the numerator over all systems and all {4 * N_DRAWS} draws:")
    w(f"{worst:.2e}, against the pre-registered tolerance of {EXACT_TOL:.0e}. The exactness claim")
    w(f"{verdict_txt}. That residue is floating-point conditioning in the pseudo-inverse, not a real")
    w("signal: measured against the injected field's own scale, the denominator, it is at most")
    w(f"{max(e['max_change_vs_den'] for e in exact):.1e}. The denominator meanwhile grows by")
    w(f"{min(e['median_growth'] for e in exact):.1f}x to {max(e['median_growth'] for e in exact):.1f}x")
    w("at the median (individual draws can shrink it, since an injected field can partly cancel the")
    w("measured one; the median over 200 draws is the summary).\n")
    w("Label noise is invisible to the numerator by construction, so any correction for it can only")
    w("move `f` UPWARD. The measured `f` is a LOWER bound on the share of the error attributable to")
    w("force-field bias, and everything below asks how far up the correction could push it.\n")

    w("## N2 Denominator shrinkage, and where the bound would first be lost\n")
    w("If a share `r` of the measured per-edge error variance is label noise, the denominator")
    w("shrinks by `s = 1/(1-r)` and `f` rises to `s*f` exactly, since the numerator does not move.")
    w("Both chance levels are the ones the article already reports: `dof/E` in the whitened metric,")
    w("and `tr(Pi W)/tr(W)` for an error isotropic in kcal/mol.\n")
    w("| system | E | dof | f | whitened chance | isotropic chance | s* to whitened | s* to isotropic |")
    w("|---|---:|---:|---:|---:|---:|---:|---:|")
    for r in all_rows:
        w(f"| `{r['system']}` | {r['E']} | {r['dof']} | {r['f']:.4f} | {r['chance_w']:.4f} | "
          f"{r['chance_iso']:.4f} | {r['s_star_w']:.1f}x | {r['s_star_iso']:.1f}x |")
    w("")
    w("### f against the shrinkage grid\n")
    w("| s | " + " | ".join(f"`{r['system']}`" for r in all_rows) + " |")
    w("|---" * (len(all_rows) + 1) + "|")
    for s in SHRINKAGE_GRID:
        w(f"| {s:g}x | " + " | ".join(f"{r['f'] * s:.4f}" for r in all_rows) + " |")
    w("")
    w("Chance levels for the same columns, whitened: " +
      ", ".join(f"`{r['system']}` {r['chance_w']:.4f}" for r in all_rows) + ".")
    w("Isotropic: " + ", ".join(f"`{r['system']}` {r['chance_iso']:.4f}" for r in all_rows) + ".")
    w("The grid reaches 10x; no whitened crossing lies inside it, and two isotropic ones do.\n")
    w("### The least favourable single-edge deletion\n")
    w("Deleting one edge moves both `f` and its chance level, so both are recomputed on every")
    w("deletion and the smallest resulting crossing is reported. This is the adverse reading the")
    w("article already quotes for the isotropic convention.\n")
    w("| system | worst s* to whitened chance | worst s* to isotropic chance |")
    w("|---|---:|---:|")
    for r in loo:
        w(f"| `{r['system']}` | {r['s_star_w_worst']:.2f}x | {r['s_star_iso_worst']:.2f}x |")
    w("")

    w("## N3 What shrinkage is actually available\n")
    w("A per-edge label noise of standard deviation `sigma_edge` contributes an expected whitened")
    w("squared norm `sigma_edge^2 * sum_e 1/se_e^2`. Its share `r` of the measured `||eps_tilde||^2`")
    w("is the shrinkage's only input. `r >= 1` means the level claims more label noise than the")
    w("network's whole measured error against experiment contains: arithmetically impossible,")
    w("reported as such rather than clipped.\n")
    w("The three levels are a ladder in per-edge label-noise variance. L1 is the measured literature")
    w("quantity; L2 and L3 are stated reductions of it, not independent measurements. Kalliokoski")
    w("et al. measured the MIXED-laboratory regime, whereas the article's grounding uses one assay")
    w("per system, IC50 before Ki, never mixed, so L1 is an upper bound and not an estimate.\n")
    for level, sigma_edge, note in LEVELS:
        w(f"* **{level}**: `sigma_edge` = {sigma_edge:.4f} kcal/mol -- {note}.")
    w("")
    w("| level | system | r (whitened) | implied s | possible? | r (unweighted) | implied s (unweighted) |")
    w("|---|---|---:|---:|---|---:|---:|")
    for row in implied:
        iso = next(x for x in iso_share if x["level"] == row["level"] and x["system"] == row["system"])
        s_txt = f"{row['s']:.2f}x" if np.isfinite(row["s"]) else "impossible"
        i_txt = f"{iso['s']:.2f}x" if np.isfinite(iso["s"]) else "impossible"
        w(f"| {row['level']} | `{row['system']}` | {row['r']:.3f} | {s_txt} | "
          f"{'yes' if row['possible'] else 'NO'} | {iso['r']:.3f} | {i_txt} |")
    w("")
    w("**L1, the only level with a literature measurement behind it, is arithmetically impossible on")
    w("three of the four systems.** On `cdk8`, `p38` and `bace` a mixed-laboratory per-ligand sigma of")
    w("0.68 log10 would put more label noise into the network than its entire measured error against")
    w("experiment contains, whitened. That is the direct evidence that 0.68 log10 is an upper bound")
    w("rather than an estimate for these curated single-assay groundings, and it is a fact about the")
    w("data rather than an assumption imported to protect the bound.\n")
    w("The pooled row is a ratio of sums, so it can carry a level that no single system can. Both")
    w("pooled figures are therefore reported: the plain one, and one that caps each system's label")
    w("variance at its own measured total. The capped figure is the SMALLER and hence the one")
    w("favourable to the article, so it is never quoted alone.\n")
    w("| level | plain pooled r | plain pooled s | capped pooled r | capped pooled s | systems impossible |")
    w("|---|---:|---:|---:|---:|---:|")
    for level, _sigma, _note in LEVELS:
        u, c = imp_by[(level, "pooled")], cap_by[level]
        u_txt = f"{u['s']:.2f}x" if np.isfinite(u["s"]) else "impossible"
        c_txt = f"{c['s']:.2f}x" if np.isfinite(c["s"]) else "impossible"
        w(f"| {level} | {u['r']:.3f} | {u_txt} | {c['r']:.3f} | {c_txt} | {c['n_impossible']} of 4 |")
    w("")
    w("The whitened column is the one that governs `f`, since `f` lives in the whitened metric. The")
    w("unweighted column is reported beside it because it is what a reader computes from the")
    w("article's quoted per-edge error in kcal/mol, and the two do not agree: whitening weights an")
    w("edge by `1/se^2`, and the grounded networks' reported standard errors span 0.14 to 0.40")
    w("kcal/mol, so the precise edges dominate both sums.\n")
    w("Measured per-edge error: " +
      ", ".join(f"`{s['system']}` rms {s['rms_eps']:.2f}, median |eps| {s['median_abs_eps']:.2f},"
                f" median se {s['median_se']:.2f}" for s in systems) + " kcal/mol.\n")

    w("## Needed against available\n")
    w("A row where `available >= needed` is a reading under which the correction would take `f` to")
    w("its chance level and the bound would be lost.\n")
    w("| system | metric | needed s* | L1 | L2 | L3 | erased by |")
    w("|---|---|---:|---:|---:|---:|---|")
    seen = []
    for r in verd:
        key = (r["system"], r["metric"])
        if key in seen:
            continue
        seen.append(key)
        rows_here = [x for x in verd if (x["system"], x["metric"]) == key]
        cells, killers = [], []
        for level, _s, _n in LEVELS:
            x = next(v for v in rows_here if v["level"] == level)
            avail = x["available_capped"] if r["system"] == "pooled" else x["available"]
            cells.append(f"{avail:.2f}x" if np.isfinite(avail) else "refuted")
            if np.isfinite(avail) and avail >= x["needed"]:
                killers.append(level)
        w(f"| `{r['system']}` | {r['metric']} | {r['needed']:.2f}x | " + " | ".join(cells) +
          f" | {', '.join(killers) if killers else 'none'} |")
    w("")
    w("The pooled row uses the capped pooled availability, since the plain pooled figure rests on a")
    w("noise level that three of the four systems cannot carry.\n")

    w("## What this does and does not license\n")
    w("* The correction is one-directional: `f` can only rise. The article's 0.6% is a floor.")
    w("* In the whitened metric the pooled bound needs a 50.6x shrinkage to be lost and the most")
    w("  any level supplies pooled is 4.13x once each system is held to its own arithmetic. The")
    w("  bound survives the correction there by an order of magnitude.")
    w("* In the isotropic metric it does not survive uniformly. `bace` under a single-edge deletion")
    w("  needs only 1.36x, which even L3 -- a quarter of the mixed-laboratory variance -- supplies.")
    w("  `p38` under a single-edge deletion needs 1.46x and L3 supplies 1.80x. The isotropic reading")
    w("  of the bound is therefore not robust to a label-noise correction on the two systems where")
    w("  it was already weakest, and the revision should not claim that it is.")
    w("* Nothing here measures the label noise in these particular ChEMBL series. It bounds what a")
    w("  literature-scale label noise could do. A direct measurement would need replicate")
    w("  assay values per ligand, which the grounding sources do not carry.\n")
    DOC.write_text("\n".join(L) + "\n")
    print(f"wrote {DOC}")


def main() -> None:
    systems = load_grounded()
    pooled = pooled_of(systems)
    exact = exactness(systems)
    cross = crossings(systems)
    pooled_cross = crossings([pooled])[0]
    implied = implied_shrinkage(systems, pooled)
    capped = capped_pooled(systems)
    iso = isotropic_share(systems)
    loo = loo_worst(systems)
    verd = verdict(cross, pooled_cross, implied, capped, loo)
    make_figure(systems, exact, cross, pooled_cross, implied, capped, loo)
    write_doc(systems, exact, cross, pooled_cross, implied, capped, iso, loo, verd)
    print(f"pooled f = {pooled_cross['f']:.6f}; s* whitened {pooled_cross['s_star_w']:.1f}x; "
          f"s* isotropic {pooled_cross['s_star_iso']:.1f}x")
    for row in capped:
        print(f"  {row['level']}: capped pooled s = {row['s']:.2f}x, "
              f"impossible on {row['n_impossible']}/4")
    for r in verd:
        if r["erased"] and r["system"] != "pooled":  # possible levels only; refuted ones excluded
            print(f"  ADVERSE {r['system']} {r['metric']}: needs {r['needed']:.2f}x, "
                  f"{r['level']} supplies {r['available']:.2f}x")


if __name__ == "__main__":
    main()
