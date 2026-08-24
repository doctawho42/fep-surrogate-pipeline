"""Fig Cal -- one heavy-tailed per-edge ratio, and the four published summaries of it.

Four numbers are scattered through the manuscript and a reader loses the thread by the third:

  * ``1.41`` -- the pooled ratio of reported standard error to across-replicate spread (Fig Arep);
  * ``0.34`` -- the median reduced chi-square of the closure test over the 48 systems (Fig L);
  * ``0.853`` -- the curl-leverage-weighted mean of ``c_e^-2``, the functional the closure
    statistic identifies (Theorem D1 / Fig OOS);
  * ``0.92`` -- the calibrated scale ``se^true / se^rep``, which is the square root of the third.

The manuscript states, correctly, that these do not reduce to one scalar. This panel shows it
instead of asserting it. All four are functionals of ONE distribution: the per-edge calibration
ratio ``c_e = se_e^rep / se_e^true``, measured on this benchmark as the reported standard error
over the across-replicate standard deviation. The panel draws that distribution over the 1143
edges of the 48 systems carrying a cycle -- the same edge set as Fig Hodge, Fig OOS and check C1
of Fig Inf -- and marks each published number on it, on the ``c`` scale.

What the panel shows, and what the record quantifies:

  * the four published numbers occupy a 1.6x window while the edges they summarise span two and a
    half decades, so which functional is quoted matters more than any of their values;
  * two of the four are the SAME functional written two ways (``0.92 = sqrt(0.853)``), so there are
    three distinct functionals, not four;
  * a large part of the spread is guaranteed by the denominator alone. With three replicates
    ``s^2/sigma^2 ~ chi^2_2/2``, so even if every reported bar were exactly right the per-edge
    ratio would be distributed as ``1/sqrt(Exp(1))`` -- median 1.20, and one edge in a hundred
    above 10. That parameter-free null is drawn behind the measured histogram, and the record
    reports how much of the functional-dependence survives it.

Nothing here is a new claim about calibration. The panel is a reading aid for numbers the article
already reports, plus the audit of whether each of them reproduces from the released generators.

Design constants, fixed before any number was computed:
  * edge set: `make_figOOS.load_records` (rows complete in all three replicates, systems with at
    least 3 such edges and at least one independent cycle) -- 1143 edges over 48 systems;
  * ``c_e`` uses the manuscript's own definition, symmetric in the three replicates: the RMS of the
    per-replicate reported se over the ddof=1 sample SD of the three ddG values;
  * each of the four published numbers is recomputed by ITS OWN released recipe, not by the
    symmetric ``c_e`` above, so the audit tests the shipped generator rather than a restatement;
  * the perfect-calibration null is analytic; the null band on the power-mean grid is 400 seeded
    draws that hold each edge's reported se and redraw its replicate spread.

Run:  PYTHONPATH=src python figs/make_figCal.py   (or `make figCal`).  Deterministic (SEED below).
Data: `data/openfe_replicates/combined_pymbar4_edge_data.csv` (OpenFE IndustryBenchmarks2024,
public, 3 independent replicates per edge). No new molecular dynamics.
"""
from __future__ import annotations

import math
import pathlib
import sys

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "figs"))

from figs.make_figA_replicates import load as load_arep_edges  # noqa: E402
from figs.make_figOOS import load_records, predicted_chi2  # noqa: E402
from figs.make_figStab import REPLICATES  # noqa: E402
from paperstyle import (  # noqa: E402
    INK, MUTED, OURS, REF, check_min_type, figsize, finish, panel, reference_line,
    tint, use_paper_style,
)

FIGDIR = pathlib.Path(__file__).resolve().parent
DOC = ROOT / "docs" / "results_figCal.md"
STEM = "figCal_calibration_functionals"

SEED = 20260810      # the seed every replicate-benchmark figure in this article uses
N_NULL = 400         # perfect-calibration draws for the power-mean null band
BINS = 44            # histogram bins in log10(c_e)
LO, HI = -1.0, 2.3   # histogram support in log10(c_e); covers the measured min and max
EXPONENTS = (-2.0, -1.0, 0.0, 1.0, 2.0)   # the power-mean grid reported in the record


# --------------------------------------------------------------------------- the distribution


def edge_table():
    """Per-edge ``(c_e, s_e, se_e, h_e, system)`` over the 1143 network-fitted benchmark edges.

    ``se_e`` is the RMS of the three per-replicate reported standard errors and ``s_e`` the ddof=1
    sample SD of the three replicate ddG values, which is `make_figA_replicates.load`'s convention
    and the manuscript's definition of the measured calibration ratio. ``h_e`` is the curl-leverage
    of replicate 0's fit (Theorem D1), the weight the closure statistic applies.
    """
    recs = load_records()
    names = sorted(recs)
    c, s, se, h, sysn = [], [], [], [], []
    for name in names:
        r = recs[name]
        y = np.vstack([r["y"][k] for k in REPLICATES])
        e = np.vstack([r["se"][k] for k in REPLICATES])
        sd = np.std(y, axis=0, ddof=1)
        rms = np.sqrt(np.mean(e ** 2, axis=0))
        c.append(rms / sd)
        s.append(sd)
        se.append(rms)
        h.append(np.asarray(r["h"][0]))
        sysn.extend([name] * r["E"])
    return {
        "c": np.concatenate(c), "s": np.concatenate(s), "se": np.concatenate(se),
        "h": np.concatenate(h), "system": np.array(sysn), "recs": recs, "names": names,
    }


def power_mean(p, values, weights):
    """``(sum w v^p / sum w)^(1/p)``, with ``p = 0`` the weighted geometric mean.

    Every published summary of ``c_e`` that is an edge functional is one of these, for some
    exponent and some weight; naming the pair is what makes two of them comparable at all.
    """
    v = np.asarray(values, dtype=float)
    w = np.asarray(weights, dtype=float)
    tot = float(np.sum(w))
    if tot <= 0:
        return math.nan
    if p == 0:
        return float(np.exp(np.sum(w * np.log(v)) / tot))
    return float((np.sum(w * v ** p) / tot) ** (1.0 / p))


def weightings(c, s, h):
    """The three weight vectors the published functionals actually use.

    ``uniform`` -- every edge counts once (the typical edge);
    ``s2`` -- weight ``s_e^2``, which is what pooling in quadrature does, and the weight under
    which the ``p = +2`` mean is exactly the reported RMS over the replicate RMS;
    ``h`` -- the curl-leverage, the weight Theorem D1 puts on the closure statistic.
    """
    return {"uniform": np.ones_like(c), "s2": np.asarray(s, float) ** 2, "h": np.asarray(h, float)}


# --------------------------------------------------------------------------- the four numbers


def reproduce_four(tab):
    """Recompute each published number by its own released recipe; report, do not adjust.

    Returns a list of records with the published value, the recomputed one, the position of the
    functional on the ``c`` scale, and whether it reproduces at the precision the article prints.
    """
    recs, names = tab["recs"], tab["names"]

    # (1) 1.41 -- the pooled ratio, on Fig Arep's OWN admission rule (1145 edges, 49 systems).
    rows = load_arep_edges()
    rep = np.array([r[0] for r in rows])
    repl = np.array([r[1] for r in rows])
    keep = repl > 1e-6
    rep, repl = rep[keep], repl[keep]
    pooled_arep = math.sqrt(float(np.sum(rep ** 2) / np.sum(repl ** 2)))
    pooled_here = math.sqrt(float(np.sum(tab["se"] ** 2) / np.sum(tab["s"] ** 2)))

    # (2) 0.34 -- the median reduced chi^2 over the 48 systems, replicate 0 (Fig L's headline).
    med_chi2 = float(np.median([recs[n]["rc"][0] for n in names]))

    # (3) 0.853 -- the curl-leverage-weighted mean of c_e^-2. The released headline predicts
    #     replicate 0 from the spread of replicates 1 and 2 ONLY (independence by rotation);
    #     the in-sample variant uses all three and `results_figOOS.md` says it is never quoted.
    dof = np.array([recs[n]["dof"] for n in names], dtype=float)
    pred_indep = np.array([predicted_chi2(recs[n], 0, (1, 2)) for n in names])
    pred_insample = np.array([predicted_chi2(recs[n], 0, REPLICATES) for n in names])
    lev_indep = float(np.sum(pred_indep * dof) / np.sum(dof))
    lev_insample = float(np.sum(pred_insample * dof) / np.sum(dof))
    # the same functional computed on this panel's symmetric c_e, for the record only
    lev_symmetric = float(np.sum(tab["h"] * tab["c"] ** -2) / np.sum(tab["h"]))

    # (4) 0.92 -- the calibrated scale, the square root of (3).
    scale_indep = math.sqrt(lev_indep)

    # The diagnosis of the one mismatch is pinned in code, not left to prose: the manuscript's
    # third decimal is the in-sample variant, which the released record excludes from quotation.
    assert f"{lev_insample:.3f}" == "0.853", lev_insample
    assert f"{pooled_arep:.2f}" == "1.41", pooled_arep
    assert f"{med_chi2:.2f}" == "0.34", med_chi2
    assert f"{scale_indep:.2f}" == "0.92", scale_indep

    return {
        "pooled_arep": pooled_arep, "pooled_here": pooled_here, "n_arep": int(rep.size),
        "med_chi2": med_chi2, "lev_indep": lev_indep, "lev_insample": lev_insample,
        "lev_symmetric": lev_symmetric, "scale_indep": scale_indep,
        "marks": [
            {"key": "lev", "published": "0.853 / 0.92", "recomputed": lev_indep,
             "c": lev_indep ** -0.5,
             "what": "curl-leverage-weighted mean of $c_e^{-2}$, and its square root",
             "ok": False},
            {"key": "pooled", "published": "1.41", "recomputed": pooled_arep, "c": pooled_arep,
             "what": "pooled reported se over replicate SD", "ok": True},
            {"key": "median", "published": "0.34", "recomputed": med_chi2,
             "c": med_chi2 ** -0.5,
             "what": "median reduced $\\chi^2$ over the 48 systems", "ok": True},
        ],
    }


# --------------------------------------------------------------------------- the null


def null_density_log10(c):
    """Density of ``log10 c_e`` when every reported bar is exactly right, at three replicates.

    With ``s^2/sigma^2 ~ chi^2_2/2`` the ratio is ``c = sigma/s = 1/sqrt(E)``, ``E ~ Exp(1)``, whose
    law does not depend on ``sigma``: the null is parameter-free. Its density in ``c`` is
    ``2 c^-3 exp(-c^-2)``, and the change of variable to ``u = log10 c`` multiplies by ``c ln 10``.
    """
    c = np.asarray(c, dtype=float)
    return math.log(10.0) * 2.0 * c ** -2 * np.exp(-c ** -2)


def null_power_grid(tab, rng):
    """Power-mean grid under the perfect-calibration null: hold ``se_e``, redraw the spread.

    Answers the only question that keeps the panel honest -- how much of the disagreement between
    the four functionals is guaranteed by having three replicates rather than by real per-edge
    calibration heterogeneity.
    """
    se, h = tab["se"], tab["h"]
    out = {k: [] for k in ("uniform", "s2", "h")}
    for _ in range(N_NULL):
        s_null = se * np.sqrt(rng.chisquare(2, size=se.size) / 2.0)
        c_null = se / s_null
        w = {"uniform": np.ones_like(c_null), "s2": s_null ** 2, "h": h}
        for k in out:
            out[k].append([power_mean(p, c_null, w[k]) for p in EXPONENTS])
    return {k: np.asarray(v) for k, v in out.items()}


# --------------------------------------------------------------------------- the figure


MARK_STYLE = {"lev": (OURS, (0, (4.5, 1.4, 1.0, 1.4))), "pooled": (OURS, (0, (4, 1.6))),
              "median": (OURS, "-")}
LABEL = {
    "lev": ("$c = 1.09$   the closure functional",
            "$\\langle c_e^{-2}\\rangle_h = 0.85$; calibrated scale $0.92$"),
    "pooled": ("$c = 1.41$   the pooled bar-to-spread ratio",
               "$\\sqrt{\\sum \\mathrm{se}^2 / \\sum s^2}$, RMS over RMS"),
    "median": ("$c = 1.71$   the median closure $\\chi^2_\\nu = 0.34$",
               "inverted as $0.34^{-1/2}$; exact only if $c_e$ were constant"),
}
LABEL_Y = {"lev": 0.93, "pooled": 0.71, "median": 0.49}
ARROW_Y = {"lev": 0.88, "pooled": 0.66, "median": 0.44}


def draw(tab, four):
    use_paper_style()
    c = tab["c"]
    u = np.log10(c)
    q1, med, q3 = (float(np.percentile(c, p)) for p in (25, 50, 75))
    p1, p99 = (float(np.percentile(c, p)) for p in (1, 99))

    fig, ax = plt.subplots(figsize=figsize(1, 3.05))

    # the perfect-calibration null, as a silhouette behind the data: a reference, not a series
    grid = np.linspace(LO, HI, 600)
    ax.fill_between(10.0 ** grid, null_density_log10(10.0 ** grid), color=tint(REF, 0.80),
                    zorder=0.6, linewidth=0)
    ax.plot(10.0 ** grid, null_density_log10(10.0 ** grid), color=REF, linestyle=":",
            linewidth=1.0, zorder=0.7)

    # the histogram is built in log10(c) -- equal-width bins in the metric the axis uses -- and
    # drawn against a log x-axis, so the bar heights are a density per decade of c_e.
    edges = np.linspace(LO, HI, BINS + 1)
    dens, _ = np.histogram(u, bins=edges, density=True)
    ax.stairs(dens, 10.0 ** edges, fill=True, facecolor=tint(OURS, 0.62), edgecolor=OURS,
              linewidth=0.6, zorder=2)
    ax.set_xscale("log")

    ax.set_xlim(10.0 ** LO, 10.0 ** HI)
    ymax = float(np.max(null_density_log10(10.0 ** grid))) * 1.10
    ax.set_ylim(0.0, ymax)

    reference_line(ax, "vline", 1.0, linewidth=0.9, zorder=2.6)
    ax.annotate("$c_e = 1$", xy=(1.0, 0.0), xycoords=("data", "axes fraction"),
                xytext=(-2.0, 3.0), textcoords="offset points", ha="right", va="bottom",
                fontsize=7.5, color=REF)

    # the window the four published numbers occupy, shaded so the sliver is visible at all
    lo_mark = min(m["c"] for m in four["marks"])
    hi_mark = max(m["c"] for m in four["marks"])
    ax.axvspan(lo_mark, hi_mark, color=tint(OURS, 0.88), zorder=0.9, linewidth=0)

    for mark in four["marks"]:
        colour, dash = MARK_STYLE[mark["key"]]
        ax.axvline(mark["c"], color=colour, linestyle=dash, linewidth=1.3, zorder=3)
        head, tail = LABEL[mark["key"]]
        ax.annotate(
            f"{head}\n{tail}",
            xy=(mark["c"], ARROW_Y[mark["key"]] * ymax), xycoords="data",
            xytext=(0.505, LABEL_Y[mark["key"]]), textcoords="axes fraction",
            ha="left", va="top", fontsize=8.0, color=INK, linespacing=1.35,
            arrowprops={"arrowstyle": "-", "linewidth": 0.7, "color": INK,
                        "shrinkA": 3, "shrinkB": 1, "connectionstyle": "arc3,rad=0.12"},
        )

    # the distribution's own centre, kept subordinate: it is not one of the four
    ax.plot([med], [0.0], marker="^", markersize=5.0, color=MUTED, clip_on=False, zorder=4)
    ax.annotate(f"median $c_e = {med:.2f}$", xy=(med, 0.0), xycoords=("data", "axes fraction"),
                xytext=(3.0, 9.0), textcoords="offset points", ha="left", va="bottom",
                fontsize=7.5, color=MUTED)

    ax.set_xlabel("per-edge calibration ratio  $c_e = \\mathrm{se}^{\\mathrm{rep}}_e\\,/\\,"
                  "\\mathrm{se}^{\\mathrm{true}}_e$   (reported bar over replicate spread)")
    ax.set_ylabel("density per decade of $c_e$")
    ax.set_xticks([0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 100.0])
    ax.set_xticklabels(["0.1", "0.3", "1", "3", "10", "30", "100"])
    ax.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())

    # the two distributions are labelled where they are, not in a legend box: the leader lines
    # for the three marks own the upper right, and a box there would have to be crossed.
    ax.annotate("what perfect bars would give\nwith three replicates",
                xy=(0.60, float(null_density_log10(np.array([0.60]))[0])), xycoords="data",
                xytext=(0.015, 0.98), textcoords="axes fraction", ha="left", va="top",
                fontsize=7.5, color=REF, linespacing=1.35,
                arrowprops={"arrowstyle": "-", "linewidth": 0.7, "color": REF,
                            "shrinkA": 4, "shrinkB": 2, "connectionstyle": "arc3,rad=-0.25"})
    ax.annotate("measured", xy=(2.3, 0.32), xycoords="data", ha="left", va="bottom",
                fontsize=8.0, color=OURS)

    panel(ax, "", "four published summaries of one per-edge ratio",
          f"{c.size} edges, {len(tab['names'])} systems: the middle 98 per cent spans a factor "
          f"of {p99 / p1:.0f}, the four summaries a factor of {hi_mark / lo_mark:.2f}")

    offenders = check_min_type(fig)
    assert offenders == [], offenders
    finish(fig, STEM)
    return {"q1": q1, "med": med, "q3": q3, "p1": p1, "p99": p99,
            "lo_mark": lo_mark, "hi_mark": hi_mark, "ymax": ymax}


# --------------------------------------------------------------------------- the record


def write_doc(tab, four, shape, obs_grid, null_grid, tails):
    c = tab["c"]
    L = []
    A = L.append
    A("# Results — Fig Cal: one heavy-tailed per-edge ratio, four published summaries")
    A("")
    A(f"**Figure:** `figs/{STEM}.{{pdf,png}}` · **Reproduce:** `make figCal`")
    A("(or `PYTHONPATH=src python figs/make_figCal.py`). Deterministic")
    A(f"(`SEED = {SEED}`, {N_NULL} null draws). Edge construction, network fits and")
    A("curl-leverages are imported from `figs/make_figOOS.py`, `figs/make_figStab.py` and")
    A("`figs/make_figA_replicates.py` rather than reimplemented, so the audit below tests the")
    A("shipped generators. Data: `data/openfe_replicates/combined_pymbar4_edge_data.csv` (OpenFE")
    A("IndustryBenchmarks2024, public, 3 independent replicates per edge). No new MD.")
    A("")
    A(f"![Fig Cal](../figs/{STEM}.png)")
    A("")
    A("## What this panel is for")
    A("")
    A("Four numbers are quoted in different sections of the manuscript and a reader loses the")
    A("thread by the third: `1.41`, `0.34`, `0.853` and `0.92`. The manuscript says, correctly,")
    A("that they do not reduce to one scalar. This panel shows that instead of asserting it: all")
    A("four are functionals of ONE distribution, the per-edge calibration ratio")
    A("`c_e = se_e^rep / se_e^true`, and the panel draws that distribution with each of them")
    A("marked on the `c` scale. **The panel makes no new claim about calibration.** It is a")
    A("reading aid for numbers the article already reports, plus an audit of whether each")
    A("reproduces from the released generators.")
    A("")
    A("## The edge set and the ratio")
    A("")
    A(f"{c.size} edges over {len(tab['names'])} systems: `make_figOOS.load_records`, i.e. rows")
    A("complete in all three replicates, in systems with at least 3 such edges and at least one")
    A("independent cycle. This is the same edge set as Fig Hodge, Fig OOS and check C1 of Fig Inf,")
    A("and it is admitted by a **different rule** from Fig Arep's 1145-edge set (which needs no")
    A("cycle), so the two counts are not expected to agree and neither replaces the other.")
    A("")
    A("Per edge, in the manuscript's own definition, symmetric in the three replicates:")
    A("")
    A("```")
    A("se_e = sqrt(mean_k se_{e,k}^2)      s_e = SD_k(ddG_{e,k}) (ddof = 1)      c_e = se_e / s_e")
    A("```")
    A("")
    A(f"Quartiles {shape['q1']:.3f} / {shape['med']:.3f} / {shape['q3']:.3f}; percentiles 1 to 99")
    A(f"run {shape['p1']:.3f} to {shape['p99']:.2f}; minimum {c.min():.3f}, maximum {c.max():.1f}.")
    A(f"{100 * float(np.mean(c < 1)):.1f}% of edges carry a bar tighter than their replicate")
    A(f"spread and {100 * float(np.mean(c > 10)):.1f}% carry one more than ten times wider. The")
    A("first of those is check C6 of Fig Inf computed on this edge set; C6 itself reports 0.279 on")
    A("Fig Arep's 1145 edges, and the two are separate computations that happen to agree.")
    A("")
    A("## Do the four published numbers reproduce?")
    A("")
    A("Each was recomputed by **its own released recipe**, not by the symmetric `c_e` above, so")
    A("what is tested is the shipped generator.")
    A("")
    A("| published | what it is | recomputed | reproduces? |")
    A("|---|---|---:|---|")
    A(f"| `1.41` | pooled reported se over replicate SD, Fig Arep's {four['n_arep']}-edge set | "
      f"{four['pooled_arep']:.4f} | **yes** (1.41) |")
    A(f"| `0.34` | median reduced χ² over the {len(tab['names'])} systems, replicate 0 | "
      f"{four['med_chi2']:.4f} | **yes** (0.34) |")
    A(f"| `0.853` | curl-leverage-weighted mean of `c_e^−2` | {four['lev_indep']:.4f} | "
      "**no** — see below |")
    A(f"| `0.92` | calibrated scale `se^true/se^rep`, = √ of the row above | "
      f"{four['scale_indep']:.4f} | **yes** (0.92) |")
    A("")
    A("The last row reproduces *because the square root absorbs the mismatch in the row above*:")
    A(f"√{four['lev_indep']:.4f} and √{four['lev_insample']:.4f} both print as `0.92`. The")
    A("audit is therefore three independent checks, not four.")
    A("")
    A("### The one that does not reproduce, stated loudly")
    A("")
    A("`0.853` appears in exactly one place in this repository outside this script and this")
    A("record: `docs/paper_body.tex`. No released record and no other generator prints it.")
    A("The released generator's headline value is")
    A(f"**{four['lev_indep']:.4f}**, which `docs/results_figOOS.md` prints as `0.85`.")
    A("")
    A("The third decimal is diagnostic, and the diagnosis is asserted in the script rather than")
    A("argued in prose: the **in-sample** variant of the same functional is")
    A(f"{four['lev_insample']:.4f}, i.e. `0.853` to three decimals. That variant predicts")
    A("replicate 0 from a spread computed on all three replicates, so the predicted replicate")
    A("contributes to its own predictor; `docs/results_figOOS.md` reports it for comparison and")
    A("states it is **never the quoted number**. The headline value predicts replicate 0 from")
    A("replicates 1 and 2 only.")
    A("")
    A("The consequence is small in size and specific in place. Both variants give the calibrated")
    A(f"scale `0.92` at the printed two decimals (√{four['lev_indep']:.4f} = "
      f"{four['scale_indep']:.4f}, √{four['lev_insample']:.4f} = "
      f"{math.sqrt(four['lev_insample']):.4f}), so")
    A("the operating point, the flag counts and the `[0.79, 1.04]` interval — which is the square")
    A("root of the **independent** `[0.63, 1.08]` interval — are all unaffected. What is affected")
    A("is the derived bar width printed beside it: the manuscript says the bars are wide by")
    A("`1.08×` (= 0.853^−1/2), while the released record says **`1.09×`**")
    A(f"(= {four['lev_indep']:.4f}^−1/2 = {four['lev_indep'] ** -0.5:.4f}) with interval")
    A("`[0.96, 1.26]`.")
    A("")
    A("**No number has been adjusted here and none should be adjusted to fit this figure.** The")
    A("finding is reported for the manuscript's owner to act on: the internally consistent pair is")
    A("either (`0.85`, `1.09`) from the independent predictor, matching `results_figOOS.md`")
    A("and the quoted interval, or (`0.853`, `1.08`) from a variant the record excludes from")
    A("quotation.")
    A("The figure is drawn at the released headline value.")
    A("")
    A("## Two of the four are one functional")
    A("")
    A("`0.92 = √0.853` by construction — the manuscript derives one from the other in the same")
    A("sentence — so on the `c` scale they are a single mark. There are **three distinct")
    A("functionals**, not four, and the panel draws three lines. Their positions on the `c` scale:")
    A("")
    A("| functional | as published | on the `c` scale |")
    A("|---|---|---:|")
    A(f"| curl-leverage mean of `c^−2` (and its root, the scale) | 0.853 / 0.92 | "
      f"{four['lev_indep'] ** -0.5:.3f} |")
    A(f"| pooled reported se over replicate SD | 1.41 | {four['pooled_arep']:.3f} |")
    A(f"| median reduced χ², inverted | 0.34 | {four['med_chi2'] ** -0.5:.3f} |")
    A("")
    A("The median-χ² entry is the one that needs a warning label, and the manuscript already")
    A("carries it: inverting a median recovers a scalar conservatism only if `c_e` were constant")
    A("across edges, which is exactly what this panel shows it is not. It is also not an edge")
    A("functional at all — it is a median over 48 **system** statistics, each of them itself a")
    A("leverage-weighted mean of `c_e^−2` within its own system — which is why it does not appear")
    A("in the power-mean grid below.")
    A("")
    A("## Why the axis is logarithmic")
    A("")
    A("The panel's job is to make the heavy tail visible, and a linear axis cannot do it here.")
    A("")
    A("1. `c_e` is a ratio, so its natural metric is multiplicative: a bar twice too wide")
    A("   (`c = 2`) and one twice too tight (`c = 0.5`) are equally miscalibrated in opposite")
    A("   directions, and only a logarithmic axis places them symmetrically about 1.")
    A(f"2. The support runs from {c.min():.3f} to {c.max():.1f}. On a linear axis over that range")
    A("   the middle half of the edges would occupy "
      f"{100 * (shape['q3'] - shape['q1']) / c.max():.1f}%")
    A("   of the axis width and the three functional marks would sit within")
    A(f"   {100 * (shape['hi_mark'] - shape['lo_mark']) / c.max():.1f}% of it — one line, not")
    A("   three.")
    A("3. The tail is not a nuisance to be clipped: it is the *reason* the functionals disagree,")
    A("   since each of them weights it differently. No edge is dropped, no axis break is used and")
    A("   the histogram support covers the measured minimum and maximum.")
    A("")
    A("On the log axis the four published summaries occupy a")
    A(f"{shape['hi_mark'] / shape['lo_mark']:.2f}× window while the edges span a factor of")
    A(f"{shape['p99'] / shape['p1']:.0f} between their 1st and 99th percentiles. The")
    A(f"distribution's own median, {shape['med']:.2f}, sits above all three marks; it is drawn")
    A("subordinate because")
    A("it is not one of the quoted four (it is published in `docs/results_figInf.md` as the median")
    A("per-edge ratio 1.81×).")
    A("")
    A("## How much of the spread is guaranteed by three replicates")
    A("")
    A("This is the caveat that keeps the panel honest, and it is drawn as the grey silhouette.")
    A("The denominator of `c_e` is a sample SD on three replicates, so `s²/σ² ~ χ²₂/2` and even if")
    A("every reported bar were exactly right the per-edge ratio would be `1/√Exp(1)`, whose law")
    A("does not depend on σ. That null is parameter-free, and heavy in its own right:")
    A("")
    A("| | measured | perfect bars, n = 3 |")
    A("|---|---:|---:|")
    A(f"| median `c_e` | {shape['med']:.2f} | {1 / math.sqrt(math.log(2)):.2f} |")
    A(f"| fraction below 1 | {tails['obs_lt1']:.3f} | {math.exp(-1.0):.3f} |")
    A(f"| fraction above 10 | {tails['obs_gt10']:.3f} | {1 - math.exp(-0.01):.3f} |")
    A(f"| percentiles 1 to 99 | {shape['p1']:.2f} to {shape['p99']:.1f} | "
      f"{(math.log(100.0)) ** -0.5:.2f} to {(-math.log(0.99)) ** -0.5:.1f} |")
    A("")
    A("So a reader who sees the tail and concludes that per-edge calibration is wildly")
    A("heterogeneous is going too far: a good part of that tail is the denominator's own small-")
    A("sample noise. The measured distribution is nonetheless clearly wider than the null on every")
    A("row, and the next section quantifies what that does to the four functionals.")
    A("")
    A("## The functional-dependence, quantified")
    A("")
    A("Every one of the published edge functionals is a weighted power mean")
    A("`M_p = (Σ w c^p / Σ w)^(1/p)` for some exponent and some weight; naming the pair is what")
    A("makes them comparable. The pooled ratio is `M_+2` under weight `s_e²` (that weighting is")
    A("what quadrature pooling does, and it is exactly `√(Σ se² / Σ s²)`); the closure functional")
    A("is `M_−2` under the curl-leverage weight `h_e`. The grid below is the same 1143 edges under")
    A("every combination, with the perfect-calibration null in brackets (median of")
    A(f"{N_NULL} draws that hold each edge's reported se and redraw its spread).")
    A("")
    A("| weight | `M_−2` | `M_−1` | `M_0` | `M_+1` | `M_+2` |")
    A("|---|---:|---:|---:|---:|---:|")
    for key, name in (("uniform", "uniform (the typical edge)"),
                      ("s2", "`s_e²` (quadrature pooling)"),
                      ("h", "`h_e` (curl-leverage, the closure weight)")):
        cells = " | ".join(f"{obs_grid[key][i]:.2f} [{np.median(null_grid[key][:, i]):.2f}]"
                           for i in range(len(EXPONENTS)))
        A(f"| {name} | {cells} |")
    A("")
    A("Read it as follows. **If the ratio were the same on every edge, every cell would hold the")
    A("same number** and the four published summaries would agree exactly. Instead the measured")
    A(f"cells range from {obs_grid['s2'][0]:.2f} to {obs_grid['h'][4]:.2f}, a factor of")
    A(f"{obs_grid['h'][4] / obs_grid['s2'][0]:.0f}, and the two cells the manuscript actually")
    A(f"quotes — `s_e²` at `M_+2` ({obs_grid['s2'][4]:.2f} = the published 1.41) and `h_e` at")
    A(f"`M_−2` ({obs_grid['h'][0]:.2f}) — differ by a third. That second cell is this panel's")
    A("symmetric-`c_e` estimator, not the released rotation recipe, which gives")
    A(f"{four['lev_indep'] ** -0.5:.3f}; its agreement to two decimals with the manuscript's")
    A("printed")
    A("`1.08` is a rounding coincidence between two different estimators and is not a reproduction")
    A("of it. The disagreement between the manuscript's numbers is therefore a")
    A("property of the ratio distribution, not a discrepancy between measurements, which is what")
    A("the manuscript says and what this panel now shows.")
    A("")
    A("Two further readings the grid supports, both stated as bounds rather than corrections:")
    A("")
    A("- Under the null the same grid already ranges from")
    A(f"  {np.median(null_grid['s2'][:, 0]):.2f} to "
      f"{np.median(null_grid['uniform'][:, 4]):.2f}, so")
    A("  **part of the functional-dependence is guaranteed by n = 3 alone** and only the excess is")
    A("  per-edge heterogeneity.")
    A("- The `M_+2` cell under uniform weight is the fragile one: under the null the")
    A("  expectation of")
    A("  `c²` does not exist (`E[1/E]` diverges for `E ~ Exp(1)`), so that cell is a finite-sample")
    A("  artefact of the denominator and is not a quantity to quote. The published pooled ratio")
    A("  avoids exactly this, because the `s_e²` weight cancels the `1/s²` blow-up — which is why")
    A("  its null cell sits at a well-behaved")
    A(f"  {np.median(null_grid['s2'][:, 4]):.2f} rather than running away.")
    A("")
    A("## What this licenses, and what it does not")
    A("")
    A("- It licenses the manuscript's sentence that the four numbers are different functionals of")
    A("  one heavy-tailed distribution and do not reduce to one scalar: that is now shown, with")
    A("  the window they occupy and the spread they summarise both measured.")
    A("- It does **not** license picking a preferred scalar. The panel is a reason not to.")
    A("- It does **not** re-measure calibration. Every value here is recomputed from the released")
    A("  generators, and the one that does not reproduce is reported rather than adjusted.")
    A("- The restriction that governs every calibration statement in this article governs this one")
    A("  too: the three OpenFE repeats share starting coordinates, so the denominator is a lower")
    A("  bound on run-to-run reproducibility and every `c_e` here is an upper bound on the")
    A("  conservatism a preparation-resampling replicate set would show.")
    A("")
    DOC.write_text("\n".join(L) + "\n")


# --------------------------------------------------------------------------- main


def main():
    tab = edge_table()
    four = reproduce_four(tab)
    shape = draw(tab, four)

    w = weightings(tab["c"], tab["s"], tab["h"])
    obs_grid = {k: np.array([power_mean(p, tab["c"], w[k]) for p in EXPONENTS]) for k in w}
    rng = np.random.default_rng(SEED)
    null_grid = null_power_grid(tab, rng)
    tails = {"obs_lt1": float(np.mean(tab["c"] < 1.0)), "obs_gt10": float(np.mean(tab["c"] > 10.0))}

    write_doc(tab, four, shape, obs_grid, null_grid, tails)
    print(f"figCal: {tab['c'].size} edges / {len(tab['names'])} systems; "
          f"c quartiles {shape['q1']:.2f}/{shape['med']:.2f}/{shape['q3']:.2f}, "
          f"1-99% {shape['p1']:.2f}-{shape['p99']:.1f}")
    print(f"  pooled ratio      published 1.41  recomputed {four['pooled_arep']:.4f}  OK")
    print(f"  median chi2_nu    published 0.34  recomputed {four['med_chi2']:.4f}  OK")
    print(f"  leverage mean c-2 published 0.853 recomputed {four['lev_indep']:.4f}  MISMATCH "
          f"(in-sample variant = {four['lev_insample']:.4f} = 0.853)")
    print(f"  calibrated scale  published 0.92  recomputed {four['scale_indep']:.4f}  OK")
    print(f"  wrote {DOC}")


if __name__ == "__main__":
    main()
