"""Fig G — calibrated stopping (the honest Fig A -> active-learning bridge).

Fig C showed the sandwich-vs-naive WEIGHTING is second-order for ranking (GLS: the
contrast mean is unbiased for any positive weights). So where does calibration pay off
in the loop? In knowing WHEN TO STOP. Both learners run the SAME acquisition with the
SAME (correct) GLS posterior mean -- identical top-k guess at every budget. They differ
only in the uncertainty they ASSUME when deciding the top-k is resolved:

  calibrated    uses the sandwich posterior covariance (correct).
  overconfident scales the posterior se by 0.15 -- the learned-variance head's
                overconfidence measured in Fig A (~7x too small). Same number, reused.

Confidence = P(current top-k is the true top-k) estimated by Monte-Carlo over the
Gaussian posterior (gauge-safe: a global shift never changes the ranking). For the
CORRECT posterior this is calibrated -- claimed confidence tracks the real correctness
frequency. The overconfident learner's claimed confidence races ahead of reality, so it
stops too early with a top-k that is actually wrong.

Run:  python figs/make_figG.py   (or `make figG`)
Message: calibration (Fig A) pays off in the loop as trustworthy STOPPING -- independent
of the (dead) weights-efficiency claim.
"""
from __future__ import annotations

import pathlib
import sys
import warnings

import matplotlib.pyplot as plt
import numpy as np

warnings.filterwarnings("ignore")
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "src"))
from bar.active import BeliefGraph, kg_scores  # noqa: E402
from make_figC import all_pairs  # noqa: E402
from paperstyle import (  # noqa: E402
    FOIL,
    INK,
    OURS,
    REF,
    figsize,
    finish,
    legend,
    panel,
    reference_line,
    use_paper_style,
)

FIGDIR = pathlib.Path(__file__).resolve().parent
K = 4
WARMUP = 6
OC_SE = 0.15          # learned-variance head se / true se, from Fig A (~7x overconfident)
# Assumed-se sweep. What the stopping experiment measures is a FUNCTION of how wrong the
# assumed se is; printing that function is what lets a reader place a given estimator on it.
# 0.94 is the 6% under-estimate a same-budget estimator pooling single-label edges by overlap
# reaches on the identical labels (Fig A); 1.06 is the same miscalibration, safe direction.
SWEEP_SCALES = (1.06, 0.94, 0.80, 0.60, 0.40, 0.20)
STOP_CONF = 0.90
N_MC = 400
# Semantic colours (paperstyle): the calibrated learner is OURS, the counterfactually
# overconfident stand-in is the FOIL, and the actual top-k correctness is not a method at
# all -- it is the measured reality the two claimed-confidence curves are judged against,
# i.e. a REFERENCE, so it takes REF and a dashed reference line style in both panels (the
# panel-B dots are the same quantity at the stopping point). INK is text and spines only.
# REF also carries the nominal stopping level, which is the other reference here.


def make_problem(n=16, seed=0):
    """Smaller, more separable graph so the top-k is resolvable within budget."""
    rng = np.random.default_rng(seed)
    xy = rng.normal(size=(n, 2))
    dG = xy @ np.array([1.7, -0.9]) + 0.25 * rng.normal(size=n)
    D = np.linalg.norm(xy[:, None] - xy[None], axis=-1)
    np.fill_diagonal(D, np.inf)
    edges = set()
    for i in range(n):
        for j in np.argsort(D[i])[:4]:
            edges.add((min(i, int(j)), max(i, int(j))))
    while len(edges) < n * 2.2:
        i, j = rng.integers(0, n, 2)
        if i != j:
            edges.add((min(i, j), max(i, j)))
    edges = sorted(edges)
    o = np.array([np.exp(-D[i, j] / 1.5) for i, j in edges])
    V = 0.03 + 0.25 * (1.0 - o) ** 2
    cost = 1.0 + 3.0 * (1.0 - o)
    return dict(n=n, dG=dG, edges=edges, V=V, cost=cost)


def topk_confidence_mc(mean, cov, rng, k=K, se_scale=1.0):
    """MC estimate of P(argsort(posterior)[:k] == current top-k guess). Gauge-safe:
    sampled via eigh; the huge gauge (all-ones) eigen-direction is a global shift that
    leaves argsort unchanged."""
    w, U = np.linalg.eigh(se_scale**2 * cov)
    w = np.clip(w, 0, None)
    guess = frozenset(np.argsort(mean)[:k])
    z = rng.standard_normal((N_MC, len(mean)))
    samples = mean + (z * np.sqrt(w)) @ U.T
    idx = np.argpartition(samples, k, axis=1)[:, :k]
    hits = sum(frozenset(row) == guess for row in idx)
    return hits / N_MC


def run(prob, budget, seed, sweep_scales=()):
    """One trajectory: actual top-k correctness, the two shipped claimed-confidence arms,
    and one claimed-confidence arm per sweep scale.

    The two shipped arms draw from one Monte-Carlo stream and the sweep scales from a
    second, so adding or removing sweep points leaves the figure's two arms bit-identical.
    Every arm reads the same posterior; they differ only in the uncertainty the learner
    assumes when deciding the top-k is resolved.
    """
    rng = np.random.default_rng(seed + 1000)
    mcrng = np.random.default_rng(seed + 5000)
    swrng = np.random.default_rng(seed + 9000)
    n, edges, dG, V = prob["n"], prob["edges"], prob["dG"], prob["V"]
    true_top = frozenset(np.argsort(dG)[:K])
    bg = BeliefGraph(n, 1e-3)
    rem = list(range(len(edges)))
    correct, conf_c, conf_oc = [], [], []
    sweep: list[list[float]] = [[] for _ in sweep_scales]
    for step in range(budget):
        if rem:
            if step < WARMUP:
                pick = rem[rng.integers(len(rem))]
            else:
                contr = all_pairs(n)
                sc = kg_scores(bg, [edges[e] for e in rem], [1.0 / V[e] for e in rem],
                               [prob["cost"][e] for e in rem], contr, np.ones(len(contr)))
                pick = rem[int(np.argmax(sc))]
            i, j = edges[pick]
            y = (dG[j] - dG[i]) + rng.normal() * np.sqrt(V[pick])
            bg.add_measurement(i, j, y, 1.0 / V[pick])  # TRUE precision -> correct GLS posterior
            rem.remove(pick)
        m, cov = bg.mean, bg.cov
        correct.append(float(frozenset(np.argsort(m)[:K]) == true_top))
        conf_c.append(topk_confidence_mc(m, cov, mcrng, se_scale=1.0))
        conf_oc.append(topk_confidence_mc(m, cov, mcrng, se_scale=OC_SE))
        for k, sc in enumerate(sweep_scales):
            sweep[k].append(topk_confidence_mc(m, cov, swrng, se_scale=sc))
    return (np.array(correct), np.array(conf_c), np.array(conf_oc),
            [np.array(c) for c in sweep])


def main() -> None:
    use_paper_style()
    budget, seeds = 40, 60
    C, Cc, Coc = [], [], []
    Csw: list[list[np.ndarray]] = [[] for _ in SWEEP_SCALES]
    for sd in range(seeds):
        c, cc, co, sw = run(make_problem(seed=sd), budget, sd, SWEEP_SCALES)
        C.append(c); Cc.append(cc); Coc.append(co)
        for k, arr in enumerate(sw):
            Csw[k].append(arr)
    C, Cc, Coc = np.array(C), np.array(Cc), np.array(Coc)
    Csw = [np.array(c) for c in Csw]

    def stop_stats(conf):
        b, corr = [], []
        for s in range(seeds):
            hit = np.where(conf[s] >= STOP_CONF)[0]
            t = int(hit[0]) if len(hit) else budget - 1
            b.append(t); corr.append(C[s][t])
        return np.array(b), np.array(corr)

    b_c, corr_c = stop_stats(Cc)
    b_oc, corr_oc = stop_stats(Coc)

    fig, (axA, axB) = plt.subplots(1, 2, figsize=figsize(2, 3.15))
    x = np.arange(budget)
    axA.plot(x, C.mean(0), color=REF, ls=(0, (4.0, 1.8)), lw=1.3, zorder=4,
             label="actual top-k correctness")
    axA.plot(x, Cc.mean(0), color=OURS, lw=1.8, zorder=3, label="calibrated claimed conf.")
    axA.plot(x, Coc.mean(0), color=FOIL, lw=1.8, zorder=3,
             label="overconfident claimed conf.")
    reference_line(axA, "hline", STOP_CONF)
    axA.annotate(f"{STOP_CONF:.2f} stopping rule", xy=(0.4, STOP_CONF), xycoords="data",
                 xytext=(0, 3), textcoords="offset points", ha="left", va="baseline",
                 fontsize=7, color=REF)
    axA.set_xlabel("oracle calls"); axA.set_ylabel("probability")
    axA.set_xlim(0, budget - 1); axA.set_ylim(0, 1.03)
    panel(axA, "A", "claimed confidence vs reality")
    legend(axA, loc="lower right")

    xb = np.arange(2); w = 0.5
    axB.bar(xb, [b_c.mean(), b_oc.mean()], w, color=[OURS, FOIL], zorder=2)
    axB2 = axB.twinx()
    # the twin carries the right-hand scale, so it keeps its own spine; only the top goes
    axB2.spines["top"].set_visible(False)
    axB2.spines["right"].set_visible(True)
    # points, not a connected line: a connector would run across the top of the taller bar
    axB2.plot(xb, [corr_c.mean(), corr_oc.mean()], color=REF, marker="o", ls="none", ms=6,
              zorder=5)
    axB2.set_ylabel("top-k correct at stop"); axB2.set_ylim(0, 1.05)
    for i, v in enumerate([corr_c.mean(), corr_oc.mean()]):
        axB2.text(xb[i], v + 0.05, f"{v:.2f}", ha="center", fontsize=8, fontweight="bold",
                  color=INK)
    axB.set_xticks(xb); axB.set_xticklabels(["calibrated", "overconfident"])
    axB.tick_params(axis="x", length=0)
    axB.set_xlim(-0.6, 1.6)
    axB.set_ylabel("oracle calls at stop"); axB.set_ylim(0, budget)
    panel(axB, "B", "when each stops, and is it right",
          subtitle="bars = oracle calls · dots = top-k correctness")

    finish(fig, "figG_calibrated_stopping")
    print(f"wrote figG_calibrated_stopping.(pdf|png) to {FIGDIR}")
    print(f"  calibrated   : stop @ {b_c.mean():.1f} calls, top-k correct {corr_c.mean():.2f}")
    print(f"  overconfident: stop @ {b_oc.mean():.1f} calls, top-k correct {corr_oc.mean():.2f}")
    print(f"  calibration of claimed-vs-actual (mean |Δ| over budget): "
          f"calib {np.mean(np.abs(Cc.mean(0)-C.mean(0))):.3f}  oc {np.mean(np.abs(Coc.mean(0)-C.mean(0))):.3f}")

    print("  assumed-se sweep (assumed se / true se -> stop budget, top-k correct at stop, "
          "mean |claimed-actual|); * = the figure's own arms:")
    arms = list(zip(SWEEP_SCALES, Csw, strict=True)) + [(1.0, Cc), (OC_SE, Coc)]
    for sc, cf in sorted(arms, key=lambda t: -t[0]):
        b, corr = stop_stats(cf)
        star = "*" if cf is Cc or cf is Coc else " "
        print(f"   {star}{sc:5.2f}x : stop @ {b.mean():5.1f} calls, correct {corr.mean():.2f}, "
              f"|Δ| {np.mean(np.abs(cf.mean(0)-C.mean(0))):.3f}")


if __name__ == "__main__":
    main()
