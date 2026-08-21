"""Fig D — gauge-aware identifiability (Theorem 3(iv)).

Two ligand clusters (two congeneric series): dense, cheap, high-overlap edges WITHIN
each cluster; a few expensive, low-overlap BRIDGES between. The decision is top-k
WITHIN cluster A; the A-B absolute offset is a gauge / nuisance direction.

Panel A: KG = 0 (exactly) on the all-ones gauge contrast for every candidate edge;
and under the within-A decision objective, bridge / cross-cluster edges score ~0
(they cannot help rank A) while within-A edges carry the value -> acquisition is
automatically gauge-invariant.

Panel B: a gauge-UNAWARE objective (rank A and B jointly -> values resolving the A-B
offset) wastes budget on the expensive bridges; the gauge-AWARE objective spends it
within A and identifies A's top-k sooner. Ablating gauge-awareness wastes budget.

Run:  python figs/make_figD.py   (or `make figD`)
Kill (plan Fig D): if gauge-awareness gives no budget saving.
"""
from __future__ import annotations

import pathlib
import sys
import warnings

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm

warnings.filterwarnings("ignore")
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "src"))
from paperstyle import (  # noqa: E402
    MUTED,
    OURS,
    REF,
    THIRD,
    figsize,
    finish,
    legend,
    panel,
    reference_line,
    use_paper_style,
)

from bar.active import BeliefGraph, kg_scores  # noqa: E402

FIGDIR = pathlib.Path(__file__).resolve().parent
K = 4
WARMUP = 6
TOL = 0.10
# House colour rule: the gauge-aware acquisition is ours-by-construction (OURS); the
# gauge-unaware one is the named policy this article corrects, not a null (THIRD); the
# all-ones gauge direction is the null quantity the theorem sends to zero (MUTED).
C_AWARE, C_UNAWARE, C_GAUGE = OURS, THIRD, MUTED


def make_problem(n_per=14, seed=0):
    rng = np.random.default_rng(seed)
    nA, nB = n_per, n_per
    n = nA + nB
    A, B = np.arange(nA), np.arange(nA, n)
    xyA = rng.normal(size=(nA, 2))
    xyB = rng.normal(size=(nB, 2))
    dG = np.empty(n)
    dG[A] = xyA @ np.array([1.2, -0.6]) + 0.3 * rng.normal(size=nA)
    dG[B] = xyB @ np.array([1.2, -0.6]) + 0.3 * rng.normal(size=nB) + 6.0  # B offset = gauge

    edges, etype = [], []  # 'A','B','bridge'

    def knn(idx, xy, k=3):
        D = np.linalg.norm(xy[:, None] - xy[None], axis=-1)
        np.fill_diagonal(D, np.inf)
        for a in range(len(idx)):
            for b in np.argsort(D[a])[:k]:
                e = (int(idx[min(a, b)]), int(idx[max(a, b)]))
                if e not in edges:
                    edges.append(e); etype.append("A" if idx is A else "B")

    knn(A, xyA); knn(B, xyB)
    for _ in range(3):  # expensive bridges
        e = (int(rng.choice(A)), int(rng.choice(B)))
        edges.append((min(e), max(e))); etype.append("bridge")
    etype = np.array(etype)
    o = np.where(etype == "bridge", 0.2, 0.85)     # bridges: low overlap
    V = 0.05 + 0.6 * (1 - o) ** 2
    cost = np.where(etype == "bridge", 6.0, 1.0)   # bridges: expensive
    return dict(n=n, A=A, B=B, dG=dG, edges=edges, etype=etype, V=V, cost=cost)


def regret_in_A(mean, dG, A, k=K):
    sel = A[np.argsort(mean[A])[:k]]
    return float(dG[sel].mean() - np.sort(dG[A])[:k].mean())


def contrasts_within(mean, nodes, k, buffer=6):
    cand = nodes[np.argsort(mean[nodes])[: k + buffer]]
    return [(int(a), int(b)) for ia, a in enumerate(cand) for b in cand[ia + 1:]]


def weighted(mean, cov, contrasts):
    w = []
    for a, b in contrasts:
        s = np.sqrt(max(cov[a, a] + cov[b, b] - 2 * cov[a, b], 1e-12))
        w.append(float(norm.cdf(-abs(mean[a] - mean[b]) / s)))
    return np.asarray(w)


def run(objective, prob, budget, seed):
    rng = np.random.default_rng(seed + 1000)
    n, A, B, dG = prob["n"], prob["A"], prob["B"], prob["dG"]
    edges, V, cost = prob["edges"], prob["V"], prob["cost"]
    bg = BeliefGraph(n, 1e-3)
    rem = list(range(len(edges)))
    reg = [regret_in_A(bg.mean, dG, A)]
    out_of_A = 0
    for step in range(budget):
        if not rem:
            reg.append(reg[-1]); continue
        if step < WARMUP:
            pick = rem[rng.integers(len(rem))]
        else:
            cand = [edges[e] for e in rem]
            if objective == "aware":                      # rank within A only
                contr = contrasts_within(bg.mean, A, K)
            else:                                          # gauge-unaware: rank A & B jointly
                contr = contrasts_within(bg.mean, np.arange(n), K)
            w = weighted(bg.mean, bg.cov, contr)
            prec = [1.0 / V[e] for e in rem]
            cst = [cost[e] for e in rem]
            scores = kg_scores(bg, cand, prec, cst, contr, w)
            pick = rem[int(np.argmax(scores))]
        i, j = edges[pick]
        if prob["etype"][pick] != "A":
            out_of_A += 1
        y = (dG[j] - dG[i]) + rng.normal() * np.sqrt(V[pick])
        bg.add_measurement(i, j, y, 1.0 / V[pick])
        rem.remove(pick)
        reg.append(regret_in_A(bg.mean, dG, A))
    return np.array(reg), out_of_A


def main() -> None:
    use_paper_style()
    budget, seeds = 34, 24

    # ---- Panel A data: gauge KG ~ 0 exactly; bridges score ~0 under within-A decision
    prob = make_problem(seed=0)
    bg = BeliefGraph(prob["n"], 1e-3)
    rng = np.random.default_rng(1)
    for e in rng.choice(len(prob["edges"]), WARMUP, replace=False):
        i, j = prob["edges"][e]
        bg.add_measurement(i, j, prob["dG"][j] - prob["dG"][i], 1.0 / prob["V"][e])
    contr = contrasts_within(bg.mean, prob["A"], K)
    w = weighted(bg.mean, bg.cov, contr)
    cand = prob["edges"]
    aware_scores = kg_scores(bg, cand, [1.0 / v for v in prob["V"]], list(prob["cost"]), contr, w)
    # gauge (all-ones) contrast variance reduction per edge (should be ~0)
    cov = bg.cov
    ones = np.ones(prob["n"])
    gauge_red = []
    for c, (i, j) in enumerate(cand):
        b = np.zeros(prob["n"]); b[j], b[i] = 1.0, -1.0
        covb = cov @ b
        g = 1.0 / prob["V"][c]
        gauge_red.append(g * (ones @ covb) ** 2 / (1 + g * (b @ covb)))
    gauge_red = np.array(gauge_red)

    # ---- Panel B data: aware vs unaware
    curves = {"aware": [], "unaware": []}
    waste = {"aware": [], "unaware": []}
    for sd in range(seeds):
        pb = make_problem(seed=sd)
        for obj in ("aware", "unaware"):
            r, oa = run(obj, pb, budget, sd)
            curves[obj].append(r); waste[obj].append(oa)
    curves = {k: np.array(v) for k, v in curves.items()}

    fig, (axA, axB) = plt.subplots(1, 2, figsize=figsize(2, 3.15))

    et = prob["etype"]
    xmap = {"A": 0, "B": 1, "bridge": 2}
    xs = np.array([xmap[t] for t in et]) + np.random.default_rng(0).normal(0, 0.06, len(et))
    axA.scatter(xs, aware_scores, c=C_AWARE, s=24, lw=0, label="within-A decision KG", zorder=3)
    # the grey series is the one the panel is about (it is zero everywhere), so it is drawn
    # ON TOP: behind the blue it disappears entirely at within-B and at the bridges
    axA.scatter(xs, gauge_red, c=C_GAUGE, s=16, marker="x", lw=1.0,
                label="gauge (all-ones) KG", zorder=4)
    axA.set_xticks([0, 1, 2]); axA.set_xticklabels(["within-A", "within-B", "bridge"])
    axA.set_xlim(-0.45, 2.45)
    axA.tick_params(axis="x", length=0)
    axA.set_ylabel("KG / variance reduction")
    panel(axA, "A", "KG routes to relevant edges")
    legend(axA, loc="upper right")
    # the reading of the grey series, placed in the empty band above it rather than on top
    # of the within-B cluster, where the old transform-coordinate text landed
    mant, expo = f"{np.abs(gauge_red).max():.0e}".split("e")
    axA.annotate(rf"max gauge KG $= {mant}\times10^{{{int(expo)}}}$",
                 xy=(0.30, 0.12), xycoords="axes fraction",
                 ha="left", va="bottom", fontsize=7.5, color=C_GAUGE)

    x = np.arange(budget + 1)
    for obj, c in [("aware", C_AWARE), ("unaware", C_UNAWARE)]:
        m = curves[obj].mean(0); se = curves[obj].std(0, ddof=1) / np.sqrt(seeds)
        lab = "gauge-aware" if obj == "aware" else "gauge-unaware"
        axB.plot(x, m, color=c, lw=1.6, label=lab, zorder=3)
        axB.fill_between(x, m - 1.96 * se, m + 1.96 * se, color=c, alpha=0.16, lw=0)
    reference_line(axB, "vline", WARMUP)
    axB.annotate("warm-up", xy=(WARMUP, 1.0), xycoords=("data", "axes fraction"),
                 xytext=(3, -3), textcoords="offset points",
                 ha="left", va="top", fontsize=7.5, color=REF)
    axB.set_xlim(0, budget)
    axB.set_xlabel("oracle calls")
    axB.set_ylabel("regret on cluster-A top-k (kcal/mol)")
    panel(axB, "B", "gauge-awareness saves budget")
    legend(axB, loc="upper right"); axB.set_ylim(0, None)

    finish(fig, "figD_gauge_identifiability")
    print(f"\n[Panel A] max |gauge (all-ones) KG| over edges = {np.abs(gauge_red).max():.2e} (exact 0)")
    print(f"          mean within-A-decision KG: within-A={aware_scores[et=='A'].mean():.4f} "
          f"bridge={aware_scores[et=='bridge'].mean():.4f}")
    print(f"[Panel B] regret@budget: aware={curves['aware'][:,-1].mean():.3f}  "
          f"unaware={curves['unaware'][:,-1].mean():.3f}")
    print(f"          wasted measurements outside cluster A (mean): "
          f"aware={np.mean(waste['aware']):.1f}  unaware={np.mean(waste['unaware']):.1f}")


if __name__ == "__main__":
    main()
