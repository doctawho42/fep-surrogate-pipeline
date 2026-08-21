"""Fig C — active-learning efficiency (Theorem 3 in action).

Controlled FEP-edge-graph simulation. L ligands with hidden ΔG; a congeneric-like
edge graph with HETEROSCEDASTIC per-edge sandwich variances V_e (similar ligands ->
high overlap -> low V_e) and non-uniform costs. Oracle reveals y_e ~ N(ΔΔG_e, V_e).
Goal: identify the top-k strongest binders (lowest ΔG) in as few oracle calls as
possible. All strategies share a random warm-up (fair cold-start), then:

  ours        gauge-aware cost-aware sandwich KG; decision contrasts = top-(k+buffer)
              pairs weighted by P(misordered) (decision-theoretic boundary focus)
  uncertainty cost-aware KG over ALL pairs, uniform (A-optimal-style baseline)
  random      random unmeasured edge
  naive-KG    like 'ours' but the WRONG naive precision 1/(1/I_e)=I_e for acquisition
              AND posterior update (invariant #4 ablation: over-conservative on the
              high-overlap edges, where it matters most -- see Fig A)

Run:  python figs/make_figC.py   (or `make figC`)
Kill (plan Fig C): if 'ours' is NOT fewer oracle-calls-to-top-k than the baselines.
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
    ALT,
    INK,
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

from bar.active import BeliefGraph, bootstrap_ratio_ci, kg_scores  # noqa: E402

FIGDIR = pathlib.Path(__file__).resolve().parent
K = 5
WARMUP = 8
TOL = 0.10  # kcal/mol: "practically identified" top-k
# House colour rule: 'ours' is the calibrated sandwich-weighted policy (OURS); 'naive-KG' is
# the naive 1/I weighting this article corrects, a named rival rather than a null (THIRD);
# 'uncertainty' is a second named acquisition rule that is neither ours nor corrected (ALT);
# 'random' is the null baseline (MUTED). Nothing here is an overconfident stand-in, so FOIL
# is deliberately unused.
STRATS = {"ours": OURS, "uncertainty": ALT, "random": MUTED, "naive-KG": THIRD}
LABELS = {"ours": "ours (sandwich KG)", "uncertainty": "uncertainty-only",
          "random": "random", "naive-KG": "naive-weight KG"}


def make_problem(n_lig=40, seed=0):
    rng = np.random.default_rng(seed)
    xy = rng.normal(size=(n_lig, 2))
    dG = xy @ np.array([1.3, -0.7]) + 0.4 * rng.normal(size=n_lig)
    D = np.linalg.norm(xy[:, None] - xy[None], axis=-1)
    np.fill_diagonal(D, np.inf)
    edges = set()
    for i in range(n_lig):
        for j in np.argsort(D[i])[:4]:
            edges.add((min(i, int(j)), max(i, int(j))))
    while len(edges) < n_lig * 2.4:
        i, j = rng.integers(0, n_lig, 2)
        if i != j:
            edges.add((min(i, j), max(i, j)))
    edges = sorted(edges)
    o = np.array([np.exp(-D[i, j] / 1.5) for i, j in edges])  # overlap (0,1)
    V = 0.05 + 0.6 * (1.0 - o) ** 2          # sandwich variance: small when similar
    rho = 1.0 + 12.0 * o**2                  # naive/sandwich inflation (Fig A direction)
    cost = 1.0 + 4.0 * (1.0 - o)             # hard (low-overlap) edges cost more
    return dict(n=n_lig, dG=dG, edges=edges, V=V, V_naive=V * rho, cost=cost)


def regret(mean, dG, k=K):
    sel = np.argsort(mean)[:k]
    return float(dG[sel].mean() - np.sort(dG)[:k].mean())


def ours_contrasts(mean, cov, k=K, buffer=10):
    cand = np.argsort(mean)[: k + buffer]
    contrasts, weights = [], []
    for ia, a in enumerate(cand):
        for b in cand[ia + 1:]:
            s = np.sqrt(max(cov[a, a] + cov[b, b] - 2 * cov[a, b], 1e-12))
            contrasts.append((int(a), int(b)))
            weights.append(float(norm.cdf(-abs(mean[a] - mean[b]) / s)))  # P(misordered)
    return contrasts, np.asarray(weights)


ALL_PAIRS_CACHE: dict[int, list] = {}


def all_pairs(n):
    if n not in ALL_PAIRS_CACHE:
        ALL_PAIRS_CACHE[n] = [(a, b) for a in range(n) for b in range(a + 1, n)]
    return ALL_PAIRS_CACHE[n]


def run(strategy, prob, budget, seed):
    rng = np.random.default_rng(seed + 1000)
    n, edges, dG = prob["n"], prob["edges"], prob["dG"]
    V, V_naive, cost = prob["V"], prob["V_naive"], prob["cost"]
    belief = BeliefGraph(n, prior_precision=1e-3)
    remaining = list(range(len(edges)))
    reg = [regret(belief.mean, dG)]

    def measure(pick, strat):
        i, j = edges[pick]
        y = (dG[j] - dG[i]) + rng.normal() * np.sqrt(V[pick])
        upd = 1.0 / (V_naive[pick] if strat == "naive-KG" else V[pick])
        belief.add_measurement(i, j, y, upd)
        remaining.remove(pick)

    for step in range(budget):
        if not remaining:
            reg.append(reg[-1]); continue
        if step < WARMUP or strategy == "random":
            pick = remaining[rng.integers(len(remaining))]
        else:
            cand = [edges[e] for e in remaining]
            cost_c = [cost[e] for e in remaining]
            if strategy == "uncertainty":
                contrasts = all_pairs(n)
                weights = np.ones(len(contrasts))
                prec = [1.0 / V[e] for e in remaining]
            else:  # ours / naive-KG
                contrasts, weights = ours_contrasts(belief.mean, belief.cov)
                prec = [1.0 / (V_naive[e] if strategy == "naive-KG" else V[e]) for e in remaining]
            scores = kg_scores(belief, cand, prec, cost_c, contrasts, weights)
            pick = remaining[int(np.argmax(scores))]
        measure(pick, strategy)
        reg.append(regret(belief.mean, dG))
    return np.array(reg)


def main() -> None:
    use_paper_style()
    budget, seeds = 60, 24
    curves = {s: [] for s in STRATS}
    calls = {s: [] for s in STRATS}
    for sd in range(seeds):
        prob = make_problem(seed=sd)
        for s in STRATS:
            r = run(s, prob, budget, sd)
            curves[s].append(r)
            hit = np.where(r <= TOL)[0]
            calls[s].append(int(hit[0]) if len(hit) else budget + 1)
    curves = {s: np.array(v) for s, v in curves.items()}

    fig, (axA, axB) = plt.subplots(1, 2, figsize=figsize(2, 3.15),
                                   gridspec_kw={"width_ratios": [1.22, 1.0]})
    x = np.arange(budget + 1)
    for s in STRATS:
        m = curves[s].mean(0)
        se = curves[s].std(0, ddof=1) / np.sqrt(seeds)
        axA.plot(x, m, color=STRATS[s], lw=1.6, label=LABELS[s], zorder=3)
        axA.fill_between(x, m - 1.96 * se, m + 1.96 * se, color=STRATS[s], alpha=0.16, lw=0)
    reference_line(axA, "vline", WARMUP)
    # anchored to the axes top, not to a y-limit read before the limit is set
    axA.annotate("warm-up", xy=(WARMUP, 1.0), xycoords=("data", "axes fraction"),
                 xytext=(3, -3), textcoords="offset points",
                 ha="left", va="top", fontsize=7.5, color=REF)
    axA.set_xlim(0, budget)
    axA.set_xlabel("oracle calls (FEP edges evaluated)")
    axA.set_ylabel("simple regret (kcal/mol)")
    panel(axA, "A", "regret vs budget")
    legend(axA, loc="upper right")
    axA.set_ylim(0, None)

    order = ["ours", "uncertainty", "naive-KG", "random"]
    med = [np.median(calls[s]) for s in order]
    lo = [np.percentile(calls[s], 25) for s in order]
    hi = [np.percentile(calls[s], 75) for s in order]
    yerr = np.array([np.array(med) - np.array(lo), np.array(hi) - np.array(med)])
    axB.bar(range(len(order)), med, yerr=yerr, capsize=3,
            color=[STRATS[s] for s in order], width=0.62, zorder=2,
            error_kw={"ecolor": INK, "elinewidth": 0.9, "capthick": 0.9})
    # three of the four medians sit at the censoring value: without this line the panel
    # reads as "reached top-k at 60 calls" rather than "never reached it inside the budget"
    reference_line(axB, "hline", budget)
    axB.annotate(f"budget exhausted ({budget})", xy=(1.0, budget),
                 xycoords=("axes fraction", "data"), xytext=(0, 6),
                 textcoords="offset points", ha="right", va="bottom",
                 fontsize=7.5, color=REF)
    axB.set_xticks(range(len(order)))
    axB.set_xticklabels([LABELS[s].split(" (")[0] for s in order],
                        rotation=30, ha="right", rotation_mode="anchor")
    axB.tick_params(axis="x", length=0)
    axB.set_xlim(-0.66, len(order) - 0.34)
    axB.set_ylim(0, budget * 1.16)
    axB.set_ylabel(f"oracle calls to top-{K}\n(regret \u2264 {TOL})")
    panel(axB, "B", "calls to top-k", subtitle=f"median and IQR, {seeds} seeds")

    finish(fig, "figC_active_learning")
    print(f"\n[Panel B] median oracle-calls to top-{K} (regret<={TOL}); budget={budget}, seeds={seeds}:")
    for s in order:
        c = np.array(calls[s])
        print(f"  {LABELS[s]:24s} median {np.median(c):5.1f}  IQR [{np.percentile(c,25):.0f},{np.percentile(c,75):.0f}]"
              f"  regret@budget {curves[s][:,-1].mean():.3f}")
    om = np.median(calls["ours"])
    print(f"\n  ours vs random   : {np.median(calls['random'])/om:.2f}x fewer calls")
    print(f"  ours vs naive-KG : {np.median(calls['naive-KG'])/om:.2f}x fewer calls (invariant #4)")
    print(f"  ours vs uncert.  : {np.median(calls['uncertainty'])/om:.2f}x")

    # --- Bootstrap CI for edge-weighting ratio (reviewer S5.7) ---
    # Per-seed paired ratio: ours_calls / naive_calls.
    # Ratio > 1 means ours needs more calls (naive wins); < 1 means ours wins.
    # We report mean(ratio) with a 95% percentile bootstrap CI across the 24 seeds.
    ours_arr = np.array(calls["ours"], dtype=float)
    naive_arr = np.array(calls["naive-KG"], dtype=float)
    per_seed_ratios = ours_arr / naive_arr
    R = float(per_seed_ratios.mean())
    ci_lo, ci_hi = bootstrap_ratio_ci(per_seed_ratios, n_boot=2000, seed=0, alpha=0.05)
    print(f"\n[Edge-weighting CI — reviewer S5.7]")
    print(f"  per-seed ratio (ours/naive-KG calls): {per_seed_ratios.round(3)}")
    print(f"  mean ratio R = {R:.3f}  95% bootstrap CI [{ci_lo:.3f}, {ci_hi:.3f}]")
    print(f"  n_seeds = {seeds}, n_boot = 2000, definition: ours_calls / naive_calls (per seed, paired)")
    print(f"  Interpretation: CI brackets 1.0 → no gain from sandwich weighting (honest null)")


if __name__ == "__main__":
    main()
