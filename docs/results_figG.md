# Results — Fig G: calibrated stopping (the Fig A → active-learning bridge)

**Figure:** `figs/figG_calibrated_stopping.{pdf,png}` · **Reproduce:** `make figG`.
Reuses `src/bar/active.py` + the Fig C problem generator. Deterministic (60 seeds).

![Fig G](../figs/figG_calibrated_stopping.png)

## Why this figure exists
Fig C found the sandwich-vs-naive **weighting** is second-order for *ranking* (GLS: the
contrast posterior mean is unbiased for any positive weights), and that this is a
ceiling even for an integrated qKG. So where does the sandwich's calibration pay off in
the loop? **In knowing when to stop.** This figure is the honest Fig A → loop bridge,
and it does **not** rely on the dead weights-efficiency claim.

## Setup
A small, separable FEP-edge graph (16 ligands, top-k = 4). Both learners run the **same**
acquisition and update with the **same correct (true-precision) GLS posterior** — so
their top-k *guess* is identical at every budget. They differ only in the uncertainty
they **assume** when deciding the top-k is resolved:
- **calibrated** — the sandwich posterior covariance (correct).
- **overconfident** — the posterior se scaled by **0.15**, the learned-variance head's
  overconfidence *measured in Fig A* (~7× too small). Same number, reused.

Confidence = `P(current top-k == true top-k)` by Monte-Carlo over the Gaussian posterior
(gauge-safe: a global shift never changes the ranking). Stop at the first budget where
claimed confidence ≥ 0.90.

## Result: calibration → trustworthy stopping

| learner | claimed-vs-actual gap (mean over budget) | stop budget | top-k correct at stop |
|---|---|---|---|
| **calibrated** | **0.022** | 29.6 | **0.82** |
| overconfident | 0.137 | 14.4 | 0.60 |

- **Panel A** — the calibrated learner's claimed confidence **tracks the actual top-k
  correctness** (mean |Δ| = 0.022); the overconfident learner's claimed confidence races
  ahead of reality and crosses the 0.90 stop line while it is still only ~50–60% right.
- **Panel B** — the overconfident learner stops **2× earlier** (14 vs 30 calls) but is
  far **less likely to be correct** (0.60 vs 0.82). Calibration buys a *trustworthy stop*.

## The point
This is calibration (Fig A) paying off in the active-learning loop — through the
**stopping decision**, which depends on the *uncertainty*, not on the edge weights. It
escapes the GLS ceiling that kills the weights-efficiency story (Fig C), and together
with **Fig D** (gauge-aware routing) it is the paper's honest AL narrative:

> the BAR-bottleneck's value in the loop is *calibrated decisions* (when to stop) and
> *gauge-aware routing* (where to spend) — not a weight-tuning efficiency win.

## Honest scope
Controlled simulation; the overconfidence factor is taken from Fig A's measured
learned-variance head, not tuned here. The MC confidence is an estimate of
`P(top-k correct)`; for the correct posterior it is calibrated by construction, which is
the very property being demonstrated. A real-FEP-network version is the natural next step.

## Gate
`make check` green (41 tests; Fig G adds no `src/` code) **and** Fig G regenerable by one
command (`make figG`). Calibrated stopping is trustworthy; overconfident stopping is
premature and wrong → the AL narrative (Fig D + calibrated stopping) is supported.
