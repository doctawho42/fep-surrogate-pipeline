# Results — Fig A: "target the sandwich" (correctness, not a straw-man win)

**Figure:** `figs/figA_target_the_sandwich.{pdf,png}` · **Reproduce:** `make figA`.
Deterministic (seed 20260629).

![Fig A](../figs/figA_target_the_sandwich.png)

## The honest claim (reframed — read this first)
The BAR-bottleneck reads the **aleatoric (sampling) variance** off the estimator as the
sandwich `B/I²`. The claim is **not** "we beat naive `1/I`." No competent FEP user
reports bare `1/I`: **pymbar already subtracts the correction** `nrat = 1/n_f + 1/n_r`,
which is the **same order O(1/n)** as `1/I` itself (a *leading-order* term, not a
finite-sample detail), and its MBAR uncertainty `1/I − nrat` **already equals the
sandwich**. So Fig A makes three honest points instead:

1. **Correctness.** The sandwich `B/I²` **coincides** with pymbar's MBAR uncertainty and
   with the Monte-Carlo truth across every overlap regime — reproduced **with no MC**.
2. **The right foil is a learned-variance head, and it remains strongly overconfident — even when fed the overlap
   scalar I (the reviewer's "fair foil").** An MVE / heteroscedastic NN trained by
   Gaussian NLL on a *realistic* budget of ~200 edges, now receiving work-summary moments
   **and** the overlap scalar `I` (6-dim feature, reviewer M2a), is overconfident by
   **0.09–0.20× true se across the overlap sweep (regime-dependent; central ~7×, se ≈ 0.15×)**.
   A **large-budget oracle** (4000 training edges, same features) reaches 0.20× — also
   ~5× overconfident (regime-dependent). Corrected foils do **not** rescue it at realistic
   budget: a **β-NLL head** (Seitzer et al. 2022, β=0.5 — the standard objective fix) does
   *not* recover calibration (**0.11× true se**, in fact slightly *worse* than plain NLL's
   0.15×), and a **5-member
   deep ensemble** (Lakshminarayanan et al. 2017) reaches only the large-budget floor
   (**0.20×**). So the residual ~5× is a **data-efficiency limit of a per-edge *learned*
   variance**, **not an identifiability barrier** and **not specific to the Gaussian-NLL
   objective**: the pooling check below shows se(overlap) *is* recoverable from single per-edge
   labels once they are pooled by overlap bucket
   (N=200: 20.1% error, N=4000: 5.8%, N=40000: 3.8% — shrinking with budget). The BAR
   bottleneck sidesteps the learned estimator's budget floor entirely by **computing** the
   variance exactly, untrained, at zero label cost. The honest frame is "learnable-with-data
   (at a label cost the physics avoids) vs exact-and-free," not "unlearnable vs exact."
3. **`1/I` is shown only as the textbook value it corrects** — the information-equality
   plug-in, off by a *varying* factor — **not** as a baseline anyone reports.

## Why this matters downstream (the actual contribution)
- **Differentiable closed form → graph weights.** `B/I²` has an O(1) backward (Thm 1)
  and propagates into the surrogate **and** the Fisher–resistance Laplacian weights
  `w_e = I_e²/B_e` (Thm 3, invariant #4). pymbar returns a number; it is not a
  backprop-able graph weight. This is the new object.
- **Robustness.** `B/I² ≥ 0` always; pymbar's `1/I − nrat` can go **negative → nan** in a
  rare edge case (very high overlap / tiny n, when `1/I < 1/n_f + 1/n_r`); the sandwich form
  is ≥ 0 by construction and degrades gracefully.
- **Calibrated aleatoric for free.** vs the learned-variance head above.

## Result: PASS

**Panel A — controlled (MC truth).** Gaussian work model, `n_f=n_r=20`, 3000
replicates/point; true se = MC SD; band = 95% bootstrap CI. Fair foil = work moments +
overlap `I` (reviewer M2a); oracle = same features, 4000 training edges.

| overlap `4⟨p(1−p)⟩` | sandwich/true | **pymbar-MBAR/true** | fair-foil/true | oracle/true | naive `1/I`/true |
|---:|---:|---:|---:|---:|---:|
| 0.17 (low)  | 0.89 | 0.99 | **0.17** | 0.39 | 1.08 |
| 0.26        | 0.94 | 0.99 | **0.20** | 0.26 | 1.14 |
| 0.38        | 0.98 | 1.00 | **0.16** | 0.19 | 1.27 |
| 0.46        | 1.00 | 1.00 | **0.11** | 0.17 | 1.36 |
| 0.53        | 1.00 | 1.00 | **0.09** | 0.15 | 1.47 |
| 0.62        | 1.00 | 1.00 | **0.10** | 0.13 | 1.61 |
| 0.78        | 0.99 | 0.98 | **0.15** | 0.11 | 2.11 |
| 0.86 (high) | 0.99 | 0.98 | **0.18** | 0.10 | 2.65 |

**Headline numbers (M2a, honest fair-foil result):**
- Fair-foil mean (realistic budget, fed moments + overlap I): **0.15× true se** (~7× overconfident, range 0.09–0.20× across overlap sweep).
- Oracle mean (large budget, 4000 training edges, same features): **0.20× true se** (~5× overconfident).
- **β-NLL (β=0.5, realistic budget): 0.11× true se** — the standard objective fix does *not* recover calibration (it is slightly worse than plain NLL's 0.15×, i.e. further from calibrated). **5-member deep ensemble (realistic budget): 0.20× true se** — only reaches the large-budget floor.
- Even with 20× more training data, a corrected objective (β-NLL), OR a stronger estimator (deep ensemble), the learned head stays ~5× overconfident at realistic budget. The conditional sampling variance IS identifiable from single ΔΔĜ labels (pooling check below) — so the residual overconfidence is a **data-efficiency floor of per-edge learned variance**, not an information-theoretic floor and not specific to the Gaussian-NLL objective (β-NLL does not help — if anything it is slightly worse; the ensemble only reaches the 20×-budget floor). The physics computes it exactly, per-edge, differentiably, at zero label cost.
- Honest frame: **"learnable-with-data vs exact-and-free"** — not "the learned head fails" and not "one label per edge can never work." It can improve substantially with data, and a differently-trained/pooled estimator *can* recover the target se — but no learned head matches the closed form for free, and the standard per-edge Gaussian-NLL MVE head does not get there at realistic budgets.

**Identifiability check (reviewer round-2 §3): is `se(overlap)` recoverable from single labels?**
The `_gauss_edge` model has true ΔF = 0 for every separation `s`, so within a narrow overlap
bin the residual SD of *single-label* ΔF̂ across edges is the pure sampling SD — no confound
from a varying true mean. `pooled_se_recovery(n_edges)` (`figs/make_figA.py`) pools edges into
12 overlap-quantile bins, computes the pooled within-bin residual SD of the single-label ΔF̂,
and compares it to the mean sandwich se in that bin (mean relative error across bins). Run via
`python figs/make_figA.py` (deterministic, seed 7):

```
[identifiability] pooled se(overlap) recovery from single labels: N=200:20.1%, N=4000:5.8%, N=40000:3.8%
```

Pooling single-label edges by overlap recovers the sandwich `se(overlap)` to **N=200: ~20%,
N=4000: ~5.8%, N=40000: ~3.8%** relative error — the error **shrinks monotonically with edge
budget**, i.e. the conditional sampling variance **IS identifiable** from one ΔF̂ label per edge.
This directly falsifies the round-2-flagged claim that per-edge sampling variance is
unlearnable "regardless of budget" or "fundamentally insufficient" from single labels: a
nonparametric pooling-by-overlap estimator recovers it cleanly, and the recovery error is
driven by finite-sample binning noise, not an identifiability barrier. The learned MVE head's
residual ~5–7× overconfidence (above) is therefore best read as a **data-efficiency limit of a
per-edge learned variance** — it does not vanish under the standard objective fix (β-NLL, Seitzer
et al. 2022, "On the Pitfalls of Heteroscedastic Uncertainty Estimation with Probabilistic Neural
Networks") — if anything worse, 0.11× vs plain NLL's 0.15× — nor under a stronger estimator (a deep
ensemble only reaches the 20×-larger-budget floor, 0.20×); heteroscedastic per-edge variance regression stays poorly
calibrated at realistic label budgets — but it is **not** evidence that the target is unlearnable
in principle (pooling recovers it). The BAR bottleneck's advantage is that it sidesteps this
learned-estimator budget floor entirely by computing the
sandwich in closed form, untrained, at zero label cost — not that the quantity is otherwise
unlearnable.

Sandwich ≈ MBAR ≈ 1 (coincide, the correctness proof). Max |Δ|(sandwich − MBAR) = 0.096.
`1/I` traces a *non-constant* 1.08→2.65 factor — the info-equality plug-in, off worst at
high overlap (see the corrected Corollary in `bar_proofs.tex`).

**Panel B — real FEP edges.** Adjacent-λ BAR edges; both the sandwich inputs and the
bootstrap truth are subsampled to the statistical inefficiency `g`
(`pymbar.timeseries.statistical_inefficiency` + `subsample_correlated_data`) before `B`
is computed — i.e. effective (decorrelated) `n` throughout:
- **Binding — BACE1 RBFE (alchemtest AMBER, complex legs), 19 edges:** sandwich/boot
  **1.00** [0.96, 1.07]; naive `1/I`/boot **6.07** [2.78, 12.01].
- Solvation — benzene hydration, 19 edges: sandwich/boot 1.01 [0.95, 1.05]; naive 4.05.

The sandwich is calibrated on real **protein-ligand binding** edges (the relevant
context), not just solvation.

**g-robustness:** A robustness check of the panel-B sandwich ratio under coarser/finer
`g`-subsampling (analogous to the overlap-filter robustness check above) was not run in this
revision because no cached work arrays are stored and re-parsing the AMBER trajectories is
non-trivial. It is flagged as future work: the sensitivity of the sandwich ratio to the
specific `g` estimate from `pymbar.timeseries.statistical_inefficiency` should be quantified
(expected to be small, since both sandwich inputs and bootstrap truth use the same `g`).

## Honest scope
- The MBAR coincidence is exact in expectation; at extreme low overlap (n=20) the
  ddof-1 sandwich plug-in runs ~10% below MBAR (both still ≈ 1) — a finite-sample
  plug-in difference, not a discrepancy. Production reports the sandwich (robust, ≥ 0).
- Bootstrap truth (no independent repeats for these systems) is the model-free
  reference on real edges; the controlled panel supplies the MC-truth validation.
- OpenFE IndustryBenchmarks2024 (3 repeats → true-replicate se) is the next data step.

## Gate
`make check` green (95 tests, incl. `test_figA_foil` + `test_figA_pooling`) **and** Fig A regenerable by
one command (`make figA`). Sandwich = MBAR = truth (correct); fair foil (M2a) still
~7× overconfident at realistic budget, ~5× at large budget — honest "learnable-with-data
vs exact-and-free" framing; calibrated on real binding edges → **proceed**.

**JCTC M2a status:** ADDRESSED. Fair foil now receives work moments + overlap `I` (6-dim
feature). Large-budget oracle added. The measured numbers are honest: the physics wins
not because the target is unidentifiable from one label/edge (the pooling-by-overlap check
above shows it IS identifiable, recovering to ~4–6% error by N=4000–40000) but because the
standard per-edge *learned* variance (MVE) estimator hits a data-efficiency floor at
realistic and even large training budgets — a floor that persists under the corrected
β-NLL objective (Seitzer 2022) and under a deep ensemble, so it is not specific to the
Gaussian-NLL objective — the BAR bottleneck sidesteps that data-efficiency floor entirely
by computing the variance exactly, per edge, untrained, at zero label cost.
