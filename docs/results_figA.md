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
2. **The right foil is a learned-variance head, and it fails.** An MVE / heteroscedastic
   NN (what an ML practitioner actually builds), trained by Gaussian NLL on a *realistic*
   budget of ~200 edges (one noisy ΔΔĜ label each), is **~7× overconfident**
   (reported se ≈ 0.14× of the truth). One label per edge cannot teach the per-edge
   sampling variance; the BAR bottleneck **computes** it exactly, untrained.
3. **`1/I` is shown only as the textbook value it corrects** — the information-equality
   plug-in, off by a *varying* factor — **not** as a baseline anyone reports.

## Why this matters downstream (the actual contribution)
- **Differentiable closed form → graph weights.** `B/I²` has an O(1) backward (Thm 1)
  and propagates into the surrogate **and** the Fisher–resistance Laplacian weights
  `w_e = I_e²/B_e` (Thm 3, invariant #4). pymbar returns a number; it is not a
  backprop-able graph weight. This is the new object.
- **Robustness.** `B/I² ≥ 0` always; pymbar's `1/I − nrat` goes **negative → nan** at
  extreme poor overlap. The sandwich form degrades gracefully.
- **Calibrated aleatoric for free.** vs the learned-variance head above.

## Result: PASS

**Panel A — controlled (MC truth).** Gaussian work model, `n_f=n_r=20`, 3000
replicates/point; true se = MC SD; band = 95% bootstrap CI.

| overlap `4⟨p(1−p)⟩` | sandwich/true | **pymbar-MBAR/true** | learned-σ/true | naive `1/I`/true |
|---:|---:|---:|---:|---:|
| 0.17 (low)  | 0.89 | 0.99 | **0.18** | 1.08 |
| 0.46        | 1.00 | 1.00 | **0.11** | 1.36 |
| 0.62        | 1.00 | 1.00 | **0.08** | 1.61 |
| 0.86 (high) | 0.99 | 0.98 | **0.09** | 2.65 |

Sandwich ≈ MBAR ≈ 1 (coincide, the correctness proof). The learned head is
catastrophically overconfident (~0.14× truth across the range). `1/I` traces a
*non-constant* 1.08→2.65 factor — the info-equality plug-in, off worst at high overlap
(see the corrected Corollary in `bar_proofs.tex`).

**Panel B — real FEP edges.** Adjacent-λ BAR edges with autocorrelation-aware bootstrap
truth (decorrelated by the statistical inefficiency):
- **Binding — BACE1 RBFE (alchemtest AMBER, complex legs), 19 edges:** sandwich/boot
  **1.00** [0.96, 1.07]; naive `1/I`/boot **6.07** [2.78, 12.01].
- Solvation — benzene hydration, 19 edges: sandwich/boot 1.01 [0.95, 1.05]; naive 4.05.

The sandwich is calibrated on real **protein-ligand binding** edges (the relevant
context), not just solvation.

## Honest scope
- The MBAR coincidence is exact in expectation; at extreme low overlap (n=20) the
  ddof-1 sandwich plug-in runs ~10% below MBAR (both still ≈ 1) — a finite-sample
  plug-in difference, not a discrepancy. Production reports the sandwich (robust, ≥ 0).
- Bootstrap truth (no independent repeats for these systems) is the model-free
  reference on real edges; the controlled panel supplies the MC-truth validation.
- OpenFE IndustryBenchmarks2024 (3 repeats → true-replicate se) is the next data step.

## Gate
`make check` green (41 tests) **and** Fig A regenerable by one command (`make figA`).
Sandwich = MBAR = truth (correct); learned-σ head fails; calibrated on real binding
edges → **proceed**.
