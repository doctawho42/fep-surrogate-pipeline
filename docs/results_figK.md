# Fig K — calibrated reward vs generation (Increment 3a): KILL (honest negative)

Spec: `docs/superpowers/specs/2026-06-30-gflownet-generator-design.md`. A minimal tabular
trajectory-balance GFlowNet generates ligands over an EVALUABLE trie (target `p38`'s
measured congeneric ligands; trunk trained target-disjoint). Three reward variants drive
three GFlowNets; metric = true-hit-rate (fraction of sampled terminals in the top-quartile
of measured affinity). The GFlowNet samples ∝ reward (toy convergence test verified); the
gate ran at full convergence (steps=1500, ~4 min via the deterministic-state rollout
speedup). `make figK`.

## PRIMARY gate — calibrated generation: KILL
True-hit-rate (5 seeds): calibrated-LCB 0.323, raw-μ̂ 0.368,
MVE-σ 0.278. cal−raw -0.045, CI [-0.087, -0.014]; cal−mve
+0.045, CI [+0.001, +0.097] -> **KILL** (calibration wins iff BOTH CI
lower bounds > 0).

## Honest negative — the GLS ceiling (third manifestation)
On p38 the calibrated reward does NOT improve generation: cal ≈ raw (point estimate
even slightly below). The cause is the same Gauss-Markov / GLS ceiling behind Fig C's
ranking null and Fig I's selection KILL — p38's trunk σ_total is near-uniform
(≈0.83–1.23) and the candidate value (−ΔΔG) spread dominates it, so the LCB (value − κσ)
preserves the raw ordering. Calibration changes the generation DISTRIBUTION only when σ
varies enough to re-rank candidates (the OOD regime of Fig J); a σ-uniform target is null.
The MVE-σ foil trends worst (overconfidence hurts), but the cal-vs-MVE margin is not
significant at convergence either.

## What this means for the contour
The contour's value is the calibrated DECISION, proven twice: **Fig I** (the
commit-to-synthesis decision is trustworthy) and **Fig J** (that trust survives amortization
to unseen molecules). The GENERATION-distribution step (Fig K) is an honest negative on
p38 — consistent with the project's recurring finding that calibration is second-order
for ranking/selection/generation, first-order for decisions.

## Honest scope
Evaluable terminal set (real measured ligands, experimental ΔΔG). The toy convergence test
confirms ∝-reward sampling, so the null is NOT an under-training artifact (the gate runs at
steps=1500). A σ-dominated target (not present here) is where calibration could help
generation — untested.

## Verdict
Primary **KILL** (honest negative). The contour rests on Fig I (commit) +
Fig J (amortization); calibration does not improve the generation distribution in the
σ-uniform regime.
