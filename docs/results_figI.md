# Fig I — calibrated reward (target-contour Increment 1): PASS

Spec: `docs/superpowers/specs/2026-06-30-target-contour-design.md` (see the PIVOT
section). Reward built on `src/bar` (no trunk, no generator). κ / z from the claimed
confidence. Regenerate with `make figI`.

## PRIMARY gate (pivoted) — calibrated commit-to-synthesis: PASS
Selection failed (diagnostic below); calibration's value is the DECISION. Commit a
candidate to synthesis iff its lower-confidence bound `risk_adjusted_reward(μ̂, σ, κ=z_(1−α)) ≥ τ` (τ = top-quartile true reward). Actual
fraction of commits with true μ ≥ τ, per claimed confidence:

| claimed 1−α | actual (calibrated σ) | actual (overconfident σ) |
|--:|--:|--:|
| 0.50 | 0.874 | 0.874 |
| 0.60 | 0.892 | 0.877 |
| 0.70 | 0.916 | 0.881 |
| 0.80 | 0.941 | 0.883 |
| 0.90 | 0.974 | 0.887 |
| 0.95 | 0.982 | 0.889 |

Calibrated σ: actual **tracks/exceeds** claimed (rises 0.874→0.982,
trustworthy). Overconfident σ: **plateaus ~0.88** regardless of claimed
level (claim 0.95, deliver 0.889 — over-claims). Mean shortfall
(claimed−actual): calibrated -0.188, overconfident -0.140;
overconfident−calibrated diff +0.048, bootstrap CI [+0.037, +0.060] → **PASS**
(trustworthy iff CI lower bound > 0). This is Fig G's calibrated-decision
logic ported to the generative commit: the same LCB that does NOT improve selection DOES
make the commit trustworthy — calibration is for decisions, not rankings.

## Diagnostic (honest negative that motivated the pivot) — selection is NULL
Calibrated reward does **not** significantly improve candidate SELECTION (the GLS /
Gauss–Markov ceiling, cf. Fig C). Regime sweep precision@K, hardest regime (0.5×σ):
cal-shrink − raw +0.042, CI [-0.021, +0.094]
(not significant; shrinkage even significantly *hurts* at the 4×σ regime). Real FEP
(38 edges, BACE1 + benzene): cal−raw regret +0.000, CI [+0.000, +0.000] (broad ΔΔG spread ⇒ selection trivial). Selection is the wrong gate.

## SECONDARY gate — 0o chirality contract: PASS
Even-only reward enantiomer collapse max|Δ| = 0.0e+00 (must be ~0); 0o
reward separation median |Δ| = 4.05. The reward representation carries the
chirality bit (necessary for per-enantiomer generation, invariant #6).

## Verdict
Primary commit-calibration **PASS**; chirality **PASS**. Proceed to Increment 2 (amortised reward / trunk).
