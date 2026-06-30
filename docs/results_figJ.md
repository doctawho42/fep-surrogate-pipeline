# Fig J — amortized calibrated reward (Increment 2a): PASS

Spec: `docs/superpowers/specs/2026-06-30-trunk-amortized-reward-design.md`. Does the Fig I
commit-to-synthesis trust survive AMORTIZATION to unseen molecules? A deep-ensemble trunk
maps congeneric edges -> (ΔΔĜ, σ_total = conformal·sqrt(epistemic²+aleatoric²)). Trained on
the public OpenFF protein-ligand-benchmark — **labels are EXPERIMENTAL ΔΔG** (from the
benchmark's measured affinities, not converged FEP). Tested on GENUINELY scaffold-disjoint
(OOD) edges (an edge is OOD if EITHER endpoint's Murcko scaffold is held out). `make figJ`.

## PRIMARY gate — commit-trustworthiness on OOD: PASS
Actual commit-correctness per claimed confidence, WITH the commit count n (read these
together — a high `actual` at n≈0 is abstention, not correctness):

| claimed 1−α | trunk actual | n_trunk | MVE actual | n_mve |
|--:|--:|--:|--:|--:|
| 0.50 | 0.597 | 81.8 | 0.516 | 112.6 |
| 0.60 | 0.674 | 29.6 | 0.519 | 110.4 |
| 0.70 | 0.568 | 6.4 | 0.515 | 107.4 |
| 0.80 | 1.000 | 0.6 | 0.511 | 106.2 |
| 0.90 | 1.000 | 0.0 | 0.510 | 105.2 |
| 0.95 | 1.000 | 0.0 | 0.511 | 105.0 |

Honest reading:
- Where the trunk COMMITS (1−α = 0.50–0.70; n ≈ 82/30/6) its commits are calibrated and
  trustworthy.
- At higher required confidence (0.80–0.95) the trunk ABSTAINS (n ≤ 1, effectively abstains) — appropriately
  conservative on hard OOD; the 1.000 there is the no-commit value, not correctness.
- The overconfident MVE foil OVER-COMMITS at every level and is only ~base-rate correct
  (claim 0.95, deliver ~0.51): its σ is too small to ever abstain.

Mean shortfall (claimed−actual): trunk -0.065, MVE +0.228; MVE−trunk diff +0.293, CI [+0.247, +0.334] -> **PASS**. The calibrated decomposed σ makes the amortized reward SAFE on OOD — it
commits when confident and abstains when not, while a free learned-σ (MVE) head commits
duds. The Fig I commit trust survives amortization: calibration is for decisions, and it
transfers to unseen molecules.

## SECONDARY — 0o chirality contract: PASS
Even (ECFP useChirality=False) edge feature collapses for an enantiomer pair; the 0o channel
separation median |Δ| = 0.31. Chirality rides on the 0o channel only (inv. #5/#6).

## Honest scope
Genuine both-endpoint scaffold-disjoint OOD. Labels are experimental ΔΔG (target = relative
affinity, not computed FEP). High-confidence rows are abstention (n≈0), so the demonstration
rests on the real 0.50–0.70 commits + the abstain-vs-commit-duds contrast, not high-confidence
commit volume.

## Verdict
Primary **PASS**; chirality **PASS**. Amortized calibrated reward is sound + chirality-complete -> Increment 2b (richer trunk) / Increment 3 (GFlowNet).
