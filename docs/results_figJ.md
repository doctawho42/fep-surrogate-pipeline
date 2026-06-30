# Fig J — amortized calibrated commit on OOD: honest re-audit

Trunk maps congeneric edges -> (ΔΔĜ, σ_total = conformal·√(epistemic²+aleatoric²)) on the
public OpenFF protein-ligand benchmark (EXPERIMENTAL ΔΔG labels); tested on both-endpoint
scaffold-disjoint (OOD) edges, 8 seeds. `make figJ`. Two honest questions the
earlier draft skipped: (A) decision quality at MATCHED commit volume (no abstention escape),
(B) calibration vs proper RECALIBRATION baselines (not just a raw overconfident MVE).

## PANEL A — matched commit volume (the abstention-proof gate)
Rank candidates by the standardized margin (μ̂−τ)/σ and commit the top-n; compare commit
precision at the SAME n. Precision at commit-fraction 25% (random base rate 0.50):
trunk-σ 0.674 · raw-μ̂ 0.680 · MVE-σ 0.504.
trunk−raw -0.006, CI [-0.016, +0.007] -> **TIE: σ-ranking does NOT beat raw-μ̂ ranking at matched volume**.
This is the honest Gauss-Markov-style ceiling: at equal commit volume the *choice* of σ does
not reliably change WHICH molecules clear the bar — calibration is not a better-ranking tool.

## PANEL B — calibration, σ isolated on the SAME (trunk) mean
All four use the trunk's ensemble mean; only σ differs, so this measures σ calibration
alone. Actual commit-correctness per claimed 1−α (0.50, 0.60, 0.70, 0.80, 0.90, 0.95):
- decomposed conformal σ (ours): 0.61 · 0.69 · 0.67 · 1.00 · 1.00 · 1.00
- epistemic-only σ (overconfident): 0.61 · 0.62 · 0.64 · 0.66 · 0.68 · 0.67
- + temperature recalibration: 0.61 · 0.68 · 0.73 · 0.85 · 0.90 · 1.00
- + split-conformal recalibration: 0.61 · 0.71 · 0.86 · 0.96 · 1.00 · 1.00

Mean shortfall (claimed−actual): ours -0.085, overconfident +0.097, temperature -0.052, conformal -0.114.
overconfident−ours +0.182 CI[+0.132,+0.225] (the raw epistemic-only σ over-claims); temperature−ours +0.034 CI[-0.018,+0.081]; conformal−ours -0.029 CI[-0.077,+0.007].

Honest reading: the overconfident epistemic-only σ over-claims, and standard recalibration (temperature / split-conformal) of the SAME mean RECOVERS the trustworthiness — so the decomposed physics σ is NOT uniquely calibrated; its advantage is that it needs no calibration set for the aleatoric term (Fig A) and is differentiable into the acquisition graph.

## SECONDARY — 0o chirality contract: PASS
Even (ECFP useChirality=False) edge feature collapses for an enantiomer pair; 0o separation
median |Δ| = 0.31 (chirality rides on the 0o channel only).

## Honest scope
Genuine both-endpoint scaffold-disjoint OOD; experimental ΔΔG labels (not computed FEP). The
matched-volume panel removes the earlier draft's abstention artifact (high-confidence cells
were n≈0). Net: calibration buys TRUST (knowing when a commit is safe), not a better commit
RANKING — consistent with the contour's decision-not-ranking finding.

## Verdict
calibration trustworthy; matched-volume ranking TIE.
