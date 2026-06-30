# Fig B (REAL) — decomposed conformal on real FEP residuals: honest audit

Replaces the synthetic 1-D Fig B headline with the REAL EnsembleTrunk on the OpenFF
protein-ligand benchmark (scaffold-disjoint OOD edges; experimental ΔΔG). At equal
~0.90 MARGINAL coverage we compare normalized conformal on the decomposed σ (**ours**)
vs flat **split**-conformal vs a fair **Mondrian**-by-target conformal. `make figBreal`.

## Marginal coverage (target 0.90)
ours 0.908 · split 0.895 · Mondrian 0.896 (8 seeds).

## Sharpness (mean half-width, kcal/mol; lower = tighter)
ours 2.16 · Mondrian 2.04 · split 1.99. **split/ours = 0.92×.**
On REAL residuals ours is **NOT sharper (wider)** than flat split — the synthetic figure's '2.2× sharper' does NOT transfer, because the
epistemic se is small relative to the large irreducible OOD error (adapt corr(half-width,
|residual|) = +0.11).

## Conditional coverage — max |bin coverage − 0.90| (lower = better)
| axis | ours | Mondrian | split |
|---|--:|--:|--:|
| per target | 0.173 | 0.224 | 0.217 |
| per \|ΔΔG\| | 0.157 | 0.174 | 0.187 |

ours−split per-target max-error diff -0.044, CI [-0.121, +0.005] (ours better iff hi < 0).

## Honest reading
- The decomposed σ buys CONDITIONAL calibration over a flat split-conformal (per-target error 0.17 vs 0.22) WITHOUT being told
  the target partition — but a fair Mondrian conformal, which IS told the partition,
  matches it (0.22).
- It is **not** sharper on real OOD data (split/ours 0.92×); the synthetic
  '2× sharper' was an artifact of an oracle DGP whose conditioning axis built σ.
- There is no real force-field-correction (Δ) head; σ is epistemic (ensemble) + a scalar
  aleatoric floor. The 3-way synthetic decomposition is illustrative only.
