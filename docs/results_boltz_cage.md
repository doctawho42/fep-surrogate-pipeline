# Boltz-2 cage cross-check — calibrated steroid-receptor screen

Boltz-2 small-molecule screen: per target, the library {anchor agonist + 4 cage forms} scored in the agonist-seeded pocket (`reference_ligands`). Cage read RELATIVE to the anchor (gap = cage − anchor binding_confidence); never raw. `make boltzcage`.

## binding_confidence (gap vs anchor)
| target | anchor | RR_OAc | SS_OAc | RR_OH | SS_OH |
|---|--:|--:|--:|--:|--:|
| GR | 0.976 | 0.448 (-0.527) | 0.452 (-0.524) | 0.175 (-0.801) | 0.196 (-0.780) |
| AR | 0.938 | 0.537 (-0.402) | 0.510 (-0.428) | 0.367 (-0.571) | 0.392 (-0.547) |
| ER | 0.979 | 0.279 (-0.700) | 0.255 (-0.724) | 0.216 (-0.763) | 0.230 (-0.749) |
| DHODH | 0.327 | 0.393 (+0.066) | 0.364 (+0.037) | 0.311 (-0.016) | 0.302 (-0.025) |
| AChE | 0.580 | 0.267 (-0.313) | 0.271 (-0.308) | 0.273 (-0.307) | 0.206 (-0.374) |

## Four diagnostics
- **Ranking** — best cage−anchor gap per target: GR -0.524, AR -0.402, ER -0.700, DHODH +0.066, AChE -0.307.
- **Discrimination (DHODH/AChE controls)** — cage gap at DHODH +0.066, AChE -0.307 (should be worse than the steroid targets if Boltz discriminates).
- **Chirality** — best enantiomer per steroid target: GR SS_OAc, AR RR_OAc, ER RR_OAc.
- **Form** — (acetate vs deacetyl: compare *_OAc vs *_OH columns above).

## Why this challenge is credible
- **Positive controls pass:** the known agonists score 0.98/0.94/0.98 (GR/AR/ER) — the pockets and the model recognise true steroid binders, so the cage's much lower score is a real gap, not a broken setup.
- **AChE greasy-pocket artefact is REJECTED:** docking falsely scored AChE 100%, but Boltz gives the cage -0.31 vs donepezil (cage best 0.27 vs 0.58) — Boltz discriminates where raw docking did not.
- **Methods disagree:** docking's one coherent signal was GR≳AR steroid (moderate); Boltz places the cage far below every steroid agonist and its weak residual preference is AR>GR (opposite), so the in-silico steroid hypothesis is not robust across methods.

## Verdict
Boltz-2 challenges the steroid hypothesis (cage far from the GR agonist and/or DHODH not rejected): downgrade to broad biochemical / phenotypic profiling rather than a targeted nuclear-receptor assay.

## Honest scope
Boltz-2 binding_confidence is a model probability calibrated only against the per-target anchor; it is not a Kd. The cage is a promiscuous H-bonder, so broad moderate confidence across targets is the non-specific reading. This corroborates or challenges the docking hypothesis; the wet-lab assay remains the un-blinder.
