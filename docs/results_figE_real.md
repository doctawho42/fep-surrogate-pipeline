# Results — Fig E (real-molecule illustration): Theorem 4 on thalidomide

**Figure:** `figs/figE_chirality_real.{pdf,png}` · **Reproduce:** `make figEreal`
(`PYTHONPATH=src python figs/make_figE_real.py`). Deterministic (RDKit embed seed `0xC0FFEE`).
SI illustration of Theorem 4 (chirality completeness) on a real chiral drug, complementing the
synthetic designed-set proof in Fig E.

## The point
Fig E proves the chirality-completeness statement on designed 4-point tetrahedra. This adds the
same demonstration on a **real molecule** — **thalidomide**, the textbook case where the two
enantiomers differ in biological activity, so a per-enantiomer binding surrogate genuinely needs
the chirality bit.

## Construction (rigorous, no free parameters)
1. Embed one enantiomer of thalidomide in 3D (RDKit `EmbedMolecule` + MMFF), giving 19 heavy-atom
   coordinates.
2. Build the exact mirror geometry by reflecting the coordinates (`x → −x`). A reflection is an
   isometry, so this is the true enantiomer, not a re-embedding.
3. Read out both clouds with `src/bar/chiral.py`.

## Result (exact, machine precision)

| readout | M (R) vs mirror M′ (S) |
|---|---|
| all **171** pairwise distances `‖pᵢ−pⱼ‖` (O(3)-invariant) | **max\|Δ\| = 0.0e+00** |
| all pairwise dot products (Gram `VᵀV`, O(3)-invariant) | **max\|Δ\| = 0.0e+00** |
| parity-odd `0o` pseudoscalar `χ = det[p₁−p₀, p₂−p₀, p₃−p₀]` | **+1.037 → −1.037** (sum = 0, exact sign flip) |
| RDKit CIP label (independent check) | **R → S** |

So a norm-and-dot-product ("even", O(3)) readout is **provably blind** to the enantiomer pair —
every invariant it can compute is identical — while the single `0o` triple-product channel
separates them (and its sign agrees with RDKit's independent CIP assignment). This is the
collapse-and-restore of the chirality bit (Theorem 4) on real chemistry.

## Honest scope
- **Illustrative, one molecule.** The load-bearing evidence is the synthetic Fig E (collapse
  `max|Δ| ~ 1e-16` over 400 random pairs; 11.7× lower enantiomer ΔΔG MAE with the `0o` channel).
  This figure makes the same fact concrete and checkable on a named drug; it is not a new claim.
- **Coordinates as readout points.** We use the heavy-atom positions as the concrete geometric
  points `pᵢ`; an equivariant trunk's per-atom vector features transform identically under
  reflection, so the argument carries over verbatim.
- The `0o` value (`χ ≈ 1.04`) is in arbitrary coordinate units; only its **sign** (and the exact
  cancellation of the even readouts) is meaningful.
