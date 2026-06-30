# NIOCH cage screen — ranked target hypotheses (PRELIMINARY · pending assays)

**This is an operational report, NOT a validated finding.** The retrospective
gate (`docs/results_figH.md`) showed structure-based target-ID does not beat a
ligand baseline and raw docking just picks the greasiest pocket; so this screen is
**heavily caveated** and exists to propose experiments, not to name a target.

## Molecule
Difluoronaphthalenone + a 1,3-dicarbonyl Michael donor (here barbituric acid),
O-acetylated; rigid fused scaffold, MW 350, two ring-fusion-coupled stereocentres
→ one diastereomer + its mirror. **Orphan scaffold** (0 ChEMBL similarity ≥40%).
Screened **per enantiomer** (R,R given / S,S mirror); racemate apparent affinity ≈
the stronger binder.

## Method (and why it is calibrated)
smina docking of each enantiomer into each pocket of an 11-target diverse panel.
We do **not** rank by raw kcal/mol (cross-pocket-incomparable). We rank each target
by the **cage's percentile among that target's known binders** docked into the same
pocket — "the cage docks better than X% of known inhibitors of this target" —
which removes the pocket-size artefact. References pooled from the gate + campaign
docking caches.

## Ranked hypotheses

| target | n_ref | R,R (given) kcal (pctl) | S,S (mirror) kcal (pctl) | best pctl |
|---|---:|---:|---:|---:|
| Acetylcholinesterase | 19 | -11.6 (100%) | -10.8 (89%) | **100%** |
| Carbonic anhydrase II | 18 | -9.4 (100%) | -8.7 (94%) | **100%** |
| Glucocorticoid receptor | 9 | -8.5 (89%) | -8.6 (89%) | **89%** |
| HIV-1 protease | 5 | -7.2 (80%) | -7.7 (80%) | **80%** |
| CDK2 (kinase) | 19 | -9.6 (74%) | -9.5 (63%) | **74%** |
| PDE5A | 12 | -9.4 (17%) | -9.9 (67%) | **67%** |
| EGFR (kinase) | 13 | -8.3 (54%) | -7.6 (31%) | **54%** |
| Estrogen receptor α | 18 | -8.7 (0%) | -8.3 (0%) | **0%** |
| Factor Xa | 8 | -8.4 (0%) | -8.0 (0%) | **0%** |
| Thrombin | 3 | -0.0 (0%) | 0.0 (0%) | **0%** |

*Calibration pending (insufficient reference docks): PDE10A (raw -8.6/-8.3).*

## How to read this (caveats — do not skip)
- **Greasy-pocket artefact persists.** Large hydrophobic pockets (e.g. AChE) score
  the rigid lipophilic cage well even after calibration; a high percentile there is
  weak evidence. Trust the **data-rich** rows (large n_ref) over thin ones.
- **Promiscuity.** The cyclic-ureide / uracil diamide is a strong, promiscuous
  H-bonder; broad moderate percentiles are consistent with a non-specific scaffold,
  not a single target. Calibration is the guard against a false promiscuous lead.
- **Docking is approximate**, the panel is **small (11 targets, not proteome-wide)**,
  and the true target may not be in it. Metal/heme/covalent targets were excluded
  (docking unreliable there).
- **Enantiomer differences** in the table are real signal worth testing (the cage is
  chiral; one enantiomer may bind preferentially).

## Broad reverse-dock (33-target pocketome) → nuclear-receptor hypothesis
To widen the search beyond the 11-target panel above, the cage (4 forms: acetate &
deacetyl, each R,R / S,S) was docked into a **~33-target diverse pocketome** and
calibrated the same way (percentile among each target's known binders). Full table:
`docs/reversedock_shortlist.md` (476 docks). The pattern that survives calibration:

- **🟠 7 targets at 100% are greasy-pocket artefacts — ignore them.** CDK2,
  VEGFR2, DPP-4, HIV integrase, AChE, MAO-B, CA-II all top out, but they span *unrelated*
  families (kinases / protease / integrase / hydrolases). A high percentile scattered
  across families is the rigid lipophilic cage fitting any roomy pocket, **not** a
  target. This is exactly the artefact the gate (`results_figH.md`) warned about.
- **🔵 The one coherent class signal is nuclear / steroid receptors:**
  **GR 80%, AR 75%, ER-α 60%** — three steroid receptors clustered in the
  moderate-high band, each data-rich (n_ref 8–10). Structural rationale: the rigid,
  lipophilic **difluoronaphthalenone is steroid-fragment-like**, so the cage plausibly
  occupies steroid-type pockets. This is the best *mechanistically-grounded* hypothesis
  in silico — **but it is moderate and non-selective** (AR≈GR≈ER), i.e. a steroid-pocket
  occupant, not a selective lead.
- **Robustness caveat (honest):** only **GR is consistent across the two independent
  screens** (panel 89% / reverse-dock 80%, different PDBs + reference sets). **ER-α is
  structure-dependent** (panel 0% vs reverse-dock 60%) → treat ER as the weakest of the
  three. So the steroid signal is real as a *class* but its strongest, most reproducible
  member is GR; rank the assay GR ≳ AR > ER accordingly.
- **Low fit (cage does not belong):** DHODH **0%** (consistent with the focused
  brequinar-tunnel test, `data/campaign/dhodh/`), Thrombin 0%, PDE10A 0%, HSP90 10%,
  COX-2 10%. Calibration correctly reports "worse than known binders" here.
- **Method-validation read:** calibration separated artefact-100% from genuine signal
  and reproduced the DHODH rejection independently — the screen behaves as designed.

## Boltz-2 cross-check (orthogonal structure model) → steroid hypothesis CHALLENGED
The docking hypothesis was tested with a stronger, orthogonal method: a Boltz-2 (AlphaFold3-class
co-folding + binding) small-molecule screen of {known agonist + the 4 cage forms} in each
receptor's **agonist-seeded pocket**, read RELATIVE to that agonist (full table + figure:
`docs/results_boltz_cage.md`, `figs/boltz_cage_crosscheck.{pdf,png}`; `make boltzcage`).
**Boltz does NOT corroborate the steroid hypothesis:**
- The known agonists score high (binding_confidence GR/AR/ER = **0.98 / 0.94 / 0.98**) — pockets
  and model recognise true steroid binders — but the cage's best form sits **far below** each:
  gaps **GR −0.52, AR −0.40, ER −0.70**. The cage is not a steroid-pocket binder at the agonist
  level; the docking GR-percentile signal looks like the "rigid lipophilic cage fits a roomy
  pocket" artefact the gate warned about.
- **The AChE greasy-pocket artefact is REJECTED by Boltz** (cage 0.21–0.27 vs donepezil 0.58):
  Boltz discriminates where raw docking falsely scored AChE 100% — evidence its negative read on
  the steroid pockets is credible.
- **The two methods disagree** (docking: GR≳AR moderate steroid; Boltz: cage far below all
  agonists, weak residual preference **AR>GR**, opposite ranking) → the in-silico steroid signal
  is **not robust across methods**. DHODH is an inconclusive control (brequinar itself scores
  low 0.33).
- Honest caveat: Boltz is not independently validated for this exact molecule, but its positive
  controls pass and it rejects the known docking artefact, so its disagreement carries weight.

## Recommended experiments (VoI-ordered, pending assays)
1. **Broad biochemical / phenotypic profiling is now the primary un-blinder** (thermal-shift /
   SPR against-panel, or Cell Painting / chemoproteomics). After the Boltz cross-check the steroid
   hypothesis is **method-dependent** (docking moderate, Boltz negative) and in-silico is
   exhausted; the orphan, promiscuous scaffold argues for breadth before depth.
2. **The nuclear-receptor reporter panel is now a cheap *bound*, not the lead.** If run, test
   **GR** (docking's one cross-screen-reproducible hit) and **AR** (Boltz's weak residual
   preference) on **both enantiomers** and **both forms** (acetate as-given + deacetyl active; on
   GR/AR the acetate scored slightly higher in Boltz). A clean negative *bounds* the orphan claim,
   but the two in-silico methods no longer agree that it will be positive.
3. A **broad biochemical/binding panel** (e.g. a kinase/against-panel + a
   thermal-shift / SPR screen) or **phenotypic profiling (Cell Painting /
   chemoproteomics)** is the real un-blinder — the orphan scaffold + promiscuity flag
   argue for breadth before depth, and in-silico has now been exhausted.
4. For any hit, **test both enantiomers separately** to confirm a chirality preference.
5. Escalate only confirmed hits up the ladder (binding → functional → co-crystal).

*Generated by `make nioch` from the docking caches; regenerates as the campaign adds
reference docks. Not a publication claim.*
