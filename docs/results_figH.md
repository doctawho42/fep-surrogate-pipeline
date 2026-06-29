# Results — Fig H: THE GATE (structure-based vs ligand-based target-ID) — **FAIL → two papers**

**This is the load-bearing gate** (`docs/target_finding_plan.md` §2). Honest verdict:
the structure-based (docking) target-finding arm **does not beat the ligand baseline**.
Per the pre-registered decision — and the advisor's explicit instruction *"не бьёт →
клетка в Paper 2, спайн шипишь как есть; не спасай мердж"* — the cage/target-finding
work defers to **Paper 2**; the methods+theory spine (Figs A–G + `paper_draft.tex`) ships
as **Paper 1 as-is**.

**Code:** `src/screen/{dock,recovery}.py` (+4 tests) · `figs/make_figH.py` · `make figH`.

## What was tested
A diverse-fold panel of 10 targets (EGFR, CDK2, HIV-1 protease, AChE, ERα, CA-II, FXa,
PDE5A, thrombin, GR), each a holo PDB pocket (co-crystal ligand auto-detected for the
docking box) + ChEMBL actives. Each scaffold-disjoint test ligand is docked (smina) into
**all 10 pockets**; the true target should rank top. Compared to the **ligand-shape**
baseline (max Tanimoto to a target's training actives — the Fig F method) and random.

## Result: structure loses, decisively
(≈37–50 fully-docked scaffold-disjoint test ligands)

| ranker | top-1 recovery | top-3 | recovery AUROC |
|---|---|---|---|
| **structure (raw docking)** | **0.16** | 0.31 | 0.62 |
| structure (per-pocket z-normalised) | 0.08 | — | — |
| **ligand-shape (Fig F)** | **0.95** | 1.00 | 0.99 |
| random | 0.10 | 0.30 | 0.50 |

## Why (the honest diagnosis)
1. **Raw docking scores are not comparable across pockets** — the known failure of
   "reverse docking". Docking simply prefers the largest/greasiest pocket: AChE recovers
   4/5 of *its* ligands but also attracts everyone else's, while EGFR/CDK2/CA-II/FXa get
   0/5. This is the same disease as raw `1/I` for BAR variance — an uncalibrated score.
2. **Cheap calibration doesn't rescue it.** Per-pocket z-normalisation made recovery
   *worse* (0.08): docking here carries little genuine target-discriminating signal.
3. **The retrospective ChEMBL benchmark cannot reach the orphan regime.** Scaffold-split
   test ligands on well-populated targets stay similar to training (median max-Tanimoto
   **0.64**; only 5% < 0.4), so the ligand-shape baseline trivially wins and there is no
   room to demonstrate structure's value *where ligand-shape fails* — which is precisely
   the cage's regime (0 ChEMBL similarity). Validating structure-based target-ID for
   genuinely **orphan** ligands retrospectively needs orphan ligand→target→structure
   triples that this public benchmark does not provide.

## Decision (ratified): two papers
- **Paper 1 (now):** the methods+theory spine — four theorems, the differentiable
  self-calibrating BAR bottleneck, Figs A–G (calibration, OOD decomposition, chirality,
  gauge identifiability, calibrated stopping; two honest negatives). `docs/paper_draft.tex`.
- **Paper 2 (matures with work/data):** the cage target-finding arm — but it needs the
  real structure-prediction stack (Boltz-ABFE / cofolding rather than raw docking;
  target-disjoint calibration with reference-active docking per pocket; a curated pocket
  library; and ideally orphan ligand→target validation). That is research, **not a gate
  patch** — do not rescue the merge now.
- **NIOCH operational debt stays separate and deliverable:** the cage can still be run
  through a docking/Boltz campaign and handed back as ranked hypotheses + recommended
  assays ("pending assays"), independent of where the publication lands — but **as a
  report, not a paper claim** (the uncalibrated docking shown here is exactly why a bare
  "we found target X" would be a false lead).

## Gate
`make check` green (45 tests incl. recovery metrics). Structure-based target-ID does
**not** beat the ligand baseline on this fair retrospective test → **second arm does not
stand → ship Paper 1 spine, defer the cage to Paper 2.**
