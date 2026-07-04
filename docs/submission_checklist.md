# Submission checklist — Paper 1 (JCTC)

Snapshot tag: **`paper1-jctc-submission-2026-07-04`** on `bf4e88e`. `make check` = 101 tests green;
both builds 0 true undefined refs; `docs/paper_draft.bbl` + `docs/paper_jctc.bbl` present (arXiv-ready).

**One open pre-submission item:** the tracked compiled PDFs (`docs/paper_jctc.pdf`, `docs/paper_draft.pdf`)
predate the CAL-1 text fix (0.14→0.15). Regenerate before uploading: `make jctc && make paper`.

## A. JCTC / ACS submission
- [ ] Regenerate PDFs: `make jctc` (achemso, 31 pp incl. SI via `\input{paper_si}`) + `make paper` (arXiv, 20 pp).
- [ ] Main document = `paper_jctc.pdf`; SI is embedded (same PDF).
- [ ] Cover letter — one paragraph of significance: a differentiable **closed-form** BAR aleatoric variance
      (sandwich `B/I²`, O(1) backward) that (i) equals pymbar MBAR/BAR and 3000-rep MC to ~1 %, (ii) gives
      Fisher–resistance graph weights a black-box estimator cannot, and (iii) an **honest audit** of where
      the calibrated σ helps (the decision — commit/stop; the positive downstream — calibrated cycle-closure
      QC) vs where it does not (sharpness/ranking/reward-re-ranking). Frame the honest-negative scope as a
      feature, not a weakness.
- [ ] Suggested reviewers / exclusions (author's call).
- [ ] Author/affiliation confirmed: **N. L. Polomoshnov**, MSU Faculty of Bioengineering and Bioinformatics
      + Institute of Biomedical Chemistry (IBMC), Moscow; `nikitapol@fbb.msu.ru`. No competing interest.
- [ ] Data & Software Availability statement present (in `paper_jctc.tex`) — paste the Zenodo DOI once minted (§C).

## B. arXiv preprint (arXiv-first, then JCTC)
- [ ] Upload **source**, not just the PDF. Minimal bundle: `paper_draft.tex`, `paper_body.tex`,
      `paper_si.tex`, `refs.bib`, and **`paper_draft.bbl`** — arXiv runs no bibtex, so the `.bbl` MUST be
      included or the bibliography vanishes.
- [ ] **`\graphicspath{{../figs/}}` gotcha:** the figures live in `../figs` relative to `docs/`, but arXiv
      flattens every uploaded file into one directory. Either copy the referenced `figs/*.pdf` alongside the
      `.tex` and set `\graphicspath{{./}}`, or list them at the upload root. Referenced figures:
      figA_target_the_sandwich, figA_replicate_validation, figB_ood_decomposition, figB_real_decomposition,
      figC_active_learning, figD_gauge_identifiability, figE_chirality_completeness, figE_chirality_real,
      figF_target_id, figG_calibrated_stopping, figL_calibrated_cycle_closure, figL_validation,
      figQC_structures, graphical_abstract (+ any others `\includegraphics`'d in body/SI).
- [ ] Categories: `physics.chem-ph` primary; cross-list `q-bio.BM`, `stat.ML`.
- [ ] License: CC-BY recommended (max reuse/reach).
- [ ] Build the arXiv bundle locally once and compile it in a clean dir to confirm it stands alone.

## C. Zenodo archive → DOI for Data & Software Availability
- [ ] Archive the `fluor_screening/` subtree at the submission tag (`git archive paper1-jctc-submission-2026-07-04`).
- [ ] Include: `src/`, `tests/`, `figs/*.py`, `Makefile`, `docs/*.tex` + `refs.bib`, `data/` (public-derived
      only — alchemtest / OpenFF protein-ligand-benchmark / OpenFE IndustryBenchmarks2024 / ChEMBL; note each
      license), `results_*.md`, `README`.
- [ ] Mint the DOI, paste into the Data & Software Availability statement, commit that one-line edit.
- [ ] README = `make check` + per-figure `make figX` reproduction table.

## D. Reproducibility sanity (pre-upload gate)
- [ ] `make check` green (101 tests).
- [ ] `make jctc` + `make paper` → 0 true undefined refs (grep `LaTeX Warning:.*undefined`, NOT `-c undefined`
      — mciteplus head-entry lines are benign).
- [ ] Spot-regenerate figures: `make figA figL figEreal qcstruct` reproduce their headline numbers.

## E. Companion deliverables (this session)
- `docs/anticipated_referee_responses.md` — pre-drafted responses to the likely JCTC referee asks.
- `docs/paper2_scoping_brief.md` — falsifiable scoping brief for the target-finding sequel (Paper 2).
