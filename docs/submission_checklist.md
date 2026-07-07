# Submission checklist — Paper 1 (JCTC)

Snapshot tag: **`paper1-jctc-submission-2026-07-08`** on `267ae0c` (supersedes the
`2026-07-04`/`bf4e88e` snapshot). `make check` = **160 tests green**; both builds regenerated
with the **D1 conservation-law theorem (Thm 4, `thm:cons`) + Fig Lev** (`make figLev`); `make jctc`
and `make paper` both **0 true undefined refs**; tracked `docs/paper_jctc.pdf` + `docs/paper_draft.pdf`
are current; `docs/paper_draft.bbl` + `docs/paper_jctc.bbl` present (arXiv-ready).

**arXiv source bundle assembled + verified standalone:** `paper1_arxiv_bundle.tar.gz` (source +
`.bbl` + 17 figs, `\graphicspath` flattened to `./`) compiles to **0 undefined** via the arXiv model
(`pdflatex` x3, no bibtex, shipped `.bbl`). Built in the session scratchpad -- copy it out before the
scratchpad is cleared.

**arXiv build gotcha (local-only):** `make paper` (article, single-pass `latexmk`) leaves NEW
`\label`/`\ref` undefined on the first run from a cold `.aux` (run `make paper` twice, or use `-f`);
`make jctc` uses `-f` and resolves in one call. arXiv's own `latexmk` runs to fixpoint.

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
      figA_replicate_validation, figA_target_the_sandwich, figB_real_decomposition, figC_active_learning,
      figD_gauge_identifiability, figE_chirality_completeness, figE_chirality_real, figF_target_id,
      figG_calibrated_stopping, figJ_amortized_reward, figK_calibrated_generation,
      figL_calibrated_cycle_closure, figL_validation, figLcausal_guided_vs_random,
      figLev_observability, figQC_structures, graphical_abstract (17 total, verified via the
      `\includegraphics` grep; all copied into `paper1_arxiv_bundle.tar.gz`).
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
- [ ] `make check` green (160 tests).
- [ ] `make jctc` + `make paper` → 0 true undefined refs (grep `LaTeX Warning:.*undefined`, NOT `-c undefined`
      — mciteplus head-entry lines are benign).
- [ ] Spot-regenerate figures: `make figA figL figEreal qcstruct` reproduce their headline numbers.

## E. Companion deliverables (this session)
- `docs/anticipated_referee_responses.md` — pre-drafted responses to the likely JCTC referee asks.
- `docs/paper2_scoping_brief.md` — falsifiable scoping brief for the target-finding sequel (Paper 2).
