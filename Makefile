# fep-surrogate-pipeline — task runner.
# Runs on the interpreter named by PY. Point it at the environment built from
# pyproject.toml, e.g. `make PY=/path/to/env/bin/python figA`, or export PY once.
PY ?= python

.PHONY: figDesign figGround figNoise figRelease help test lint type check verify figA figAseeds figE figC figD figB figF figG figH figI figJ figK nioch reversedock figM figN figOOS figInf figCal figP4 figP4b figP8 figSelf graphical gabs all

help:  # list every target with its description
	@echo "Targets. Run with PY=/path/to/python if the default interpreter is not the right one."
	@echo ""
	@awk 'BEGIN{FS=":"} /^[a-zA-Z0-9_]+:/ { \
		target=$$1; desc=""; \
		if (match($$0, /#[ \t]*/)) desc=substr($$0, RSTART+RLENGTH); \
		printf "  %-12s %s\n", target, desc }' $(MAKEFILE_LIST)

test:  # the unit tests: the correctness invariants of the estimator and the instrument
	$(PY) -m pytest

lint:  # ruff over src and tests
	$(PY) -m ruff check src tests

type:  # mypy over src/bar and src/gen
	$(PY) -m mypy src/bar src/gen

check: lint type test  # lint, type check and tests; the gate

verify:  # check the theorem statements numerically, outside the test suite
	$(PY) tests/verify_proofs.py

figA:  # Fig 2: sandwich variance against MBAR, truth and the learned foils
	$(PY) figs/make_figA.py

figAseeds:  # P6c: multi-seed spread of the learned-variance foils (frozen 5 seeds)
	$(PY) figs/make_figA.py seeds

figE:  # SI: chirality completeness, the parity-odd readout ablation
	$(PY) figs/make_figE.py

figEreal:  # SI illustration: chirality completeness on a real drug (thalidomide)
	$(PY) figs/make_figE_real.py

qcstruct:  # SI: structural context of the Fig L QC flags (BACE1 catalytic dyad; BRD4 buried waters)
	$(PY) figs/analyze_qc_structures.py

figC:  # SI: gauge-aware cost-aware active learning
	$(PY) figs/make_figC.py

figD:  # SI: the knowledge gradient vanishes on gauge-redundant directions
	$(PY) figs/make_figD.py

figB:  # SI: decomposed uncertainty against conformal baselines, synthetic
	$(PY) figs/make_figB.py

figF:  # SI: ligand-similarity reverse screening on ChEMBL
	$(PY) figs/make_figF.py

figG:  # SI: calibrated stopping, and the assumed-error sweep behind it
	$(PY) figs/make_figG.py

figH:  # the target-finding gate (structure vs ligand); docking, cached to data/figH
	PYTHONPATH=src $(PY) figs/make_figH.py

figI:  # target-contour Increment 1: calibrated commit-to-synthesis gate + 0o chirality
	PYTHONPATH=src $(PY) figs/make_figI.py

figJ:  # target-contour Increment 2a: amortized reward commit-on-OOD gate
	PYTHONPATH=src $(PY) figs/make_figJ.py

figK:  # target-contour Increment 3a: calibrated-generation GFlowNet gate (honest negative)
	PYTHONPATH=src $(PY) figs/make_figK.py

figBreal:  # honest audit: Fig B on REAL OOD FEP residuals (not sharper; fair Mondrian foil)
	PYTHONPATH=src $(PY) figs/make_figB_real.py

figArep:  # independent-replicate validation of the reported se (OpenFE IndustryBenchmarks2024)
	PYTHONPATH=src $(PY) figs/make_figA_replicates.py

gsweep:  # SI: autocorrelation robustness of the sandwich calibration (raw n vs n_eff = n/g, BACE1)
	PYTHONPATH=src $(PY) figs/make_figAC_gsweep.py

figL:  # IMPACT: calibrated cycle-closure QC (separates systematic from sampling; FPR = f(calibration))
	PYTHONPATH=src $(PY) figs/make_figL.py

figLval:  # IMPACT validation: flags are causal (repair test) + reproduce out-of-sample (held-out replicates)
	PYTHONPATH=src $(PY) figs/make_figL_validation.py

figLcausal:  # close-the-loop: does acting on the QC flag improve accuracy vs experiment?
	PYTHONPATH=src $(PY) figs/make_figLcausal.py

figLev:  # D1: per-edge observability map (sum_h==dof, bridges) + pre-registered predictive falsifier
	PYTHONPATH=src $(PY) figs/make_figLev.py

figHodge:  # Theorem 5: auditability map + where the error against experiment lives + repair race
	PYTHONPATH=src $(PY) figs/make_figHodge.py

figDesign:  # the observability map as a design rule: edges to add for no bridge + a delta* target
	PYTHONPATH=src $(PY) figs/make_figDesign.py

figRelease:  # the cross-release variance split over every multi-release system, all three replicates
	PYTHONPATH=src $(PY) figs/make_figRelease.py

figNoise:  # the label-noise floor under the visible fraction: exactness, shrinkage, crossings
	PYTHONPATH=src $(PY) figs/make_figNoise.py

figGround:  # the visible fraction on the benchmark's CURATED per-edge labels, 14 systems
	PYTHONPATH=src $(PY) figs/make_figGround.py

figCut:  # peer-review P1+P2: fixed-cutoff head-to-head + chi^2 reconciliation
	PYTHONPATH=src $(PY) figs/make_figCut.py

figStab:  # replicate stability of the QC flag set (rep 0 vs 1 vs 2 vs pooled); honest negative
	PYTHONPATH=src $(PY) figs/make_figStab.py

figOOS:  # out-of-sample QC checks: predicted-vs-observed closure chi2 + held-out localization
	PYTHONPATH=src $(PY) figs/make_figOOS.py

figInf:  # referee round 2: inference for six claims stated without it (cluster CIs, permutation, nulls)
	PYTHONPATH=src $(PY) figs/make_figInf.py

figCal:  # the per-edge calibration ratio, with the four published summaries of it marked on it
	PYTHONPATH=src $(PY) figs/make_figCal.py

figP4:  # peer-review P4: QC sweep under the real heterogeneous learned-sigma profile
	PYTHONPATH=src $(PY) figs/make_figP4.py

figP4b:  # peer-review P4b: dose-response of the QC to sigma miscalibration
	PYTHONPATH=src $(PY) figs/make_figP4b.py

figP8:  # peer-review P8: per-system flag robustness to each system's own calibration error
	PYTHONPATH=src $(PY) figs/make_figP8.py

figSelf:  # exploratory TWO-SIDED self-calibration of the QC null (P8 = the pre-registered one-sided check)
	PYTHONPATH=src $(PY) figs/make_figSelf.py

boltzinputs:  # build the Boltz cage cross-check inputs (manifest + 5 screen payloads)
	PYTHONPATH=scripts $(PY) scripts/boltz_cage_inputs.py

boltzcage:  # analyze downloaded Boltz cage screen runs -> calibrated results doc + figure
	PYTHONPATH=figs $(PY) figs/make_boltz_cage.py

graphical:  # graphical abstract (Figma-style pastel redraw: BAR bottleneck -> cycle-closure QC)
	$(PY) figs/make_graphical_abstract.py

gabs: graphical  # alias for `make graphical`

paper:  # generic / arXiv build: manuscript and Supporting Information as separate files
	# Same two-file shape as `make jctc`, in the article class: docs/paper_draft.pdf and
	# docs/paper_draft_si.pdf, the second S-numbered. The pair references each other's labels
	# through xr-hyper, so neither settles until both have been built after the other; -g
	# forces the last two passes, which latexmk would skip as up to date while their external
	# numbers are still stale.
	cd docs && latexmk -pdf -interaction=nonstopmode paper_draft.tex
	cd docs && latexmk -pdf -interaction=nonstopmode paper_draft_si.tex
	cd docs && latexmk -pdf -g -interaction=nonstopmode paper_draft.tex
	cd docs && latexmk -pdf -g -interaction=nonstopmode paper_draft_si.tex

jctc:  # JCTC (ACS) submission build: manuscript and Supporting Information as separate files
	# ACS requires the Supporting Information as its own file, so this target builds two:
	# docs/paper_jctc.pdf (title to references) and docs/paper_jctc_si.pdf (S-numbered).
	# They reference each other's labels through xr, which reads the other document's .aux,
	# so neither is settled until both have been built after the other. Four alternating
	# passes reach that fixed point from a cold tree; -g forces the last two, which latexmk
	# would otherwise skip as up to date while their external numbers are still stale.
	# achemso/mciteplus emits a benign head-entry PackageError under nonstopmode, so
	# pdflatex returns 1 even when the PDF is complete and every citation resolves. -f
	# forces latexmk through all passes (bibtex + reruns); the leading `-` lets make ignore
	# the benign nonzero exit, and the `test -f` guard still fails the build if no PDF was
	# produced. (`make paper`, the arXiv/article build, has no mciteplus.)
	-cd docs && latexmk -pdf -f -interaction=nonstopmode paper_jctc.tex
	-cd docs && latexmk -pdf -f -interaction=nonstopmode paper_jctc_si.tex
	-cd docs && latexmk -pdf -f -g -interaction=nonstopmode paper_jctc.tex
	-cd docs && latexmk -pdf -f -g -interaction=nonstopmode paper_jctc_si.tex
	cd docs && test -f paper_jctc.pdf && test -f paper_jctc_si.pdf

all: check figA figE figC figD figB figF figG  # the gate plus the figures it covers

nioch:  # NIOCH cage screen report (operational deliverable, not a paper claim)
	PYTHONPATH=src $(PY) figs/make_nioch.py

reversedock:  # broad calibrated reverse docking -> cage hypothesis shortlist (orphan)
	PYTHONPATH=src $(PY) figs/make_reversedock.py

figM:  # Paper-2 Increment 1: orphan-benchmark validity gate (go/no-go; no scoring)
	PYTHONPATH=src $(PY) figs/make_figM.py

figN:  # Paper-2 Increment 2 Step-0: collapse-stratum validity gate (go/no-go; no scoring)
	PYTHONPATH=src $(PY) figs/make_figN.py
