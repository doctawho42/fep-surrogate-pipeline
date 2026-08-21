# fep-surrogate-pipeline — task runner.
# Runs on the interpreter named by PY. Point it at the environment built from
# pyproject.toml, e.g. `make PY=/path/to/env/bin/python figA`, or export PY once.
PY ?= python

.PHONY: help test lint type check verify figA figAseeds figE figC figD figB figF figG figH figI figJ figK nioch reversedock figM figN figOOS figInf figP4 figP4b figP8 figSelf graphical gabs all

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

figCut:  # peer-review P1+P2: fixed-cutoff head-to-head + chi^2 reconciliation
	PYTHONPATH=src $(PY) figs/make_figCut.py

figStab:  # replicate stability of the QC flag set (rep 0 vs 1 vs 2 vs pooled); honest negative
	PYTHONPATH=src $(PY) figs/make_figStab.py

figOOS:  # out-of-sample QC checks: predicted-vs-observed closure chi2 + held-out localization
	PYTHONPATH=src $(PY) figs/make_figOOS.py

figInf:  # referee round 2: inference for six claims stated without it (cluster CIs, permutation, nulls)
	PYTHONPATH=src $(PY) figs/make_figInf.py

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

paper:  # generic / arXiv build (article class); shared body in docs/paper_body.tex
	cd docs && latexmk -pdf -interaction=nonstopmode paper_draft.tex

jctc:  # JCTC (ACS) submission build (achemso class); same shared body
	# achemso/mciteplus emits a benign head-entry PackageError under nonstopmode, so
	# pdflatex returns 1 even when the PDF is complete and every citation resolves. -f
	# forces latexmk through all passes (bibtex + reruns) to build the resolved PDF from a
	# cold tree in one call; the leading `-` lets make ignore the benign nonzero exit, and
	# the `test -f` guard still fails the build if no PDF was produced. (`make paper`, the
	# arXiv/article build, has no mciteplus and stays a single clean pass.)
	-cd docs && latexmk -pdf -f -interaction=nonstopmode paper_jctc.tex
	cd docs && test -f paper_jctc.pdf

all: check figA figE figC figD figB figF figG  # the gate plus the figures it covers

nioch:  # NIOCH cage screen report (operational deliverable, not a paper claim)
	PYTHONPATH=src $(PY) figs/make_nioch.py

reversedock:  # broad calibrated reverse docking -> cage hypothesis shortlist (orphan)
	PYTHONPATH=src $(PY) figs/make_reversedock.py

figM:  # Paper-2 Increment 1: orphan-benchmark validity gate (go/no-go; no scoring)
	PYTHONPATH=src $(PY) figs/make_figM.py

figN:  # Paper-2 Increment 2 Step-0: collapse-stratum validity gate (go/no-go; no scoring)
	PYTHONPATH=src $(PY) figs/make_figN.py
