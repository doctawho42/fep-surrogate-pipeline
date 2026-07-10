# fep-surrogate-pipeline — task runner.
# Uses the project conda env `fluor_screening` by default; override with `PY=...`.
PY ?= /Users/nikitapolomosnov/anaconda3/envs/fluor_screening/bin/python

.PHONY: help test lint type check verify figA figE figC figD figB figF figG figH figI figJ figK nioch reversedock figM figN all

help:
	@echo "make test    - run the pytest suite (theorem invariants)"
	@echo "make lint    - ruff lint"
	@echo "make type    - mypy type check"
	@echo "make check   - lint + type + test (the CI gate)"
	@echo "make verify  - run the standalone Phase-0 proof verification"
	@echo "make figA    - regenerate Fig A (sandwich calibration) end-to-end"
	@echo "make figE    - regenerate Fig E (chirality completeness) end-to-end"

test:
	$(PY) -m pytest

lint:
	$(PY) -m ruff check src tests

type:
	$(PY) -m mypy src/bar src/gen

check: lint type test

verify:
	$(PY) tests/verify_proofs.py

figA:
	$(PY) figs/make_figA.py

figE:
	$(PY) figs/make_figE.py

figEreal:  # SI illustration: Theorem 4 (chirality completeness) on a real drug (thalidomide)
	$(PY) figs/make_figE_real.py

qcstruct:  # SI: structural context of the Fig L QC flags (BACE1 catalytic dyad; BRD4 buried waters)
	$(PY) figs/analyze_qc_structures.py

figC:
	$(PY) figs/make_figC.py

figD:
	$(PY) figs/make_figD.py

figB:
	$(PY) figs/make_figB.py

figF:
	$(PY) figs/make_figF.py

figG:
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

boltzinputs:  # build the Boltz cage cross-check inputs (manifest + 5 screen payloads)
	PYTHONPATH=scripts $(PY) scripts/boltz_cage_inputs.py

boltzcage:  # analyze downloaded Boltz cage screen runs -> calibrated results doc + figure
	PYTHONPATH=figs $(PY) figs/make_boltz_cage.py

graphical:  # graphical abstract (vector schematic, every label exact)
	$(PY) figs/make_graphical_abstract.py

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

all: check figA figE figC figD figB figF figG

nioch:  # NIOCH cage screen report (operational deliverable, not a paper claim)
	PYTHONPATH=src $(PY) figs/make_nioch.py

reversedock:  # broad calibrated reverse docking -> cage hypothesis shortlist (orphan)
	PYTHONPATH=src $(PY) figs/make_reversedock.py

figM:  # Paper-2 Increment 1: orphan-benchmark validity gate (go/no-go; no scoring)
	PYTHONPATH=src $(PY) figs/make_figM.py

figN:  # Paper-2 Increment 2 Step-0: collapse-stratum validity gate (go/no-go; no scoring)
	PYTHONPATH=src $(PY) figs/make_figN.py
