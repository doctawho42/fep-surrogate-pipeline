# fep-surrogate-pipeline — task runner.
# Uses the project conda env `fluor_screening` by default; override with `PY=...`.
PY ?= /Users/nikitapolomosnov/anaconda3/envs/fluor_screening/bin/python

.PHONY: help test lint type check verify figA figE figC figD figB figF figG nioch reversedock all

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
	$(PY) -m mypy src/bar

check: lint type test

verify:
	$(PY) tests/verify_proofs.py

figA:
	$(PY) figs/make_figA.py

figE:
	$(PY) figs/make_figE.py

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

all: check figA figE figC figD figB figF figG

nioch:  # NIOCH cage screen report (operational deliverable, not a paper claim)
	PYTHONPATH=src $(PY) figs/make_nioch.py

reversedock:  # broad calibrated reverse docking -> cage hypothesis shortlist (orphan)
	PYTHONPATH=src $(PY) figs/make_reversedock.py
