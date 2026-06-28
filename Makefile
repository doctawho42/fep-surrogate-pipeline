# fep-surrogate-pipeline — task runner.
# Uses the project conda env `fluor_screening` by default; override with `PY=...`.
PY ?= /Users/nikitapolomosnov/anaconda3/envs/fluor_screening/bin/python

.PHONY: help test lint type check verify figA all

help:
	@echo "make test    - run the pytest suite (theorem invariants)"
	@echo "make lint    - ruff lint"
	@echo "make type    - mypy type check"
	@echo "make check   - lint + type + test (the CI gate)"
	@echo "make verify  - run the standalone Phase-0 proof verification"
	@echo "make figA    - regenerate Fig A end-to-end (one command)"

test:
	$(PY) -m pytest

lint:
	$(PY) -m ruff check src tests

type:
	$(PY) -m mypy src

check: lint type test

verify:
	$(PY) tests/verify_proofs.py

# Fig A is the Phase-1 gate; target wired once src/bar + figs/make_figA.py land.
figA:
	$(PY) figs/make_figA.py

all: check figA
