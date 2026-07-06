"""Close-the-loop accuracy race: does acting on the calibrated QC flag improve accuracy?

Guided (remove top-|z| QC-flagged edges) vs random edge removal, MUE-vs-experiment, on flagged
systems grounded on public ChEMBL affinity. Reuses src/bar/qc.py (GLS network, repair_order) and
figs/analyze_eg5_accuracy.py (mean-aligned MUE). NO new MD.

See docs/superpowers/specs/2026-07-05-closeloop-accuracy-race-design.md.
"""
from __future__ import annotations

import hashlib
from dataclasses import dataclass
from pathlib import Path

import yaml

CLOSELOOP_PREREG_PATH = "data/openfe_replicates/closeloop_prereg.yaml"
# frozen anchor; must equal tests/test_closeloop.py::PINNED_SHA256
CLOSELOOP_PREREG_SHA256 = "35cb3b5a81c6878fddaa2f8cfc0107596fbdd9f1c6de3ab132524632c97ff86e"


@dataclass(frozen=True)
class Prereg:
    systems: list
    coverage_min: float
    n_perm: int
    assay_order: list
    target_reduced_chi2: float
    min_grounded: int
    sha256: str


def sha256_of(path: str) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_prereg(path: str = CLOSELOOP_PREREG_PATH,
                expected_sha256: str = CLOSELOOP_PREREG_SHA256) -> Prereg:
    """Parse the frozen pre-registration read-only, asserting its SHA-256 equals an EXTERNAL anchor
    (pinned in the test suite, never edited) — defeats coordinated post-hoc tuning. Raises
    ValueError on mismatch."""
    actual = sha256_of(path)
    if actual != expected_sha256:
        raise ValueError(
            f"prereg immutability violation: {path} sha256={actual} != {expected_sha256}"
        )
    d = yaml.safe_load(Path(path).read_text())
    return Prereg(systems=d["systems"], coverage_min=d["coverage_min"], n_perm=d["n_perm"],
                  assay_order=d["assay_order"], target_reduced_chi2=d["target_reduced_chi2"],
                  min_grounded=d["min_grounded"], sha256=actual)
