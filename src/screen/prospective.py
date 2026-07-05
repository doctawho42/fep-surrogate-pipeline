"""Prospective-validation decision scaffold for the cage (Paper-2 prospective loop).

Thin, assay-agnostic decision layer over the summary statistics a partner lab returns
(effect + MEASURED replicate/Hill-fit sigma per (species, arm)). Makes NO target claim and NO
per-target probability (K1 + TERMINAL C forbid it): it only donates the Fig-I/G/L decision
CALCULUS applied to MEASURED assay sigma. Every sigma here is a measured sigma; no primitive
accepts an opaque pre-computed confidence.

See docs/superpowers/specs/2026-07-05-cage-prospective-loop-design.md and
docs/cage_prospective_protocol.md.
"""
from __future__ import annotations

import hashlib
import math
from dataclasses import dataclass
from pathlib import Path

import yaml


def decision_lcb(effect: float, sigma: float, z: float, tau: float) -> dict:
    """Fig-I risk-adjusted commit: LCB = effect - z*sigma; commit iff LCB >= tau.
    `sigma` is a MEASURED replicate/Hill-fit sigma."""
    lcb = effect - z * sigma
    return {"lcb": lcb, "commit": lcb >= tau}


def enantiopreference(rr: float, ss: float, sigma_rr: float, sigma_ss: float, z: float) -> dict:
    """Co-primary chirality statistic (invariants #5-6). delta = rr - ss; independent arms ->
    sigma_delta = sqrt(sigma_rr**2 + sigma_ss**2); CI = delta +/- z*sigma_delta. `discordant`
    (specificity-consistent) iff 0 is OUTSIDE the CI; 0-covering is aggregation-consistent."""
    delta = rr - ss
    sigma_delta = math.sqrt(sigma_rr ** 2 + sigma_ss ** 2)
    lo, hi = delta - z * sigma_delta, delta + z * sigma_delta
    return {"delta": delta, "ci": (lo, hi), "discordant": not (lo <= 0.0 <= hi)}


def aggregation_guard(signal: float, detergent_signal: float, sigma_detergent: float,
                      frac: float, z: float) -> dict:
    """Fig-L systematic-vs-sampling detector (sigma-aware). A signal SURVIVES detergent (is NOT a
    colloidal-aggregation artifact) iff its detergent-arm LCB stays above `frac` of the
    no-detergent signal: (detergent_signal - z*sigma_detergent) >= frac*signal."""
    survives = (detergent_signal - z * sigma_detergent) >= frac * signal
    return {"survives": survives, "artifact": not survives}


def _normal_cdf(x: float) -> float:
    return 0.5 * (1.0 + math.erf(x / math.sqrt(2.0)))


def stop_rule(effect: float, sigma_assay: float, tau: float, bound: float) -> dict:
    """Fig-G calibrated stopping: conf = Phi((effect - tau)/sigma_assay); stop iff conf >= bound.

    Takes MEASURED-sigma only; no opaque pre-computed confidence argument."""
    if sigma_assay <= 0:
        raise ValueError("sigma_assay must be > 0 (measured replicate/Hill-fit sigma)")
    conf = _normal_cdf((effect - tau) / sigma_assay)
    return {"conf": conf, "stop": conf >= bound}


PREREG_PATH = "data/cage/prospective_prereg.yaml"
# frozen anchor; must equal tests/test_prospective.py::PINNED_SHA256
PREREG_SHA256 = "66ef0bc848ab0b515b2e5b32be44230b0f528504b0b018caccc6884a0c73c488"


@dataclass(frozen=True)
class Prereg:
    """Parsed, read-only view of data/cage/prospective_prereg.yaml."""
    version: int
    species: list
    forecast: dict
    decision: dict
    disposition_table: list
    sha256: str


def sha256_of(path: str) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_prereg(path: str = PREREG_PATH, expected_sha256: str = PREREG_SHA256) -> Prereg:
    """Parse the frozen pre-registration read-only, asserting its SHA-256 equals an EXTERNAL
    anchor (`expected_sha256`, pinned once in the test suite, never edited). Comparing against an
    anchor the file-editor cannot re-derive in the same commit is what defeats coordinated
    post-hoc tuning; a self-recomputed hash would not. Raises ValueError on hash mismatch."""
    actual = sha256_of(path)
    if actual != expected_sha256:
        raise ValueError(
            f"prereg immutability violation: {path} sha256={actual} != anchor {expected_sha256}")
    d = yaml.safe_load(Path(path).read_text())
    return Prereg(version=d["version"], species=d["species"], forecast=d["forecast"],
                  decision=d["decision"], disposition_table=d["disposition_table"], sha256=actual)
