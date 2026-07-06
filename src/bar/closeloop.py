"""Close-the-loop accuracy race: does acting on the calibrated QC flag improve accuracy?

Guided (remove top-|z| QC-flagged edges) vs random edge removal, MUE-vs-experiment, on flagged
systems grounded on public ChEMBL affinity. Reuses src/bar/qc.py (GLS network, repair_order) and
figs/analyze_eg5_accuracy.py (mean-aligned MUE). NO new MD.

See docs/superpowers/specs/2026-07-05-closeloop-accuracy-race-design.md.
"""
from __future__ import annotations

import csv
import hashlib
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import yaml
from scipy.stats import norm

from bar.qc import gls_network, repair_order

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


def _f(x) -> float:
    try:
        return float(x)
    except (TypeError, ValueError):
        return math.nan


def load_system_edges(system: str, exp: dict) -> list:
    """combined-csv edges for `system` (3-replicate mean; se=sqrt(mean se_k^2 / n)); keep edges
    whose both endpoints are grounded (in `exp`)."""
    p = Path("data/openfe_replicates/combined_pymbar4_edge_data.csv")
    edges = []
    for r in csv.DictReader(p.open()):
        if r["system name"] != system:
            continue
        dd, se = [], []
        for k in (0, 1, 2):
            cD = _f(r[f"complex_repeat_{k}_DG (kcal/mol)"])
            cd = _f(r[f"complex_repeat_{k}_dDG (kcal/mol)"])
            sD = _f(r[f"solvent_repeat_{k}_DG (kcal/mol)"])
            sd = _f(r[f"solvent_repeat_{k}_dDG (kcal/mol)"])
            if not any(math.isnan(v) for v in (cD, cd, sD, sd)):
                dd.append(cD - sD)
                se.append(math.sqrt(cd * cd + sd * sd))
        if len(dd) >= 2 and r["ligand_A"] in exp and r["ligand_B"] in exp:
            edges.append((r["ligand_A"], r["ligand_B"], float(np.mean(dd)),
                          math.sqrt(float(np.mean(np.array(se) ** 2)) / len(dd))))
    return edges


def _mue_on(fit, exp: dict, nodes: list) -> float:
    """Mean-aligned MUE of `fit` potentials vs exp, restricted to `nodes` (gauge re-aligned on
    set)."""
    pred = dict(zip(fit.nodes, fit.potentials, strict=True))
    keys = [k for k in nodes if k in pred and k in exp]
    if len(keys) < 2:
        return math.nan
    p = np.array([pred[k] for k in keys])
    e = np.array([exp[k] for k in keys])
    p = p - p.mean() + e.mean()
    return float(np.mean(np.abs(p - e)))


def delta_mue(edges: list, removed: list, exp: dict) -> float:
    """MUE(full) - MUE(after removing `removed`), BOTH on the node set surviving the removal."""
    keep = [e for i, e in enumerate(edges) if i not in set(removed)]
    if len(keep) < 1:
        return math.nan
    fit_after = gls_network(keep)
    surviving = [n for n in fit_after.nodes if n in exp]
    fit_full = gls_network(edges)
    mue_full = _mue_on(fit_full, exp, surviving)
    mue_after = _mue_on(fit_after, exp, surviving)
    return mue_full - mue_after


def system_effect(edges: list, exp: dict, n_perm: int, target_rchi2: float, seed: int = 0) -> dict:
    """guided (repair_order top-|z|) vs random ΔMUE at the trajectory endpoint K; per-system p."""
    removed, _ = repair_order(edges, target_reduced_chi2=target_rchi2)
    K = len(removed)
    guided = delta_mue(edges, removed, exp)
    rng = np.random.default_rng(seed)
    E = len(edges)
    rnd = np.array([delta_mue(edges, list(rng.choice(E, K, replace=False)), exp)
                    for _ in range(n_perm)])
    rnd = rnd[~np.isnan(rnd)]
    p = float(np.mean(rnd >= guided)) if rnd.size else math.nan
    below_5pct = bool(guided > np.percentile(rnd, 95)) if rnd.size else False
    random_mean = float(np.mean(rnd)) if rnd.size else math.nan
    return {"k": K, "guided": float(guided), "random_mean": random_mean,
            "p": p, "below_5pct": below_5pct,
            "curve": {"guided_removed": K, "n_edges": E}}


def combine(effects: list) -> dict:
    """One-sided Stouffer on per-system p + a sign test on the per-system effect (guided-random).
    Systems with a NaN p (degenerate networks where every permutation draw was NaN) are dropped
    from ALL three sub-statistics consistently, so Stouffer, the sign test, and the majority-below
    count share one effective sample size."""
    valid = [e for e in effects if not math.isnan(e["p"])]
    n = len(valid)
    if n == 0:
        return {"stouffer_p": math.nan, "sign_p": math.nan, "verdict": "NULL"}
    z = norm.isf(np.array([e["p"] for e in valid]))  # one-sided: small p -> large z
    stouffer = float(norm.sf(z.sum() / math.sqrt(n)))
    signs = sum(1 for e in valid if e["guided"] > e["random_mean"])
    sign_p = float(sum(math.comb(n, i) for i in range(signs, n + 1)) / (2 ** n))
    majority_below = sum(bool(e.get("below_5pct")) for e in valid) > n / 2
    verdict = "SUCCESS" if (stouffer < 0.05 and majority_below) else "NULL"
    return {"stouffer_p": stouffer, "sign_p": sign_p, "verdict": verdict}
