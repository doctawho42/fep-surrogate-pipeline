"""Step-0 validity gate for the Paper-2 orphan benchmark (Increment 1).

Does the ligand-shape null collapse to near-random in the orphan stratum? No structure scoring;
the "score" is max Tanimoto to each pocket's actives (higher = better -> lower_better=False).
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from numpy.typing import NDArray

from screen.recovery import recovery_at_k, recovery_auroc
from screen.stratify import ecfp, max_tanimoto

SEED = 0
# pre-registered thresholds (spec §2)
REC1_MAX, AUROC_MAX, CI_HI_MAX = 0.10, 0.60, 0.15
MIN_ORPHAN_QUERIES, MIN_FOLD_CLUSTERS = 30, 8


def shape_score_matrix(queries: pd.DataFrame, pocket_actives: dict[str, list[str]],
                       pocket_order: list[str]) -> tuple[NDArray, NDArray]:
    """scores[q, t] = max Tanimoto of query q to pocket_order[t]'s actives (self-excluded).
    true_idx[q] = index in pocket_order of the query's true target."""
    pocket_fps = {t: [fp for fp in (ecfp(s) for s in pocket_actives.get(t, [])) if fp is not None]
                  for t in pocket_order}
    Q, N = len(queries), len(pocket_order)
    scores = np.zeros((Q, N))
    true_idx = np.zeros(Q, dtype=int)
    tpos = {t: i for i, t in enumerate(pocket_order)}
    for qi, (_, row) in enumerate(queries.iterrows()):
        qfp = ecfp(row["smiles"])
        for t in pocket_order:
            scores[qi, tpos[t]] = 0.0 if qfp is None else max_tanimoto(qfp, pocket_fps[t])
        true_idx[qi] = tpos[row["target"]]
    return scores, true_idx


def bootstrap_recovery1(scores: NDArray, true_idx: NDArray,
                        n_seed: int = 5) -> tuple[float, float, float]:
    """Point recovery@1 + 95% bootstrap CI (resample queries) across n_seed*200 resamples."""
    point = float(recovery_at_k(scores, true_idx, lower_better=False)[0])
    boots = []
    Q = len(true_idx)
    for seed in range(n_seed):
        rng = np.random.default_rng(SEED + seed)
        for _ in range(200):
            idx = rng.integers(0, Q, Q)
            boots.append(float(recovery_at_k(scores[idx], true_idx[idx], lower_better=False)[0]))
    lo, hi = np.percentile(boots, [2.5, 97.5])
    return point, float(lo), float(hi)


def _stratum_result(queries: pd.DataFrame, scores: NDArray, true_idx: NDArray,
                    stratum: str) -> dict:
    mask = (queries["stratum"] == stratum).to_numpy()
    if mask.sum() == 0:
        return {"n": 0, "recovery1": float("nan"), "auroc": float("nan"), "ci": (float("nan"),) * 2}
    s, t = scores[mask], true_idx[mask]
    pt, lo, hi = bootstrap_recovery1(s, t)
    return {"n": int(mask.sum()), "recovery1": pt,
            "auroc": recovery_auroc(s, t, lower_better=False), "ci": (lo, hi)}


def verdict(queries: pd.DataFrame, scores: NDArray, true_idx: NDArray,
            n_pockets: int, n_fold_clusters: int) -> dict:
    """Apply P1/P2/P3 -> {'verdict': 'PASS'|'VALIDITY_KILL', per-stratum results, reasons}."""
    orphan = _stratum_result(queries, scores, true_idx, "orphan")
    high = _stratum_result(queries, scores, true_idx, "high")
    reasons: list[str] = []
    p1 = orphan["n"] >= MIN_ORPHAN_QUERIES and n_fold_clusters >= MIN_FOLD_CLUSTERS
    if not p1:
        reasons.append(f"P1 fail: {orphan['n']} orphan queries (need {MIN_ORPHAN_QUERIES}), "
                       f"{n_fold_clusters} fold clusters (need {MIN_FOLD_CLUSTERS})")
    p2 = (orphan["n"] > 0 and orphan["recovery1"] <= REC1_MAX
          and orphan["auroc"] <= AUROC_MAX and orphan["ci"][1] <= CI_HI_MAX)
    if not p2:
        reasons.append(f"P2 fail: orphan shape-null recovery@1={orphan['recovery1']:.3f} "
                       f"(<= {REC1_MAX}), AUROC={orphan['auroc']:.3f} (<= {AUROC_MAX}), "
                       f"CI_hi={orphan['ci'][1]:.3f} (<= {CI_HI_MAX})")
    v = "PASS" if (p1 and p2) else "VALIDITY_KILL"
    return {"verdict": v, "orphan": orphan, "high": high, "n_pockets": n_pockets,
            "n_fold_clusters": n_fold_clusters, "reasons": reasons}
