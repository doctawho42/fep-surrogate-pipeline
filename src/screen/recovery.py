"""Reverse-screening recovery metrics (the gate's honest core).

A query is one ligand with a known true target; a ranker assigns each candidate target a
score. We measure whether the true target is recovered top-k, and the per-query AUROC of
the true target against the decoys. Structure-based (docking) is compared to the
ligand-shape (Tanimoto) baseline on the SAME queries (docs/target_finding_plan.md §4).
"""
from __future__ import annotations

import numpy as np
from numpy.typing import NDArray


def _ranks(scores: NDArray, true_idx: NDArray, lower_better: bool) -> NDArray:
    """0-based rank of the true target per query (0 = best)."""
    s = scores if lower_better else -scores
    order = np.argsort(s, axis=1)
    ranks = np.empty(len(true_idx), dtype=int)
    for q in range(len(true_idx)):
        ranks[q] = int(np.where(order[q] == true_idx[q])[0][0])
    return ranks


def recovery_at_k(scores: NDArray, true_idx: NDArray, lower_better: bool = True) -> NDArray:
    """Fraction of queries whose true target is within the top-k, for k=1..n_targets."""
    ranks = _ranks(scores, true_idx, lower_better)
    n = scores.shape[1]
    return np.array([(ranks < k).mean() for k in range(1, n + 1)])


def recovery_auroc(scores: NDArray, true_idx: NDArray, lower_better: bool = True) -> float:
    """Mean per-query AUROC of the true target vs the decoy targets."""
    from sklearn.metrics import roc_auc_score

    s = scores if lower_better else -scores
    aucs = []
    for q in range(len(true_idx)):
        y = np.zeros(scores.shape[1])
        y[true_idx[q]] = 1
        aucs.append(roc_auc_score(y, -s[q]))  # lower score -> higher rank -> larger -s
    return float(np.mean(aucs))
