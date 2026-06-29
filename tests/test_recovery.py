"""Recovery-metric invariants for the target-finding gate."""
from __future__ import annotations

import numpy as np

from screen.recovery import recovery_at_k, recovery_auroc


def test_perfect_ranker_recovers_at_1():
    # true target has the lowest (best) score for every query
    scores = np.array([[0.1, 5.0, 5.0], [5.0, 0.2, 5.0], [5.0, 5.0, 0.3]])
    true = np.array([0, 1, 2])
    rec = recovery_at_k(scores, true, lower_better=True)
    assert rec[0] == 1.0
    assert recovery_auroc(scores, true) == 1.0


def test_worst_ranker_recovers_only_at_k_equals_n():
    # true target has the WORST score every time -> recovered only at k=n
    scores = np.array([[9.0, 0.0, 1.0], [1.0, 9.0, 0.0]])
    true = np.array([0, 1])
    rec = recovery_at_k(scores, true, lower_better=True)
    assert rec[0] == 0.0
    assert rec[-1] == 1.0          # everyone is within top-n
    assert recovery_auroc(scores, true) == 0.0


def test_recovery_monotone_nondecreasing_in_k():
    rng = np.random.default_rng(0)
    scores = rng.normal(size=(20, 6))
    true = rng.integers(0, 6, 20)
    rec = recovery_at_k(scores, true)
    assert np.all(np.diff(rec) >= 0)
    assert rec[-1] == 1.0


def test_random_recovery_near_uniform():
    rng = np.random.default_rng(1)
    n_t = 8
    scores = rng.normal(size=(4000, n_t))
    true = rng.integers(0, n_t, 4000)
    rec = recovery_at_k(scores, true)
    assert abs(rec[0] - 1 / n_t) < 0.03       # top-1 ~ 1/n
    assert abs(recovery_auroc(scores, true) - 0.5) < 0.03
