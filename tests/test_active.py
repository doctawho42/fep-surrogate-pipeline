"""Active-learning core over the FEP-edge graph (Theorem 3 in action).

Gaussian belief over node potentials (ΔG) from relative measurements; gauge-aware,
cost-aware KG = Sherman-Morrison reduction of decision-contrast variance.
"""
from __future__ import annotations

import numpy as np
import pytest

from bar.active import BeliefGraph, kg_scores


def test_measurement_reduces_own_contrast_variance():
    bg = BeliefGraph(4, prior_precision=1e-3)
    before = bg.contrast_variance(0, 1)
    bg.add_measurement(0, 1, y=2.0, precision=5.0)
    assert bg.contrast_variance(0, 1) < before


def test_posterior_recovers_contrasts_up_to_gauge():
    # true potentials; measure a connected path with low noise -> contrasts ~ truth
    true = np.array([0.0, -1.0, 2.0, 0.5])
    bg = BeliefGraph(4, prior_precision=1e-6)
    for i, j in [(0, 1), (1, 2), (2, 3)]:
        bg.add_measurement(i, j, y=true[j] - true[i], precision=1e4)
    m = bg.mean
    for a, b in [(0, 1), (0, 2), (1, 3)]:
        assert (m[b] - m[a]) == pytest.approx(true[b] - true[a], abs=1e-2)


def test_variance_reduction_matches_direct_recompute():
    bg = BeliefGraph(5, prior_precision=1e-2)
    for i, j, y in [(0, 1, 1.0), (1, 2, 0.5), (3, 4, -1.0)]:
        bg.add_measurement(i, j, y, precision=2.0)
    contrasts = [(0, 2), (1, 4)]
    g = 3.0
    sm = bg.variance_reduction(0, 3, precision=g, contrasts=contrasts)
    # direct: clone, add the measurement, recompute contrast variances
    before = np.array([bg.contrast_variance(a, b) for a, b in contrasts])
    bg2 = bg.copy()
    bg2.add_measurement(0, 3, y=0.0, precision=g)
    after = np.array([bg2.contrast_variance(a, b) for a, b in contrasts])
    np.testing.assert_allclose(sm, before - after, rtol=1e-9, atol=1e-12)


def test_kg_zero_for_edge_irrelevant_to_decision_contrast():
    # decision cares about contrast (0,1); a candidate edge in the disjoint {2,3}
    # part cannot reduce var(phi_1 - phi_0) -> KG ~ 0 for that edge.
    bg = BeliefGraph(4, prior_precision=1e-3)
    bg.add_measurement(0, 1, 1.0, precision=1.0)
    bg.add_measurement(2, 3, 1.0, precision=1.0)
    cand = [(0, 1), (2, 3)]
    scores = kg_scores(bg, cand, cand_precision=[5.0, 5.0], costs=[1.0, 1.0],
                       contrasts=[(0, 1)], weights=[1.0])
    assert scores[1] == pytest.approx(0.0, abs=1e-12)
    assert scores[0] > 1e-6


def test_cost_awareness_scales_score():
    bg = BeliefGraph(3, prior_precision=1e-2)
    cand = [(0, 1)]
    s_cheap = kg_scores(bg, cand, [4.0], costs=[1.0], contrasts=[(0, 1)], weights=[1.0])[0]
    s_dear = kg_scores(bg, cand, [4.0], costs=[4.0], contrasts=[(0, 1)], weights=[1.0])[0]
    assert s_cheap == pytest.approx(4.0 * s_dear, rel=1e-9)


def test_global_level_is_gauge_free():
    # connected measured graph + tiny prior: a contrast has finite variance while the
    # all-ones (gauge) level is nearly unconstrained (variance ~ 1/tau, much larger).
    bg = BeliefGraph(4, prior_precision=1e-8)
    for i, j in [(0, 1), (1, 2), (2, 3)]:
        bg.add_measurement(i, j, 1.0, precision=10.0)
    cov = bg.cov
    ones = np.ones(4) / 2.0
    gauge_var = float(ones @ cov @ ones)
    assert bg.contrast_variance(0, 3) < 1e-3 * gauge_var
