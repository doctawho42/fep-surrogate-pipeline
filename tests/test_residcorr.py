"""Invariants for the residual-correlation / effective-dof diagnostic (peer-review item P2)."""
from __future__ import annotations

import numpy as np
import pytest

from bar.qc import gls_network
from bar.residcorr import (
    effective_dof,
    empirical_pair_correlation,
    null_pair_correlation,
    pair_masks,
    residual_maker,
)


def _net():
    """Two triangles sharing node C, plus a chord: 5 nodes, 7 edges, dof = 3."""
    return [("A", "B", 0.3, 0.2), ("B", "C", 0.4, 0.25), ("A", "C", 0.7, 0.3),
            ("C", "D", 0.2, 0.15), ("D", "E", 0.5, 0.2), ("C", "E", 0.7, 0.35),
            ("A", "D", 0.9, 0.4)]


def test_residual_maker_is_a_projector_with_trace_dof():
    M = residual_maker(_net())
    assert np.allclose(M, M @ M, atol=1e-9)          # idempotent
    assert np.allclose(M, M.T, atol=1e-12)           # symmetric
    assert abs(np.trace(M) - gls_network(_net()).dof) < 1e-9


def test_pair_masks_partition_the_upper_triangle():
    edges = _net()
    shared, disjoint = pair_masks(edges)
    E = len(edges)
    iu = np.triu(np.ones((E, E), dtype=bool), 1)
    assert not np.any(shared & disjoint)
    assert np.array_equal(shared | disjoint, iu)
    # ("A","B") and ("B","C") share node B -> shared; ("A","B") and ("D","E") share none
    assert shared[0, 1] and disjoint[0, 4]


def test_null_correlation_is_nonzero_by_construction():
    """The residual-maker induces correlation even under a perfect null -- this is exactly the
    confound the empirical estimate must be compared against, not against zero."""
    shared, _ = pair_masks(_net())
    assert abs(null_pair_correlation(residual_maker(_net()), shared)) > 1e-3


def test_empirical_matches_null_when_errors_are_independent():
    """Simulate the exact null (independent N(0, V_e) edge errors) and confirm the empirical
    pair correlation converges to the M-induced null value."""
    edges = _net()
    M = residual_maker(edges)
    shared, _ = pair_masks(edges)
    rng = np.random.default_rng(0)
    se = np.array([e[3] for e in edges])
    zs = []
    for _ in range(4000):
        eps = rng.standard_normal(len(edges)) * se
        pert = [(a, b, y + d, s) for (a, b, y, s), d in zip(edges, eps, strict=True)]
        zs.append(gls_network(pert).z)
    emp = empirical_pair_correlation(np.vstack(zs), shared)
    assert abs(emp - null_pair_correlation(M, shared)) < 0.03


def test_effective_dof_equals_dof_when_errors_are_uncorrelated():
    edges = _net()
    M = residual_maker(edges)
    shared, disjoint = pair_masks(edges)
    assert abs(effective_dof(M, shared, disjoint, 0.0, 0.0) - np.trace(M)) < 1e-9


def test_effective_dof_moves_with_excess_correlation():
    edges = _net()
    M = residual_maker(edges)
    shared, disjoint = pair_masks(edges)
    base = effective_dof(M, shared, disjoint, 0.0, 0.0)
    moved = effective_dof(M, shared, disjoint, 0.3, 0.0)
    assert moved != base
    # exact value: tr(M) + 2 * rho_shared * sum(M[shared]), rho_disjoint=0 contributes nothing --
    # catches a sign flip or a factor-of-2 regression in the off-diagonal symmetrization.
    assert moved == pytest.approx(float(np.trace(M)) + 2 * 0.3 * float(M[shared].sum()))
