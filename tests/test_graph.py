"""Fisher-resistance graph utilities (Theorem 3). Invariant #4: Laplacian edge
weights are inverse sandwich variances w_e = I_e^2 / B_e, NOT raw overlaps."""
from __future__ import annotations

import numpy as np
import pytest

from bar.graph import (
    edge_weight_from_sandwich,
    effective_resistance,
    laplacian_nullity,
    sherman_morrison_resistance,
    weighted_laplacian,
)

# triangle nodes 0,1,2; conductances w(0,1)=1, w(0,2)=0.5, w(1,2)=2
TRIANGLE = {(0, 1): 1.0, (0, 2): 0.5, (1, 2): 2.0}


def test_effective_resistance_series_parallel():
    # Omega_12 = direct(1) parallel [series 1/0.5 + 1/2] = 1/(1 + 0.4) = 1/1.4
    L = weighted_laplacian(TRIANGLE, 3)
    assert effective_resistance(L, 0, 1) == pytest.approx(0.714286, abs=1e-5)
    assert effective_resistance(L, 0, 1) == pytest.approx(1 / 1.4, rel=1e-12)


def test_single_edge_resistance_equals_variance():
    # one edge with conductance w -> Omega = 1/w = V_e (Theorem 3(ii))
    L = weighted_laplacian({(0, 1): 4.0}, 2)
    assert effective_resistance(L, 0, 1) == pytest.approx(0.25, rel=1e-12)


def test_kernel_is_constant_vector():
    # connected graph: nullity 1, null vector proportional to 1 (Theorem 3(iv))
    L = weighted_laplacian(TRIANGLE, 3)
    assert laplacian_nullity(L) == 1
    evals, evecs = np.linalg.eigh(L)
    null = evecs[:, 0]
    cos = abs(null @ np.ones(3)) / (np.linalg.norm(null) * np.sqrt(3))
    assert cos == pytest.approx(1.0, abs=1e-9)


def test_disconnected_graph_nullity():
    # two components {0,1} and {2,3}: nullity 2 (ker L = span of component indicators)
    L = weighted_laplacian({(0, 1): 1.0, (2, 3): 1.0}, 4)
    assert laplacian_nullity(L) == 2


def test_sherman_morrison_update():
    # add measurement along (0,1) with precision g=0.7 -> Omega' = Omega/(1+g Omega)
    L = weighted_laplacian(TRIANGLE, 3)
    omega = effective_resistance(L, 0, 1)
    omega_sm = sherman_morrison_resistance(omega, 0.7)
    assert omega_sm == pytest.approx(0.476190, abs=1e-5)
    # exact agreement with a direct Laplacian recompute (add g to that edge)
    L2 = weighted_laplacian({(0, 1): 1.0 + 0.7, (0, 2): 0.5, (1, 2): 2.0}, 3)
    assert effective_resistance(L2, 0, 1) == pytest.approx(omega_sm, rel=1e-9)


def test_edge_weight_is_inverse_sandwich():
    # invariant #4: w_e = I_e^2 / B_e = 1 / V_e  (V_e = B_e/I_e^2 the sandwich)
    I_e, B_e = 0.6, 0.15
    V_e = B_e / I_e**2
    assert edge_weight_from_sandwich(I_e, B_e) == pytest.approx(1.0 / V_e, rel=1e-12)
    assert edge_weight_from_sandwich(I_e, B_e) == pytest.approx(I_e**2 / B_e, rel=1e-12)
