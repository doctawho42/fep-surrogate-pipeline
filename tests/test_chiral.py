"""Chirality completeness (Theorem 4, invariants #5/#6).

O(3)-invariant ("even") readouts collapse enantiomers; the parity-odd 0o
pseudoscalar (signed volume / triple product) separates them. Pairwise distances +
sgn(chi) is a complete SO(3) invariant — exactly one lost bit, restored by 0o.
"""
from __future__ import annotations

import numpy as np
import pytest
from scipy.spatial.transform import Rotation

from bar.chiral import (
    chiral_readout,
    even_features,
    signed_volume,
)

REFLECT = np.diag([-1.0, 1.0, 1.0])


def _tetra(seed: int) -> np.ndarray:
    """A random chiral 4-point configuration (generic -> non-coplanar)."""
    return np.random.default_rng(seed).normal(size=(4, 3))


def _mirror(coords: np.ndarray) -> np.ndarray:
    return coords @ REFLECT.T


# --- even (O(3)-invariant) features --------------------------------------
def test_even_features_invariant_under_reflection():
    M = _tetra(0)
    np.testing.assert_allclose(even_features(M), even_features(_mirror(M)), atol=1e-12)


def test_even_features_invariant_under_rotation():
    M = _tetra(1)
    R = Rotation.random(random_state=7).as_matrix()
    np.testing.assert_allclose(even_features(M), even_features(M @ R.T), atol=1e-10)


# --- 0o pseudoscalar ------------------------------------------------------
def test_signed_volume_equals_triple_product():
    M = _tetra(2)
    v1, v2, v3 = M[1] - M[0], M[2] - M[0], M[3] - M[0]
    assert signed_volume(M) == pytest.approx(v1 @ np.cross(v2, v3), rel=1e-12)


def test_signed_volume_flips_under_reflection():
    M = _tetra(3)
    assert signed_volume(_mirror(M)) == pytest.approx(-signed_volume(M), rel=1e-12)


def test_signed_volume_invariant_under_rotation():
    M = _tetra(4)
    R = Rotation.random(random_state=11).as_matrix()
    assert signed_volume(M @ R.T) == pytest.approx(signed_volume(M), rel=1e-9)


def test_chi_of_linear_map_scales_by_det():
    M = _tetra(5)
    A = np.random.default_rng(9).normal(size=(3, 3))
    assert signed_volume(M @ A.T) == pytest.approx(np.linalg.det(A) * signed_volume(M), rel=1e-9)


# --- the theorem: even collapses, 0o separates ----------------------------
def test_even_readout_collapses_enantiomers():
    M = _tetra(6)
    np.testing.assert_allclose(
        chiral_readout(M, include_0o=False),
        chiral_readout(_mirror(M), include_0o=False),
        atol=1e-12,
    )


def test_0o_readout_separates_enantiomers():
    M = _tetra(7)
    rM = chiral_readout(M, include_0o=True)
    rMp = chiral_readout(_mirror(M), include_0o=True)
    assert not np.allclose(rM, rMp)
    # they differ ONLY in the 0o channel, which flips sign
    np.testing.assert_allclose(rM[:-1], rMp[:-1], atol=1e-12)
    assert rM[-1] == pytest.approx(-rMp[-1], rel=1e-12)


def test_coplanar_pseudoscalar_is_zero():
    # achiral (coplanar) config: chi = 0, no enantiomer to separate
    flat = np.array([[0.0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0]])
    assert signed_volume(flat) == pytest.approx(0.0, abs=1e-12)
