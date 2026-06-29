"""Chirality completeness readout (Theorem 4; invariants #5/#6).

``even_features`` are O(3)-invariant (pairwise distances) and therefore *constant on
an enantiomer pair* {M, sigma M} — a distance-only / norm-and-dot-product readout is
provably enantiomer-blind. ``signed_volume`` is the parity-odd ``0o`` pseudoscalar
chi(M)=det[p1-p0,p2-p0,p3-p0]=v1.(v2 x v3): SO(3)-invariant, sign-flipping under
reflection. ``(pairwise distances, sgn chi)`` is a complete SO(3) invariant — O(3)
loses exactly the chirality bit, and the 0o channel restores it.
"""
from __future__ import annotations

from itertools import combinations

import numpy as np
from numpy.typing import ArrayLike, NDArray


def even_features(coords: ArrayLike) -> NDArray:
    """O(3)-invariant features: all pairwise distances (fixed i<j order)."""
    c = np.asarray(coords, dtype=float)
    return np.array([np.linalg.norm(c[i] - c[j]) for i, j in combinations(range(c.shape[0]), 2)])


def signed_volume(coords: ArrayLike) -> float:
    """Parity-odd 0o pseudoscalar chi = det[p1-p0, p2-p0, p3-p0] (first 4 points)."""
    c = np.asarray(coords, dtype=float)
    v = (c[1:4] - c[0]).T  # columns v1, v2, v3
    return float(np.linalg.det(v))


def chiral_readout(coords: ArrayLike, include_0o: bool = True) -> NDArray:
    """Readout = [even O(3) features] (+ [0o pseudoscalar] if include_0o).
    Without 0o the readout collapses enantiomers; with it, they separate."""
    feats = even_features(coords)
    if include_0o:
        feats = np.append(feats, signed_volume(coords))
    return feats
