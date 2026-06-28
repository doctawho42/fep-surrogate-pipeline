"""Fisher-resistance graph utilities (Theorem 3).

Relative free-energy measurements ``y_e = (phi_j - phi_i) + eps_e`` with
``eps_e ~ N(0, V_e)`` and ``V_e = B_e/I_e^2`` the per-edge sandwich variance. The
Gaussian posterior precision of the node potentials is the weighted Laplacian
``L = sum_e w_e (e_i - e_j)(e_i - e_j)^T`` with conductances ``w_e = 1/V_e =
I_e^2/B_e`` (invariant #4 — inverse sandwich variances, NOT raw overlaps). Then
contrast variance = effective resistance, single edge -> V_e, Sherman-Morrison
gives O(1) updates, and ``ker L = span{component indicators}`` so gauge directions
have zero acquisition value.
"""
from __future__ import annotations

from collections.abc import Iterable, Mapping

import numpy as np
from numpy.typing import NDArray


def weighted_laplacian(edges, n: int) -> NDArray:
    """Build the n-by-n weighted Laplacian. ``edges`` is either a mapping
    ``{(i, j): w}`` or an iterable of ``(i, j, w)`` triples."""
    items: Iterable
    if isinstance(edges, Mapping):
        items = ((i, j, w) for (i, j), w in edges.items())
    else:
        items = edges
    L = np.zeros((n, n), dtype=float)
    for i, j, w in items:
        L[i, i] += w
        L[j, j] += w
        L[i, j] -= w
        L[j, i] -= w
    return L


def effective_resistance(L: NDArray, a: int, b: int) -> float:
    """Omega_ab = (e_a - e_b)^T L^+ (e_a - e_b), the variance of the contrast
    phi_b - phi_a (Theorem 3(i)). Uses the Moore-Penrose pseudoinverse (gauge-safe)."""
    Lp = np.linalg.pinv(L)
    e = np.zeros(L.shape[0])
    e[a], e[b] = 1.0, -1.0
    return float(e @ Lp @ e)


def laplacian_nullity(L: NDArray, tol: float = 1e-9) -> int:
    """Dimension of ker L = number of connected components (Theorem 3(iv))."""
    evals = np.linalg.eigvalsh(L)
    scale = max(1.0, float(np.max(np.abs(evals))))
    return int(np.sum(np.abs(evals) < tol * scale))


def sherman_morrison_resistance(omega: float, g: float) -> float:
    """O(1) update after adding a measurement of precision g along the same
    contrast: Omega' = Omega / (1 + g*Omega) (Theorem 3(iii))."""
    return omega / (1.0 + g * omega)


def edge_weight_from_sandwich(I_e: float, B_e: float) -> float:
    """Laplacian conductance w_e = I_e^2 / B_e = 1 / V_e (invariant #4)."""
    return I_e**2 / B_e
