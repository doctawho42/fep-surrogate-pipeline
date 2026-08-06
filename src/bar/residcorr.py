"""Residual-correlation diagnostic and effective dof (peer-review item P2).

The manuscript reports two conservatism factors that do not agree: the closure-implied factor
(median reduced chi^2 = 0.34, i.e. ~1.71x in se) and the replicate-validated factor (1.25-1.41x).
The referees' proposed explanation is that the GLS fit assumes a DIAGONAL V, while residuals of
edges sharing a ligand endpoint are correlated; correlated errors change E[X^2] away from the
nominal dof, deflating the reduced chi^2 independently of how wide the bars are.

The load-bearing subtlety: residuals are correlated even under a PERFECT null, because
``r = M eps`` with ``M = I - H`` the residual-maker projector. So the empirical correlation must
be compared against the correlation M itself induces (``null_pair_correlation``), never against
zero. The excess over that null is the evidence for genuine error correlation.

Effective dof uses ``E[X^2] = tr(M C)`` with ``C`` the correlation of the *whitened errors*: under
independence ``C = I`` and ``E[X^2] = tr(M) = dof``. We plug in the structured estimate
``C = I + rho_shared * S + rho_disjoint * D`` with the EXCESS correlations, which approximates the
error correlation by the excess residual correlation. That approximation is first-order (exact
only when the M-induced coupling for those pairs is small) and must be stated wherever the number
is reported.
"""
from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from bar.qc import Edge, _incidence


def residual_maker(edges: list[Edge]) -> NDArray:
    """The whitened residual-maker ``M = I - H``, an orthogonal projector onto the cycle space
    (``bar.leverage.curl_leverage`` returns its diagonal; here we need the full matrix)."""
    if len(edges) < 1:
        raise ValueError("need at least one edge")
    _nodes, B, _y, V = _incidence(edges)
    if np.any(V <= 0):
        raise ValueError("edge variances must be positive")
    xt = (1.0 / np.sqrt(V))[:, None] * B
    u, s, _ = np.linalg.svd(xt, full_matrices=False)
    tol = float(s.max()) * max(xt.shape) * np.finfo(float).eps if s.size else 0.0
    r = int(np.sum(s > tol))
    return np.eye(len(edges)) - u[:, :r] @ u[:, :r].T


def pair_masks(edges: list[Edge]) -> tuple[NDArray, NDArray]:
    """Upper-triangular boolean masks ``(shared_node, disjoint)`` over edge pairs."""
    E = len(edges)
    shared = np.zeros((E, E), dtype=bool)
    ends = [{a, b} for a, b, _y, _se in edges]
    for i in range(E):
        for j in range(i + 1, E):
            shared[i, j] = bool(ends[i] & ends[j])
    iu = np.triu(np.ones((E, E), dtype=bool), 1)
    return shared & iu, (~shared) & iu


def null_pair_correlation(M: NDArray, mask: NDArray) -> float:
    """Residual correlation the projector induces under a perfect null, over ``mask`` pairs.

    Uses the SAME ratio-of-means form as :func:`empirical_pair_correlation`, namely
    ``mean_{(e,f)} M_ef / mean_e M_ee``, because under the null ``Cov(r) = M`` gives
    ``E[z_e z_f] = M_ef`` and ``E[z_e^2] = M_ee``, so that is exactly what the empirical
    estimator converges to. A mean-of-ratios form (``mean of M_ef/sqrt(M_ee M_ff)``) is a
    DIFFERENT estimator and disagrees whenever the per-edge leverage is heterogeneous, which
    would make the empirical-vs-null comparison invalid.
    """
    ii, jj = np.where(mask)
    if ii.size == 0:
        return 0.0
    denom = float(np.mean(np.diag(M)))
    if abs(denom) < 1e-12:
        return 0.0
    return float(np.mean(M[ii, jj]) / denom)


def empirical_pair_correlation(z_reps: NDArray, mask: NDArray) -> float:
    """Pooled empirical residual correlation over ``mask`` pairs.

    ``z_reps`` is ``(n_reps, E)`` standardized residuals from independently fitted replicates.
    Uses the ratio-of-means form ``mean_{(e,f),k} z_e z_f / mean_{e,k} z_e^2``, which is far more
    stable at small ``n_reps`` than averaging per-pair correlations.
    """
    z = np.asarray(z_reps, dtype=float)
    ii, jj = np.where(mask)
    if ii.size == 0:
        return 0.0
    var = float(np.mean(z * z))
    if var <= 0:
        return 0.0
    return float(np.mean(z[:, ii] * z[:, jj]) / var)


def effective_dof(
    M: NDArray, shared: NDArray, disjoint: NDArray, rho_shared: float, rho_disjoint: float
) -> float:
    """``tr(M C)`` with ``C = I + rho_shared * S + rho_disjoint * D`` (S, D symmetrized).

    Equals ``tr(M) = dof`` when both excess correlations are zero. Feed the EXCESS (empirical
    minus null) correlations, not the raw empirical ones.
    """
    off = 2.0 * (rho_shared * float(np.sum(M[shared])) + rho_disjoint * float(np.sum(M[disjoint])))
    return float(np.trace(M)) + off
