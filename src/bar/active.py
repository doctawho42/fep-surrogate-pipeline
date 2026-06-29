"""Active learning over the FEP-edge graph (Theorem 3 in action).

A Gaussian belief over node potentials phi (the per-ligand ΔG) built from relative
measurements ``y_e ≈ phi_j - phi_i`` with precision ``w_e = 1/V_e = I_e^2/B_e`` (the
inverse sandwich variance, invariant #4). The precision matrix is the weighted
Laplacian plus a vague prior; contrasts ``phi_b - phi_a`` are gauge-invariant.

Acquisition = gauge-aware, cost-aware Knowledge Gradient: the decision-weighted
reduction of decision-contrast variance from a candidate measurement, computed in
O(1) per (contrast, candidate) by Sherman-Morrison (Theorem 3(iii)). Edges that
cannot reduce any decision contrast (gauge / redundant directions) score 0 — the
acquisition is automatically gauge-invariant (Theorem 3(iv)).
"""
from __future__ import annotations

import numpy as np
from numpy.typing import NDArray


class BeliefGraph:
    """Gaussian posterior over node potentials from relative measurements.

    Edge ``(i, j)`` measures ``y ≈ phi_j - phi_i`` (contrast vector ``e_j - e_i``).
    Precision ``L = tau*I + sum_e w_e b_e b_e^T``; information vector ``h``; mean
    solves ``L mean = h``; covariance ``L^{-1}``.
    """

    def __init__(self, n_nodes: int, prior_precision: float = 1e-3):
        self.n = int(n_nodes)
        self.tau = float(prior_precision)
        self.L = self.tau * np.eye(self.n)
        self.h = np.zeros(self.n)

    def _b(self, i: int, j: int) -> NDArray:
        b = np.zeros(self.n)
        b[j], b[i] = 1.0, -1.0
        return b

    def add_measurement(self, i: int, j: int, y: float, precision: float) -> None:
        b = self._b(i, j)
        self.L += precision * np.outer(b, b)
        self.h += precision * y * b

    @property
    def mean(self) -> NDArray:
        return np.linalg.solve(self.L, self.h)

    @property
    def cov(self) -> NDArray:
        return np.linalg.inv(self.L)

    def contrast_variance(self, a: int, b: int) -> float:
        d = np.zeros(self.n)
        d[b], d[a] = 1.0, -1.0
        cov = self.cov
        return float(d @ cov @ d)

    def variance_reduction(self, i: int, j: int, precision: float, contrasts) -> NDArray:
        """Sherman-Morrison reduction in each contrast's variance from adding the
        measurement (i, j) with the given precision. O(1) per contrast."""
        cov = self.cov
        bc = self._b(i, j)
        cov_bc = cov @ bc
        denom = 1.0 + precision * float(bc @ cov_bc)
        out = np.empty(len(contrasts))
        for k, (a, b) in enumerate(contrasts):
            d = np.zeros(self.n)
            d[b], d[a] = 1.0, -1.0
            num = float(d @ cov_bc) ** 2
            out[k] = precision * num / denom
        return out

    def copy(self) -> BeliefGraph:
        new = BeliefGraph(self.n, self.tau)
        new.L = self.L.copy()
        new.h = self.h.copy()
        return new


def kg_scores(belief, candidates, cand_precision, costs, contrasts, weights) -> NDArray:
    """Decision-weighted, cost-aware KG score for each candidate edge:
    ``sum_d weight_d * variance_reduction_d(c) / cost_c``. Vectorised: one covariance
    solve, then batch Sherman-Morrison over all (candidate, contrast) pairs."""
    n = belief.n
    cov = belief.cov
    w = np.asarray(weights, dtype=float)
    g = np.asarray(cand_precision, dtype=float)
    cst = np.asarray(costs, dtype=float)

    # incidence matrices: candidates Bc (m x n), contrasts Dd (p x n)
    Bc = np.zeros((len(candidates), n))
    for c, (i, j) in enumerate(candidates):
        Bc[c, j], Bc[c, i] = 1.0, -1.0
    Dd = np.zeros((len(contrasts), n))
    for k, (a, b) in enumerate(contrasts):
        Dd[k, b], Dd[k, a] = 1.0, -1.0

    covBc = cov @ Bc.T                       # n x m
    denom = 1.0 + g * np.einsum("cn,nc->c", Bc, covBc)   # m
    cross = Dd @ covBc                       # p x m  = (d . cov . b_c)
    red = g[None, :] * cross**2 / denom[None, :]         # p x m
    return (w @ red) / cst                    # m
