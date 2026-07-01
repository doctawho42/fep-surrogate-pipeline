"""Calibrated cycle-closure quality control (Fig L; the impact result).

A perturbation network must close its cycles: for any loop the signed sum of ΔΔG is zero in
expectation. Given the differentiable BAR bottleneck's calibrated per-edge aleatoric variance
``V_e = B_e / I_e**2`` (Theorem 2), cycle closure becomes a *properly-calibrated* statistical
test that separates **edge-level systematic error** (non-convergence, hysteresis, charge / water
/ protonation artifacts) from **sampling noise**:

    fit node potentials by GLS with weights ``1/V_e``; residual ``r_e``;
    ``X2 = sum_e r_e**2 / V_e ~ chi2(dof)``,  ``dof = #independent cycles = E - rank(incidence)``.

``reduced_chi2 >> 1`` ⇒ the network is over-dispersed relative to sampling ⇒ systematic error.
The test's usefulness is entirely a function of σ-calibration: an overconfident σ makes every
system fail (FPR → 1); only a calibrated σ (the sandwich) yields a selective test. Cycle closure
is (correctly) blind to node-consistent force-field bias, which cancels around a loop.

Pure NumPy; no SciPy/Torch dependency, so it imports fast and unit-tests run standalone.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

Edge = tuple[object, object, float, float]  # (node_a, node_b, ddg, se)


@dataclass
class NetworkFit:
    """Result of a GLS network fit + cycle-closure goodness-of-fit test."""

    nodes: list
    potentials: NDArray          # GLS node potentials (up to an additive gauge)
    residuals: NDArray           # r_e = y_e - (phi_b - phi_a)
    z: NDArray                   # standardized residuals r_e / sqrt(V_e)
    chi2: float                  # X^2 = sum_e r_e^2 / V_e
    dof: int                     # independent cycles = E - rank(incidence)

    @property
    def reduced_chi2(self) -> float:
        return self.chi2 / self.dof if self.dof > 0 else math.nan

    @property
    def p_value(self) -> float:
        """Upper-tail chi^2 p-value (Wilson--Hilferty; accurate in the tail we test)."""
        return chi2_sf(self.chi2, self.dof)


def _incidence(edges: list[Edge]):
    nodes = sorted({e[0] for e in edges} | {e[1] for e in edges}, key=str)
    idx = {n: i for i, n in enumerate(nodes)}
    E = len(edges)
    B = np.zeros((E, len(nodes)))
    y = np.zeros(E)
    V = np.zeros(E)
    for e, (a, b, ddg, se) in enumerate(edges):
        B[e, idx[b]] += 1.0
        B[e, idx[a]] -= 1.0
        y[e] = float(ddg)
        V[e] = float(se) ** 2
    return nodes, B, y, V


def gls_network(edges: list[Edge]) -> NetworkFit:
    """GLS fit of node potentials with inverse-sandwich-variance weights ``1/V_e``.

    The minimum-norm least-squares solution absorbs the (rank-1 per component) gauge freedom;
    the residuals and χ² are gauge-invariant. ``dof`` is the number of independent cycles.
    """
    if len(edges) < 1:
        raise ValueError("need at least one edge")
    nodes, B, y, V = _incidence(edges)
    if np.any(V <= 0):
        raise ValueError("edge variances must be positive (pass decorrelated, positive se)")
    Wh = 1.0 / np.sqrt(V)
    phi, *_ = np.linalg.lstsq(Wh[:, None] * B, Wh * y, rcond=None)
    r = y - B @ phi
    chi2 = float(np.sum(r * r / V))
    dof = int(len(edges) - np.linalg.matrix_rank(B))
    return NetworkFit(nodes=nodes, potentials=phi, residuals=r, z=r / np.sqrt(V),
                      chi2=chi2, dof=dof)


def chi2_sf(x: float, k: int) -> float:
    """Upper-tail survival function of chi^2(k) via the Wilson--Hilferty approximation."""
    if k <= 0:
        return math.nan
    z = ((x / k) ** (1.0 / 3.0) - (1.0 - 2.0 / (9 * k))) / math.sqrt(2.0 / (9 * k))
    return 0.5 * math.erfc(z / math.sqrt(2.0))


def cycle_closure_test(edges: list[Edge]) -> NetworkFit:
    """Calibrated cycle-closure goodness-of-fit test. Alias of :func:`gls_network` returning the
    full :class:`NetworkFit` (``reduced_chi2``, ``p_value``, per-edge ``z``)."""
    return gls_network(edges)


def benjamini_hochberg(pvals, alpha: float = 0.05) -> NDArray:
    """Boolean flags for BH-FDR control at level ``alpha``."""
    p = np.asarray(pvals, dtype=float)
    m = p.size
    order = np.argsort(p)
    thr = alpha * np.arange(1, m + 1) / m
    passed = np.where(p[order] <= thr)[0]
    cut = thr[passed.max()] if passed.size else 0.0
    return p <= cut


def repair_order(
    edges: list[Edge], target_reduced_chi2: float = 1.0, max_remove: int | None = None
):
    """Greedily remove the highest-|z| edge until ``reduced_chi2 <= target`` (or the graph is
    exhausted). Returns ``(removed_indices, trace)`` where ``trace`` is the reduced-χ² after each
    removal. Used to show that the *flagged* edges are the ones that cause the over-dispersion
    (guided removal reaches consistency far faster than random)."""
    remaining = list(range(len(edges)))
    removed: list[int] = []
    trace = [gls_network(edges).reduced_chi2]
    while len(remaining) > 1:
        sub = [edges[i] for i in remaining]
        fit = gls_network(sub)
        # A graph with no remaining cycles (dof == 0) is trivially consistent -> stop.
        if fit.dof == 0 or fit.reduced_chi2 <= target_reduced_chi2:
            break
        if max_remove is not None and len(removed) >= max_remove:
            break
        worst_local = int(np.argmax(np.abs(fit.z)))
        removed.append(remaining.pop(worst_local))
        trace.append(gls_network([edges[i] for i in remaining]).reduced_chi2)
    return removed, trace
