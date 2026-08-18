"""Theorem 5: the Hodge split of a perturbation network's error field.

A network's per-edge systematic error ``mu`` lives in ``R^E``. Whitening by the reported standard
errors and projecting onto the column space of the incidence matrix splits it, orthogonally and
uniquely, into

* a **gradient** part, ``mu_e = b_j - b_i`` for some per-ligand offsets ``b``, of dimension
  ``N - c``; and
* a **cycle** part, of dimension ``dof = E - N + c``.

Two facts make the split worth naming. First, the gradient part is *not identifiable*: for any
``delta``, the parameter pairs ``(phi, mu)`` and ``(phi + delta, mu - A delta)`` induce the same
distribution of the observed edge values, so no statistic computed from a network's own edges can
tell a per-ligand bias from no bias at all, at any magnitude. Second, the bias such an error
induces in the fitted node potentials depends *only* on the gradient part, while the closure
noncentrality depends *only* on the cycle part -- what a cycle-closure test can see and what
corrupts the answer are exact orthogonal complements.

The practical consequence is the ranking. Closure ranks edges by ``|z_e|``, which is powered by the
curl-leverage ``h_e``; the bias an edge's error induces scales with ``1 - h_e``. Ranking by
``|z_e| * sqrt((1 - h_e) / h_e)`` -- the classical DFFITS influence statistic, which the
conservation law ``h_e + w_e Omega_e = 1`` identifies with the un-auditable half of the budget --
targets harm rather than visibility.

Pure NumPy. Reuses ``bar.qc`` for the incidence and the GLS fit and ``bar.leverage`` for ``h_e``.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from bar.leverage import curl_leverage
from bar.qc import Edge, _incidence, gls_network

_H_TOL = 1e-9


@dataclass
class HodgeSplit:
    """An edge field decomposed into its gradient and cycle parts, in the original units."""

    gradient: NDArray            # the node-consistent part: invisible to closure, biases the fit
    cycle: NDArray               # the auditable part: visible to closure, biases nothing
    visible_fraction: float      # ||cycle||^2 / ||field||^2, both whitened
    dof: int
    rank: int


def _whitened(edges: list[Edge]):
    nodes, B, _y, V = _incidence(edges)
    if np.any(V <= 0):
        raise ValueError("edge variances must be positive")
    root_w = 1.0 / np.sqrt(V)
    return nodes, B, root_w, root_w[:, None] * B


def gradient_field(edges: list[Edge], potentials: dict) -> NDArray:
    """The edge field a set of per-ligand offsets induces: ``mu_e = potentials[b] - potentials[a]``.

    This is the shape every per-ligand systematic error takes, whether it comes from the force
    field, from the experimental assay, or from an annotation join.
    """
    return np.array([float(potentials[b]) - float(potentials[a]) for a, b, _, _ in edges])


def hodge_split(edges: list[Edge], field: NDArray) -> HodgeSplit:
    """Split ``field`` into its gradient and cycle parts, orthogonally in the whitened metric."""
    field = np.asarray(field, dtype=float)
    if field.shape != (len(edges),):
        raise ValueError(f"field must have one entry per edge, got {field.shape}")
    _nodes, B, root_w, xt = _whitened(edges)
    whitened = root_w * field
    hat = xt @ np.linalg.pinv(xt.T @ xt) @ xt.T
    grad_w = hat @ whitened
    cycle_w = whitened - grad_w
    total = float(whitened @ whitened)
    rank = int(np.linalg.matrix_rank(B))
    return HodgeSplit(
        gradient=grad_w / root_w,
        cycle=cycle_w / root_w,
        visible_fraction=float(cycle_w @ cycle_w) / total if total > 0 else 0.0,
        dof=len(edges) - rank,
        rank=rank,
    )


def gradient_r2(edges: list[Edge], field: NDArray) -> tuple[float, float, float]:
    """How much of ``field`` a per-ligand offset model explains: ``(R2, adjusted R2, chance)``.

    The model spends ``rank`` parameters on ``E`` observations and carries no intercept, so the
    no-intercept adjustment ``1 - (1 - R2) * E / (E - rank)`` has expectation zero for a field with
    no preferential alignment. ``chance`` is ``dof / E``, the visible fraction such a field shows.
    """
    split = hodge_split(edges, field)
    n_edges = len(edges)
    r2 = 1.0 - split.visible_fraction
    if n_edges <= split.rank:
        return r2, float("nan"), split.dof / n_edges
    adjusted = 1.0 - (1.0 - r2) * n_edges / (n_edges - split.rank)
    return r2, adjusted, split.dof / n_edges


def influence_rank(edges: list[Edge], z: NDArray | None = None) -> NDArray:
    """Per-edge influence ``|z_e| * sqrt((1 - h_e) / h_e)``; zero on un-auditable (bridge) edges.

    ``z`` defaults to the standardized residuals of the network's own GLS fit. An edge with
    ``h_e = 0`` carries no evidence at any magnitude, so it is given influence zero rather than the
    infinity the formula would otherwise produce: there is nothing to act on, not everything.
    """
    h = curl_leverage(edges)
    resid = gls_network(edges).z if z is None else np.asarray(z, dtype=float)
    auditable = h > _H_TOL
    out = np.zeros(len(edges))
    out[auditable] = np.abs(resid[auditable]) * np.sqrt((1.0 - h[auditable]) / h[auditable])
    return out


def influence_repair_order(edges: list[Edge], k: int) -> list[int]:
    """Indices of the ``k`` edges to remove, greedily by influence, refitting after each removal.

    Mirrors ``bar.qc.repair_order`` step for step and differs only in the statistic ranked, so a
    race between the two isolates the ranking. Stops early if no auditable edge remains.
    """
    remaining = list(range(len(edges)))
    removed: list[int] = []
    while len(removed) < k and len(remaining) > 1:
        sub = [edges[i] for i in remaining]
        infl = influence_rank(sub)
        if not np.any(infl > 0.0):
            break
        removed.append(remaining.pop(int(np.argmax(infl))))
    return removed


def unidentifiable_dimension(edges: list[Edge], anchored: set) -> int:
    """Dimension of the error subspace no measurement can resolve, given absolute anchors.

    Equals ``(N - |S|) - c_free(S)`` with ``S`` the anchored ligands and ``c_free`` the number of
    connected components holding none of them. The first anchor in a component removes only that
    component's arbitrary offset and buys no identifiability; each further anchor buys exactly one
    dimension; resolving a per-ligand bias in full needs an absolute measurement on every ligand.
    """
    nodes, _B, _rw, _xt = _whitened(edges)
    anchored = {n for n in anchored if n in set(nodes)}
    adjacency: dict = {n: set() for n in nodes}
    for a, b, _, _ in edges:
        adjacency[a].add(b)
        adjacency[b].add(a)
    seen: set = set()
    free_components = 0
    for start in nodes:
        if start in seen:
            continue
        stack, component = [start], set()
        while stack:
            node = stack.pop()
            if node in component:
                continue
            component.add(node)
            stack.extend(adjacency[node] - component)
        seen |= component
        if not (component & anchored):
            free_components += 1
    return (len(nodes) - len(anchored)) - free_components
