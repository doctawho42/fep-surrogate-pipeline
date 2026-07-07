"""Estimation-detection conservation law and per-edge observability map (D1).

The GLS cycle-closure fit has a residual-maker projector ``M = I - H`` onto the whitened cycle
space. Its diagonal ``h_e = M_ee`` (the *curl-leverage*) is how much of edge ``e``'s error survives
into the closure residual, i.e. how observable edge ``e`` is to the QC. Theorem D1 (conservation
law): ``h_e + w_e * Omega_e = 1`` with ``w_e = 1/V_e`` the sandwich conductance and ``Omega_e`` the
effective resistance of the edge's endpoints (Theorem 3's object) -- estimation self-influence and
detection leverage are exactly complementary. Corollaries: ``sum_e h_e = dof`` (independent cycles),
and the closure-chi^2 noncentrality against a systematic shift ``mu`` on edge ``e`` is
``h_e * mu^2 / V_e`` (so the minimax-detectable shift is ``delta*_e = sqrt(V_e / h_e)``). A bridge
edge lies in no cycle, so ``h_e = 0`` and it is structurally un-auditable; node-consistent
force-field bias is a gradient in ``range(H)``, annihilated by ``M`` -- which *derives* the QC's
node-bias blindness rather than asserting it.

Pure NumPy; SciPy is used lazily only for the detectability constant. Reuses ``bar.qc`` and
``bar.graph``; ``h_e`` is the *residual-maker* diagonal, never the fit-hat diagonal.
"""
from __future__ import annotations

import math
from collections import defaultdict

import numpy as np
from numpy.typing import NDArray

from bar.graph import effective_resistance, weighted_laplacian
from bar.qc import Edge, _incidence


def curl_leverage(edges: list[Edge], tol: float = 1e-9) -> NDArray:
    """Per-edge curl-leverage ``h_e = (I - H)_ee``, computed two independent ways and cross-checked.

    Path A (projector): SVD hat matrix of the whitened incidence ``X_tilde = sqrt(w) * B``.
    Path B (dual, Theorem D1): ``h_e = 1 - w_e * Omega_e`` via the weighted Laplacian and effective
    resistance. Their agreement to ``tol`` is the numerical proof of the conservation law. Raises
    ``ValueError`` on mismatch (a wrong orientation/weight/index would break it).
    """
    if len(edges) < 1:
        raise ValueError("need at least one edge")
    nodes, B, _y, V = _incidence(edges)
    if np.any(V <= 0):
        raise ValueError("edge variances must be positive")
    w = 1.0 / V
    idx = {n: i for i, n in enumerate(nodes)}

    # Path A: residual-maker diagonal via an SVD projector onto col(X_tilde).
    xt = np.sqrt(w)[:, None] * B
    u, s, _ = np.linalg.svd(xt, full_matrices=False)
    s_tol = float(s.max()) * max(xt.shape) * np.finfo(float).eps if s.size else 0.0
    r = int(np.sum(s > s_tol))
    hat = u[:, :r] @ u[:, :r].T
    h_proj = 1.0 - np.diag(hat)

    # Path B: 1 - w_e * Omega_e via the weighted Laplacian (Theorem 3 object).
    lap = weighted_laplacian(
        [(idx[a], idx[b], float(w[e])) for e, (a, b, _, _) in enumerate(edges)], len(nodes)
    )
    h_dual = np.array(
        [1.0 - float(w[e]) * effective_resistance(lap, idx[a], idx[b])
         for e, (a, b, _, _) in enumerate(edges)]
    )

    dev = float(np.max(np.abs(h_proj - h_dual))) if len(edges) else 0.0
    if dev > tol:
        raise ValueError(f"conservation law h + w*Omega = 1 violated: max deviation {dev:.2e}")
    return h_proj


def bridges(edges: list[Edge]) -> set[int]:
    """Indices of bridge edges (edges in no cycle) via Tarjan DFS, tracking the parent *edge* id so
    parallel edges are correctly not-bridges."""
    adj: dict[object, list[tuple[object, int]]] = defaultdict(list)
    for i, (a, b, _, _) in enumerate(edges):
        adj[a].append((b, i))
        adj[b].append((a, i))
    disc: dict[object, int] = {}
    low: dict[object, int] = {}
    timer = [0]
    out: set[int] = set()

    def dfs(u: object, parent_edge: int) -> None:
        disc[u] = low[u] = timer[0]
        timer[0] += 1
        for v, ei in adj[u]:
            if ei == parent_edge:
                continue
            if v not in disc:
                dfs(v, ei)
                low[u] = min(low[u], low[v])
                if low[v] > disc[u]:
                    out.add(ei)
            else:
                low[u] = min(low[u], disc[v])

    for node in list(adj):
        if node not in disc:
            dfs(node, -1)
    return out


def _lambda_star(alpha: float, power: float) -> float:
    """Noncentrality of a 1-dof noncentral-chi^2 test achieving ``power`` at level ``alpha``."""
    from scipy.optimize import brentq
    from scipy.stats import chi2, ncx2

    crit = float(chi2.ppf(1.0 - alpha, 1))
    return float(brentq(lambda lam: float(ncx2.sf(crit, 1, lam)) - power, 0.0, 1000.0))


def observability_certificate(
    edges: list[Edge], alpha: float = 0.05, power: float = 0.8, h_min: float = 0.05
) -> list[dict]:
    """Per-edge observability record: leverage ``h``, estimation share ``w*Omega`` (= ``1 - h``),
    detectable-shift resolution ``delta_star = sqrt(lambda* * V / h)`` (inf for bridges / ``h<=0``),
    bridge flag, and auditability (``h >= h_min``)."""
    h = curl_leverage(edges)
    br = bridges(edges)
    _nodes, _B, _y, V = _incidence(edges)
    lam = _lambda_star(alpha, power)
    out: list[dict] = []
    for e, (a, b, _ddg, _se) in enumerate(edges):
        he = float(h[e])
        is_bridge = e in br
        delta_star = math.inf if (is_bridge or he <= 0.0) else math.sqrt(lam * float(V[e]) / he)
        out.append({
            "index": e, "node_a": a, "node_b": b, "h": he, "w_times_Omega": 1.0 - he,
            "V": float(V[e]), "delta_star": delta_star, "is_bridge": is_bridge,
            "auditable": he >= h_min,
        })
    return out
