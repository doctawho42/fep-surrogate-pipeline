"""Invariants for the estimation-detection conservation law + observability map (D1)."""
from __future__ import annotations

import numpy as np

from bar.leverage import bridges, curl_leverage, observability_certificate
from bar.qc import gls_network


def _rng_network(n=8, n_chords=6, seed=0):
    """Connected network: spanning tree over n nodes + n_chords chords; positive random se."""
    rng = np.random.default_rng(seed)
    edges = []
    for k in range(1, n):
        edges.append((str(k), str(rng.integers(0, k)), float(rng.normal()),
                      float(rng.uniform(0.15, 0.6))))
    seen = {frozenset((a, b)) for a, b, _, _ in edges}
    while len(edges) < (n - 1) + n_chords:
        a, b = int(rng.integers(0, n)), int(rng.integers(0, n))
        if a != b and frozenset((str(a), str(b))) not in seen:
            edges.append((str(a), str(b), float(rng.normal()), float(rng.uniform(0.15, 0.6))))
            seen.add(frozenset((str(a), str(b))))
    return edges


def test_conservation_law_two_paths_agree():
    # curl_leverage internally checks h_M (projector) == 1 - w*Omega (dual) to tol;
    # a mismatch raises.
    h = curl_leverage(_rng_network(seed=1))
    assert np.all(h >= -1e-9) and np.all(h <= 1 + 1e-9)


def test_sum_h_equals_dof():
    e = _rng_network(seed=2)
    assert abs(float(curl_leverage(e).sum()) - gls_network(e).dof) < 1e-9


def test_triangle_plus_bridge():
    # triangle A-B-C (all V=1): each edge h=1/3, dof=1, sum h=1; pendant bridge C-D: h=0.
    e = [("A", "B", 0.0, 1.0), ("B", "C", 0.0, 1.0), ("C", "A", 0.0, 1.0), ("C", "D", 0.0, 1.0)]
    h = curl_leverage(e)
    assert np.allclose(h[:3], 1.0 / 3.0, atol=1e-9)
    assert abs(float(h[3])) < 1e-9
    assert abs(float(h.sum()) - gls_network(e).dof) < 1e-9
    assert bridges(e) == {3}


def test_bridge_edges_have_zero_leverage():
    # a chorded component joined by one bridge to a pendant node
    e = [("0", "1", 0.0, 0.7), ("1", "2", 0.0, 0.5), ("2", "0", 0.0, 0.9),  # triangle
         ("2", "3", 0.0, 0.4), ("3", "4", 0.0, 0.6)]                        # 2 chained bridges
    h = curl_leverage(e)
    br = bridges(e)
    assert br == {3, 4}
    assert np.all(np.abs(h[list(br)]) < 1e-9)


def test_gradient_bias_is_annihilated():
    # adding a node-consistent (gradient) bias to every edge leaves the closure chi^2 unchanged
    e = _rng_network(seed=3)
    nodes = sorted({a for a, _, _, _ in e} | {b for _, b, _, _ in e}, key=str)
    rng = np.random.default_rng(9)
    phi = {nd: float(rng.normal()) for nd in nodes}
    e_biased = [(a, b, ddg + (phi[b] - phi[a]), se) for a, b, ddg, se in e]
    assert abs(gls_network(e).chi2 - gls_network(e_biased).chi2) < 1e-7


def test_leverage_orientation_and_relabel_invariant():
    e = _rng_network(seed=4)
    h = curl_leverage(e)
    flipped = [(b, a, -ddg, se) for a, b, ddg, se in e]      # orientation flip
    assert np.allclose(h, curl_leverage(flipped), atol=1e-9)
    relabel = [("N" + a, "N" + b, ddg, se) for a, b, ddg, se in e]  # node relabel
    assert np.allclose(np.sort(h), np.sort(curl_leverage(relabel)), atol=1e-9)


def test_observability_certificate_fields():
    e = [("A", "B", 0.0, 1.0), ("B", "C", 0.0, 1.0), ("C", "A", 0.0, 1.0), ("C", "D", 0.0, 1.0)]
    cert = observability_certificate(e)
    by_idx = {c["index"]: c for c in cert}
    assert by_idx[3]["is_bridge"] and not by_idx[3]["auditable"]
    assert np.isinf(by_idx[3]["delta_star"])
    assert by_idx[0]["auditable"] and by_idx[0]["delta_star"] < np.inf
    # conservation law field: h + w*Omega == 1
    assert all(abs(c["h"] + c["w_times_Omega"] - 1.0) < 1e-9 for c in cert)
