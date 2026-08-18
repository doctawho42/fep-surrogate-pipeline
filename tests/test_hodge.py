"""Theorem 5: the Hodge split of a perturbation network's error field.

Every invariant of the theorem is an assertion here, and each is checked on a network whose
answer is known independently of the implementation under test.
"""
from __future__ import annotations

import hashlib
import pathlib

import numpy as np
import pytest

from bar.hodge import (
    gradient_field,
    gradient_r2,
    hodge_split,
    influence_rank,
    influence_repair_order,
    unidentifiable_dimension,
)
from bar.qc import _incidence, gls_network

ROOT = pathlib.Path(__file__).resolve().parent.parent
PREREG = ROOT / "data" / "openfe_replicates" / "hodge_prereg.yaml"
# Anchors the pre-registration against silent edits (same device as closeloop/prospective).
PREREG_SHA256 = "16279af73d90e1e9b39f33046ba6ee432b80ad750a6e5687a0e05b3e54b01585"


def triangle(se: float = 1.0):
    """3 nodes, 3 edges, one independent cycle: dof = 1, rank = 2."""
    return [("a", "b", 0.0, se), ("b", "c", 0.0, se), ("a", "c", 0.0, se)]


def two_triangles_and_a_tail():
    """A network with a bridge, so h_e = 0 is exercised."""
    return [("a", "b", 0.0, 1.0), ("b", "c", 0.0, 1.0), ("a", "c", 0.0, 1.0),
            ("c", "d", 0.0, 1.0),                                   # bridge
            ("d", "e", 0.0, 1.0), ("e", "f", 0.0, 1.0), ("d", "f", 0.0, 1.0)]


def two_components():
    return [("a", "b", 0.0, 1.0), ("b", "c", 0.0, 1.0), ("a", "c", 0.0, 1.0),
            ("x", "y", 0.0, 1.0), ("y", "z", 0.0, 1.0), ("x", "z", 0.0, 1.0)]


def test_prereg_is_immutable():
    got = hashlib.sha256(PREREG.read_bytes()).hexdigest()
    assert got == PREREG_SHA256, "hodge_prereg.yaml changed after freezing"


def test_gradient_field_has_no_cycle_part():
    """Theorem 5(i): a per-ligand bias is annihilated by the residual maker, exactly."""
    edges = two_triangles_and_a_tail()
    rng = np.random.default_rng(0)
    nodes, _, _, _ = _incidence(edges)
    mu = gradient_field(edges, dict(zip(nodes, rng.normal(size=len(nodes)), strict=True)))
    split = hodge_split(edges, mu)
    assert np.abs(split.cycle).max() < 1e-10
    assert split.visible_fraction < 1e-12


def test_cycle_field_has_no_gradient_part_and_the_split_is_orthogonal():
    edges = two_triangles_and_a_tail()
    rng = np.random.default_rng(1)
    v = rng.normal(size=len(edges))
    split = hodge_split(edges, v)
    # re-splitting the cycle part returns it unchanged, and its gradient part vanishes
    again = hodge_split(edges, split.cycle)
    assert np.abs(again.gradient).max() < 1e-10
    assert np.allclose(again.cycle, split.cycle, atol=1e-10)
    # Pythagoras in the whitened metric
    w = 1.0 / np.array([e[3] for e in edges]) ** 2
    tot = float(np.sum(w * v * v))
    assert np.sum(w * split.gradient ** 2) + np.sum(w * split.cycle ** 2) == pytest.approx(tot)


def test_visible_fraction_of_an_unaligned_field_is_dof_over_E():
    """Theorem 5(iii): the chance level, in the whitened metric."""
    edges = two_triangles_and_a_tail()
    fit = gls_network(edges)
    se = np.array([e[3] for e in edges])
    rng = np.random.default_rng(2)
    fr = [hodge_split(edges, rng.normal(size=len(edges)) * se).visible_fraction
          for _ in range(4000)]
    assert np.mean(fr) == pytest.approx(fit.dof / len(edges), abs=0.02)


def test_adjusted_r2_of_an_unaligned_field_is_zero_on_average():
    edges = two_triangles_and_a_tail()
    se = np.array([e[3] for e in edges])
    rng = np.random.default_rng(3)
    adj = [gradient_r2(edges, rng.normal(size=len(edges)) * se)[1] for _ in range(4000)]
    assert np.mean(adj) == pytest.approx(0.0, abs=0.03)


def test_adjusted_r2_is_one_for_a_pure_gradient():
    edges = two_triangles_and_a_tail()
    nodes, _, _, _ = _incidence(edges)
    mu = gradient_field(edges, dict(zip(nodes, [0.0, 1.0, 2.0, 3.0, 4.0, 5.0], strict=True)))
    r2, adj, chance = gradient_r2(edges, mu)
    assert r2 == pytest.approx(1.0)
    assert adj == pytest.approx(1.0)
    assert chance == pytest.approx(gls_network(edges).dof / len(edges))


def test_cycle_part_induces_no_bias_and_gradient_part_no_detectability():
    """Theorem 5(i'): the two halves are exact orthogonal complements.

    The bias a systematic edge error induces in the fitted potentials is computed by re-fitting a
    noiseless network carrying that error, so this checks the claim through the shipped GLS path
    rather than through a re-derivation of the same algebra.
    """
    edges = two_triangles_and_a_tail()
    rng = np.random.default_rng(4)
    split = hodge_split(edges, rng.normal(size=len(edges)))
    base = gls_network(edges).potentials

    def refit_with(mu):
        return gls_network([(a, b, mu[i], se) for i, (a, b, _, se) in enumerate(edges)]).potentials

    # the cycle part shifts no potential ...
    moved = refit_with(split.cycle) - base
    assert np.abs(moved - moved.mean()).max() < 1e-10
    # ... and the gradient part contributes nothing to the closure statistic
    chi2 = gls_network([(a, b, split.gradient[i], se)
                        for i, (a, b, _, se) in enumerate(edges)]).chi2
    assert chi2 < 1e-16


def test_influence_rank_is_zero_on_bridges_and_grows_as_leverage_falls():
    edges = two_triangles_and_a_tail()
    infl = influence_rank(edges, z=np.ones(len(edges)))
    assert infl[3] == 0.0, "a bridge carries no auditable evidence, so it cannot be ranked"
    # at equal |z|, the sqrt((1-h)/h) factor orders edges by how much bias their error induces
    from bar.leverage import curl_leverage
    h = curl_leverage(edges)
    inner = [i for i in range(len(edges)) if h[i] > 1e-9]
    assert np.corrcoef(h[inner], infl[inner])[0, 1] < 0


def inconsistent_network():
    """The same shape, but the cycles do not close, so there is something to repair."""
    return [("a", "b", 1.0, 1.0), ("b", "c", 1.0, 1.0), ("a", "c", 0.0, 1.0),
            ("c", "d", 0.0, 1.0),                                   # bridge
            ("d", "e", 1.0, 1.0), ("e", "f", 1.0, 1.0), ("d", "f", 0.0, 1.0)]


def test_influence_repair_order_respects_k_and_skips_bridges():
    edges = inconsistent_network()
    order = influence_repair_order(edges, k=2)
    assert len(order) == 2
    assert 3 not in order


def test_influence_repair_order_stops_when_there_is_nothing_to_repair():
    """A network whose cycles close exactly offers no edge to act on, so the order is empty."""
    assert influence_repair_order(two_triangles_and_a_tail(), k=3) == []


def test_unidentifiable_dimension_counts_components():
    """Theorem 5(iv): the first anchor in each component buys nothing."""
    tri = triangle()
    assert unidentifiable_dimension(tri, anchored=set()) == 2          # N - c = 3 - 1
    assert unidentifiable_dimension(tri, anchored={"a"}) == 2          # first anchor: no gain
    assert unidentifiable_dimension(tri, anchored={"a", "b"}) == 1
    assert unidentifiable_dimension(tri, anchored={"a", "b", "c"}) == 0

    two = two_components()
    assert unidentifiable_dimension(two, anchored=set()) == 4          # N - c = 6 - 2
    assert unidentifiable_dimension(two, anchored={"a"}) == 4          # only one component pinned
    assert unidentifiable_dimension(two, anchored={"a", "x"}) == 4     # both offsets pinned
    assert unidentifiable_dimension(two, anchored={"a", "b", "x"}) == 3


def test_unidentifiable_dimension_matches_a_brute_force_rank():
    edges = two_triangles_and_a_tail()
    nodes, B, _, _ = _incidence(edges)
    for anchored in [set(), {"a"}, {"a", "d"}, {"a", "b", "d"}, set(nodes)]:
        keep = [i for i, n in enumerate(nodes) if n not in anchored]
        brute = int(np.linalg.matrix_rank(B[:, keep])) if keep else 0
        assert unidentifiable_dimension(edges, anchored=anchored) == brute
