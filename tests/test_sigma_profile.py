"""Invariants for the rank-transferred learned-sigma profile (peer-review item P4)."""
from __future__ import annotations

import numpy as np

from bar.sigma_profile import PROFILE_POINTS, rank_transfer, shuffled


def test_output_length_matches_input():
    ov = np.array([0.01, 0.05, 0.10, 0.20, 0.23])
    assert rank_transfer(ov).shape == ov.shape


def test_ratios_stay_inside_the_measured_profile_range():
    ov = np.linspace(0.0001, 0.233, 50)
    r = rank_transfer(ov)
    lo = min(p[1] for p in PROFILE_POINTS)
    hi = max(p[1] for p in PROFILE_POINTS)
    assert r.min() >= lo - 1e-12 and r.max() <= hi + 1e-12


def test_rank_not_raw_value_drives_the_transfer():
    """Two overlap sets with the same ORDERING but wildly different scales must get the same
    ratios -- that is what makes the transfer a rank transfer rather than a value lookup."""
    a = np.array([0.001, 0.002, 0.003, 0.004])
    b = np.array([0.10, 0.20, 0.30, 0.40])
    assert np.allclose(rank_transfer(a), rank_transfer(b))


def test_endpoints_hit_the_profile_endpoints():
    ov = np.array([0.01, 0.05, 0.10, 0.23])
    r = rank_transfer(ov)
    pts = sorted(PROFILE_POINTS)
    assert abs(float(r[0]) - pts[0][1]) < 1e-12    # lowest-overlap edge -> lowest-overlap ratio
    assert abs(float(r[-1]) - pts[-1][1]) < 1e-12  # highest-overlap edge -> highest-overlap ratio


def test_all_equal_overlaps_do_not_crash_and_are_constant():
    r = rank_transfer(np.full(6, 0.12))
    assert np.all(np.isfinite(r)) and np.allclose(r, r[0])


def test_single_edge_is_finite():
    r = rank_transfer(np.array([0.12]))
    assert r.shape == (1,) and np.isfinite(r[0])


def test_shuffle_preserves_the_multiset_and_is_seeded():
    r = rank_transfer(np.linspace(0.01, 0.23, 20))
    s1 = shuffled(r, seed=20260808)
    s2 = shuffled(r, seed=20260808)
    assert np.allclose(np.sort(s1), np.sort(r))   # same multiset
    assert np.allclose(s1, s2)                    # deterministic
    assert not np.allclose(s1, r)                 # actually permuted
