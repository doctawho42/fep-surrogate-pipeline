"""Invariants for the comparative cycle-closure detectors (peer-review item P1)."""
from __future__ import annotations

import numpy as np

from bar.detectors import (
    anchor_score,
    auc_flag_vs_anchor,
    flag_calibrated,
    flag_fixed_cutoff,
    flag_fixed_se,
    fundamental_cycle_sums,
    paired_auc_bootstrap,
)
from bar.qc import gls_network


def _triangle(imbalance: float, se: float = 0.1):
    """A-B-C triangle whose signed cycle sum is exactly `imbalance`."""
    return [("A", "B", 1.0, se), ("B", "C", 1.0, se), ("A", "C", 2.0 - imbalance, se)]


def test_cycle_sum_matches_injected_imbalance():
    # cycle A->B->C then back C->A: (+1) + (+1) - (2 - imb) = imb
    sums = fundamental_cycle_sums(_triangle(0.7))
    assert len(sums) == 1
    assert abs(abs(sums[0]) - 0.7) < 1e-12


def test_cycle_count_equals_dof():
    edges = _triangle(0.0) + [("C", "D", 0.5, 0.1), ("A", "D", 1.5, 0.1)]
    assert len(fundamental_cycle_sums(edges)) == gls_network(edges).dof


def test_bridge_only_network_has_no_cycles():
    edges = [("A", "B", 1.0, 0.1), ("B", "C", 1.0, 0.1)]
    assert fundamental_cycle_sums(edges) == []


def test_fixed_cutoff_flags_only_above_threshold():
    systems = {"big": _triangle(1.5), "small": _triangle(0.2)}
    flags = flag_fixed_cutoff(systems, cutoff=1.0)
    assert flags["big"] and not flags["small"]


def test_fixed_se_equals_calibrated_when_all_se_equal():
    # with a constant per-edge se the pooled-se detector IS the calibrated one
    systems = {"s1": _triangle(1.2), "s2": _triangle(0.05), "s3": _triangle(0.6)}
    assert flag_fixed_se(systems) == flag_calibrated(systems)


def test_calibrated_flags_the_imbalanced_system():
    # tight bars + a large imbalance must be flagged; a clean system must not
    systems = {"bad": _triangle(2.0, se=0.05), "ok": _triangle(0.0, se=0.05)}
    flags = flag_calibrated(systems)
    assert flags["bad"] and not flags["ok"]


def test_anchor_is_high_for_reproducible_systematic_error():
    """The anchor is |mean_k z| per edge: a systematic offset repeated across replicates
    scores high; independent noise of the same size averages down."""
    rng = np.random.default_rng(0)
    syst = [[("A", "B", 1.0, 0.1), ("B", "C", 1.0, 0.1), ("A", "C", 3.0, 0.1)] for _ in range(3)]
    noisy = [
        [("A", "B", 1.0, 0.1), ("B", "C", 1.0, 0.1),
         ("A", "C", 2.0 + 1.0 * rng.standard_normal(), 0.1)]
        for _ in range(3)
    ]
    assert anchor_score(syst) > anchor_score(noisy)


def test_auc_is_one_when_flags_perfectly_order_the_anchor():
    flags = np.array([True, True, False, False])
    anchor = np.array([9.0, 8.0, 2.0, 1.0])
    assert auc_flag_vs_anchor(flags, anchor) == 1.0


def test_auc_is_half_when_a_group_is_empty():
    assert auc_flag_vs_anchor(np.array([False, False]), np.array([1.0, 2.0])) == 0.5


def test_paired_bootstrap_zero_difference_for_identical_flags():
    f = np.array([True, False, True, False, True, False])
    anchor = np.array([5.0, 1.0, 6.0, 2.0, 7.0, 0.5])
    out = paired_auc_bootstrap(f, f, f, anchor, n_boot=200, seed=1)
    assert abs(out["diff"]) < 1e-12
    assert out["verdict"] == "TIE"
    assert set(out) >= {"auc_a", "auc_b", "auc_c", "diff", "ci_lo", "ci_hi", "verdict"}
