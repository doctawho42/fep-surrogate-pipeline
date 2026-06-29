from __future__ import annotations

import numpy as np

from bar.reward import (
    edge_reward,
    realized_hitrate,
    regret,
    regret_difference_ci,
    risk_adjusted_reward,
    select_topk,
)


def test_sigma_zero_collapses_to_raw_value():
    # σ=0 ⇒ LCB penalty vanishes ⇒ reward == raw value (machine precision)
    assert risk_adjusted_reward(1.234, 0.0, kappa=5.0) == 1.234


def test_reward_strictly_decreases_in_sigma():
    base = risk_adjusted_reward(0.0, 0.0, kappa=1.0)
    more = risk_adjusted_reward(0.0, 0.5, kappa=1.0)
    assert more < base


def test_reward_strictly_decreases_in_kappa():
    a = risk_adjusted_reward(0.0, 0.3, kappa=1.0)
    b = risk_adjusted_reward(0.0, 0.3, kappa=2.0)
    assert b < a


def test_edge_reward_uses_sandwich_and_negates_delta():
    rng = np.random.default_rng(0)
    s = 1.5
    xf = rng.normal(s**2 / 2, s, 200)
    xr = rng.normal(-(s**2) / 2, s, 200)
    # with kappa=0 the reward is exactly -delta_f
    from bar.estimator import bar_estimate
    r = bar_estimate(xf, xr)
    assert abs(edge_reward(xf, xr, kappa=0.0) - (-r.delta_f)) < 1e-12
    # a positive kappa lowers it (sandwich σ > 0 at finite overlap)
    assert edge_reward(xf, xr, kappa=1.0) < edge_reward(xf, xr, kappa=0.0)


def test_sigma_se_overrides_sandwich():
    rng = np.random.default_rng(1)
    xf = rng.normal(0.0, 1.0, 100)
    xr = rng.normal(0.0, 1.0, 100)
    base = edge_reward(xf, xr, kappa=1.0)
    overridden = edge_reward(xf, xr, kappa=1.0, sigma_se=999.0)
    assert overridden < base  # a large override σ makes the reward much smaller
    # negative override is clamped to 0 -> equals the kappa=0 (no-penalty) reward
    clamped = edge_reward(xf, xr, kappa=1.0, sigma_se=-5.0)
    assert abs(clamped - edge_reward(xf, xr, kappa=0.0)) < 1e-12


def test_select_topk_picks_highest_rewards_in_order():
    rewards = [0.1, 0.9, 0.3, 0.7]
    sel = select_topk(rewards, 2)
    assert list(sel) == [1, 3]


def test_hitrate_is_fraction_above_threshold():
    truth = np.array([0.0, 1.0, 2.0, 3.0])
    assert realized_hitrate(truth, [2, 3], threshold=1.5) == 1.0
    assert realized_hitrate(truth, [0, 1], threshold=1.5) == 0.0
    assert realized_hitrate(truth, [1, 2], threshold=1.5) == 0.5
    assert realized_hitrate(truth, [], threshold=1.5) == 0.0


def test_regret_zero_for_true_topk_and_positive_otherwise():
    truth = np.array([5.0, 1.0, 4.0, 2.0])
    assert regret(truth, [0, 2]) == 0.0           # true top-2 = {5,4}
    assert regret(truth, [1, 3]) > 0.0            # worst-2 selected


def test_sigma_zero_makes_calibrated_selection_equal_raw():
    # with σ=0 the calibrated reward == raw value ⇒ identical selection
    values = np.array([0.2, 0.5, 0.1, 0.9])
    raw = values
    cal = np.array([risk_adjusted_reward(v, 0.0, kappa=3.0) for v in values])
    assert list(select_topk(raw, 2)) == list(select_topk(cal, 2))


def test_regret_difference_ci_detects_calibrated_win():
    # calibrated regret consistently below raw ⇒ upper CI bound < 0
    cal = np.full(8, 0.10)
    raw = np.full(8, 0.40)
    mean_diff, lo, hi = regret_difference_ci(cal, raw)
    assert mean_diff < 0
    assert hi < 0  # KILL verdict: calibrated strictly better


def test_regret_difference_ci_null_when_equal():
    same = np.array([0.2, 0.25, 0.18, 0.22, 0.2, 0.21, 0.19, 0.23])
    mean_diff, lo, hi = regret_difference_ci(same, same)
    assert mean_diff == 0.0
    assert lo <= 0.0 <= hi
