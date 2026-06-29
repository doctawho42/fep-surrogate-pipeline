from __future__ import annotations

import numpy as np

from bar.reward import edge_reward, risk_adjusted_reward


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
