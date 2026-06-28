"""BAR estimator invariants (Theorems 1 & 2). All numbers cross-checked against
docs/bar_proofs.tex and reproduced in scratchpad/verify_proofs.py.

Correctness invariants exercised here:
  #1 variance = sandwich B/I^2 (never 1/I); matches pymbar to machine precision
  #2 calibration measured against the sandwich
  (T1) information-share gradient + curvature = -I
"""
from __future__ import annotations

import numpy as np
import pytest

from bar.estimator import (
    bar_estimate,
    bar_score,
    bennett_variance,
    default_log_count_ratio,
    fisher_info,
    information_shares,
    mbar_variance,
    naive_variance,
    sandwich_variance,
    solve_bar,
)

# --- proof-sheet instance: x_f={0,1}, x_r={-0.5}, M=0 ----------------------
XF0 = np.array([0.0, 1.0])
XR0 = np.array([-0.5])


def gaussian_works(s: float, n: int, rng: np.random.Generator):
    """Correctly-specified Gaussian work model: forward x~N(mf,s^2), reverse
    x~N(mr,s^2) with mf-mr=s^2, dF=(mf+mr)/2=0. Separation in sigma = s.
    Returns unified coordinates (n_f=n_r=n so M=0)."""
    mf, mr = s**2 / 2, -(s**2) / 2
    return rng.normal(mf, s, n), rng.normal(mr, s, n)


# --------------------------------------------------------------------------
# T1 — root, score, curvature, information share
# --------------------------------------------------------------------------
def test_solve_bar_known_instance():
    assert solve_bar(XF0, XR0, M=0.0) == pytest.approx(-0.59571, abs=1e-4)


def test_default_log_count_ratio():
    assert default_log_count_ratio(20, 10) == pytest.approx(np.log(2.0))
    assert default_log_count_ratio(7, 7) == 0.0


def test_score_is_zero_at_root():
    d = solve_bar(XF0, XR0, M=0.0)
    assert bar_score(d, XF0, XR0, M=0.0) == pytest.approx(0.0, abs=1e-10)


def test_information_share_closed_form():
    d = solve_bar(XF0, XR0, M=0.0)
    If, _Ir, I = fisher_info(d, XF0, XR0, M=0.0)
    assert If / I == pytest.approx(0.596825, abs=1e-5)


def test_information_shares_sum_to_one():
    share_f, share_r = information_shares(XF0, XR0, M=0.0)
    assert abs(share_f) + abs(share_r) == pytest.approx(1.0, abs=1e-12)


def test_curvature_equals_minus_I():
    d = solve_bar(XF0, XR0, M=0.0)
    eps = 1e-6
    dS = (bar_score(d + eps, XF0, XR0, 0.0) - bar_score(d - eps, XF0, XR0, 0.0)) / (2 * eps)
    _, _, I = fisher_info(d, XF0, XR0, M=0.0)
    assert dS == pytest.approx(-I, abs=1e-5)


# --------------------------------------------------------------------------
# #1 — variance forms: sandwich, pymbar match, never 1/I
# --------------------------------------------------------------------------
def test_naive_variance_is_one_over_I():
    d = solve_bar(XF0, XR0, M=0.0)
    _, _, I = fisher_info(d, XF0, XR0, M=0.0)
    assert naive_variance(XF0, XR0, M=0.0) == pytest.approx(1.0 / I, rel=1e-12)


def test_mbar_variance_matches_pymbar_to_machine_precision():
    """Invariant #1: our closed form 1/I-(1/nf+1/nr) IS pymbar's 'MBAR' uncertainty."""
    pmbar = pytest.importorskip("pymbar.other_estimators").bar
    rng = np.random.default_rng(0)
    for s in (1.0, 1.7, 2.4):
        xf, xr = gaussian_works(s, 40, rng)
        ours = mbar_variance(xf, xr)
        # mapping: unified x_i=W^f_i, x_j=-W^r_j  =>  bar(w_F=xf, w_R=-xr)
        pm = pmbar(xf, -xr, uncertainty_method="MBAR")["dDelta_f"] ** 2
        assert ours == pytest.approx(pm, rel=1e-9)


def test_bennett_variance_matches_pymbar_to_machine_precision():
    pmbar = pytest.importorskip("pymbar.other_estimators").bar
    rng = np.random.default_rng(1)
    for s in (1.0, 1.7, 2.4):
        xf, xr = gaussian_works(s, 40, rng)
        ours = bennett_variance(xf, xr)
        pm = pmbar(xf, -xr, uncertainty_method="BAR")["dDelta_f"] ** 2
        assert ours == pytest.approx(pm, rel=1e-9)


def test_sandwich_matches_mbar_asymptotically():
    """B/I^2 and 1/I-nrat are asymptotically equal IN EXPECTATION (verified ~1% at
    N>=50). On a single draw they differ by plug-in noise, so average over reps."""
    rng = np.random.default_rng(7)
    sand, mbar = [], []
    for _ in range(300):
        xf, xr = gaussian_works(1.5, 400, rng)
        sand.append(sandwich_variance(xf, xr))
        mbar.append(mbar_variance(xf, xr))
    assert np.mean(sand) == pytest.approx(np.mean(mbar), rel=0.03)


# --------------------------------------------------------------------------
# #2 — calibration: sandwich tracks MC truth; naive 1/I does NOT (varying factor)
# --------------------------------------------------------------------------
def _mc_ratios(s: float, n: int, reps: int, seed: int):
    rng = np.random.default_rng(seed)
    dhats, sand_se, naive_se = [], [], []
    for _ in range(reps):
        xf, xr = gaussian_works(s, n, rng)
        res = bar_estimate(xf, xr)
        dhats.append(res.delta_f)
        sand_se.append(np.sqrt(res.var_sandwich))
        naive_se.append(np.sqrt(res.var_naive))
    emp = float(np.std(dhats, ddof=1))
    return emp, np.mean(sand_se) / emp, np.mean(naive_se) / emp


def test_sandwich_calibrated_at_high_overlap():
    # n=20, >=2000 reps, fixed seed -> deterministic, no flakiness
    emp, sand_ratio, _ = _mc_ratios(s=1.0, n=20, reps=2000, seed=12345)
    assert emp == pytest.approx(0.160, abs=0.02)
    assert 0.95 <= sand_ratio <= 1.05


def test_naive_overestimates_at_high_overlap():
    _, _, naive_ratio = _mc_ratios(s=1.0, n=20, reps=2000, seed=12345)
    assert naive_ratio == pytest.approx(2.2, abs=0.25)  # ~2.2x too large in se


def test_naive_error_is_a_varying_factor():
    """The Fig A thesis: naive/emp shrinks 2.3 -> 1.2 across overlap; no constant
    rescale fixes 1/I. Sandwich stays ~1 everywhere."""
    _, s1_sand, s1_naive = _mc_ratios(s=1.0, n=20, reps=2000, seed=777)
    _, s2_sand, s2_naive = _mc_ratios(s=2.4, n=20, reps=2000, seed=778)
    assert s1_naive > 1.9              # poor overlap: badly over-conservative
    assert s2_naive < 1.45             # good overlap: closer to 1
    assert s1_naive / s2_naive > 1.4   # the factor genuinely varies
    assert 0.93 <= s1_sand <= 1.07     # sandwich calibrated at both ends
    assert 0.93 <= s2_sand <= 1.07


def test_bar_result_fields_consistent():
    res = bar_estimate(XF0, XR0, M=0.0)
    assert res.delta_f == pytest.approx(solve_bar(XF0, XR0, M=0.0))
    assert res.var_sandwich > 0
    assert res.overlap > 0
    assert abs(res.info_share_f) + abs(res.info_share_r) == pytest.approx(1.0)
