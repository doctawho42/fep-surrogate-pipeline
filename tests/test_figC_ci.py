"""TDD test for bootstrap_ratio_ci (reviewer S5.7 — edge-weighting row CI)."""
import numpy as np
import pytest


def test_bootstrap_ratio_ci_ordered_finite():
    """CI helper returns a finite, ordered pair from a seeded ratio array."""
    from src.bar.active import bootstrap_ratio_ci  # noqa: PLC0415

    rng = np.random.default_rng(0)
    ratios = rng.normal(1.0, 0.03, 40)
    lo, hi = bootstrap_ratio_ci(ratios, n_boot=500, seed=1)
    assert np.isfinite(lo) and np.isfinite(hi) and lo <= hi


def test_bootstrap_ratio_ci_covers_mean():
    """At large n the 95% CI should contain the population mean."""
    from src.bar.active import bootstrap_ratio_ci  # noqa: PLC0415

    rng = np.random.default_rng(42)
    ratios = rng.normal(1.05, 0.10, 200)
    lo, hi = bootstrap_ratio_ci(ratios, n_boot=2000, seed=0)
    # With n=200 the CI should contain 1.05 with high probability
    assert lo <= 1.10 and hi >= 1.00


def test_bootstrap_ratio_ci_deterministic():
    """Same seed produces same result."""
    from src.bar.active import bootstrap_ratio_ci  # noqa: PLC0415

    rng = np.random.default_rng(7)
    ratios = rng.normal(1.0, 0.05, 30)
    r1 = bootstrap_ratio_ci(ratios, n_boot=300, seed=99)
    r2 = bootstrap_ratio_ci(ratios, n_boot=300, seed=99)
    assert r1 == r2
