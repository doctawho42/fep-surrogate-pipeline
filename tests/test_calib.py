from __future__ import annotations

import numpy as np

from bar.calib import conformal_q, coverage, ensemble_predict, train_ensemble


def test_ensemble_predicts_and_spreads():
    rng = np.random.default_rng(0)
    X = rng.normal(size=(200, 4))
    y = X[:, 0] * 2.0 - X[:, 1] + rng.normal(0, 0.1, 200)
    nets = train_ensemble(X, y, n_members=5, epochs=120, seed0=0)
    mu, sd = ensemble_predict(nets, X)
    assert mu.shape == (200,) and sd.shape == (200,)
    assert np.all(sd >= 0)
    assert sd.mean() > 1e-3  # members actually disagree (real epistemic spread)
    # in-sample fit is decent
    assert np.corrcoef(mu, y)[0, 1] > 0.9


def test_normalized_conformal_hits_marginal_coverage():
    rng = np.random.default_rng(1)
    n = 4000
    mu = np.zeros(n)
    sigma = rng.uniform(0.5, 2.0, n)            # heteroscedastic
    y = rng.normal(mu, sigma)                   # true noise = sigma
    half = n // 2
    q = conformal_q(np.abs(y[:half] - mu[:half]), sigma[:half], alpha=0.10)
    cov = coverage(y[half:], mu[half:], sigma[half:], q)
    assert 0.86 <= cov <= 0.94                  # ~90% marginal coverage
