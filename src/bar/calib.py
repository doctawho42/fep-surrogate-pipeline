"""Reusable uncertainty machinery for the amortized reward: a bootstrap deep ensemble
(epistemic spread) and normalized/Mondrian split-conformal (marginal-coverage guarantee
with per-edge width). Generalises the 1-D demonstration in figs/make_figB.py to an
arbitrary feature dimension.
"""
from __future__ import annotations

import numpy as np
import torch
from numpy.typing import NDArray


def _make_net(in_dim: int, seed: int) -> torch.nn.Module:
    torch.manual_seed(seed)
    return torch.nn.Sequential(
        torch.nn.Linear(in_dim, 64), torch.nn.SiLU(),
        torch.nn.Linear(64, 64), torch.nn.SiLU(), torch.nn.Linear(64, 1),
    )


def train_ensemble(X: NDArray, y: NDArray, n_members: int = 8, epochs: int = 300,
                   seed0: int = 0) -> list[tuple[torch.nn.Module, NDArray, NDArray]]:
    """Bootstrap deep ensemble of small MLPs (epistemic uncertainty = member spread)."""
    Xa = np.asarray(X, dtype=np.float32)
    ya = np.asarray(y, dtype=np.float32)
    x_mean, x_std = Xa.mean(0), Xa.std(0) + 1e-6
    nets: list[tuple[torch.nn.Module, NDArray, NDArray]] = []
    for m in range(n_members):
        seed = seed0 + m
        net = _make_net(Xa.shape[1], seed)
        opt = torch.optim.Adam(net.parameters(), lr=5e-3, weight_decay=1e-4)
        rng = np.random.default_rng(seed)
        idx = rng.integers(0, len(Xa), len(Xa))  # bootstrap resample
        xt = torch.tensor((Xa[idx] - x_mean) / x_std)
        yt = torch.tensor(ya[idx, None])
        lossf = torch.nn.MSELoss()
        for _ in range(epochs):
            opt.zero_grad()
            loss = lossf(net(xt), yt)
            loss.backward()
            opt.step()
        nets.append((net, x_mean, x_std))
    return nets


def ensemble_predict(
    nets: list[tuple[torch.nn.Module, NDArray, NDArray]], X: NDArray
) -> tuple[NDArray, NDArray]:
    """Return (mean, std) over ensemble members for inputs X."""
    Xa = np.asarray(X, dtype=np.float32)
    preds: list[NDArray] = []
    for net, x_mean, x_std in nets:
        xt = torch.tensor((Xa - x_mean) / x_std)
        with torch.no_grad():
            preds.append(net(xt).numpy().ravel())
    P = np.stack(preds)
    return P.mean(0), P.std(0, ddof=1)


def conformal_q(residuals: NDArray, sigma: NDArray, alpha: float) -> float:
    """Normalized split-conformal quantile q s.t. |y-mu| <= q*sigma holds at marginal
    (1-alpha) coverage. residuals = |y-mu| on a calibration fold; sigma = sigma_total there."""
    r = np.asarray(residuals, dtype=float)
    s = np.asarray(sigma, dtype=float)
    n = r.size
    k = min(np.ceil((n + 1) * (1.0 - alpha)) / n, 1.0)
    return float(np.quantile(r / np.maximum(s, 1e-9), k))


def coverage(y: NDArray, mu: NDArray, sigma: NDArray, q: float) -> float:
    """Fraction of points within mu +/- q*sigma."""
    return float(np.mean(np.abs(np.asarray(y) - np.asarray(mu)) <= q * np.asarray(sigma)))


# --- Recalibration baselines (referee foils) ------------------------------------------
# Each fits on a held-out calibration fold and returns a *new sigma for the SAME mean*.
# They are the standard, cheap recalibrators a practitioner would apply to ANY uncertainty
# head (incl. an overconfident MVE), so the honest question is whether the physics-derived
# decomposed sigma beats a properly-recalibrated learned sigma, not just a raw one.

def temperature_scale(residuals: NDArray, sigma_cal: NDArray) -> float:
    """Single multiplicative factor T (variance-temperature) minimizing Gaussian NLL on the
    cal fold: closed form T = sqrt(mean((res/sigma)^2)). Recalibrated sigma_test = T*sigma_test."""
    r = np.asarray(residuals, dtype=float)
    s = np.maximum(np.asarray(sigma_cal, dtype=float), 1e-9)
    return float(np.sqrt(np.mean((r / s) ** 2)))


def platt_scale(residuals: NDArray, sigma_cal: NDArray) -> tuple[float, float]:
    """Affine-in-log-sigma recalibration: fit log(sigma') = a*log(sigma)+b by Gaussian NLL.
    Returns (a, b); apply as sigma' = exp(b) * sigma**a. Generalizes temperature (a=1)."""
    from scipy.optimize import minimize  # local import; scipy is a core dep

    r = np.asarray(residuals, dtype=float)
    ls = np.log(np.maximum(np.asarray(sigma_cal, dtype=float), 1e-9))
    r2 = r ** 2

    def nll(theta: NDArray) -> float:
        a, b = float(theta[0]), float(theta[1])
        logv = 2.0 * (a * ls + b)  # log(sigma'^2)
        return float(np.mean(0.5 * logv + 0.5 * r2 / np.exp(logv)))

    res = minimize(nll, x0=np.array([1.0, 0.0]), method="Nelder-Mead")
    a, b = float(res.x[0]), float(res.x[1])
    return a, b


def isotonic_sigma(residuals: NDArray, sigma_cal: NDArray):
    """Monotone (isotonic) map from base sigma to empirical error magnitude, fit on the cal
    fold. Returns a callable sigma_test -> recalibrated sigma. Non-parametric analogue of
    temperature/Platt (lets the sigma->error relation bend), using sklearn (already a dep)."""
    from sklearn.isotonic import IsotonicRegression

    s = np.asarray(sigma_cal, dtype=float)
    r = np.abs(np.asarray(residuals, dtype=float))
    iso = IsotonicRegression(out_of_bounds="clip", y_min=1e-6).fit(s, r)

    def apply(sigma_test: NDArray) -> NDArray:
        return np.maximum(iso.predict(np.asarray(sigma_test, dtype=float)), 1e-6)

    return apply
