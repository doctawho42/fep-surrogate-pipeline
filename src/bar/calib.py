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
    mu, sd = Xa.mean(0), Xa.std(0) + 1e-6
    nets: list[tuple[torch.nn.Module, NDArray, NDArray]] = []
    for m in range(n_members):
        seed = seed0 + m
        net = _make_net(Xa.shape[1], seed)
        opt = torch.optim.Adam(net.parameters(), lr=5e-3, weight_decay=1e-4)
        rng = np.random.default_rng(seed)
        idx = rng.integers(0, len(Xa), len(Xa))  # bootstrap resample
        xt = torch.tensor((Xa[idx] - mu) / sd)
        yt = torch.tensor(ya[idx, None])
        lossf = torch.nn.MSELoss()
        for _ in range(epochs):
            opt.zero_grad()
            loss = lossf(net(xt), yt)
            loss.backward()
            opt.step()
        nets.append((net, mu, sd))
    return nets


def ensemble_predict(
    nets: list[tuple[torch.nn.Module, NDArray, NDArray]], X: NDArray
) -> tuple[NDArray, NDArray]:
    """Return (mean, std) over ensemble members for inputs X."""
    Xa = np.asarray(X, dtype=np.float32)
    preds = []
    for net, mu, sd in nets:
        xt = torch.tensor((Xa - mu) / sd)
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
