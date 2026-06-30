"""Risk-adjusted reward from the calibrated BAR estimate (target-contour Increment 1).

The reward turns the BAR-bottleneck's calibrated SANDWICH uncertainty into a
decision: a higher-is-better value (``-ΔΔĜ``) penalised by its standard error
(a lower-confidence bound, ``value - κ·σ``). Spec:
docs/superpowers/specs/2026-06-30-target-contour-design.md.
"""
from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray

from bar.chiral import chiral_readout
from bar.estimator import bar_estimate


def risk_adjusted_reward(value: float, sigma: float, kappa: float = 1.0) -> float:
    """Lower-confidence-bound reward on a higher-is-better ``value``:
    ``r = value - kappa*sigma``. At ``sigma=0`` it is the raw value; it decreases
    monotonically in ``sigma`` (``kappa>=0``) and in ``kappa`` (``sigma>=0``)."""
    return float(value - kappa * sigma)


def edge_reward(x_f: ArrayLike, x_r: ArrayLike, kappa: float = 1.0,
                sigma_se: float | None = None) -> float:
    """Risk-adjusted reward for a BAR edge. value = ``-delta_f`` (more negative
    ΔΔG = stronger binder = higher reward); sigma = sqrt(max(var_sandwich,0)) unless
    an overriding standard error ``sigma_se`` (clamped to >= 0) is given."""
    r = bar_estimate(x_f, x_r)
    if sigma_se is None:
        sigma = float(np.sqrt(max(r.var_sandwich, 0.0)))
    else:
        sigma = float(max(sigma_se, 0.0))
    return risk_adjusted_reward(-r.delta_f, sigma, kappa)


def select_topk(rewards: ArrayLike, k: int) -> NDArray:
    """Indices of the top-``k`` by reward (descending; ties broken by index)."""
    r = np.asarray(rewards, dtype=float)
    return np.argsort(-r, kind="stable")[:k]


def realized_hitrate(truth: ArrayLike, selected: ArrayLike, threshold: float) -> float:
    """Fraction of the ``selected`` whose true value is >= ``threshold`` (0 if none)."""
    t = np.asarray(truth, dtype=float)
    sel = np.asarray(selected, dtype=int)
    if sel.size == 0:
        return 0.0
    return float(np.mean(t[sel] >= threshold))


def regret(truth: ArrayLike, selected: ArrayLike) -> float:
    """Simple regret = mean(true top-k) - mean(true value of selected k). >= 0; 0 iff
    the selected set is a true top-k set."""
    t = np.asarray(truth, dtype=float)
    sel = np.asarray(selected, dtype=int)
    k = sel.size
    if k == 0:
        return 0.0
    best = float(np.sort(t)[::-1][:k].mean())
    got = float(t[sel].mean())
    return best - got


def regret_difference_ci(regret_a: ArrayLike, regret_b: ArrayLike, alpha: float = 0.05,
                         n_boot: int = 2000, seed: int = 0) -> tuple[float, float, float]:
    """Bootstrap CI on the paired regret difference ``a - b`` across seeds. Returns
    ``(mean_diff, lo, hi)``. For the gate (a=calibrated, b=raw), calibrated beats raw
    iff ``hi < 0`` (its regret is strictly lower). A NEGATIVE ``mean_diff``/``hi``
    means ``a`` has LOWER regret than ``b`` (i.e. a is better). Generic paired-bootstrap
    CI on per-seed difference a-b; callers may pass any paired metric (regret: lower
    better -> winner iff hi<0; precision/coverage: higher better -> winner iff lo>0)."""
    a = np.asarray(regret_a, dtype=float)
    b = np.asarray(regret_b, dtype=float)
    if a.shape != b.shape:
        raise ValueError(f"regret arrays must match shape: {a.shape} != {b.shape}")
    d = a - b
    rng = np.random.default_rng(seed)
    boots = np.array([d[rng.integers(0, d.size, d.size)].mean() for _ in range(n_boot)])
    lo, hi = np.percentile(boots, [100 * alpha / 2, 100 * (1 - alpha / 2)])
    return float(d.mean()), float(lo), float(hi)


def linear_readout_reward(coords: ArrayLike, weights: ArrayLike,
                          include_0o: bool = True) -> float:
    """Reward = ``weights · chiral_readout(coords, include_0o)``. The representation
    contract for the future amortised reward: WITH the parity-odd ``0o`` channel the
    reward separates enantiomers (mirror images); WITHOUT it the readout is
    O(3)-invariant and the reward is identical for an enantiomer pair (chirality-blind,
    Thm 4). ``weights`` must match the readout length (6 even, or 7 with ``0o``)."""
    feats = chiral_readout(coords, include_0o=include_0o)
    w = np.asarray(weights, dtype=float)
    if w.ndim != 1:
        raise ValueError(f"weights must be 1-D, got shape {w.shape}")
    if w.shape[0] != feats.shape[0]:
        raise ValueError(f"weights length {w.shape[0]} != readout length {feats.shape[0]}")
    return float(w @ feats)
