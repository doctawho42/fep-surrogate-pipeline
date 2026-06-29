"""Risk-adjusted reward from the calibrated BAR estimate (target-contour Increment 1).

The reward turns the BAR-bottleneck's calibrated SANDWICH uncertainty into a
decision: a higher-is-better value (``-ΔΔĜ``) penalised by its standard error
(a lower-confidence bound, ``value - κ·σ``). Spec:
docs/superpowers/specs/2026-06-30-target-contour-design.md.
"""
from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike

from bar.estimator import bar_estimate


def risk_adjusted_reward(value: float, sigma: float, kappa: float = 1.0) -> float:
    """Lower-confidence-bound reward on a higher-is-better ``value``:
    ``r = value - kappa*sigma``. At ``sigma=0`` it is the raw value; it decreases
    monotonically in ``sigma`` (``kappa>=0``) and in ``kappa`` (``sigma>=0``)."""
    return float(value - kappa * sigma)


def edge_reward(x_f: ArrayLike, x_r: ArrayLike, kappa: float = 1.0,
                sigma_se: float | None = None) -> float:
    """Risk-adjusted reward for a BAR edge. value = ``-delta_f`` (more negative
    ΔΔG = stronger binder = higher reward); ``sigma`` = sandwich se unless an
    overriding standard error ``sigma_se`` (e.g. a learned-variance head) is given."""
    r = bar_estimate(x_f, x_r)
    sigma = float(np.sqrt(max(r.var_sandwich, 0.0))) if sigma_se is None else float(sigma_se)
    return risk_adjusted_reward(-r.delta_f, sigma, kappa)
