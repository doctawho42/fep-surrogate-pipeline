"""BAR estimator core (Theorems 1 & 2).

Unified coordinates (per the proofs sheet): ``x_i = W^f_i`` (forward),
``x_j = -W^r_j`` (reverse), ``M = ln(n_f/n_r)``, ``p(x;d) = sigma(x - d + M)``.
The Bennett estimator ``dF_hat`` is the root of the score ``S(d)=0``.

Variance forms (all reported in :class:`BARResult`):
  * ``sandwich``  = ``B / I^2`` with ``B = n_f Var_f[p] + n_r Var_r[p]``  — Theorem 2,
    the headline self-calibrating uncertainty. Always >= 0.
  * ``mbar``      = ``1/I - (1/n_f + 1/n_r)``  — equals pymbar ``bar(..,'MBAR')`` exactly.
  * ``bennett``   = Bennett 1976 Eq. 10a       — equals pymbar ``bar(..,'BAR')`` exactly.
  * ``naive``     = ``1/I``  — the information-equality limit; **only** for the Fig A
    ablation. NEVER report this as the BAR variance (invariant #1).
The three non-naive forms coincide asymptotically (verified ~1% at N>=50).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.optimize import brentq
from scipy.special import expit


def default_log_count_ratio(n_f: int, n_r: int) -> float:
    """M = ln(n_f / n_r)."""
    return float(np.log(n_f / n_r))


def _prep(x_f: ArrayLike, x_r: ArrayLike, M: float | None):
    xf = np.asarray(x_f, dtype=float)
    xr = np.asarray(x_r, dtype=float)
    if M is None:
        M = default_log_count_ratio(xf.size, xr.size)
    return xf, xr, float(M)


def bar_score(d: float, x_f: ArrayLike, x_r: ArrayLike, M: float) -> float:
    """S(d) = sum_r p(x_j) - sum_f (1 - p(x_i)).  Strictly decreasing in d."""
    xf, xr, M = _prep(x_f, x_r, M)
    pf = expit(xf - d + M)
    pr = expit(xr - d + M)
    return float(np.sum(pr) - np.sum(1.0 - pf))


def solve_bar(x_f: ArrayLike, x_r: ArrayLike, M: float | None = None) -> float:
    """Root of S(d)=0 (Brent). S decreases monotonically from +n_r to -n_f, so the
    root is unique (Lemma: strict concavity of the log-likelihood)."""
    xf, xr, M = _prep(x_f, x_r, M)

    def f(d: float) -> float:
        return bar_score(d, xf, xr, M)

    allx = np.concatenate([xf, xr])
    lo, hi = float(allx.min()) - 1.0, float(allx.max()) + 1.0
    step = 1.0
    while f(lo) <= 0.0:  # ensure S(lo) > 0
        lo -= step
        step *= 2.0
    step = 1.0
    while f(hi) >= 0.0:  # ensure S(hi) < 0
        hi += step
        step *= 2.0
    return float(brentq(f, lo, hi, xtol=1e-13))


def fisher_info(d: float, x_f: ArrayLike, x_r: ArrayLike, M: float) -> tuple[float, float, float]:
    """(I_f, I_r, I) where I = sum p(1-p) is the curvature / overlap."""
    xf, xr, M = _prep(x_f, x_r, M)
    pf = expit(xf - d + M)
    pr = expit(xr - d + M)
    If = float(np.sum(pf * (1.0 - pf)))
    Ir = float(np.sum(pr * (1.0 - pr)))
    return If, Ir, If + Ir


def _solve_if_needed(xf: NDArray, xr: NDArray, M: float, d: float | None) -> float:
    return solve_bar(xf, xr, M) if d is None else d


def information_shares(
    x_f: ArrayLike, x_r: ArrayLike, M: float | None = None, d: float | None = None
) -> tuple[float, float]:
    """Signed information-share gradients (Theorem 1): (dF/dmu_f, dF/dmu_r) =
    (+I_f/I, -I_r/I). Absolute values sum to 1."""
    xf, xr, M = _prep(x_f, x_r, M)
    d = _solve_if_needed(xf, xr, M, d)
    If, Ir, I = fisher_info(d, xf, xr, M)
    return If / I, -Ir / I


def _sample_var(p: NDArray) -> float:
    """Unbiased (ddof=1) variance of p; 0 for n<2 (spread unestimable)."""
    return float(np.var(p, ddof=1)) if p.size > 1 else 0.0


def sandwich_variance(
    x_f: ArrayLike, x_r: ArrayLike, M: float | None = None, d: float | None = None
) -> float:
    """Theorem 2 sandwich Var = B / I^2, B = n_f Var_f[p] + n_r Var_r[p].

    Assumes independent samples: ``B`` uses the raw array sizes ``n_f`` and ``n_r``,
    so callers must pass works that have already been decorrelated to the statistical
    inefficiency ``g`` (e.g. via ``pymbar.timeseries.subsample_correlated_data``).
    If the input arrays still contain autocorrelated MD samples, the effective sample
    sizes are smaller than the raw sizes, and ``B/I^2`` will be biased
    (under-estimated variance).
    """
    xf, xr, M = _prep(x_f, x_r, M)
    d = _solve_if_needed(xf, xr, M, d)
    pf = expit(xf - d + M)
    pr = expit(xr - d + M)
    _, _, I = fisher_info(d, xf, xr, M)
    B = xf.size * _sample_var(pf) + xr.size * _sample_var(pr)
    return B / I**2


def mbar_variance(
    x_f: ArrayLike, x_r: ArrayLike, M: float | None = None, d: float | None = None
) -> float:
    """1/I - (1/n_f + 1/n_r). Identical to pymbar bar(..,'MBAR')."""
    xf, xr, M = _prep(x_f, x_r, M)
    d = _solve_if_needed(xf, xr, M, d)
    _, _, I = fisher_info(d, xf, xr, M)
    return 1.0 / I - (1.0 / xf.size + 1.0 / xr.size)


def bennett_variance(
    x_f: ArrayLike, x_r: ArrayLike, M: float | None = None, d: float | None = None
) -> float:
    """Bennett 1976 Eq. 10a. Identical to pymbar bar(..,'BAR'). In our notation
    f_F = 1 - p over forward, f_R = p over reverse."""
    xf, xr, M = _prep(x_f, x_r, M)
    d = _solve_if_needed(xf, xr, M, d)
    nf, nr = xf.size, xr.size
    fF = 1.0 - expit(xf - d + M)
    fR = expit(xr - d + M)
    afF, afF2 = float(np.mean(fF)), float(np.mean(fF**2))
    afR, afR2 = float(np.mean(fR)), float(np.mean(fR**2))
    nrat = 1.0 / nf + 1.0 / nr
    return (afF2 / afF**2) / nf + (afR2 / afR**2) / nr - nrat


def naive_variance(
    x_f: ArrayLike, x_r: ArrayLike, M: float | None = None, d: float | None = None
) -> float:
    """1/I — information-equality limit. ABLATION ONLY (Fig A); never the BAR variance."""
    xf, xr, M = _prep(x_f, x_r, M)
    d = _solve_if_needed(xf, xr, M, d)
    _, _, I = fisher_info(d, xf, xr, M)
    return 1.0 / I


@dataclass
class BARResult:
    delta_f: float
    var_sandwich: float
    var_mbar: float
    var_naive: float
    overlap: float  # I
    info_share_f: float  # +I_f/I
    info_share_r: float  # -I_r/I


def bar_estimate(x_f: ArrayLike, x_r: ArrayLike, M: float | None = None) -> BARResult:
    """Full BAR readout: dF_hat + sandwich/mbar/naive variances + overlap + shares."""
    xf, xr, M = _prep(x_f, x_r, M)
    d = solve_bar(xf, xr, M)
    If, Ir, I = fisher_info(d, xf, xr, M)
    return BARResult(
        delta_f=d,
        var_sandwich=sandwich_variance(xf, xr, M, d),
        var_mbar=mbar_variance(xf, xr, M, d),
        var_naive=1.0 / I,
        overlap=I,
        info_share_f=If / I,
        info_share_r=-Ir / I,
    )
