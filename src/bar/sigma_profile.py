"""Rank-transferred learned-sigma profile for the QC calibration sweep (peer-review item P4).

Figure L's panel B shrinks every edge's se by a single uniform factor to stand in for an
overconfident learned sigma. Referees objected that a uniform ``x0.15`` is near-mechanically forced
(it inflates every chi-square by ``1/0.15^2 ~ 44``), so it says little about a real learned head,
whose miscalibration is heterogeneous and overlap-dependent (measured at ``0.09``-``0.20x`` across
the Fig A sweep).

This module applies that measured profile PER EDGE. The catch: the two overlap quantities are
not the same measurement. Fig A's controlled sweep uses a normalized BAR overlap spanning
``0.26``-``0.78``; the OpenFE benchmark reports pymbar's ``smallest_overlap`` spanning
``0.0001``-``0.233``. Their values are not interchangeable, so the transfer is by **rank
(percentile), not by raw value**: an edge at the p-th percentile of the real overlap
distribution receives the ratio the learned head showed at the p-th percentile of the Fig A
sweep. That preserves the ordering and spread of the head's miscalibration, which is what
the objection is about, without pretending the scales are commensurable. It is an
approximation and must be reported as one.
"""
from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

# Verbatim from the committed docs/results_figA.md Panel-A table: (overlap, reported/true se).
# Frozen by the P4 pre-registration; never re-fit or re-tuned.
PROFILE_POINTS: list[tuple[float, float]] = [(0.26, 0.20), (0.46, 0.11), (0.53, 0.09), (0.78, 0.15)]


def rank_transfer(
    overlaps: NDArray, profile: list[tuple[float, float]] = PROFILE_POINTS
) -> NDArray:
    """Per-edge learned-sigma ratio, transferred by percentile rank.

    ``overlaps`` are the real per-edge overlaps (any scale). The returned array has the same shape:
    element ``e`` is the profile ratio at edge ``e``'s percentile within ``overlaps``. Ties share a
    percentile; a single edge or an all-equal input maps to the profile's midpoint by convention
    (``np.interp`` at 0.0 would silently pin every such edge to the lowest-overlap ratio, which
    would be an artifact of degeneracy rather than a measurement).
    """
    ov = np.asarray(overlaps, dtype=float)
    if ov.size == 0:
        return ov
    pts = sorted(profile)
    prof_pct = np.linspace(0.0, 1.0, len(pts))
    prof_ratio = np.array([p[1] for p in pts], dtype=float)

    order = np.argsort(ov, kind="stable")
    ranks = np.empty(ov.size, dtype=float)
    ranks[order] = np.arange(ov.size, dtype=float)
    if ov.size == 1 or np.allclose(ov, ov[0]):
        pct = np.full(ov.size, 0.5)
    else:
        pct = ranks / (ov.size - 1.0)
    return np.interp(pct, prof_pct, prof_ratio)


def shuffled(ratios: NDArray, seed: int) -> NDArray:
    """The same ratio multiset, randomly permuted (seeded).

    The control arm: it holds the marginal distribution of per-edge ratios fixed while destroying
    the association with overlap, separating "the heterogeneity matters" from "only the average
    shrink matters".
    """
    r = np.asarray(ratios, dtype=float)
    return np.random.default_rng(seed).permutation(r)
