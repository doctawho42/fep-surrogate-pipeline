"""Comparative cycle-closure detectors and the cross-replicate anchor (peer-review item P1).

The manuscript claims the calibrated per-edge null improves on the field's fixed hysteresis
cutoffs. This module makes that claim testable by flagging the SAME systems three ways:

  A ``flag_calibrated``  -- GLS closure chi^2 with the per-edge sandwich V_e, BH-FDR (ours).
  B ``flag_fixed_cutoff`` -- the field-standard rule: flag if ANY independent cycle closes
                             worse than 1.0 kcal/mol. No variance model at all.
  C ``flag_fixed_se``    -- the same chi^2 machinery but with ONE pooled se per system, which
                             isolates the value of *per-edge adaptivity* from the chi^2 frame.

They are scored against ``anchor_score``: the cross-replicate reproducible-systematic magnitude,
median_e |mean_k z_e^(k)| over independently fitted replicate networks. The anchor is computed
only from replicate residuals, so it does not depend on which detector is being scored (though
A and C do share the standardized-residual machinery with it, which the results doc must state).
Higher anchor = more reproducible systematic error, so a good detector flags high-anchor systems.

Pre-registered (never tuned after the fact): the 1.0 kcal/mol cutoff, BH-FDR 0.05, and the WIN
rule (A must beat BOTH B and C with the paired bootstrap CI excluding zero; otherwise TIE).
"""
from __future__ import annotations

import math
from collections import defaultdict

import numpy as np
from numpy.typing import NDArray

from bar.qc import Edge, benjamini_hochberg, chi2_sf, gls_network

HYSTERESIS_CUTOFF = 1.0  # kcal/mol; pre-registered field-standard fixed hysteresis cutoff


def _tree_path_sum(tree_adj: dict, a: object, b: object) -> float | None:
    """Signed ddG sum along the unique spanning-tree path ``a -> b`` (None if disconnected)."""
    stack: list[tuple[object, float]] = [(a, 0.0)]
    seen = {a}
    while stack:
        node, acc = stack.pop()
        if node == b:
            return acc
        for nbr, w in tree_adj.get(node, ()):
            if nbr not in seen:
                seen.add(nbr)
                stack.append((nbr, acc + w))
    return None


def fundamental_cycle_sums(edges: list[Edge]) -> list[float]:
    """Signed ddG sum around each independent cycle.

    Builds a spanning forest by union-find; every non-tree edge closes exactly one fundamental
    cycle, so the count equals ``gls_network(edges).dof``. Parallel edges are handled naturally
    (the second copy is a non-tree edge). Sign convention matches ``bar.qc._incidence``:
    ``y_e = phi_b - phi_a``, so traversing an edge from ``a`` to ``b`` contributes ``+y_e``.
    """
    parent: dict = {}

    def find(x: object) -> object:
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    tree_adj: dict = defaultdict(list)
    cotree: list[tuple[object, object, float]] = []
    for a, b, y, _se in edges:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb
            tree_adj[a].append((b, float(y)))
            tree_adj[b].append((a, -float(y)))
        else:
            cotree.append((a, b, float(y)))

    sums: list[float] = []
    for a, b, y in cotree:
        path = _tree_path_sum(tree_adj, a, b)
        if path is not None:
            sums.append(y - path)
    return sums


def flag_calibrated(systems: dict, alpha: float = 0.05) -> dict:
    """Detector A (ours): closure chi^2 with the per-edge sandwich V_e, BH-FDR across systems."""
    names: list[str] = []
    pvals: list[float] = []
    for name in sorted(systems):
        fit = gls_network(systems[name])
        if fit.dof < 1:
            continue
        names.append(name)
        pvals.append(chi2_sf(fit.chi2, fit.dof))
    if not names:
        return {}
    flags = benjamini_hochberg(pvals, alpha)
    return {n: bool(f) for n, f in zip(names, flags, strict=True)}


def flag_fixed_cutoff(systems: dict, cutoff: float = HYSTERESIS_CUTOFF) -> dict:
    """Detector B (field standard): flag if any independent cycle closes worse than ``cutoff``."""
    out: dict = {}
    for name in sorted(systems):
        edges = systems[name]
        if gls_network(edges).dof < 1:
            continue
        sums = fundamental_cycle_sums(edges)
        out[name] = bool(sums) and max(abs(s) for s in sums) > cutoff
    return out


def flag_fixed_se(systems: dict, alpha: float = 0.05) -> dict:
    """Detector C: the chi^2 test with ONE pooled se per system (its median per-edge se)."""
    names: list[str] = []
    pvals: list[float] = []
    for name in sorted(systems):
        edges = systems[name]
        if gls_network(edges).dof < 1:
            continue
        se_med = float(np.median([float(se) for _a, _b, _y, se in edges]))
        fit = gls_network([(a, b, y, se_med) for a, b, y, _se in edges])
        names.append(name)
        pvals.append(chi2_sf(fit.chi2, fit.dof))
    if not names:
        return {}
    flags = benjamini_hochberg(pvals, alpha)
    return {n: bool(f) for n, f in zip(names, flags, strict=True)}


def anchor_score(rep_edges: list) -> float:
    """Cross-replicate reproducible-systematic magnitude: ``median_e |mean_k z_e^(k)|``.

    ``rep_edges`` is one edge list per replicate, all listing the SAME edges in the SAME order
    (so the per-edge z vectors align). Systematic error reproduces across replicates and
    survives the mean; sampling noise averages toward zero. Returns NaN if any replicate
    network is acyclic.
    """
    zs: list[NDArray] = []
    for edges in rep_edges:
        fit = gls_network(edges)
        if fit.dof < 1:
            return math.nan
        zs.append(fit.z)
    return float(np.median(np.abs(np.mean(np.vstack(zs), axis=0))))


def auc_flag_vs_anchor(flags: NDArray, anchor: NDArray) -> float:
    """Mann-Whitney AUC of the anchor score of flagged vs unflagged systems.

    0.5 = no discrimination, 1.0 = every flagged system out-scores every unflagged one. A
    detector that flags nothing or everything is undefined here and returns 0.5 by convention.
    """
    from scipy.stats import mannwhitneyu

    f = np.asarray(flags, dtype=bool)
    a = np.asarray(anchor, dtype=float)
    ok = ~np.isnan(a)
    hi, lo = a[ok & f], a[ok & ~f]
    if hi.size == 0 or lo.size == 0:
        return 0.5
    u = float(mannwhitneyu(hi, lo, alternative="greater").statistic)
    return u / (hi.size * lo.size)


def paired_auc_bootstrap(
    flags_a: NDArray, flags_b: NDArray, flags_c: NDArray, anchor: NDArray,
    n_boot: int = 2000, seed: int = 0,
) -> dict:
    """Bootstrap over systems of ``AUC(A) - max(AUC(B), AUC(C))`` (paired: one resample scores
    all three detectors). Pre-registered verdict: WIN iff the point difference is positive AND
    the 95% CI lower bound excludes zero; otherwise TIE."""
    rng = np.random.default_rng(seed)
    A = np.asarray(flags_a, dtype=bool)
    B = np.asarray(flags_b, dtype=bool)
    C = np.asarray(flags_c, dtype=bool)
    s = np.asarray(anchor, dtype=float)

    def diff(idx: NDArray) -> tuple[float, float, float, float]:
        aa = auc_flag_vs_anchor(A[idx], s[idx])
        bb = auc_flag_vs_anchor(B[idx], s[idx])
        cc = auc_flag_vs_anchor(C[idx], s[idx])
        return aa - max(bb, cc), aa, bb, cc

    n = s.size
    point, auc_a, auc_b, auc_c = diff(np.arange(n))
    boot = np.array([diff(rng.integers(0, n, n))[0] for _ in range(n_boot)])
    lo, hi = (float(v) for v in np.percentile(boot, [2.5, 97.5]))
    return {"auc_a": auc_a, "auc_b": auc_b, "auc_c": auc_c, "diff": point,
            "ci_lo": lo, "ci_hi": hi,
            "verdict": "WIN" if (point > 0 and lo > 0) else "TIE"}
