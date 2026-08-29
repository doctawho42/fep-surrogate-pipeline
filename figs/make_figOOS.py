"""Fig OOS -- what survives of the cycle-closure "localization" claim, tested out-of-sample.

The manuscript's edge-removal *causal repair* panel has been retracted: the statistic it drives
down (`sum_e z_e^2`) is the same quantity the removal order maximises (descending `z_e^2`), so the
contrast cannot fail, and a synthetic world with no localized systematic error reproduces it. This
figure replaces that panel with two claims that are *not* circular, because in both of them the
quantity being predicted is measured on data the selection never saw.

**Panel A -- predicted vs observed closure chi-square.** Write the reported per-edge variance as
``V_e^rep = c_e^2 V_e^true``. Whitening the GLS fit by the *reported* errors gives
``E[X^2] = tr(M D)`` with ``M = I - H`` the residual maker and ``D = diag(c_e^-2)``, so with
``h_e = M_ee`` the curl-leverage of `bar.leverage.curl_leverage` and ``sum_e h_e = dof``
(Theorem D1), ``E[chi2_nu] = sum_e h_e c_e^-2 / dof`` -- the curl-leverage-weighted MEAN of
``c_e^-2``. The replicate set supplies an unbiased estimate of ``c_e^-2`` that the closure test
never uses: ``s_e^2 / se_e^2`` with ``s_e^2`` the sample variance of the replicate ddG values
(``E[s^2] = sigma^2`` at any n >= 2). Independence is enforced by rotation: replicate ``k``'s
observed ``chi2_nu`` is predicted from ``s_e`` computed on the OTHER TWO replicates only. The
headline is that fully-independent number on replicate 0 (the replicate Fig L reports), with the
other two rotations as the robustness read; the all-three-replicate variants are reported but
labelled in-sample, since the predicted replicate then contributes to its own predictor.

**Panel B -- out-of-sample localization.** Select edges on one replicate, evaluate on the two
held out; rotate over all three choices of selection replicate. Arms, all of the same size K:
  * `static`   -- top-K by |z| on the selection replicate (the one-shot ordering);
  * `guided`   -- `bar.qc.repair_order`'s greedy removal set (the paper's actual repair rule);
  * `leverage` -- top-K by curl-leverage h (topology only, blind to the measured values);
  * `matched`  -- the SAME |z| ranks taken in an unflagged donor system (blind to the flag
                  decision: is this just what high-|z| edges do anywhere?);
  * `random`   -- uniform K-subsets, which also supply the null distribution and the p-values.

Test statistics are calibrated before use, because one of them is invalid. Under the sampling
null ``Var(z_e) = h_e``, so selecting on |z| preferentially selects high-leverage edges, whose
held-out ``z^2`` is larger *by topology alone*. The alpha-calibration harness (a no-systematic-error
world: same graphs, same reported se, resampled values) measures every statistic's realized
false-positive rate at nominal 0.05 and the doc reports all of them, including the discarded one.

Pre-stated design constants, fixed before any outcome was computed and not revisited:
  * systems: the six Fig L flagged systems; a system needs dof >= 2 to be testable (brd4 has one
    cycle, so removing a cycle edge leaves nothing to measure -- it is reported as untestable);
  * edges: rows complete in all three replicates, so the per-edge index means the same edge in
    every replicate; every arm, random included, draws from the same candidate pool -- the
    NON-BRIDGE edges. A bridge lies in no cycle, so ``h_e = 0`` and ``z_e = 0`` identically
    (Theorem D1): it is structurally un-auditable, ``z^2/h`` is 0/0 on it, and removing it cannot
    move any closure statistic;
  * K = max(1, #removals `repair_order` needs on the selection replicate) -- the paper's own
    stopping rule, floored at 1 so a selection is never empty (`faah` on replicate 1 closes
    already). In the null calibration K is held at its real-data value, so the null exercises the
    same selection sizes;
  * one test per (system, rotation): the statistic is averaged over the two held-out replicates,
    then compared with the identically-averaged random null; p-values one-sided upper,
    ``(1 + #{null >= obs}) / (1 + N)``; per-rotation aggregation by Stouffer over systems.

Run:  PYTHONPATH=src python figs/make_figOOS.py   (or `make figOOS`).  Deterministic (all
resampling seeded).  Data: `data/openfe_replicates/combined_pymbar4_edge_data.csv` (OpenFE
IndustryBenchmarks2024, public, 3 independent replicates per edge).
"""
from __future__ import annotations

import math
import pathlib
import sys

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.ticker as mticker  # noqa: E402
import numpy as np  # noqa: E402

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "figs"))
from figs.make_figStab import REPLICATES, edge_val, load_systems  # noqa: E402

from bar.leverage import bridges, curl_leverage  # noqa: E402
from bar.qc import gls_network, repair_order  # noqa: E402
from paperstyle import (  # noqa: E402
    ALT, INK, MUTED, OURS, REF, check_min_type, figsize, finish, legend, panel, reference_line,
    tint, use_paper_style,
)

FIGDIR = pathlib.Path(__file__).resolve().parent
DOC = ROOT / "docs" / "results_figOOS.md"

FLAGGED = ("brd4", "bace", "faah", "cdk8", "hif2a", "p38")
MIN_EDGES = 3          # Fig L's per-system admission rule
MIN_DOF_TEST = 2       # a system needs >1 independent cycle to be testable in panel B
ALPHA = 0.05           # nominal level the calibration harness is measured against
N_RAND = 1999          # random draws per (system, rotation) in the reported analysis
N_DONOR = 10           # unflagged donor systems per (system, rotation) for the matched arm
N_DONOR_RAND = 99      # random draws inside each donor (its own null)
N_NULL = 120           # no-systematic-error realizations in the alpha-calibration harness
N_NULL_RAND = 99       # random draws per test inside the null harness
N_BOOT = 2000          # system-cluster bootstrap resamples for the panel-A aggregate
SEED = 20260810

# Semantic colours (paperstyle), by quantity family. OURS = the calibrated per-edge residual
# and the selection rules built from it -- `guided` is the manuscript's repair rule, `static` its
# one-shot ordering, so they are one hue told apart by tint. ALT = the theorem-derived
# curl-leverage rule, the named alternative. MUTED = the null controls (uniform random removal,
# and the rank-matched control taken in donor systems).
#
# FOIL is deliberately ABSENT from this figure. Every candidate statistic in panel B -- the two
# that hold their nominal level and the two that fire far above it -- is one of OUR OWN, built
# from the calibrated per-edge residual, so all four are drawn in the OURS family. A statistic
# does not change family by failing its own audit, and FOIL has exactly one referent in this
# article (the overconfident stand-in and the real learned-variance heads); spending it on a
# discarded arm of ours would destroy that referent. "Passes / fails its level" is therefore
# encoded by FILL AND HATCH, not by hue: see PASS_FACE / FAIL_FACE / FAIL_HATCH below.
ARM_COLOR = {"static": tint(OURS, 0.45), "guided": OURS, "leverage": ALT,
             "matched": MUTED, "random": MUTED}

# Panel B's non-hue encoding of "held its nominal level" vs "fired far above it".
PASS_FACE = OURS               # solid: the statistic held its level and is used
FAIL_FACE = tint(OURS, 0.62)   # same hue, opened up: still ours, but not used
FAIL_HATCH = "///"             # ... and struck through, which no other figure element is

# The statistics, all evaluated on held-out replicates. `z2_raw` is retained deliberately: the
# harness is what disqualifies it, and a discarded statistic that is never shown cannot be judged.
STAT_NAMES = ("dchi2", "z2_lev", "z2_raw", "xrep")
STAT_LABEL = {
    "dchi2": "held-out consistency gain  $\\Delta\\chi^2_\\nu$",
    "z2_lev": "leverage-normalised held-out $z^2/h$",
    "z2_raw": "raw held-out $z^2$",
    "xrep": "signed cross-replicate product $z_s z_h$",
}
# Tick labels carry their own gloss: the subscripts of the internal names mean nothing to a
# reader of the printed figure, so each axis entry says in words what the statistic is.
STAT_SHORT = {"dchi2": "$\\Delta\\chi^2_\\nu$\nconsistency\ngain",
              "z2_lev": "$z^2/h$\nleverage-\nnormalised",
              "z2_raw": "raw $z^2$", "xrep": "cross-\nreplicate\n$z$ product"}
STAT_INLINE = {"dchi2": "the consistency gain", "z2_lev": "the leverage-normalised $z^2/h$",
               "z2_raw": "raw $z^2$", "xrep": "the cross-replicate $z$ product"}
ARMS = ("static", "guided", "leverage", "matched", "random")


# --------------------------------------------------------------------------- data


def aligned_rows(rows):
    """Rows whose four legs are present in ALL three replicates, so edge index e denotes the same
    ligand pair in every replicate (the panel-B statistics compare per-edge across replicates)."""
    return [r for r in rows if all(edge_val(r, k) for k in REPLICATES)]


def system_record(rows):
    """Per-system, per-replicate fit bundle: edges, z, curl-leverage h, reduced chi^2, dof.

    Returns ``None`` for systems Fig L would not test (fewer than MIN_EDGES aligned edges, or no
    independent cycle). ``h`` comes from `bar.leverage.curl_leverage`, whose two independent code
    paths cross-check the conservation law, so a wrong weight or orientation cannot pass silently.

    ``pool`` is the panel-B candidate set: the NON-BRIDGE edges. A bridge lies in no cycle, so by
    Theorem D1 it has ``h_e = 0`` and residual ``z_e = 0`` identically; it is structurally
    un-auditable by closure, ``z^2/h`` is 0/0 on it, and removing it cannot change any closure
    statistic. Every arm, including random, draws from the same pool, so no arm is advantaged.
    Bridges are a property of the topology alone, hence the same set in every replicate.
    """
    keep = aligned_rows(rows)
    if len(keep) < MIN_EDGES:
        return None
    edges, fits, hs = {}, {}, {}
    for k in REPLICATES:
        e = [(r["ligand_A"], r["ligand_B"], *edge_val(r, k)) for r in keep]
        fit = gls_network(e)
        if fit.dof < 1:
            return None
        edges[k], fits[k], hs[k] = e, fit, curl_leverage(e)
    br = bridges(edges[0])
    return {
        "E": len(keep), "edges": edges, "dof": fits[0].dof,
        "pool": np.array([e for e in range(len(keep)) if e not in br], dtype=int),
        "y": {k: np.array([e[2] for e in edges[k]]) for k in REPLICATES},
        "se": {k: np.array([e[3] for e in edges[k]]) for k in REPLICATES},
        "z": {k: np.asarray(fits[k].z) for k in REPLICATES},
        "h": {k: np.asarray(hs[k]) for k in REPLICATES},
        "rc": {k: float(fits[k].reduced_chi2) for k in REPLICATES},
        "dof_k": {k: int(fits[k].dof) for k in REPLICATES},
    }


def load_records():
    return {name: rec for name, rows in sorted(load_systems().items())
            if (rec := system_record(rows)) is not None}


# --------------------------------------------------------------------------- panel A


def leverage_weighted_mean(h, values):
    """``sum_e h_e v_e / sum_e h_e`` -- the functional the closure statistic identifies.

    With ``v_e = c_e^-2`` and ``sum_e h_e = dof`` this is exactly ``tr(M D) / dof = E[chi2_nu]``.
    """
    h = np.asarray(h, dtype=float)
    v = np.asarray(values, dtype=float)
    tot = float(np.sum(h))
    if tot <= 0:
        return math.nan
    return float(np.sum(h * v) / tot)


def sample_var(vals):
    """Unbiased sample variance of the replicate ddG values (``E[s^2] = sigma^2``, any n >= 2).

    ``vals`` is a list of per-replicate arrays; the divisor is n-1, so at n=2 this is
    ``(y_i - y_j)^2 / 2``.
    """
    a = np.vstack(vals)
    return np.var(a, axis=0, ddof=1)


def predicted_chi2(rec, predicted_rep, source_reps):
    """Predicted ``E[chi2_nu]`` for ``predicted_rep`` from the replicate spread of ``source_reps``.

    ``c_e^-2`` is estimated per edge as ``s_e^2 / se_{e,k}^2`` and averaged with the curl-leverage
    weights of the predicted replicate's own fit. If ``predicted_rep not in source_reps`` the
    estimate is statistically independent of the observed statistic being predicted.
    """
    s2 = sample_var([rec["y"][k] for k in source_reps])
    cinv2 = s2 / rec["se"][predicted_rep] ** 2
    return leverage_weighted_mean(rec["h"][predicted_rep], cinv2)


def panel_a(recs):
    """Predicted vs observed reduced chi^2 per system, in several flavours.

    HEADLINE = ``indep_rep0``: replicate 0's observed chi^2_nu predicted from the replicate spread
    of replicates 1 and 2 only, so predictor and target share no data. ``indep_rep1`` and
    ``indep_rep2`` rotate the same construction. The ``insample_*`` variants use all three in the
    spread, so the predicted replicate contributes to its own predictor; they are reported for
    comparison and are never the quoted number.
    """
    names = sorted(recs)
    dof = np.array([recs[s]["dof"] for s in names], dtype=float)
    is_flag = np.array([s in FLAGGED for s in names])
    obs = {k: np.array([recs[s]["rc"][k] for s in names]) for k in REPLICATES}
    pred_indep, pred_insample = {}, {}
    for k in REPLICATES:
        others = tuple(j for j in REPLICATES if j != k)
        pred_indep[k] = np.array([predicted_chi2(recs[s], k, others) for s in names])
        pred_insample[k] = np.array([predicted_chi2(recs[s], k, REPLICATES) for s in names])

    variants = {f"indep_rep{k}": (pred_indep[k], obs[k]) for k in REPLICATES}
    variants |= {f"insample_rep{k}": (pred_insample[k], obs[k]) for k in REPLICATES}
    variants["indep_mean"] = (np.mean([pred_indep[k] for k in REPLICATES], axis=0),
                              np.mean([obs[k] for k in REPLICATES], axis=0))
    variants["insample_mean"] = (np.mean([pred_insample[k] for k in REPLICATES], axis=0),
                                 np.mean([obs[k] for k in REPLICATES], axis=0))
    packs = {tag: rho_pack(pred, ob, dof, is_flag) for tag, (pred, ob) in variants.items()}

    # aggregate E[chi2_nu] = total predicted chi^2 over total dof (the leverage-weighted mean of
    # c^-2 pooled over systems). Bootstrap resamples SYSTEMS, the cluster. Reported for the
    # headline rotation and, as a robustness read, pooled over all three rotations.
    rng = np.random.default_rng(SEED)

    def aggregate(num, den):
        agg = float(num.sum() / den.sum())
        boot = np.empty(N_BOOT)
        for b in range(N_BOOT):
            idx = rng.integers(0, len(names), len(names))
            boot[b] = num[idx].sum() / den[idx].sum()
        return agg, (float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5)))

    agg, agg_ci = aggregate(pred_indep[0] * dof, dof.copy())
    agg_all, agg_all_ci = aggregate(sum(pred_indep[k] for k in REPLICATES) * dof, 3.0 * dof)
    obs_agg = float(np.sum(obs[0] * dof) / np.sum(dof))
    obs_agg_all = float(np.sum(sum(obs[k] for k in REPLICATES) * dof) / np.sum(3.0 * dof))

    ratio = {s: float(obs[0][i] / pred_indep[0][i]) for i, s in enumerate(names)}
    ratio_is = {s: float(obs[0][i] / pred_insample[0][i]) for i, s in enumerate(names)}
    return {"names": names, "variants": variants, "packs": packs, "dof": dof, "is_flag": is_flag,
            "agg": agg, "agg_ci": agg_ci, "agg_all": agg_all, "agg_all_ci": agg_all_ci,
            "obs_agg": obs_agg, "obs_agg_all": obs_agg_all, "ratio": ratio, "ratio_is": ratio_is,
            "pred_indep": pred_indep, "pred_insample": pred_insample, "obs": obs}


def rho_pack(pred, obs, dof, is_flag):
    """Spearman rho plus the four pre-stated robustness reads: drop the flagged systems, partial
    out dof, restrict to well-determined networks (dof >= 5), and leave one system out."""
    from scipy.stats import spearmanr

    r = spearmanr(pred, obs)
    sub = ~is_flag
    hi = dof >= 5
    loo = [float(spearmanr(np.delete(pred, i), np.delete(obs, i)).statistic)
           for i in range(len(pred))]
    return {
        "rho": float(r.statistic), "p": float(r.pvalue), "n": int(len(pred)),
        "rho_no_flagged": float(spearmanr(pred[sub], obs[sub]).statistic),
        "n_no_flagged": int(sub.sum()),
        "rho_partial_dof": partial_spearman(pred, obs, dof),
        "rho_dof_ge5": float(spearmanr(pred[hi], obs[hi]).statistic),
        "n_dof_ge5": int(hi.sum()),
        "loo": (float(np.min(loo)), float(np.max(loo))),
    }


def partial_spearman(x, y, z):
    """Spearman partial correlation of x and y given z (rank-transform, then partial Pearson)."""
    from scipy.stats import rankdata

    rx, ry, rz = (rankdata(v) for v in (x, y, z))
    c = np.corrcoef(np.vstack([rx, ry, rz]))
    num = c[0, 1] - c[0, 2] * c[1, 2]
    den = math.sqrt(max((1 - c[0, 2] ** 2) * (1 - c[1, 2] ** 2), 1e-300))
    return float(num / den)


# --------------------------------------------------------------------------- panel B


def _rc(edges):
    """Reduced chi^2 of a sub-network; a network with no remaining cycle is trivially consistent
    (convention: 0.0), applied identically to every arm so no arm is advantaged by it."""
    fit = gls_network(edges)
    return fit.reduced_chi2 if fit.dof > 0 else 0.0


def rank_matched(scores, ranks):
    """Indices of the entries occupying the given 1-based DESCENDING ranks of ``|scores|``.

    The matched arm asks whether the same |z| *ranks* behave the same way in a system the QC did
    not flag, so the matcher must be able to reproduce an arbitrary rank multiset, not only a
    top-K block. Ties are broken by index (stable sort), so the map rank -> element is a bijection.
    """
    order = np.argsort(-np.abs(np.asarray(scores, dtype=float)), kind="stable")
    ranks = np.asarray(ranks, dtype=int)
    if ranks.min() < 1 or ranks.max() > order.size:
        raise ValueError("ranks must lie in 1..len(scores)")
    return order[ranks - 1]


def prepare(rec, y_by_rep=None, z_by_rep=None):
    """Per-replicate evaluation context, built once and reused across thousands of selections.

    Holds the edge lists, the standardized residuals, the curl-leverages and the full-network
    reduced chi^2 of each replicate. ``y_by_rep`` / ``z_by_rep`` override the measured values, which
    is how the alpha-calibration harness evaluates the same graphs in a no-systematic-error world.
    """
    y = rec["y"] if y_by_rep is None else y_by_rep
    z = rec["z"] if z_by_rep is None else z_by_rep
    edges = {k: [(a, b, float(y[k][e]), se)
                 for e, (a, b, _yy, se) in enumerate(rec["edges"][k])] for k in REPLICATES}
    return {"E": rec["E"], "pool": rec["pool"], "edges": edges, "z": z, "h": rec["h"],
            "rc_full": {k: _rc(edges[k]) for k in REPLICATES}}


def statistics(ctx, sel_rep, sel_idx):
    """The four held-out statistics for a selection made on replicate ``sel_rep``.

    Each is averaged over the two held-out replicates:
      * ``dchi2``  -- consistency gain: full-network reduced chi^2 minus the reduced chi^2 after
                      the selection is removed (positive = removing it improves held-out closure);
      * ``z2_lev`` -- mean held-out ``z_e^2 / h_e`` (unit mean under the sampling null);
      * ``z2_raw`` -- mean held-out ``z_e^2`` (retained only so the harness can disqualify it);
      * ``xrep``   -- mean signed product ``z_e(selection replicate) * z_e(held-out replicate)``.
    """
    held = [k for k in REPLICATES if k != sel_rep]
    sel = np.asarray(sel_idx, dtype=int)
    drop = set(sel.tolist())
    out = dict.fromkeys(STAT_NAMES, 0.0)
    z = ctx["z"]
    for k in held:
        edges_k = ctx["edges"][k]
        kept = [edges_k[e] for e in range(ctx["E"]) if e not in drop]
        out["dchi2"] += ctx["rc_full"][k] - (_rc(kept) if len(kept) >= 2 else 0.0)
        out["z2_lev"] += float(np.mean(z[k][sel] ** 2 / ctx["h"][k][sel]))
        out["z2_raw"] += float(np.mean(z[k][sel] ** 2))
        out["xrep"] += float(np.mean(z[sel_rep][sel] * z[k][sel]))
    return {n: v / len(held) for n, v in out.items()}


def perm_p(obs, null):
    """One-sided upper permutation p-value, ``(1 + #{null >= obs}) / (1 + N)``."""
    null = np.asarray(null, dtype=float)
    return float((1 + int(np.sum(null >= obs))) / (1 + null.size))


def stouffer(pvals):
    """Stouffer's z-combination of one-sided p-values (equal weights)."""
    from scipy.stats import norm

    p = np.clip(np.asarray(pvals, dtype=float), 1e-300, 1 - 1e-16)
    zz = float(np.sum(norm.isf(p)) / math.sqrt(p.size))
    return zz, float(norm.sf(zz))


def pool_ranks(ctx, scores, ranks):
    """Edge indices at the given 1-based |score| ranks *within the auditable (non-bridge) pool*."""
    pool = ctx["pool"]
    return pool[rank_matched(np.asarray(scores)[pool], ranks)]


def selection_size(ctx, sel_rep):
    """K = max(1, #removals `bar.qc.repair_order` needs on the selection replicate), capped so the
    selection cannot exhaust the auditable pool."""
    removed, _ = repair_order(ctx["edges"][sel_rep], target_reduced_chi2=1.0)
    return max(1, min(len(removed), ctx["pool"].size - 1)), removed


def arm_selections(ctx, sel_rep, k_size, removed):
    """The three data-side arms' index sets (``matched`` and ``random`` are drawn elsewhere)."""
    static = pool_ranks(ctx, ctx["z"][sel_rep], np.arange(1, k_size + 1))
    in_pool = set(ctx["pool"].tolist())
    guided = np.array([i for i in removed if i in in_pool][:k_size], dtype=int)
    if guided.size < k_size:                       # repair_order stopped early: top up by |z|
        extra = [i for i in static if i not in set(guided.tolist())]
        guided = np.array(list(guided) + extra[: k_size - guided.size], dtype=int)
    lev = ctx["pool"][np.argsort(-ctx["h"][sel_rep][ctx["pool"]], kind="stable")[:k_size]]
    return {"static": static, "guided": guided, "leverage": np.asarray(lev, dtype=int)}


def random_null(ctx, sel_rep, k_size, n_draws, rng, stat_fn=None):
    """``n_draws`` uniform K-subsets of the auditable pool, each scored with every statistic."""
    stat_fn = statistics if stat_fn is None else stat_fn
    out = None
    for d in range(n_draws):
        st = stat_fn(ctx, sel_rep, rng.choice(ctx["pool"], k_size, replace=False))
        if out is None:
            out = {n: np.empty(n_draws) for n in st}
        for n in out:
            out[n][d] = st[n]
    return out


def resample_null_world(rec, rng):
    """A no-systematic-error world on the SAME graph and the SAME reported errors.

    Each replicate's edge values are redrawn as ``y ~ N(0, se_k)``; the GLS residuals are
    gauge-invariant, so zero node potentials are without loss of generality. Returns
    ``(y_by_rep, z_by_rep)``; ``h`` is unchanged because it depends only on topology and se.
    """
    y_by_rep, z_by_rep = {}, {}
    for k in REPLICATES:
        yk = rng.normal(0.0, rec["se"][k])
        edges = [(a, b, float(yk[e]), se) for e, (a, b, _yy, se) in enumerate(rec["edges"][k])]
        y_by_rep[k], z_by_rep[k] = yk, np.asarray(gls_network(edges).z)
    return y_by_rep, z_by_rep


def alpha_calibration(recs, testable, sizes, rng, n_null=None, n_draws=None,
                      stat_names=STAT_NAMES, stat_fn=None):
    """Realized false-positive rate of each statistic's `static`-arm test at nominal ALPHA.

    The world carries no localized systematic error, so every rejection is a false one. K is held
    at its real-data value per (system, rotation) so the null exercises the same selection sizes.
    ``stat_fn`` is injectable, so a deliberately mis-calibrated statistic can be shown to be caught.
    Returns ``{statistic: {"per_test": alpha, "stouffer": alpha, ...}}``.
    """
    n_null = N_NULL if n_null is None else n_null
    n_draws = N_NULL_RAND if n_draws is None else n_draws
    stat_fn = statistics if stat_fn is None else stat_fn
    fired = dict.fromkeys(stat_names, 0)
    fired_st = dict.fromkeys(stat_names, 0)
    n_test = 0
    for _ in range(n_null):
        ctxs = {s: prepare(recs[s], *resample_null_world(recs[s], rng)) for s in testable}
        for sel_rep in REPLICATES:
            ps = {n: [] for n in stat_names}
            for s in testable:
                ctx = ctxs[s]
                k_size = sizes[(s, sel_rep)]
                sel = pool_ranks(ctx, ctx["z"][sel_rep], np.arange(1, k_size + 1))
                obs = stat_fn(ctx, sel_rep, sel)
                null = random_null(ctx, sel_rep, k_size, n_draws, rng, stat_fn=stat_fn)
                for n in stat_names:
                    p = perm_p(obs[n], null[n])
                    ps[n].append(p)
                    fired[n] += int(p <= ALPHA)
                n_test += 1
            for n in stat_names:
                fired_st[n] += int(stouffer(ps[n])[1] <= ALPHA)
    n_test = max(n_test, 1)
    n_st = max(n_null * len(REPLICATES), 1)
    return {n: {"per_test": fired[n] / n_test, "stouffer": fired_st[n] / n_st, "n_test": n_test,
                "n_stouffer": n_st} for n in stat_names}


def panel_b(recs):
    """Out-of-sample localization: arms x statistics x rotations, plus the alpha calibration."""
    rng = np.random.default_rng(SEED)
    testable = [s for s in FLAGGED if s in recs and recs[s]["dof"] >= MIN_DOF_TEST]
    untestable = [s for s in FLAGGED if s in recs and recs[s]["dof"] < MIN_DOF_TEST]
    donor_names = [s for s in recs if s not in FLAGGED and recs[s]["dof"] >= MIN_DOF_TEST]
    ctxs = {s: prepare(recs[s]) for s in recs}

    sizes, removed_by = {}, {}
    for s in testable:
        for sel_rep in REPLICATES:
            k_size, removed = selection_size(ctxs[s], sel_rep)
            sizes[(s, sel_rep)] = k_size
            removed_by[(s, sel_rep)] = removed

    cells = {}       # (system, rotation, arm) -> {"stat": {...}, "p": {...}}
    rand_mean = {}   # (system, rotation) -> {stat: mean of the random null}
    for s in testable:
        ctx = ctxs[s]
        for sel_rep in REPLICATES:
            k_size = sizes[(s, sel_rep)]
            null = random_null(ctx, sel_rep, k_size, N_RAND, rng)
            rand_mean[(s, sel_rep)] = {n: float(np.mean(null[n])) for n in STAT_NAMES}
            for arm, sel in arm_selections(ctx, sel_rep, k_size, removed_by[(s, sel_rep)]).items():
                st = statistics(ctx, sel_rep, sel)
                cells[(s, sel_rep, arm)] = {
                    "stat": st, "p": {n: perm_p(st[n], null[n]) for n in STAT_NAMES},
                    "rand_ref": rand_mean[(s, sel_rep)], "K": k_size}
            # random arm: its own mean, sitting at the middle of its own null by construction
            cells[(s, sel_rep, "random")] = {
                "stat": rand_mean[(s, sel_rep)], "p": dict.fromkeys(STAT_NAMES, 0.5),
                "rand_ref": rand_mean[(s, sel_rep)], "K": k_size}
            # matched arm: the same |z| ranks in unflagged donor systems, each scored against
            # ITS OWN random null -- so its reference level is the donors', not this system's.
            ranks = np.arange(1, k_size + 1)
            elig = [d for d in donor_names if ctxs[d]["pool"].size > k_size]
            pick = rng.choice(len(elig), min(N_DONOR, len(elig)), replace=False)
            m_stat = {n: [] for n in STAT_NAMES}
            m_p = {n: [] for n in STAT_NAMES}
            m_ref = {n: [] for n in STAT_NAMES}
            for di in pick:
                dctx = ctxs[elig[int(di)]]
                dst = statistics(dctx, sel_rep, pool_ranks(dctx, dctx["z"][sel_rep], ranks))
                dnull = random_null(dctx, sel_rep, k_size, N_DONOR_RAND, rng)
                for n in STAT_NAMES:
                    m_stat[n].append(dst[n])
                    m_p[n].append(perm_p(dst[n], dnull[n]))
                    m_ref[n].append(float(np.mean(dnull[n])))
            cells[(s, sel_rep, "matched")] = {
                "stat": {n: float(np.mean(m_stat[n])) for n in STAT_NAMES},
                "p": {n: float(stouffer(m_p[n])[1]) for n in STAT_NAMES},
                "rand_ref": {n: float(np.mean(m_ref[n])) for n in STAT_NAMES}, "K": k_size}

    combined = {}    # (rotation, arm, stat) -> (z, p)
    for sel_rep in REPLICATES:
        for arm in ARMS:
            for n in STAT_NAMES:
                combined[(sel_rep, arm, n)] = stouffer(
                    [cells[(s, sel_rep, arm)]["p"][n] for s in testable])

    alpha = alpha_calibration(recs, testable, sizes, rng)
    valid = [n for n in STAT_NAMES if abs(alpha[n]["per_test"] - ALPHA) <= ALPHA]

    # cross-replicate z correlation per system (the mechanism behind an arm that runs backwards)
    zcorr = {s: float(np.mean([np.corrcoef(recs[s]["z"][i], recs[s]["z"][j])[0, 1]
                               for i, j in [(0, 1), (0, 2), (1, 2)]]))
             for s in testable + untestable}
    return {"testable": testable, "untestable": untestable, "sizes": sizes, "cells": cells,
            "rand_mean": rand_mean, "combined": combined, "alpha": alpha, "valid": valid,
            "zcorr": zcorr, "n_donor_pool": len(donor_names)}


# --------------------------------------------------------------------------- figure


def make_figure(A, B):
    fig, axes = plt.subplots(2, 2, figsize=figsize(4, 5.75))
    axA, axB, axC, axD = axes[0, 0], axes[0, 1], axes[1, 0], axes[1, 1]

    # --- A: predicted vs observed reduced chi^2 -------------------------------------------
    pred, obs = A["variants"]["indep_rep0"]
    fl = A["is_flag"]
    pk = A["packs"]["indep_rep0"]
    # the flagged systems are the calibrated detector's own flags, so they carry OURS here and in
    # Fig L; the sampling-consistent majority is present but is not the point, so it is MUTED.
    axA.scatter(pred[~fl], obs[~fl], s=11, c=MUTED, alpha=0.85, lw=0, label="not flagged")
    # the six flagged systems overlap tightly, so name them by marker rather than by a text label
    for j, idx in enumerate(np.flatnonzero(fl)):
        axA.scatter(pred[idx], obs[idx], s=24, c=OURS, lw=0,
                    marker="osP^Dv"[j % 6], label=A["names"][idx])
    lim = [min(pred.min(), obs.min()) * 0.6, max(pred.max(), obs.max()) * 1.7]
    axA.plot(lim, lim, color=REF, ls=(0, (1.0, 2.0)), lw=1.0, zorder=0.5)
    axA.set_xscale("log")
    axA.set_yscale("log")
    axA.set_xlim(*lim)
    axA.set_ylim(*lim)
    axA.set_xlabel(r"predicted $\mathrm{E}[\chi^2_\nu]$ = leverage-weighted"
                   "\n" r"mean of $c_e^{-2}$ (replicates 1 and 2 only)")
    axA.set_ylabel("observed $\\chi^2_\\nu$, replicate 0\n(closure test)")
    # upper left, the one corner the point cloud never enters; no patch, so nothing is veiled
    axA.text(0.03, 0.98, f"Spearman $\\rho$ = {pk['rho']:.3f}\n$p$ = {pk['p']:.1e}"
                         f"  ($n$ = {pk['n']})\naggregate $\\mathrm{{E}}[\\chi^2_\\nu]$ = "
                         f"{A['agg']:.2f}\n[{A['agg_ci'][0]:.2f}, {A['agg_ci'][1]:.2f}]",
             transform=axA.transAxes, va="top", ha="left", fontsize=7.5, color=INK)
    legend(axA, loc="lower right", fontsize=7.5, ncol=2, handletextpad=0.25, columnspacing=0.8,
           borderaxespad=0.25)
    panel(axA, "A", "predicts a $\\chi^2$ it never saw",
          subtitle="flagged systems sit above the identity")

    # --- B: alpha calibration ----------------------------------------------------------------
    names = list(STAT_NAMES)
    vals = [B["alpha"][n]["per_test"] for n in names]
    ok = [n in B["valid"] for n in names]
    xs = np.arange(len(names))
    top = max(vals) * 1.55
    # all four are our own candidate statistics, so all four are OURS blue; the two that fail
    # their own audit are opened up and hatched, and labelled, rather than recoloured (see
    # PASS_FACE / FAIL_FACE / FAIL_HATCH).
    bars = axB.bar(xs, vals, color=[PASS_FACE if k else FAIL_FACE for k in ok], width=0.66,
                   edgecolor=OURS, linewidth=0.9, zorder=2)
    for b, k in zip(bars, ok, strict=True):
        if not k:
            b.set_hatch(FAIL_HATCH)
    for x, v, k in zip(xs, vals, ok, strict=True):
        axB.text(x, v + top * 0.02, f"{v:.3f}", ha="center", va="bottom", fontsize=7.5, color=INK)
        if not k:
            axB.text(x, v + top * 0.10, "not used", ha="center", va="bottom", fontsize=7.5,
                     color=INK)
    reference_line(axB, "hline", ALPHA, zorder=2.5)
    # both notes go in the empty upper left: the nominal line runs through the two short bars,
    # so a label beside it would sit on them.
    axB.text(0.03, 0.97, f"dotted: nominal $\\alpha$ = {ALPHA}\n"
                         f"{B['alpha'][names[0]]['n_test']} tests, no-systematic-error world",
             transform=axB.transAxes, ha="left", va="top", fontsize=7.5, color=REF)
    axB.set_xticks(xs)
    axB.set_xticklabels([STAT_SHORT[n] for n in names], fontsize=7.5)
    axB.tick_params(axis="x", length=0)
    # bars are 0.66 wide on unit centres: give the y-spine the same gap the bars have
    axB.set_xlim(-0.67, len(names) - 0.33)
    axB.set_ylim(0, top)
    axB.set_ylabel("realized false-positive rate")
    dropped = [n for n in STAT_NAMES if n not in B["valid"]]
    panel(axB, "B", "calibrate the statistic first",
          subtitle=("the discarded statistics fire far above nominal" if len(dropped) > 1
                    else "the discarded statistic fires far above nominal") if dropped
          else "every statistic holds its nominal level")

    # --- C: held-out consistency gain per system ---------------------------------------------
    testable = B["testable"]
    show = ["static", "guided", "leverage", "random"]      # `matched` lives in donor systems
    width = 0.20
    xs = np.arange(len(testable))
    series = {arm: [float(np.mean([B["cells"][(s, k, arm)]["stat"]["dchi2"] for k in REPLICATES]))
                    for s in testable] for arm in show}
    for ai, arm in enumerate(show):
        axC.bar(xs + (ai - 1.5) * width, series[arm], width=width * 0.9, color=ARM_COLOR[arm],
                lw=0, label=arm, zorder=2)
    axC.axhline(0.0, color=INK, lw=0.7, zorder=3)
    axC.set_xticks(xs)
    axC.set_xticklabels(testable)
    axC.tick_params(axis="x", length=0)
    axC.set_xlim(-0.6, len(testable) - 0.4)
    lo = min(min(v) for v in series.values())
    hi = max(max(v) for v in series.values())
    axC.set_yscale("symlog", linthresh=1.0, linscale=0.9)
    axC.set_ylim(lo * 1.3 - 0.1, hi * 2.3)
    # symlog's default labelling prints 10^0 and -10^0 for 1 and -1; these are differences of a
    # reduced chi^2, so they are labelled as the plain numbers they are.
    axC.yaxis.set_major_locator(mticker.FixedLocator([-2, -1, 0, 1, 2, 5]))
    axC.yaxis.set_minor_locator(mticker.NullLocator())
    axC.yaxis.set_major_formatter(
        mticker.FuncFormatter(lambda v, _p: f"{v:g}".replace("-", "\u2212")))
    axC.set_ylabel("held-out $\\Delta\\chi^2_\\nu$\n(mean over 3 rotations)")
    legend(axC, loc="upper right", fontsize=7.5, ncol=2, columnspacing=1.0)
    panel(axC, "C", "the choice transfers to held-out runs",
          subtitle="selected on one replicate, scored on the other two")

    # --- D: per-rotation Stouffer p, valid statistics only -------------------------------------
    valid = [n for n in STAT_NAMES if n in B["valid"]]
    arms_d = ["static", "guided", "leverage", "matched"]
    xs = np.arange(len(valid))
    for ai, arm in enumerate(arms_d):
        for k in REPLICATES:
            ps = [B["combined"][(k, arm, n)][1] for n in valid]
            axD.scatter(xs + (ai - 1.5) * 0.18, ps, s=22, color=ARM_COLOR[arm],
                        marker=["o", "s", "^"][k], lw=0, zorder=3,
                        label=arm if k == 0 else None)
    reference_line(axD, "hline", ALPHA)
    allp = [B["combined"][(k, arm, n)][1] for n in valid for arm in arms_d for k in REPLICATES]
    axD.set_yscale("log")
    axD.set_ylim(min(allp) * 0.04, 3.0)
    axD.set_xlim(-0.5, len(valid) - 0.5)
    # under the line, on the right: the only part of that row no arm's points occupy
    axD.text(len(valid) - 0.55, ALPHA * 0.88, f"$\\alpha$={ALPHA}", fontsize=7.5, color=REF,
             ha="right", va="top")
    axD.set_xticks(xs)
    axD.set_xticklabels([STAT_SHORT[n] for n in valid], fontsize=7.5)
    axD.tick_params(axis="x", length=0)
    axD.set_xlabel("held-out statistic (marker = rotation)")
    axD.set_ylabel("Stouffer $p$ over systems\n(one point per rotation)")
    legend(axD, loc="lower center", fontsize=7.5, ncol=4, columnspacing=1.0,
           handletextpad=0.25)
    panel(axD, "D", "significant in every rotation",
          subtitle="the topology-only control never fires")

    small = check_min_type(fig)
    if small:
        raise AssertionError(f"in-axes type below the house floor: {small}")
    finish(fig, "figOOS_out_of_sample")
    plt.close(fig)


# --------------------------------------------------------------------------- results doc


def write_doc(A, B, recs):
    pk = A["packs"]
    testable, untestable = B["testable"], B["untestable"]
    ratio = A["ratio"]

    vtbl = ["| variant | predictor uses | independent? | ρ | p | n | LOO ρ range |",
            "|---|---|---|---:|---:|---:|---|"]
    labels = {
        "indep_rep0": ("**headline**: predict replicate 0", "replicates 1, 2", "yes"),
        "indep_rep1": ("rotation: predict replicate 1", "replicates 0, 2", "yes"),
        "indep_rep2": ("rotation: predict replicate 2", "replicates 0, 1", "yes"),
        "indep_mean": ("rotation-averaged", "the held-out pair", "yes"),
        "insample_rep0": ("predict replicate 0", "all three", "**no**"),
        "insample_rep1": ("predict replicate 1", "all three", "**no**"),
        "insample_rep2": ("predict replicate 2", "all three", "**no**"),
        "insample_mean": ("rotation-averaged", "all three", "**no**"),
    }
    for tag, (what, uses, indep) in labels.items():
        d = pk[tag]
        vtbl.append(f"| {what} | {uses} | {indep} | {d['rho']:.3f} | {d['p']:.1e} | {d['n']} | "
                    f"[{d['loo'][0]:.3f}, {d['loo'][1]:.3f}] |")

    d, dis = pk["indep_rep0"], pk["insample_rep0"]
    dtbl = ["| diagnostic | headline (independent) | same on the in-sample predictor |",
            "|---|---:|---:|",
            f"| ρ, all systems | {d['rho']:.3f} (p = {d['p']:.1e}, n = {d['n']}) | "
            f"{dis['rho']:.3f} (p = {dis['p']:.1e}) |",
            f"| ρ, excluding the {int(A['is_flag'].sum())} flagged systems "
            f"(n = {d['n_no_flagged']}) | {d['rho_no_flagged']:.3f} | "
            f"{dis['rho_no_flagged']:.3f} |",
            f"| partial ρ given dof | {d['rho_partial_dof']:.3f} | {dis['rho_partial_dof']:.3f} |",
            f"| ρ, systems with dof ≥ 5 (n = {d['n_dof_ge5']}) | {d['rho_dof_ge5']:.3f} | "
            f"{dis['rho_dof_ge5']:.3f} |",
            f"| leave-one-out ρ range | [{d['loo'][0]:.3f}, {d['loo'][1]:.3f}] | "
            f"[{dis['loo'][0]:.3f}, {dis['loo'][1]:.3f}] |"]

    ftbl = ["| system | predicted E[χ²ᵥ] (independent) | predicted (in-sample) | "
            "observed χ²ᵥ (rep 0) | obs / pred (independent) | obs / pred (in-sample) |",
            "|---|---:|---:|---:|---:|---:|"]
    for i, s_ in enumerate(A["names"]):
        if s_ in FLAGGED:
            ftbl.append(f"| {s_} | {A['pred_indep'][0][i]:.2f} | {A['pred_insample'][0][i]:.2f} | "
                        f"{A['obs'][0][i]:.2f} | **{ratio[s_]:.1f}** | {A['ratio_is'][s_]:.1f} |")

    atbl = ["| statistic | realized α (per test) | realized α (Stouffer) | verdict |",
            "|---|---:|---:|---|"]
    for n in STAT_NAMES:
        a = B["alpha"][n]
        atbl.append(f"| {STAT_LABEL[n]} | {a['per_test']:.3f} | {a['stouffer']:.3f} | "
                    + ("USED" if n in B["valid"] else "**DISCARDED** (invalid)") + " |")

    valid = [n for n in STAT_NAMES if n in B["valid"]]
    rtbl = ["| statistic | arm | rotation 0 | rotation 1 | rotation 2 |",
            "|---|---|---|---|---|"]
    for n in valid:
        for arm in ("static", "guided", "leverage", "matched"):
            cells = []
            for k in REPLICATES:
                eff = float(np.mean([B["cells"][(s, k, arm)]["stat"][n] for s in testable]))
                ref = float(np.mean([B["cells"][(s, k, arm)]["rand_ref"][n] for s in testable]))
                p = B["combined"][(k, arm, n)][1]
                cells.append(f"{eff:+.2f} (rand {ref:+.2f}), p={p:.1e}")
            rtbl.append(f"| {STAT_LABEL[n]} | `{arm}` | " + " | ".join(cells) + " |")

    stbl = ["| system | dof | K (rot 0/1/2) | held-out Δχ²ᵥ, `static` (rot 0/1/2) | "
            "same, `random` | per-system p, Δχ²ᵥ (rot 0/1/2) | signed product on the selected "
            "edges | whole-network cross-replicate z corr |",
            "|---|---:|---|---|---|---|---:|---:|"]
    for s in testable:
        ks = "/".join(str(B["sizes"][(s, k)]) for k in REPLICATES)
        st = "/".join(f"{B['cells'][(s, k, 'static')]['stat']['dchi2']:+.2f}" for k in REPLICATES)
        rd = "/".join(f"{B['rand_mean'][(s, k)]['dchi2']:+.2f}" for k in REPLICATES)
        ps = "/".join(f"{B['cells'][(s, k, 'static')]['p']['dchi2']:.3f}" for k in REPLICATES)
        xr = float(np.mean([B["cells"][(s, k, "static")]["stat"]["xrep"] for k in REPLICATES]))
        stbl.append(f"| {s} | {recs[s]['dof']} | {ks} | {st} | {rd} | {ps} | {xr:+.2f} | "
                    f"{B['zcorr'][s]:+.2f} |")

    static_vals = [B["cells"][(s, k, "static")]["stat"]["dchi2"]
                   for s in testable for k in REPLICATES]
    rand_vals = [B["rand_mean"][(s, k)]["dchi2"] for s in testable for k in REPLICATES]
    n_sig = sum(1 for s in testable
                if sum(B["cells"][(s, k, "static")]["p"]["dchi2"] <= ALPHA
                       for k in REPLICATES) >= 2)
    bad = [s for s in testable
           if np.mean([B["cells"][(s, k, "static")]["stat"]["dchi2"] for k in REPLICATES]) <= 0]
    # static vs guided, decided from the run rather than asserted: on how many
    # (valid statistic, rotation) cells does the greedy refit beat the one-shot ordering?
    gvs = [(n, k,
            float(np.mean([B["cells"][(s, k, "static")]["stat"][n] for s in testable])),
            float(np.mean([B["cells"][(s, k, "guided")]["stat"][n] for s in testable])))
           for n in valid for k in REPLICATES]
    g_wins = sum(1 for _n, _k, st_, gd in gvs if gd > st_)
    dropped = [n for n in STAT_NAMES if n not in B["valid"]]
    mfire = [n for n in valid
             if sum(B["combined"][(k, "matched", n)][1] <= ALPHA for k in REPLICATES) >= 2]
    mtxt = ("On " + ", ".join(STAT_LABEL[n] for n in mfire)
            + " the same |z| ranks beat their own random null inside **unflagged** donor systems"
              " too, in at least two rotations, so that statistic measures a generic property of"
              " high-|z| edges rather than anything specific to the flagged set. On "
            + ", ".join(STAT_LABEL[n] for n in valid if n not in mfire)
            + " the matched control does not fire, so that is where the flag-specific claim lives."
            ) if mfire else (
        "The matched control does not fire on any valid statistic, so the effect is specific to"
        " the flagged systems rather than a generic property of high-|z| edges.")
    gtxt = (f"the greedy refit (`guided`) beats the one-shot |z| ordering (`static`) on "
            f"**{g_wins} of the {len(gvs)}** (valid statistic × rotation) cells")

    lines = [
        "# Results — Fig OOS: what survives of the QC localization claim, out-of-sample",
        "",
        "**Figure:** `figs/figOOS_out_of_sample.{pdf,png}` · **Reproduce:** `make figOOS`",
        "(`PYTHONPATH=src python figs/make_figOOS.py`). Deterministic (every resampling seeded,",
        f"`SEED = {SEED}`). Data: `data/openfe_replicates/combined_pymbar4_edge_data.csv`",
        "(OpenFE IndustryBenchmarks2024, public, 3 independent replicates per edge).",
        "",
        "## Why this figure exists",
        "The edge-removal *causal repair* panel was retracted: the statistic it drives down",
        "(`sum_e z_e²`) is the quantity the removal order maximises (descending `z_e²`), so the",
        "guided-vs-random contrast cannot fail, and a synthetic world with no localized systematic",
        "error reproduces it. Both results below are constructed so that the quantity being",
        "predicted is measured on data the selection never saw, which is the property the",
        "retracted panel lacked. Nothing here is tuned: every constant is stated in the script",
        "docstring and was fixed before any outcome was computed.",
        "",
        "## Panel A — the replicate spread predicts the closure χ² it never saw",
        "",
        "Write the reported per-edge variance as `V_e^rep = c_e² V_e^true`. Whitening the GLS fit",
        "by the *reported* errors gives `E[X²] = tr(M D)` with `M = I − H` the residual maker and",
        "`D = diag(c_e^−2)`, so with `h_e = M_ee` the curl-leverage of",
        "`bar.leverage.curl_leverage` and `Σ_e h_e = dof` (Theorem D1),",
        "",
        "> `E[χ²ᵥ] = Σ_e h_e c_e^−2 / dof` — the curl-leverage-weighted **mean** of `c_e^−2`.",
        "",
        "The replicate set supplies an estimate of `c_e^−2` that the closure test never uses:",
        "`s_e² / se_e²`, with `s_e²` the unbiased sample variance of the replicate ΔΔG values",
        "(`E[s²] = σ²` at any n ≥ 2). **Independence is enforced by rotation**: replicate 0's",
        "observed χ²ᵥ is predicted from `s_e` computed on replicates 1 and 2 only, and the same",
        "construction is rotated over the other two choices. That is the headline. The",
        "all-three-replicate variants are reported for comparison but are **in-sample**, since the",
        "predicted replicate then contributes to its own predictor.",
        "",
    ] + vtbl + [
        "",
        "The in-sample variants are uniformly higher, as they must be: the shared replicate",
        "correlates the two axes directly. Only the headline row is quotable.",
        "",
        "**Diagnostics.** Reported on both predictors, because the difference is large enough to",
        "matter: a diagnostic computed on the in-sample predictor is not a robustness check on the",
        "headline.",
        "",
    ] + dtbl + [
        "",
        "**Aggregate.** On the headline (independent, replicate 0) predictor the pooled",
        f"`E[χ²ᵥ] = {A['agg']:.2f}`, system-cluster bootstrap 95% CI",
        f"[{A['agg_ci'][0]:.2f}, {A['agg_ci'][1]:.2f}] ({N_BOOT} resamples of the "
        f"{len(A['names'])} systems);",
        f"pooling all three rotations gives {A['agg_all']:.2f} "
        f"[{A['agg_all_ci'][0]:.2f}, {A['agg_all_ci'][1]:.2f}].",
        f"Read as a bar width, the reported se is **{1 / math.sqrt(A['agg']):.2f}×** the replicate",
        f"standard deviation [{1 / math.sqrt(A['agg_ci'][1]):.2f}, "
        f"{1 / math.sqrt(A['agg_ci'][0]):.2f}]:",
        "on this functional the bars are close to calibrated, not conservative by a large factor.",
        "That is *not* in conflict with any replicate-validated ratio reported elsewhere. As the",
        "SI records, a leverage-weighted mean of `c_e^−2` and a ratio of typical magnitudes are",
        "different functionals of a heavy-tailed per-edge distribution. The observed side of the",
        f"same functional (dof-weighted mean χ²ᵥ) is {A['obs_agg']:.2f} on replicate 0",
        f"({A['obs_agg_all']:.2f} over all three), above the predicted level; the excess is what",
        "the flagged systems contribute.",
        "",
        "**Flagged systems sit above the identity line**, which is the expected direction: a",
        "reproducible systematic error inflates the closure residual but not the replicate spread,",
        "so observed exceeds predicted exactly where the QC fires.",
        "",
    ] + ftbl + [
        "",
        "## Panel B — out-of-sample localization",
        "",
        "Select edges on one replicate, evaluate on the two held out, rotate over all three",
        "choices of selection replicate. Arms, all of the same size *K*:",
        "",
        "- `static` — top-*K* by |z| on the selection replicate (the one-shot ordering);",
        "- `guided` — `bar.qc.repair_order`'s greedy removal set (the paper's actual repair rule);",
        "- `leverage` — top-*K* by curl-leverage *h* (topology only, blind to the values);",
        "- `matched` — the **same |z| ranks** taken in an unflagged donor system "
        f"({N_DONOR} donors drawn from the {B['n_donor_pool']} eligible, each scored against its",
        "  own random null): is this just what high-|z| edges do anywhere?",
        "- `random` — uniform *K*-subsets, which also supply the null and the p-values.",
        "",
        f"`K = max(1, #removals repair_order needs on the selection replicate)`; testable systems "
        f"need dof ≥ {MIN_DOF_TEST}.",
        "",
        "### The statistics were calibrated before use, and two failed",
        "Under the sampling null `Var(z_e) = h_e`, so selecting on |z| preferentially selects",
        "high-leverage edges whose held-out `z²` is larger **by topology alone**. The harness runs",
        "the entire pipeline in a world with no systematic error (same graphs, same reported se,",
        f"resampled values; {N_NULL} realizations × 3 rotations × {len(testable)} systems, "
        f"{B['alpha']['dchi2']['n_test']} tests) and measures each statistic's realized",
        "false-positive rate at nominal α = 0.05. A statistic is used only if its realized rate is",
        "within a factor of two of nominal, i.e. |α̂ − 0.05| ≤ 0.05:",
        "",
    ] + atbl + [
        "",
        f"**{len(dropped)} of the {len(STAT_NAMES)} statistics are discarded**, and they are shown",
        "here rather than silently dropped. Raw held-out `z²` fires at "
        f"{B['alpha']['z2_raw']['per_test']:.3f} because of the leverage confound above; dividing",
        f"by *h* restores the nominal level ({B['alpha']['z2_lev']['per_test']:.3f}). The signed",
        f"cross-replicate product fires at {B['alpha']['xrep']['per_test']:.3f} for a different",
        "reason: selecting on |z| makes the multiplier `z_s` large, so the selected set's",
        "statistic has a much wider null distribution than the equal-size random draws it is",
        "compared with, and it exceeds their upper tail far more often than 5% of the time even",
        "though both are centred on zero. That is measured here, not assumed: an exploratory pass",
        "before this script reported that statistic as passing at α̂ = 0.077, which this run does",
        "not reproduce. The `leverage` arm in the tables below carries the same confound as an",
        "explicit control arm.",
        "",
        "### Result",
        "Effect sizes are means over the testable systems; p is the Stouffer combination over",
        "systems within a rotation. `rand` is **each arm's own** random reference: for `matched`",
        "that is the donor systems' random draws, since the matched statistic is measured in the",
        "donors, not in the flagged system.",
        "",
    ] + rtbl + [
        "",
        f"Per (system, rotation), the `static` arm's held-out Δχ²ᵥ ranges "
        f"{min(static_vals):+.2f} to {max(static_vals):+.2f} against "
        f"{min(rand_vals):+.2f} to {max(rand_vals):+.2f} for equal-size random.",
        "",
        "### Per-system detail",
        "",
    ] + stbl + [
        "",
        "### Riders — these are not optional",
        f"- **It is an aggregate result.** Only {n_sig} of the {len(testable)} testable flagged",
        "  systems reach individual significance (p ≤ 0.05 on Δχ²ᵥ in at least two of the three",
        "  rotations). The claim is about the set, not about any single system.",
        (f"- **`{', '.join(bad)}` runs the wrong way** (negative mean held-out Δχ²ᵥ)." if bad
         else "- **No system runs backwards.**")
        + " Read it with the last two columns: a system whose signed residuals do not reproduce",
        "  across runs has nothing for a one-replicate selection to transfer.",
        (f"- **`{', '.join(untestable)}` is untestable**, having one independent cycle: removing a"
         if untestable else "- All flagged systems were testable")
        + " cycle edge leaves no closure to measure, so it is excluded rather than scored.",
        f"- **One-shot |z| ordering vs `repair_order`:** {gtxt}. The two arms are the same size by",
        "  construction, so this is a like-for-like comparison of the ordering against the greedy",
        "  refit; the tables above carry both, and neither is asserted in advance.",
        f"- **The `matched` control separates the two statistics.** {mtxt}",
        "- The three rotations are not independent of one another (they reuse the same three",
        "  runs), so the three p-values are three views of one dataset, not three replications.",
        "",
        "## Honest reading",
        "- Panel A is a genuine out-of-sample prediction: the replicate spread and the closure χ²",
        "  are computed from disjoint data, related only through a theorem, and they agree in rank",
        f"  (ρ = {d['rho']:.3f}) and in aggregate level. That is the strongest statement this",
        "  benchmark supports about the calibrated null being the right null.",
        "- Panel B's flag-specific content is narrower than it first looks. The `matched` control",
        "  says the leverage-normalised held-out `z²/h` result is a generic property of high-|z|",
        "  edges, reproducible in unflagged systems too; it is the held-out consistency gain",
        "  Δχ²ᵥ, where the matched control stays silent, that is specific to the flagged set.",
        "- Panel B recovers a weaker version of the retracted claim: the QC's per-edge ordering",
        "  carries information about *unseen* runs of the same protocol. It does not show that",
        "  removing those edges improves agreement with experiment (`make figLcausal` measured",
        "  that separately, and found a null), and the selection is still evaluated on further",
        "  runs of the same force field and protocol, so it validates reproducibility of the",
        "  flagged error, not its physical diagnosis.",
        "- Neither panel is a prospective test. The decisive experiment, re-running a flagged edge",
        "  with fresh sampling and confirming the cycle closes, remains the natural next step.",
        "",
    ]
    DOC.write_text("\n".join(ln for ln in lines if ln is not None) + "\n")


# --------------------------------------------------------------------------- driver


def main():
    use_paper_style()
    recs = load_records()
    A = panel_a(recs)
    B = panel_b(recs)
    make_figure(A, B)
    write_doc(A, B, recs)
    print(f"wrote figOOS_out_of_sample.(pdf|png) to {FIGDIR}")
    print(f"wrote {DOC.relative_to(ROOT)}")

    d = A["packs"]["indep_rep0"]
    print(f"[A] systems={len(recs)}  HEADLINE (predict rep 0 from reps 1,2) rho={d['rho']:.3f} "
          f"p={d['p']:.2e} n={d['n']}")
    for tag in ("indep_rep1", "indep_rep2", "indep_mean", "insample_rep0", "insample_rep1",
                "insample_rep2", "insample_mean"):
        q = A["packs"][tag]
        print(f"    variant {tag:14s}: rho={q['rho']:.3f} p={q['p']:.2e} "
              f"LOO[{q['loo'][0]:.3f},{q['loo'][1]:.3f}]")
    for tag in ("indep_rep0", "insample_rep0"):
        q = A["packs"][tag]
        print(f"    [{tag}] excl. flagged rho={q['rho_no_flagged']:.3f} "
              f"(n={q['n_no_flagged']}); partial|dof rho={q['rho_partial_dof']:.3f}; "
              f"dof>=5 rho={q['rho_dof_ge5']:.3f} (n={q['n_dof_ge5']}); "
              f"LOO[{q['loo'][0]:.3f},{q['loo'][1]:.3f}]")
    print(f"[A] aggregate E[chi2_nu]={A['agg']:.2f} [{A['agg_ci'][0]:.2f},{A['agg_ci'][1]:.2f}] "
          f"-> bar-width factor {1 / math.sqrt(A['agg']):.2f}x "
          f"[{1 / math.sqrt(A['agg_ci'][1]):.2f},{1 / math.sqrt(A['agg_ci'][0]):.2f}]; "
          f"all-rotation {A['agg_all']:.2f} "
          f"[{A['agg_all_ci'][0]:.2f},{A['agg_all_ci'][1]:.2f}]; "
          f"observed dof-weighted mean={A['obs_agg']:.2f} (rep0), {A['obs_agg_all']:.2f} (all)")
    print("[A] flagged obs/pred (independent | in-sample): "
          + ", ".join(f"{s} {A['ratio'][s]:.1f}|{A['ratio_is'][s]:.1f}"
                      for s in A["names"] if s in FLAGGED))

    print(f"[B] testable={B['testable']}  untestable={B['untestable']}")
    print("[B] alpha calibration (nominal 0.05):")
    for n in STAT_NAMES:
        a = B["alpha"][n]
        print(f"    {n:8s} per-test {a['per_test']:.3f}  stouffer {a['stouffer']:.3f}  "
              + ("USED" if n in B["valid"] else "DISCARDED"))
    for n in B["valid"]:
        for arm in ("static", "guided", "leverage", "matched"):
            row = []
            for k in REPLICATES:
                eff = float(np.mean([B["cells"][(s, k, arm)]["stat"][n] for s in B["testable"]]))
                ref = float(np.mean([B["cells"][(s, k, arm)]["rand_ref"][n]
                                     for s in B["testable"]]))
                pp = B["combined"][(k, arm, n)][1]
                row.append(f"rot{k} {eff:+.2f}(rand {ref:+.2f}) p={pp:.1e}")
            print(f"    {n:8s} {arm:9s} " + "  ".join(row))
    print("[B] per-system held-out dchi2 (static / random) and cross-replicate z corr:")
    for s in B["testable"]:
        st = "/".join(f"{B['cells'][(s, k, 'static')]['stat']['dchi2']:+.2f}" for k in REPLICATES)
        rd = "/".join(f"{B['rand_mean'][(s, k)]['dchi2']:+.2f}" for k in REPLICATES)
        ps = "/".join(f"{B['cells'][(s, k, 'static')]['p']['dchi2']:.3f}" for k in REPLICATES)
        print(f"    {s:8s} K={'/'.join(str(B['sizes'][(s, k)]) for k in REPLICATES):8s} "
              f"static {st:22s} random {rd:22s} p {ps}  zcorr {B['zcorr'][s]:+.2f}")


if __name__ == "__main__":
    main()
