"""Unit tests for the out-of-sample QC checks (figs/make_figOOS.py, Fig OOS).

Behavioural invariants of the added analysis, each one a property the figure's claims rest on:
  * the panel-A identity: the curl-leverage-weighted mean of ``c_e^-2`` really is ``tr(M D)/dof``
    for the residual maker ``M`` of the whitened GLS fit, checked against an explicit projector;
  * the panel-A estimator is unbiased in the direction that matters: on a network whose reported
    bars are inflated by a known factor ``c``, the replicate-spread prediction recovers ``c^-2``
    and tracks the observed reduced chi^2;
  * bridges never enter a selection (they have ``h = 0`` and ``z = 0``, so ``z^2/h`` is 0/0);
  * the rank matcher returns exactly the requested rank multiset, which is what makes the
    ``matched`` control |z|-comparable rather than merely equal-size;
  * a synthetic network carrying a known injected localized error is recovered OUT OF SAMPLE:
    selecting on one replicate improves closure on the two held out, and a no-error twin does not;
  * the alpha-calibration harness catches a deliberately mis-calibrated statistic (one that reads
    the selection replicate it was selected on), which is the guard that disqualified raw
    held-out ``z^2`` in the figure itself.
"""
import numpy as np
import pytest
from figs.make_figOOS import (
    alpha_calibration,
    curl_leverage,
    leverage_weighted_mean,
    perm_p,
    pool_ranks,
    predicted_chi2,
    prepare,
    random_null,
    rank_matched,
    sample_var,
    statistics,
    stouffer,
    system_record,
)

from bar.qc import _incidence, gls_network

REPS = (0, 1, 2)


# --------------------------------------------------------------------------- builders

def _row(a, b, vals):
    """A CSV-shaped row in the paper's edge construction: the complex leg carries ``(ddg, se)``
    and the solvent leg is zero, so ``edge_val`` returns exactly ``(ddg, se)``."""
    r = {"ligand_A": a, "ligand_B": b}
    for k, (ddg, se) in enumerate(vals):
        r[f"complex_repeat_{k}_DG (kcal/mol)"] = ddg
        r[f"complex_repeat_{k}_dDG (kcal/mol)"] = se
        r[f"solvent_repeat_{k}_DG (kcal/mol)"] = 0.0
        r[f"solvent_repeat_{k}_dDG (kcal/mol)"] = 0.0
    return r


def _cycle_rows(n_nodes, rng, se=0.3, bias=None, c=1.0, extra=()):
    """Rows for a complete graph on ``n_nodes`` with three independent replicates.

    ``y = phi_b - phi_a + N(0, se) [+ bias]``; the REPORTED bar is ``c * se``, so ``c`` is the
    known miscalibration factor of the panel-A error model. ``bias`` maps an edge ``(i, j)`` to a
    systematic offset present in every replicate. ``extra`` adds pendant (bridge) edges.
    """
    phi = rng.normal(0.0, 2.0, n_nodes)
    rows = []
    for i in range(n_nodes):
        for j in range(i + 1, n_nodes):
            vals = []
            for _k in REPS:
                y = phi[j] - phi[i] + rng.normal(0.0, se) + (bias or {}).get((i, j), 0.0)
                vals.append((y, c * se))
            rows.append(_row(str(i), str(j), vals))
    for a, b in extra:
        rows.append(_row(str(a), str(b), [(rng.normal(0.0, se), c * se) for _k in REPS]))
    return rows


# --------------------------------------------------------------------------- panel A identity

def test_leverage_weighted_mean_is_tr_MD_on_a_small_case():
    """``sum_e h_e d_e == tr(M D)`` with ``M = I - H`` built explicitly from the whitened
    incidence, and ``sum_e h_e == dof`` (so the weighted mean is exactly ``E[chi2_nu]``)."""
    edges = [("a", "b", 0.4, 0.2), ("b", "c", -0.3, 0.5), ("a", "c", 0.2, 0.3),
             ("c", "d", 1.1, 0.25), ("b", "d", 0.7, 0.4)]
    _nodes, B, _y, V = _incidence(edges)
    xt = B / np.sqrt(V)[:, None]
    hat = xt @ np.linalg.pinv(xt)
    resid_maker = np.eye(len(edges)) - hat

    h = curl_leverage(edges)
    assert np.allclose(h, np.diag(resid_maker), atol=1e-9)
    dof = gls_network(edges).dof
    assert abs(float(np.sum(h)) - dof) < 1e-9

    for d in ([1.0] * 5, [0.5, 2.0, 1.0, 0.25, 4.0], [3.0, 0.1, 0.1, 0.1, 0.1]):
        d = np.asarray(d)
        assert abs(float(np.sum(h * d)) - float(np.trace(resid_maker @ np.diag(d)))) < 1e-9
        # ...hence the leverage-weighted mean IS the predicted reduced chi^2
        assert abs(leverage_weighted_mean(h, d)
                   - float(np.trace(resid_maker @ np.diag(d))) / dof) < 1e-9


def test_sample_var_is_the_unbiased_replicate_variance():
    a = [np.array([1.0, 0.0]), np.array([3.0, 0.0]), np.array([5.0, 6.0])]
    assert np.allclose(sample_var(a[:2]), [(1.0 - 3.0) ** 2 / 2, 0.0])   # n=2 -> (x-y)^2/2
    assert np.allclose(sample_var(a), np.var(np.vstack(a), axis=0, ddof=1))


def test_predicted_chi2_recovers_a_known_bar_miscalibration():
    """Reported bars inflated by a known ``c`` -> predicted ``E[chi2_nu] ~ c^-2``, and the
    prediction tracks the observed closure chi^2 of the replicate it never saw."""
    for c in (1.0, 1.5):
        preds, obss = [], []
        for seed in range(30):
            rng = np.random.default_rng(1000 + seed)
            rec = system_record(_cycle_rows(6, rng, se=0.3, c=c))
            preds.append(predicted_chi2(rec, 0, (1, 2)))
            obss.append(rec["rc"][0])
        assert abs(float(np.mean(preds)) - c ** -2) < 0.15 * c ** -2
        assert abs(float(np.mean(obss)) - c ** -2) < 0.25 * c ** -2
        # the predictor is genuinely independent of its target: it uses replicates 1 and 2 only
        rng = np.random.default_rng(7)
        rec = system_record(_cycle_rows(6, rng, se=0.3, c=c))
        moved = dict(rec, y={**rec["y"], 0: rec["y"][0] + 5.0})
        assert predicted_chi2(moved, 0, (1, 2)) == pytest.approx(predicted_chi2(rec, 0, (1, 2)))


# --------------------------------------------------------------------------- selection machinery

def test_bridges_are_excluded_from_the_selection_pool():
    rng = np.random.default_rng(3)
    rec = system_record(_cycle_rows(4, rng, extra=[(0, "x"), ("x", "w")]))
    assert rec["E"] == 8 and rec["pool"].size == 6          # 6 cycle edges + 2 pendant bridges
    h = rec["h"][0]
    assert np.all(h[rec["pool"]] > 1e-9)
    off = [e for e in range(rec["E"]) if e not in set(rec["pool"].tolist())]
    assert np.allclose(h[off], 0.0, atol=1e-9)
    assert np.allclose(rec["z"][0][off], 0.0, atol=1e-9)    # a bridge is fitted exactly


def test_rank_matched_selects_the_requested_rank_multiset():
    scores = np.array([0.1, -5.0, 2.0, -3.0, 0.5])          # |.| ranks: 5 > 3 > 2 > 4 > 1
    assert list(rank_matched(scores, [1])) == [1]
    assert list(rank_matched(scores, [1, 2, 3])) == [1, 3, 2]
    assert list(rank_matched(scores, [1, 3, 5])) == [1, 2, 0]
    assert list(rank_matched(scores, [5, 4])) == [0, 4]      # order follows the ranks asked for
    assert sorted(rank_matched(scores, range(1, 6))) == [0, 1, 2, 3, 4]
    assert list(rank_matched([1.0, 1.0, 1.0], [1, 2])) == [0, 1]   # ties break by index
    for bad in ([0], [6], [1, 9]):
        with pytest.raises(ValueError):
            rank_matched(scores, bad)


def test_pool_ranks_matches_inside_the_pool_only():
    rng = np.random.default_rng(4)
    rec = system_record(_cycle_rows(4, rng, extra=[(0, "x"), ("x", "w")]))
    sel = pool_ranks(rec, rec["z"][0], [1, 2])
    assert set(sel.tolist()) <= set(rec["pool"].tolist())
    # the same ranks computed by hand on the pool subset
    pool = rec["pool"]
    order = pool[np.argsort(-np.abs(rec["z"][0][pool]), kind="stable")]
    assert list(sel) == list(order[:2])


# --------------------------------------------------------------------------- out-of-sample power

def _oos_test(bias, seed, k_size=2, n_draws=199):
    """Select the top-K |z| edges on replicate 0, score the held-out replicates against a
    uniform-random null of the same size. Returns (statistics, p-values, selection)."""
    rng = np.random.default_rng(seed)
    rec = system_record(_cycle_rows(6, rng, se=0.3, bias=bias))
    ctx = prepare(rec)
    sel = pool_ranks(ctx, ctx["z"][0], np.arange(1, k_size + 1))
    obs = statistics(ctx, 0, sel)
    null = random_null(ctx, 0, k_size, n_draws, np.random.default_rng(seed + 1))
    return obs, {n: perm_p(obs[n], null[n]) for n in obs}, sel, rec


def test_injected_localized_error_is_recovered_out_of_sample():
    """A systematic offset on two edges, present in all three replicates: selecting on replicate 0
    localizes it, and removing the selection improves closure on the replicates never seen."""
    bias = {(1, 3): 2.5, (0, 4): -2.5}
    obs, p, sel, rec = _oos_test(bias, seed=11)
    # the two injected edges are the two selected ones
    injected = {i for i, (a, b, _y, _se) in enumerate(rec["edges"][0])
                if (int(a), int(b)) in bias}
    assert set(sel.tolist()) == injected
    assert obs["dchi2"] > 0 and p["dchi2"] <= 0.05          # held-out closure genuinely improves
    assert obs["xrep"] > 0                                   # the same signed error, both runs
    assert obs["z2_lev"] > 1.0                               # and it is large where predicted


def test_a_no_error_twin_shows_no_out_of_sample_gain():
    """The same construction with no injected error: the held-out gain is not significant.

    Run over several seeds and require the p-values to behave like a valid test (not all small),
    which is the property the alpha harness measures at scale in the figure.
    """
    ps = [_oos_test(None, seed=s)[1]["dchi2"] for s in range(20, 32)]
    assert sum(x <= 0.05 for x in ps) <= 3, ps
    assert float(np.median(ps)) > 0.1, ps


# --------------------------------------------------------------------------- alpha calibration

def _bogus_stat(ctx, sel_rep, sel):
    """A deliberately mis-calibrated statistic: it reads the SELECTION replicate, so the |z|-top-K
    arm maximises it by construction and every random draw loses. Realized alpha must be ~1."""
    out = statistics(ctx, sel_rep, sel)
    out["bogus"] = float(np.mean(np.abs(ctx["z"][sel_rep][np.asarray(sel, dtype=int)])))
    return out


def test_alpha_harness_catches_a_miscalibrated_statistic():
    rng = np.random.default_rng(5)
    recs = {f"s{i}": system_record(_cycle_rows(5, rng, se=0.3)) for i in range(3)}
    sizes = {(s, k): 2 for s in recs for k in REPS}
    alpha = alpha_calibration(recs, list(recs), sizes, np.random.default_rng(6),
                              n_null=12, n_draws=99, stat_names=("dchi2", "bogus"),
                              stat_fn=_bogus_stat)
    assert alpha["bogus"]["per_test"] > 0.9                  # fires on essentially every null test
    assert alpha["bogus"]["stouffer"] > 0.9
    assert alpha["dchi2"]["per_test"] < 0.25                 # the honest statistic stays near 0.05
    assert alpha["dchi2"]["n_test"] == 12 * len(REPS) * len(recs)


def test_perm_p_and_stouffer_are_the_stated_formulas():
    null = np.arange(100, dtype=float)
    assert perm_p(1e9, null) == pytest.approx(1 / 101)       # nothing in the null beats it
    assert perm_p(-1e9, null) == pytest.approx(1.0)          # everything does
    assert perm_p(50.0, null) == pytest.approx(51 / 101)
    z_half, p_half = stouffer([0.5, 0.5, 0.5])
    assert z_half == pytest.approx(0.0) and p_half == pytest.approx(0.5)
    z_lo, p_lo = stouffer([0.01, 0.01])
    assert z_lo > 0 and p_lo < 0.01                          # combining evidence sharpens it
    assert stouffer([0.9, 0.9])[1] > 0.9
