"""Unit tests for the replicate-stability analysis (figs/make_figStab.py, Fig Stab).

Behavioural invariants of the added analysis functions:
  * Jaccard is the standard set overlap on known sets, and the closed-form random-draw
    reference matches a Monte-Carlo draw of the same sizes;
  * the pooling rule is Fig L panel C's: mean ddG and se -> se/sqrt(n) for equal bars, with
    edges below the required replicate count dropped;
  * the BH q-values reported in the results table reproduce the paper's own
    ``benjamini_hochberg`` flags exactly;
  * a synthetic network carrying a known injected systematic error is flagged in EVERY
    replicate (so the observed instability is not an artefact of the machinery), while the
    flag rate has a transition band in effect size where replicates disagree -- the mechanism
    the figure reports on the real data.
"""
from itertools import pairwise

import numpy as np
from figs.make_figStab import (
    bh_qvalues,
    edge_val,
    expected_random_jaccard,
    flag_set,
    jaccard,
    pool_replicates,
    pooled_edges,
    stability_summary,
    system_flags,
)

from bar.qc import benjamini_hochberg

# --------------------------------------------------------------------------- helpers

def _row(a, b, vals):
    """A CSV-shaped row: complex leg carries (ddg, se), solvent leg is zero, so that
    ``edge_val`` returns exactly ``(ddg, se)``. ``vals[k] is None`` = replicate k missing."""
    r = {"ligand_A": a, "ligand_B": b}
    for k, v in enumerate(vals):
        ddg, se = ("", "") if v is None else v
        r[f"complex_repeat_{k}_DG (kcal/mol)"] = ddg
        r[f"complex_repeat_{k}_dDG (kcal/mol)"] = se
        r[f"solvent_repeat_{k}_DG (kcal/mol)"] = "" if v is None else 0.0
        r[f"solvent_repeat_{k}_dDG (kcal/mol)"] = "" if v is None else 0.0
    return r


def _network(phi, se, rng, bias=None):
    """Complete graph on len(phi) nodes: y = phi_j - phi_i + N(0, se) + optional edge bias."""
    edges = []
    for i in range(len(phi)):
        for j in range(i + 1, len(phi)):
            y = phi[j] - phi[i] + rng.normal(0.0, se)
            if bias and (i, j) in bias:
                y += bias[(i, j)]
            edges.append((i, j, y, se))
    return edges


def _three_replicate_systems(n_clean, bias, seed=0, n_nodes=6, se=0.3):
    """``n_clean`` sampling-only systems plus one ('bad') carrying the SAME injected edge bias
    in all three replicates. Returns one ``{system: edges}`` dict per replicate."""
    rng = np.random.default_rng(seed)
    phis = {f"clean{i}": rng.normal(0, 2, n_nodes) for i in range(n_clean)}
    phis["bad"] = rng.normal(0, 2, n_nodes)
    per_rep = []
    for _k in range(3):
        per_rep.append({name: _network(phi, se, rng, bias=bias if name == "bad" else None)
                        for name, phi in phis.items()})
    return per_rep


# --------------------------------------------------------------------------- set overlap

def test_jaccard_on_known_sets():
    assert jaccard({"a", "b", "c"}, {"b", "c", "d"}) == 0.5      # |and|=2, |or|=4
    assert jaccard({"a"}, {"a"}) == 1.0
    assert jaccard({"a"}, {"b"}) == 0.0
    assert jaccard(set(), set()) == 1.0                          # stated convention
    assert jaccard({"a", "b"}, set()) == 0.0


def test_expected_random_jaccard_matches_monte_carlo():
    rng = np.random.default_rng(11)
    for n_a, n_b in [(6, 6), (6, 3), (3, 3)]:
        draws = []
        for _ in range(20000):
            a = set(rng.choice(48, n_a, replace=False))
            b = set(rng.choice(48, n_b, replace=False))
            draws.append(len(a & b) / len(a | b))
        assert abs(expected_random_jaccard(n_a, n_b, 48) - float(np.mean(draws))) < 0.005


def test_expected_random_jaccard_degenerate_cases():
    assert expected_random_jaccard(48, 48, 48) == 1.0            # both draws take everything
    assert expected_random_jaccard(0, 0, 48) == 1.0
    assert expected_random_jaccard(1, 1, 48) < 0.03              # near-disjoint by chance


def test_stability_summary_reports_intersection_union_and_reference():
    sets = [{"a", "b", "c"}, {"a", "b", "d"}, {"a", "e"}]
    out = stability_summary(sets, n_total=48)
    assert out["intersection"] == ["a"]
    assert out["union"] == ["a", "b", "c", "d", "e"]
    assert out["jaccard"][0] == 0.5                              # {a,b,c} vs {a,b,d}
    assert all(o > r for o, r in zip(out["jaccard"], out["jaccard_random"], strict=True))


# --------------------------------------------------------------------------- pooling rule

def test_pool_replicates_averages_and_shrinks_se():
    ddg, se = pool_replicates([(1.0, 0.2), (2.0, 0.2), (3.0, 0.2)])
    assert abs(ddg - 2.0) < 1e-12                                # replicate mean
    assert abs(se - 0.2 / np.sqrt(3)) < 1e-12                    # se -> se/sqrt(3)
    # unequal bars: root-mean-square of the se, then /sqrt(n)
    _d, se2 = pool_replicates([(0.0, 0.3), (0.0, 0.5)])
    assert abs(se2 - np.sqrt((0.09 + 0.25) / 2 / 2)) < 1e-12
    # a single replicate is a no-op on both quantities
    assert pool_replicates([(1.5, 0.4)]) == (1.5, 0.4)


def test_pooled_edges_requires_the_declared_replicate_count():
    rows = [_row("x", "y", [(1.0, 0.2), (1.0, 0.2), (1.0, 0.2)]),   # 3 replicates
            _row("y", "z", [(2.0, 0.2), (2.0, 0.2), None]),          # 2 replicates
            _row("x", "z", [(3.0, 0.2), None, None])]                # 1 replicate
    assert len(pooled_edges(rows, min_reps=1)) == 3
    assert len(pooled_edges(rows, min_reps=2)) == 2
    assert len(pooled_edges(rows, min_reps=3)) == 1
    # the surviving 2-replicate edge is pooled over the replicates it HAS (n=2, not n=3)
    kept = {(a, b): (y, se) for a, b, y, se in pooled_edges(rows, min_reps=2)}
    assert abs(kept[("y", "z")][1] - 0.2 / np.sqrt(2)) < 1e-12
    assert abs(kept[("x", "y")][1] - 0.2 / np.sqrt(3)) < 1e-12


def test_edge_val_is_the_paper_edge_construction():
    r = _row("x", "y", [(1.25, 0.4), None, None])
    r["complex_repeat_0_DG (kcal/mol)"] = 3.0
    r["solvent_repeat_0_DG (kcal/mol)"] = 1.0
    r["complex_repeat_0_dDG (kcal/mol)"] = 0.3
    r["solvent_repeat_0_dDG (kcal/mol)"] = 0.4
    ddg, se = edge_val(r, 0)
    assert abs(ddg - 2.0) < 1e-12                                # complex - solvent leg
    assert abs(se - 0.5) < 1e-12                                 # quadrature of the two legs
    assert edge_val(r, 1) is None                                # incomplete replicate dropped


# --------------------------------------------------------------------------- the test itself

def test_bh_qvalues_reproduce_the_paper_flag_function():
    rng = np.random.default_rng(5)
    for alpha in (0.01, 0.05, 0.20):
        for _ in range(20):
            p = np.concatenate([rng.uniform(0, 1, 40), rng.uniform(0, 1e-4, 5)])
            q = bh_qvalues(p)
            assert np.array_equal(q <= alpha, benjamini_hochberg(p, alpha=alpha))
    # step-up monotonicity in the p-ordering, and the [0,1] range
    p = np.array([0.001, 0.01, 0.02, 0.9])
    q = bh_qvalues(p)
    assert np.all(np.diff(q[np.argsort(p)]) >= -1e-12)
    assert q.max() <= 1.0


def test_system_flags_drops_systems_without_a_cycle():
    tree = [(0, 1, 0.0, 0.2), (1, 2, 0.0, 0.2), (2, 3, 0.0, 0.2)]   # 3 edges but no cycle
    tiny = [(0, 1, 0.0, 0.2), (1, 2, 0.0, 0.2)]                     # below MIN_EDGES
    rng = np.random.default_rng(1)
    ok = _network(rng.normal(0, 2, 5), 0.3, rng)
    rows = system_flags({"tree": tree, "tiny": tiny, "ok": ok})
    assert [d["sys"] for d in rows] == ["ok"]
    assert rows[0]["dof"] == len(ok) - 4                            # cyclomatic number


def test_injected_systematic_error_is_flagged_in_every_replicate():
    per_rep = _three_replicate_systems(n_clean=11, bias={(1, 3): 3.0}, seed=7)  # ~10 sigma
    sets = [flag_set(system_flags(reps)) for reps in per_rep]
    assert all("bad" in s for s in sets), sets
    out = stability_summary(sets, n_total=12)
    assert out["intersection"] == ["bad"]                           # survives all three runs
    assert out["union"] == ["bad"]                                  # and nothing else appears
    assert out["jaccard"] == [1.0, 1.0, 1.0]


def test_pooling_also_flags_the_injected_systematic_error():
    per_rep = _three_replicate_systems(n_clean=11, bias={(1, 3): 3.0}, seed=7)
    pooled = {}
    for name in per_rep[0]:
        vals = [[(y, se) for _a, _b, y, se in per_rep[k][name]] for k in range(3)]
        keys = [(a, b) for a, b, _y, _se in per_rep[0][name]]
        pooled[name] = [(a, b, *pool_replicates([vals[k][e] for k in range(3)]))
                        for e, (a, b) in enumerate(keys)]
    assert flag_set(system_flags(pooled)) == {"bad"}

    # ...and the pooled bars are 1/sqrt(3) of the single-replicate ones, so the reduced chi^2
    # of the systematic system rises toward the 3x variance-components ceiling (Fig L panel C).
    rc_single = next(d["rc"] for d in system_flags(per_rep[0]) if d["sys"] == "bad")
    rc_pooled = next(d["rc"] for d in system_flags(pooled) if d["sys"] == "bad")
    assert rc_pooled > 1.5 * rc_single


def test_flag_stability_has_a_transition_region_in_effect_size():
    """The observed instability is a property of thresholding a low-power statistic, not a bug in
    the machinery: as the injected systematic error grows the flag goes from never-fires to
    always-fires, and in between there is a band where the three replicates DISAGREE. No single
    bias value is asserted -- only that the three regimes exist."""
    counts = {}
    for bias in (0.0, 0.4, 0.8, 1.2, 1.5, 2.0, 3.0):
        per_rep = _three_replicate_systems(n_clean=11, bias={(1, 3): bias}, seed=3)
        counts[bias] = sum("bad" in flag_set(system_flags(reps)) for reps in per_rep)
    assert counts[0.0] == 0                                        # no error -> never flagged
    assert counts[3.0] == 3                                        # large error -> always flagged
    assert any(0 < c < 3 for c in counts.values()), counts         # a disagreement band exists
    assert all(counts[a] <= counts[b] for a, b in pairwise(sorted(counts)))  # monotone in effect
