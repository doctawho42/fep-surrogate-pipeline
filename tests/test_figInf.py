"""Unit tests for the inference added in figs/make_figInf.py (Fig Inf).

Each test is a property one of the six reported numbers rests on, and several of them are the
guard that stops the analysis being wrong in the *favourable* direction:

  * the cluster bootstrap really resamples clusters -- on data whose rows are perfectly
    correlated inside a cluster it must be about sqrt(rows-per-cluster) wider than a naive row
    bootstrap, which is the whole reason C1's interval was recomputed;
  * the within-system permutation moves rows only inside their own system, is centred on zero,
    and the resulting test has power (small p on a genuinely reproducing signal, large p on
    independent replicates);
  * the C2 contrast machinery reproduces `bar.detectors.paired_auc_bootstrap` exactly on its own
    contrast, so the extra rows are the same bootstrap and not a second, differently seeded one;
  * the C3 aggregation rules: the matched rule is `sum X^2 / sum dof`, and the three rules really
    do disagree on a skewed network set (the substantive point of that section);
  * the C4 stratum summary restricts to the stratum and exposes the empty-set Jaccard convention
    rather than hiding it;
  * the C6 reference level `exp(-1)` is the right one -- checked against a simulation, both raw
    and c4-corrected, since the whole Wade comparison is read against it;
  * the C5 circularity harness is null when `c_e == 1` despite the shared denominator, and is NOT
    null when a real per-system miscalibration is present (so its silence is informative).
"""
import math

import numpy as np
import pytest
from figs.make_figInf import (
    aggregate_rules,
    cluster_blocks,
    cluster_resample,
    fraction_below,
    null_rho,
    paired_auc_contrasts,
    pearson,
    perfect_calibration_fraction_below,
    perm_p_upper,
    stratum_summary,
    within_cluster_permutation,
)
from figs.make_figOOS import system_record

from bar.detectors import paired_auc_bootstrap

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


def _cycle_rows(n_nodes, rng, se=0.3, c=1.0):
    """Complete graph on ``n_nodes`` with three replicates; reported bar is ``c * se``."""
    phi = rng.normal(0.0, 2.0, n_nodes)
    rows = []
    for i in range(n_nodes):
        for j in range(i + 1, n_nodes):
            vals = [(phi[j] - phi[i] + rng.normal(0.0, se), c * se) for _k in REPS]
            rows.append(_row(str(i), str(j), vals))
    return rows


# --------------------------------------------------------------------------- clustering

def test_cluster_bootstrap_resamples_clusters_not_rows():
    """Rows inside a cluster carry no independent information; the cluster bootstrap must say so.

    20 clusters x 10 identical rows: an honest interval sees n = 20, a row bootstrap sees
    n = 200, so the cluster SD has to be roughly sqrt(10) times the naive one.
    """
    rng = np.random.default_rng(0)
    per, n_clust = 10, 20
    values = np.repeat(rng.normal(0.0, 1.0, n_clust), per)
    labels = np.repeat([f"s{i}" for i in range(n_clust)], per)
    names, blocks = cluster_blocks(labels)
    assert len(names) == n_clust
    assert sorted(np.concatenate(blocks).tolist()) == list(range(values.size))

    boot_rng = np.random.default_rng(1)
    clustered = np.array([values[cluster_resample(blocks, boot_rng)].mean() for _ in range(600)])
    naive = np.array([values[boot_rng.integers(0, values.size, values.size)].mean()
                      for _ in range(600)])
    ratio = clustered.std() / naive.std()
    assert 2.0 < ratio < 4.5, ratio          # sqrt(10) = 3.16


def test_within_cluster_permutation_never_leaves_its_cluster():
    labels = np.array(["a"] * 4 + ["b"] * 3 + ["c"] * 5)
    _names, blocks = cluster_blocks(labels)
    rng = np.random.default_rng(2)
    for _ in range(50):
        perm = within_cluster_permutation(blocks, labels.size, rng)
        assert sorted(perm.tolist()) == list(range(labels.size))    # a bijection
        assert list(labels[perm]) == list(labels)                   # membership preserved


def test_pearson_and_perm_p_are_the_stated_formulas():
    rng = np.random.default_rng(3)
    a, b = rng.normal(size=200), rng.normal(size=200)
    assert pearson(a, b) == pytest.approx(float(np.corrcoef(a, b)[0, 1]))
    assert math.isnan(pearson(np.ones(5), a[:5]))                   # zero variance -> undefined
    null = np.arange(100, dtype=float)
    assert perm_p_upper(1e9, null) == pytest.approx(1 / 101)
    assert perm_p_upper(-1e9, null) == pytest.approx(1.0)
    assert perm_p_upper(50.0, null) == pytest.approx(51 / 101)


def test_permutation_null_is_centred_on_zero_and_the_test_has_power():
    """Independent replicates -> a valid p (not uniformly small); a shared systematic -> it fires.

    Run over several independent data sets rather than one, since a single draw's p is itself
    random: a valid test is allowed the occasional small p and must not be asserted away.
    """
    rng = np.random.default_rng(4)
    labels = np.repeat([f"s{i}" for i in range(12)], 25)
    n = labels.size
    _names, blocks = cluster_blocks(labels)

    def run(shared_sd):
        systematic = rng.normal(0.0, shared_sd, n)
        zi = systematic + rng.normal(0.0, 1.0, n)
        zj = systematic + rng.normal(0.0, 1.0, n)
        obs = pearson(zi, zj)
        null = np.array([pearson(zi, zj[within_cluster_permutation(blocks, n, rng)])
                         for _ in range(199)])
        return obs, float(np.mean(null)), perm_p_upper(obs, null)

    quiet = [run(0.0) for _ in range(12)]
    assert max(abs(m) for _o, m, _p in quiet) < 0.06        # the null is centred on zero
    assert float(np.median([p for _o, _m, p in quiet])) > 0.1
    assert sum(p <= 0.05 for _o, _m, p in quiet) <= 3       # valid level, not conservative-by-luck

    loud = [run(1.0) for _ in range(4)]
    assert max(abs(m) for _o, m, _p in loud) < 0.06         # the null does not move with the signal
    assert all(o > 0.35 and p <= 1 / 200 + 1e-12 for o, _m, p in loud)


# --------------------------------------------------------------------------- C2 contrasts

def test_paired_auc_contrasts_reproduces_the_repo_bootstrap():
    """The published `A - max(B, C)` row must come out of the SAME resamples, bit for bit."""
    rng = np.random.default_rng(5)
    n = 40
    anchor = rng.normal(size=n)
    A = anchor > np.percentile(anchor, 82)               # a good detector
    B = rng.random(n) < 0.5                              # noise
    C = anchor > np.percentile(anchor, 60)               # a decent detector
    ours = paired_auc_contrasts(A, B, C, anchor, n_boot=300, seed=0)
    ref = paired_auc_bootstrap(A, B, C, anchor, n_boot=300, seed=0)
    for key in ("auc_a", "auc_b", "auc_c"):
        assert ours["auc"][{"auc_a": "calibrated", "auc_b": "cutoff",
                            "auc_c": "fixed_se"}[key]] == pytest.approx(ref[key])
    assert ours["vs_max"]["diff"] == pytest.approx(ref["diff"])
    assert ours["vs_max"]["ci_lo"] == pytest.approx(ref["ci_lo"])
    assert ours["vs_max"]["ci_hi"] == pytest.approx(ref["ci_hi"])
    assert ours["vs_max"]["verdict"] == ref["verdict"]
    # and the individual contrasts are the individual differences, which the max can hide
    assert ours["vs_cutoff"]["diff"] == pytest.approx(ref["auc_a"] - ref["auc_b"])
    assert ours["vs_fixed_se"]["diff"] == pytest.approx(ref["auc_a"] - ref["auc_c"])
    assert ours["vs_max"]["diff"] == pytest.approx(min(ours["vs_cutoff"]["diff"],
                                                       ours["vs_fixed_se"]["diff"]))


# --------------------------------------------------------------------------- C3 aggregation

def test_aggregate_rules_and_their_disagreement():
    obs = np.array([0.2, 0.3, 0.25, 6.0])
    dof = np.array([1.0, 1.0, 1.0, 20.0])
    r = aggregate_rules(obs, dof)
    assert r["matched"] == pytest.approx(float(np.sum(obs * dof) / np.sum(dof)))
    assert r["matched"] == pytest.approx(float(np.sum(obs * dof) / np.sum(dof)))  # = sum X2/sum nu
    assert r["median"] == pytest.approx(0.275)
    assert r["unweighted_mean"] == pytest.approx(float(np.mean(obs)))
    # the three disagree badly exactly when a few high-dof systems dominate -- C3's whole point
    assert r["matched"] > 4 * r["median"]
    # with equal dof the matched rule collapses onto the unweighted mean
    eq = aggregate_rules(obs, np.ones_like(dof))
    assert eq["matched"] == pytest.approx(eq["unweighted_mean"])


# --------------------------------------------------------------------------- C4 stratification

def test_stratum_summary_restricts_and_exposes_the_empty_convention():
    sets = [{"a", "b", "x"}, {"a", "b"}, {"a", "x"}]
    swings = {"a": 2.0, "b": 4.0, "x": 10.0, "q": 1.5}
    hi = stratum_summary(sets, {"a", "b", "q"}, swings, "hi")
    assert hi["n"] == 3 and hi["sizes"] == [2, 2, 1]
    assert hi["ever"] == ["a", "b"] and hi["inter"] == ["a"]
    assert hi["jaccard"] == pytest.approx([1.0, 0.5, 0.5])
    assert hi["median_swing"] == pytest.approx(2.0)               # median of 2.0, 4.0, 1.5

    empty = stratum_summary(sets, {"q"}, swings, "none flagged")
    assert empty["sizes"] == [0, 0, 0]
    assert empty["mean_jaccard"] == pytest.approx(1.0)            # the empty-set convention...
    assert empty["ever"] == [] and empty["inter"] == []           # ...with nothing behind it


# --------------------------------------------------------------------------- C6 reference

def test_perfect_calibration_reference_matches_a_simulation():
    """`exp(-1)` really is the fraction of edges a correct bar puts below a 3-replicate SD."""
    assert perfect_calibration_fraction_below(3) == pytest.approx(math.exp(-1.0), rel=1e-9)
    c4 = math.sqrt(2.0 / 2.0) * math.gamma(3 / 2) / math.gamma(2 / 2)
    assert perfect_calibration_fraction_below(3, c4) == pytest.approx(math.exp(-c4 ** 2), rel=1e-9)

    rng = np.random.default_rng(6)
    sigma = np.exp(rng.normal(0.0, 0.8, 60000))          # heterogeneous, correctly reported bars
    draws = rng.normal(0.0, 1.0, (3, sigma.size)) * sigma
    s = draws.std(axis=0, ddof=1)
    assert fraction_below(sigma, s) == pytest.approx(math.exp(-1.0), abs=0.01)
    assert fraction_below(sigma, s, c4) == pytest.approx(math.exp(-c4 ** 2), abs=0.01)
    # 5 replicates has a different reference, so the number is n-specific and not a constant
    assert perfect_calibration_fraction_below(5) != pytest.approx(math.exp(-1.0), abs=1e-3)


# --------------------------------------------------------------------------- C5 circularity

def _synthetic_recs(n_sys=8, seed=7, c_true=None):
    rng = np.random.default_rng(seed)
    return {f"s{i}": system_record(_cycle_rows(5, rng, se=0.3,
                                               c=1.0 if c_true is None else c_true[i]))
            for i in range(n_sys)}


def test_circularity_harness_is_null_but_can_see_a_real_signal():
    """The shared denominator alone does not manufacture rank correlation; a real one shows."""
    recs = _synthetic_recs()
    rng = np.random.default_rng(8)
    null = np.array([null_rho(recs, rng) for _ in range(40)])
    assert abs(float(np.mean(null))) < 0.25, float(np.mean(null))

    names = sorted(recs)
    cmap = dict(zip(names, np.exp(np.linspace(math.log(0.4), math.log(2.5), len(names))),
                    strict=True))
    rng_pc = np.random.default_rng(9)
    ctrl = np.array([null_rho(recs, rng_pc, cmap) for _ in range(20)])
    assert float(np.median(ctrl)) > 0.6, float(np.median(ctrl))
