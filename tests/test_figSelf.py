"""Unit tests for the two-sided self-calibration analysis (figs/make_figSelf.py, Fig Self).

Behavioural invariants of the added analysis, each one a property a claim in
`docs/results_figSelf.md` rests on:

  * the ratio table is `make_figA_replicates.load()`'s, edge for edge -- no second ratio
    definition was invented (D1);
  * the direction convention is the stated one: ``ratio > 1`` tightens a bar, ``ratio < 1``
    loosens it, and the one-sided mode keeps only the loosening branch, i.e. P8's rule (D3);
  * an unmeasurable ratio leaves the bar exactly as reported, so it can neither add nor remove
    a flag (D2), and the cap clips both tails (D4);
  * the nominal arm reproduces Fig L's published population, median reduced chi^2 and flag set,
    which is the script's own STOP condition;
  * a neutral lookup reproduces the nominal arm exactly, so the self-calibration machinery is
    a no-op when it measures nothing;
  * the mechanism the figure claims is real on synthetic data with a KNOWN answer: a system
    carrying a genuine systematic error but reported with over-wide bars is missed by the
    nominal test and recovered once the bars are corrected to their measured width;
  * the synthetic-null harness reproduces the ``nu/(nu-2)`` ordering it is quoted for -- the
    per-edge (n=3) arm inflates its own false-positive rate above the aggregate arm's, and the
    nominal arm stays near alpha.
"""
import math

import numpy as np
import pytest
from figs.make_figSelf import (
    ALPHA,
    CAP,
    FLAGGED_ANCHOR,
    MEDIAN_RC_ANCHOR,
    _neutral,
    analyze,
    corrected_se,
    edge_ratios,
    flag_set,
    median_rc,
    null_calibration,
    per_edge_lookup,
    per_system_lookup,
    run_test,
    system_ratios,
    transitions,
)

# --------------------------------------------------------------------------- helpers

def _row(a, b, vals, system="s"):
    """A CSV-shaped row: the complex leg carries ``(ddg, se)`` and the solvent leg is zero, so
    ``make_figL.edge_val`` returns exactly ``(ddg, se)`` for replicate ``k``."""
    r = {"system name": system, "ligand_A": a, "ligand_B": b}
    for k, v in enumerate(vals):
        ddg, se = ("", "") if v is None else v
        r[f"complex_repeat_{k}_DG (kcal/mol)"] = ddg
        r[f"complex_repeat_{k}_dDG (kcal/mol)"] = se
        r[f"solvent_repeat_{k}_DG (kcal/mol)"] = "" if v is None else 0.0
        r[f"solvent_repeat_{k}_dDG (kcal/mol)"] = "" if v is None else 0.0
    return r


def _system(name, phi, sigma, reported, rng, bias=None):
    """A complete graph on ``len(phi)`` nodes, 3 replicates. The TRUE per-replicate noise is
    ``sigma``; the REPORTED bar is ``reported`` (so the measured ratio is ``reported/sigma``).
    ``bias`` maps an edge to a systematic offset present in every replicate."""
    rows = []
    for i in range(len(phi)):
        for j in range(i + 1, len(phi)):
            off = (bias or {}).get((i, j), 0.0)
            vals = [(phi[j] - phi[i] + off + rng.normal(0.0, sigma), reported) for _ in range(3)]
            rows.append(_row(f"L{i}", f"L{j}", vals, system=name))
    return rows


def _world(n_clean, bias, seed=0, n_nodes=6, sigma=0.10, reported=0.25):
    rng = np.random.default_rng(seed)
    by = {}
    for i in range(n_clean):
        by[f"clean{i}"] = _system(f"clean{i}", rng.normal(0, 2, n_nodes), sigma, reported, rng)
    by["bad"] = _system("bad", rng.normal(0, 2, n_nodes), sigma, reported, rng, bias=bias)
    return by


def _ratios(by):
    """Per-edge ratio table for a synthetic ``by``, using the module's own formula."""
    from figs.make_figL import edge_val

    table = {}
    for name, rows in by.items():
        for r in rows:
            vals = [edge_val(r, k) for k in (0, 1, 2)]
            rep = float(np.sqrt(np.mean([v[1] ** 2 for v in vals])))
            repl = float(np.std([v[0] for v in vals], ddof=1))
            table[(name, r["ligand_A"], r["ligand_B"])] = {
                "rep": rep, "repl": repl, "ratio": rep / repl if repl > 1e-6 else None}
    return table


# --------------------------------------------------------------------------- D1: the ratio table

def test_ratio_table_is_figA_replicates_edge_for_edge():
    """The per-edge (rep, repl) pair is `make_figA_replicates.load()`'s, not a re-derivation."""
    from figs.make_figA_replicates import load

    mine = sorted((round(d["rep"], 12), round(d["repl"], 12))
                  for d in edge_ratios().values())
    theirs = sorted((round(r[0], 12), round(r[1], 12)) for r in load())
    assert mine == theirs


def test_system_ratio_is_the_pooled_rms_formula():
    """The aggregate ratio is Fig A-rep's / P8's pooled RMS, over non-degenerate edges only."""
    table = {("s", "a", "b"): {"rep": 2.0, "repl": 1.0, "ratio": 2.0},
             ("s", "b", "c"): {"rep": 1.0, "repl": 2.0, "ratio": 0.5},
             ("s", "c", "d"): {"rep": 9.0, "repl": 0.0, "ratio": None}}   # degenerate: excluded
    got = system_ratios(table)["s"]
    assert got == pytest.approx(math.sqrt((4 + 1) / 2) / math.sqrt((1 + 4) / 2))


# --------------------------------------------------------------------------- D2/D3/D4: the rule

def test_direction_convention_is_the_stated_one():
    assert corrected_se(1.0, 2.0) == pytest.approx(0.5)     # ratio > 1 -> bar TIGHTER
    assert corrected_se(1.0, 0.5) == pytest.approx(2.0)     # ratio < 1 -> bar LOOSER
    assert corrected_se(1.0, 1.0) == pytest.approx(1.0)


def test_one_sided_mode_is_p8s_direction():
    """P8's rule: loosen where the bar was too tight, never tighten."""
    assert corrected_se(1.0, 0.5, mode="one-sided") == pytest.approx(2.0)
    assert corrected_se(1.0, 2.0, mode="one-sided") == 1.0
    assert corrected_se(1.0, 1.0, mode="one-sided") == 1.0


@pytest.mark.parametrize("bad", [None, 0.0, -1.0, math.nan, math.inf])
def test_unmeasurable_ratio_leaves_the_bar_untouched(bad):
    assert corrected_se(0.37, bad) == 0.37


def test_degenerate_replicate_spread_is_marked_unmeasurable():
    from figs.make_figSelf import _edge_ratio
    _rep, _repl, ratio = _edge_ratio([0.2, 0.2, 0.2], [1.0, 1.0, 1.0])   # zero spread
    assert ratio is None


def test_cap_clips_both_tails():
    assert corrected_se(1.0, 100.0, cap=CAP) == pytest.approx(1.0 / CAP)
    assert corrected_se(1.0, 0.001, cap=CAP) == pytest.approx(CAP)
    assert corrected_se(1.0, 2.0, cap=CAP) == pytest.approx(0.5)         # inside: untouched


def test_a_missing_edge_key_is_left_unchanged():
    look = per_edge_lookup({("s", "a", "b"): {"ratio": 4.0}})
    assert look("s", "a", "b") == 4.0
    assert look("s", "x", "y") is None
    assert corrected_se(1.0, look("s", "x", "y")) == 1.0


# --------------------------------------------------------------------------- the pipeline

def test_neutral_lookup_reproduces_the_nominal_arm_exactly():
    by = _world(n_clean=3, bias={(0, 1): 0.4})
    a = run_test(by, _neutral)
    b = run_test(by, per_edge_lookup({}))          # empty table: nothing measurable anywhere
    assert [d["rc"] for d in a] == [d["rc"] for d in b]
    assert flag_set(a) == flag_set(b)


def test_bh_q_reproduces_the_flags():
    by = _world(n_clean=5, bias={(0, 1): 0.6, (2, 3): -0.5})
    rows = run_test(by, per_edge_lookup(_ratios(by)))
    assert all(d["flag"] == (d["q"] <= ALPHA) for d in rows)


def test_scaling_a_whole_network_scales_chi2_by_the_square():
    by = _world(n_clean=1, bias={})
    base = run_test(by, _neutral)
    scaled = run_test(by, per_system_lookup({d["sys"]: 2.0 for d in base}))
    for a, b in zip(base, scaled, strict=True):
        assert b["rc"] == pytest.approx(4.0 * a["rc"])


def test_transitions_report_kept_gained_lost():
    nom = [{"sys": "a", "flag": True}, {"sys": "b", "flag": True}, {"sys": "c", "flag": False}]
    new = [{"sys": "a", "flag": True}, {"sys": "b", "flag": False}, {"sys": "c", "flag": True}]
    assert transitions(nom, new) == (["a"], ["c"], ["b"])


# --------------------------------------------------------------------------- the mechanism

def test_selfcalibration_recovers_a_detection_the_wide_null_suppresses():
    """The figure's claim, on synthetic data with a known answer: a system that really does carry
    a systematic error, but whose bars are reported ~2.5x too wide, is MISSED by the nominal test
    and recovered once the bars are corrected to the width its own replicates measure."""
    by = _world(n_clean=8, bias={(0, 1): 0.45, (2, 3): -0.40, (1, 4): 0.35}, seed=7)
    nom = run_test(by, _neutral)
    agg = run_test(by, per_system_lookup(system_ratios(_ratios(by))))
    assert "bad" not in flag_set(nom)                      # suppressed by the over-wide null
    assert "bad" in flag_set(agg)                          # recovered by self-calibration
    assert flag_set(nom) <= flag_set(agg)                  # nothing is lost on the way


def test_selfcalibration_does_not_manufacture_a_detection_from_nothing():
    """A clean world (no injected error anywhere, bars reported 2.5x too wide) must not turn into
    a wall of flags: the correction restores chi^2 to ~1, it does not push it past the threshold
    for every system."""
    by = _world(n_clean=9, bias={}, seed=11)
    agg = run_test(by, per_system_lookup(system_ratios(_ratios(by))))
    assert median_rc(agg) == pytest.approx(1.0, abs=0.45)
    assert len(flag_set(agg)) <= 1


# --------------------------------------------------------------------------- the null harness

def test_null_harness_reproduces_the_nu_over_nu_minus_two_ordering():
    """The quoted diagnostic: on a perfectly-calibrated null the nominal arm stays near alpha,
    the per-edge (n=3, nu=2) arm inflates its own false-positive rate above the aggregate arm's,
    and the held-out (n=2, nu=1) arm is worst of all."""
    by = _world(n_clean=9, bias={}, seed=3)
    out = null_calibration(by, ["nominal", "aggregate", "per-edge", "held-out"], n_draws=40,
                           seed=1234)
    assert out["nominal"]["p_any"] <= 0.25
    assert out["per-edge"]["mean_flagged"] > out["aggregate"]["mean_flagged"]
    assert out["held-out"]["mean_flagged"] > out["per-edge"]["mean_flagged"]
    assert out["nominal"]["median_rc"] == pytest.approx(1.0, abs=0.3)


def test_null_harness_is_deterministic():
    by = _world(n_clean=3, bias={}, seed=5)
    a = null_calibration(by, ["nominal", "aggregate"], n_draws=8, seed=99)
    b = null_calibration(by, ["nominal", "aggregate"], n_draws=8, seed=99)
    assert a == b


# --------------------------------------------------------------------------- the real-data anchor

def test_real_data_anchor_and_headline():
    """Two things at once, on one pass over the real data (the pass is the expensive part):

    * the script's own STOP condition -- the nominal arm must BE Fig L (48 systems, median
      reduced chi^2 0.34, the published six), otherwise nothing downstream may be trusted;
    * the load-bearing half of the headline -- the two-sided correction removes no published
      flag (an anti-conservative arm cannot lose a flag by accident).
    """
    res = analyze(n_null=2)
    nominal = res["nominal"]
    assert len(nominal) == 48
    assert median_rc(nominal) == pytest.approx(MEDIAN_RC_ANCHOR, abs=0.01)
    assert flag_set(nominal) == set(FLAGGED_ANCHOR)

    _kept, _gained, lost = transitions(nominal, res["aggregate"])
    assert lost == []
    assert flag_set(nominal) <= flag_set(res["aggregate"])
