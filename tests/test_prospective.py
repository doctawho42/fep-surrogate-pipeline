"""Unit tests for the cage prospective-loop decision scaffold (src/screen/prospective.py)."""
from __future__ import annotations

import inspect
import math
import pathlib

import pytest

from screen import prospective as P

# external witness; independent literal, never edited
PINNED_SHA256 = "66ef0bc848ab0b515b2e5b32be44230b0f528504b0b018caccc6884a0c73c488"


def test_decision_lcb_commit_and_boundary():
    r = P.decision_lcb(effect=3.0, sigma=0.5, z=1.645, tau=1.0)
    assert math.isclose(r["lcb"], 3.0 - 1.645 * 0.5)
    assert r["commit"] is True
    # boundary: lcb exactly == tau commits; a hair above tau does not
    assert P.decision_lcb(2.0, 1.0, 1.0, 1.0)["commit"] is True
    assert P.decision_lcb(2.0, 1.0, 1.0, 1.01)["commit"] is False


def test_enantiopreference_sqrt_combination_is_pinned():
    # sqrt(0.3**2 + 0.4**2) = 0.5 -> CI = 1.0 +/- 1.645*0.5 EXCLUDES 0 (discordant).
    # a wrong linear sum 0.3+0.4=0.7 -> CI = 1.0 +/- 1.645*0.7 would COVER 0. This case
    # discriminates the correct combination from the common wrong one.
    r = P.enantiopreference(rr=1.0, ss=0.0, sigma_rr=0.3, sigma_ss=0.4, z=1.645)
    assert math.isclose(r["delta"], 1.0)
    assert math.isclose(r["ci"][1] - r["ci"][0], 2 * 1.645 * 0.5)  # half-width from sigma_delta=0.5
    assert r["discordant"] is True
    # null pair: delta small vs sigma -> CI covers 0 -> not discordant
    assert P.enantiopreference(1.0, 0.9, 0.5, 0.5, 1.645)["discordant"] is False


def test_aggregation_guard_is_sigma_aware():
    # survives: detergent LCB (0.9 - 1.645*0.05 = 0.818) >= 0.5*1.0
    assert P.aggregation_guard(1.0, 0.9, 0.05, 0.5, 1.645)["survives"] is True
    # artifact: detergent abolishes signal (0.2 - 0.082 = 0.118) < 0.5
    assert P.aggregation_guard(1.0, 0.2, 0.05, 0.5, 1.645)["artifact"] is True
    # sigma matters: point 0.55 > 0.5 but LCB (0.55 - 0.1645 = 0.3855) < 0.5 -> artifact
    assert P.aggregation_guard(1.0, 0.55, 0.1, 0.5, 1.645)["artifact"] is True


def test_stop_rule_is_measured_sigma_only():
    # signature exposes sigma_assay, NOT a confidence argument
    params = list(inspect.signature(P.stop_rule).parameters)
    assert params == ["effect", "sigma_assay", "tau", "bound"]
    assert "conf" not in params and "claimed_conf" not in params
    # strong effect vs small sigma -> high confidence -> stop
    assert P.stop_rule(3.0, 0.5, 1.0, 0.9)["stop"] is True
    assert P.stop_rule(1.2, 1.0, 1.0, 0.9)["stop"] is False
    # conf is monotone increasing as sigma shrinks (effect > tau)
    assert P.stop_rule(2.0, 0.5, 1.0, 0.9)["conf"] > P.stop_rule(2.0, 1.0, 1.0, 0.9)["conf"]
    with pytest.raises(ValueError):
        P.stop_rule(2.0, 0.0, 1.0, 0.9)


def test_prereg_loads_and_matches_pinned_anchor():
    # the module constant and the test's independent witness must agree with the file on disk
    assert P.PREREG_SHA256 == PINNED_SHA256
    assert P.sha256_of(P.PREREG_PATH) == PINNED_SHA256
    pr = P.load_prereg()
    assert pr.version == 1
    assert [s["id"] for s in pr.species] == ["RR-OAc", "SS-OAc", "RR-OH", "SS-OH"]
    assert set(pr.forecast) == {"F1", "F2", "F3", "F4", "F5"}
    assert pr.decision["z"] == 1.645
    assert pr.decision["frac"] == 0.5
    assert pr.decision["stop_bound"] == 0.90


def test_prereg_immutable_against_external_anchor(tmp_path):
    # mutate the frozen file's content -> its hash changes -> load must raise vs the external anchor
    tampered = tmp_path / "tampered.yaml"
    tampered.write_bytes(pathlib.Path(P.PREREG_PATH).read_bytes() + b"\n# post-hoc edit\n")
    with pytest.raises(ValueError):
        P.load_prereg(str(tampered), expected_sha256=PINNED_SHA256)


def test_prereg_f3_not_scored_and_power_justified():
    pr = P.load_prereg()
    assert pr.forecast["F3"]["scored"] is False
    mr = pr.decision["min_replicates"]
    assert mr["value"] == 3
    assert mr["power_justification"].strip()  # non-empty derivation present


def test_disposition_table_is_exhaustive():
    pr = P.load_prereg()
    dispositions = {row["disposition"] for row in pr.disposition_table}
    assert len(pr.disposition_table) == 6
    assert "HIT" in dispositions and "F5-confirmed-promiscuous" in dispositions
    assert "ambiguous-engagement" in dispositions and "inconclusive-but-suggestive" in dispositions


def test_score_forecast_all_negative_outcomes():
    pr = P.load_prereg()
    obs = {"specific_engagement": False, "nr_signal": None,
           "steroid_rescue": False, "detergent_surviving": False, "enantiodiscordant": False}
    rep = P.score_forecast(pr, obs)
    by = {o.fid: o.outcome for o in rep.outcomes}
    assert by == {"F1": "confirmed", "F2": "confirmed", "F4": "not_testable", "F5": "confirmed"}
    assert "F3" not in by  # F3 is context, never scored
    assert rep.scorecard == {"confirmed": 3, "refuted": 0, "not_testable": 1}


def test_score_forecast_all_refuted():
    pr = P.load_prereg()
    obs = {"specific_engagement": True,
           "nr_signal": {"target": "ER", "agonist_comparable": True},
           "steroid_rescue": True, "enantiomer_call": {"GR": "RR-OAc"},  # prereg says GR->SS-OAc
           "detergent_surviving": True, "enantiodiscordant": True}
    rep = P.score_forecast(pr, obs)
    by = {o.fid: o.outcome for o in rep.outcomes}
    assert by == {"F1": "refuted", "F2": "refuted", "F4": "refuted", "F5": "refuted"}
    assert rep.scorecard["refuted"] == 4


def test_score_forecast_f4_confirmed_on_matching_rescue():
    pr = P.load_prereg()
    obs = {"steroid_rescue": True, "enantiomer_call": {"GR": "SS-OAc", "AR": "RR-OAc"}}
    rep = P.score_forecast(pr, obs)
    f4 = next(o for o in rep.outcomes if o.fid == "F4")
    assert f4.outcome == "confirmed"
