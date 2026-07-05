"""Unit tests for the cage prospective-loop decision scaffold (src/screen/prospective.py)."""
from __future__ import annotations

import inspect
import math

import pytest

from screen import prospective as P


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
