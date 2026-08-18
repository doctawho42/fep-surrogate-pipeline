"""Unit tests for the close-the-loop accuracy race (src/bar/closeloop.py)."""
from __future__ import annotations

import math
import pathlib

import numpy as np
import pytest

from bar import closeloop as C

# external witness; independent literal, never edited
PINNED_SHA256 = "35cb3b5a81c6878fddaa2f8cfc0107596fbdd9f1c6de3ab132524632c97ff86e"


def test_prereg_loads_and_matches_pinned_anchor():
    assert C.CLOSELOOP_PREREG_SHA256 == PINNED_SHA256
    assert C.sha256_of(C.CLOSELOOP_PREREG_PATH) == PINNED_SHA256
    pr = C.load_prereg()
    assert pr.systems == ["cdk8", "hif2a", "p38", "bace"]
    assert pr.coverage_min == 0.60
    assert pr.n_perm == 1000
    assert pr.min_grounded == 2
    assert pr.assay_order == ["IC50", "Ki"]


def test_prereg_immutable_against_external_anchor(tmp_path):
    tampered = tmp_path / "t.yaml"
    tampered.write_bytes(pathlib.Path(C.CLOSELOOP_PREREG_PATH).read_bytes() + b"\n# edit\n")
    with pytest.raises(ValueError):
        C.load_prereg(str(tampered), expected_sha256=PINNED_SHA256)


def _line_network(n=8, bad_edge_err=3.0, seed=1):
    """A path graph 0-1-...-n with true node dG = i; all edges se=1, one edge given a large
    systematic error. exp = true dG. Removing the bad edge should improve MUE.

    NOTE (deviation from task-3-brief.md, flagged in task-3-report.md): the brief's fixture put
    the bad edge at i == 3, which is edge (3, 4) -- a bridge, not part of either triangle cycle
    (0,1,2) or (5,6,7). Cycle-closure chi^2 is a sum over *cycle* residuals only, so a bridge's
    systematic error is invisible to gls_network/repair_order by construction (this is the same
    "closure is blind to node-consistent bias" property documented in bar/qc.py) -- repair_order
    then removes nothing and the test's `> 0` assertion is unreachable. Moving the bad edge to
    i == 0 (edge (0, 1), which IS on the first triangle cycle) preserves every other stated
    property of the fixture (path + 2 chords, n=8, deterministic seed) while making it a real
    stress test of repair_order.
    """
    rng = np.random.default_rng(seed)
    exp = {i: float(i) for i in range(n)}
    edges = []
    for i in range(n - 1):
        err = bad_edge_err if i == 0 else 0.0
        edges.append((i, i + 1, 1.0 + err + 0.01 * rng.standard_normal(), 1.0))
    # add a couple of chords so removal leaves cycles/connectivity
    edges.append((0, 2, 2.0, 1.0))
    edges.append((5, 7, 2.0, 1.0))
    return edges, exp


def test_delta_mue_positive_when_bad_edge_removed():
    edges, exp = _line_network()
    from bar.qc import repair_order
    removed, _ = repair_order(edges, target_reduced_chi2=1.0)
    # guided removal improves accuracy on this synthetic
    assert C.delta_mue(edges, removed, exp) > 0


def test_system_effect_guided_beats_random_on_synthetic():
    edges, exp = _line_network()
    eff = C.system_effect(edges, exp, n_perm=500, target_rchi2=1.0, seed=0)
    assert eff["guided"] > eff["random_mean"]
    # guided significantly better than random on the seeded network
    assert eff["p"] < 0.2
    assert set(eff) >= {"k", "guided", "random_mean", "p", "below_5pct", "curve"}


def test_combine_stouffer_and_sign():
    # NOTE (deviation from task-3-brief.md, flagged in task-3-report.md): the brief's `effs`
    # fixtures omitted "below_5pct". combine()'s SUCCESS branch requires a majority of systems
    # with below_5pct=True; without the key, e.get("below_5pct") is always None -> False -> the
    # verdict is unreachably "NULL" even though stouffer_p < 0.05. Added below_5pct=True,
    # consistent with p=0.01/0.03 (both would clear a random-null 5th-percentile cut).
    effs = [{"p": 0.01, "guided": 0.5, "random_mean": 0.1, "below_5pct": True},
            {"p": 0.03, "guided": 0.4, "random_mean": 0.1, "below_5pct": True}]
    out = C.combine(effs)
    assert out["stouffer_p"] < 0.05 and out["verdict"] == "SUCCESS"
    null = C.combine([{"p": 0.6, "guided": 0.0, "random_mean": 0.1, "below_5pct": False}] * 2)
    assert null["verdict"] == "NULL"


def test_combine_drops_nan_p_consistently():
    # a NaN-p system must be dropped from ALL sub-statistics (n_effective=2), not just Stouffer
    effs = [{"p": 0.01, "guided": 0.5, "random_mean": 0.1, "below_5pct": True},
            {"p": 0.02, "guided": 0.4, "random_mean": 0.1, "below_5pct": True},
            {"p": float("nan"), "guided": 0.0, "random_mean": 0.1, "below_5pct": False}]
    out = C.combine(effs)
    assert out["stouffer_p"] < 0.05       # NaN dropped -> Stouffer over the 2 valid systems
    assert out["sign_p"] == 0.25          # 2/2 successes over n=2 -> comb(2,2)/2**2 = 0.25
    assert out["verdict"] == "SUCCESS"
    assert C.combine([{"p": float("nan"), "guided": 0.0, "random_mean": 0.1}])["verdict"] == "NULL"


def test_combine_refuses_to_report_a_saturated_stouffer():
    """A permutation p of exactly 1 drives Stouffer to exactly 1 whatever the other systems say.

    That is the resolution of the permutation grid speaking, not a combined p-value, so the
    function reports NaN and names the saturated input instead of emitting the artefact.
    """
    effs = [
        {"p": 0.937, "guided": -0.05, "random_mean": 0.0, "below_5pct": False},
        {"p": 0.819, "guided": -0.10, "random_mean": 0.0, "below_5pct": False},
        {"p": 1.000, "guided": -0.61, "random_mean": 0.0, "below_5pct": False},
        {"p": 0.023, "guided": +0.14, "random_mean": 0.0, "below_5pct": True},
    ]
    out = C.combine(effs, n_perm=1000)
    assert math.isnan(out["stouffer_p"])
    assert out["stouffer_saturated"] == [1.0]
    # the clipped sensitivity is finite, and says the same thing the raw artefact hid: not extreme
    assert 0.9 < out["stouffer_p_clipped"] < 1.0
    # the pre-registered decision is unchanged by the NaN, and the sign test still works
    assert out["verdict"] == "NULL"
    assert out["sign_p"] == pytest.approx(15 / 16)


def test_combine_is_unchanged_when_no_p_saturates():
    """The published close-the-loop numbers must not move: none of their p-values saturate."""
    effs = [{"p": p, "guided": 0.0, "random_mean": 0.1, "below_5pct": False}
            for p in (0.905, 0.136, 0.880, 0.071)]
    out = C.combine(effs, n_perm=1000)
    assert out["stouffer_saturated"] == []
    assert out["stouffer_p"] == pytest.approx(0.4838, abs=1e-3)
    assert out["stouffer_p_clipped"] == pytest.approx(out["stouffer_p"], abs=1e-9)
