"""The NIOCH cage report must be reproducible from its caches, byte for byte.

This report is a hand-off to an experimental group, and it was edited by hand once: a Boltz-2
cross-check section was typed straight into the markdown, and the recommendation order was
rewritten to demote the nuclear-receptor panel after that cross-check contradicted the docking
hypothesis. Neither change reached the generator, so the next `make nioch` would have deleted the
section and restored a recommendation the evidence had retracted, telling a wet lab to run the
wrong assay first. The test below makes that impossible: the committed report must equal a fresh
generation, so the generator and the document cannot disagree without the suite failing.
"""
from __future__ import annotations

import pathlib
import sys

import pytest

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "figs"))

REPORT = ROOT / "docs" / "nioch_cage_report.md"


@pytest.fixture(scope="module")
def generated(tmp_path_factory):
    import make_nioch as M
    out = tmp_path_factory.mktemp("nioch") / "report.md"
    M.report(M.rank(), out)
    return out.read_text()


def test_committed_report_is_exactly_what_the_generator_produces(generated):
    assert generated == REPORT.read_text(), (
        "docs/nioch_cage_report.md has drifted from figs/make_nioch.py. Re-run `make nioch` and "
        "commit the result, or move the hand-written change into the generator."
    )


def test_the_recommendation_order_follows_the_boltz_verdict(generated):
    """The order is the part a hand edit got right and a regeneration would have undone."""
    import make_nioch as M
    _lines, challenged, _ranked = M.boltz_cross_check()
    body = generated.split("## Recommended experiments")[1]
    first = body.split("\n2.")[0]
    if challenged:
        assert "Broad biochemical / phenotypic profiling" in first
        assert "CHALLENGED" in generated
    else:
        assert "Nuclear-receptor reporter panel" in first


def test_the_broad_profiling_item_appears_once(generated):
    """The hand-edited version promoted broad profiling to the lead and left the old item behind."""
    body = generated.split("## Recommended experiments")[1]
    assert body.count("is the real un-blinder") + body.count("the primary un-blinder") == 1


def test_the_boltz_numbers_come_from_the_committed_screen_results(generated):
    """Spot-check that the section is derived, not transcribed: the anchors must match the data."""
    import make_boltz_cage as mbc
    for target in ("GR", "AR", "ER"):
        anchor = mbc.load_target_conf(target)["anchor"]
        assert f"{anchor:.2f}" in generated
