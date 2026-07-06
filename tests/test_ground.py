"""Unit tests for grounding logic (src/bar/ground.py).

Network calls are integration-tested separately (cached CSV).
"""
from __future__ import annotations

import math

from bar import ground as G


def test_exp_dg_conversion():
    # exp dG = -RT ln10 * pChEMBL; pIC50 7.0 -> -1.364*7 = -9.548 kcal/mol
    assert math.isclose(G.exp_dg(7.0), -1.364 * 7.0)


def test_coverage_and_grounded_gate():
    # a system maps 7/10 ligands -> coverage 0.7 >= 0.60 -> grounded; 5/10 -> not grounded
    assert G._coverage({f"l{i}": -9.0 for i in range(7)}, n_ligands=10) == 0.7
    assert G._coverage({f"l{i}": -9.0 for i in range(5)}, n_ligands=10) == 0.5


def test_median_pchembl_single_assay():
    # multiple IC50 records -> median pChEMBL; Ki ignored when assay=IC50
    recs = [{"standard_type": "IC50", "pchembl_value": "7.0"},
            {"standard_type": "IC50", "pchembl_value": "7.5"},
            {"standard_type": "Ki", "pchembl_value": "9.9"}]
    assert G._median_pchembl(recs, assay="IC50") == 7.25
    assert G._median_pchembl(recs, assay="Ki") == 9.9
    null_rec = [{"standard_type": "IC50", "pchembl_value": None}]
    assert G._median_pchembl(null_rec, assay="IC50") is None


def test_assay_fallback_prefers_first_available():
    # IC50 has matches -> IC50 chosen even though Ki also has matches (never mix, first wins)
    both = {"a": {"IC50": 7.0, "Ki": 6.9}, "b": {"IC50": 6.5}}
    assert G.select_assay(both, ["IC50", "Ki"]) == "IC50"

    # IC50 has zero matches for this system, Ki has matches -> fall back to Ki
    ic50_empty = {"a": {"Ki": 7.0}, "b": {"Ki": 6.5}}
    assert G.select_assay(ic50_empty, ["IC50", "Ki"]) == "Ki"

    # neither assay has any match -> falls back to the last assay_order entry (definite, no raise)
    assert G.select_assay({"a": {}}, ["IC50", "Ki"]) == "Ki"
    assert G.select_assay({}, ["IC50", "Ki"]) == "Ki"


def test_collapsed_enantiomer_pairs_hif2a():
    # hif2a ligands 237/15 are ChEMBL-flexmatch-collapsed enantiomers (same pchembl=8.0, 2D graph
    # identical, stored SMILES differ only in one stereocenter) -- reads the committed CSV.
    pairs = G.collapsed_enantiomer_pairs("hif2a")
    normalized = {frozenset(p) for p in pairs}
    assert frozenset({"237", "15"}) in normalized
