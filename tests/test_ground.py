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
