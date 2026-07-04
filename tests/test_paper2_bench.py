"""Test suite for Paper-2 benchmark infrastructure (Increment 1)."""
import pathlib

import pandas as pd

from screen.bench_sources import (
    DEFAULT_LIG_RESNAME,
    _resname_from_ligand_mol2,
    _resname_from_substructure_line,
    _resolve_fold,
    assign_folds,
    audit_sources,
    parse_litpcba_target,
    triples_from_records,
)
from screen.stratify import assign_stratum, ecfp, max_tanimoto, stratify

FIXTURE_DIR = pathlib.Path(__file__).parent / "fixtures" / "litpcba_mini"


def test_audit_reports_unreachable_without_crashing():
    """Verify audit_sources handles unreachable URLs gracefully."""
    fake = {"litpcba": "http://127.0.0.1:9/nope", "chembl": "http://127.0.0.1:9/nope"}
    report = audit_sources(fake, timeout=1)
    assert set(report) == {"litpcba", "chembl"}
    assert report["litpcba"]["reachable"] is False
    assert isinstance(report["litpcba"]["note"], str)


def test_triples_from_records_shapes_and_dedupes():
    recs = [
        {"mol_id": "a", "smiles": "CCO", "target": "T1", "pdb_id": "1ABC", "lig_resname": "LIG",
         "affinity_nm": 10.0, "source": "litpcba", "fold": "PF00001"},
        {"mol_id": "a", "smiles": "CCO", "target": "T1", "pdb_id": "1ABC", "lig_resname": "LIG",
         "affinity_nm": 12.0, "source": "chembl", "fold": "PF00001"},   # dup (mol_id,target)
        {"mol_id": "b", "smiles": "c1ccccc1", "target": "T2", "pdb_id": "2XYZ",
         "lig_resname": "BEN", "affinity_nm": 50.0, "source": "litpcba", "fold": "PF00002"},
    ]
    df = triples_from_records(recs)
    assert list(df.columns) == ["mol_id", "smiles", "target", "pdb_id", "lig_resname",
                                "affinity_nm", "source", "fold"]
    assert len(df) == 2                                  # (a,T1) deduped
    assert set(df["target"]) == {"T1", "T2"}


def test_triples_from_records_drops_empty_smiles():
    """Fix 3: dropna() alone misses empty strings; triples_from_records must drop them too."""
    recs = [
        {"mol_id": "a", "smiles": "CCO", "target": "T1", "pdb_id": "1ABC", "lig_resname": "LIG",
         "affinity_nm": 10.0, "source": "litpcba", "fold": "PF00001"},
        {"mol_id": "b", "smiles": "", "target": "T1", "pdb_id": "1ABC", "lig_resname": "LIG",
         "affinity_nm": 10.0, "source": "litpcba", "fold": "PF00001"},
    ]
    df = triples_from_records(recs)
    assert len(df) == 1
    assert set(df["mol_id"]) == {"a"}


# --- Fixture-based parsing tests (offline, deterministic) ----------------------------------


def test_parse_litpcba_target_reads_fixture_smiles():
    recs = parse_litpcba_target(str(FIXTURE_DIR / "TARGET_A"))
    assert len(recs) == 3
    assert {r["smiles"] for r in recs} == {"CCO", "c1ccccc1", "CC(=O)Oc1ccccc1C(=O)O"}
    for r in recs:
        assert r["target"] == "TARGET_A"
        assert r["source"] == "litpcba"
        assert r["mol_id"].startswith("TARGET_A:")


def test_parse_litpcba_target_picks_representative_structure_and_resname():
    recs = parse_litpcba_target(str(FIXTURE_DIR / "TARGET_A"))
    assert recs[0]["pdb_id"] == "1abc"
    assert recs[0]["lig_resname"] == "3XG"


def test_resname_from_substructure_line_extracts_real_code():
    line = "     1 Q3XG603    55 GROUP             4 A     ****    0 ROOT 3XG A 603"
    assert _resname_from_substructure_line(line) == "3XG"


def test_resname_from_ligand_mol2_with_substructure_block():
    mol2 = FIXTURE_DIR / "TARGET_A" / "1abc_ligand.mol2"
    assert _resname_from_ligand_mol2(mol2) == "3XG"


def test_resname_from_ligand_mol2_falls_back_without_substructure_block():
    mol2 = FIXTURE_DIR / "TARGET_B" / "9xyz_ligand.mol2"
    assert _resname_from_ligand_mol2(mol2) == DEFAULT_LIG_RESNAME


# --- Fold resolution / fallback chain -------------------------------------------------------


def test_resolve_fold_falls_through_to_target_name_when_pfam_missing(monkeypatch):
    """Monkeypatch _pfam_for_pdb to return None -> _resolve_fold falls through to the
    conservative target-name fallback (InterPro fallback is internal to _pfam_for_pdb)."""
    monkeypatch.setattr("screen.bench_sources._pfam_for_pdb", lambda pdb_id, **kw: None)
    assert _resolve_fold("SOME_TARGET", "1abc") == "SOME_TARGET"


def test_resolve_fold_uses_pfam_when_available(monkeypatch):
    monkeypatch.setattr("screen.bench_sources._pfam_for_pdb", lambda pdb_id, **kw: "PF00001")
    assert _resolve_fold("SOME_TARGET", "1abc") == "PF00001"


def test_resolve_fold_with_no_pdb_id_uses_target_name():
    assert _resolve_fold("SOME_TARGET", None) == "SOME_TARGET"


# --- Fix 1 regression: fold-cache None-key bug ----------------------------------------------


def test_assign_folds_gives_distinct_labels_to_structure_less_targets(monkeypatch):
    """Regression for Fix 1: two structure-less targets (pdb_id=None) must each get their OWN
    fold label (their target name), never a value cached under a shared `None` key."""
    monkeypatch.setattr("screen.bench_sources._pfam_for_pdb", lambda pdb_id, **kw: "PF00001")
    records_by_target = {
        "TARGET_NOSTRUCT_1": [{"pdb_id": None, "target": "TARGET_NOSTRUCT_1"}],
        "TARGET_NOSTRUCT_2": [{"pdb_id": None, "target": "TARGET_NOSTRUCT_2"}],
        "TARGET_WITH_STRUCT": [{"pdb_id": "1abc", "target": "TARGET_WITH_STRUCT"}],
    }
    folds = assign_folds(records_by_target)
    assert folds["TARGET_NOSTRUCT_1"] == "TARGET_NOSTRUCT_1"
    assert folds["TARGET_NOSTRUCT_2"] == "TARGET_NOSTRUCT_2"
    assert folds["TARGET_NOSTRUCT_1"] != folds["TARGET_NOSTRUCT_2"]
    assert folds["TARGET_WITH_STRUCT"] == "PF00001"


def test_assign_folds_shares_cache_across_targets_with_same_pdb(monkeypatch):
    """Targets that DO share a representative pdb_id may share the resolved fold (this is the
    intended caching behavior Fix 1 preserves for the real, PDB-keyed case)."""
    calls = []

    def fake_pfam(pdb_id, **kw):
        calls.append(pdb_id)
        return "PF00099"

    monkeypatch.setattr("screen.bench_sources._pfam_for_pdb", fake_pfam)
    records_by_target = {
        "TARGET_X": [{"pdb_id": "1abc", "target": "TARGET_X"}],
        "TARGET_Y": [{"pdb_id": "1abc", "target": "TARGET_Y"}],
    }
    folds = assign_folds(records_by_target)
    assert folds["TARGET_X"] == folds["TARGET_Y"] == "PF00099"
    assert calls == ["1abc"]  # resolved once, cached by pdb_id for the second target


# --- Stratification tests (orphan benchmark) --


def test_tanimoto_self_is_one_and_stratum_cuts():
    fp = ecfp("CCO")
    assert abs(max_tanimoto(fp, [ecfp("CCO")]) - 1.0) < 1e-9
    assert max_tanimoto(fp, []) == 0.0
    assert assign_stratum(0.7) == "high"
    assert assign_stratum(0.4) == "mid"
    assert assign_stratum(0.30) == "orphan"
    assert assign_stratum(0.15) == "deep_orphan"


def test_stratify_self_excludes_query_from_active_pool():
    # a query identical to an active of ANOTHER pocket must NOT count itself -> orphan if unique
    queries = pd.DataFrame({"mol_id": ["q1"], "smiles": ["CCO"], "target": ["T1"]})
    pocket_actives = {"T1": ["CCO"], "T2": ["c1ccccc1C(=O)O"]}   # T1's only active IS the query
    out = stratify(queries, pocket_actives)
    assert out.loc[0, "s"] < 0.35            # self excluded -> dissimilar to benzoic acid -> orphan
    assert out.loc[0, "stratum"] in ("orphan", "deep_orphan")
