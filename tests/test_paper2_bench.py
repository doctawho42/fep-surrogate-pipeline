"""Test suite for Paper-2 benchmark infrastructure (Increment 1)."""
from screen.bench_sources import audit_sources, triples_from_records


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
