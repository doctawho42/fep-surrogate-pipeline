"""Test suite for Paper-2 benchmark infrastructure (Increment 1)."""
from screen.bench_sources import audit_sources


def test_audit_reports_unreachable_without_crashing():
    """Verify audit_sources handles unreachable URLs gracefully."""
    fake = {"litpcba": "http://127.0.0.1:9/nope", "chembl": "http://127.0.0.1:9/nope"}
    report = audit_sources(fake, timeout=1)
    assert set(report) == {"litpcba", "chembl"}
    assert report["litpcba"]["reachable"] is False
    assert isinstance(report["litpcba"]["note"], str)
