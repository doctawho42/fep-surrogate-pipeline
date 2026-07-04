"""Paper-2 orphan benchmark — public data acquisition + accessibility audit (Increment 1).

NO structure scoring here. This module only pulls ligand->target->holo-pocket triples and
caches them under data/paper2_bench/. audit_sources() is the cheapest-fail-first probe: if the
public sources are unreachable in this environment we stop before any curation effort.
"""
from __future__ import annotations

import urllib.request

# Candidate public sources. HEAD/GET probe only in the audit; real layout confirmed in Task 2.
SOURCE_URLS: dict[str, str] = {
    # LIT-PCBA (AVE-debiased actives + a bundled structure per target) — the orphan-honest anchor.
    "litpcba": "https://drugdesign.unistra.fr/LIT-PCBA/",
    # ChEMBL EBI REST (supplementary actives) — the repo already uses this API elsewhere.
    "chembl": "https://www.ebi.ac.uk/chembl/api/data/status.json",
    # RCSB (holo pockets) — used by dock.py already.
    "rcsb": "https://files.rcsb.org/",
    # ECOD (PDB-chain -> fold labels) — domain-classification flat file.
    "ecod": "http://prodata.swmed.edu/ecod/distributions/",
}


def audit_sources(urls: dict[str, str] = SOURCE_URLS, timeout: int = 20) -> dict[str, dict]:
    """Probe each source with a lightweight GET; return reachability + a short note. No parsing."""
    report: dict[str, dict] = {}
    for name, url in urls.items():
        try:
            req = urllib.request.Request(
                url,
                method="GET",
                headers={"User-Agent": "fep-paper2/1.0"},
            )
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                report[name] = {
                    "reachable": 200 <= resp.status < 400,
                    "note": f"HTTP {resp.status}",
                }
        except Exception as exc:  # noqa: BLE001 - audit must never crash; unreachable is a valid outcome
            report[name] = {"reachable": False, "note": f"{type(exc).__name__}: {exc}"[:200]}
    return report
