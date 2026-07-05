"""BindingDB breadth supplement -> active-ligand records for the Paper-2 collapse-count
re-test (Increment 2). NO structure scoring here.

Access-path note (Increment-2 empirical test, 2026-07-05): BindingDB's legacy
`axis2/services/BDBService/getLigandsByUniprot` endpoint is gone (404), and the documented
`/rest/getLigandsByUniprot` JSON API (advertised on `BindingDBRESTfulAPI.jsp`) times out from
this environment on every host/scheme combination tried (www/bare host, http/https). The full
`BindingDB_All` bulk dump is a multi-GB download, ruled out by the task brief. The one
tractable, reachable path is the **HTML** by-target search page
`rwd/bind/ByUniProtids.jsp?uniprotids=<UniProt>&cutoff=<nM>`, which returns HTTP 200 with a
real per-ligand results table (SMILES via `setClipboard('<smiles>')`, target name, `monomerid`,
IC50/Ki/Kd/EC50). This module scrapes that page with a plain regex (no HTML parser dependency
needed — the SMILES/InChI live inside single-quoted `onclick` JS calls, not markup). Verified
against ByUniProtids.jsp for EGFR/CA9/BTK/renin/ADORA2A: all 200, all with real ligand rows.
Pagination (`startPg`) is NOT honored by ByUniProtids.jsp directly (returns the same first
~50-155 rows regardless); the only way to page is a second request to
`PrimarySearch_ki.jsp` with the resolved `polymerid`/`complexid` list scraped from page 1 -- NOT
implemented here because the collapse-count test only needs chemotype-diverse breadth
(dozens of ligands per target), not exhaustive per-target coverage, matching the depth already
used by `chembl_diverse.fetch_actives`'s `max_records` cap.
"""
from __future__ import annotations

import re
import urllib.request

# 15 chemotype/fold-diverse targets NOT in the existing 29-target aggregate (see
# docs/ or `triples_aggregate.parquet`'s `target` column for the existing set: kinases
# EGFR/CDK2/P38, proteases HIVPR/FXA/THROMB, GPCR ADRB2, nuclear receptors ESR1/GR/PPARG/VDR,
# hydrolases ACHE/CAII/DHFR/HSP90/MMP9/BACE1, etc.). Each entry pairs a UniProt accession with
# a holo PDB (representative bound structure) so `_pfam_for_pdb` can resolve a real fold label,
# spanning kinases (RTK + non-RTK), proteases (aspartic + cysteine), an ADP-ribosyltransferase,
# a deacetylase, a glycosyltransferase, a phosphatase, an oxidoreductase, and four aminergic
# GPCRs -- deliberately overlapping fold classes already in the pool (kinase, GPCR) as well as
# NEW ones (PARP, HDAC, NAMPT, PTP, NOS), consistent with a "breadth" addition.
DIVERSE_TARGETS: list[dict] = [
    {"target": "JAK2",    "uniprot": "O60674", "pdb_id": "2B7A", "lig_resname": "STI"},
    {"target": "BTK",     "uniprot": "Q06187", "pdb_id": "3GEN", "lig_resname": "B96"},
    {"target": "ALK",     "uniprot": "Q9UM73", "pdb_id": "2XP2", "lig_resname": "VX6"},
    {"target": "MET",     "uniprot": "P08581", "pdb_id": "3LQ8", "lig_resname": "03P"},
    {"target": "REN",     "uniprot": "P00797", "pdb_id": "2REN", "lig_resname": "CGP"},
    {"target": "CATL",    "uniprot": "P07711", "pdb_id": "2XU3", "lig_resname": "07C"},
    {"target": "PARP1",   "uniprot": "P09874", "pdb_id": "1UK0", "lig_resname": "L3P"},
    {"target": "HDAC2",   "uniprot": "Q92769", "pdb_id": "3MAX", "lig_resname": "SHH"},
    {"target": "NAMPT",   "uniprot": "P43490", "pdb_id": "2GVJ", "lig_resname": "APO"},
    {"target": "PTP1B",   "uniprot": "P18031", "pdb_id": "1PTY", "lig_resname": "788"},
    {"target": "ADORA2A", "uniprot": "P29274", "pdb_id": "3EML", "lig_resname": "ZMA"},
    {"target": "DRD3",    "uniprot": "P35462", "pdb_id": "3PBL", "lig_resname": "ETQ"},
    {"target": "HTR2A",   "uniprot": "P28223", "pdb_id": "6A93", "lig_resname": "8XC"},
    {"target": "CB1",     "uniprot": "P21554", "pdb_id": "5TGZ", "lig_resname": "8HW"},
    {"target": "NOS3",    "uniprot": "P29474", "pdb_id": "1P6L", "lig_resname": "P6L"},
]

TRIPLE_KEYS = [
    "mol_id", "smiles", "target", "pdb_id", "lig_resname", "affinity_nm", "source", "fold",
]

# The results page embeds each ligand's SMILES in a `Copy SMILES` button's onclick JS call,
# e.g. onclick="setClipboard('COc1cc(...)cc1Nc1nccc(n1)-c1cn(C)c2ccccc12')" -- single-quoted,
# no HTML entities inside (unlike the surrounding markup). Each ligand row also carries a
# `monomerid=<int>` BindingDB compound id we use to dedupe (mirrors ChEMBL's molecule_chembl_id).
_SMILES_RE = re.compile(r"setClipboard\('([^'\"<>]+)'\)\s*\">\s*Copy&nbsp;SMILES")
_MONOMERID_RE = re.compile(r"monomerid=(\d+)")


def _extract_ligand_rows(html: str) -> list[tuple[str, str]]:
    """Pull (monomerid, smiles) pairs out of a ByUniProtids.jsp results page, in document
    order. Each ligand row's `monomerid=<id>` anchor appears immediately before that ligand's
    `Copy SMILES` button, so pairing the two regexes' matches positionally (by scanning once
    left-to-right) recovers the association without a full HTML parse."""
    pairs: list[tuple[str, str]] = []
    last_monomerid: str | None = None
    # Single pass over interleaved matches, ordered by position in the document.
    events = [(m.start(), "id", m.group(1)) for m in _MONOMERID_RE.finditer(html)]
    events += [(m.start(), "smiles", m.group(1)) for m in _SMILES_RE.finditer(html)]
    events.sort(key=lambda e: e[0])
    for _, kind, val in events:
        if kind == "id":
            last_monomerid = val
        elif kind == "smiles" and last_monomerid is not None:
            pairs.append((last_monomerid, val))
            last_monomerid = None
    return pairs


def records_from_html(html: str, target: str, pdb_id: str, lig_resname: str) -> list[dict]:
    """Shape a ByUniProtids.jsp results page into benchmark records; drop empty/dup molecules."""
    seen, recs = set(), []
    for monomerid, smi in _extract_ligand_rows(html):
        if not smi:
            continue
        mid = f"{target}:BDBM{monomerid}"
        if mid in seen:
            continue
        seen.add(mid)
        recs.append({"mol_id": mid, "smiles": smi, "target": target, "pdb_id": pdb_id,
                     "lig_resname": lig_resname, "affinity_nm": None, "source": "bindingdb",
                     "fold": None})
    return recs


def fetch_html(uniprot: str, cutoff_nm: int = 10000, timeout: int = 60) -> str:
    """Fetch the ByUniProtids.jsp HTML results page for one UniProt accession."""
    url = (f"https://www.bindingdb.org/rwd/bind/ByUniProtids.jsp"
           f"?uniprotids={uniprot}&cutoff={cutoff_nm}")
    req = urllib.request.Request(url, headers={"User-Agent": "paper2/1.0"})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return resp.read().decode("utf-8", "replace")


def load_bindingdb_records(cache_dir: str = "data/paper2_bench") -> list[dict]:
    """All BindingDB breadth-target active records (fold left None; the aggregator assigns
    Pfam folds). `cache_dir` is accepted for signature symmetry with the other loaders but
    unused (no per-source cache file; the aggregate parquet is the cache)."""
    out: list[dict] = []
    for t in DIVERSE_TARGETS:
        html = fetch_html(t["uniprot"])
        out.extend(records_from_html(html, t["target"], t["pdb_id"], t["lig_resname"]))
    return out
