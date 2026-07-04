"""Paper-2 orphan benchmark — public data acquisition + accessibility audit.

NO structure scoring here. This module only pulls ligand->target->holo-pocket triples and
caches them under data/paper2_bench/. audit_sources() is the cheapest-fail-first probe: if the
public sources are unreachable in this environment we stop before any curation effort.

Task 2 scope: LIT-PCBA (AVE-debiased set) ONLY. The ChEMBL supplement and BindingDB breadth
fallback mentioned in the spec are deferred to a later increment gated on orphan-thinness.

Increment-2 prep: LIT-PCBA-specific parsing/download code has moved to `screen.sources.litpcba`.
This module now keeps the shared, source-agnostic helpers and `build_triples`, the thin
aggregator over all sources.
"""
from __future__ import annotations

import json
import pathlib
import urllib.request

import pandas as pd

from screen.sources import chembl_diverse, litpcba

# Candidate public sources. HEAD/GET probe only in the audit; real layout confirmed in Task 2.
SOURCE_URLS: dict[str, str] = {
    # LIT-PCBA (AVE-debiased actives + a bundled structure per target) — the orphan-honest anchor.
    # NOTE: the site's own https://lab.<host> redirect is broken (404); the plain https host
    # (no "lab." subdomain, no http->https redirect hop) is what actually serves the tarball.
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


# --- Triple assembly (Task 2) ------------------------------------------------------------

TRIPLE_COLS = [
    "mol_id", "smiles", "target", "pdb_id", "lig_resname", "affinity_nm", "source", "fold",
]


def triples_from_records(records: list[dict]) -> pd.DataFrame:
    """Normalize raw records into the canonical triple table; dedupe on (mol_id, target).

    `dropna` only catches None/NaN, so empty-string identifiers are dropped explicitly too —
    an empty `smiles` (or `mol_id`/`target`) is not a valid triple and must not survive to the
    cached parquet.
    """
    df = pd.DataFrame(records, columns=TRIPLE_COLS)
    df = df.dropna(subset=["mol_id", "smiles", "target"])
    for col in ("mol_id", "smiles", "target"):
        df = df[df[col].astype(str).str.strip() != ""]
    df = df.drop_duplicates(subset=["mol_id", "target"], keep="first").reset_index(drop=True)
    return df[TRIPLE_COLS]


def _pfam_for_pdb(pdb_id: str, timeout: int = 30, max_entities: int = 8) -> str | None:
    """Resolve a Pfam id for a PDB entry via RCSB's data API (ECOD returns 403 here).
    Tries polymer entities 1..max_entities in order and returns the first Pfam annotation
    found; falls back to the entry's first InterPro id if no entity carries a Pfam."""
    interpro_fallback: str | None = None
    for entity in range(1, max_entities + 1):
        url = f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id.upper()}/{entity}"
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "fep-paper2/1.0"})
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                data = json.loads(resp.read().decode("utf-8", "replace"))
        except Exception:  # noqa: BLE001 - missing entity / network hiccup -> stop trying
            break
        anns = data.get("rcsb_polymer_entity_annotation") or []
        for ann in anns:
            if ann.get("type") == "Pfam" and ann.get("annotation_id"):
                return ann["annotation_id"]
            if interpro_fallback is None and ann.get("type") == "InterPro" and ann.get(
                "annotation_id"
            ):
                interpro_fallback = ann["annotation_id"]
    return interpro_fallback


def _resolve_fold(target: str, pdb_id: str | None) -> str:
    """Fold label for a target: Pfam id of its representative PDB, else InterPro, else the
    conservative fallback of the target name itself (per Task-2 resolution #2)."""
    if pdb_id:
        pfam = _pfam_for_pdb(pdb_id)
        if pfam:
            return pfam
    return target


def assign_folds(records_by_target: dict[str, list[dict]]) -> dict[str, str]:
    """Resolve one fold label per target, given each target's parsed records (for `pdb_id`).

    Pure and network-free except for the `_resolve_fold` call it delegates to. Real PDB-keyed
    folds are cached by `pdb_id` (so targets sharing a representative structure share a Pfam
    lookup), but a structure-less target (`pdb_id is None`) NEVER shares that cache: it gets
    its own fold label, namely its target name. This is the Fix-1 invariant — without it, the
    first structure-less target's fold would get cached under a shared `None` key and silently
    reused for every later structure-less target.
    """
    fold_by_pdb: dict[str, str] = {}
    folds: dict[str, str] = {}
    for target, recs in records_by_target.items():
        pdb_id = recs[0]["pdb_id"] if recs else None
        if pdb_id:
            if pdb_id not in fold_by_pdb:
                fold_by_pdb[pdb_id] = _resolve_fold(target, pdb_id)
            folds[target] = fold_by_pdb[pdb_id]
        else:
            folds[target] = target
    return folds


def aggregate_records(records: list[dict]) -> pd.DataFrame:
    """Assign one Pfam fold per target across ALL sources, then shape+dedupe into the triple
    table. Fold resolution is done ONCE over the union of all sources' records for a given
    target, so a target seen from two sources doesn't get resolved twice (or inconsistently)."""
    by_target: dict[str, list[dict]] = {}
    for r in records:
        by_target.setdefault(r["target"], []).append(r)
    folds = assign_folds(by_target)
    for r in records:
        r["fold"] = folds[r["target"]]
    return triples_from_records(records)


def build_triples(
    cache_dir: str = "data/paper2_bench",
    limit: int | None = None,
    include_chembl_diverse: bool = True,
) -> pd.DataFrame:
    """Aggregate active->target->holo-pocket triples from all sources, assign folds, cache
    parquet. Idempotent (returns cached parquet if present). `limit` caps the number of targets
    kept, for a fast smoke run.

    `include_chembl_diverse=True` (default) merges LIT-PCBA + the ChEMBL-diverse targets into
    `triples_aggregate.parquet` (Increment 2). `include_chembl_diverse=False` reproduces the
    Increment-1 LIT-PCBA-only `triples.parquet` (kept for that increment's own tests/Fig M)."""
    cache = pathlib.Path(cache_dir)
    cache.mkdir(parents=True, exist_ok=True)
    name = "triples_aggregate.parquet" if include_chembl_diverse else "triples.parquet"
    out = cache / name
    if out.exists():
        return pd.read_parquet(out)

    records = litpcba.load_litpcba_records(cache_dir)
    if include_chembl_diverse:
        records += chembl_diverse.load_chembl_diverse_records(cache_dir)

    df = aggregate_records(records)
    if limit is not None:
        keep = set(sorted(df["target"].unique())[:limit])
        df = df[df["target"].isin(keep)].reset_index(drop=True)
    df.to_parquet(out)
    return df
