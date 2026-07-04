"""Paper-2 orphan benchmark — public data acquisition + accessibility audit.

NO structure scoring here. This module only pulls ligand->target->holo-pocket triples and
caches them under data/paper2_bench/. audit_sources() is the cheapest-fail-first probe: if the
public sources are unreachable in this environment we stop before any curation effort.

Task 2 scope: LIT-PCBA (AVE-debiased set) ONLY. The ChEMBL supplement and BindingDB breadth
fallback mentioned in the spec are deferred to a later increment gated on orphan-thinness.
"""
from __future__ import annotations

import json
import pathlib
import re
import tarfile
import urllib.request

import pandas as pd

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

# Real download URL for the AVE-debiased (orphan-honest) LIT-PCBA tarball, confirmed by probe:
# HTTP 200, Content-Length 57396648 (~57 MB, well under the 2 GB stop threshold).
LITPCBA_TGZ_URL = "https://drugdesign.unistra.fr/LIT-PCBA/Files/AVE_unbiased.tgz"

# Confirmed against the extracted tarball layout (Task 2 probe): 15 targets, each a folder of
# <pdbid>_ligand.mol2 / <pdbid>_protein.mol2 structure pairs plus active_T.smi + active_V.smi
# (train/validation split of actives; no combined "actives.smi", no bundled PDB file).
LITPCBA_TARGETS: list[str] = [
    "ADRB2", "ALDH1", "ESR1_ago", "ESR1_ant", "FEN1", "GBA", "IDH1", "KAT2A",
    "MAPK1", "MTORC1", "OPRK1", "PKM2", "PPARG", "TP53", "VDR",
]
ACTIVES_FILES = ["active_T.smi", "active_V.smi"]  # supersedes the brief's single "actives.smi"

# Fallback resname when a ligand .mol2 carries no parseable PDB het code.
DEFAULT_LIG_RESNAME = "LIG"

_CHAIN_RESID_RE = re.compile(r"^[A-Za-z0-9]{1,4}\d*$")
_CHAIN_ONLY_RE = re.compile(r"^[A-Za-z]$")
_CHAIN_GLUED_RESID_RE = re.compile(r"^[A-Za-z]\d+$")


def _resname_from_substructure_line(line: str) -> str | None:
    """Parse the PDB het code out of a Tripos SUBSTRUCTURE record's data row, e.g.:
    '     1 Q3XG603    55 GROUP             4 A     ****    0 ROOT 3XG A 603' -> '3XG'.
    Not every ligand .mol2 carries this block (30/333 in the real AVE_unbiased tarball do);
    return None when absent or unparseable so the caller can fall back to DEFAULT_LIG_RESNAME.
    """
    parts = line.split()
    while parts and parts[-1] in ("ROOT", "****", "0"):
        parts.pop()
    if len(parts) >= 3 and _CHAIN_RESID_RE.match(parts[-1]) and _CHAIN_ONLY_RE.match(parts[-2]):
        return parts[-3]
    if len(parts) >= 2 and _CHAIN_GLUED_RESID_RE.match(parts[-1]):
        return parts[-2]
    return None


def _resname_from_ligand_mol2(path: pathlib.Path) -> str:
    """Best-effort PDB het code for a LIT-PCBA ligand .mol2; DEFAULT_LIG_RESNAME if unavailable."""
    try:
        text = path.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return DEFAULT_LIG_RESNAME
    lines = text.splitlines()
    for i, line in enumerate(lines):
        if line.strip() == "@<TRIPOS>SUBSTRUCTURE":
            if i + 1 < len(lines):
                resname = _resname_from_substructure_line(lines[i + 1])
                if resname:
                    return resname
            break
    return DEFAULT_LIG_RESNAME


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


def _pick_representative_structure(target_dir: pathlib.Path) -> tuple[str | None, str]:
    """Pick one representative PDB id + ligand resname for a LIT-PCBA target folder.
    Picks the ligand .mol2 with the smallest file (fewest atoms) as a cheap, deterministic
    tie-break among the several bundled holo structures per target."""
    ligand_files = sorted(target_dir.glob("*_ligand.mol2"))
    if not ligand_files:
        return None, DEFAULT_LIG_RESNAME
    chosen = min(ligand_files, key=lambda p: p.stat().st_size)
    pdb_id = chosen.name.removesuffix("_ligand.mol2")
    resname = _resname_from_ligand_mol2(chosen)
    return pdb_id, resname


def parse_litpcba_target(target_dir: str) -> list[dict]:
    """Return active-ligand records for one extracted LIT-PCBA target folder.

    Supersedes the brief's HTTP-index approach: the real AVE-debiased tarball ships each
    target as a folder of <pdbid>_ligand.mol2 / <pdbid>_protein.mol2 structure pairs plus
    active_T.smi + active_V.smi (train/val split; no single "actives.smi", no bundled PDB
    file). `target` is the folder's basename (the LIT-PCBA target key, e.g. "ADRB2").
    Affinity is left None: LIT-PCBA is an actives set, not a Kd table (see brief Step 5).
    `fold` is left None here; build_triples fills it once per target via the RCSB Pfam lookup.
    """
    tdir = pathlib.Path(target_dir)
    target = tdir.name
    pdb_id, lig_resname = _pick_representative_structure(tdir)

    recs: list[dict] = []
    seen_mids: set[str] = set()
    for fname in ACTIVES_FILES:
        fpath = tdir / fname
        if not fpath.exists():
            continue
        text = fpath.read_text(encoding="utf-8", errors="replace")
        for line in text.splitlines():
            parts = line.split()
            if len(parts) >= 2:
                smi, mid_raw = parts[0], parts[1]
                mol_id = f"{target}:{mid_raw}"
                if mol_id in seen_mids:
                    continue
                seen_mids.add(mol_id)
                recs.append({
                    "mol_id": mol_id, "smiles": smi, "target": target,
                    "pdb_id": pdb_id, "lig_resname": lig_resname,
                    "affinity_nm": None, "source": "litpcba", "fold": None,
                })
    return recs


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


def build_triples(cache_dir: str = "data/paper2_bench", limit: int | None = None) -> pd.DataFrame:
    """Assemble the triple table from LIT-PCBA (Task-2 scope; see module docstring) and cache
    it. Idempotent: returns the cached parquet if present. `limit` caps the number of targets
    walked, for a fast smoke run. Downloads + extracts the AVE-debiased tarball on first use."""
    cache = pathlib.Path(cache_dir)
    cache.mkdir(parents=True, exist_ok=True)
    out = cache / "triples.parquet"
    if out.exists():
        return pd.read_parquet(out)

    litpcba_dir = cache / "litpcba"
    tgz_path = cache / "AVE_unbiased.tgz"
    if not litpcba_dir.exists() or not any(litpcba_dir.iterdir()):
        if not tgz_path.exists():
            urllib.request.urlretrieve(LITPCBA_TGZ_URL, tgz_path)
        litpcba_dir.mkdir(parents=True, exist_ok=True)
        with tarfile.open(tgz_path, "r:gz") as tf:
            tf.extractall(litpcba_dir)  # noqa: S202 - trusted first-party academic mirror

    targets = LITPCBA_TARGETS[:limit] if limit is not None else LITPCBA_TARGETS
    records_by_target: dict[str, list[dict]] = {}
    for target in targets:
        target_dir = litpcba_dir / target
        if not target_dir.is_dir():
            continue
        records_by_target[target] = parse_litpcba_target(str(target_dir))

    folds = assign_folds(records_by_target)
    records: list[dict] = []
    for target, recs in records_by_target.items():
        fold = folds[target]
        for rec in recs:
            rec["fold"] = fold
        records.extend(recs)

    df = triples_from_records(records)
    df.to_parquet(out)
    return df
