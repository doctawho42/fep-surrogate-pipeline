"""LIT-PCBA (AVE-debiased) source — download/extract the tarball and parse per-target actives.

Moved out of `screen.bench_sources` (Paper-2 increment-2 prep): this module owns everything
specific to the LIT-PCBA layout. `screen.bench_sources` keeps the source-agnostic helpers
(`triples_from_records`, `assign_folds`, `_pfam_for_pdb`, `_resolve_fold`, `audit_sources`,
`TRIPLE_COLS`) and aggregates records from this (and, later, other) sources.
"""
from __future__ import annotations

import pathlib
import re
import tarfile
import urllib.request

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


def load_litpcba_records(cache_dir: str = "data/paper2_bench") -> list[dict]:
    """Download+extract AVE_unbiased.tgz (if needed) and return active-ligand records for all
    LIT-PCBA targets (TRIPLE_COLS keys; `fold` left None — the aggregator assigns folds)."""
    cache = pathlib.Path(cache_dir)
    cache.mkdir(parents=True, exist_ok=True)

    litpcba_dir = cache / "litpcba"
    tgz_path = cache / "AVE_unbiased.tgz"
    if not litpcba_dir.exists() or not any(litpcba_dir.iterdir()):
        if not tgz_path.exists():
            urllib.request.urlretrieve(LITPCBA_TGZ_URL, tgz_path)
        litpcba_dir.mkdir(parents=True, exist_ok=True)
        with tarfile.open(tgz_path, "r:gz") as tf:
            tf.extractall(litpcba_dir)  # noqa: S202 - trusted first-party academic mirror

    records: list[dict] = []
    for target in LITPCBA_TARGETS:
        target_dir = litpcba_dir / target
        if not target_dir.is_dir():
            continue
        records.extend(parse_litpcba_target(str(target_dir)))
    return records
