"""Grounding flagged QC systems against public ChEMBL affinity (for the close-the-loop race).

network ligand name -> SMILES (public name-keyed structures) -> ChEMBL single-assay pChEMBL ->
experimental dG. Results are cached to data/openfe_replicates/affinity_<system>.csv and committed;
the race reads the cache (deterministic). NO new MD.

network_ligands(system) source: OpenFE IndustryBenchmarks2024 input-structure SDFs
(https://github.com/OpenFreeEnergy/IndustryBenchmarks2024) — the SAME public benchmark the
combined-csv (data/openfe_replicates/combined_pymbar4_edge_data.csv) is built from. Verified by
direct join against that csv's `system name` + ligand_A/ligand_B columns:
  hif2a (merck/hif2a)                      -> 41/41 combined-csv ligands matched verbatim
  p38   (jacs_set/p38 + fragments/p38)     -> 40/40 matched verbatim
  bace  (jacs_set/bace)                    -> 36/36 matched verbatim
  cdk8  (merck/cdk8 + miscellaneous_set/cdk8) -> 29/35 verbatim; the remaining 6 are the SDF's
        "<n> flipped" alt-pose ligands vs the csv's "<n>-flipped" -> after space->hyphen
        normalization (documented, not fabricated -- same molecule, same SMILES), 35/35 matched.
No other renaming/normalization was needed for any system.
"""
from __future__ import annotations

import csv
import json
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path
from statistics import median

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")
RT_LN10 = 1.364  # kcal/mol at 298 K

# Verified human single-protein ChEMBL targets (Step 2), via
# https://www.ebi.ac.uk/chembl/api/data/target/search.json?q=<GENE> , cross-checked against
# https://www.ebi.ac.uk/chembl/api/data/target/<id>.json for organism=Homo sapiens,
# target_type=SINGLE PROTEIN, and (where ambiguous) the UniProt gene-symbol component.
TARGETS: dict[str, str] = {
    "hif2a": "CHEMBL1744522",  # EPAS1 / Endothelial PAS domain-containing protein 1 (verified)
    "cdk8": "CHEMBL5719",      # CDK8 / Cyclin-dependent kinase 8 (UniProt P49336, gene CDK8;
                               # search also returns CHEMBL6002 = CDK19, a paralog: excluded)
    "p38": "CHEMBL260",        # MAPK14 / Mitogen-activated protein kinase 14
    "bace": "CHEMBL4822",      # BACE1 / Beta-secretase 1
}

_GH_RAW = (
    "https://raw.githubusercontent.com/OpenFreeEnergy/IndustryBenchmarks2024/main/"
    "industry_benchmarks/input_structures/original_structures"
)

# system -> list of (sdf url suffix, name normalizer) pairs; a system's ligands may be split
# across more than one input SDF (e.g. cdk8 = merck set + a separate miscellaneous set).
_NAME_SOURCES: dict[str, list[str]] = {
    "hif2a": ["merck/hif2a/hif2a_automap_ligands.sdf"],
    "cdk8": [
        "merck/cdk8/cdk8_5cei_new_helix_loop_extra_ligands.sdf",
        "misc/cdk8_koehler_ligands.sdf",
    ],
    "p38": [
        "jacs_set/p38/p38_ligands.sdf",
        "fragments/p38/frag_p38_ligands.sdf",
    ],
    "bace": ["jacs_set/bace/bace_ligands.sdf"],
}

COMBINED_CSV = Path("data/openfe_replicates/combined_pymbar4_edge_data.csv")


def _network_ligand_names(system: str, csv_path: Path = COMBINED_CSV) -> set[str]:
    """The set of ligand names actually used in this system's combined-csv network (the
    denominator for coverage) -- read from `ligand_A`/`ligand_B` where `system name` == system."""
    names: set[str] = set()
    with csv_path.open(newline="") as f:
        for row in csv.DictReader(f):
            if row.get("system name") == system:
                names.add(row["ligand_A"])
                names.add(row["ligand_B"])
    return names


def system_dataset(system: str) -> str:
    """Combined-csv `system name` for a flag. Identity for cdk8/hif2a/p38/bace (no renaming in
    the combined csv itself -- the flag names ARE the `system name` values)."""
    return system


def _normalize_name(name: str) -> str:
    """cdk8's flipped-pose ligands: SDF title 'NN flipped' (space) vs combined-csv 'NN-flipped'
    (hyphen) -- same molecule, same SMILES, cosmetic naming difference only."""
    return name.strip().replace(" flipped", "-flipped")


def _fetch_sdf_text(url: str) -> str | None:
    try:
        with urllib.request.urlopen(url, timeout=60) as r:
            return r.read().decode("utf-8", errors="replace")
    except (urllib.error.URLError, TimeoutError, OSError):
        return None


def _parse_sdf_names_to_smiles(sdf_text: str) -> dict[str, str]:
    """Parse an SDF blob's per-record title line ('_Name') -> RDKit-canonical SMILES."""
    out: dict[str, str] = {}
    supp = Chem.SDMolSupplier()
    supp.SetData(sdf_text, removeHs=False)
    for mol in supp:
        if mol is None or not mol.HasProp("_Name"):
            continue
        name = _normalize_name(mol.GetProp("_Name"))
        if not name:
            continue
        out[name] = Chem.MolToSmiles(mol)
    return out


def network_ligands(system: str) -> dict[str, str]:
    """Combined-csv ligand name -> canonical SMILES, for one system's network.

    Source: public name-keyed OpenFE IndustryBenchmarks2024 input-structure SDFs (see module
    docstring for the verified per-system join coverage). The SDFs can carry ligands beyond the
    system's combined-csv network (e.g. alternate input sets); the result is restricted to names
    that actually appear in `data/openfe_replicates/combined_pymbar4_edge_data.csv` for this
    system, so `len(network_ligands(system))` is the correct coverage denominator. If a system
    has no known source or a fetch fails, returns whatever resolved -- never fabricates a SMILES
    for an unresolved name.
    """
    urls = _NAME_SOURCES.get(system, [])
    resolved: dict[str, str] = {}
    for suffix in urls:
        text = _fetch_sdf_text(f"{_GH_RAW}/{suffix}")
        if text is None:
            continue
        resolved.update(_parse_sdf_names_to_smiles(text))
    wanted = _network_ligand_names(system)
    if not wanted:
        return resolved
    return {name: smi for name, smi in resolved.items() if name in wanted}


def exp_dg(pchembl: float) -> float:
    return -RT_LN10 * float(pchembl)


def _coverage(mapped: dict, n_ligands: int) -> float:
    return len(mapped) / n_ligands if n_ligands else 0.0


def _median_pchembl(records: list, assay: str) -> float | None:
    vals = [float(r["pchembl_value"]) for r in records
            if r.get("standard_type") == assay and r.get("pchembl_value") not in (None, "", "None")]
    return float(median(vals)) if vals else None


def _get(url: str) -> dict:
    try:
        with urllib.request.urlopen(url, timeout=30) as r:
            return json.load(r)
    except Exception as e:  # network/parse failure -> treated as no match
        return {"_err": str(e)}


def _chembl_id_and_pchembl(smiles: str, target_chembl_id: str, assay: str
                            ) -> tuple[str | None, float | None]:
    """RDKit-canonical -> ChEMBL flexmatch -> (matched ChEMBL molecule id, median pChEMBL for
    `assay` against the target). `(None, None)` if no structure or no activity match."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, None
    q = urllib.parse.quote(Chem.MolToSmiles(mol))
    mol_j = _get(f"https://www.ebi.ac.uk/chembl/api/data/molecule.json"
                 f"?molecule_structures__canonical_smiles__flexmatch={q}&limit=1")
    ms = mol_j.get("molecules", [])
    if not ms:
        return None, None
    cid = ms[0]["molecule_chembl_id"]
    act = _get(f"https://www.ebi.ac.uk/chembl/api/data/activity.json"
               f"?molecule_chembl_id={cid}&target_chembl_id={target_chembl_id}"
               f"&pchembl_value__isnull=false&limit=50")
    p = _median_pchembl(act.get("activities", []), assay=assay)
    return cid, p


def chembl_affinity(smiles: str, target_chembl_id: str, assay: str) -> float | None:
    """RDKit-canonical -> ChEMBL flexmatch -> median-pChEMBL-derived experimental dG for `assay`
    against the target; `None` if no match."""
    _, p = _chembl_id_and_pchembl(smiles, target_chembl_id, assay)
    return exp_dg(p) if p is not None else None


def ground_system(system: str, prereg) -> dict:
    """Ground one system: name->SMILES->ChEMBL single-assay dG; cache + return coverage report."""
    lig = network_ligands(system)  # name -> canonical SMILES (public IndustryBenchmarks2024 SDFs)
    target = TARGETS.get(system, "")
    mapped: dict[str, float] = {}
    rows: list[dict] = []
    failures: list[str] = []
    assay_used = prereg.assay_order[0]
    for name, smi in lig.items():
        cid, p = _chembl_id_and_pchembl(smi, target, assay_used) if target else (None, None)
        if p is None:
            failures.append(name)
            continue
        dg = exp_dg(p)
        mapped[name] = dg
        rows.append({"ligand": name, "smiles": smi, "chembl_id": cid, "assay": assay_used,
                     "pchembl": p, "exp_dg": dg})
    out = Path(f"data/openfe_replicates/affinity_{system}.csv")
    out.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["ligand", "smiles", "chembl_id", "assay", "pchembl", "exp_dg"]
    with out.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)
    return {"map": mapped, "coverage": _coverage(mapped, len(lig)), "assay": assay_used,
            "n_ligands": len(lig), "failures": failures}


def grounded_systems(prereg) -> list:
    return [s for s in prereg.systems
            if ground_system(s, prereg)["coverage"] >= prereg.coverage_min]
