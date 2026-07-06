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

KNOWN LIMITATION -- stereo-blind enantiomer collapse (disclosed, not fabrication):
`chembl_affinity` matches structures via ChEMBL's `flexmatch`, which is connectivity-based and
stereo-blind. Consequently enantiomer pairs present in a network (same 2D graph, opposite
stereocenter) resolve to the SAME ChEMBL molecule and therefore share one experimental affinity
value (empirically observed: hif2a 237/15, hif2a 165/164, hif2a 7b/7a; cdk8 30/31). The match
itself is real -- flexmatch legitimately found that ChEMBL record -- so this is NOT a fabricated
value. It DOES inflate the ABSOLUTE MUE-vs-experiment for enantiomer edges (one arm of the pair is
necessarily "wrong" relative to its true unmeasured enantiomer affinity), but it CANCELS in the
guided-vs-random close-the-loop CONTRAST, because both arms of that comparison read the identical
cached affinity for the identical ligand -- the shared value is a constant offset, not a
differential bias. Given this project's chirality invariants (#5 "predict per enantiomer", #6
"racemate apparent affinity ~= stronger binder"), the limitation must be disclosed rather than
silently absorbed. `collapsed_enantiomer_pairs(system)` reports the affected ligand-name pairs
for a given system's committed affinity CSV, for transparency in any accuracy claim built on it.
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


def collapsed_enantiomer_pairs(system: str) -> list[tuple[str, str]]:
    """Ligand-name pairs in the committed `data/openfe_replicates/affinity_<system>.csv` that
    share a ChEMBL flexmatch experimental affinity SOLELY because flexmatch is stereo-blind (see
    module docstring). Pure logic, no network: reads only the committed cache.

    Two ligands collapse iff (a) their stored `pchembl` values are identical, AND (b) their 2D
    graphs match (`Chem.MolToSmiles(mol, isomericSmiles=False)` equal for both), AND (c) their
    stored `smiles` strings differ (i.e. they are written as distinct stereoisomers, not the same
    string twice). Order-insensitive; each unordered pair reported once.
    """
    path = Path(f"data/openfe_replicates/affinity_{system}.csv")
    with path.open(newline="") as f:
        rows = list(csv.DictReader(f))

    by_2d: dict[str, list[tuple[str, str, float]]] = {}
    for row in rows:
        smi = row["smiles"]
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        key2d = Chem.MolToSmiles(mol, isomericSmiles=False)
        by_2d.setdefault(key2d, []).append((row["ligand"], smi, float(row["pchembl"])))

    pairs: list[tuple[str, str]] = []
    for group in by_2d.values():
        if len(group) < 2:
            continue
        for i in range(len(group)):
            name_i, smi_i, p_i = group[i]
            for j in range(i + 1, len(group)):
                name_j, smi_j, p_j = group[j]
                if p_i == p_j and smi_i != smi_j:
                    pairs.append((name_i, name_j))
    return pairs


def _map_ligands_for_assay(lig: dict[str, str], target: str, assay: str
                            ) -> dict[str, tuple[str | None, float]]:
    """name -> (chembl_id, pchembl) for every ligand that has a `pchembl` match for `assay`
    against `target`. Pure wrapper around `_chembl_id_and_pchembl`; the per-system assay
    fallback in `ground_system` calls this once per candidate assay."""
    out: dict[str, tuple[str | None, float]] = {}
    for name, smi in lig.items():
        cid, p = _chembl_id_and_pchembl(smi, target, assay) if target else (None, None)
        if p is not None:
            out[name] = (cid, p)
    return out


def select_assay(name_to_assay_pchembl: dict[str, dict[str, float]], assay_order: list[str]
                  ) -> str:
    """Per-system assay fallback (plan-mandated): pick the FIRST assay in `assay_order` that
    yields >=1 mapped ligand, given a `ligand name -> {assay -> pchembl}` mapping. Never mixes
    assays within a system -- the returned assay is used for every ligand of that system. Falls
    back to the last entry of `assay_order` if none has any match (so a system with zero matches
    for every assay still reports a definite `assay_used` rather than raising)."""
    for assay in assay_order:
        if any(assay in per_ligand for per_ligand in name_to_assay_pchembl.values()):
            return assay
    return assay_order[-1]


def ground_system(system: str, prereg) -> dict:
    """Ground one system: name->SMILES->ChEMBL dG; cache + return coverage report.

    Assay selection (Finding 2): tries each assay in `prereg.assay_order` in turn and uses the
    FIRST one that yields >=1 mapped ligand for this system (never mixes assays within a system).
    For the 4 currently-flagged systems IC50 always has matches, so `assay_used` stays "IC50" and
    the committed `data/openfe_replicates/affinity_*.csv` are unchanged by this fallback logic.
    """
    lig = network_ligands(system)  # name -> canonical SMILES (public IndustryBenchmarks2024 SDFs)
    target = TARGETS.get(system, "")
    per_assay: dict[str, dict[str, tuple[str | None, float]]] = {
        assay: _map_ligands_for_assay(lig, target, assay) for assay in prereg.assay_order
    }
    name_to_assay_pchembl = {
        name: {assay: p for assay, m in per_assay.items() if name in m for _, p in [m[name]]}
        for name in lig
    }
    assay_used = select_assay(name_to_assay_pchembl, prereg.assay_order)
    chosen = per_assay[assay_used]

    mapped: dict[str, float] = {}
    rows: list[dict] = []
    failures: list[str] = []
    for name, smi in lig.items():
        if name not in chosen:
            failures.append(name)
            continue
        cid, p = chosen[name]
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
