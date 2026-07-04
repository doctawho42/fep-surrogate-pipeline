"""Diverse ChEMBL targets -> active-ligand records for the Paper-2 cross-target-confusable
benchmark. Adds chemotype-diverse targets (each with a holo PDB + real co-crystal ligand box +
Pfam fold) so more query ligands are near-zero similar to ALL actives. NO structure scoring here.

Schema-probe (Task-2 Step 1) findings: RCSB's entry-level `rcsb_entry_info` does NOT expose
non-polymer comp_ids directly; the reliable fields are `rcsb_binding_affinity[].comp_id` (when
present) and each `nonpolymer_entity/<pdb>/<entity_id>.pdbx_entity_nonpoly.comp_id`. Six of the
brief's ten seeded `lig_resname` values did not match their PDB's actual bound ligand and were
corrected below (CDK2, HIVPR, CAII, FXA); four targets (HSP90, MMP9, BACE1, DHFR) were added to
push distinct-Pfam-fold coverage safely past the >=8 gate once the known FXA/THROMB collapse
(both PF00089, trypsin-like serine proteases) and a second, brief-unanticipated collapse
(CDK2/P38, both PF00069 "Pkinase") are accounted for. A candidate ADRB2 GPCR entry (2RH1) was
tried and dropped: `_pfam_for_pdb` returns the FIRST Pfam hit per polymer entity, and 2RH1's
receptor is expressed as a T4-lysozyme fusion whose entity 1 lists the lysozyme's Pfam
(PF00959) before the receptor's own 7TM Pfam (PF00001) -- using it would silently mislabel the
fold. See docs / task-2-report.md for the full per-target verification table.
"""
from __future__ import annotations

import json
import urllib.request

# Verified in Task-2 Step 1 (RCSB nonpolymer_entity + rcsb_binding_affinity probes, 2026-07-05):
# each pdb_id is a holo structure of that ChEMBL target with a real drug-like co-crystal ligand
# matching lig_resname. 14 targets -> 12 distinct Pfam/InterPro fold labels (>= 8 gate, with
# margin for the two known kinase/protease-fold collapses).
DIVERSE_TARGETS: list[dict] = [
    {"target": "EGFR",   "chembl_id": "CHEMBL203",  "pdb_id": "1M17", "lig_resname": "AQ4"},
    {"target": "CDK2",   "chembl_id": "CHEMBL301",  "pdb_id": "1H1S", "lig_resname": "4SP"},
    {"target": "HIVPR",  "chembl_id": "CHEMBL243",  "pdb_id": "1HXW", "lig_resname": "RIT"},
    {"target": "ACHE",   "chembl_id": "CHEMBL220",  "pdb_id": "1EVE", "lig_resname": "E20"},
    {"target": "CAII",   "chembl_id": "CHEMBL205",  "pdb_id": "1BN1", "lig_resname": "AL5"},
    {"target": "FXA",    "chembl_id": "CHEMBL244",  "pdb_id": "1F0R", "lig_resname": "815"},
    {"target": "THROMB", "chembl_id": "CHEMBL204",  "pdb_id": "1DWD", "lig_resname": "MID"},
    {"target": "PDE5",   "chembl_id": "CHEMBL1827", "pdb_id": "1UDT", "lig_resname": "VIA"},
    {"target": "P38",    "chembl_id": "CHEMBL260",  "pdb_id": "1A9U", "lig_resname": "SB2"},
    {"target": "GR",     "chembl_id": "CHEMBL2034", "pdb_id": "1M2Z", "lig_resname": "DEX"},
    {"target": "HSP90",  "chembl_id": "CHEMBL3880", "pdb_id": "1YET", "lig_resname": "GDM"},
    {"target": "MMP9",   "chembl_id": "CHEMBL321",  "pdb_id": "1GKC", "lig_resname": "NFH"},
    {"target": "BACE1",  "chembl_id": "CHEMBL4822", "pdb_id": "2QMG", "lig_resname": "SC6"},
    {"target": "DHFR",   "chembl_id": "CHEMBL202",  "pdb_id": "1DRF", "lig_resname": "FOL"},
]

TRIPLE_KEYS = [
    "mol_id", "smiles", "target", "pdb_id", "lig_resname", "affinity_nm", "source", "fold",
]


def records_from_activities(
    acts: list[dict], target: str, pdb_id: str, lig_resname: str
) -> list[dict]:
    """Shape ChEMBL activity rows into benchmark records; drop empty/dup molecules."""
    seen, recs = set(), []
    for a in acts:
        mid, smi = a.get("molecule_chembl_id"), a.get("canonical_smiles")
        if not mid or not smi:
            continue
        key = f"{target}:{mid}"
        if key in seen:
            continue
        seen.add(key)
        recs.append({"mol_id": key, "smiles": smi, "target": target, "pdb_id": pdb_id,
                     "lig_resname": lig_resname, "affinity_nm": None, "source": "chembl_diverse",
                     "fold": None})
    return recs


def fetch_actives(chembl_id: str, min_pchembl: float = 6.5, max_records: int = 1500) -> list[dict]:
    """ChEMBL EBI REST actives for a target (pattern from figs/make_figF.py:fetch_target)."""
    rows, offset = [], 0
    while offset < max_records:
        url = (f"https://www.ebi.ac.uk/chembl/api/data/activity.json?target_chembl_id={chembl_id}"
               f"&pchembl_value__gte={min_pchembl}&limit=1000&offset={offset}")
        req = urllib.request.Request(url, headers={"User-Agent": "paper2/1.0"})
        d = json.load(urllib.request.urlopen(req, timeout=60))
        acts = d.get("activities", [])
        rows.extend(acts)
        if len(acts) < 1000:
            break
        offset += 1000
    return rows


def load_chembl_diverse_records(cache_dir: str = "data/paper2_bench") -> list[dict]:
    """All diverse-target active records (fold left None; the aggregator assigns Pfam folds)."""
    out: list[dict] = []
    for t in DIVERSE_TARGETS:
        acts = fetch_actives(t["chembl_id"])
        out.extend(records_from_activities(acts, t["target"], t["pdb_id"], t["lig_resname"]))
    return out
