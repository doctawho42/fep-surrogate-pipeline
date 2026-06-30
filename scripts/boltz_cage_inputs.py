"""Build the Boltz cage cross-check input manifest: cage ligand SMILES (4 forms),
anchor SMILES (PubChem), protein LBD sequences (UniProt), and the 20-complex matrix.

Run: PYTHONPATH=scripts python scripts/boltz_cage_inputs.py  (or `make boltzinputs`)
"""
from __future__ import annotations

import json
import pathlib
import urllib.request

from rdkit import Chem
from rdkit.Chem import AllChem

ROOT = pathlib.Path(__file__).resolve().parent.parent
CAGE = ROOT / "data" / "cage"

# Deacetylation: hydrolyse a single O-acetyl ester (R-O-C(=O)CH3 -> R-OH). Matches both
# alcohol and phenol acetates; the O keeps its R bond and gains an implicit H.
_DEACETYL = AllChem.ReactionFromSmarts("[#6:1][O:2][CX3](=O)[CH3]>>[#6:1][O:2]")


def _sdf_smiles(path: pathlib.Path) -> str:
    """Isomeric SMILES from a 3D SDF, with stereo assigned from the coordinates."""
    mol = next(Chem.SDMolSupplier(str(path), removeHs=False))
    if mol is None:
        raise ValueError(f"could not read {path}")
    Chem.AssignStereochemistryFrom3D(mol)
    return Chem.MolToSmiles(Chem.RemoveHs(mol))


def deacetylate(smiles: str) -> str:
    """Remove one O-acetyl ester, returning the free alcohol/phenol SMILES."""
    mol = Chem.MolFromSmiles(smiles)
    prods = _DEACETYL.RunReactants((mol,))
    if not prods:
        return smiles  # no acetate present
    out = prods[0][0]
    Chem.SanitizeMol(out)
    return Chem.MolToSmiles(out)


def cage_smiles() -> dict[str, str]:
    rr = _sdf_smiles(CAGE / "cage_given.sdf")  # R,R given
    ss = _sdf_smiles(CAGE / "cage_enantiomer.sdf")  # S,S mirror
    return {"RR_OAc": rr, "SS_OAc": ss,
            "RR_OH": deacetylate(rr), "SS_OH": deacetylate(ss)}


# Target order for the panel (AChE = greasy-pocket control; docking falsely scored it 100%).
_TARGETS = ("GR", "AR", "ER", "DHODH", "AChE")
# UniProt accession + 1-based LBD / catalytic-domain ranges (confirm vs UniProt features).
_LBD = {"GR": ("P04150", 521, 777), "AR": ("P10275", 669, 919),
        "ER": ("P03372", 305, 550), "DHODH": ("Q02127", 78, 395),
        "AChE": ("P22303", 32, 574)}

# Per-target anchor (GR/AR/ER agonists; DHODH inhibitor), as curated isomeric SMILES from
# ChEMBL (PubChem is unreachable from this environment). Source ChEMBL IDs are recorded so the
# provenance is explicit and the values are verifiable; each is RDKit-parse-checked in tests.
_ANCHOR_SMILES = {
    # dexamethasone, GR agonist (CHEMBL384467)
    "GR": "C[C@@H]1C[C@H]2[C@@H]3CCC4=CC(=O)C=C[C@]4(C)[C@@]3(F)[C@@H](O)C[C@]2(C)[C@@]1(O)C(=O)CO",
    # 5alpha-dihydrotestosterone / stanolone, AR agonist (CHEMBL27769)
    "AR": "C[C@]12CCC(=O)C[C@@H]1CC[C@@H]1[C@@H]2CC[C@]2(C)[C@@H](O)CC[C@@H]12",
    # 17beta-estradiol, ER agonist (CHEMBL135)
    "ER": "C[C@]12CC[C@@H]3c4ccc(O)cc4CC[C@H]3[C@@H]1CC[C@@H]2O",
    # brequinar, DHODH inhibitor (CHEMBL38434)
    "DHODH": "Cc1c(-c2ccc(-c3ccccc3F)cc2)nc2ccc(F)cc2c1C(=O)O",
    # donepezil, AChE inhibitor (CHEMBL502)
    "AChE": "COc1cc2c(cc1OC)C(=O)C(CC1CCN(Cc3ccccc3)CC1)C2",
}


def _uniprot_fasta(acc: str) -> str:
    url = f"https://rest.uniprot.org/uniprotkb/{acc}.fasta"
    raw = urllib.request.urlopen(url, timeout=30).read().decode()
    return "".join(line for line in raw.splitlines() if not line.startswith(">"))


def lbd_sequences() -> dict[str, str]:
    return {k: _uniprot_fasta(acc)[lo - 1:hi] for k, (acc, lo, hi) in _LBD.items()}


def anchor_smiles() -> dict[str, str]:
    return dict(_ANCHOR_SMILES)


def build_manifest() -> dict:
    """The 20-complex matrix: each target x {its anchor + the 4 cage forms}."""
    targets, anchors, cage = lbd_sequences(), anchor_smiles(), cage_smiles()
    complexes = []
    for tgt in _TARGETS:
        complexes.append({"name": f"cage-{tgt}-anchor", "target": tgt,
                          "ligand_label": "anchor", "ligand_smiles": anchors[tgt],
                          "is_anchor": True})
        for form, smi in cage.items():
            complexes.append({"name": f"cage-{tgt}-{form}", "target": tgt,
                              "ligand_label": form, "ligand_smiles": smi,
                              "is_anchor": False})
    return {"targets": targets, "anchors": anchors, "cage": cage, "complexes": complexes}


def screen_payload(target: str, m: dict) -> dict:
    """One small-molecule-screen payload: the {anchor + 4 cage forms} library docked into
    `target`, with the agonist seeding pocket detection. Filters disabled so no library
    molecule is dropped before scoring (we need the cage scored)."""
    anchor = m["anchors"][target]
    molecules = [{"smiles": anchor, "id": "anchor"}]
    molecules += [{"smiles": smi, "id": form} for form, smi in m["cage"].items()]
    return {
        "molecules": molecules,
        "target": {
            "entities": [{"type": "protein", "chain_ids": ["A"], "value": m["targets"][target]}],
            "reference_ligands": [anchor],
        },
        "molecule_filters": {"boltz_smarts_catalog_filter_level": "disabled"},
    }


def write_screen_payloads(m: dict) -> list[pathlib.Path]:
    out_dir = CAGE / "boltz_screen_payloads"
    out_dir.mkdir(exist_ok=True)
    paths = []
    for tgt in _TARGETS:
        p = out_dir / f"cage-screen-{tgt}.json"
        p.write_text(json.dumps(screen_payload(tgt, m), indent=2))
        paths.append(p)
    return paths


def main() -> None:
    m = build_manifest()
    out = CAGE / "boltz_inputs.json"
    out.write_text(json.dumps(m, indent=2))
    paths = write_screen_payloads(m)
    print(f"wrote {out} ({len(m['complexes'])} complexes) + {len(paths)} screen payloads")


if __name__ == "__main__":
    main()
