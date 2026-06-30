"""Build the Boltz cage cross-check input manifest: cage ligand SMILES (4 forms),
anchor SMILES (PubChem), protein LBD sequences (UniProt), and the 20-complex matrix.

Run: PYTHONPATH=scripts python scripts/boltz_cage_inputs.py  (or `make boltzinputs`)
"""
from __future__ import annotations

import pathlib

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
