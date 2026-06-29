"""The cage molecule (NIOCH target) — structure preparation, per enantiomer.

Synthesis (see docs): a (substituted) difluoronaphthalenone + a 1,3-dicarbonyl Michael
donor (here barbituric acid) cyclise; the product is O-acylated (here acetate). The two
stereocentres are coupled by the ring fusion, so the product is a single diastereomer;
we prepare it AND its mirror image (invariant #6: predict per enantiomer; the racemate's
apparent affinity is the stronger binder).

Note: in the product the barbituric acid has cyclised into a fused
uracil/pyrimidine-2,4-dione (`[nH]C(=O)[nH]C(=O)` diamide + enol-ether bridge), not a
free barbiturate. That cyclic-ureide diamide is a strong, promiscuous H-bond donor/acceptor
array and will dominate the complementarity signal -- the reason calibration is a guard
against a false promiscuous lead, not decoration (docs/target_finding_plan.md §1).
"""
from __future__ import annotations

import pathlib

# Reference reaction components (vary-able for the analogue library).
PRODUCT_SMILES = "CC(=O)O[C@]12C[C@H](C3=CC=CC=C3C1(F)F)C4=C(NC(=O)NC4=O)O2"
SYNTHONS = {
    "difluoronaphthalenone": "C1=CC=C2C(=C1)C=CC(=O)C2(F)F",
    "michael_donor_barbituric_acid": "O=C1CC(=O)NC(=O)N1",
    "acylating_acetic_acid": "CC(=O)O",
}
UREIDE_SMARTS = "[#7][#6](=O)[#7][#6](=O)"  # cyclic-ureide / uracil diamide (arom or aliph)


def enantiomer_smiles(smiles: str) -> str:
    """Mirror image: invert every tetrahedral stereocentre."""
    from rdkit import Chem

    m = Chem.MolFromSmiles(smiles)
    for a in m.GetAtoms():
        t = a.GetChiralTag()
        if t == Chem.ChiralType.CHI_TETRAHEDRAL_CW:
            a.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CCW)
        elif t == Chem.ChiralType.CHI_TETRAHEDRAL_CCW:
            a.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CW)
    return Chem.MolToSmiles(m)


def embed_3d(smiles: str, out_sdf: str, seed: int = 1) -> str:
    """ETKDGv3 embed + MMFF optimise; write a single-conformer SDF."""
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    AllChem.EmbedMolecule(mol, params)
    AllChem.MMFFOptimizeMolecule(mol)
    pathlib.Path(out_sdf).parent.mkdir(parents=True, exist_ok=True)
    w = Chem.SDWriter(out_sdf)
    w.write(mol)
    w.close()
    return out_sdf


def prepare(out_dir: str = "data/cage") -> dict:
    """Write the given enantiomer and its mirror as 3D SDFs."""
    enant = enantiomer_smiles(PRODUCT_SMILES)
    return {
        "given": embed_3d(PRODUCT_SMILES, f"{out_dir}/cage_given.sdf"),
        "enantiomer": embed_3d(enant, f"{out_dir}/cage_enantiomer.sdf"),
        "given_smiles": PRODUCT_SMILES,
        "enantiomer_smiles": enant,
    }


if __name__ == "__main__":
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdMolDescriptors

    m = Chem.MolFromSmiles(PRODUCT_SMILES)
    print("cage:", rdMolDescriptors.CalcMolFormula(m), f"MW {Descriptors.MolWt(m):.1f}",
          f"logP {Descriptors.MolLogP(m):.2f} RotB {Descriptors.NumRotatableBonds(m)}")
    print("stereocentres:", Chem.FindMolChiralCenters(m, useLegacyImplementation=False))
    print("ureide motif present:", m.HasSubstructMatch(Chem.MolFromSmarts(UREIDE_SMARTS)))
    out = prepare()
    print("enantiomer SMILES:", out["enantiomer_smiles"])
    print("wrote:", out["given"], out["enantiomer"])
