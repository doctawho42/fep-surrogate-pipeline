"""Amortized reward trunk: a minimal molecular featurizer + deep-ensemble predictor with a
full Fig-B decomposed sigma. The featurizer's even part (ECFP useChirality=False + achiral
descriptors) is enantiomer-blind; the parity-odd 0o channel (chiral.signed_volume on embedded
stereocentres) is the sole chirality carrier (invariant #5/#6).
Spec: docs/superpowers/specs/2026-06-30-trunk-amortized-reward-design.md.
"""
from __future__ import annotations

from functools import cache

import numpy as np
from numpy.typing import NDArray
from rdkit import Chem
from rdkit.Chem import AllChem, Crippen, Descriptors, rdMolDescriptors

from bar.calib import ensemble_predict, train_ensemble
from bar.chiral import signed_volume

_DESC = [Descriptors.MolWt, Crippen.MolLogP, rdMolDescriptors.CalcTPSA,
         rdMolDescriptors.CalcNumHBD, rdMolDescriptors.CalcNumHBA,
         rdMolDescriptors.CalcNumRotatableBonds]


@cache
def chir_pseudoscalar(smiles: str, seed: int = 1) -> float:
    """Parity-odd molecule descriptor: sum of signed volumes over assigned stereocentres of a
    3D embedding. Enantiomer (all centres inverted) -> negation; achiral -> 0."""
    m0 = Chem.MolFromSmiles(smiles)
    if m0 is None:
        raise ValueError(f"bad SMILES: {smiles}")
    mol = Chem.AddHs(m0)
    centres = Chem.FindMolChiralCenters(mol, includeUnassigned=False, useLegacyImplementation=False)
    if not centres:
        return 0.0
    pp = AllChem.ETKDGv3()
    pp.randomSeed = seed
    if AllChem.EmbedMolecule(mol, pp) != 0:
        return 0.0
    conf = mol.GetConformer()
    total = 0.0
    for idx, _label in centres:
        # [:4] captures all 4 neighbours of a stereocentre (incl. H after AddHs)
        nbrs = [a.GetIdx() for a in mol.GetAtomWithIdx(idx).GetNeighbors()][:4]
        if len(nbrs) < 4:
            continue
        coords = np.array([[conf.GetAtomPosition(j).x, conf.GetAtomPosition(j).y,
                            conf.GetAtomPosition(j).z] for j in nbrs])
        total += signed_volume(coords)
    return float(total)


@cache
def ligand_features(smiles: str) -> NDArray:
    """ECFP4 (useChirality=False, 1024 bits) + achiral descriptors. Enantiomer-blind."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"bad SMILES: {smiles}")
    fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024, useChirality=False)
    bits = np.zeros(1024, dtype=float)
    for b in fp.GetOnBits():
        bits[b] = 1.0
    desc = np.array([fn(mol) for fn in _DESC], dtype=float)
    return np.concatenate([bits, desc])


def featurize_edge(smiles_a: str, smiles_b: str, include_0o: bool = True) -> NDArray:
    """Delta-feature of a congeneric edge: f(B) - f(A), with the 0o pseudoscalar delta
    appended when include_0o (the sole chirality carrier)."""
    edge = ligand_features(smiles_b) - ligand_features(smiles_a)
    if include_0o:
        edge = np.append(edge, chir_pseudoscalar(smiles_b) - chir_pseudoscalar(smiles_a))
    return edge


class EnsembleTrunk:
    """Deep-ensemble amortized predictor over edge Delta-features -> (DDG_hat, sigma_epistemic)."""

    def __init__(self) -> None:
        self._nets: list = []
        self._include_0o = True

    def fit(self, edges_smiles: list[tuple[str, str]], ddg: NDArray, include_0o: bool = True,
            n_members: int = 8) -> EnsembleTrunk:
        self._include_0o = include_0o
        X = np.stack([featurize_edge(a, b, include_0o) for a, b in edges_smiles])
        self._nets = train_ensemble(X, np.asarray(ddg, dtype=float), n_members=n_members)
        return self

    def predict(self, edges_smiles: list[tuple[str, str]]) -> tuple[NDArray, NDArray]:
        X = np.stack([featurize_edge(a, b, self._include_0o) for a, b in edges_smiles])
        return ensemble_predict(self._nets, X)


def amortized_sigma(sigma_epistemic: NDArray, sigma_aleatoric: float | NDArray,
                    q: float) -> NDArray:
    """Conformal-scaled total sigma for the commit LCB: q * sqrt(sigma_epi^2 + sigma_ale^2)."""
    if q < 0:
        raise ValueError(f"conformal scale q must be >= 0, got {q}")
    se = np.asarray(sigma_epistemic, dtype=float)
    sa = np.asarray(sigma_aleatoric, dtype=float)
    return q * np.sqrt(se ** 2 + sa ** 2)
