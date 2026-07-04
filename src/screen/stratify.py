"""Similarity stratification for the Paper-2 orphan benchmark (Increment 1).

A query's stratum similarity s = max Tanimoto (ECFP4, achiral) to the union of ALL pockets'
reference actives, EXCLUDING the query's own fingerprint. Cuts: high >=0.5, mid [0.35,0.5),
orphan [0.2,0.35), deep_orphan <0.2. No model, no scoring.
"""
from __future__ import annotations

import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator as fpg

_GEN = fpg.GetMorganGenerator(radius=2, fpSize=2048)   # useChirality defaults False -> achiral

HIGH, ORPHAN, DEEP = 0.5, 0.35, 0.2


def ecfp(smiles: str):
    """Achiral ECFP4 bit vector, or None if the SMILES does not parse."""
    mol = Chem.MolFromSmiles(smiles)
    return _GEN.GetFingerprint(mol) if mol is not None else None


def max_tanimoto(query_fp, active_fps: list) -> float:
    """Max Tanimoto of the query to a pool of active fingerprints (0.0 if the pool is empty).
    Does NOT self-exclude — the caller is responsible for removing the query from the pool
    (see `stratify`, which excludes by value)."""
    pool = [fp for fp in active_fps if fp is not None]
    if not pool:
        return 0.0
    return float(max(DataStructs.BulkTanimotoSimilarity(query_fp, pool)))


def assign_stratum(s: float) -> str:
    if s >= HIGH:
        return "high"
    if s >= ORPHAN:
        return "mid"
    if s >= DEEP:
        return "orphan"
    return "deep_orphan"


def stratify(queries: pd.DataFrame, pocket_actives: dict[str, list[str]]) -> pd.DataFrame:
    """Add columns `s` and `stratum` to `queries` (needs `smiles`). `pocket_actives` maps each
    pocket/target to its list of active SMILES; the query is compared to the union of all."""
    pool = [ecfp(smi) for smis in pocket_actives.values() for smi in smis]
    pool = [fp for fp in pool if fp is not None]
    out = queries.copy()
    s_vals, strata = [], []
    for smi in out["smiles"]:
        qfp = ecfp(smi)
        # Exclude the query's own fingerprint from the pool (self-exclusion)
        filtered_pool = [fp for fp in pool if fp is not None and not (fp == qfp)]
        s = 0.0 if qfp is None else max_tanimoto(qfp, filtered_pool)
        s_vals.append(s)
        strata.append(assign_stratum(s))
    out["s"] = s_vals
    out["stratum"] = strata
    return out
