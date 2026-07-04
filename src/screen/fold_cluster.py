"""Fold-disjoint clustering + near-duplicate pocket removal for the Paper-2 orphan benchmark.

Cluster pockets by fold label (ECOD primary; Pfam via SIFTS as fallback). dedupe_pockets removes
library-internal near-duplicate sites (binding-site-residue Jaccard >= thr) so the >=8 clusters are
genuinely distinct folds. TM-align is NOT required here (residue-set Jaccard is the dedupe signal).
"""
from __future__ import annotations


def fold_of(pdb_id: str, chain: str, ecod_table: dict[tuple[str, str], str]) -> str | None:
    """Fold label for a (pdb, chain) from a preloaded ECOD table, or None if absent."""
    return ecod_table.get((pdb_id.upper(), chain))


def cluster_pockets(pocket_folds: dict[str, str]) -> dict[str, str]:
    """Map each pocket/target to its cluster id (= its fold label). Pockets sharing a fold share a
    cluster (they are NOT fold-disjoint from each other)."""
    return dict(pocket_folds)


def n_disjoint_clusters(clusters: dict[str, str]) -> int:
    """Number of distinct fold clusters."""
    return len(set(clusters.values()))


def dedupe_pockets(residue_sets: dict[str, set], jaccard_thr: float = 0.5) -> list[str]:
    """Return the pocket keys to DROP: for any pair whose binding-site residue sets have Jaccard
    >= thr, keep the first (sorted) and drop the rest. Deterministic."""
    keys = sorted(residue_sets)
    kept: list[str] = []
    dropped: list[str] = []
    for k in keys:
        s = residue_sets[k]
        redundant = False
        for kk in kept:
            a, b = s, residue_sets[kk]
            j = len(a & b) / len(a | b) if (a or b) else 0.0
            if j >= jaccard_thr:
                redundant = True
                break
        (dropped if redundant else kept).append(k)
    return dropped
