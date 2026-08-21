"""Fig M — Paper-2 Increment-1 orphan-benchmark validity gate (the go/no-go).

Assembles the benchmark (bench_sources) -> fold labels (bench_sources.build_triples's `fold`
column) + disjoint-cluster count (fold_cluster.n_disjoint_clusters) -> similarity strata
(stratify) -> Step-0 shape-null validity verdict (validity_gate). Prints PASS or VALIDITY_KILL
and draws the strata / shape-null-vs-random figure. NO structure scoring. Deterministic.

Run:  python figs/make_figM.py  (or `make figM`). Triples are cached in data/paper2_bench/
(built once by `build_triples`); re-runs are fast (no LIT-PCBA re-download).
"""
from __future__ import annotations

import pathlib
import sys

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from screen.bench_sources import build_triples  # noqa: E402
from screen.fold_cluster import n_disjoint_clusters  # noqa: E402
from screen.stratify import stratify  # noqa: E402
from screen.validity_gate import shape_score_matrix, verdict  # noqa: E402

FIGDIR = ROOT / "figs"


def build_pocket_actives(triples):
    return {t: g["smiles"].tolist() for t, g in triples.groupby("target")}


def main() -> None:
    triples = build_triples()
    pocket_actives = build_pocket_actives(triples)
    pocket_order = sorted(pocket_actives)
    # one query row per ligand (mol_id is already unique per row in the cached triples);
    # true target = its `target`.
    queries = triples[["mol_id", "smiles", "target"]].drop_duplicates("mol_id").reset_index(drop=True)
    queries = stratify(queries, pocket_actives)
    scores, true_idx = shape_score_matrix(queries, pocket_actives, pocket_order)

    # fold clusters: one fold label per target, from the `fold` column build_triples attaches
    # (RCSB Pfam/InterPro lookup on the representative pocket PDB; falls back to the target
    # name itself when no structure is available). n_fold_clusters = number of DISTINCT fold
    # labels among the true-target pockets of the ORPHAN-stratum queries.
    fold_by_target = triples.drop_duplicates("target").set_index("target")["fold"].to_dict()
    orphan_targets = set(queries.loc[queries["stratum"] == "orphan", "target"])
    clusters = {t: fold_by_target[t] for t in orphan_targets}
    v = verdict(queries, scores, true_idx, n_pockets=len(pocket_order),
                n_fold_clusters=n_disjoint_clusters(clusters))

    print(f"[Fig M] {v['verdict']}  | pockets={v['n_pockets']} fold-clusters={v['n_fold_clusters']}")
    print(f"  orphan: n={v['orphan']['n']} recovery@1={v['orphan']['recovery1']:.3f} "
          f"AUROC={v['orphan']['auroc']:.3f} CI={v['orphan']['ci']}")
    print(f"  high  : n={v['high']['n']} recovery@1={v['high']['recovery1']:.3f}")
    for r in v["reasons"]:
        print("  reason:", r)

    strata_order = ["high", "mid", "orphan", "deep_orphan"]
    counts = [int((queries["stratum"] == s).sum()) for s in strata_order]
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(8.2, 3.4))
    axL.bar(strata_order, counts, color="#4C72B0")
    axL.set_title("A  query counts by similarity stratum", loc="left", fontweight="bold", fontsize=9)
    axL.set_ylabel("queries")
    rec = [v["high"]["recovery1"], v["orphan"]["recovery1"]]
    rnd = 1.0 / max(v["n_pockets"], 1)
    axR.bar(["high", "orphan"], rec, color=["#C44E52", "#55A868"], label="shape-null recovery@1")
    axR.axhline(rnd, ls="--", c="k", lw=1, label=f"random 1/N={rnd:.3f}")
    axR.set_title(f"B  shape-null collapse -> {v['verdict']}", loc="left", fontweight="bold", fontsize=9)
    axR.legend(frameon=False, fontsize=7)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(FIGDIR / f"figM_orphan_validity.{ext}", dpi=300, bbox_inches="tight")
    print("wrote figs/figM_orphan_validity.(pdf|png)")


if __name__ == "__main__":
    main()
