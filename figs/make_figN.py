"""Fig N — Paper-2 Increment-2 Step-0 validity gate on the cross-target-confusable benchmark.

Aggregate library (LIT-PCBA + ChEMBL-diverse + BindingDB, the fullest breadth attempt) ->
stratify -> the s<0.15 collapse stratum -> does the shape-null collapse there, well-powered
+ monotone gradient? Prints PASS (green-light Plan-2 smina gate) or terminal C. NO scoring.

Run:  python figs/make_figN.py  (or `make figN`). Triples are cached in data/paper2_bench/
(built once by `build_triples(include_bindingdb=True)`, ~27.5k ligands / 44 targets); re-runs
are fast (no re-download). The similarity stratification over ~27k ligands is the slow step
(~minutes) on a fresh cache.
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
from screen.recovery import recovery_at_k  # noqa: E402
from screen.stratify import stratify  # noqa: E402
from screen.validity_gate import shape_score_matrix, verdict  # noqa: E402

FIGDIR = ROOT / "figs"

STRATA = ["high", "mid", "orphan", "collapse"]
COLLAPSE_TAU = 0.15


def is_monotone_gradient(rec_by_stratum: dict[str, float]) -> bool:
    """recovery@1 must not increase as similarity decreases (high>=mid>=orphan>=collapse)."""
    seq = [rec_by_stratum[s] for s in STRATA if s in rec_by_stratum]
    return all(a + 1e-9 >= b for a, b in zip(seq, seq[1:], strict=False))


def _collapse_label(s: float) -> str:
    if s >= 0.50:
        return "high"
    if s >= 0.35:
        return "mid"
    if s >= COLLAPSE_TAU:
        return "orphan"
    return "collapse"


def main() -> None:
    # fullest aggregate: LIT-PCBA + ChEMBL-diverse + BindingDB
    df = build_triples(include_bindingdb=True)
    pocket_actives = {t: g["smiles"].tolist() for t, g in df.groupby("target")}
    pocket_order = sorted(pocket_actives)
    queries = df[["mol_id", "smiles", "target"]].drop_duplicates("mol_id").reset_index(drop=True)
    queries = stratify(queries, pocket_actives)                     # adds leave-one-out s
    queries["stratum"] = queries["s"].map(_collapse_label)
    scores, true_idx = shape_score_matrix(queries, pocket_actives, pocket_order)

    # per-stratum recovery for the gradient; the gate itself runs on the collapse stratum
    fold_by_target = df.drop_duplicates("target").set_index("target")["fold"]
    rec = {}
    for st in STRATA:
        mask = (queries["stratum"] == st).to_numpy()
        if mask.sum() == 0:
            continue
        rec[st] = float(recovery_at_k(scores[mask], true_idx[mask], lower_better=False)[0])

    collapse_targets = set(queries.loc[queries["stratum"] == "collapse", "target"])
    n_folds = n_disjoint_clusters({t: fold_by_target[t] for t in collapse_targets})
    # Reuse the amended-P2 verdict on the collapse stratum: verdict() reads whichever rows are
    # labeled "orphan". _collapse_label ALSO produces a stratum literally named "orphan"
    # (0.15<=s<0.35), so naively relabeling collapse->"orphan" would silently MERGE the two
    # strata into one group (782 rows instead of the true collapse n). Relabel the pre-existing
    # "orphan" rows out of the way first so only the true collapse stratum lands under "orphan".
    qv = queries.copy()
    qv["stratum"] = qv["stratum"].map(lambda s: "_mid_orphan" if s == "orphan" else s)
    qv["stratum"] = qv["stratum"].map(lambda s: "orphan" if s == "collapse" else s)
    v = verdict(qv, scores, true_idx, n_pockets=len(pocket_order), n_fold_clusters=n_folds)

    mono = is_monotone_gradient(rec)
    step0_pass = v["verdict"] == "PASS" and mono
    verdict_str = "PASS" if step0_pass else "FAIL/TERMINAL-C"
    print(f"[Fig N Step-0] {verdict_str} | pockets={len(pocket_order)} "
          f"collapse folds={n_folds} monotone={mono}")
    grad_str = ", ".join(f"{s}={rec.get(s, float('nan')):.3f}" for s in STRATA)
    print(f"  gradient recovery@1: {grad_str}")
    print(f"  collapse stratum: n={v['orphan']['n']} recovery@1={v['orphan']['recovery1']:.3f} "
          f"AUROC={v['orphan']['auroc']:.3f} CI={v['orphan']['ci']}")
    for r in v["reasons"]:
        print("  reason:", r)

    fig, ax = plt.subplots(figsize=(5.2, 3.4))
    xs = [s for s in STRATA if s in rec]
    ax.plot(xs, [rec[s] for s in xs], "o-", color="#4C72B0")
    ax.axhline(1.0 / len(pocket_order), ls="--", c="k", lw=1, label="random 1/N")
    ax.set_ylabel("shape-null recovery@1")
    ax.set_ylim(0, 1)
    ax.set_title(f"Fig N Step-0: {'PASS' if step0_pass else 'terminal C'} "
                 f"(collapse n={v['orphan']['n']}, folds={n_folds})", fontsize=9, loc="left")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(FIGDIR / f"figN_collapse_validity.{ext}", dpi=300, bbox_inches="tight")
    print("wrote figs/figN_collapse_validity.(pdf|png)")


if __name__ == "__main__":
    main()
