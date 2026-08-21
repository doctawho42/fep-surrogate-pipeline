"""Fig Lcausal: does acting on the calibrated QC flag improve accuracy vs experiment?

Guided (remove top-|z| QC edges) vs random removal, MUE-vs-experiment, over grounded flagged
systems. Run: PYTHONPATH=src python figs/make_figLcausal.py   (or `make figLcausal`)
"""
from __future__ import annotations

import csv
import pathlib
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
from paperstyle import (  # noqa: E402
    INK,
    MUTED,
    NARROW,
    OURS,
    finish,
    legend,
    panel,
    use_paper_style,
)

from bar.closeloop import combine, load_prereg, load_system_edges, system_effect  # noqa: E402


def _exp(system: str) -> dict:
    p = ROOT / "data" / "openfe_replicates" / f"affinity_{system}.csv"
    if not p.exists():
        return {}
    return {r["ligand"]: float(r["exp_dg"]) for r in csv.DictReader(p.open())}


def main() -> None:
    pr = load_prereg()
    # a system is grounded iff Task 2 wrote its affinity cache (the coverage gate is enforced there)
    grounded = [
        s for s in pr.systems
        if (ROOT / "data" / "openfe_replicates" / f"affinity_{s}.csv").exists() and len(_exp(s))
    ]
    effects = []
    for s in grounded:
        exp = _exp(s)
        edges = load_system_edges(s, exp)
        eff = system_effect(edges, exp, n_perm=pr.n_perm, target_rchi2=pr.target_reduced_chi2)
        eff["system"] = s
        effects.append(eff)
        print(f"{s}: K={eff['k']} guided ΔMUE={eff['guided']:.3f} "
              f"random={eff['random_mean']:.3f} p={eff['p']:.3f}")
    verdict = combine(effects) if len(effects) >= pr.min_grounded else {"verdict": "ABORT",
                                                                         "stouffer_p": float("nan")}
    stouffer_p = verdict.get("stouffer_p", float("nan"))
    # The pre-registration records the outcome as a code word; the printed figure spells it out,
    # since a reader of the article has no key to the code.
    said = {"SUCCESS": "guided removal beats random removal",
            "NULL": "guided removal ties random removal",
            "ABORT": "too few grounded systems to compare"}.get(verdict["verdict"],
                                                                verdict["verdict"].lower())

    # --- presentation -------------------------------------------------------------------
    # One quantity family, one colour: the |z|-guided removal rule is built from the
    # calibrated per-edge residual, so it is OURS; the size-matched random draw is a null
    # baseline, so it is MUTED. This is the same pairing as Fig Hodge panel C.
    use_paper_style()
    fig, ax = plt.subplots(figsize=(NARROW, 3.1))
    xs = list(range(len(effects)))
    guided = [e["guided"] for e in effects]
    random_mean = [e["random_mean"] for e in effects]
    bars = ax.bar(xs, guided, width=0.62, color=OURS, lw=0, zorder=2, label="guided removal")
    pts = ax.scatter(xs, random_mean, s=26, color=MUTED, lw=0, zorder=3,
                     label="random removal (mean)")
    ax.axhline(0, color=INK, lw=0.7, zorder=1)
    lo, hi = min(guided + random_mean), max(guided + random_mean)
    span = hi - lo
    ax.set_ylim(lo - 0.12 * span, hi + 0.34 * span)
    ax.set_xlim(-0.62, len(effects) - 0.38)
    ax.set_xticks(xs)
    ax.set_xticklabels([e["system"] for e in effects])
    ax.tick_params(axis="x", length=0)
    ax.set_ylabel("$\\Delta$MUE vs experiment\n(positive is an improvement)")
    legend(ax, loc="upper left", handles=[bars, pts])
    panel(ax, "", said, subtitle=f"close-the-loop race, Stouffer $p={stouffer_p:.3f}$")

    finish(fig, "figLcausal_guided_vs_random")
    print(f"\nVERDICT: {verdict}")


if __name__ == "__main__":
    main()
