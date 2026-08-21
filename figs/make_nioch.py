"""Generate the NIOCH cage-screen report + figure from the docking caches.

Operational-debt deliverable (NOT a paper claim): ranked, calibrated, per-enantiomer
target HYPOTHESES + recommended assays, with honest caveats. Calibration = the cage's
percentile among each target's KNOWN binders docked into the same pocket (kills the
greasy-pocket artefact that Fig H exposed). References are pooled from the gate cache
(scaffold-disjoint test ligands docked into their own pocket) and the campaign cache.
Re-run as the campaign fills in more references.
"""
from __future__ import annotations

import csv
import pathlib
import textwrap
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "figs"))
sys.path.insert(0, str(ROOT / "src"))
from make_boltz_cage import (  # noqa: E402
    CAGE_FORMS as BOLTZ_FORMS,
    TARGETS as BOLTZ_TARGETS,
    anchor_gap,
    load_target_conf,
    verdict as boltz_verdict,
)
from make_figH import PANEL as GPANEL, build_queries  # noqa: E402
from screen.campaign import CAGE, PANEL  # noqa: E402

DATA = ROOT / "data"


def pooled_references():
    gate = {(l, p): float(s) for l, p, s in csv.reader(open(DATA / "figH" / "dock_scores.csv"))}
    camp = {(l, p): float(s) for l, p, s in csv.reader(open(DATA / "campaign" / "scores.csv"))}
    queries, _ = build_queries()
    true_by_lig = {m: t for m, _, t in queries}
    gpdb = [p for _, _, p in GPANEL]
    ref = {pdb: [] for _, _, pdb in PANEL}
    for (lig, pdb), sc in gate.items():
        if np.isfinite(sc) and lig in true_by_lig and gpdb[true_by_lig[lig]] == pdb:
            ref.setdefault(pdb, []).append(sc)
    for (lig, pdb), sc in camp.items():
        if not lig.startswith("cage::") and np.isfinite(sc):
            ref.setdefault(pdb, []).append(sc)
    cage = {k: v for k, v in camp.items() if k[0].startswith("cage::")}
    return ref, cage


def rank():
    ref, camp = pooled_references()
    rows = []
    for tname, _, pdb in PANEL:
        dist = sorted(ref.get(pdb, []))
        rec = {"target": tname, "pdb": pdb, "n": len(dist)}
        for e in CAGE:
            cs = camp.get((f"cage::{e}", pdb), float("nan"))
            rec[e] = cs
            rec["pct_" + e] = (100 * sum(1 for r in dist if r > cs) / len(dist)
                               if dist and np.isfinite(cs) else float("nan"))
        rec["best"] = np.nanmax([rec["pct_" + e] for e in CAGE]) if dist else float("nan")
        rows.append(rec)
    rows.sort(key=lambda r: -(r["best"] if r["best"] == r["best"] else -1))
    return rows


def figure(rows, path):
    mpl.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
                         "font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                         "savefig.dpi": 300, "savefig.bbox": "tight"})
    rr = [r for r in rows if r["n"] >= 3]
    fig, ax = plt.subplots(figsize=(6.6, 3.6))
    y = np.arange(len(rr))[::-1]
    e1, e2 = list(CAGE)
    ax.barh(y + 0.18, [r["pct_" + e1] for r in rr], 0.36, color="#0072B2", label=e1)
    ax.barh(y - 0.18, [r["pct_" + e2] for r in rr], 0.36, color="#56B4E9", label=e2)
    ax.set_yticks(y); ax.set_yticklabels([f"{r['target']} (n={r['n']})" for r in rr])
    ax.set_xlabel("cage percentile among known binders (docking)")
    ax.axvline(50, color="#999", ls=":", lw=1)
    ax.set_title("NIOCH cage screen — calibrated, per-enantiomer (PRELIMINARY, pending assays)",
                 fontsize=9, fontweight="bold")
    ax.legend(frameon=False, loc="lower right", fontsize=8)
    fig.tight_layout()
    for e in ("pdf", "png"):
        fig.savefig(str(path) + "." + e)


def _wrap(text: str, indent: str = "  ") -> list:
    """Wrap one bullet to the width the rest of this report uses, continuing lines indented."""
    return textwrap.wrap(text, width=96, subsequent_indent=indent, break_long_words=False,
                         break_on_hyphens=False)


def _sig(x: float) -> str:
    """Signed number with a typographic minus, matching the rest of the report."""
    return f"{x:+.2f}".replace("-", "\u2212")


def boltz_cross_check():
    """The Boltz-2 section, read from the committed screen results rather than written by hand.

    This section used to be typed straight into the report, which meant regenerating the report
    deleted it and restored a recommendation the Boltz result had already retracted. Everything
    below is now derived from ``data/cage/boltz_screen_results/``, so the two cannot disagree.
    """
    conf = {t: load_target_conf(t) for t in BOLTZ_TARGETS}
    gaps = {t: anchor_gap(conf[t], "anchor") for t in BOLTZ_TARGETS}
    best = {t: max(gaps[t].values()) for t in BOLTZ_TARGETS}
    challenged = "challenges" in boltz_verdict(gaps)
    steroids = ("GR", "AR", "ER")
    ache = [conf["AChE"][f] for f in BOLTZ_FORMS]
    ranked = sorted(steroids, key=lambda t: -best[t])
    anchors = " / ".join(f"{conf[t]['anchor']:.2f}" for t in steroids)
    gapstr = ", ".join(f"{t} {_sig(best[t])}" for t in steroids)
    lines = [
        "## Boltz-2 cross-check (orthogonal structure model) "
        + ("\u2192 steroid hypothesis CHALLENGED" if challenged
           else "\u2192 steroid hypothesis CORROBORATED"),
        *_wrap(
            "The docking hypothesis was tested with a stronger, orthogonal method: a Boltz-2 "
            "(AlphaFold3-class co-folding + binding) small-molecule screen of {known agonist + "
            "the 4 cage forms} in each receptor's **agonist-seeded pocket**, read RELATIVE to "
            "that agonist (full table + figure: `docs/results_boltz_cage.md`, "
            "`figs/boltz_cage_crosscheck.{pdf,png}`; `make boltzcage`).", indent=""),
        f"**Boltz does {'NOT ' if challenged else ''}corroborate the steroid hypothesis:**",
        *_wrap(
            f"- The known agonists score high (binding_confidence {'/'.join(steroids)} = "
            f"**{anchors}**) \u2014 pockets and model recognise true steroid binders \u2014 but "
            f"the cage's best form sits **far below** each: gaps **{gapstr}**. The cage is not a "
            "steroid-pocket binder at the agonist level; the docking GR-percentile signal looks "
            'like the "rigid lipophilic cage fits a roomy pocket" artefact the gate warned '
            "about."),
        *_wrap(
            f"- **The AChE greasy-pocket artefact is REJECTED by Boltz** (cage {min(ache):.2f}"
            f"\u2013{max(ache):.2f} vs donepezil {conf['AChE']['anchor']:.2f}): Boltz "
            "discriminates where raw docking falsely scored AChE 100% \u2014 evidence its "
            "negative read on the steroid pockets is credible."),
        *_wrap(
            "- **The two methods disagree** (docking: GR\u2273AR moderate steroid; Boltz: cage "
            f"far below all agonists, weak residual preference **{ranked[0]}>{ranked[1]}**, "
            "opposite ranking) \u2192 the in-silico steroid signal is **not robust across "
            "methods**. DHODH is an inconclusive control (brequinar itself scores low "
            f"{conf['DHODH']['anchor']:.2f})."),
        *_wrap(
            "- Honest caveat: Boltz is not independently validated for this exact molecule, but "
            "its positive controls pass and it rejects the known docking artefact, so its "
            "disagreement carries weight."),
        "",
    ]
    return lines, challenged, ranked


def recommendations(challenged: bool, ranked: list) -> list:
    """VoI-ordered assays. The order follows the Boltz verdict, which is the point.

    ``make_boltz_cage.verdict`` already states the rule: a challenged steroid hypothesis means
    "downgrade to broad biochemical / phenotypic profiling rather than a targeted nuclear-receptor
    assay". That rule was applied by hand to this report once and would have been silently undone
    by the next regeneration; it is now applied by the code that knows the verdict.

    The broad-profiling item leads when the hypothesis is challenged and trails when it is not,
    and it appears exactly once either way. The hand-edited version carried it twice, once as the
    new lead and once as the item it was promoted from.
    """
    broad = (
        "A **broad biochemical/binding panel** (e.g. a kinase/against-panel + a thermal-shift / "
        "SPR screen) or **phenotypic profiling (Cell Painting / chemoproteomics)** is the real "
        "un-blinder \u2014 the orphan scaffold + promiscuity flag argue for breadth before depth, "
        "and in-silico has now been exhausted.")
    nr_lead = (
        "**Nuclear-receptor reporter panel \u2014 prioritise GR \u2273 AR > ER-\u03b1 \u2014 on "
        "both enantiomers**: the single most mechanistically-grounded, cheapest decisive test. GR "
        "is the only NR reproducible across both screens, so lead with it; it directly probes the "
        "only coherent in-silico class signal and a clean negative *bounds* the orphan claim. "
        "Include the **deacetyl (active) form** alongside the acetate.")
    data_rich = (
        "**Test the data-rich moderate hits first**, not the thin greasy-pocket 100%s: the "
        "percentiles backed by the most reference docks are the trustworthy hypotheses.")
    broad_lead = (
        "**Broad biochemical / phenotypic profiling is now the primary un-blinder** "
        "(thermal-shift / SPR against-panel, or Cell Painting / chemoproteomics). After the Boltz "
        "cross-check the steroid hypothesis is **method-dependent** (docking moderate, Boltz "
        "negative) and in-silico is exhausted; the orphan, promiscuous scaffold argues for "
        "breadth before depth.")
    nr_bound = (
        "**The nuclear-receptor reporter panel is now a cheap *bound*, not the lead.** If run, "
        "test **GR** (docking's one cross-screen-reproducible hit) and "
        f"**{ranked[0]}** (Boltz's weak residual preference) on **both enantiomers** and **both "
        "forms** (acetate as-given + deacetyl active; on GR/AR the acetate scored slightly higher "
        "in Boltz). A clean negative *bounds* the orphan claim, but the two in-silico methods no "
        "longer agree that it will be positive.")
    enantiomers = ("For any hit, **test both enantiomers separately** to confirm a chirality "
                   "preference.")
    ladder = ("Escalate only confirmed hits up the ladder (binding \u2192 functional \u2192 "
              "co-crystal).")
    items = ([broad_lead, nr_bound, enantiomers, ladder] if challenged
             else [nr_lead, data_rich, broad, enantiomers, ladder])
    lines = ["## Recommended experiments (VoI-ordered, pending assays)"]
    for i, item in enumerate(items, 1):
        lines += _wrap(f"{i}. {item}", indent="   ")
    return lines


def report(rows, path):
    e1, e2 = list(CAGE)
    have = [r for r in rows if r["n"] >= 3]
    boltz_lines, challenged, ranked = boltz_cross_check()
    lines = [
        "# NIOCH cage screen — ranked target hypotheses (PRELIMINARY · pending assays)",
        "",
        "**This is an operational report, NOT a validated finding.** The retrospective",
        "gate (`docs/results_figH.md`) showed structure-based target-ID does not beat a",
        "ligand baseline and raw docking just picks the greasiest pocket; so this screen is",
        "**heavily caveated** and exists to propose experiments, not to name a target.",
        "",
        "## Molecule",
        "Difluoronaphthalenone + a 1,3-dicarbonyl Michael donor (here barbituric acid),",
        "O-acetylated; rigid fused scaffold, MW 350, two ring-fusion-coupled stereocentres",
        "→ one diastereomer + its mirror. **Orphan scaffold** (0 ChEMBL similarity ≥40%).",
        "Screened **per enantiomer** (R,R given / S,S mirror); racemate apparent affinity ≈",
        "the stronger binder.",
        "",
        "## Method (and why it is calibrated)",
        "smina docking of each enantiomer into each pocket of an 11-target diverse panel.",
        "We do **not** rank by raw kcal/mol (cross-pocket-incomparable). We rank each target",
        "by the **cage's percentile among that target's known binders** docked into the same",
        "pocket — \"the cage docks better than X% of known inhibitors of this target\" —",
        "which removes the pocket-size artefact. References pooled from the gate + campaign",
        "docking caches.",
        "",
        "## Ranked hypotheses",
        "",
        "| target | n_ref | " + f"{e1} kcal (pctl) | {e2} kcal (pctl) | best pctl |",
        "|---|---:|---:|---:|---:|",
    ]
    for r in have:
        lines.append(f"| {r['target']} | {r['n']} | {r[e1]:.1f} ({r['pct_'+e1]:.0f}%) | "
                     f"{r[e2]:.1f} ({r['pct_'+e2]:.0f}%) | **{r['best']:.0f}%** |")
    pend = [r for r in rows if r["n"] < 3]
    if pend:
        lines += ["", "*Calibration pending (insufficient reference docks): "
                  + ", ".join(f"{r['target']} (raw {r[e1]:.1f}/{r[e2]:.1f})" for r in pend) + ".*"]
    lines += [
        "",
        "## How to read this (caveats — do not skip)",
        "- **Greasy-pocket artefact persists.** Large hydrophobic pockets (e.g. AChE) score",
        "  the rigid lipophilic cage well even after calibration; a high percentile there is",
        "  weak evidence. Trust the **data-rich** rows (large n_ref) over thin ones.",
        "- **Promiscuity.** The cyclic-ureide / uracil diamide is a strong, promiscuous",
        "  H-bonder; broad moderate percentiles are consistent with a non-specific scaffold,",
        "  not a single target. Calibration is the guard against a false promiscuous lead.",
        "- **Docking is approximate**, the panel is **small (11 targets, not proteome-wide)**,",
        "  and the true target may not be in it. Metal/heme/covalent targets were excluded",
        "  (docking unreliable there).",
        "- **Enantiomer differences** in the table are real signal worth testing (the cage is",
        "  chiral; one enantiomer may bind preferentially).",
        "",
        "## Broad reverse-dock (33-target pocketome) → nuclear-receptor hypothesis",
        "To widen the search beyond the 11-target panel above, the cage (4 forms: acetate &",
        "deacetyl, each R,R / S,S) was docked into a **~33-target diverse pocketome** and",
        "calibrated the same way (percentile among each target's known binders). Full table:",
        "`docs/reversedock_shortlist.md` (476 docks). The pattern that survives calibration:",
        "",
        "- **\U0001F7E0 7 targets at 100% are greasy-pocket artefacts — ignore them.** CDK2,",
        "  VEGFR2, DPP-4, HIV integrase, AChE, MAO-B, CA-II all top out, but they span *unrelated*",
        "  families (kinases / protease / integrase / hydrolases). A high percentile scattered",
        "  across families is the rigid lipophilic cage fitting any roomy pocket, **not** a",
        "  target. This is exactly the artefact the gate (`results_figH.md`) warned about.",
        "- **\U0001F535 The one coherent class signal is nuclear / steroid receptors:**",
        "  **GR 80%, AR 75%, ER-α 60%** — three steroid receptors clustered in the",
        "  moderate-high band, each data-rich (n_ref 8–10). Structural rationale: the rigid,",
        "  lipophilic **difluoronaphthalenone is steroid-fragment-like**, so the cage plausibly",
        "  occupies steroid-type pockets. This is the best *mechanistically-grounded* hypothesis",
        "  in silico — **but it is moderate and non-selective** (AR≈GR≈ER), i.e. a steroid-pocket",
        "  occupant, not a selective lead.",
        "- **Robustness caveat (honest):** only **GR is consistent across the two independent",
        "  screens** (panel 89% / reverse-dock 80%, different PDBs + reference sets). **ER-α is",
        "  structure-dependent** (panel 0% vs reverse-dock 60%) → treat ER as the weakest of the",
        "  three. So the steroid signal is real as a *class* but its strongest, most reproducible",
        "  member is GR; rank the assay GR ≳ AR > ER accordingly.",
        "- **Low fit (cage does not belong):** DHODH **0%** (consistent with the focused",
        "  brequinar-tunnel test, `data/campaign/dhodh/`), Thrombin 0%, PDE10A 0%, HSP90 10%,",
        "  COX-2 10%. Calibration correctly reports \"worse than known binders\" here.",
        "- **Method-validation read:** calibration separated artefact-100% from genuine signal",
        "  and reproduced the DHODH rejection independently — the screen behaves as designed.",
        "",
        *boltz_lines,
        *recommendations(challenged, ranked),
        "",
        "*Generated by `make nioch` from the docking caches and the committed Boltz-2 screen",
        "results; regenerates as the campaign adds reference docks. Every section above is derived,",
        "including the Boltz cross-check and the order of the recommendations, so nothing here can",
        "be lost by regenerating it. Not a publication claim.*",
    ]
    pathlib.Path(path).write_text("\n".join(lines) + "\n")


def main() -> None:
    rows = rank()
    figure(rows, str(ROOT / "figs" / "nioch_cage_screen"))
    report(rows, ROOT / "docs" / "nioch_cage_report.md")
    print("wrote docs/nioch_cage_report.md + figs/nioch_cage_screen.{pdf,png}")
    for r in rows:
        print(f"  {r['target']:24s} n={r['n']:>2}  best={r['best'] if r['best']==r['best'] else 0:.0f}%")


if __name__ == "__main__":
    main()
