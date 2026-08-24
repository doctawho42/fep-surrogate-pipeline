"""Fig Ground: the visible fraction on curated labels, over fourteen benchmark systems.

The article's headline -- that a network's error against experiment lives almost entirely in the
direction cycle closure cannot see -- is measured on four systems, all four selected by the closure
test itself, and on experimental labels taken from a stereo-blind, target-wide ChEMBL search. This
script runs the same measurement on the OpenFF benchmark's own CURATED per-edge experimental
``ddG``, for every system whose name the OpenFE replicate benchmark shares (fourteen of fifteen;
``pde2`` is not in the replicate benchmark), and pools the result over the systems the closure test
flags and, separately, over those it does not.

The join rule, the recovery threshold and the pooling rule were fixed before any of the eleven
systems beyond the Supporting Information's three was measured; they are stated in
``docs/results_figGround.md`` and implemented in ``src/bar/curated.py``.

Run: PYTHONPATH=src python figs/make_figGround.py   (or `make figGround`)
"""
from __future__ import annotations

import csv
import math
import pathlib
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
from paperstyle import (  # noqa: E402
    FULL,
    INK,
    MUTED,
    OURS,
    REF,
    check_min_type,
    finish,
    legend,
    panel,
    tint,
    use_paper_style,
)

from bar.curated import (  # noqa: E402
    Reading,
    curated_edges,
    match_system,
    name_key,
    network_edges,
    pool,
    read_ligand_cache,
    read_system,
    resolve_ligands,
    system_groups,
    visible_fraction,
    write_ligand_cache,
)
from bar.qc import benjamini_hochberg, chi2_sf, gls_network  # noqa: E402

DOC = ROOT / "docs" / "results_figGround.md"
TABLE = ROOT / "docs" / "tab_ground.tex"
COMBINED = ROOT / "data" / "openfe_replicates" / "combined_pymbar4_edge_data.csv"
CACHE = ROOT / "data" / "curated" / "openfe_ligand_smiles.csv"
CACHE_NEUTRAL = ROOT / "data" / "curated" / "openfe_ligand_smiles_neutral.csv"

#: The fifteen targets carrying curated per-edge labels, less the one the replicate benchmark does
#: not contain. Fixed by the two data sets, not chosen.
SYSTEMS = ["cdk2", "cdk8", "cmet", "eg5", "hif2a", "mcl1", "p38", "pfkfb3", "ptp1b", "shp2",
           "syk", "thrombin", "tnks2", "tyk2"]

#: Declared before any of the eleven new systems was measured: a system whose curated edges recover
#: below this rate into the network is reported but excluded from the thresholded pooled figure.
RECOVERY_MIN = 0.50

#: The published four-system ChEMBL-grounded figures, quoted from docs/results_figHodge.md for
#: comparison only; nothing here recomputes them.
PUBLISHED = {"f": 0.0060, "chance": 0.325, "E": 167, "dof": 51, "systems": 4}


# --------------------------------------------------------------------------------------
# The article's flag, recomputed rather than copied.
# --------------------------------------------------------------------------------------
def _replicate_edge(row: dict, k: int):
    def value(key: str) -> float:
        try:
            return float(row[key])
        except (TypeError, ValueError, KeyError):
            return math.nan
    complex_dg = value(f"complex_repeat_{k}_DG (kcal/mol)")
    complex_se = value(f"complex_repeat_{k}_dDG (kcal/mol)")
    solvent_dg = value(f"solvent_repeat_{k}_DG (kcal/mol)")
    solvent_se = value(f"solvent_repeat_{k}_dDG (kcal/mol)")
    if any(math.isnan(v) for v in (complex_dg, complex_se, solvent_dg, solvent_se)):
        return None
    return complex_dg - solvent_dg, math.sqrt(complex_se ** 2 + solvent_se ** 2)


def ever_flagged() -> set[str]:
    """Systems the Benjamini--Hochberg closure test flags on at least one of the three replicates.

    The article's own rule, over the whole 49-system benchmark, so a system's flag here means what
    it means there. Independently implemented from ``make_figHodge.py``'s copy, which is a check.
    """
    by_system: dict[str, list[dict]] = {}
    with COMBINED.open(newline="") as handle:
        for row in csv.DictReader(handle):
            by_system.setdefault(row["system name"], []).append(row)
    flagged: set[str] = set()
    for replicate in (0, 1, 2):
        names, pvals = [], []
        for name, rows in sorted(by_system.items()):
            edges = [(r["ligand_A"], r["ligand_B"], *_replicate_edge(r, replicate))
                     for r in rows if _replicate_edge(r, replicate)]
            if len(edges) < 3:
                continue
            fit = gls_network(edges)
            if fit.dof < 1:
                continue
            names.append(name)
            pvals.append(chi2_sf(fit.chi2, fit.dof))
        hits = benjamini_hochberg(pvals, alpha=0.05)
        flagged |= {n for n, on in zip(names, hits, strict=True) if on}
    return flagged


# --------------------------------------------------------------------------------------
# The measurement.
# --------------------------------------------------------------------------------------
def _ligand_structures(neutralize: bool, memo: dict[str, str | None]) -> dict:
    """The per-release ligand structures for every system, from the committed cache or the network.

    One HTTP fetch per SDF at most, shared between the exact and the neutralized readings.
    """
    path = CACHE_NEUTRAL if neutralize else CACHE
    cached = read_ligand_cache(path)
    if all(s in cached and cached[s] for s in SYSTEMS):
        return cached

    def fetch(url: str) -> str | None:
        if url not in memo:
            from bar.curated import _fetch
            memo[url] = _fetch(url)
        return memo[url]

    resolved = {s: resolve_ligands(s, system_groups(s, COMBINED), fetch, neutralize)
                for s in SYSTEMS}
    write_ligand_cache(resolved, path)
    return resolved


def readings(neutralize: bool, flagged: set[str], memo: dict[str, str | None],
             by_structure: bool = True) -> list[Reading]:
    """One reading per system, in the fixed system order.

    ``by_structure`` is the node key of ``match_system``; the article quotes the structure-keyed
    reading and reports the name-keyed one beside it.
    """
    structures = _ligand_structures(neutralize, memo)
    out = []
    for system in SYSTEMS:
        match = match_system(
            network_edges(system, COMBINED),
            structures.get(system, {}),
            curated_edges(system, ROOT / "data" / "fep_edges", neutralize=neutralize),
            by_structure=by_structure,
        )
        out.append(read_system(system, match, system in flagged))
    return out


def pools(rows: list[Reading]) -> dict[str, dict[str, float]]:
    """Every pooled figure the record reports, from one set of readings."""
    return {
        "all": pool(rows),
        "thresholded": pool([r for r in rows if r.recovery >= RECOVERY_MIN]),
        "flagged": pool([r for r in rows if r.flagged]),
        "unflagged": pool([r for r in rows if not r.flagged]),
    }


def chembl_labels(system: str) -> dict[str, float]:
    """The committed ChEMBL-grounded per-ligand experimental dG the article's headline uses."""
    path = ROOT / "data" / "openfe_replicates" / f"affinity_{system}.csv"
    if not path.exists():
        return {}
    with path.open(newline="") as handle:
        return {r["ligand"]: float(r["exp_dg"]) for r in csv.DictReader(handle)}


def same_edge_comparison(rows: list[Reading], structures: dict) -> list[dict]:
    """The two label sources on THE SAME edges, for the systems that carry both.

    The pooled curated figure and the article's pooled ChEMBL figure are measured on different
    edge sets, so their ratio confounds the label source with the sub-network. This holds the edge
    set fixed: for each system with a committed ChEMBL affinity table, restrict to the matched
    sub-network's edges whose two ligands both carry a ChEMBL value, and measure the visible
    fraction of both error fields there.
    """
    out = []
    for system in SYSTEMS:
        exp = chembl_labels(system)
        if not exp:
            continue
        curated = curated_edges(system, ROOT / "data" / "fep_edges")
        by_pair = {}
        for first, second, ddg in curated:
            by_pair.setdefault(frozenset((first, second)), (first, second, ddg))
        edges, eps_cur, eps_chembl = [], [], []
        for edge in network_edges(system, COMBINED):
            head = structures.get(system, {}).get(edge.group, {}).get(name_key(edge.a))
            tail = structures.get(system, {}).get(edge.group, {}).get(name_key(edge.b))
            if head is None or tail is None or head == tail:
                continue
            hit = by_pair.get(frozenset((head, tail)))
            if hit is None or edge.a not in exp or edge.b not in exp:
                continue
            curated_head, _tail, ddg = hit
            edges.append((edge.a, edge.b, edge.ddg, edge.se))
            eps_cur.append(edge.ddg - (ddg if curated_head == head else -ddg))
            eps_chembl.append(edge.ddg - (exp[edge.b] - exp[edge.a]))
        if len(edges) < 2:
            continue
        f_cur, chance, iso, dof = visible_fraction(edges, np.array(eps_cur))
        f_chembl, _c, _i, _d = visible_fraction(edges, np.array(eps_chembl))
        out.append({"system": system, "E": len(edges), "dof": dof, "chance": chance,
                    "chance_iso": iso, "curated": f_cur, "chembl": f_chembl,
                    "median_cur": float(np.median(np.abs(eps_cur))),
                    "median_chembl": float(np.median(np.abs(eps_chembl)))})
    return out


# --------------------------------------------------------------------------------------
# The figure.
# --------------------------------------------------------------------------------------
def draw(rows: list[Reading], pooled: dict[str, dict[str, float]],
         same_edge: list[dict]) -> None:
    use_paper_style()
    fig = plt.figure(figsize=(FULL, 5.5), layout="constrained")
    grid = fig.add_gridspec(2, 3, height_ratios=(1.2, 1.0))
    ax_a = fig.add_subplot(grid[0, :])
    ax_b = fig.add_subplot(grid[1, 0])
    ax_c = fig.add_subplot(grid[1, 1])
    ax_d = fig.add_subplot(grid[1, 2])

    # A: per-system visible fraction against its two chance levels.
    live = sorted([r for r in rows if r.measurable], key=lambda r: r.visible)
    y = np.arange(len(live))
    for i, r in enumerate(live):
        ax_a.plot([max(r.loo_visible[0], 1e-6), r.loo_visible[1]], [i, i],
                  color=tint(OURS, 0.55), linewidth=3.0, solid_capstyle="butt", zorder=2)
    ax_a.scatter([r.visible for r in live], y, s=26, color=OURS, zorder=3,
                 label="$f$, curated labels")
    ax_a.scatter([r.chance for r in live], y, s=30, facecolor="none", edgecolor=REF,
                 linewidth=1.0, zorder=3, label=r"chance $\mathrm{dof}/E$")
    ax_a.scatter([r.chance_iso for r in live], y, s=26, marker="|", color=REF, zorder=3,
                 label=r"chance $\mathrm{tr}(\Pi W)/\mathrm{tr}(W)$")
    ax_a.set_yticks(y)
    ax_a.set_yticklabels([f"{r.system}{'*' if r.flagged else ''}" for r in live])
    ax_a.set_xscale("log")
    ax_a.set_xlim(1e-5, 1.0)
    ax_a.set_xlabel("share of the squared standardized error that closure can see")
    ax_a.set_ylim(-0.7, len(live) - 0.3)
    panel(ax_a, "A", "every system with a cycle sits far below its own chance level",
          subtitle="bar: single-edge-deletion range of $f$; * = flagged by the closure test")
    legend(ax_a, loc="upper center", bbox_to_anchor=(0.5, -0.20), ncol=3)

    # B: what survives the join, per system.
    order = sorted(rows, key=lambda r: -r.n_curated)
    pos = np.arange(len(order))
    ax_b.barh(pos, [r.n_curated for r in order], color=tint(MUTED, 0.55), height=0.72,
              label="curated")
    ax_b.barh(pos, [r.n_recovered for r in order], color=tint(OURS, 0.55), height=0.72,
              label="recovered")
    ax_b.barh(pos, [r.n_edges for r in order], color=OURS, height=0.42,
              label="matched")
    for i, r in enumerate(order):
        if not r.measurable:
            ax_b.annotate("no cycle" if r.n_edges else "no match", xy=(r.n_curated + 1.5, i),
                          va="center", ha="left", fontsize=7.5, color=REF)
    ax_b.set_yticks(pos)
    ax_b.set_yticklabels([r.system for r in order])
    ax_b.invert_yaxis()
    ax_b.set_xlim(0, 96)
    ax_b.set_xlabel("edges")
    panel(ax_b, "B", "what the join keeps")
    legend(ax_b, loc="upper center", bbox_to_anchor=(0.5, -0.26), ncol=3, handlelength=1.0,
           columnspacing=1.0)

    # C: the pooled figures, against the published four-system one.
    labels = ["ChEMBL, 4 sys.", "curated, all 8", "curated, flagged", "curated, not fl."]
    values = [PUBLISHED["f"], pooled["all"]["f"], pooled["flagged"]["f"],
              pooled["unflagged"]["f"]]
    chances = [PUBLISHED["chance"], pooled["all"]["chance"], pooled["flagged"]["chance"],
               pooled["unflagged"]["chance"]]
    isos = [math.nan, pooled["all"]["chance_iso"], pooled["flagged"]["chance_iso"],
            pooled["unflagged"]["chance_iso"]]
    x = np.arange(len(labels))
    ax_c.bar(x, values, width=0.56, color=[tint(OURS, 0.45)] + [OURS] * 3,
             label="pooled $f$")
    ax_c.scatter(x, chances, s=30, facecolor="none", edgecolor=REF, linewidth=1.0, zorder=3,
                 label=r"chance $\mathrm{dof}/E$")
    ax_c.scatter(x, isos, s=26, marker="|", color=REF, zorder=3, label="chance iso.")
    ax_c.set_yscale("log")
    ax_c.set_ylim(1e-4, 30.0)
    ax_c.set_yticks([1e-4, 1e-3, 1e-2, 1e-1, 1e0])
    ax_c.set_xticks(x)
    ax_c.set_xticklabels(labels, rotation=30, ha="right")
    ax_c.set_ylabel("pooled fraction")
    panel(ax_c, "C", "pooled, split by the flag")
    legend(ax_c, loc="upper left")

    # D: the two label sources on the same edges, which holds the sub-network fixed.
    live_same = [s for s in same_edge if s["dof"] >= 1]
    xs = np.arange(len(live_same))
    ax_d.bar(xs - 0.19, [s["curated"] for s in live_same], width=0.34, color=OURS,
             label="curated")
    ax_d.bar(xs + 0.19, [s["chembl"] for s in live_same], width=0.34, color=tint(OURS, 0.55),
             label="ChEMBL")
    ax_d.scatter(xs, [s["chance_iso"] for s in live_same], s=26, marker="|", color=REF, zorder=3,
                 label="chance iso.")
    ax_d.set_yscale("log")
    ax_d.set_ylim(1e-4, 30.0)
    ax_d.set_yticks([1e-4, 1e-3, 1e-2, 1e-1, 1e0])
    ax_d.set_xticks(xs)
    ax_d.set_xticklabels([f"{s['system']}\nE = {s['E']}" for s in live_same])
    ax_d.set_ylabel("visible fraction")
    panel(ax_d, "D", "same edges")
    legend(ax_d, loc="upper left")

    for spine in (ax_a, ax_b, ax_c, ax_d):
        spine.tick_params(color=INK)
    offenders = check_min_type(fig)
    if offenders:
        raise AssertionError(f"type below the floor: {offenders}")
    finish(fig, "figGround_curated_grounding")


# --------------------------------------------------------------------------------------
# The records.
# --------------------------------------------------------------------------------------
def _fmt(x: float, digits: int = 4) -> str:
    return "n/a" if x is None or (isinstance(x, float) and math.isnan(x)) else f"{x:.{digits}f}"


def write_table(rows: list[Reading], pooled: dict[str, dict[str, float]]) -> None:
    live = [r for r in rows if r.measurable]
    lines = [
        r"\begin{table}[tp]\centering",
        r"\caption{\textbf{The visible fraction on curated labels, over the benchmark.} For each "
        r"target whose curated per-edge experimental $\ddG$ the replicate benchmark can be joined "
        r"to: curated edges $E_{\mathrm{cur}}$, the share of them whose two ligands both recover "
        r"into the network ($\rho$), and the matched sub-network's ligands $N$, edges $E$, "
        r"components $c$ and independent cycles $\mathrm{dof}$. Then the visible fraction $f$ of "
        r"the error against those labels, its $\mathrm{dof}/E$ chance level, the chance level "
        r"$\mathrm{tr}(\Pi W)/\mathrm{tr}(W)$ an error isotropic in kcal/mol would produce under "
        r"the same projector, the range of $f$ over single-edge deletions, the worst margin "
        r"$\mathrm{chance}/f$ any deletion leaves under each chance level, and the median "
        r"$|\epsilon|$ against experiment. Both the fraction and its chance level are recomputed "
        r"on every deletion, since deleting an edge moves the projector as well. $f$ and its "
        r"deletion range are given in units of $10^{-3}$. A dagger marks the systems the closure "
        r"test flags. The systems whose matched sub-network carries no cycle, and so no "
        r"measurement, are named below the rule.}",
        r"\label{tab:ground}",
        r"\footnotesize",
        r"\setlength{\tabcolsep}{3.4pt}",
        r"\begin{tabular}{@{}lrrrrrrrrrlrrr@{}}",
        r"\toprule",
        r"system & $E_{\mathrm{cur}}$ & $\rho$ & $N$ & $E$ & $c$ & $\mathrm{dof}$ & "
        r"$f/10^{-3}$ & $\mathrm{dof}/E$ & iso & $f$ range$/10^{-3}$ & "
        r"\multicolumn{2}{c}{worst margin} & $|\epsilon|$ \\",
        r" &  &  &  &  &  &  &  &  &  &  & dof/$E$ & iso &  \\",
        r"\midrule",
    ]
    for r in sorted(live, key=lambda r: r.system):
        dagger = r"$^\dagger$" if r.flagged else ""
        lines.append(
            rf"\texttt{{{r.system}}}{dagger} & ${r.n_curated}$ & ${r.recovery:.2f}$ & "
            rf"${r.n_nodes}$ & ${r.n_edges}$ & ${r.components}$ & ${r.dof}$ & "
            rf"${1e3 * r.visible:.2f}$ & ${r.chance:.3f}$ & ${r.chance_iso:.3f}$ & "
            rf"$[{1e3 * r.loo_visible[0]:.2f}, {1e3 * r.loo_visible[1]:.2f}]$ & "
            rf"${r.loo_worst_margin:.0f}\times$ & ${r.loo_worst_margin_iso:.1f}\times$ & "
            rf"${r.median_abs_eps:.2f}$ \\")
    lines.append(r"\midrule")
    for key, label in (("all", "pooled, all"), ("flagged", "pooled, flagged"),
                       ("unflagged", "pooled, not flagged")):
        p = pooled[key]
        lines.append(
            rf"\multicolumn{{4}}{{@{{}}l}}{{{label} ({int(p['systems'])} systems)}} & "
            rf"${int(p['E'])}$ &  & ${int(p['dof'])}$ & ${1e3 * p['f']:.2f}$ & "
            rf"${p['chance']:.3f}$ & ${p['chance_iso']:.3f}$ &  &  &  &  \\")
    dead = [r for r in rows if not r.measurable]
    no_cycle = ", ".join(rf"\texttt{{{r.system}}} ($E={r.n_edges}$, $\mathrm{{dof}}=0$)"
                         for r in dead if r.n_edges)
    no_match = ", ".join(rf"\texttt{{{r.system}}} ($\rho={r.recovery:.2f}$)"
                         for r in dead if not r.n_edges)
    lines += [
        r"\bottomrule",
        r"\end{tabular}",
        r"\vspace{2pt}",
        r"{\footnotesize No measurement, matched sub-network acyclic: " + no_cycle + r". "
        r"No measurement, no curated edge joined: " + no_match + r".\par}",
        r"\end{table}",
    ]
    TABLE.write_text("\n".join(lines) + "\n")
    print(f"wrote {TABLE}")


def write_doc(grid: dict[tuple[bool, bool], tuple[list[Reading], dict[str, dict[str, float]]]],
              same_edge: list[dict]) -> None:
    """The record, from every (neutralize, by_structure) reading.

    Keyed ``(neutralize, by_structure)``. ``(False, True)`` is the measurement; the other three
    are the two declared alternatives and their combination, which is where the two keys
    interact.
    """
    rows, pooled = grid[(False, True)]
    name_rows, name_pooled = grid[(False, False)]
    neutral_rows, neutral_pooled = grid[(True, True)]
    neutral_name_rows, neutral_name_pooled = grid[(True, False)]
    live = [r for r in rows if r.measurable]
    dead = [r for r in rows if not r.measurable]
    all_pool, flag_pool, unflag_pool = pooled["all"], pooled["flagged"], pooled["unflagged"]
    lines = [
        "# Fig Ground: the visible fraction on curated labels, over fourteen systems",
        "",
        "Generated by `figs/make_figGround.py` (`make figGround`); the join and the measurement "
        "live in",
        "`src/bar/curated.py`.",
        "",
        "## The headline",
        "",
        f"Pooled over the {int(all_pool['systems'])} systems whose matched sub-network carries a "
        f"cycle -- {int(all_pool['E'])} edges and {int(all_pool['dof'])} residual degrees of "
        f"freedom -- the visible fraction on curated labels is "
        f"**{all_pool['f']:.5f}**, against a `dof/E` chance level of {all_pool['chance']:.3f} and "
        f"an isotropic-in-kcal/mol chance level of {all_pool['chance_iso']:.4f}. That is a factor "
        f"of {all_pool['chance'] / all_pool['f']:.0f} below the first and "
        f"{all_pool['chance_iso'] / all_pool['f']:.0f} below the second.",
        "",
        f"The article reports {PUBLISHED['f']:.4f} over {PUBLISHED['systems']} systems, "
        f"{PUBLISHED['E']} edges and {PUBLISHED['dof']} residual degrees of freedom. **The "
        f"expanded measurement moves the pooled number DOWN**, by a factor of "
        f"{PUBLISHED['f'] / all_pool['f']:.1f}. Three things change at once: eleven new systems "
        "enter, the matched sub-networks are smaller and sparser in cycles than the full ones, "
        "and the label source changes. Only the third is isolated below, and it pushes the OTHER "
        "way -- on edges held fixed the curated labels give a HIGHER visible fraction than the "
        "ChEMBL join, by a factor of three to four (see 'Same edges, two label sources'). The "
        "pooled figures below are therefore two grounded measurements of different things, and "
        "neither corrects the other; what carries across is that every reading, on every system "
        "that carries a cycle, sits one to three orders of magnitude below its own chance level.",
        "",
        "Split by the flag, which is the point of the exercise:",
        "",
        "| pool | systems | E | dof | f | chance dof/E | chance isotropic |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for key, label in (("all", "all matched"), ("thresholded", f"recovery >= {RECOVERY_MIN:.2f}"),
                       ("flagged", "flagged by closure"), ("unflagged", "not flagged")):
        p = pooled[key]
        lines.append(
            f"| {label} | {int(p['systems'])} | {int(p['E'])} | {int(p['dof'])} | "
            f"{p['f']:.5f} | {p['chance']:.3f} | {p['chance_iso']:.4f} |")
    lines += [
        "",
        f"The systems the closure test does NOT flag show the effect at least as strongly as the "
        f"ones it does: {unflag_pool['f']:.5f} against {flag_pool['f']:.5f}, or "
        f"{unflag_pool['chance_iso'] / unflag_pool['f']:.0f}x below the isotropic chance level "
        f"against {flag_pool['chance_iso'] / flag_pool['f']:.0f}x. The measurement is therefore a "
        "statement about the benchmark, not a property of the networks the test selected.",
        "",
        "## The protocol, frozen before the eleven new systems were measured",
        "",
        "1. **Systems.** The fifteen targets carrying curated per-edge experimental `ddG` in",
        "   `data/fep_edges/`, less `pde2`, which the OpenFE replicate benchmark does not "
        "contain: fourteen.",
        "2. **Structures.** A network ligand name resolves through the benchmark's own prepared",
        "   input SDF for the release that edge belongs to",
        "   (`prepared_structures/<system group>/<system name>/ligands.sdf`), so two releases that",
        "   reuse a name for different molecules cannot be mixed. Names are compared",
        "   separator-insensitively (runs of whitespace, `-` and `_` collapse to one space);",
        "   nothing else about a name is rewritten and an unresolved name is never given a",
        "   structure.",
        "3. **Match.** RDKit isomeric canonical SMILES on both sides; a curated edge matches a",
        "   network edge when the unordered canonical-SMILES pairs are equal, and the curated",
        "   `ddG` is signed to the network edge's direction.",
        "4. **Recovery threshold.** `recovery = (curated edges whose two ligands both resolve into",
        "   the network) / (curated edges)`, the Supporting Information's convention, with the",
        f"   threshold set at {RECOVERY_MIN:.2f} before any of the eleven new systems was",
        "   measured. Both the thresholded and the unthresholded pooled figures are reported.",
        "5. **Pooling.** Ratios of summed squared norms, as the article pools: "
        "`f = sum ||Pi eps~||^2 / sum ||eps~||^2`, `chance = sum dof / sum E`, "
        "`chance_iso = sum tr(Pi W) / sum tr(W)`. A system whose matched sub-network carries no",
        "   cycle contributes to neither side.",
        "",
        "## Reproduction of the Supporting Information's three-system cross-check",
        "",
        "The Supporting Information reports recovering 57%, 98% and 91% of the curated edges for",
        "`cdk8`, `hif2a` and `p38`, a visible fraction of 0.0007 for `cdk8` against a chance level",
        "of 0.071, and 0.0006 for `hif2a` against 0.147. That check's code was never committed, so",
        "the rule was reconstructed from its reported numbers. It reproduces `hif2a` and `p38`",
        "exactly -- 0.980 and 0.911 recovery, `hif2a` f = 0.00063 against dof/E = 0.1471, `p38`",
        "matched sub-network with dof = 0, the 'retains no cycle' the text records.",
        "",
        "**It does not reproduce `cdk8`, and the difference is a fix rather than a discrepancy.**",
        "`cdk8` is the one system here whose two releases (`merck`, `miscellaneous_set`) reuse",
        "seven ligand names for different molecules. The reconstruction that matches the",
        "Supporting Information's 57% builds one name-to-structure table per system, in which the",
        "second release silently overwrites the first on those seven names; resolving each edge in",
        "its own release instead recovers 44 of 46 curated edges (95.7%) and a 29-edge matched",
        "sub-network rather than a 15-edge one. `cdk8`'s visible fraction rises from 0.0007 to",
        "0.0023 as a result, against chance levels of 0.172 and 0.077 -- the release-local rule is",
        "the less favourable of the two for the article, and it is the one used throughout here.",
        "",
        "## Per system",
        "",
        "| system | flag | curated | recovered | recovery | N | E | c | dof | f | chance dof/E |"
        " chance iso | f range (1-edge deletion) | worst margin dof/E | worst margin iso |"
        " median abs eps |",
        "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---:|---:|---:|",
    ]
    for r in rows:
        flag = "yes" if r.flagged else "no"
        rng = (f"[{r.loo_visible[0]:.5f}, {r.loo_visible[1]:.5f}]"
               if r.measurable else "-")
        lines.append(
            f"| `{r.system}` | {flag} | {r.n_curated} | {r.n_recovered} | {r.recovery:.3f} | "
            f"{r.n_nodes} | {r.n_edges} | {r.components} | {r.dof} | "
            f"{_fmt(r.visible, 5)} | {_fmt(r.chance, 3)} | {_fmt(r.chance_iso, 4)} | {rng} | "
            f"{_fmt(r.loo_worst_margin, 1)} | {_fmt(r.loo_worst_margin_iso, 1)} | "
            f"{_fmt(r.median_abs_eps, 2)} |")
    no_cycle = [r for r in dead if r.n_edges]
    no_match = [r for r in dead if not r.n_edges]
    lines += [
        "",
        "## What carries no measurement, and why",
        "",
        f"- **Acyclic matched sub-network (dof = 0), {len(no_cycle)} systems**: "
        + ", ".join(f"`{r.system}` (E = {r.n_edges}, N = {r.n_nodes}, c = {r.components})"
                    for r in no_cycle)
        + ". A sub-network with no independent cycle has nothing for a closure diagnostic to see "
          "and nothing for the projector to keep; its visible fraction is identically zero by "
          "construction, which is not evidence, so these contribute to neither side of the pool.",
        f"- **No matched edge at all, {len(no_match)} systems**: "
        + ", ".join(f"`{r.system}` (recovery {r.recovery:.3f})" for r in no_match) + ".",
    ]
    thrombin = next((r for r in rows if r.system == "thrombin"), None)
    pfkfb3 = next((r for r in rows if r.system == "pfkfb3"), None)
    if thrombin is not None and pfkfb3 is not None:
        lines += [
            f"  `thrombin` recovers {thrombin.n_recovered} of {thrombin.n_curated} curated edges "
            "even though every one of its network ligand names resolves to a structure: the "
            "benchmark's prepared thrombin ligands are protonated (`[NH3+]`) where its curated "
            "table is neutral, and the exact-SMILES rule refuses the join (the sensitivity "
            "analysis below lifts exactly this). `pfkfb3` fails for a different reason -- the "
            f"replicate benchmark carries only {len(network_edges('pfkfb3', COMBINED))} pfkfb3 "
            f"edges in total, against the curated table's {pfkfb3.n_curated}.",
        ]
    n_thresh_excluded = [r for r in rows if r.recovery < RECOVERY_MIN]
    lines += [
        "",
        f"The recovery threshold of {RECOVERY_MIN:.2f} excludes "
        + (", ".join(f"`{r.system}` ({r.recovery:.3f})" for r in n_thresh_excluded) or "nothing")
        + ". Both excluded systems already contribute no matched edge, so the thresholded and "
          "unthresholded pooled figures are identical to every digit reported above; the "
          "threshold changes nothing here and is recorded because it was declared.",
        "",
        "## The alternative node key (both readings, declared)",
        "",
        "`match_system(by_structure=...)` decides how the matched sub-network's nodes are "
        "identified. The article quotes the structure-keyed reading, which is the supported one "
        "because a ligand name is not an identifier in this benchmark, and the less favourable "
        "of the two because merging nodes adds cycles. The name-keyed reading is:",
        "",
        "| pool | systems | E | dof | f | chance dof/E | chance isotropic |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for key, label in (("all", "all matched"), ("flagged", "flagged by closure"),
                       ("unflagged", "not flagged")):
        p = name_pooled[key]
        lines.append(
            f"| {label} | {int(p['systems'])} | {int(p['E'])} | {int(p['dof'])} | "
            f"{p['f']:.5f} | {p['chance']:.3f} | {p['chance_iso']:.4f} |")
    key_moved = [(a, b) for a, b in zip(name_rows, rows, strict=True)
                 if a.measurable and b.measurable and abs(a.visible - b.visible) > 1e-12]
    lines += [
        "",
        "Only "
        + (", ".join(
            f"`{b.system}` moves at all: {a.n_nodes} nodes and {a.components} components keyed by "
            f"name against {b.n_nodes} and {b.components} keyed by structure, so dof rises "
            f"{a.dof} -> {b.dof} and f rises {a.visible:.5f} -> {b.visible:.5f}; its worst "
            f"single-edge-deletion margin under the isotropic chance level falls "
            f"{a.loo_worst_margin_iso:.1f}x -> {b.loo_worst_margin_iso:.1f}x"
            for a, b in key_moved) or "no system moves")
        + ".",
        "",
        "## Sensitivity: a protonation-insensitive key (post-hoc, labelled)",
        "",
        "The frozen rule loses `thrombin` entirely to a protonation-state difference. Repeating "
        "the whole measurement with formal charges stripped from both sides before "
        "canonicalization -- a rule chosen AFTER seeing that result, and reported for that reason "
        "as a sensitivity analysis and not as the measurement --",
        "gives:",
        "",
        "| pool | systems | E | dof | f | chance dof/E | chance isotropic |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for key, label in (("all", "all matched"), ("flagged", "flagged by closure"),
                       ("unflagged", "not flagged")):
        p = neutral_pooled[key]
        lines.append(
            f"| {label} | {int(p['systems'])} | {int(p['E'])} | {int(p['dof'])} | "
            f"{p['f']:.5f} | {p['chance']:.3f} | {p['chance_iso']:.4f} |")
    neutral_live = [r for r in neutral_rows if r.measurable]
    entered = [b for a, b in zip(rows, neutral_rows, strict=True)
               if b.measurable and not a.measurable]
    lines += [
        "",
        "and, under the name key instead, the same protonation-insensitive rule gives:",
        "",
        "| pool | systems | E | dof | f | chance dof/E | chance isotropic |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for key, label in (("all", "all matched"), ("flagged", "flagged by closure"),
                       ("unflagged", "not flagged")):
        p = neutral_name_pooled[key]
        lines.append(
            f"| {label} | {int(p['systems'])} | {int(p['E'])} | {int(p['dof'])} | "
            f"{p['f']:.5f} | {p['chance']:.3f} | {p['chance_iso']:.4f} |")
    lines += [
        "",
        f"**The two keys disagree here, and only here.** Under the name key the rule changes the "
        f"pooled figure by nothing at all ({neutral_name_pooled['all']['f']:.5f} against the "
        f"name-keyed exact reading's {name_pooled['all']['f']:.5f}), because the system it lets "
        f"in comes in acyclic. Under the structure key it moves the pool to "
        f"{neutral_pooled['all']['f']:.5f} from {pooled['all']['f']:.5f}, a factor of "
        f"{neutral_pooled['all']['f'] / pooled['all']['f']:.1f}, and the whole of that comes from "
        + (", ".join(f"`{r.system}` entering with E = {r.n_edges} over {r.n_nodes} nodes in "
                     f"{r.components} components, dof = {r.dof} and f = {r.visible:.5f} against "
                     f"its own isotropic chance level of {r.chance_iso:.4f}"
                     for r in entered) or "no system entering")
        + ".",
        "",
        "That is an artefact of the key and not a measurement. The releases of `thrombin` name "
        "one molecule in several poses (`3b` in the JACS set against `3b pose 1` in the water "
        "set), so a 2D structure key merges nodes the release deliberately keeps apart, and the "
        "cycles it creates assert that two poses of two preparations have the same free energy. "
        "The mirror of the `cdk8` collision: keying by name merges distinct molecules, keying by "
        "2D structure merges distinct poses, and a name not being an identifier does not make a "
        "2D structure one. **The measurement never meets this**, because the frozen exact rule "
        "refuses every `thrombin` join on protonation state and the system contributes nothing "
        "to it. All four readings are reported; none replaces another.",
        "",
        "## Same edges, two label sources",
        "",
        "The pooled curated figure and the article's pooled ChEMBL figure are measured on "
        "different edge sets, so their ratio confounds the label source with the sub-network. "
        "Holding the edge set fixed -- the matched sub-network's edges whose two ligands also "
        "carry a committed ChEMBL value -- separates them:",
        "",
        "| system | E | dof | chance dof/E | chance iso | f, curated | f, ChEMBL |"
        " median abs eps, curated | median abs eps, ChEMBL |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for s in same_edge:
        lines.append(
            f"| `{s['system']}` | {s['E']} | {s['dof']} | {s['chance']:.3f} | "
            f"{s['chance_iso']:.4f} | {s['curated']:.5f} | {s['chembl']:.5f} | "
            f"{s['median_cur']:.2f} | {s['median_chembl']:.2f} |")
    live_same = [s for s in same_edge if s["dof"] >= 1]
    if live_same:
        lower = [s for s in live_same if s["curated"] < s["chembl"]]
        ratios = ", ".join(f"`{s['system']}` {s['curated'] / s['chembl']:.1f}x"
                           for s in live_same if s["chembl"] > 0)
        lines += [
            "",
            f"**On genuinely identical edges the curated labels give the HIGHER visible fraction, "
            f"not the lower one**: {len(lower)} of the {len(live_same)} systems that carry a "
            f"cycle come out lower on curated labels, and the curated-to-ChEMBL ratio is "
            f"{ratios}. The direction is the one the decomposition predicts, and it is the "
            "opposite of what the article's Supporting Information asserts. Per-ligand label "
            "noise is a gradient field: the whitened projector annihilates it, so it lands "
            "entirely in the denominator of `f` and pushes the visible fraction DOWN. The "
            "stereo-blind, target-wide ChEMBL join is the noisier label source -- its median "
            "absolute error against experiment is the larger one on these same edges -- so it "
            "should give, and does give, the smaller visible fraction.",
            "",
            "The Supporting Information's sentence that the curated check runs \"at or below what "
            "the ChEMBL labels give on the same edges\" does not survive this. Its curated "
            "figures are the ones in the table above (the `cdk8` row reproduces its 0.0007 "
            "against 0.071 and the `hif2a` row its 0.0006 against 0.147, to every digit quoted), "
            "but the ChEMBL numbers they were compared against are the FULL systems' 0.0065 and "
            "0.0017, measured on 31 and 59 edges rather than on these 14 and 34. Restricted to "
            "the same edges the comparison reverses. What the reversal costs the article is one "
            "sentence of the label-noise rebuttal, not the rebuttal: the curated reading is the "
            "cleaner label source, it is the one that should be believed, and it still puts the "
            "visible fraction two orders of magnitude below its own chance level. It is a "
            "consistency check between two label sources on shared calculated values, not an "
            "independent confirmation of either.",
        ]
    lines += [
        "",
        "## Caveats a referee should hold against this",
        "",
        "- The matched sub-networks are sub-networks. Restricting a network to the edges a "
        "curated table happens to cover removes cycles, and the pooled `dof/E` here "
        f"({all_pool['chance']:.3f}) is well below the whole benchmark's 0.325. A smaller "
        "auditable dimension makes the chance level smaller too, which is why both chance levels "
        "are carried in every row rather than one.",
        "- Curated per-ligand experimental error is a gradient field and lands entirely in the "
        "denominator of `f`, so any label noise biases this measurement toward the article's "
        "conclusion. That is an argument for the curated reading being the cleaner one -- it is "
        "the less noisy label source -- and not an argument that either reading is unbiased.",
        "- **The label-noise floor is of the same size as the errors being measured.** Mixed "
        "public affinity data reproduce across laboratories at 0.68 log10 units per ligand "
        "(Kalliokoski et al., PLOS ONE 2013), which at 1.3642 kcal/mol per log10 unit is 0.93 "
        "kcal/mol per ligand and 1.31 kcal/mol per edge. The median per-edge error against the "
        f"curated labels here runs {min(r.median_abs_eps for r in live):.2f} to "
        f"{max(r.median_abs_eps for r in live):.2f} kcal/mol, at or below that floor. Much of "
        "the denominator of `f` may therefore be label noise rather than force-field error. "
        "Because a per-ligand label error is exactly a gradient field, the whitened projector "
        "annihilates it and it cannot enter the numerator at all -- so it cannot manufacture a "
        "small `f`, but it does mean `f` is a lower bound on the visible share of the "
        "CALCULATION's own error rather than an estimate of it. Bounding that share from above "
        "needs a label source whose noise is known, which public affinity data is not.",
        "- Four systems with a real matched sub-network carry no cycle in it. Their error is not "
        "measured here at all, in either direction.",
        f"- Per-edge error against experiment runs "
        f"{min(r.median_abs_eps for r in live):.2f} to {max(r.median_abs_eps for r in live):.2f} "
        "kcal/mol at the median across these systems, against reported standard errors of "
        f"{min(r.median_se for r in live):.2f} to {max(r.median_se for r in live):.2f}. The "
        "networks are wrong against experiment by several standard errors while closing their "
        "cycles, which is the same picture the four-system measurement gives.",
    ]
    DOC.write_text("\n".join(lines) + "\n")
    print(f"wrote {DOC}")


def main() -> None:
    flagged = ever_flagged()
    memo: dict[str, str | None] = {}
    grid = {(neutralize, by_structure):
            (rs := readings(neutralize, flagged, memo, by_structure=by_structure), pools(rs))
            for neutralize in (False, True) for by_structure in (True, False)}
    rows, pooled = grid[(False, True)]
    same_edge = same_edge_comparison(rows, _ligand_structures(False, memo))
    for r in rows:
        print(f"{r.system:9s} rec {r.recovery:.3f}  E {r.n_edges:3d}  dof {r.dof:3d}  "
              f"f {r.visible:.5f}  chance {r.chance:.4f}  iso {r.chance_iso:.4f}")
    print("pooled:", {k: round(v["f"], 6) for k, v in pooled.items()})
    draw(rows, pooled, same_edge)
    write_table(rows, pooled)
    write_doc(grid, same_edge)


if __name__ == "__main__":
    main()
