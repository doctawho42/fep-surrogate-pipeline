"""Fig Hodge: where a perturbation network's error against experiment actually lives.

Theorem 5 splits the per-edge error field into a gradient part (dimension N-c, invisible to any
diagnostic computed from the network's own edges) and a cycle part (dimension dof, all that closure
can ever see). This script runs the three pre-registered items of
``data/openfe_replicates/hodge_prereg.yaml``:

  H3  the auditability map over the whole benchmark (descriptive, no threshold);
  H1  the measured visible fraction of the error against experiment, against its dof/E chance level;
  H2  the influence-ranked repair race, the rule the theorem derives, against random removal.

Run: PYTHONPATH=src python figs/make_figHodge.py   (or `make figHodge`)
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
    ALT,
    INK,
    MUTED,
    OURS,
    REF,
    figsize,
    finish,
    legend,
    panel,
    tint,
    use_paper_style,
)

from bar.closeloop import combine, delta_mue, load_prereg, load_system_edges  # noqa: E402
from bar.hodge import gradient_r2, hodge_split, influence_repair_order  # noqa: E402
from bar.leverage import curl_leverage  # noqa: E402
from bar.qc import _incidence, benjamini_hochberg, chi2_sf, gls_network, repair_order  # noqa: E402

DOC = ROOT / "docs" / "results_figHodge.md"
TABLE = ROOT / "docs" / "tab_auditability.tex"
METRIC_TABLE = ROOT / "docs" / "tab_metrics.tex"
COMBINED = ROOT / "data" / "openfe_replicates" / "combined_pymbar4_edge_data.csv"
GROUNDED = ["cdk8", "hif2a", "p38", "bace"]


def _f(x: str) -> float:
    try:
        return float(x)
    except (TypeError, ValueError):
        return math.nan


def _edge(row: dict, k: int):
    cD, cd = _f(row[f"complex_repeat_{k}_DG (kcal/mol)"]), _f(row[f"complex_repeat_{k}_dDG (kcal/mol)"])
    sD, sd = _f(row[f"solvent_repeat_{k}_DG (kcal/mol)"]), _f(row[f"solvent_repeat_{k}_dDG (kcal/mol)"])
    if any(math.isnan(v) for v in (cD, cd, sD, sd)):
        return None
    return cD - sD, math.sqrt(cd * cd + sd * sd)


def experimental(system: str) -> dict:
    path = ROOT / "data" / "openfe_replicates" / f"affinity_{system}.csv"
    if not path.exists():
        return {}
    return {r["ligand"]: float(r["exp_dg"]) for r in csv.DictReader(path.open())}


def _ever_flagged(by: dict) -> set:
    """Systems the Benjamini--Hochberg test flags on at least one of the three replicates.

    Recomputed here rather than copied from the text, so the table's dagger cannot drift away
    from the flag sets the rest of the analysis produces.
    """
    ever: set = set()
    for replicate in (0, 1, 2):
        names, pvals = [], []
        for name, rows in sorted(by.items()):
            edges = [(r["ligand_A"], r["ligand_B"], *_edge(r, replicate))
                     for r in rows if _edge(r, replicate)]
            if len(edges) < 3:
                continue
            fit = gls_network(edges)
            if fit.dof < 1:
                continue
            names.append(name)
            pvals.append(chi2_sf(fit.chi2, fit.dof))
        flags = benjamini_hochberg(pvals, alpha=0.05)
        ever |= {n for n, on in zip(names, flags, strict=True) if on}
    return ever


def auditability_map() -> dict:
    """H3: dof/E per system and the per-edge detectable-shift resolution over the benchmark."""
    by: dict[str, list] = {}
    for row in csv.DictReader(COMBINED.open()):
        by.setdefault(row["system name"], []).append(row)
    flagged = _ever_flagged(by)
    per_system, resolutions, n_bridges, n_edges = [], [], 0, 0
    for name, rows in sorted(by.items()):
        edges = [(r["ligand_A"], r["ligand_B"], *_edge(r, 0)) for r in rows if _edge(r, 0)]
        if len(edges) < 3:
            continue
        fit = gls_network(edges)
        if fit.dof < 1:
            continue
        h = curl_leverage(edges)
        variance = np.array([e[3] for e in edges]) ** 2
        auditable = h > 1e-9
        delta_star = np.full(len(edges), np.inf)
        delta_star[auditable] = np.sqrt(variance[auditable] / h[auditable])
        finite = delta_star[np.isfinite(delta_star)]
        per_system.append({"system": name, "E": len(edges), "N": len(fit.nodes), "dof": fit.dof,
                           "c": len(fit.nodes) - (len(edges) - fit.dof),
                           "visible_dim_frac": fit.dof / len(edges),
                           "bridges": int((~auditable).sum()),
                           "median_delta": float(np.median(finite)) if finite.size else float("nan"),
                           "flagged": name in flagged})
        resolutions.extend(delta_star.tolist())
        n_bridges += int((~auditable).sum())
        n_edges += len(edges)
    res = np.array(resolutions)
    return {"per_system": per_system, "delta_star": res, "n_edges": n_edges,
            "n_bridges": n_bridges,
            "frac_at": {d: float(np.mean(res <= d)) for d in (0.5, 1.0, 2.0)},
            "quartiles": [float(q) for q in np.nanpercentile(res[np.isfinite(res)], [25, 50, 75])],
            "pooled_visible_dim": sum(s["dof"] for s in per_system) / n_edges}


def release_decomposition() -> list[dict]:
    """Split each system's closure statistic into within-release and across-release cycles.

    The benchmark names a system by target, and eight of its systems union edges from more than one
    separately released data set. Where those releases share ligands, the union carries cycles that
    exist in neither release alone. The cycle space of the union contains the two releases' cycle
    spaces as orthogonal subspaces, so the split below is exact: the chi-squared and the degrees of
    freedom both add.
    """
    by: dict[str, list] = {}
    for row in csv.DictReader(COMBINED.open()):
        by.setdefault(row["system name"], []).append(row)
    out = []
    for name, rows in sorted(by.items()):
        per_release: dict[str, list] = {}
        edges = []
        for r in rows:
            val = _edge(r, 0)
            if not val:
                continue
            edge = (r["ligand_A"], r["ligand_B"], *val)
            edges.append(edge)
            per_release.setdefault(r["dataset_name"], []).append(edge)
        if len(per_release) < 2 or len(edges) < 3:
            continue
        union = gls_network(edges)
        if union.dof < 1:
            continue
        within_chi2 = within_dof = 0.0
        for block in per_release.values():
            if len(block) < 2:
                continue
            fit = gls_network(block)
            if fit.dof < 1:
                continue
            within_chi2 += fit.chi2
            within_dof += fit.dof
        across_dof = union.dof - int(within_dof)
        across_chi2 = union.chi2 - within_chi2
        out.append({"system": name, "releases": len(per_release),
                    "union_chi2": union.chi2, "union_dof": union.dof,
                    "union_reduced": union.reduced_chi2,
                    "within_chi2": within_chi2, "within_dof": int(within_dof),
                    "across_chi2": across_chi2, "across_dof": across_dof,
                    "across_reduced": across_chi2 / across_dof if across_dof > 0 else float("nan"),
                    "across_p": chi2_sf(across_chi2, across_dof) if across_dof > 0 else float("nan")})
    return out


def where_the_error_lives() -> list[dict]:
    """H1: the visible fraction of the measured error against experiment, per grounded system."""
    out = []
    for system in GROUNDED:
        exp = experimental(system)
        edges = load_system_edges(system, exp)
        if len(edges) < 3:
            continue
        eps = np.array([ddg - (exp[b] - exp[a]) for a, b, ddg, _ in edges])
        se = np.array([e[3] for e in edges])
        split = hodge_split(edges, eps)
        r2, adjusted, chance = gradient_r2(edges, eps)
        # a median-based variant, to show the fraction is not driven by a few large residuals
        robust = float(np.median((split.cycle / se) ** 2) / np.median((eps / se) ** 2))
        weights = np.diag(1.0 / se ** 2)
        whitened = np.diag(1.0 / se) @ _incidence(edges)[1]
        proj = np.eye(len(edges)) - whitened @ np.linalg.pinv(whitened.T @ whitened) @ whitened.T
        chance_iso = float(np.trace(proj @ weights) / np.trace(weights))
        out.append({"system": system, "E": len(edges), "dof": split.dof,
                    "visible": split.visible_fraction, "chance": chance,
                    "chance_iso": chance_iso, "robust": robust,
                    "r2": r2, "r2_adj": adjusted,
                    "rms_eps": float(np.sqrt((eps ** 2).mean())),
                    "median_abs_eps": float(np.median(np.abs(eps))),
                    "median_se": float(np.median(se)),
                    "cycle_sq": float(np.sum((split.cycle / se) ** 2)),
                    "total_sq": float(np.sum((eps / se) ** 2))})
    return out


def _projectors(edges):
    """(whitened projector, unweighted projector, weight matrix) for one edge list."""
    se = np.array([e[3] for e in edges])
    inc = _incidence(edges)[1]
    whitened = np.diag(1.0 / se) @ inc
    proj_w = np.eye(len(edges)) - whitened @ np.linalg.pinv(whitened.T @ whitened) @ whitened.T
    proj_u = np.eye(len(edges)) - inc @ np.linalg.pinv(inc.T @ inc) @ inc.T
    return proj_w, proj_u, np.diag(1.0 / se ** 2)


def _three_readings(edges, eps):
    """The visible fraction and its chance level under the three metrics the article quotes.

    whitened   fraction and chance both in the variance-weighted metric the closure statistic
               lives in; the chance level is dof/E.
    isotropic  the same whitened fraction against the chance level an error isotropic in
               kcal/mol would produce under the test's own projector, tr(Pi W)/tr(W).
    unweighted the projector RECOMPUTED in the unweighted metric, so both the fraction and its
               dof/E chance level are taken with no variance weighting at all.
    """
    se = np.array([e[3] for e in edges])
    proj_w, proj_u, weights = _projectors(edges)
    split = hodge_split(edges, eps)
    tilde = eps / se
    vis_w = float((tilde @ proj_w @ tilde) / (tilde @ tilde))
    vis_u = float((eps @ proj_u @ eps) / (eps @ eps))
    chance = split.dof / len(edges)
    chance_iso = float(np.trace(proj_w @ weights) / np.trace(weights))
    return {"dof": split.dof, "E": len(edges),
            "whitened": (vis_w, chance), "isotropic": (vis_w, chance_iso),
            "unweighted": (vis_u, chance)}


def metric_readings() -> list[dict]:
    """Item-8 bookkeeping: four grounded systems by three metric conventions, each with its own
    chance level and its leave-one-edge-out range. Both the fraction and the chance level are
    recomputed on every deletion, since deleting an edge moves the projector as well."""
    out = []
    for system in GROUNDED:
        exp = experimental(system)
        edges = load_system_edges(system, exp)
        if len(edges) < 3:
            continue
        eps = np.array([ddg - (exp[b] - exp[a]) for a, b, ddg, _ in edges])
        full = _three_readings(edges, eps)
        loo: dict[str, list] = {k: [] for k in ("whitened", "isotropic", "unweighted")}
        for k in range(len(edges)):
            sub = [e for i, e in enumerate(edges) if i != k]
            sub_eps = np.array([eps[i] for i in range(len(edges)) if i != k])
            r = _three_readings(sub, sub_eps)
            if r["dof"] < 1:
                continue
            for key in loo:
                v, c = r[key]
                if v > 0:
                    loo[key].append((v, c))
        row = {"system": system, "E": full["E"], "dof": full["dof"]}
        for key in ("whitened", "isotropic", "unweighted"):
            v, c = full[key]
            vs = [x[0] for x in loo[key]]
            margins = [c_ / v_ for v_, c_ in loo[key]]
            row[key] = {"visible": v, "chance": c, "margin": c / v,
                        "loo_lo": min(vs), "loo_hi": max(vs), "loo_worst_margin": min(margins)}
        out.append(row)
    return out


def pooled_readings() -> dict:
    """Pooled visible fraction in the whitened and unweighted metrics, with and without cdk8.

    cdk8's flag is a cross-release join rather than a statement about either release's
    simulations (Section 4), so the pooled figure is reported both ways at every point of use.
    """
    acc = {}
    for system in GROUNDED:
        exp = experimental(system)
        edges = load_system_edges(system, exp)
        if len(edges) < 3:
            continue
        eps = np.array([ddg - (exp[b] - exp[a]) for a, b, ddg, _ in edges])
        se = np.array([e[3] for e in edges])
        proj_w, proj_u, _ = _projectors(edges)
        tilde = eps / se
        acc[system] = {"cw": float(tilde @ proj_w @ tilde), "tw": float(tilde @ tilde),
                       "cu": float(eps @ proj_u @ eps), "tu": float(eps @ eps)}
    keys = list(acc)
    no_cdk8 = [k for k in keys if k != "cdk8"]

    def pool(ks, a, b):
        return sum(acc[k][a] for k in ks) / sum(acc[k][b] for k in ks)

    return {"whitened": pool(keys, "cw", "tw"), "whitened_no_cdk8": pool(no_cdk8, "cw", "tw"),
            "unweighted": pool(keys, "cu", "tu"), "unweighted_no_cdk8": pool(no_cdk8, "cu", "tu")}


def pooled_r2_adj(lives: list[dict]) -> dict:
    """The pooled degrees-of-freedom-adjusted gradient-model R^2 the pre-registration froze.

    H1's success rule is ``pooled R2_adj >= 0.50``, with R2_adj = 1 - (1 - R^2) * E / (E - rank)
    and R^2 = 1 - f. Pooling is over the four grounded sub-networks: the visible fraction is the
    ratio of summed squared norms, and E and rank are summed likewise, so the adjustment carries
    the pooled residual degrees of freedom rather than an average of per-system adjustments. The
    adjudicator released with the first version of this script tested ``1 - f_pooled >= 0.5``,
    the UNADJUSTED pooled criterion, which is a weaker rule than the one frozen.
    """
    cycle = sum(r["cycle_sq"] for r in lives)
    total = sum(r["total_sq"] for r in lives)
    n_edges = sum(r["E"] for r in lives)
    rank = sum(r["E"] - r["dof"] for r in lives)
    f = cycle / total
    r2 = 1.0 - f
    adj = 1.0 - (1.0 - r2) * n_edges / (n_edges - rank)
    return {"f": f, "r2": r2, "r2_adj": adj, "E": n_edges, "rank": rank,
            "resid_dof": n_edges - rank}


def write_metric_table(readings: list[dict], pooled: dict) -> None:
    """Emit the three-convention bookkeeping as LaTeX, for the supplement to input."""
    body = []
    for r in readings:
        name = r["system"].replace("_", r"\_")
        first = True
        for key, label in (("whitened", "whitened"), ("isotropic", "isotropic, kcal/mol"),
                           ("unweighted", "unweighted")):
            d = r[key]
            lead = rf"\texttt{{{name}}} & ${r['E']}$ & ${r['dof']}$" if first else "& &"
            body.append(
                rf"{lead} & {label} & ${d['visible']:.4f}$ & ${d['chance']:.4f}$ & "
                rf"${d['margin']:.1f}\times$ & $[{d['loo_lo']:.4f},{d['loo_hi']:.4f}]$ & "
                rf"${d['loo_worst_margin']:.2f}\times$ \\")
            first = False
        body.append(r"\addlinespace[2pt]")
    METRIC_TABLE.write_text("\n".join([
        r"\begin{table}[tp]\centering",
        r"\caption{\textbf{The visible fraction under three metric conventions.} For each grounded "
        r"system: edges $E$, independent cycles $\nu$, the auditable share of that system's error "
        r"against experiment, the chance level it is measured against, their ratio, the range the "
        r"fraction takes when one edge is deleted at a time, and the worst margin over those "
        r"deletions with both the fraction and the chance level recomputed on each. "
        r"\emph{Whitened} takes both in the variance-weighted metric the closure statistic lives "
        r"in, against $\nu/E$. \emph{Isotropic} keeps that fraction and measures it against the "
        r"chance level an error isotropic in kcal/mol would produce under the same projector, "
        r"$\mathrm{tr}(\Pi W)/\mathrm{tr}(W)$. \emph{Unweighted} recomputes the projector in the "
        r"unweighted metric, so neither the fraction nor its $\nu/E$ chance level carries any "
        r"variance weighting. Every cell of every convention, and every single-edge deletion, "
        r"sits below its own chance level; the weakest is \texttt{bace} isotropic under deletion. "
        rf"Pooled over the four systems the fraction is ${pooled['whitened']:.4f}$ whitened and "
        rf"${pooled['unweighted']:.4f}$ unweighted, and ${pooled['whitened_no_cdk8']:.4f}$ and "
        rf"${pooled['unweighted_no_cdk8']:.4f}$ with \texttt{{cdk8}} dropped. "
        r"Produced by \texttt{make figHodge}.}",
        r"\label{tab:metrics}",
        r"\small",
        r"\begin{tabular}{@{}lrrlrrrlr@{}}",
        r"\toprule",
        r"system & $E$ & $\nu$ & metric & visible & chance & margin & leave-one-out & worst \\",
        r"\midrule",
        *body,
        r"\bottomrule",
        r"\end{tabular}",
        r"\end{table}",
    ]) + "\n")
    print(f"wrote {METRIC_TABLE}")


def repair_race(n_perm: int, target_rchi2: float) -> dict:
    """H2: influence-ranked removal against random, at the removal count the |z| rule chooses."""
    influence, residual = [], []
    for system in GROUNDED:
        exp = experimental(system)
        edges = load_system_edges(system, exp)
        removed_z, _ = repair_order(edges, target_reduced_chi2=target_rchi2)
        k = len(removed_z)
        removed_influence = influence_repair_order(edges, k=k)
        gain_z = delta_mue(edges, removed_z, exp)
        gain_influence = delta_mue(edges, removed_influence, exp)
        rng = np.random.default_rng(0)
        draws = np.array([delta_mue(edges, list(rng.choice(len(edges), k, replace=False)), exp)
                          for _ in range(n_perm)])
        draws = draws[~np.isnan(draws)]
        cut = float(np.percentile(draws, 95))
        for gain, bucket in ((gain_influence, influence), (gain_z, residual)):
            bucket.append({"system": system, "k": k, "guided": float(gain),
                           "random_mean": float(draws.mean()),
                           "p": float(np.mean(draws >= gain)), "below_5pct": bool(gain > cut)})
    wins = sum(1 for a, b in zip(influence, residual, strict=True) if a["guided"] > b["guided"])
    n = len(influence)
    return {"influence": influence, "residual": residual, "n_perm": n_perm,
            "combined_influence": combine(influence, n_perm),
            "combined_residual": combine(residual, n_perm),
            "head_to_head_wins": wins,
            "head_to_head_p": sum(math.comb(n, i) for i in range(wins, n + 1)) / 2 ** n}


def main() -> None:
    prereg = load_prereg()
    audit = auditability_map()
    releases = release_decomposition()
    lives = where_the_error_lives()
    race = repair_race(prereg.n_perm, prereg.target_reduced_chi2)

    use_paper_style()
    fig, (axA, axB, axC) = plt.subplots(
        1, 3, figsize=figsize(3, 2.85), gridspec_kw={"width_ratios": [1.18, 1.0, 1.12]})
    fig.subplots_adjust(left=0.075, right=0.995, bottom=0.155, top=0.83, wspace=0.44)

    # A: how fine an error each edge could detect, over the whole benchmark
    finite = np.sort(audit["delta_star"][np.isfinite(audit["delta_star"])])
    frac = np.arange(1, len(finite) + 1) / audit["n_edges"]
    ceiling = 1 - audit["n_bridges"] / audit["n_edges"]
    axA.set_xscale("log")
    axA.set_xlim(finite[0] * 0.75, finite[-1] * 1.4)
    axA.set_ylim(0.0, 1.0)
    axA.axhspan(ceiling, 1.0, color=tint(REF, 0.86), lw=0, zorder=0)
    axA.axhline(ceiling, color=REF, lw=0.8, ls="--", zorder=1)
    # two lines, so the label stops well short of the curve instead of being crossed by it
    axA.annotate(f"{audit['n_bridges']} bridge edges:\nnever auditable",
                 xy=(0.03, ceiling), xycoords=("axes fraction", "data"),
                 xytext=(0, -3), textcoords="offset points",
                 ha="left", va="top", fontsize=7.5, color=REF, linespacing=1.3)
    # Each guide stops at the crossing it marks, and each label sits on the open side of
    # that crossing -- 0.5 to the left of its guide and above the curve, 1.0 to the right
    # of its guide and below it. Set two lines wide, no label steps over any line.
    for d, ha, dx, dy, va in ((0.5, "right", -4, 3, "bottom"), (1.0, "left", 4, -3, "top")):
        y_at = audit["frac_at"][d]
        axA.plot([d, d], [0.0, y_at], color=REF, lw=0.7, ls=":", zorder=1)
        axA.annotate(f"{y_at:.0%} reach\n$\\delta^\\ast_e\\leq{d}$",
                     (d, y_at), xytext=(dx, dy), textcoords="offset points",
                     ha=ha, ma=ha, va=va, fontsize=7.5, color=REF, linespacing=1.3)
    axA.step(finite, frac, color=OURS, lw=1.6, where="post", zorder=3)
    # upper left, the one quadrant the curve never enters: anchored bottom-right this block
    # reached back to the x = 1 guide and touched it, the collision class this pass removes.
    axA.annotate("80% power:\n2.8–4.7×",
                 xy=(0.03, 0.72), xycoords="axes fraction",
                 ha="left", ma="left", va="top", fontsize=7.5, color=REF,
                 linespacing=1.3)
    axA.set_xlabel(r"$\delta^\ast_e$ (kcal/mol)")
    axA.set_ylabel(f"fraction of the {audit['n_edges']} edges")
    panel(axA, "A", "what the audit can resolve",
          subtitle="unit noncentrality, not a threshold")

    # B: the measured error against experiment, against the chance level
    xs = np.arange(len(lives))
    w = 0.26
    # both chance levels are null baselines (MUTED); they are told apart by fill weight
    # and a hatch, not by two greys a print reader would have to discriminate
    axB.bar(xs - 0.28, [r["chance"] for r in lives], w, color=MUTED, lw=0,
            label=r"chance, $\nu/E$")
    axB.bar(xs, [r["chance_iso"] for r in lives], w, facecolor=tint(MUTED, 0.72),
            edgecolor=MUTED, lw=0.6, hatch="///", label="chance, isotropic")
    axB.bar(xs + 0.28, [r["visible"] for r in lives], w, color=OURS, lw=0,
            label="error vs experiment")
    axB.set_yscale("log")
    axB.set_ylim(4e-4, 6.0)
    axB.set_yticks([1e-3, 1e-2, 1e-1, 1e0])
    axB.set_xlim(-0.6, len(lives) - 0.4)
    axB.set_xticks(xs)
    axB.set_xticklabels([r["system"] for r in lives])
    axB.set_ylabel("auditable fraction\nof the error")
    pooled = sum(r["cycle_sq"] for r in lives) / sum(r["total_sq"] for r in lives)
    panel(axB, "B", "where the error lives",
          subtitle=f"{pooled:.1%} visible, variance-weighted")
    legend(axB, loc="upper left", ncol=1)

    # C: the race the theorem motivated, and the rule it was meant to beat
    # colour by quantity family: |z|-ranked is built from the calibrated residual (OURS),
    # influence-ranked is the theorem-derived alternative (ALT), random is a null (MUTED)
    axC.bar(xs - 0.28, [r["guided"] for r in race["influence"]], w, color=ALT, lw=0,
            label=f"influence-ranked, $p={race['combined_influence']['sign_p']:.2f}$")
    axC.bar(xs, [r["guided"] for r in race["residual"]], w, color=OURS, lw=0,
            label=f"$|z|$-ranked, $p={race['combined_residual']['sign_p']:.2f}$")
    axC.bar(xs + 0.28, [r["random_mean"] for r in race["influence"]], w, color=MUTED, lw=0,
            label="random mean")
    axC.axhline(0, color=INK, lw=0.7, zorder=2)
    axC.set_ylim(-0.68, 0.48)
    axC.set_xlim(-0.6, len(race["influence"]) - 0.4)
    axC.set_xticks(xs)
    axC.set_xticklabels([r["system"] for r in race["influence"]])
    axC.set_ylabel("$\\Delta$MUE vs experiment\n(higher is better)")
    panel(axC, "C", "repair by removal", subtitle="neither rule beats random")
    legend(axC, loc="upper left")

    finish(fig, "figHodge_where_the_error_lives")
    write_table(audit)
    readings = metric_readings()
    pooled_all = pooled_readings()
    write_metric_table(readings, pooled_all)
    write_doc(audit, lives, race, pooled, releases, readings, pooled_all)


def write_table(audit: dict) -> None:
    """Emit the per-system auditability table as LaTeX, for the Supporting Information to input.

    Generated rather than transcribed: a hand-copied table is the one artefact in this repository
    that can silently disagree with the analysis that produced it.
    """
    rows = sorted(audit["per_system"], key=lambda r: r["visible_dim_frac"])
    body = []
    for r in rows:
        name = r["system"].replace("_", r"\_")
        dagger = r"$^\dagger$" if r["flagged"] else ""
        body.append(
            rf"\texttt{{{name}}}{dagger} & ${r['N']}$ & ${r['E']}$ & ${r['c']}$ & ${r['dof']}$ & "
            rf"${r['visible_dim_frac']:.2f}$ & ${r['bridges']}$ & ${r['median_delta']:.2f}$ \\")
    half = (len(rows) + 1) // 2
    left, right = body[:half], body[half:]
    right += [r"& & & & & & & \\"] * (len(left) - len(right))
    paired = [f"{a[:-3]} & {b}" for a, b in zip(left, right, strict=True)]
    head = (r"system & $N$ & $E$ & $c$ & $\nu$ & $\nu/E$ & br & $\tilde\delta^\ast$")
    TABLE.write_text("\n".join([
        r"\begin{table}[tp]\centering",
        r"\caption{\textbf{What each benchmark system can be audited for.} For the $48$ systems "
        r"carrying at least one independent cycle: $N$ ligands, $E$ edges, $c$ connected "
        r"components, $\nu$ independent cycles, the auditable fraction $\nu/E$ of that system's "
        r"error space, the number of bridge edges (br), which lie in no cycle and carry no "
        r"evidence at any magnitude, and the median shift $\tilde\delta^\ast$ in kcal/mol at "
        r"which the closure statistic carries unit noncentrality, over that system's auditable "
        r"edges. Rows are ordered by auditable fraction ascending. A dagger marks the eight "
        r"systems the Benjamini--Hochberg test flags on at least one of the three replicates. "
        r"Values are replicate-0 and are produced by the same run of the released analysis code "
        r"that draws Figure~\ref{fig:Hodge}. Note that $\tilde\delta^\ast$ is not a detection "
        r"threshold: reaching $80\%$ power at $\alpha=0.05$ needs a shift of $2.8$ to $4.7$ "
        r"times it, and it is carried on the raw reported standard errors.}",
        r"\label{tab:auditability}",
        r"\footnotesize",
        r"\setlength{\tabcolsep}{3.2pt}",
        r"\begin{tabular}{@{}lrrrrrrr@{\hspace{7pt}}lrrrrrrr@{}}",
        r"\toprule",
        head + " & " + head + r" \\",
        r"\midrule",
        *paired,
        r"\bottomrule",
        r"\end{tabular}",
        r"\end{table}",
    ]) + "\n")
    print(f"wrote {TABLE}")


def write_doc(audit: dict, lives: list, race: dict, pooled: float, releases: list,
              readings: list, pooled_all: dict) -> None:
    lines = [
        "# Fig Hodge: where a perturbation network's error against experiment lives",
        "",
        "Generated by `figs/make_figHodge.py` (`make figHodge`). Pre-registered in",
        "`data/openfe_replicates/hodge_prereg.yaml`, frozen before any number below was computed.",
        "",
        "## H3 auditability map (descriptive)",
        "",
        f"- {audit['n_edges']} edges over {len(audit['per_system'])} systems carrying at least one cycle.",
        f"- {audit['n_bridges']} bridge edges ({audit['n_bridges'] / audit['n_edges']:.1%}) lie in no",
        "  cycle and are un-auditable at any magnitude.",
        f"- Auditable dimensions: {sum(s['dof'] for s in audit['per_system'])} of {audit['n_edges']}",
        f"  ({audit['pooled_visible_dim']:.1%}); the rest of the error space is node-consistent and",
        "  invisible to every closure-based diagnostic.",
        f"- Detectable shift `delta*` quartiles: {audit['quartiles'][0]:.2f} / "
        f"{audit['quartiles'][1]:.2f} / {audit['quartiles'][2]:.2f} kcal/mol.",
        f"- Fraction of edges able to detect an error of 0.5 kcal/mol: {audit['frac_at'][0.5]:.1%};",
        f"  of 1.0 kcal/mol: {audit['frac_at'][1.0]:.1%}.",
        "",
        "## H1 where the measured error lives (pre-registered quantification)",
        "",
        "| system | edges | dof | chance dof/E | visible | median-based | adjusted R2 | median abs eps |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for r in lives:
        lines.append(
            f"| `{r['system']}` | {r['E']} | {r['dof']} | {r['chance']:.3f} | {r['visible']:.4f} "
            f"| {r['robust']:.4f} | {r['r2_adj']:.3f} | {r['median_abs_eps']:.2f} |")
    # H1's frozen rule is on the ADJUSTED pooled R^2, not on 1 - f_pooled.
    pool_adj = pooled_r2_adj(lives)
    passed = pool_adj["r2_adj"] >= 0.50 and sum(1 for r in lives if r["r2_adj"] > 0) >= 3
    lines += [
        "",
        "",
        "## The visible fraction under three metric conventions",
        "",
        "Both the fraction and its chance level are recomputed on every single-edge deletion.",
        "",
        "| system | E | dof | metric | visible | chance | margin | leave-one-out | worst |",
        "|---|---:|---:|---|---:|---:|---:|---|---:|",
    ] + [
        f"| `{r['system']}` | {r['E']} | {r['dof']} | {k} | {r[k]['visible']:.4f} | "
        f"{r[k]['chance']:.4f} | {r[k]['margin']:.1f}x | "
        f"[{r[k]['loo_lo']:.4f}, {r[k]['loo_hi']:.4f}] | {r[k]['loo_worst_margin']:.2f}x |"
        for r in readings for k in ("whitened", "isotropic", "unweighted")
    ] + [
        "",
        f"Pooled adjusted R2 (H1's frozen criterion): {pool_adj['r2_adj']:.4f} from pooled "
        f"f = {pool_adj['f']:.6f} over E = {pool_adj['E']} edges and rank {pool_adj['rank']}, "
        f"residual dof {pool_adj['resid_dof']}; the criterion is >= 0.50.",
        "",
        f"Pooled: whitened {pooled_all['whitened']:.4f} "
        f"({pooled_all['whitened_no_cdk8']:.4f} without cdk8), unweighted "
        f"{pooled_all['unweighted']:.4f} ({pooled_all['unweighted_no_cdk8']:.4f} without cdk8).",
        "",
        f"Pooled visible fraction {pooled:.4f}; the pre-registered criterion (pooled adjusted R2",
        f"at least 0.50 and positive on at least three systems) is {'met' if passed else 'not met'}.",
        "",
        f"The visible fraction runs below its own chance level by a factor of "
        f"{min(r['chance'] / r['visible'] for r in lives):.0f} to "
        f"{max(r['chance'] / r['visible'] for r in lives):.0f}. The median-based column repeats",
        "the calculation on medians rather than sums, so the result is not carried by a few",
        f"large residuals. Per-edge error against experiment runs "
        f"{min(r['median_abs_eps'] for r in lives):.2f} to "
        f"{max(r['median_abs_eps'] for r in lives):.2f} kcal/mol at the median, against reported",
        f"standard errors of {min(r['median_se'] for r in lives):.2f} to "
        f"{max(r['median_se'] for r in lives):.2f}. These networks are wrong against experiment",
        "by several standard errors while closing their cycles, and the discrepancy sits almost",
        "entirely in the direction closure cannot see.",
        "",
        "Declared confound, from the pre-registration: experimental measurement error and the",
        "stereo-blind ChEMBL join are per-ligand and therefore gradient fields too, so this design",
        "shows the dominant part of the error is node-consistent, not that it is force-field bias.",
        "",
        "## H2 influence-ranked repair (pre-registered, second look at the same data)",
        "",
        "| system | K | influence-ranked | \\|z\\|-ranked | random mean | p (influence) |",
        "|---|---|---|---|---|---|",
    ]
    for a, b in zip(race["influence"], race["residual"], strict=True):
        lines.append(f"| `{a['system']}` | {a['k']} | {a['guided']:+.4f} | {b['guided']:+.4f} "
                     f"| {a['random_mean']:+.4f} | {a['p']:.3f} |")
    inf_c, res_c = race["combined_influence"], race["combined_residual"]
    lines += [
        "",
        f"Pre-registered sign test: p = {inf_c['sign_p']:.3f} for the influence-ranked rule and "
        f"{res_c['sign_p']:.3f} for the",
        f"residual-ranked one. The influence rule wins on {race['head_to_head_wins']} of "
        f"{len(race['influence'])} systems (paired sign p = {race['head_to_head_p']:.3f}). "
        f"Verdict: {inf_c['verdict']}.",
        "",
        "The pre-registration also names a one-sided Stouffer combination, and for the",
        "influence-ranked arm it does not return a usable number. A permutation p is bounded by the",
        f"number of draws, and `p38` lands on exactly 1.000: all {race['n_perm']} random removals did at",
        "least as well as the guided one. Stouffer maps that to an infinite z, so the combination",
        "returns exactly 1.000 whatever the other three systems say. That value is the resolution of",
        "the permutation grid speaking, not a combined p-value, and it is reported here as not",
        "computable rather than as a result. Clipping the inputs to the grid's own range gives",
        f"{inf_c['stouffer_p_clipped']:.3f} as a sensitivity, which is not extreme in either",
        "direction. The residual-ranked arm has no saturated input and its Stouffer combination",
        f"stands unchanged at {res_c['stouffer_p']:.3f}. The decision rule is unaffected, since a",
        "value that cannot be computed fails the same threshold the degenerate one did. What the",
        "saturation hides is read off directly below rather than through a combined p-value.",
        "",
        f"On `p38` the influence rule is not merely null: its guided change of "
        f"{race['influence'][2]['guided']:+.4f} falls below every one of the "
        f"{race['n_perm']} random draws, so the opposite tail is at most 1/{race['n_perm']+1}. That",
        "tail is post hoc, it is one system of four, and the magnitude is convention-dependent",
        "(-0.612 greedy, -0.232 one-shot); the sign survives both. The two edges removed there",
        "carry curl-leverage 0.099 and 0.388, so this is not the near-bridge degenerate branch.",
        "",
        "This is a null overall, and no better than the rule it was meant to improve on. H1",
        "explains it.",
        "A repair rule reads the visible part of the error, and the visible part carries a few per",
        "cent of what is wrong; dividing the residual by a small curl-leverage, which is what the",
        "influence form does, amplifies noise wherever the error is not concentrated on a single",
        "edge. H1 says it is not. The theorem's reach on these data is explanation, not repair.",
        "",
        "## Released-data-set composition (descriptive, found while checking a referee report)",
        "",
        "The benchmark names a system by target, and eight systems union edges from more than one",
        "separately released data set. Where the releases share ligands the union carries cycles",
        "present in neither release alone; the split below is exact, since the releases' cycle",
        "spaces sit orthogonally inside the union's.",
        "",
        "| system | releases | union chi2/dof | within-release | across-release | reduced | p |",
        "|---|---|---|---|---|---|---|",
    ]
    for r in releases:
        lines.append(
            f"| `{r['system']}` | {r['releases']} | {r['union_chi2']:.1f}/{r['union_dof']} "
            f"| {r['within_chi2']:.1f}/{r['within_dof']} | {r['across_chi2']:.1f}/{r['across_dof']} "
            f"| {r['across_reduced']:.2f} | {r['across_p']:.1e} |")
    hot = [r for r in releases if r["across_dof"] > 0 and r["across_p"] < 0.05]
    lines += [
        "",
        "Only `" + "`, `".join(r["system"] for r in hot) + "` carries across-release cycles that",
        "matter. In six of the eight the releases share no ligand and form disjoint components, so",
        "there are no such cycles at all, and in `tyk2` the single one is unremarkable.",
    ]
    DOC.write_text("\n".join(lines) + "\n")
    print(f"wrote {DOC}")


if __name__ == "__main__":
    main()
