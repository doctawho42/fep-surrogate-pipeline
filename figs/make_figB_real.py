"""Fig B (REAL) — decomposed/normalized conformal on REAL FEP residuals, vs a FAIR
feature-conditional (Mondrian-by-target) conformal and flat split-conformal.

The original figs/make_figB.py is a controlled 1-D synthetic where the conditioning axis is
the same variable that builds the sigma decomposition (an oracle DGP). This script is the
honest audit: it runs the *real* EnsembleTrunk on the OpenFF protein-ligand benchmark
(scaffold-disjoint OOD edges, experimental ddG labels) and compares, at equal ~0.90 marginal
coverage:
  - ours  : normalized conformal, half-width = q_norm * sqrt(se_epistemic^2 + sigma_ale^2)
  - split : flat split-conformal, constant half-width
  - mondrian : per-target (feature-conditional) split-conformal — the FAIR adaptive baseline
The metrics: marginal coverage, sharpness (mean half-width), and CONDITIONAL coverage on two
axes — per target (Mondrian's home turf) and per |ddG|-magnitude bin (orthogonal to all).

Honest expectation (and finding): on real OOD residuals the epistemic se is weak relative to
the large irreducible error, so ours is NOT sharper than flat split (the synthetic "2x" does
not transfer); the decomposed sigma buys CONDITIONAL calibration comparable to a Mondrian
conformal that must be told the target partition, at the cost of wider intervals.

Run:  PYTHONPATH=src python figs/make_figB_real.py
"""
from __future__ import annotations

import csv
import pathlib
import sys
import warnings

import matplotlib.pyplot as plt
import numpy as np

warnings.filterwarnings("ignore")
ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from bar.calib import conformal_q  # noqa: E402
from bar.reward import regret_difference_ci  # noqa: E402
from bar.trunk import EnsembleTrunk  # noqa: E402
from paperstyle import (  # noqa: E402
    ALT, INK, OURS, THIRD, figsize, finish, legend, panel, use_paper_style,
)
from rdkit import Chem  # noqa: E402
from rdkit.Chem.Scaffolds import MurckoScaffold  # noqa: E402

DATA = ROOT / "data" / "fep_edges"
ALPHA = 0.10            # target marginal coverage 1 - ALPHA = 0.90
SIGMA_ALE = 0.4         # benchmark label noise ~0.4 kcal/mol
COV = 1.0 - ALPHA
MIN_TARGET_N = 4        # min OOD test edges to report a per-target coverage
MONDRIAN_MIN_CAL = 5    # min cal edges in a target group to use its own conformal q


def load_edges():
    rows = []
    with open(DATA / "all_edges.csv") as f:
        for r in csv.DictReader(f):
            rows.append((r["smiles_a"], r["smiles_b"], float(r["ddg"]), r["target"]))
    return rows


def _scaffold(smi):
    m = Chem.MolFromSmiles(smi)
    return MurckoScaffold.MurckoScaffoldSmiles(mol=m) if m else ""


def scaffold_split(rows, seed):
    scaffs = sorted({_scaffold(a) for a, _b, _d, _t in rows}
                    | {_scaffold(b) for _a, b, _d, _t in rows})
    rng = np.random.default_rng(seed)
    rng.shuffle(scaffs)
    ood = set(scaffs[: max(1, len(scaffs) // 4)])
    train, test = [], []
    for a, b, d, t in rows:
        (test if (_scaffold(a) in ood or _scaffold(b) in ood) else train).append((a, b, d, t))
    return train, test


def _flat_q(abs_res):
    """Flat split-conformal half-width = normalized conformal with sigma=1 (a constant)."""
    return conformal_q(abs_res, np.ones_like(abs_res), ALPHA)


def _cond_cov_error(covered, key):
    """Per-group |coverage - COV|: returns (max, mean) over groups with >= MIN_TARGET_N."""
    errs = []
    for g in sorted(set(key)):
        m = key == g
        if m.sum() >= MIN_TARGET_N:
            errs.append(abs(float(np.mean(covered[m])) - COV))
    return (float(np.max(errs)), float(np.mean(errs))) if errs else (np.nan, np.nan)


def run(seeds=range(8)):
    rows = load_edges()
    acc = {k: [] for k in (
        "marg_ours", "marg_split", "marg_mond",
        "width_ours", "width_split", "width_mond",
        "tgt_max_ours", "tgt_max_split", "tgt_max_mond",
        "tgt_mean_ours", "tgt_mean_split", "tgt_mean_mond",
        "mag_max_ours", "mag_max_split", "mag_max_mond",
        "adapt_corr",  # corr(ours half-width, |residual|) — does ours adapt on real data?
    )}
    for seed in seeds:
        train, test = scaffold_split(rows, seed)
        if len(test) < 12 or len(train) < 40:
            continue
        cut = int(0.75 * len(train))
        fit_e = [(a, b) for a, b, _d, _t in train[:cut]]
        fit_y = np.array([d for _a, _b, d, _t in train[:cut]])
        cal_e = [(a, b) for a, b, _d, _t in train[cut:]]
        cal_y = np.array([d for _a, _b, d, _t in train[cut:]])
        cal_t = np.array([t for _a, _b, _d, t in train[cut:]])
        te_e = [(a, b) for a, b, _d, _t in test]
        te_y = np.array([d for _a, _b, d, _t in test])
        te_t = np.array([t for _a, _b, _d, t in test])

        trunk = EnsembleTrunk().fit(fit_e, fit_y, n_members=8)
        mu_cal, se_cal = trunk.predict(cal_e)
        mu_te, se_te = trunk.predict(te_e)
        res_cal = np.abs(cal_y - mu_cal)
        res_te = np.abs(te_y - mu_te)

        # ours: normalized conformal on decomposed sigma_total
        st_cal = np.sqrt(se_cal ** 2 + SIGMA_ALE ** 2)
        st_te = np.sqrt(se_te ** 2 + SIGMA_ALE ** 2)
        q_norm = conformal_q(res_cal, st_cal, ALPHA)
        hw_ours = q_norm * st_te
        # split: flat conformal
        hw_split = np.full(te_y.shape, _flat_q(res_cal))
        # mondrian: per-target flat conformal (fair feature-conditional baseline)
        q_global = _flat_q(res_cal)
        hw_mond = np.empty(te_y.shape)
        for i, t in enumerate(te_t):
            m = cal_t == t
            hw_mond[i] = _flat_q(res_cal[m]) if m.sum() >= MONDRIAN_MIN_CAL else q_global

        cov_ours = (res_te <= hw_ours).astype(float)
        cov_split = (res_te <= hw_split).astype(float)
        cov_mond = (res_te <= hw_mond).astype(float)

        acc["marg_ours"].append(cov_ours.mean())
        acc["marg_split"].append(cov_split.mean())
        acc["marg_mond"].append(cov_mond.mean())
        acc["width_ours"].append(float(hw_ours.mean()))
        acc["width_split"].append(float(hw_split.mean()))
        acc["width_mond"].append(float(hw_mond.mean()))

        for tag, cov in (("ours", cov_ours), ("split", cov_split), ("mond", cov_mond)):
            mx, mn = _cond_cov_error(cov, te_t)
            acc[f"tgt_max_{tag}"].append(mx)
            acc[f"tgt_mean_{tag}"].append(mn)
            # orthogonal axis: |ddG| magnitude tertiles
            mag = np.digitize(np.abs(te_y), np.quantile(np.abs(te_y), [1 / 3, 2 / 3]))
            mmx, _ = _cond_cov_error(cov, mag)
            acc[f"mag_max_{tag}"].append(mmx)
        if res_te.std() > 1e-9:
            acc["adapt_corr"].append(float(np.corrcoef(hw_ours, res_te)[0, 1]))
    return {k: np.array([v for v in vals if np.isfinite(v)]) for k, vals in acc.items()}


# Semantic colours (paperstyle): OURS = the decomposed-sigma interval, this article's own
# quantity; ALT and THIRD = the two NAMED conformal baselines it is measured against. Neither
# is a random or null baseline, so neither is grey; Mondrian, the fair adaptive rival, takes
# ALT and flat split-conformal takes THIRD.
C_OURS, C_MOND, C_SPLIT = OURS, ALT, THIRD


def main():
    use_paper_style()
    r = run()

    def m(k):
        return float(np.mean(r[k])) if r[k].size else float("nan")

    # sharpness ratio split/ours (>1 => ours sharper; the synthetic claimed ~2.2)
    sharp_ratio = m("width_split") / m("width_ours") if m("width_ours") else float("nan")
    # ours-vs-split per-target conditional-error advantage (lower better => winner iff hi<0)
    md, lo, hi = regret_difference_ci(r["tgt_max_ours"], r["tgt_max_split"])

    fig, (axA, axB) = plt.subplots(1, 2, figsize=figsize(2, 3.1))
    # Panel A: per-target & per-|ddG| max conditional-coverage error (lower = better)
    groups = ["per target", "per |ΔΔG|"]
    ours = [m("tgt_max_ours"), m("mag_max_ours")]
    split = [m("tgt_max_split"), m("mag_max_split")]
    mond = [m("tgt_max_mond"), m("mag_max_mond")]
    x = np.arange(2)
    axA.bar(x - 0.25, ours, 0.25, color=C_OURS, label="ours (decomposed σ)", zorder=2)
    axA.bar(x, mond, 0.25, color=C_MOND, label="Mondrian-by-target", zorder=2)
    axA.bar(x + 0.25, split, 0.25, color=C_SPLIT, label="flat split-conformal", zorder=2)
    axA.set_xticks(x)
    axA.set_xticklabels(groups)
    axA.set_xlim(-0.62, 1.62)
    # headroom for the key, which is the one thing that has to sit inside the panel: the
    # three-bar groups start at zero, so above them is the only clear space there is.
    axA.set_ylim(0.0, max(ours + mond + split) * 1.32)
    axA.tick_params(axis="x", length=0)
    axA.set_ylabel("max |conditional cov − 0.90|")
    panel(axA, "A", "conditional calibration (real OOD)")
    legend(axA, loc="upper right", fontsize=7, borderaxespad=0.2)

    # Panel B: sharpness (mean half-width) — the honest correction to "2x sharper"
    widths = [m("width_ours"), m("width_mond"), m("width_split")]
    axB.bar([0, 1, 2], widths, color=[C_OURS, C_MOND, C_SPLIT], width=0.6, zorder=2)
    axB.set_xticks([0, 1, 2])
    axB.set_xticklabels(["ours", "Mondrian", "split"])
    # bars 0.6 wide on unit centres: the y-spine gets the same 0.4 gap as the bars do
    axB.set_xlim(-0.7, 2.7)
    axB.set_ylim(0.0, max(widths) * 1.22)
    axB.tick_params(axis="x", length=0)
    axB.set_ylabel("mean interval half-width (kcal/mol)")
    panel(axB, "B", "sharpness (lower = tighter)")
    # upper right, clear of the three bars, opaque and in the text ink rather than a wash
    axB.text(0.98, 0.99,
             f"split/ours = {sharp_ratio:.2f}×\n"
             f"(ours {'tighter' if sharp_ratio > 1 else 'wider'})",
             transform=axB.transAxes, ha="right", va="top", fontsize=7, color=INK,
             linespacing=1.45)

    finish(fig, "figB_real_decomposition")

    out = (
        f"# Fig B (REAL) — decomposed conformal on real FEP residuals: honest audit\n\n"
        f"Replaces the synthetic 1-D Fig B headline with the REAL EnsembleTrunk on the OpenFF\n"
        f"protein-ligand benchmark (scaffold-disjoint OOD edges; experimental ΔΔG). At equal\n"
        f"~0.90 MARGINAL coverage we compare normalized conformal on the decomposed σ (**ours**)\n"
        f"vs flat **split**-conformal vs a fair **Mondrian**-by-target conformal. `make figBreal`.\n\n"
        f"## Marginal coverage (target 0.90)\n"
        f"ours {m('marg_ours'):.3f} · split {m('marg_split'):.3f} · Mondrian {m('marg_mond'):.3f} "
        f"({r['marg_ours'].size} seeds).\n\n"
        f"## Sharpness (mean half-width, kcal/mol; lower = tighter)\n"
        f"ours {m('width_ours'):.2f} · Mondrian {m('width_mond'):.2f} · split {m('width_split'):.2f}. "
        f"**split/ours = {sharp_ratio:.2f}×.**\n"
        f"On REAL residuals ours is **{'sharper' if sharp_ratio > 1.02 else 'NOT sharper (wider)'}** "
        f"than flat split — the synthetic figure's 'up to 2.2× sharper vs split-conformal "
        f"(2.0× vs CQR)' does NOT transfer, because the\n"
        f"epistemic se is small relative to the large irreducible OOD error (adapt corr(half-width,\n"
        f"|residual|) = {m('adapt_corr'):+.2f}).\n\n"
        f"## Conditional coverage — max |bin coverage − 0.90| (lower = better)\n"
        f"| axis | ours | Mondrian | split |\n|---|--:|--:|--:|\n"
        f"| per target | {m('tgt_max_ours'):.3f} | {m('tgt_max_mond'):.3f} | {m('tgt_max_split'):.3f} |\n"
        f"| per \\|ΔΔG\\| | {m('mag_max_ours'):.3f} | {m('mag_max_mond'):.3f} | {m('mag_max_split'):.3f} |\n\n"
        f"ours−split per-target max-error diff {md:+.3f}, CI [{lo:+.3f}, {hi:+.3f}] "
        f"(ours better iff hi < 0).\n\n"
        f"## Honest reading\n"
        f"- The decomposed σ buys CONDITIONAL calibration over a flat split-conformal "
        f"(per-target error {m('tgt_max_ours'):.2f} vs {m('tgt_max_split'):.2f}) WITHOUT being told\n"
        f"  the target partition — but a fair Mondrian conformal, which IS told the partition,\n"
        f"  matches it ({m('tgt_max_mond'):.2f}).\n"
        f"- It is **not** sharper on real OOD data (split/ours {sharp_ratio:.2f}×); the synthetic\n"
        f"  'up to 2.2× sharper vs split-conformal (2.0× vs CQR)' was an artifact of an oracle "
        f"DGP whose conditioning axis built σ.\n"
        f"- There is no real force-field-correction (Δ) head; σ is epistemic (ensemble) + a scalar\n"
        f"  aleatoric floor. The 3-way synthetic decomposition is illustrative only.\n"
    )
    (ROOT / "docs" / "results_figB_real.md").write_text(out)
    print("wrote figB_real_decomposition.(pdf|png) + docs/results_figB_real.md")
    print(f"[marginal] ours {m('marg_ours'):.3f} split {m('marg_split'):.3f} mond {m('marg_mond'):.3f}")
    print(f"[sharpness] split/ours {sharp_ratio:.2f}x  (ours width {m('width_ours'):.2f})")
    print(f"[cond per-target max-err] ours {m('tgt_max_ours'):.3f} mond {m('tgt_max_mond'):.3f} "
          f"split {m('tgt_max_split'):.3f}")


if __name__ == "__main__":
    main()
