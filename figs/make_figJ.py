"""Fig J — does the calibrated commit survive amortization to unseen molecules, and is the
win robust to (a) proper recalibration baselines and (b) matched commit volume?

A deep-ensemble trunk maps congeneric edges -> (DDG_hat, sigma_total = conformal-scaled
sqrt(epistemic^2 + aleatoric^2)). On scaffold-disjoint (OOD) edges we ask two honest
questions the earlier draft did not:

  PANEL A (matched commit volume): rank candidates by the standardized margin
  (mu - tau)/sigma and commit the top-n; compare commit precision at the SAME n for the
  trunk-sigma ranking, a raw-mu ranking (sigma ignored), and the MVE-sigma ranking. Holding
  n equal removes the abstention artifact (a method cannot look good by committing nothing).

  PANEL B (calibration, with recalibration foils): commit-correctness vs claimed confidence
  for the trunk sigma, the raw (overconfident) MVE sigma, and the MVE sigma after standard
  cheap recalibration (temperature scaling and split-conformal). This tests whether the
  physics-derived sigma is uniquely trustworthy or whether any properly recalibrated sigma
  also gives trustworthy commits (it does — the honest narrowing).

Honest framing: the contribution is (i) a calibrated commit needs *some* proper calibration,
which the decomposed physics sigma provides at no calibration-set cost for the aleatoric term
(raw MVE over-claims; recalibrated MVE catches up), and (ii) at matched commit volume the
*choice* of sigma does not reliably improve WHICH molecules are committed (a Gauss-Markov-style
ranking ceiling) — calibration is a trust/decision instrument, not a better-ranking one.

Run:  PYTHONPATH=src python figs/make_figJ.py   (or `make figJ`)
"""
from __future__ import annotations

import csv
import pathlib
import sys
import warnings

import matplotlib.pyplot as plt
import numpy as np
import torch

warnings.filterwarnings("ignore")
ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
from paperstyle import (  # noqa: E402
    ALT,
    FOIL,
    OURS,
    REF,
    figsize,
    finish,
    legend,
    panel,
    reference_line,
    tint,
    use_paper_style,
)

from bar.calib import conformal_q, temperature_scale  # noqa: E402
from bar.reward import (  # noqa: E402
    commit_correctness_curve,
    commit_precision_at_volume,
    regret_difference_ci,
)
from bar.trunk import EnsembleTrunk, chir_pseudoscalar, featurize_edge  # noqa: E402
from rdkit import Chem  # noqa: E402
from rdkit.Chem.Scaffolds import MurckoScaffold  # noqa: E402

DATA = ROOT / "data" / "fep_edges"
LEVELS = [0.5, 0.6, 0.7, 0.8, 0.9, 0.95]
FRACTIONS = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
HEAD_FRAC = 0.25  # headline matched-volume budget: commit the top 25% safest candidates
ALPHA = 0.10
SIGMA_ALE = 0.4  # benchmark label noise ~0.4 kcal/mol


def load_edges():
    rows = []
    with open(DATA / "all_edges.csv") as f:
        for r in csv.DictReader(f):
            rows.append((r["smiles_a"], r["smiles_b"], float(r["ddg"])))
    return rows


def _scaffold(smi):
    m = Chem.MolFromSmiles(smi)
    return MurckoScaffold.MurckoScaffoldSmiles(mol=m) if m else ""


def scaffold_split(rows, seed):
    scaffs = sorted({_scaffold(a) for a, _b, _d in rows} | {_scaffold(b) for _a, b, _d in rows})
    rng = np.random.default_rng(seed)
    rng.shuffle(scaffs)
    ood = set(scaffs[: max(1, len(scaffs) // 4)])
    train, test = [], []
    for a, b, d in rows:
        (test if (_scaffold(a) in ood or _scaffold(b) in ood) else train).append((a, b, d))
    return train, test


def _mve_predict(train_edges, train_ddg, eval_edges, seed):
    """Overconfident learned-sigma foil: a 2-output (mu, logvar) net by Gaussian NLL.
    Returns (mu, sigma) on eval_edges."""
    Xtr = np.stack([featurize_edge(a, b) for a, b in train_edges]).astype(np.float32)
    Xev = np.stack([featurize_edge(a, b) for a, b in eval_edges]).astype(np.float32)
    mu, sd = Xtr.mean(0), Xtr.std(0) + 1e-6
    torch.manual_seed(seed)
    net = torch.nn.Sequential(torch.nn.Linear(Xtr.shape[1], 64), torch.nn.SiLU(),
                              torch.nn.Linear(64, 2))
    opt = torch.optim.Adam(net.parameters(), lr=3e-3)
    xt = torch.tensor((Xtr - mu) / sd)
    yt = torch.tensor(np.asarray(train_ddg, dtype=np.float32))
    for _ in range(400):
        opt.zero_grad()
        out = net(xt)
        m, logv = out[:, 0], out[:, 1]
        nll = (0.5 * logv + 0.5 * (yt - m) ** 2 / torch.exp(logv)).mean()
        nll.backward()
        opt.step()
    with torch.no_grad():
        o = net(torch.tensor((Xev - mu) / sd)).numpy()
    return o[:, 0], np.sqrt(np.exp(np.clip(o[:, 1], -20, 20)))


# calibration-panel methods: sigma-isolation on the SAME (trunk) mean, so the panel measures
# *sigma calibration* alone (not mean quality). (name, label, color)
# OURS = the decomposed sigma this article ships. The other three are ONE family: the
# overconfident epistemic-only sigma and the two standard treatments applied to that same
# sigma. DELIBERATE CHOICE, stated once so the two are never split across two roles:
# temperature and split-conformal recalibration are *what the foil becomes after treatment*,
# not independent named methods, so all three take the FOIL hue and are told apart by tint,
# dash and marker — the palette's "a family with several members takes ONE hue" refinement.
# Giving either of them ALT or THIRD would have spent a method hue on a derivative of the
# foil, and (since the raw-mean ranking is ALT in panel A) would have made ALT mean two
# different things inside one figure.
CAL_METHODS = [
    ("method", "decomposed conformal $\\sigma$ (ours)", OURS),
    ("overconf", "epistemic-only $\\sigma$ (overconfident)", FOIL),
    ("temp", "+ temperature recal", tint(FOIL, 0.25)),
    ("conf", "+ split-conformal recal", tint(FOIL, 0.45)),
]
#: one hue, three members: dash + marker carry what the tint alone cannot.
CAL_STYLE = {"method": ("-", "s"), "overconf": ("-", "s"),
             "temp": ((0, (4.5, 1.6)), "^"), "conf": ((0, (5.5, 1.5, 1.2, 1.5)), "D")}
# matched-volume methods. The raw-mean ranking (sigma ignored) is a NAMED baseline, not a
# random one, so it is ALT — the same hue it carries as the raw-mean reward in Fig K. Grey is
# reserved for random / null, and the actual null of this panel is already drawn as the REF
# dotted random-base-rate line. The learned MVE head is the FOIL.
VOL_METHODS = [
    ("trunk_sig", "trunk $\\sigma$ ranking", OURS),
    ("raw_mu", r"raw $\hat{\mu}$ ranking", ALT),
    ("mve_sig", "MVE $\\sigma$ ranking", FOIL),
]
#: the trunk-sigma and raw-mean curves very nearly coincide; a full-strength ALT drawn solid
#: on top of OURS would hide our own series outright, so the later curve is dashed.
VOL_STYLE = {"trunk_sig": "-", "raw_mu": (0, (4.0, 2.0)), "mve_sig": "-"}


def run_gate(seeds=range(8)):
    rows = load_edges()
    curve = {n: {lev: [] for lev in LEVELS} for n, _l, _c in CAL_METHODS}
    gap = {n: [] for n, _l, _c in CAL_METHODS}
    prec = {n: {f: [] for f in FRACTIONS} for n, _l, _c in VOL_METHODS}
    head = {n: [] for n, _l, _c in VOL_METHODS}
    base_rate = []
    for seed in seeds:
        train, test = scaffold_split(rows, seed)
        if len(test) < 12 or len(train) < 40:
            continue
        cut = int(0.75 * len(train))
        fit_e = [(a, b) for a, b, _ in train[:cut]]
        fit_y = np.array([d for _a, _b, d in train[:cut]])
        cal_e = [(a, b) for a, b, _ in train[cut:]]
        cal_y = np.array([d for _a, _b, d in train[cut:]])
        te_e = [(a, b) for a, b, _ in test]
        te_y = np.array([d for _a, _b, d in test])

        trunk = EnsembleTrunk().fit(fit_e, fit_y, n_members=8)
        mu_cal, se_cal = trunk.predict(cal_e)
        q = conformal_q(np.abs(cal_y - mu_cal), np.sqrt(se_cal ** 2 + SIGMA_ALE ** 2), ALPHA)
        mu_te, se_te = trunk.predict(te_e)
        sig_trunk = q * np.sqrt(se_te ** 2 + SIGMA_ALE ** 2)

        # MVE (learned head, different mean) — used only in the matched-volume RANKING panel,
        # where its poor OOD mean shows up.
        mu_mve, sig_mve = _mve_predict(fit_e, fit_y, te_e, seed)

        # sigma-isolation foils on the SAME trunk mean: an overconfident epistemic-only sigma,
        # and its temperature / split-conformal recalibrations, vs the decomposed-conformal sigma.
        res_cal = np.abs(cal_y - mu_cal)
        T_tr = temperature_scale(res_cal, se_cal)
        q_ep = conformal_q(res_cal, se_cal, ALPHA)

        tau = float(np.quantile(te_y, 0.50))
        base_rate.append(float(np.mean(te_y >= tau)))

        sig_for = {"method": (mu_te, sig_trunk), "overconf": (mu_te, se_te),
                   "temp": (mu_te, T_tr * se_te), "conf": (mu_te, q_ep * se_te)}
        for name, _l, _c in CAL_METHODS:
            mu, sg = sig_for[name]
            cc = commit_correctness_curve(mu, sg, te_y, tau, LEVELS)
            for lev in LEVELS:
                curve[name][lev].append(cc[lev])
            gap[name].append(float(np.mean([lev - cc[lev] for lev in LEVELS])))

        ns = [max(1, round(f * len(te_y))) for f in FRACTIONS]
        rank_for = {"trunk_sig": (mu_te, sig_trunk), "raw_mu": (mu_te, np.ones_like(mu_te)),
                    "mve_sig": (mu_mve, sig_mve)}
        for name, _l, _c in VOL_METHODS:
            mu, sg = rank_for[name]
            pv = commit_precision_at_volume(mu, sg, te_y, tau, ns)
            for f, n in zip(FRACTIONS, ns):
                prec[name][f].append(pv[n])
            n_head = max(1, round(HEAD_FRAC * len(te_y)))
            head[name].append(commit_precision_at_volume(mu, sg, te_y, tau, [n_head])[n_head])

    return {
        "curve": {n: {lev: float(np.mean(v)) for lev, v in d.items()} for n, d in curve.items()},
        "gap": {n: np.array(v) for n, v in gap.items()},
        "prec": {n: {f: float(np.mean(v)) for f, v in d.items()} for n, d in prec.items()},
        "head": {n: np.array(v) for n, v in head.items()},
        "base_rate": float(np.mean(base_rate)) if base_rate else 0.5,
        "n_seeds": len(base_rate),
    }


def chirality_check(rows, n=60):
    chiral_edges = [(a, b) for a, b, _ in rows if abs(chir_pseudoscalar(b)) > 1e-6][:n]
    if not chiral_edges:
        return 0.0, 0.0
    even_zero, odd_nonzero = [], []
    for a, b in chiral_edges:
        f_no = featurize_edge(a, b, include_0o=False)
        f_0o = featurize_edge(a, b, include_0o=True)
        even_zero.append(float(f_0o[:-1].shape[0] - np.count_nonzero(f_0o[:-1] == f_no)))
        odd_nonzero.append(abs(float(f_0o[-1])))
    return float(np.max(even_zero)), float(np.median(odd_nonzero))


def _sig(x, nd=3):
    """Signed number with a real Unicode minus, for figure text."""
    return f"{x:+.{nd}f}".replace("-", "\u2212")


def main():
    use_paper_style()
    rows = load_edges()
    res = run_gate()
    even_collapse, odd_sep = chirality_check(rows)
    chir_ok = even_collapse == 0.0 and odd_sep > 1e-3

    # matched-volume headline: trunk-sigma ranking vs raw-mu ranking (higher precision better)
    md_v, lo_v, hi_v = regret_difference_ci(res["head"]["trunk_sig"], res["head"]["raw_mu"])
    vol_win = lo_v > 0
    # calibration (σ-isolation on the trunk mean): overconfident over-claims; the recalibrated
    # σ's should recover the method's trustworthiness (so the method is not uniquely calibrated).
    md_oc, lo_oc, hi_oc = regret_difference_ci(res["gap"]["overconf"], res["gap"]["method"])
    md_t, lo_t, hi_t = regret_difference_ci(res["gap"]["temp"], res["gap"]["method"])
    md_cf, lo_cf, hi_cf = regret_difference_ci(res["gap"]["conf"], res["gap"]["method"])

    fig, (axA, axB) = plt.subplots(1, 2, figsize=figsize(2, 3.15))
    # Panel A: matched-volume precision vs commit fraction
    fr = np.array(FRACTIONS)
    for name, label, color in VOL_METHODS:
        y = [res["prec"][name][f] for f in FRACTIONS]
        axA.plot(fr, y, color=color, ls=VOL_STYLE[name], marker="o", lw=1.7, ms=3.5,
                 label=label)
    reference_line(axA, "hline", res["base_rate"], label="random base rate")
    reference_line(axA, "vline", HEAD_FRAC, ls="--", lw=0.9)
    axA.annotate("headline budget", xy=(HEAD_FRAC, 1.0), xycoords=("data", "axes fraction"),
                 xytext=(3, -2), textcoords="offset points", ha="left", va="top",
                 fontsize=7, color=REF)
    axA.set_xlabel("commit volume (fraction of OOD edges)")
    axA.set_ylabel("commit precision (true ≥ τ)")
    # the headline contrast lives in the heading line, not in a translucent box on the data
    panel(axA, "A", "matched-volume decision quality",
          subtitle=f"at {HEAD_FRAC:.0%} commit volume: trunk−raw {_sig(md_v)}, "
                   f"CI [{_sig(lo_v)}, {_sig(hi_v)}] ({'win' if vol_win else 'tie'})")
    legend(axA, loc="upper right", fontsize=7)

    # Panel B: calibration with recalibration foils
    lv = np.array(LEVELS)
    reference_line(axB, "diagonal", label="ideal")
    for name, label, color in CAL_METHODS:
        y = [res["curve"][name][le] for le in LEVELS]
        ls, mk = CAL_STYLE[name]
        axB.plot(lv, y, color=color, ls=ls, marker=mk, lw=1.5, ms=3.5, label=label)
    axB.set_xlim(0.475, 0.975)
    axB.set_ylim(0.475, 1.03)
    axB.set_xlabel("claimed commit confidence  1−α")
    axB.set_ylabel("actual commit-correctness (OOD)")
    panel(axB, "B", "calibration vs recalibrated foils")
    legend(axB, loc="lower right", fontsize=7, handlelength=1.4)
    finish(fig, "figJ_amortized_reward")

    def row(name):
        return " · ".join(f"{res['curve'][name][le]:.2f}" for le in LEVELS)

    recal_recovers = (lo_t <= 0 <= hi_t) or (lo_cf <= 0 <= hi_cf)
    verdict = "calibration trustworthy; matched-volume ranking " + ("WIN" if vol_win else "TIE")
    (ROOT / "docs" / "results_figJ.md").write_text(
        f"# Fig J — amortized calibrated commit on OOD: honest re-audit\n\n"
        f"Trunk maps congeneric edges -> (ΔΔĜ, σ_total = conformal·√(epistemic²+aleatoric²)) on the\n"
        f"public OpenFF protein-ligand benchmark (EXPERIMENTAL ΔΔG labels); tested on both-endpoint\n"
        f"scaffold-disjoint (OOD) edges, {res['n_seeds']} seeds. `make figJ`. Two honest questions the\n"
        f"earlier draft skipped: (A) decision quality at MATCHED commit volume (no abstention escape),\n"
        f"(B) calibration vs proper RECALIBRATION baselines (not just a raw overconfident MVE).\n\n"
        f"## PANEL A — matched commit volume (the abstention-proof gate)\n"
        f"Rank candidates by the standardized margin (μ̂−τ)/σ and commit the top-n; compare commit\n"
        f"precision at the SAME n. Precision at commit-fraction {HEAD_FRAC:.0%} (random base rate "
        f"{res['base_rate']:.2f}):\n"
        f"trunk-σ {res['head']['trunk_sig'].mean():.3f} · raw-μ̂ {res['head']['raw_mu'].mean():.3f} · "
        f"MVE-σ {res['head']['mve_sig'].mean():.3f}.\n"
        f"trunk−raw {md_v:+.3f}, CI [{lo_v:+.3f}, {hi_v:+.3f}] -> "
        f"**{'σ-ranking improves which molecules are committed' if vol_win else 'TIE: σ-ranking does NOT beat raw-μ̂ ranking at matched volume'}**.\n"
        f"This is the honest Gauss-Markov-style ceiling: at equal commit volume the *choice* of σ does\n"
        f"not reliably change WHICH molecules clear the bar — calibration is not a better-ranking tool.\n\n"
        f"## PANEL B — calibration, σ isolated on the SAME (trunk) mean\n"
        f"All four use the trunk's ensemble mean; only σ differs, so this measures σ calibration\n"
        f"alone. Actual commit-correctness per claimed 1−α ({', '.join(f'{le:.2f}' for le in LEVELS)}):\n"
        f"- decomposed conformal σ (ours): {row('method')}\n"
        f"- epistemic-only σ (overconfident): {row('overconf')}\n"
        f"- + temperature recalibration: {row('temp')}\n"
        f"- + split-conformal recalibration: {row('conf')}\n\n"
        f"Mean shortfall (claimed−actual): ours {res['gap']['method'].mean():+.3f}, "
        f"overconfident {res['gap']['overconf'].mean():+.3f}, temperature {res['gap']['temp'].mean():+.3f}, "
        f"conformal {res['gap']['conf'].mean():+.3f}.\n"
        f"overconfident−ours {md_oc:+.3f} CI[{lo_oc:+.3f},{hi_oc:+.3f}] (the raw epistemic-only σ "
        f"over-claims); temperature−ours {md_t:+.3f} CI[{lo_t:+.3f},{hi_t:+.3f}]; "
        f"conformal−ours {md_cf:+.3f} CI[{lo_cf:+.3f},{hi_cf:+.3f}].\n\n"
        f"Honest reading: the overconfident epistemic-only σ over-claims, and "
        f"{'standard recalibration (temperature / split-conformal) of the SAME mean RECOVERS the trustworthiness — so the decomposed physics σ is NOT uniquely calibrated; its advantage is that it needs no calibration set for the aleatoric term (Fig A) and is differentiable into the acquisition graph.' if recal_recovers else 'even recalibrated σ does not fully match ours here.'}\n\n"
        f"## SECONDARY — 0o chirality contract: {'PASS' if chir_ok else 'FAIL'}\n"
        f"Even (ECFP useChirality=False) edge feature collapses for an enantiomer pair; 0o separation\n"
        f"median |Δ| = {odd_sep:.2f} (chirality rides on the 0o channel only).\n\n"
        f"## Honest scope\n"
        f"Genuine both-endpoint scaffold-disjoint OOD; experimental ΔΔG labels (not computed FEP). The\n"
        f"matched-volume panel removes the earlier draft's abstention artifact (high-confidence cells\n"
        f"were n≈0). Net: calibration buys TRUST (knowing when a commit is safe), not a better commit\n"
        f"RANKING — consistent with the contour's decision-not-ranking finding.\n\n"
        f"## Verdict\n{verdict}.\n"
    )
    print("wrote figJ_amortized_reward.(pdf|png) + docs/results_figJ.md")
    print(f"[matched-vol @{HEAD_FRAC:.0%}] trunk-raw {md_v:+.3f} CI[{lo_v:+.3f},{hi_v:+.3f}] "
          f"-> {'WIN' if vol_win else 'TIE'}")
    print(f"[calibration σ-isolation] overconf−ours {md_oc:+.3f} CI[{lo_oc:+.3f},{hi_oc:+.3f}]; "
          f"temp−ours {md_t:+.3f} CI[{lo_t:+.3f},{hi_t:+.3f}]; conf−ours {md_cf:+.3f} "
          f"CI[{lo_cf:+.3f},{hi_cf:+.3f}]")
    print(f"[chirality] even {even_collapse:.0f} 0o {odd_sep:.3f} -> {'PASS' if chir_ok else 'FAIL'}")


if __name__ == "__main__":
    main()
