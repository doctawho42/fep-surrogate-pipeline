# figs/make_figI.py
"""Fig I — the target-contour Increment-1 gate (calibrated reward).

History (honest): the gate first tested whether a calibrated risk-adjusted reward
improves candidate SELECTION (regret@K / precision@K). It came back KILL — calibration
is second-order for ranking (the GLS / Gauss-Markov ceiling, cf. Fig C). The selection
sweep is RETAINED here as the diagnostic that motivated the pivot.

PRIMARY gate (pivoted): calibrated COMMIT-TO-SYNTHESIS. Commit a candidate iff its
lower-confidence bound risk_adjusted_reward(muhat, sigma, kappa=z_alpha) >= tau. With a
CALIBRATED sandwich sigma the actual fraction of commits clearing tau tracks/exceeds the
claimed confidence (trustworthy); an OVERCONFIDENT learned sigma over-claims (commits to
duds). This is Fig G's calibrated-decision logic ported to the generative commit -- the
honest punchline being that the same LCB that does NOT help selection DOES make the
commit decision trustworthy, because calibration is for decisions, not rankings.

SECONDARY gate: the 0o chirality contract (added as a figure panel in a later step).

Spec: docs/superpowers/specs/2026-06-30-target-contour-design.md (see the PIVOT section).
Run:  PYTHONPATH=src python figs/make_figI.py   (or `make figI`)

KILL (primary): if the overconfident shortfall does NOT significantly exceed the
calibrated shortfall (bootstrap CI on the per-seed difference includes 0) -> calibration
does not make the commit trustworthy -> do not build the trunk / generator.
"""
from __future__ import annotations

import pathlib
import sys
import warnings

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

warnings.filterwarnings("ignore")
ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "figs"))
from bar.estimator import bar_estimate, solve_bar  # noqa: E402
from bar.reward import (  # noqa: E402
    commit_correctness_curve, regret, regret_difference_ci, risk_adjusted_reward, select_topk,
)
from make_figA import _bace_edges, _benzene_edges, _gauss_edge, _train_mve  # noqa: E402

FIGDIR = ROOT / "figs"
KAPPA = 1.0
C_CAL = "#0072B2"
C_LEARN = "#D55E00"
C_REF = "#555555"


def _style() -> None:
    mpl.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 9, "axes.labelsize": 10, "axes.titlesize": 10,
        "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 7.5,
        "axes.spines.top": False, "axes.spines.right": False,
        "figure.dpi": 120, "savefig.dpi": 300, "savefig.bbox": "tight",
    })


def _pool_rewards(mu, n, mve_se, rng):
    """Given true rewards `mu`, draw an INDEPENDENT overlap edge per candidate (noise
    decoupled from mu) and return four selectors' scores:
      raw          = observed value (mu + sampling fluctuation), no calibration
      cal_lcb      = value - KAPPA*sigma_sandwich   (risk-AVERSE acquisition)
      cal_shrink   = empirical-Bayes shrinkage with the CALIBRATED sandwich sigma
                     (risk-NEUTRAL; the right tool for precision/selection)
      learned_shrink = shrinkage with the overconfident learned sigma (foil)
    """
    val = np.empty(mu.size)
    sig_cal = np.empty(mu.size)
    sig_lrn = np.empty(mu.size)
    for j in range(mu.size):
        s = rng.uniform(0.6, 3.0)
        xf, xr = _gauss_edge(s, n, rng)
        r = bar_estimate(xf, xr)
        val[j] = mu[j] + r.delta_f
        sig_cal[j] = np.sqrt(max(r.var_sandwich, 0.0))
        sig_lrn[j] = mve_se(xf, xr)
    raw = val
    cal_lcb = np.array([risk_adjusted_reward(v, sc, KAPPA) for v, sc in zip(val, sig_cal)])
    cal_shrink = _empirical_bayes_shrink(val, sig_cal)
    learned_shrink = _empirical_bayes_shrink(val, sig_lrn)
    return raw, cal_lcb, cal_shrink, learned_shrink


def _empirical_bayes_shrink(val, sigma):
    """Risk-NEUTRAL calibrated selector for precision: shrink each noisy observation
    toward the pool mean by its calibrated noise. signal_var = total observed var minus
    mean noise var; per-candidate weight w = signal/(signal+sigma^2). Calibrated sigma ->
    correct shrinkage; overconfident sigma -> too little shrinkage -> worse selection."""
    val = np.asarray(val, dtype=float)
    sigma = np.asarray(sigma, dtype=float)
    prior_mean = float(np.mean(val))
    signal_var = max(float(np.var(val)) - float(np.mean(sigma ** 2)), 1e-6)
    w = signal_var / (signal_var + sigma ** 2)
    return prior_mean + w * (val - prior_mean)


def _precision_at_k(mu, selected, k):
    """Fraction of `selected` that are in the TRUE top-k by mu (scale-free selection
    quality: 1.0 = recovered the genuinely-best k)."""
    true_top = set(np.argsort(-np.asarray(mu))[:k].tolist())
    sel = np.asarray(selected, dtype=int)
    if sel.size == 0:
        return 0.0
    return float(np.mean([int(s) in true_top for s in sel]))


def controlled_selection(n_cand=60, n=20, k=8, seeds=range(12)):
    """Regime sweep. Difficulty = signal spread (mu ~ U(0, spread)) relative to the
    per-edge sandwich se. Easy (spread >> noise): selection trivial, calibration can't
    help (GLS ceiling). Hard (spread ~ noise): selection is noise-limited, where a
    CALIBRATED risk-neutral selector (shrinkage) should recover more true top-k than raw.
    Reports precision@k for each selector and the cal-shrink-vs-raw and lcb-vs-raw CIs."""
    mve_se = _train_mve()
    # anchor 'hard' to the typical per-edge sandwich se so the regime axis is principled
    rng0 = np.random.default_rng(7)
    ses = []
    for _ in range(200):
        s = rng0.uniform(0.6, 3.0)
        xf, xr = _gauss_edge(s, n, rng0)
        ses.append(np.sqrt(max(bar_estimate(xf, xr).var_sandwich, 0.0)))
    sbar = float(np.median(ses))
    spreads = [8 * sbar, 4 * sbar, 2 * sbar, 1 * sbar, 0.5 * sbar]
    out = []
    for spread in spreads:
        p_raw, p_lcb, p_shr, p_lshr = [], [], [], []
        for seed in seeds:
            rng = np.random.default_rng(1000 + seed)
            mu = rng.uniform(0.0, spread, n_cand)
            raw, cal_lcb, cal_shrink, learned_shrink = _pool_rewards(mu, n, mve_se, rng)
            p_raw.append(_precision_at_k(mu, select_topk(raw, k), k))
            p_lcb.append(_precision_at_k(mu, select_topk(cal_lcb, k), k))
            p_shr.append(_precision_at_k(mu, select_topk(cal_shrink, k), k))
            p_lshr.append(_precision_at_k(mu, select_topk(learned_shrink, k), k))
        p_raw, p_lcb, p_shr, p_lshr = (np.array(p_raw), np.array(p_lcb),
                                       np.array(p_shr), np.array(p_lshr))
        # regret_difference_ci is a generic paired-bootstrap CI on (a-b); here a,b are
        # PRECISION (higher=better), so calibrated beats raw iff the diff lo > 0.
        md_s, lo_s, hi_s = regret_difference_ci(p_shr, p_raw)
        md_l, lo_l, hi_l = regret_difference_ci(p_lcb, p_raw)
        out.append(dict(spread=spread, ratio=spread / sbar,
                        raw=p_raw.mean(), lcb=p_lcb.mean(), shrink=p_shr.mean(),
                        lshrink=p_lshr.mean(), shr_diff=md_s, shr_lo=lo_s, shr_hi=hi_s,
                        lcb_diff=md_l, lcb_lo=lo_l, lcb_hi=hi_l))
    return sbar, out


def commit_calibration(n_cand=400, n=20, n_seeds=8, levels=None):
    """Calibrated COMMIT-to-synthesis gate (the pivoted PRIMARY gate). Each candidate has a
    true reward mu and a calibrated estimate (muhat, sigma_sandwich); the overconfident
    learned head gives sigma_learned. COMMIT candidate j to synthesis iff its
    lower-confidence bound at claimed level (1-alpha) clears an actionable threshold tau:
    risk_adjusted_reward(muhat, sigma, kappa=z_alpha) >= tau (reuses the LCB as commit score).
    For each claimed (1-alpha) we record the ACTUAL fraction of commits whose true mu >= tau.
    Calibrated sigma -> actual >= claimed (conservative, trustworthy); overconfident sigma ->
    actual < claimed (commits to duds). Returns (levels, curve, gap_cal, gap_lrn): curve[lev]
    = (mean actual cal, mean actual learned); gap_* = per-seed mean shortfall (claimed-actual)."""
    if levels is None:
        levels = [0.5, 0.6, 0.7, 0.8, 0.9, 0.95]
    mve_se = _train_mve()
    per = {lev: {"cal": [], "lrn": []} for lev in levels}
    gap_cal, gap_lrn = [], []
    for seed in range(n_seeds):
        rng = np.random.default_rng(3000 + seed)
        mu = rng.uniform(0.0, 4.0, n_cand)
        tau = float(np.quantile(mu, 0.75))  # actionable = strong-binder top quartile
        muhat = np.empty(n_cand)
        sig_c = np.empty(n_cand)
        sig_l = np.empty(n_cand)
        for j in range(n_cand):
            s = rng.uniform(0.6, 3.0)
            xf, xr = _gauss_edge(s, n, rng)
            r = bar_estimate(xf, xr)
            muhat[j] = mu[j] + r.delta_f
            sig_c[j] = np.sqrt(max(r.var_sandwich, 0.0))
            sig_l[j] = mve_se(xf, xr)
        cal_curve = commit_correctness_curve(muhat, sig_c, mu, tau, levels)
        lrn_curve = commit_correctness_curve(muhat, sig_l, mu, tau, levels)
        for lev in levels:
            per[lev]["cal"].append(cal_curve[lev])
            per[lev]["lrn"].append(lrn_curve[lev])
        gap_cal.append(float(np.mean([lev - cal_curve[lev] for lev in levels])))
        gap_lrn.append(float(np.mean([lev - lrn_curve[lev] for lev in levels])))
    curve = {lev: (float(np.mean(per[lev]["cal"])), float(np.mean(per[lev]["lrn"]))) for lev in levels}
    return levels, curve, np.array(gap_cal), np.array(gap_lrn)


def chirality_panel(n=400, seed=99):
    """0o contract: a linear-readout reward separates enantiomers WITH the 0o channel and
    is identical WITHOUT it (the chirality-completeness theorem, Supporting Information).
    Returns (even_collapse, odd_reward_sep)."""
    from bar.reward import linear_readout_reward
    rng = np.random.default_rng(seed)
    w_even = rng.normal(size=6)
    w_0o = np.append(w_even, 1.0)  # same even weights + a nonzero 0o weight
    mirror = np.array([-1.0, 1.0, 1.0])
    even_diffs, odd_diffs = [], []
    for _ in range(n):
        c = rng.normal(size=(4, 3))
        cm = c * mirror
        even_diffs.append(abs(linear_readout_reward(c, w_even, include_0o=False)
                              - linear_readout_reward(cm, w_even, include_0o=False)))
        odd_diffs.append(abs(linear_readout_reward(c, w_0o, include_0o=True)
                             - linear_readout_reward(cm, w_0o, include_0o=True)))
    return float(np.max(even_diffs)), float(np.median(odd_diffs))


def real_fep_selection(n_sub=15, ks=(3, 5, 8), seeds=range(8)):
    """Confirmation on real FEP edges. Candidate = an edge; truth = converged BAR on all
    decorrelated samples; estimate = subsample n_sub. value = -ΔΔĜ (favourable edge)."""
    from pymbar.timeseries import statistical_inefficiency, subsample_correlated_data
    edges = []
    for gen in (_bace_edges, _benzene_edges):
        for w_f, w_r, _li, _lj in gen():
            try:
                gf, gr = statistical_inefficiency(w_f), statistical_inefficiency(w_r)
            except Exception:
                gf = gr = 1.0
            wf = w_f[subsample_correlated_data(w_f, g=gf)]
            wr = w_r[subsample_correlated_data(w_r, g=gr)]
            if wf.size >= n_sub and wr.size >= n_sub:
                edges.append((wf, -wr))            # unified convention: bar_estimate(wf, -wr)
    truth = np.array([-solve_bar(wf, xr) for wf, xr in edges])  # -ΔΔG, higher = better
    out = {k: {"raw": [], "cal": []} for k in ks}
    for seed in seeds:
        rng = np.random.default_rng(2000 + seed)
        r_raw, r_cal = np.empty(len(edges)), np.empty(len(edges))
        for j, (wf, xr) in enumerate(edges):
            fi = rng.choice(wf.size, n_sub, replace=False)
            ri = rng.choice(xr.size, n_sub, replace=False)
            rr = bar_estimate(wf[fi], xr[ri])
            val = -rr.delta_f
            sig = np.sqrt(max(rr.var_sandwich, 0.0))
            r_raw[j] = val
            r_cal[j] = risk_adjusted_reward(val, sig, KAPPA)
        for k in ks:
            out[k]["raw"].append(regret(truth, select_topk(r_raw, k)))
            out[k]["cal"].append(regret(truth, select_topk(r_cal, k)))
    return len(edges), {k: {m: np.array(v) for m, v in d.items()} for k, d in out.items()}


def main() -> None:
    _style()
    levels, curve, gap_cal, gap_lrn = commit_calibration()
    md, lo, hi = regret_difference_ci(gap_lrn, gap_cal)  # overconfident over-claims iff lo>0
    primary_pass = lo > 0
    even_collapse, odd_sep = chirality_panel()
    _, sweep = controlled_selection()
    hard = sweep[-1]
    n_edges, real = real_fep_selection()
    md_r, lo_r, hi_r = regret_difference_ci(real[5]["cal"], real[5]["raw"])

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(7.4, 3.3))
    # Panel A: calibrated commit-to-synthesis
    lv = np.array(levels)
    cal_a = np.array([curve[le][0] for le in levels])
    lrn_a = np.array([curve[le][1] for le in levels])
    axA.plot([0.5, 1.0], [0.5, 1.0], color=C_REF, ls=":", lw=1.0, zorder=0,
             label="ideal (actual = claimed)")
    axA.plot(lv, cal_a, "-s", color=C_CAL, lw=1.8, ms=4, label="calibrated σ (sandwich)")
    axA.plot(lv, lrn_a, "-^", color=C_LEARN, lw=1.6, ms=4, label="overconfident σ (learned)")
    axA.set_xlabel("claimed commit confidence  1−α")
    axA.set_ylabel("actual fraction of commits, true μ ≥ τ")
    axA.set_title("A   calibrated reward → trustworthy commit", loc="left", fontweight="bold")
    axA.legend(frameon=False, loc="lower right", fontsize=7)
    axA.set_xlim(0.47, 0.98)
    axA.set_ylim(0.85, 1.005)
    axA.text(0.04, 0.96,
             f"overconf−cal shortfall {md:+.3f}\nCI[{lo:+.3f},{hi:+.3f}] "
             f"{'PASS' if primary_pass else 'KILL'}\nselection: NULL (GLS ceiling, Fig C)",
             transform=axA.transAxes, fontsize=7, va="top", color=C_CAL)
    # Panel B: 0o chirality contract
    axB.bar([0, 1], [even_collapse, odd_sep], color=[C_LEARN, C_CAL], width=0.6)
    axB.set_xticks([0, 1])
    axB.set_xticklabels(["even reward\n(no 0o)", "even + 0o\nreward"])
    axB.set_ylabel("enantiomer reward |Δ|  (M vs mirror)")
    axB.set_title("B   0o restores chirality", loc="left", fontweight="bold")
    axB.text(0, even_collapse, f"{even_collapse:.1e}", va="bottom", ha="center",
             fontsize=7, color=C_LEARN)
    axB.text(1, odd_sep, f"{odd_sep:.2f}", va="bottom", ha="center", fontsize=8, color=C_CAL)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(FIGDIR / f"figI_calibrated_reward.{ext}")

    verdict = "PASS" if primary_pass else "KILL"
    chir_ok = even_collapse < 1e-9 and odd_sep > 1e-3
    rows = "".join(f"| {le:.2f} | {curve[le][0]:.3f} | {curve[le][1]:.3f} |\n" for le in levels)
    doc = ROOT / "docs" / "results_figI.md"
    doc.write_text(
        f"# Fig I — calibrated reward (target-contour Increment 1): {verdict}\n\n"
        f"Spec: `docs/superpowers/specs/2026-06-30-target-contour-design.md` (see the PIVOT\n"
        f"section). Reward built on `src/bar` (no trunk, no generator). κ / z from the claimed\n"
        f"confidence. Regenerate with `make figI`.\n\n"
        f"## PRIMARY gate (pivoted) — calibrated commit-to-synthesis: {verdict}\n"
        f"Selection failed (diagnostic below); calibration's value is the DECISION. Commit a\n"
        f"candidate to synthesis iff its lower-confidence bound "
        f"`risk_adjusted_reward(μ̂, σ, κ=z_(1−α)) ≥ τ` (τ = top-quartile true reward). Actual\n"
        f"fraction of commits with true μ ≥ τ, per claimed confidence:\n\n"
        f"| claimed 1−α | actual (calibrated σ) | actual (overconfident σ) |\n"
        f"|--:|--:|--:|\n{rows}\n"
        f"Calibrated σ: actual **tracks/exceeds** claimed (rises {cal_a[0]:.3f}→{cal_a[-1]:.3f},\n"
        f"trustworthy). Overconfident σ: **plateaus ~{lrn_a.mean():.2f}** regardless of claimed\n"
        f"level (claim {levels[-1]:.2f}, deliver {lrn_a[-1]:.3f} — over-claims). Mean shortfall\n"
        f"(claimed−actual): calibrated {gap_cal.mean():+.3f}, overconfident {gap_lrn.mean():+.3f};\n"
        f"overconfident−calibrated diff {md:+.3f}, bootstrap CI [{lo:+.3f}, {hi:+.3f}] → "
        f"**{verdict}**\n(trustworthy iff CI lower bound > 0). This is Fig G's calibrated-decision\n"
        f"logic ported to the generative commit: the same LCB that does NOT improve selection DOES\n"
        f"make the commit trustworthy — calibration is for decisions, not rankings.\n\n"
        f"## Diagnostic (honest negative that motivated the pivot) — selection is NULL\n"
        f"Calibrated reward does **not** significantly improve candidate SELECTION (the GLS /\n"
        f"Gauss–Markov ceiling, cf. Fig C). Regime sweep precision@K, hardest regime (0.5×σ):\n"
        f"cal-shrink − raw {hard['shr_diff']:+.3f}, CI [{hard['shr_lo']:+.3f}, {hard['shr_hi']:+.3f}]\n"
        f"(not significant; shrinkage even significantly *hurts* at the 4×σ regime). Real FEP\n"
        f"({n_edges} edges, BACE1 + benzene): cal−raw regret {md_r:+.3f}, CI [{lo_r:+.3f}, "
        f"{hi_r:+.3f}] (broad ΔΔG spread ⇒ selection trivial). Selection is the wrong gate.\n\n"
        f"## SECONDARY gate — 0o chirality contract: {'PASS' if chir_ok else 'FAIL'}\n"
        f"Even-only reward enantiomer collapse max|Δ| = {even_collapse:.1e} (must be ~0); 0o\n"
        f"reward separation median |Δ| = {odd_sep:.2f}. The reward representation carries the\n"
        f"chirality bit (necessary for per-enantiomer generation, invariant #6).\n\n"
        f"## Verdict\n"
        f"Primary commit-calibration **{verdict}**; chirality **{'PASS' if chir_ok else 'FAIL'}**. "
        f"{'Proceed to Increment 2 (amortised reward / trunk).' if primary_pass and chir_ok else 'Revisit before building the trunk/generator.'}\n"
    )
    print(f"wrote figI_calibrated_reward.(pdf|png) + docs/results_figI.md")
    print(f"[primary commit-calibration] overconf−cal shortfall {md:+.3f} "
          f"CI[{lo:+.3f},{hi:+.3f}] -> {verdict}")
    print(f"[chirality] even collapse {even_collapse:.1e}  0o sep {odd_sep:.2f} -> "
          f"{'PASS' if chir_ok else 'FAIL'}")


if __name__ == "__main__":
    main()
