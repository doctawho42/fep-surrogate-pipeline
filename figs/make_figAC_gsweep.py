"""SI Fig (autocorrelation g-sweep) — the sandwich calibration is robust to real MD
autocorrelation *provided the works are decorrelated to the statistical inefficiency g*.

Referee ask #1 [MUST] (docs/anticipated_referee_responses.md): report the sandwich se at raw n
vs n_eff = n/g and the panel-B calibration ratio under both. NO new MD — this reuses the same
alchemtest BACE1 RBFE complex-leg works as Fig A panel B.

Why this is the one item a fresh FEP referee insists on: the sandwich B = n_f Var_f[p] +
n_r Var_r[p] (Theorem 2) assumes *independent* samples, but raw MD works are autocorrelated. The
estimator docstring is explicit that feeding raw correlated samples under-estimates the variance.
This panel quantifies exactly that on real data:

  * RAW n (naive, treat correlated samples as iid): the sandwich se is OVERCONFIDENT vs the
    correlation-aware (decorrelated-bootstrap) truth, by ~sqrt(g_bar).
  * n_eff = n/g (the paper's protocol, stated at paper_body.tex:172): the sandwich se MATCHES the
    truth (ratio ~ 1), reproducing the calibrated Fig A panel-B result.
  * The sandwich also tracks its OWN matched bootstrap under both n (ratio ~ 1 either way): the
    formula is internally consistent — the only thing that must be right is the sample count n,
    which the paper decorrelates to g throughout. So the calibration VERDICT is invariant to
    autocorrelation once n_eff is used.

Run:  python figs/make_figAC_gsweep.py   (or `make gsweep`)
"""
from __future__ import annotations

import pathlib
import sys
import warnings

import matplotlib.pyplot as plt
import numpy as np

warnings.filterwarnings("ignore")
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "src"))
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from bar.estimator import sandwich_variance, solve_bar  # noqa: E402
from paperstyle import (  # noqa: E402
    OURS, REF, THIRD, figsize, finish, legend, panel, reference_line, tint, use_paper_style,
)

FIGDIR = pathlib.Path(__file__).resolve().parent
RNG_SEED = 20260710

# Semantic colours (paperstyle). The n_eff sandwich is the paper's own protocol -> OURS. The
# raw-n sandwich is the naive practice this section corrects, exactly as naive 1/I is in Fig A,
# so it takes THIRD; its self-consistency check against its own raw bootstrap is the same
# (naive) quantity family and takes a tint of it. Grey is reserved for reference lines.
C_RAW = THIRD
C_NEFF = OURS
C_RAW_SELF = tint(THIRD, 0.45)


def _edges_from_unk(frames):
    """frames: list of (lambda, u_nk DataFrame). Yields (w_f, w_r) per adjacent lambda pair."""
    frames = sorted(frames, key=lambda t: t[0])
    cols = sorted(frames[0][1].columns)
    for i in range(len(frames) - 1):
        ci, cj = cols[i], cols[i + 1]
        u_i, u_j = frames[i][1], frames[i + 1][1]
        yield (u_i[cj] - u_i[ci]).to_numpy(), (u_j[ci] - u_j[cj]).to_numpy()


def _bace_edges():
    """Real protein-ligand BINDING FEP edges (BACE1 RBFE, alchemtest AMBER, complex legs)."""
    from alchemlyb.parsing.amber import extract_u_nk
    from alchemtest.amber import load_bace_example
    data = load_bace_example()["data"]["complex"]
    for files in data.values():
        frames = []
        for f in files:
            u = extract_u_nk(f, T=298.0)
            frames.append((float(u.index.get_level_values(-1)[0]), u))
        yield from _edges_from_unk(frames)


def _boot_se(wf, wr, rng, n_boot=400):
    """iid bootstrap SD of the BAR estimate on the given (already-chosen) work arrays."""
    d = np.array([solve_bar(wf[rng.integers(0, wf.size, wf.size)],
                            -wr[rng.integers(0, wr.size, wr.size)]) for _ in range(n_boot)])
    return float(d.std(ddof=1))


def collect(seed=RNG_SEED, n_boot=400):
    from pymbar.timeseries import statistical_inefficiency, subsample_correlated_data
    rng = np.random.default_rng(seed)
    rows = []
    for w_f, w_r in _bace_edges():
        # statistical inefficiency per leg
        try:
            gf, gr = statistical_inefficiency(w_f), statistical_inefficiency(w_r)
        except Exception:
            gf = gr = 1.0
        idx_f = subsample_correlated_data(w_f, g=gf)
        idx_r = subsample_correlated_data(w_r, g=gr)
        wf_s, wr_s = w_f[idx_f], w_r[idx_r]
        if wf_s.size < 8 or wr_s.size < 8 or w_f.size < 8 or w_r.size < 8:
            continue
        # sandwich se at raw n and at n_eff = n/g
        se_raw = float(np.sqrt(max(sandwich_variance(w_f, -w_r), 0.0)))
        se_neff = float(np.sqrt(max(sandwich_variance(wf_s, -wr_s), 0.0)))
        # correlation-aware truth = decorrelated (n_eff) iid bootstrap SD (same truth as Fig A-B)
        truth_neff = _boot_se(wf_s, wr_s, rng, n_boot)
        # matched (wrong) truth = raw iid bootstrap SD (treats correlated points as iid)
        truth_raw = _boot_se(w_f, w_r, rng, n_boot)
        if truth_neff <= 0 or truth_raw <= 0:
            continue
        rows.append(dict(
            g=float(0.5 * (gf + gr)), n_raw=int(w_f.size), n_eff=int(wf_s.size),
            se_raw=se_raw, se_neff=se_neff, truth_neff=truth_neff, truth_raw=truth_raw,
            r_raw=se_raw / truth_neff,      # raw sandwich vs correlation-aware truth (overconfident)
            r_neff=se_neff / truth_neff,    # n_eff sandwich vs truth (calibrated, the paper's number)
            r_raw_matched=se_raw / truth_raw,   # sandwich vs its own (naive) bootstrap (self-consistent)
        ))
    return rows


def main() -> None:
    use_paper_style()
    rows = collect()
    n = len(rows)
    g = np.array([d["g"] for d in rows])
    r_raw = np.array([d["r_raw"] for d in rows])
    r_neff = np.array([d["r_neff"] for d in rows])
    r_raw_matched = np.array([d["r_raw_matched"] for d in rows])
    se_raw = np.array([d["se_raw"] for d in rows])
    se_neff = np.array([d["se_neff"] for d in rows])
    truth = np.array([d["truth_neff"] for d in rows])
    gbar = float(np.exp(np.mean(np.log(g))))  # geometric mean g

    fig, (axA, axB) = plt.subplots(1, 2, figsize=figsize(2, 3.1))

    # Panel A: sandwich se vs correlation-aware (decorrelated bootstrap) truth, raw vs n_eff
    lim = [0.9 * min(se_raw.min(), se_neff.min(), truth.min()),
           1.1 * max(se_raw.max(), se_neff.max(), truth.max())]
    axA.plot(lim, lim, ls=":", color=REF, lw=1.0, zorder=0)
    axA.scatter(truth, se_raw, s=26, facecolors="none", edgecolors=C_RAW, marker="o",
                lw=1.1, zorder=3, label=f"raw $n$ (naive): {r_raw.mean():.2f}$\\times$ truth")
    axA.scatter(truth, se_neff, s=26, c=C_NEFF, marker="D", zorder=4,
                label=f"$n_{{\\rm eff}}=n/g$ (used): {r_neff.mean():.2f}$\\times$ truth")
    axA.set_xscale("log"); axA.set_yscale("log")
    axA.set_xlabel("correlation-aware truth se (decorrelated bootstrap)")
    axA.set_ylabel("sandwich $(B/I^2)^{1/2}$ se")
    legend(axA, loc="upper left")
    panel(axA, "A", "raw $n$ overconfident, $n_{\\rm eff}$ calibrated")

    # Panel B: aggregate calibration ratios (mean +/- bootstrap CI over edges)
    def _mean_ci(x, nb=2000):
        rr = np.random.default_rng(0)
        bs = np.array([x[rr.integers(0, x.size, x.size)].mean() for _ in range(nb)])
        return x.mean(), np.percentile(bs, 2.5), np.percentile(bs, 97.5)

    labels = ["sandwich(raw)\n/ truth", "sandwich($n_{\\rm eff}$)\n/ truth",
              "sandwich(raw)\n/ raw-boot"]
    series = [r_raw, r_neff, r_raw_matched]
    colors = [C_RAW, C_NEFF, C_RAW_SELF]
    xs = np.arange(len(series))
    bounds = []
    for x0, xr_, c in zip(xs, series, colors):
        m, lo, hi = _mean_ci(xr_)
        axB.errorbar([x0], [m], yerr=[[m - lo], [hi - m]], fmt="o", color=c, ms=6, capsize=4)
        bounds += [lo, hi]
    g_line = 1.0 / np.sqrt(gbar)
    # two reference lines, so both are REF grey and are told apart by dash pattern, not by hue.
    reference_line(axB, "hline", 1.0, label="calibrated (=1)")
    axB.axhline(g_line, color=REF, ls=(0, (4.5, 2.0)), lw=1.0, zorder=0.5,
                label=f"$\\bar g^{{-1/2}}={g_line:.2f}$ ($\\bar g={gbar:.1f}$)")
    axB.set_xticks(xs); axB.set_xticklabels(labels, fontsize=8)
    axB.set_xlim(-0.6, len(series) - 0.4)
    axB.set_ylabel("reported se / truth se")
    # the ratios and their intervals live in a band around 1; a 0-based axis spent three
    # quarters of its height on empty space. The legend goes in the band below the lowest
    # interval, which no marker and no error bar reaches.
    lo_all, hi_all = min(bounds + [g_line]), max(bounds + [1.0])
    axB.set_ylim(lo_all - 0.13, hi_all + 0.04)
    axB.tick_params(axis="x", length=0)
    legend(axB, loc="lower right", fontsize=7)
    panel(axB, "B", "the fix is $n_{\\rm eff}$, not the formula")

    finish(fig, "figAC_gsweep")

    # --- machine-readable summary for the SI text ---
    print(f"\n[g-sweep] BACE1 binding edges: n={n}")
    print(f"  statistical inefficiency g: geo-mean {gbar:.2f}, "
          f"range [{g.min():.2f},{g.max():.2f}]; median n_raw={int(np.median([d['n_raw'] for d in rows]))}, "
          f"median n_eff={int(np.median([d['n_eff'] for d in rows]))}")
    print(f"  sandwich(n_eff)/truth  = {r_neff.mean():.2f}  [{r_neff.min():.2f},{r_neff.max():.2f}]   (calibrated; the paper's protocol)")
    print(f"  sandwich(raw)  /truth  = {r_raw.mean():.2f}  [{r_raw.min():.2f},{r_raw.max():.2f}]   (overconfident ~ 1/sqrt(g_bar)={1/np.sqrt(gbar):.2f})")
    print(f"  sandwich(raw)/raw-boot = {r_raw_matched.mean():.2f}  [{r_raw_matched.min():.2f},{r_raw_matched.max():.2f}]   (formula self-consistent under raw n too)")


if __name__ == "__main__":
    main()
