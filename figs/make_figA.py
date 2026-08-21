"""Fig A — "target the sandwich" (reframed: correctness + the right foil).

The BAR-bottleneck reads the aleatoric (sampling) variance off the estimator as the
SANDWICH B/I^2. The honest message is NOT "we beat naive 1/I" (no competent FEP user
reports bare 1/I — pymbar already subtracts the nrat correction). It is:

  * CORRECTNESS — the sandwich COINCIDES with pymbar's MBAR uncertainty and with the
    Monte-Carlo truth across every overlap regime, with NO Monte-Carlo needed.
  * THE RIGHT FOIL — a learned-variance head (MVE / heteroscedastic NN, what an ML
    practitioner actually builds) trained on a realistic label budget is badly
    OVERCONFIDENT (se ~ 5-10x too small), because one noisy ΔΔĜ per edge cannot teach
    the per-edge sampling variance. The BAR bottleneck COMPUTES it for free, untrained.
  * The naive `1/I` (information-equality plug-in) is shown only as the textbook value
    it corrects, off by a *varying* factor — NOT as a baseline anyone reports.

Why it matters downstream: this calibrated aleatoric variance is a differentiable
closed form (O(1) backward) that propagates into the surrogate and the Laplacian edge
weights I_e^2/B_e — pymbar gives a number, not a backprop-able graph weight. And B/I^2
is robust (>= 0 always; pymbar's 1/I - nrat can go negative -> nan in a rare edge
case: very high overlap / tiny n).

Panel A: controlled MC-truth sweep. Panel B: real FEP edges — protein-ligand BINDING
(BACE1 RBFE, alchemtest AMBER) plus benzene solvation, sandwich vs bootstrap truth.

Run:  python figs/make_figA.py   (or `make figA`)
"""
from __future__ import annotations

import pathlib
import sys
import warnings

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedLocator, LogLocator, NullFormatter

warnings.filterwarnings("ignore")
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "src"))
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from bar.estimator import bar_estimate, mbar_variance, sandwich_variance, solve_bar  # noqa: E402
from paperstyle import (  # noqa: E402
    ALT,
    FOIL,
    INK,
    MUTED,
    OURS,
    REF,
    THIRD,
    figsize,
    finish,
    legend,
    panel,
    reference_line,
    tint,
    use_paper_style,
)

FIGDIR = pathlib.Path(__file__).resolve().parent
RNG_SEED = 20260629


def _norm_overlap(I, n_f, n_r):
    return 4.0 * I / (n_f + n_r)


def _gauss_edge(s, n, rng):
    mf, mr = s**2 / 2, -(s**2) / 2
    return rng.normal(mf, s, n), rng.normal(mr, s, n)


def pooled_se_recovery(n_edges, n=20, nbins=12, seed=7):
    """Nonparametric identifiability check (reviewer round-2, §3): from n_edges edges with ONE
    ΔF-hat label each, can the conditional sampling se(overlap) be recovered? The _gauss_edge model
    has true ΔF = 0 for every s, so within a narrow overlap bin the residual SD of the single-label
    ΔF-hat IS the pure sampling SD. We pool by overlap and compare to the mean sandwich se per bin.
    Returns the mean relative error; it shrinks with n_edges -> the conditional variance is
    identifiable from single labels (so a learned head's residual overconfidence is an objective/
    optimization artifact, not an identifiability barrier)."""
    rng = np.random.default_rng(seed)
    ov, dhat, sand = [], [], []
    for _ in range(int(n_edges)):
        s = rng.uniform(0.6, 3.0)
        xf, xr = _gauss_edge(s, n, rng)
        r = bar_estimate(xf, xr)
        ov.append(_norm_overlap(r.overlap, n, n))
        dhat.append(r.delta_f)
        sand.append(np.sqrt(r.var_sandwich))
    ov, dhat, sand = np.asarray(ov), np.asarray(dhat), np.asarray(sand)
    edges = np.quantile(ov, np.linspace(0.0, 1.0, nbins + 1))
    idx = np.clip(np.digitize(ov, edges[1:-1]), 0, nbins - 1)
    rel = []
    for b in range(nbins):
        m = idx == b
        if m.sum() < 8:
            continue
        pooled = np.sqrt(((dhat[m] - dhat[m].mean()) ** 2).mean())
        target = sand[m].mean()
        if target > 0:
            rel.append(abs(pooled - target) / target)
    return float(np.mean(rel))


def _train_mve_core(n_train=200, n=20, seed=RNG_SEED, use_overlap=True, beta=0.0):
    """Shared trainer: small heteroscedastic NN trained by (beta-)Gaussian NLL on a budget of
    edges (one noisy ΔΔĜ label each). With use_overlap it also receives the overlap scalar I —
    the 'fair foil' that is fed work moments AND overlap (the reviewer's request).

    beta=0.0 reproduces plain Gaussian NLL exactly (byte-for-byte with the original
    implementation). beta>0 reweights each sample's NLL by the stop-gradient predictive
    variance to the beta power (Seitzer et al. 2022, arXiv:2203.09168), which decouples the
    mean-fit/variance-starvation pathology blamed for Gaussian-NLL MVE overconfidence.

    Returns (net, mu, sd, predict) where predict(xf, xr) -> (mean, se).
    """
    import torch

    rng = np.random.default_rng(seed)
    seps = rng.uniform(0.6, 3.0, n_train)
    feats, labels = [], []
    for s in seps:
        xf, xr = _gauss_edge(s, n, rng)
        r = bar_estimate(xf, xr)
        row = [xf.mean(), xr.mean(), xf.std(), xr.std(), abs(xf.mean() - xr.mean())]
        if use_overlap:
            row.append(_norm_overlap(r.overlap, n, n))
        feats.append(row)
        labels.append(r.delta_f)
    X = torch.tensor(np.array(feats), dtype=torch.float32)
    y = torch.tensor(labels, dtype=torch.float32)
    mu, sd = X.mean(0), X.std(0) + 1e-6
    Xn = (X - mu) / sd
    in_dim = X.shape[1]
    torch.manual_seed(seed)
    net = torch.nn.Sequential(
        torch.nn.Linear(in_dim, 32), torch.nn.SiLU(),
        torch.nn.Linear(32, 32), torch.nn.SiLU(), torch.nn.Linear(32, 2),
    )
    opt = torch.optim.Adam(net.parameters(), lr=5e-3)
    for _ in range(1500):
        opt.zero_grad()
        out = net(Xn)
        m, logv = out[:, 0], out[:, 1]
        per = 0.5 * logv + 0.5 * (y - m) ** 2 / torch.exp(logv)
        if beta == 0.0:
            nll = per.mean()
        else:
            weight = torch.exp(logv * beta).detach()  # sigma^{2*beta}, stop-gradient
            nll = (weight * per).mean()
        nll.backward()
        opt.step()

    def predict(xf, xr):
        row = [xf.mean(), xr.mean(), xf.std(), xr.std(), abs(xf.mean() - xr.mean())]
        if use_overlap:
            r = bar_estimate(xf, xr)
            row.append(_norm_overlap(r.overlap, xf.size, xr.size))
        fe = torch.tensor([row], dtype=torch.float32)
        with torch.no_grad():
            out = net((fe - mu) / sd)[0]
            return float(out[0]), float(torch.exp(0.5 * out[1]))

    return predict


def _train_mve(n_train=200, n=20, seed=RNG_SEED, use_overlap=True, beta=0.0):
    """Public entry point: Gaussian-NLL (beta=0, default, unchanged) or beta-NLL (beta>0,
    Seitzer et al. 2022) learned-variance head. Returns predict_se(xf, xr) -> se, matching the
    original signature/behavior exactly at beta=0.0."""
    predict = _train_mve_core(n_train=n_train, n=n, seed=seed, use_overlap=use_overlap, beta=beta)

    def predict_se(xf, xr):
        return predict(xf, xr)[1]

    return predict_se


def _train_mve_ensemble(n_train=200, n=20, seed=RNG_SEED, use_overlap=True, n_members=5):
    """Deep ensemble of plain Gaussian-NLL (beta=0) MVE heads (Lakshminarayanan et al. 2017),
    seeds seed, seed+1, ..., seed+n_members-1. Returns predict_se(xf, xr) -> se using the
    mixture SD: se = sqrt(mean_m(sigma_m^2) + var_m(mu_m)) — aleatoric mean + epistemic
    spread of the member means."""
    members = [
        _train_mve_core(n_train=n_train, n=n, seed=seed + i, use_overlap=use_overlap, beta=0.0)
        for i in range(n_members)
    ]

    def predict_se(xf, xr):
        mus = np.empty(n_members)
        sigmas = np.empty(n_members)
        for i, predict in enumerate(members):
            m, s = predict(xf, xr)
            mus[i], sigmas[i] = m, s
        aleatoric = float(np.mean(sigmas ** 2))
        epistemic = float(np.var(mus))
        return float(np.sqrt(aleatoric + epistemic))

    return predict_se


# --- multi-seed foil spread (peer-review item P6c) ------------------------------------------
# The learned-foil ratios reported in Fig A come from one training seed. These constants are
# FROZEN: the spread is reported over exactly these seeds and this sweep, whatever it turns out
# to be. The Monte-Carlo evaluation stream is deliberately held FIXED across training seeds
# (`eval_seed`), so the across-seed spread isolates training-seed variability rather than
# re-sampling noise in the Monte-Carlo truth.
N_FOIL_SEEDS = 5
FOIL_SEPS = np.linspace(0.8, 3.2, 11)
# Single source of truth for the ensemble's member count. `foil_seed_spread` below binds BOTH
# the disjoint-block seed stride and the `n_members=` argument to this one name, so raising the
# member count can never silently desynchronize the two and reintroduce overlapping seed blocks.
N_ENSEMBLE_MEMBERS = 5


def foil_seed_spread(seps=FOIL_SEPS, n=20, reps=1500, n_seeds=N_FOIL_SEEDS,
                     seed0=RNG_SEED, eval_seed=RNG_SEED):
    """Reported/true-se ratio of each learned foil, per training seed.

    The row is indexed by training-seed index ``i = 0 .. n_seeds-1``. For ``plain``,
    ``oracle``, and ``betanll`` (each a single net) row ``i`` trains at seed ``seed0 + i`` —
    5 independent retrains, one seed apart, as for any single-net foil.

    ``ensemble`` is a 5-member deep ensemble, so a single row cannot reuse the same
    ``seed0 + i`` scheme without seed collisions ACROSS rows: `_train_mve_ensemble(seed=sd,
    n_members=5)` internally trains members at `sd, sd+1, ..., sd+4`, so adjacent rows under
    the single-net scheme would share 4 of their 5 member seeds (a sliding window, not 5
    independent retrains). Instead, row ``i``'s ensemble consumes a DISJOINT block of
    ``n_members`` seeds: base seed ``seed0 + n_members * i``, i.e. member seeds
    ``seed0 + n_members*i .. seed0 + n_members*i + n_members - 1``. With ``n_members=5`` and
    ``seed0=RNG_SEED`` the 5 rows' member blocks are ``{629-633}, {634-638}, {639-643},
    {644-648}, {649-653}`` — no overlap, so the ensemble column's cross-row spread is a
    genuine 5-independent-retrain spread (25 distinct member networks total, not 9).

    For each training seed we retrain all four foils (plain Gaussian-NLL, the large-budget
    oracle, beta-NLL, and the 5-member deep ensemble), then sweep the separations and take
    ``mean(reported se) / MonteCarlo-true se`` per separation, averaging over the sweep. A ratio
    of 1.0 is perfect calibration; BELOW 1.0 is overconfident, and lower is worse.

    Returns one dict per seed with keys ``seed, plain, oracle, betanll, ensemble``.
    """
    rows = []
    for i in range(n_seeds):
        sd = seed0 + i
        foils = {
            "plain": _train_mve(seed=sd, use_overlap=True),
            "oracle": _train_mve(seed=sd, n_train=4000, use_overlap=True),
            "betanll": _train_mve(seed=sd, use_overlap=True, beta=0.5),
            "ensemble": _train_mve_ensemble(seed=seed0 + N_ENSEMBLE_MEMBERS * i, n_train=200,
                                            use_overlap=True, n_members=N_ENSEMBLE_MEMBERS),
        }
        rng = np.random.default_rng(eval_seed)   # fixed across training seeds by design
        acc = {k: [] for k in foils}
        for s in seps:
            dhat = np.empty(reps)
            pred = {k: np.empty(reps) for k in foils}
            for k in range(reps):
                xf, xr = _gauss_edge(s, n, rng)
                dhat[k] = bar_estimate(xf, xr).delta_f
                for name, f in foils.items():
                    pred[name][k] = f(xf, xr)
            true_se = float(dhat.std(ddof=1))
            for name in foils:
                acc[name].append(float(pred[name].mean()) / true_se)
        rows.append({"seed": sd, **{k: float(np.mean(v)) for k, v in acc.items()}})
    return rows


def print_foil_seed_spread():
    """Print the frozen multi-seed table plus mean and min-max interval per foil."""
    rows = foil_seed_spread()
    names = ("plain", "oracle", "betanll", "ensemble")
    print(f"{'seed':>10} " + " ".join(f"{k:>10}" for k in names))
    for r in rows:
        print(f"{r['seed']:>10} " + " ".join(f"{r[k]:>10.4f}" for k in names))
    print()
    for k in names:
        v = np.array([r[k] for r in rows])
        print(f"{k:>10}: mean {v.mean():.4f}  min {v.min():.4f}  max {v.max():.4f}  "
              f"sd {v.std(ddof=1):.4f}")


def controlled_panel(seps, n=20, reps=3000, n_boot=400, seed=RNG_SEED):
    import torch  # noqa: F401  (ensure torch import error surfaces early)

    rng = np.random.default_rng(seed)
    # fair foil, realistic budget
    mve_se = _train_mve(seed=seed, use_overlap=True)
    # fair foil, large budget
    mve_oracle_se = _train_mve(seed=seed, n_train=4000, use_overlap=True)
    # corrected foil #1 (beta-NLL, Seitzer 2022)
    mve_betanll_se = _train_mve(seed=seed, use_overlap=True, beta=0.5)
    # corrected foil #2 (deep ensemble of MVEs)
    mve_ens_se = _train_mve_ensemble(seed=seed, n_train=200, use_overlap=True, n_members=5)
    out = []
    for s in seps:
        dhat = np.empty(reps)
        sand = np.empty(reps)
        mbar = np.empty(reps)
        naive = np.empty(reps)
        mve = np.empty(reps)
        mve_o = np.empty(reps)
        mve_beta = np.empty(reps)
        mve_ens = np.empty(reps)
        ovl = np.empty(reps)
        for k in range(reps):
            xf, xr = _gauss_edge(s, n, rng)
            r = bar_estimate(xf, xr)
            dhat[k] = r.delta_f
            sand[k] = np.sqrt(r.var_sandwich)
            mbar[k] = np.sqrt(max(r.var_mbar, 0.0))
            naive[k] = np.sqrt(r.var_naive)
            mve[k] = mve_se(xf, xr)
            mve_o[k] = mve_oracle_se(xf, xr)
            mve_beta[k] = mve_betanll_se(xf, xr)
            mve_ens[k] = mve_ens_se(xf, xr)
            ovl[k] = _norm_overlap(r.overlap, n, n)
        true_se = dhat.std(ddof=1)
        row = dict(overlap=float(ovl.mean()),
                   sand=float(sand.mean() / true_se), mbar=float(mbar.mean() / true_se),
                   naive=float(naive.mean() / true_se), mve=float(mve.mean() / true_se),
                   mve_oracle=float(mve_o.mean() / true_se),
                   mve_betanll=float(mve_beta.mean() / true_se),
                   mve_ensemble=float(mve_ens.mean() / true_se))
        sr = np.empty(n_boot)
        for b in range(n_boot):
            idx = rng.integers(0, reps, reps)
            sr[b] = sand[idx].mean() / dhat[idx].std(ddof=1)
        row["sand_lo"], row["sand_hi"] = np.percentile(sr, 2.5), np.percentile(sr, 97.5)
        out.append(row)
    out.sort(key=lambda d: d["overlap"])
    return out


def _edges_from_unk(frames):
    """frames: list of (lambda, u_nk DataFrame). Yields (w_f, w_r, lam_i, lam_j)."""
    frames = sorted(frames, key=lambda t: t[0])
    cols = sorted(frames[0][1].columns)
    for i in range(len(frames) - 1):
        ci, cj = cols[i], cols[i + 1]
        u_i, u_j = frames[i][1], frames[i + 1][1]
        yield (u_i[cj] - u_i[ci]).to_numpy(), (u_j[ci] - u_j[cj]).to_numpy(), frames[i][0], frames[i + 1][0]


def _benzene_edges():
    from alchemlyb.parsing.gmx import extract_u_nk
    from alchemtest.gmx import load_benzene
    for leg, files in load_benzene()["data"].items():
        frames = [(float(extract_u_nk(f, T=300.0).index.get_level_values(-1)[0]),
                   extract_u_nk(f, T=300.0)) for f in files]
        yield from _edges_from_unk(frames)


def _bace_edges():
    """Real protein-ligand BINDING FEP edges (BACE1 RBFE, complex legs)."""
    from alchemlyb.parsing.amber import extract_u_nk
    from alchemtest.amber import load_bace_example
    data = load_bace_example()["data"]["complex"]
    for leg, files in data.items():
        frames = []
        for f in files:
            u = extract_u_nk(f, T=298.0)
            frames.append((float(u.index.get_level_values(-1)[0]), u))
        yield from _edges_from_unk(frames)


def real_panel(n_boot=400, seed=RNG_SEED):
    from pymbar.timeseries import statistical_inefficiency, subsample_correlated_data
    rng = np.random.default_rng(seed)
    rows = []
    for kind, gen in [("binding", _bace_edges), ("solvation", _benzene_edges)]:
        for w_f, w_r, lam_i, lam_j in gen():
            try:
                gf, gr = statistical_inefficiency(w_f), statistical_inefficiency(w_r)
            except Exception:
                gf = gr = 1.0
            wf = w_f[subsample_correlated_data(w_f, g=gf)]
            wr = w_r[subsample_correlated_data(w_r, g=gr)]
            if wf.size < 8 or wr.size < 8:
                continue
            r = bar_estimate(wf, -wr)
            dboot = np.array([solve_bar(wf[rng.integers(0, wf.size, wf.size)],
                                        -wr[rng.integers(0, wr.size, wr.size)]) for _ in range(n_boot)])
            true_se = dboot.std(ddof=1)
            if true_se <= 0:
                continue
            rows.append(dict(kind=kind, overlap=_norm_overlap(r.overlap, wf.size, wr.size),
                             sand=np.sqrt(max(r.var_sandwich, 0)) / true_se,
                             naive=np.sqrt(r.var_naive) / true_se))
    return rows


def bace_table_rows(n_boot=400, seed=RNG_SEED):
    """Same computation as `real_panel()`'s ``binding`` branch, but also carries the
    (lambda_i, lambda_j) BAR-window labels through -- so Supplementary Table
    ``tab:bace1`` (I4, whole-branch review) has a programmatic, re-derivable source
    instead of hand-paired rows. Numerically identical sand/naive values to
    `real_panel()`'s binding rows: same RNG_SEED, same edge order from `_bace_edges()`,
    and the binding branch is consumed from the shared rng before the solvation branch
    in `real_panel()`, so a fresh rng(seed) here reproduces the exact same draws."""
    from pymbar.timeseries import statistical_inefficiency, subsample_correlated_data
    rng = np.random.default_rng(seed)
    rows = []
    for w_f, w_r, lam_i, lam_j in _bace_edges():
        try:
            gf, gr = statistical_inefficiency(w_f), statistical_inefficiency(w_r)
        except Exception:
            gf = gr = 1.0
        wf = w_f[subsample_correlated_data(w_f, g=gf)]
        wr = w_r[subsample_correlated_data(w_r, g=gr)]
        if wf.size < 8 or wr.size < 8:
            continue
        r = bar_estimate(wf, -wr)
        dboot = np.array([solve_bar(wf[rng.integers(0, wf.size, wf.size)],
                                    -wr[rng.integers(0, wr.size, wr.size)]) for _ in range(n_boot)])
        true_se = dboot.std(ddof=1)
        if true_se <= 0:
            continue
        rows.append(dict(lam_i=lam_i, lam_j=lam_j,
                         overlap=_norm_overlap(r.overlap, wf.size, wr.size),
                         sand=np.sqrt(max(r.var_sandwich, 0)) / true_se,
                         naive=np.sqrt(r.var_naive) / true_se))
    rows.sort(key=lambda d: d["overlap"])
    return rows


def print_bace_table() -> None:
    """Emit the 19-row BACE1 table (with lambda labels) that Supplementary
    Table `tab:bace1` (`docs/paper_si.tex`) cites; committed verbatim to
    `docs/results_figA.md` so the SI table has a re-derivable provenance."""
    rows = bace_table_rows()
    print(f"[BACE1 table] {len(rows)} edges, sorted by overlap ascending "
          f"(lambda_i, lambda_j, 4<p(1-p)>, sandwich/boot, naive/boot):")
    for d in rows:
        print(f"  {d['lam_i']:.3f}  {d['lam_j']:.3f}  {d['overlap']:.3f}  "
              f"{d['sand']:.3f}  {d['naive']:.3f}")


def main() -> None:
    use_paper_style()
    ctrl = controlled_panel(seps=np.linspace(0.8, 3.2, 11))
    real = real_panel()

    # Panel A carries seven series over a wide ratio range, panel B a narrow scatter, so A
    # gets the wider column. Both panels plot the same quantity -- a reported/true se ratio
    # against overlap -- so both use the same log ratio axis, the same reference line at 1,
    # the same colours, and a key in the same place under the axes.
    fig = plt.figure(figsize=figsize(2, height=3.4))
    gs = fig.add_gridspec(1, 2, width_ratios=(1.3, 1.0),
                          left=0.070, right=0.995, bottom=0.325, top=0.905, wspace=0.30)
    axA = fig.add_subplot(gs[0, 0])
    axB = fig.add_subplot(gs[0, 1])

    ox = np.array([d["overlap"] for d in ctrl])

    def col(key):
        return np.array([d[key] for d in ctrl])

    # --- panel A -------------------------------------------------------------------------
    # The four learned foils are one family (a learned per-edge variance, four ways), so they
    # are drawn inside their own envelope and told apart by dash and marker, not by hue alone.
    foils = (
        ("mve", "MVE head, realistic budget", FOIL, "-", "^", 1.5, 3.8),
        ("mve_oracle", "MVE head, large budget", FOIL, (0, (3.2, 1.4)), "v", 1.1, 3.4),
        # every learned head is the FOIL family by the palette rule; they are told apart by
        # dash pattern and marker, not by hue, so amber is left to mean ALT and nothing else.
        ("mve_betanll", r"$\beta$-NLL head", tint(FOIL, 0.30), (0, (1.1, 1.3)), "s", 1.1, 3.0),
        ("mve_ensemble", "deep ensemble (5 members)", tint(FOIL, 0.55), (0, (4.5, 1.3, 1.0, 1.3)), "D",
         1.1, 2.7),
    )
    stack = np.vstack([col(k) for k, *_ in foils])
    axA.fill_between(ox, stack.min(0), stack.max(0), color=tint(FOIL, 0.86), lw=0, zorder=1.2)
    h_foil = []
    for key, lab, c, ls, mk, lw, ms in foils:
        h_foil += axA.plot(ox, col(key), ls=ls, marker=mk, color=c, lw=lw, ms=ms,
                           mew=0.0, zorder=2.5, label=lab)

    # naive 1/I: the textbook plug-in this corrects, a de-emphasised baseline rather than a
    # rival. It is a NAMED estimator this article corrects, not a null and not de-emphasised
    # data, so it takes THIRD rather than a grey: grey on the facing pages means "random".
    h_naive, = axA.plot(ox, col("naive"), ls=(0, (5, 2)), color=THIRD, lw=1.4, zorder=2.2,
                        label="naive $1/I$")
    # sandwich = MBAR = truth
    axA.fill_between(ox, [d["sand_lo"] for d in ctrl], [d["sand_hi"] for d in ctrl],
                     color=tint(OURS, 0.55), lw=0, zorder=3.5)
    h_sand, = axA.plot(ox, col("sand"), "-", color=OURS, lw=2.0, zorder=4,
                       label="sandwich $B/I^2$")
    h_mbar, = axA.plot(ox, col("mbar"), "o", color=INK, ms=2.8, mew=0.0, zorder=5,
                       label="pymbar-MBAR")
    reference_line(axA, "hline", 1.0)

    axA.set_yscale("log")
    axA.set_xlim(0.13, 0.91)
    axA.set_ylim(0.058, 3.5)
    axA.yaxis.set_major_locator(FixedLocator([0.1, 0.2, 0.5, 1.0, 2.0]))
    axA.yaxis.set_major_formatter(lambda v, _: f"{v:g}")
    axA.yaxis.set_minor_locator(LogLocator(subs=tuple(np.arange(2, 10) * 0.1)))
    axA.yaxis.set_minor_formatter(NullFormatter())
    axA.tick_params(axis="y", which="minor", length=1.6)
    axA.set_xlabel(r"overlap   $4\langle p(1-p)\rangle$")
    axA.set_ylabel("reported se / true se")
    axA.annotate("learned-variance foils", xy=(0.60, 0.235), xycoords="data",
                 ha="center", va="bottom", fontsize=7.5, color=FOIL)
    # Both panels label the reference line the same way: a fixed 9 pt above it in points, on
    # whichever side the data has left free (right here, left in B), so the clearance from the
    # nearest marker does not depend on the axis limits or on the log scale.
    axA.annotate("1 = calibrated", xy=(0.905, 1.0), xycoords="data",
                 xytext=(0, 9), textcoords="offset points",
                 ha="right", va="baseline", fontsize=7, color=REF)
    panel(axA, "A", "controlled Monte-Carlo sweep")

    # One key per panel, under its own axes, so nothing sits on the data. Column 1 is the
    # three estimators, column 2 the four learned foils; the blank pads column 1 to length.
    blank = Line2D([], [], color="none", label=" ")
    legend(axA, handles=[h_sand, h_mbar, h_naive, blank] + h_foil,
           loc="upper left", bbox_to_anchor=(-0.01, -0.185), ncol=2, columnspacing=1.4,
           handlelength=2.0)

    # --- panel B -------------------------------------------------------------------------
    C_LEG = {"binding": OURS, "solvation": tint(OURS, 0.45)}
    for kind, mk in (("binding", "o"), ("solvation", "D")):
        pts = [d for d in real if d["kind"] == kind]
        if not pts:
            continue
        ob = [d["overlap"] for d in pts]
        axB.scatter(ob, [d["naive"] for d in pts], c=THIRD, s=20, marker="^", zorder=2,
                    linewidths=0.0, alpha=0.85)
        # A point is one adjacent-lambda BAR window, not a ligand-to-ligand free energy; the
        # legend says which leg the window comes from rather than calling it an edge.
        leg = {"binding": "BACE1 complex leg", "solvation": "benzene solvation leg"}[kind]
        axB.scatter(ob, [d["sand"] for d in pts], c=C_LEG[kind], s=26, marker=mk, zorder=3,
                    linewidths=0.4, edgecolors="white", label=f"sandwich \u2014 {leg}")
    axB.scatter([], [], c=THIRD, marker="^", s=20, linewidths=0.0, label="naive $1/I$")
    reference_line(axB, "hline", 1.0)

    axB.set_yscale("log")
    axB.set_ylim(0.75, 16.0)
    axB.yaxis.set_major_locator(FixedLocator([1.0, 2.0, 5.0, 10.0]))
    axB.yaxis.set_major_formatter(lambda v, _: f"{v:g}")
    axB.yaxis.set_minor_locator(LogLocator(subs=tuple(np.arange(2, 10) * 0.1)))
    axB.yaxis.set_minor_formatter(NullFormatter())
    axB.tick_params(axis="y", which="minor", length=1.6)
    axB.set_xlabel(r"overlap   $4\langle p(1-p)\rangle$")
    axB.set_ylabel("reported se / bootstrap se")
    # Panel B's marker cluster runs the full width of the y = 1 line, so the label cannot sit
    # at axes-right the way panel A's does: it is anchored at axes-left instead, a fixed 9 pt
    # above the line in points, which clears the tallest marker there whatever the data limits.
    axB.annotate("1 = calibrated", xy=(0.02, 1.0), xycoords=("axes fraction", "data"),
                 xytext=(0, 9), textcoords="offset points",
                 ha="left", va="baseline", fontsize=7, color=REF)
    panel(axB, "B", "real adjacent-\u03bb BAR windows")
    legend(axB, loc="upper left", bbox_to_anchor=(-0.01, -0.185), handlelength=1.0)

    finish(fig, "figA_target_the_sandwich")

    print("\n[Panel A] overlap  sand/true  MBAR/true  fair-foil/true  oracle/true  "
          "beta-NLL/true  ensemble/true  naive/true")
    for d in ctrl:
        print(f"  {d['overlap']:.3f}    {d['sand']:.3f}     {d['mbar']:.3f}      "
              f"{d['mve']:.3f}          {d['mve_oracle']:.3f}       "
              f"{d['mve_betanll']:.3f}         {d['mve_ensemble']:.3f}         {d['naive']:.3f}")
    print(f"  -> sandwich == MBAR (max |Δ| = {max(abs(d['sand']-d['mbar']) for d in ctrl):.3f}); "
          f"fair foil {np.mean([d['mve'] for d in ctrl]):.2f}x of true se "
          f"(large-budget {np.mean([d['mve_oracle'] for d in ctrl]):.2f}x)")
    print(f"  -> corrected foils at realistic budget: beta-NLL "
          f"{np.mean([d['mve_betanll'] for d in ctrl]):.2f}x of true se; "
          f"ensemble-of-MVEs {np.mean([d['mve_ensemble'] for d in ctrl]):.2f}x of true se "
          f"(sandwich {np.mean([d['sand'] for d in ctrl]):.2f}x = reference)")
    for kind in ("binding", "solvation"):
        pts = [d for d in real if d["kind"] == kind]
        if pts:
            sr = np.array([d["sand"] for d in pts])
            nr = np.array([d["naive"] for d in pts])
            print(f"[Panel B] {kind:9s}: {len(pts)} edges  sandwich/boot {sr.mean():.2f} "
                  f"[{sr.min():.2f},{sr.max():.2f}]  naive/boot {nr.mean():.2f} [{nr.min():.2f},{nr.max():.2f}]")

    print("[identifiability] pooled se(overlap) recovery from single labels: "
          + ", ".join(f"N={N}:{100*pooled_se_recovery(N):.1f}%" for N in (200, 4000, 40000)))


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "seeds":
        print_foil_seed_spread()
    elif len(sys.argv) > 1 and sys.argv[1] == "bace_table":
        print_bace_table()
    else:
        main()
