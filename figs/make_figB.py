"""Fig B — uncertainty decomposition & OOD calibration (the differentiator).

Total Belief variance = aleatoric (sandwich B/I^2) + epistemic (ensemble) + Delta
(force-field). On test edges sorted by trunk-space domain distance, the physically
DECOMPOSED per-edge uncertainty holds conditional coverage, while a marginal split-
conformal band (flat width) does not — even though BOTH are conformal-calibrated to
90% MARGINAL coverage. The decomposition also shows *why* you are uncertain
(epistemic + Delta dominate OOD); conformal gives a black-box flat band.

Fair comparison (both have the conformal marginal guarantee):
  ours          mu +/- q_norm * sigma_total   (normalized / Mondrian conformal; per-edge)
  split-conf    mu +/- q                       (standard split conformal; flat width)
  const-sigma   mu +/- 1.645 * rms             (single-number heuristic)

Run:  python figs/make_figB.py   (or `make figB`)
Falsifier (plan Fig B): split-conformal matches per-edge conditional coverage.
Honest scope (plan §9): the sandwich is aleatoric-only; OOD coverage is the ensemble's
and Delta-head's job — this is the TEST, not a guarantee.
"""
from __future__ import annotations

import pathlib
import sys
import warnings

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import torch

warnings.filterwarnings("ignore")
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "src"))

FIGDIR = pathlib.Path(__file__).resolve().parent
C_OURS, C_CONF, C_HEUR, C_REF = "#0072B2", "#D55E00", "#999999", "#555555"
ALPHA = 0.10  # 90% prediction interval


def _style() -> None:
    mpl.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 9, "axes.labelsize": 10, "axes.titlesize": 10,
        "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 7.5,
        "axes.spines.top": False, "axes.spines.right": False,
        "figure.dpi": 120, "savefig.dpi": 300, "savefig.bbox": "tight",
    })


def domain_distance(x):
    """Trunk-space distance from the training support [-2, 2]."""
    return np.maximum(np.abs(x) - 2.0, 0.0)


def sigma_aleatoric(x):
    """Known heteroscedastic aleatoric (the BAR sandwich), grows with |x|."""
    return 0.12 + 0.10 * np.abs(x)


def sigma_delta(x):
    """Delta-head force-field uncertainty: grows with domain distance (novel chemistry)."""
    return 0.55 * domain_distance(x) ** 1.3


def true_mean(x):
    return np.sin(1.5 * x) + 0.3 * x


def ff_bias(x):
    """Systematic force-field error, OOD only (zero on training support)."""
    return 0.45 * np.sign(x) * domain_distance(x) ** 1.3


def sample(x, rng):
    return true_mean(x) + ff_bias(x) + rng.normal(0, sigma_aleatoric(x))


def make_net(seed):
    torch.manual_seed(seed)
    return torch.nn.Sequential(
        torch.nn.Linear(1, 24), torch.nn.Tanh(),
        torch.nn.Linear(24, 24), torch.nn.Tanh(), torch.nn.Linear(24, 1),
    )


def train_member(x, y, seed, epochs=300):
    net = make_net(seed)
    opt = torch.optim.Adam(net.parameters(), lr=8e-3, weight_decay=1e-4)
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, len(x), len(x))  # bootstrap
    xt = torch.tensor(x[idx, None], dtype=torch.float32)
    yt = torch.tensor(y[idx, None], dtype=torch.float32)
    lossf = torch.nn.MSELoss()
    for _ in range(epochs):
        opt.zero_grad(); loss = lossf(net(xt), yt); loss.backward(); opt.step()
    return net


def ensemble_predict(nets, x):
    xt = torch.tensor(x[:, None], dtype=torch.float32)
    with torch.no_grad():
        preds = np.stack([n(xt).numpy().ravel() for n in nets])  # (M, N)
    return preds.mean(0), preds.std(0, ddof=1)


def main() -> None:
    _style()
    rng = np.random.default_rng(7)
    # ensemble trained on the IN-DISTRIBUTION support only -> real OOD extrapolation
    x_tr = rng.uniform(-2, 2, 200); y_tr = sample(x_tr, rng)
    # calibration and test are EXCHANGEABLE (both span the full range), so split
    # conformal gets its ~90% MARGINAL guarantee fairly; we then probe CONDITIONAL
    # coverage across domain distance.
    x_cal = rng.uniform(-5.5, 5.5, 500); y_cal = sample(x_cal, rng)
    x_te = rng.uniform(-5.5, 5.5, 4000); y_te = sample(x_te, rng)

    nets = [train_member(x_tr, y_tr, s) for s in range(12)]

    def components(x):
        mu, s_ens = ensemble_predict(nets, x)
        return mu, sigma_aleatoric(x), s_ens, sigma_delta(x)

    mu_cal, sa_cal, se_cal, sd_cal = components(x_cal)
    st_cal = np.sqrt(sa_cal**2 + se_cal**2 + sd_cal**2)
    mu_te, sa_te, se_te, sd_te = components(x_te)
    st_te = np.sqrt(sa_te**2 + se_te**2 + sd_te**2)

    # conformal calibration (all give ~90% MARGINAL coverage by construction)
    res_cal = np.abs(y_cal - mu_cal)
    qn = 1.0 - ALPHA
    k = min(int(np.ceil((len(x_cal) + 1) * qn)) / len(x_cal), 1.0)
    q_norm = np.quantile(res_cal / st_cal, k)             # ours: normalized conformal (per-edge σ_total)
    q_flat = np.quantile(res_cal, k)                      # split conformal: flat
    # CQR: base quantiles use the ALEATORIC σ only (what an adaptive-conformal user has);
    # the conformal adjustment is a flat scalar -> adapts to aleatoric, NOT to epistemic/Δ.
    z = 1.645
    e_cal = np.maximum((mu_cal - z * sa_cal) - y_cal, y_cal - (mu_cal + z * sa_cal))
    q_cqr = np.quantile(e_cal, k)

    covered = {
        "ours": np.abs(y_te - mu_te) <= q_norm * st_te,
        "CQR": (y_te >= mu_te - z * sa_te - q_cqr) & (y_te <= mu_te + z * sa_te + q_cqr),
        "split-conformal": np.abs(y_te - mu_te) <= q_flat,
    }
    width = {
        "ours": 2 * q_norm * st_te,
        "CQR": 2 * (z * sa_te + q_cqr),
        "split-conformal": np.full_like(x_te, 2 * q_flat),
    }

    d_te = domain_distance(x_te)
    bins = np.linspace(0, d_te.max(), 9)
    bc = 0.5 * (bins[:-1] + bins[1:])
    which = np.digitize(d_te, bins) - 1
    which = np.clip(which, 0, len(bc) - 1)

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(7.6, 3.3))

    # Panel A: conditional coverage vs domain distance
    for name, c, mk in [("ours", C_OURS, "o"), ("CQR", "#CC79A7", "D"), ("split-conformal", C_CONF, "s")]:
        cov = np.array([covered[name][which == b].mean() if np.any(which == b) else np.nan for b in range(len(bc))])
        axA.plot(bc, cov, "-" + mk, color=c, ms=4, lw=1.6, label=name)
    axA.axhline(0.90, color=C_REF, ls="--", lw=1.0)
    axA.set_xlabel("trunk-space domain distance (OOD →)")
    axA.set_ylabel("conditional coverage (90% PI)")
    axA.set_title("A   coverage holds OOD vs conformal", loc="left", fontweight="bold")
    axA.legend(frameon=False, loc="lower left")
    axA.set_ylim(0.3, 1.02)

    # Panel B: variance decomposition vs domain distance + conformal flat half-width
    order = np.argsort(d_te)
    dd = d_te[order]
    axB.plot(dd, sa_te[order], color="#56B4E9", lw=1.6, label="aleatoric (sandwich)")
    axB.plot(dd, se_te[order], color="#009E73", lw=1.6, label="epistemic (ensemble)")
    axB.plot(dd, sd_te[order], color="#CC79A7", lw=1.6, label="Δ (force-field)")
    axB.plot(dd, (q_norm * st_te)[order], color=C_OURS, lw=1.8, label="ours half-width")
    axB.axhline(q_flat, color=C_CONF, ls="--", lw=1.4, label="split-conformal (flat)")
    axB.set_xlabel("trunk-space domain distance (OOD →)")
    axB.set_ylabel("uncertainty (kcal/mol)")
    axB.set_title("B   decomposition: epistemic+Δ grow OOD", loc="left", fontweight="bold")
    axB.legend(frameon=False, loc="upper left", ncol=1)
    sharp = width["split-conformal"].mean() / width["ours"].mean()
    axB.text(0.97, 0.05, f"same 90% marginal cov;\nours {sharp:.1f}× sharper",
             transform=axB.transAxes, fontsize=7.5, ha="right", color=C_OURS)

    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(FIGDIR / f"figB_ood_decomposition.{ext}")
    print(f"wrote figB_ood_decomposition.(pdf|png) to {FIGDIR}")

    print("\nmarginal coverage (target 0.90):")
    for name in covered:
        print(f"  {name:16s} {covered[name].mean():.3f}   mean width {width[name].mean():.2f}")
    far = d_te > 1.5
    print("\nFAR-OOD (domain distance > 1.5) conditional coverage:")
    for name in covered:
        print(f"  {name:16s} {covered[name][far].mean():.3f}")
    # conditional-coverage error: max |bin coverage - 0.90|
    print("\nmax |bin coverage - 0.90| (lower = better conditional calibration):")
    for name in covered:
        cov = np.array([covered[name][which == b].mean() for b in range(len(bc)) if np.any(which == b)])
        print(f"  {name:16s} {np.max(np.abs(cov - 0.90)):.3f}")


if __name__ == "__main__":
    main()
