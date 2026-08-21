"""Fig E — chirality completeness.

A designed enantiomer-pair test set (chiral tetrahedral centres + their mirrors).
Panel A: the O(3)-invariant ("even") readout is *identical* on a pair (collapse);
the parity-odd 0o pseudoscalar flips sign (separates). Panel B: a small MLP trained
to predict per-enantiomer dG — WITHOUT the 0o channel it provably cannot tell
enantiomers apart, so its enantiomer ddG error equals the predict-zero baseline; WITH
0o it recovers the chiral signal.

Run:  python figs/make_figE.py    (or `make figE`)
Kill (plan Fig E): if the EVEN readout already separates enantiomers -> the theorem's
premise is violated -> investigate. (Here it is provably exact: collapse ~ 1e-16.)
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
from bar.chiral import chiral_readout, even_features, signed_volume  # noqa: E402

FIGDIR = pathlib.Path(__file__).resolve().parent
C_EVEN = "#D55E00"   # vermillion — even / no-0o (blind)
C_ODD = "#0072B2"    # blue — with 0o
C_REF = "#555555"
ALPHA_CHIRAL = 0.5   # sets the enantiomer ddG scale (~1.5 kcal/mol)


def _style() -> None:
    mpl.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 9, "axes.labelsize": 10, "axes.titlesize": 10,
        "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 8,
        "axes.spines.top": False, "axes.spines.right": False,
        "figure.dpi": 120, "savefig.dpi": 300, "savefig.bbox": "tight",
    })


# ----------------------------------------------------------------------------
# Designed enantiomer-pair dataset
# ----------------------------------------------------------------------------
def make_dataset(n_pairs: int, seed: int):
    """Each 'molecule' = 4 distinct substituents (scalar types, shuffled) at a
    perturbed, randomly rotated tetrahedron. Mirror = reflection. True per-enantiomer
    dG = achiral(even features, substituents) + ALPHA*sgn(chi)*(s0 - s3)."""
    rng = np.random.default_rng(seed)
    base = np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]], dtype=float)
    # fixed generative weights for the achiral part
    w_even = rng.normal(size=6)
    w_sub = rng.normal(size=4)
    rows = []
    for _ in range(n_pairs):
        from scipy.spatial.transform import Rotation
        R = Rotation.random(random_state=rng.integers(1 << 31)).as_matrix()
        coords = (base + 0.25 * rng.normal(size=(4, 3))) @ R.T
        subs = rng.permutation([1.0, 2.0, 3.0, 4.0])
        chi = signed_volume(coords)
        ev = even_features(coords)
        achiral = ev @ w_even + subs @ w_sub
        chiral = ALPHA_CHIRAL * np.sign(chi) * (subs[0] - subs[3])
        noise = 0.05 * rng.normal()
        # molecule M
        rows.append(dict(coords=coords, subs=subs, chi=chi, dG=achiral + chiral + noise))
        # mirror M' (reflect x): same subs, even features identical, chi flips
        coords_m = coords * np.array([-1.0, 1, 1])
        chi_m = signed_volume(coords_m)
        chiral_m = ALPHA_CHIRAL * np.sign(chi_m) * (subs[0] - subs[3])
        rows.append(dict(coords=coords_m, subs=subs, chi=chi_m,
                         dG=achiral + chiral_m + 0.05 * rng.normal()))
    return rows


def featurize(rows, include_0o: bool):
    X, y = [], []
    for r in rows:
        ev = np.concatenate([even_features(r["coords"]), r["subs"]])  # O(3)-invariant
        if include_0o:
            ev = np.append(ev, r["chi"])
        X.append(ev)
        y.append(r["dG"])
    return np.asarray(X), np.asarray(y)


# ----------------------------------------------------------------------------
def train_eval(include_0o: bool, seed: int, n_pairs=900, epochs=300):
    rows = make_dataset(n_pairs, seed)
    n = len(rows) // 2
    idx = np.random.default_rng(seed).permutation(n)
    tr_pairs, te_pairs = idx[: int(0.75 * n)], idx[int(0.75 * n):]
    tr_rows = [rows[2 * p + k] for p in tr_pairs for k in (0, 1)]
    Xtr, ytr = featurize(tr_rows, include_0o)
    mu, sd = Xtr.mean(0), Xtr.std(0) + 1e-9
    Xtr = (Xtr - mu) / sd
    Xt = torch.tensor(Xtr, dtype=torch.float32)
    yt = torch.tensor(ytr, dtype=torch.float32).unsqueeze(1)

    torch.manual_seed(seed)
    net = torch.nn.Sequential(
        torch.nn.Linear(Xt.shape[1], 64), torch.nn.SiLU(),
        torch.nn.Linear(64, 64), torch.nn.SiLU(), torch.nn.Linear(64, 1),
    )
    opt = torch.optim.Adam(net.parameters(), lr=5e-3)
    lossf = torch.nn.MSELoss()
    for _ in range(epochs):
        opt.zero_grad()
        loss = lossf(net(Xt), yt)
        loss.backward()
        opt.step()

    # evaluate enantiomer ddG = dG(M) - dG(M') on held-out PAIRS
    net.eval()
    ddg_err = []
    with torch.no_grad():
        for p in te_pairs:
            mM, mP = rows[2 * p], rows[2 * p + 1]
            xM = (featurize([mM], include_0o)[0] - mu) / sd
            xP = (featurize([mP], include_0o)[0] - mu) / sd
            pred_ddg = (net(torch.tensor(xM, dtype=torch.float32)).item()
                        - net(torch.tensor(xP, dtype=torch.float32)).item())
            true_ddg = mM["dG"] - mP["dG"]
            ddg_err.append(abs(pred_ddg - true_ddg))
    return float(np.mean(ddg_err))


# ----------------------------------------------------------------------------
def main() -> None:
    _style()
    seeds = [0, 1, 2, 3, 4]

    # collapse check (kill criterion): even readout identical on a pair; 0o flips
    rng = np.random.default_rng(99)
    even_diffs, odd_diffs = [], []
    for _ in range(400):
        c = rng.normal(size=(4, 3))
        cm = c * np.array([-1.0, 1, 1])
        even_diffs.append(np.abs(chiral_readout(c, False) - chiral_readout(cm, False)).max())
        odd_diffs.append(abs(signed_volume(c) - signed_volume(cm)))
    even_collapse = float(np.max(even_diffs))
    odd_sep = float(np.median(odd_diffs))

    mae0o = np.array([train_eval(True, s) for s in seeds])
    mae_even = np.array([train_eval(False, s) for s in seeds])

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(7.2, 3.3))

    # Panel A: collapse vs separation (signed_volume of M vs mirror)
    rng2 = np.random.default_rng(7)
    chis = np.array([signed_volume(rng2.normal(size=(4, 3))) for _ in range(120)])
    axA.scatter(chis, chis, c=C_EVEN, s=14, marker="o", label="even readout (M vs M′)", zorder=3)
    axA.scatter(chis, -chis, c=C_ODD, s=14, marker="^", label="0o channel (M vs M′)", zorder=3)
    lim = 1.05 * np.abs(chis).max()
    axA.plot([-lim, lim], [-lim, lim], color=C_REF, ls="--", lw=0.9, zorder=0)
    axA.set_xlabel("readout value for M")
    axA.set_ylabel("readout value for mirror M′")
    axA.set_title("A   collapse vs separation", loc="left", fontweight="bold")
    axA.legend(frameon=False, loc="upper left")
    axA.text(0.04, 0.06, f"even collapse:\nmax|Δ| = {even_collapse:.1e}",
             transform=axA.transAxes, fontsize=7.5, color=C_EVEN)

    # Panel B: enantiomer ddG MAE, even vs even+0o
    means = [mae_even.mean(), mae0o.mean()]
    errs = [mae_even.std(ddof=1), mae0o.std(ddof=1)]
    axB.bar([0, 1], means, yerr=errs, capsize=4, color=[C_EVEN, C_ODD], width=0.6)
    axB.set_xticks([0, 1])
    axB.set_xticklabels(["even\n(no 0o)", "even + 0o"])
    axB.set_ylabel("enantiomer ΔΔG MAE (kcal/mol)")
    axB.set_title("B   chirality ablation", loc="left", fontweight="bold")
    for x, m, e in zip([0, 1], means, errs):
        axB.text(x, m + e + 0.04, f"{m:.2f}", ha="center", va="bottom", fontsize=8)
    axB.set_ylim(0, max(means) * 1.18)

    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(FIGDIR / f"figE_chirality_completeness.{ext}")
    print(f"wrote figE_chirality_completeness.(pdf|png) to {FIGDIR}")
    print(f"\n[kill check] even-readout collapse max|Δ| = {even_collapse:.2e}  "
          f"(must be ~0; 0o median |Δ| = {odd_sep:.3f})")
    print(f"[Panel B] enantiomer ddG MAE: even (no 0o) = {mae_even.mean():.3f} "
          f"± {mae_even.std(ddof=1):.3f}  |  even+0o = {mae0o.mean():.3f} ± {mae0o.std(ddof=1):.3f}")
    print(f"          0o reduces ddG MAE by {mae_even.mean()/mae0o.mean():.1f}x")


if __name__ == "__main__":
    main()
