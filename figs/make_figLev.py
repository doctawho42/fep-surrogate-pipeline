"""Fig Lev -- per-edge observability map + the pre-registered predictive falsifier (D1).

Structural (must pass exactly, per system): sum_e h_e == dof; bridge edges have h==0.
Predictive (pre-registered): does the curl-leverage h_e predict where *reproducible* systematic
error concentrates? Confound: under the sampling null Var(z_e)=h_e, so a raw Spearman(h,|z|) is
positive by construction. Fix: studentize z_tilde = z/sqrt(h) (unit null variance; systematic
signal still scales as sqrt(h)), take the cross-replicate reproducible magnitude
S_e = |mean_k z_tilde|, and test Spearman(h_e, S_e) against a WITHIN-SYSTEM block-permutation
null. Pre-registered:
  CONFIRMED if rho>0 with permutation p<0.05  -> observability predicts detectable systematic error;
  KILL if p>=0.05                             -> drop only the predictive claim; the structural
                                                 theorem + bridge/auditability map stand regardless.

Run:  PYTHONPATH=src python figs/make_figLev.py   (or `make figLev`).  Deterministic.
"""
from __future__ import annotations

import csv
import math
import pathlib
import sys
from collections import defaultdict

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from bar.leverage import bridges, curl_leverage  # noqa: E402
from bar.qc import gls_network  # noqa: E402
from paperstyle import (  # noqa: E402
    MUTED, OURS, figsize, finish, panel, use_paper_style,
)

DATA = ROOT / "data" / "openfe_replicates" / "combined_pymbar4_edge_data.csv"
FIGDIR = pathlib.Path(__file__).resolve().parent
FLAGGED = {"brd4", "bace", "faah", "cdk8", "hif2a", "p38"}
H_MIN = 0.05
N_PERM = 10000
# Semantic colours (paperstyle). The six systems the calibrated cycle-closure detector flags are
# the detector's own output, so they carry OURS -- the same blue that marks them in Fig L; the
# sampling-consistent majority is present but not the point, so it is MUTED. The binned
# reproducible-residual series in B is likewise computed from the calibrated per-edge bar and is
# OURS whatever the verdict turns out to be: the hue names the quantity, the heading names the
# verdict. FOIL is reserved for the overconfident stand-in and does not appear here.


def _f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return math.nan


def edge_val(r, k):
    cD = _f(r[f"complex_repeat_{k}_DG (kcal/mol)"])
    cd = _f(r[f"complex_repeat_{k}_dDG (kcal/mol)"])
    sD = _f(r[f"solvent_repeat_{k}_DG (kcal/mol)"])
    sd = _f(r[f"solvent_repeat_{k}_dDG (kcal/mol)"])
    if any(math.isnan(v) for v in (cD, cd, sD, sd)):
        return None
    return cD - sD, math.sqrt(cd * cd + sd * sd)


def load_systems():
    by = defaultdict(list)
    for r in csv.DictReader(open(DATA)):
        by[r["system name"]].append(r)
    return by


def _edges_rep(rows, k):
    return [(r["ligand_A"], r["ligand_B"], *edge_val(r, k)) for r in rows if edge_val(r, k)]


def structural(by):
    """Per-system sum_h==dof and bridge stats on replicate 0."""
    rows = []
    for sysname, rs in sorted(by.items()):
        e = _edges_rep(rs, 0)
        if len(e) < 3:
            continue
        fit = gls_network(e)
        if fit.dof < 1:
            continue
        h = curl_leverage(e)
        br = bridges(e)
        rows.append({
            "sys": sysname, "E": len(e), "dof": fit.dof, "sum_h": float(h.sum()),
            "n_bridge": len(br), "frac_audit": float(np.mean(h >= H_MIN)),
            "h_med": float(np.median(h)),
        })
    return rows


def predictive(by):
    """Studentized reproducibility S_e vs h_e, pooled with per-system block info for permutation."""
    h_all, s_all, sys_all = [], [], []
    for sysname, rs in sorted(by.items()):
        reps = [_edges_rep(rs, k) for k in (0, 1, 2)]
        if any(len(e) < 3 for e in reps):
            continue
        # key edges present in all three replicates (undirected key)
        maps = []
        for e in reps:
            m = {frozenset((a, b)): (a, b, dd, se) for a, b, dd, se in e}
            maps.append(m)
        common = set(maps[0]) & set(maps[1]) & set(maps[2])
        if len(common) < 3:
            continue
        # per replicate: fit, z, h aligned to the common key order
        keys = sorted(common, key=lambda fs: sorted(map(str, fs)))
        per_rep_z, per_rep_h, ok = [], [], True
        for k in range(3):
            e = [maps[k][key] for key in keys]
            fit = gls_network(e)
            if fit.dof < 1:
                ok = False
                break
            per_rep_z.append(fit.z)
            per_rep_h.append(curl_leverage(e))
        if not ok:
            continue
        z = np.vstack(per_rep_z)            # 3 x m
        h = np.mean(np.vstack(per_rep_h), axis=0)   # m; per-rep h varies ~9%, same KILL
        good = h >= H_MIN
        if good.sum() < 3:
            continue
        zt = z[:, good] / np.sqrt(h[good])          # studentize -> unit null variance
        s = np.abs(np.mean(zt, axis=0))             # reproducible systematic magnitude
        h_all.extend(h[good].tolist())
        s_all.extend(s.tolist())
        sys_all.extend([sysname] * int(good.sum()))
    return np.array(h_all), np.array(s_all), np.array(sys_all, dtype=object)


def _spearman(a, b):
    ra = np.argsort(np.argsort(a))
    rb = np.argsort(np.argsort(b))
    ra = ra - ra.mean()
    rb = rb - rb.mean()
    denom = math.sqrt(float(np.sum(ra * ra)) * float(np.sum(rb * rb)))
    return float(np.sum(ra * rb) / denom) if denom > 0 else 0.0


def block_permutation_p(h, s, sysnames, n_perm=N_PERM, seed=0):
    """One-sided p for rho(h,s) under within-system permutation of s against h."""
    rng = np.random.default_rng(seed)
    obs = _spearman(h, s)
    groups = defaultdict(list)
    for i, g in enumerate(sysnames):
        groups[g].append(i)
    idx_groups = [np.array(v) for v in groups.values()]
    ge = 0
    for _ in range(n_perm):
        sp = s.copy()
        for g in idx_groups:
            sp[g] = s[rng.permutation(g)]
        if _spearman(h, sp) >= obs:
            ge += 1
    return obs, (ge + 1) / (n_perm + 1)


def main():
    use_paper_style()
    by = load_systems()
    srows = structural(by)
    # structural assertion: sum_h == dof to machine precision
    max_dev = max(abs(r["sum_h"] - r["dof"]) for r in srows)
    assert max_dev < 1e-6, f"sum_h != dof somewhere (max dev {max_dev:.2e})"

    h, s, sysnames = predictive(by)
    rho, pval = block_permutation_p(h, s, sysnames)
    verdict = "CONFIRMED" if (rho > 0 and pval < 0.05) else "KILL"

    fig, (axA, axB) = plt.subplots(1, 2, figsize=figsize(2, 3.0))
    # A: leverage distribution + bridge fraction per system (flagged highlighted)
    order = sorted(srows, key=lambda r: r["h_med"])
    ys = np.arange(len(order), dtype=float)
    hm = np.array([r["h_med"] for r in order])
    fl = np.array([r["sys"] in FLAGGED for r in order])
    axA.scatter(hm[~fl], ys[~fl], c=MUTED, s=10, lw=0, alpha=0.8, zorder=3)
    axA.scatter(hm[fl], ys[fl], c=OURS, s=18, lw=0, zorder=4)
    for r, yv in zip(order, ys, strict=True):
        if r["sys"] in FLAGGED:
            lab = f"{r['sys']} ({r['n_bridge']} bridge{'' if r['n_bridge'] == 1 else 's'})"
            axA.annotate(lab, (r["h_med"], yv), xytext=(5.5, -0.5),
                         textcoords="offset points", fontsize=7, va="center", color=OURS,
                         zorder=5)
    # the names run to the right of their own points and the rightmost one is the longest, so
    # the right limit is opened enough to hold it inside the axes.
    axA.set_xlim(hm.min() - 0.018, hm.max() + 0.075)
    axA.set_ylim(-1.5, len(order) + 0.5)
    axA.set_xlabel(r"median curl-leverage $h_e$ (replicate 0)")
    axA.set_yticks([])
    axA.set_ylabel(f"{len(order)} systems (ranked)")
    panel(axA, "A", "per-system observability",
          r"bridges are un-auditable ($h_e{=}0$)")
    # B: predictive falsifier -- binned S_e vs h_e + verdict
    nb = 8
    qs = np.quantile(h, np.linspace(0, 1, nb + 1))
    qs[-1] += 1e-9
    xs, ms = [], []
    for i in range(nb):
        m = (h >= qs[i]) & (h < qs[i + 1])
        if m.sum():
            xs.append(float(np.mean(h[m])))
            ms.append(float(np.mean(s[m])))
    axB.plot(xs, ms, marker="o", ms=4.0, color=OURS, lw=1.4, zorder=3)
    axB.set_xlabel(r"curl-leverage $h_e$")
    axB.set_ylabel(r"reproducible $|{\rm mean}_k \tilde z_e|$")
    panel(axB, "B",
          f"predictive falsifier: {'supported' if verdict == 'CONFIRMED' else 'null result'}",
          rf"Spearman $\rho={rho:+.2f}$, perm $p={pval:.3f}$")

    finish(fig, "figLev_observability")
    print(f"[structural] {len(srows)} systems; sum_h==dof max dev {max_dev:.2e}; "
          f"total bridges {sum(r['n_bridge'] for r in srows)}")
    print(f"[predictive] n_edges={len(h)}; Spearman rho={rho:+.3f}; perm p={pval:.4f}")
    print(f"VERDICT: {verdict}")


if __name__ == "__main__":
    main()
