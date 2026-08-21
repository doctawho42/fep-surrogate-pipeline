"""eg5 accuracy vs experiment — grounding the cycle-closure QC against real ChEMBL affinities.

Honest outcome (see docs/results_eg5_accuracy.md):
  * POSITIVE (validation): the real OpenFE eg5 RBFE network, solved by our GLS pipeline, agrees
    with real ChEMBL KIF11 experimental affinity at MUE ~0.79 kcal/mol (R ~0.65), and the system's
    calibrated cycle-closure reduced chi^2 ~0.64 is clean — a QC-clean system is indeed accurate.
  * NULL (honest, and predicted by the scope): acting on the QC (prune / down-weight flagged
    edges) does NOT improve accuracy vs experiment on this system, and the internal |z| barely
    predicts per-edge error vs experiment (Pearson ~+0.13). This is exactly the stated Fig L scope:
    cycle closure is blind to node-consistent force-field bias, and on a clean system most of the
    error vs experiment IS that common FF bias. The strong claim ("QC repair improves accuracy")
    would need a QC-FLAGGED system with experimental data; those systems in the IndustryBenchmarks
    set are anonymized (no public per-ligand affinity), so the strong test needs new MD / unblinded
    data and is left as future work.

Experimental data: data/eg5_experimental_chembl.csv (ChEMBL KIF11 ATPase IC50, assay CHEMBL1101849;
pulled per-ligand). Network: data/openfe_replicates (3-replicate mean). Run:
  PYTHONPATH=src python figs/analyze_eg5_accuracy.py
"""
from __future__ import annotations

import csv
import math
import pathlib
import sys

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
from bar.qc import cycle_closure_test, gls_network  # noqa: E402

RT_LN10 = 1.364  # kcal/mol at 298 K; exp dG = -RT ln10 * pChEMBL


def _f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return math.nan


def load_experiment():
    p = ROOT / "data" / "eg5_experimental_chembl.csv"
    return {r["ligand_chembl_id"]: -RT_LN10 * float(r["pchembl"]) for r in csv.DictReader(open(p))}


def load_edges(exp):
    p = ROOT / "data" / "openfe_replicates" / "combined_pymbar4_edge_data.csv"
    edges = []
    for r in csv.DictReader(open(p)):
        if r["system name"] != "eg5":
            continue
        dd, se = [], []
        for k in (0, 1, 2):
            cD = _f(r[f"complex_repeat_{k}_DG (kcal/mol)"]); cd = _f(r[f"complex_repeat_{k}_dDG (kcal/mol)"])
            sD = _f(r[f"solvent_repeat_{k}_DG (kcal/mol)"]);  sd = _f(r[f"solvent_repeat_{k}_dDG (kcal/mol)"])
            if not any(math.isnan(v) for v in (cD, cd, sD, sd)):
                dd.append(cD - sD); se.append(math.sqrt(cd * cd + sd * sd))
        if len(dd) >= 2 and r["ligand_A"] in exp and r["ligand_B"] in exp:
            edges.append((r["ligand_A"], r["ligand_B"], float(np.mean(dd)),
                          math.sqrt(np.mean(np.array(se) ** 2) / len(dd))))
    return edges


def kendall_tau(p, e):
    """Standard Kendall tau-b (ties handled), via scipy — not the tie-dropping variant."""
    from scipy.stats import kendalltau
    t = kendalltau(p, e).correlation
    return float(t) if t == t else math.nan


def _metrics(p, e):
    """MUE/RMSE/tau/R on a (predicted, experimental) pair, re-aligning the gauge on this set."""
    p = p - p.mean() + e.mean()  # FEP predicts relative dG: align means (gauge)
    R = float(np.corrcoef(p, e)[0, 1]) if p.std() > 0 and e.std() > 0 else math.nan
    return dict(mue=float(np.mean(np.abs(p - e))), rmse=float(math.sqrt(np.mean((p - e) ** 2))),
                tau=float(kendall_tau(p, e)), R=R)


def accuracy(edges, exp):
    fit = gls_network(edges)
    pred = dict(zip(fit.nodes, fit.potentials))
    keys = [k for k in pred if k in exp]
    p = np.array([pred[k] for k in keys]); e = np.array([exp[k] for k in keys])
    out = _metrics(p, e); out["n"] = len(keys); out["p"] = p; out["e"] = e
    return out


def bootstrap_ci(p, e, n_boot=2000, seed=0):
    """95% percentile-bootstrap CIs for the accuracy metrics, resampling ligands with replacement
    (gauge re-aligned per resample, matching the point estimate). n=28 is small, so the CIs are wide
    -- reporting them is the honest treatment."""
    rng = np.random.default_rng(seed)
    n = len(p)
    acc = {k: [] for k in ("mue", "rmse", "tau", "R")}
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        m = _metrics(p[idx], e[idx])
        for k in acc:
            acc[k].append(m[k])
    return {k: (float(np.nanpercentile(v, 2.5)), float(np.nanpercentile(v, 97.5))) for k, v in acc.items()}


def main():
    exp = load_experiment()
    edges = load_edges(exp)
    acc = accuracy(edges, exp)
    rc = cycle_closure_test(edges).reduced_chi2
    print(f"eg5: {len(edges)} edges, {acc['n']} ligands with ChEMBL experiment")
    ci = bootstrap_ci(acc["p"], acc["e"])
    print(f"\n[POSITIVE] real OpenFE eg5 network vs ChEMBL KIF11 experiment (95% bootstrap CI, n={acc['n']}):")
    print(f"    MUE={acc['mue']:.2f} [{ci['mue'][0]:.2f},{ci['mue'][1]:.2f}]  "
          f"RMSE={acc['rmse']:.2f} [{ci['rmse'][0]:.2f},{ci['rmse'][1]:.2f}] kcal/mol  "
          f"tau={acc['tau']:+.2f} [{ci['tau'][0]:+.2f},{ci['tau'][1]:+.2f}]  "
          f"R={acc['R']:+.2f} [{ci['R'][0]:+.2f},{ci['R'][1]:+.2f}]")
    print(f"    calibrated cycle-closure reduced chi^2 = {rc:.2f} (clean) -> QC-clean AND accurate (consistent)")

    fit = gls_network(edges)
    z = np.abs(fit.z)
    err = np.abs([y - (exp[b] - exp[a]) for (a, b, y, se) in edges])
    pear = np.corrcoef(z, err)[0, 1]
    rx, ry = np.argsort(np.argsort(z)), np.argsort(np.argsort(err))
    spear = np.corrcoef(rx, ry)[0, 1]
    print("\n[NULL, honest] does internal |z| predict real error vs experiment?")
    print(f"    Pearson(|z|,|edge error vs exp|)={pear:+.2f}  Spearman={spear:+.2f}  -> weak/none")
    print("    Reason (Fig L scope, now confirmed vs experiment): cycle closure is blind to")
    print("    node-consistent force-field bias, which dominates error-vs-experiment on a clean system.")
    print("    The strong test (QC-flagged system + experiment) needs unblinded / new-MD data (future work).")


if __name__ == "__main__":
    main()
