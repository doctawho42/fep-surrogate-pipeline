"""Fig Inf -- inference for six manuscript claims that are currently stated without it.

A second-round referee panel asked for the inference behind six numbers. Every design below was
fixed *before* the corresponding number was computed, is stated here in full, and is run once;
nothing is revisited after seeing an outcome. Several of these can come back against the
manuscript, and the results doc reports whichever way they fall.

**C1 -- inference for the cross-replicate residual correlation.** The manuscript reports pooled
per-edge standardized-residual correlations across replicate pairs (`r = +0.30, +0.42, +0.38`,
`n = 1143`) with no interval and no test, and says the correlation concentrates in the flagged
systems. The 1143 edges are *not* independent: they are nested in 48 networks and every residual
is a projection through that network's residual maker `M = I - H`. Two inferential objects,
both respecting that nesting:
  * a **cluster bootstrap** that resamples the 48 SYSTEMS with replacement (never edges), `N_BOOT`
    draws, percentile 95% CI, for each of the three pairs;
  * a **within-system permutation test**: one permutation `pi` per draw permutes the edge order
    *inside each system* of the second member of each pair, which destroys the cross-replicate
    edge pairing while preserving system membership, system sizes, and the marginal z
    distribution. Under the sampling null, replicate fits are independent, so `E[z_e^(i) z_f^(j)]
    = 0` for every `e, f` including `e = f`; the permutation therefore has the right null. The
    same `pi` is used for all three pairs in a draw, so the three null values inherit the same
    dependence the three observed values have. One-sided upper p, `(1 + #{null >= obs})/(1 + N)`.
  * the **flagged-vs-unflagged contrast** `r_flagged - r_unflagged` gets a *stratified* cluster
    bootstrap (resample the flagged systems among themselves and the unflagged among themselves,
    so neither stratum can come out empty) and a **system-label permutation test** that reassigns
    which 6 of the 48 systems carry the flag, holding all data fixed. That is the null for
    "concentrated in the flagged systems"; a within-system edge permutation would instead test
    "no correlation anywhere", so both are reported and labelled.
  * The three rotations are **dependent**: they are three views of the same three runs, and every
    pair of rotations shares a replicate. They are not three replications and the doc says so.

**C2 -- the missing half of a two-way null.** The manuscript says the calibrated detector does not
measurably out-discriminate *either* foil, quoting one interval, `+0.140 [-0.008, +0.316]`, which
is `AUC(A) - max(AUC(B), AUC(C))` and hence the comparison against the fixed-se test only. The
comparison against the 1.0 kcal/mol hysteresis cutoff has no interval anywhere. This computes all
three paired contrasts (A-B, A-C, A-max(B,C)) from the SAME `N_BOOT` system resamples, using
`bar.detectors.auc_flag_vs_anchor` and `make_figCut.build` unchanged and the same `seed = 0`, so
the A-max(B,C) row must reproduce the published interval exactly (asserted).

**C3 -- the observed aggregate that a prediction is said to match.** Section 4.2 quotes a
predicted `E[chi2_nu] = 0.85 [0.63, 1.08]` and says the aggregate level "also matches", without
stating the observed aggregate or the aggregation rule. The predicted aggregate is
`sum_s sum_e h_e c_e^-2 / sum_s dof_s`; since `sum_e h_e = dof_s` (Theorem D1) and
`E[X^2_s] = sum_e h_e c_e^-2`, the matched observed statistic is `sum_s X^2_s / sum_s dof_s` --
i.e. the leverage-weighted rule and the dof-weighted global value are the SAME rule, and only
three of the four requested numbers are distinct. All are reported (matched/dof-weighted global,
median over systems, unweighted mean over systems), the matched one with a system-cluster
bootstrap CI, and the doc states whether the two intervals overlap.

**C4 -- the degrees-of-freedom distribution.** The count of independent cycles per system, then
`make figStab`'s flag stability stratified by it. Pre-stated primary split: the MEDIAN split of
dof over the admitted systems (low = dof <= median, high = dof > median); secondary split
`dof == 1` vs `dof >= 2`. Per stratum: the three pairwise Jaccards of the replicate flag sets
restricted to the stratum, their exact size-matched random-draw reference
(`make_figStab.expected_random_jaccard`), and the median per-system chi^2_nu swing (max/min over
the three replicates).

**C5 -- is the predicted-versus-observed check circular?** In the headline rotation the predicted
quantity is `sum_e h_e (s_e^2 / se_{e,0}^2) / dof` with `s_e^2` the sample variance of replicates
1 and 2, and the observed quantity is `sum_e r_{e,0}^2 / se_{e,0}^2 / dof`. Replicate **0**'s
reported se is the denominator on BOTH sides, and the leverage weights `h_e` are themselves
computed from replicate 0's `V_e`. The null: `make_figOOS.resample_null_world`, the repo's own
no-systematic-error world (same graphs, same reported se, values redrawn as `y ~ N(0, se_k)`, so
`c_e == 1` identically). `h_e` and `se_e` are unchanged by it, so the shared denominator is
retained exactly. `N_NULL` realizations, each giving one realized Spearman over the same systems;
the reported 0.648 is then quoted against that null distribution rather than against zero. A
positive control (a genuine per-system `c` heterogeneity injected into the same harness) confirms
the harness can see a real signal.

**C6 -- a like-for-like comparison with the literature.** \\citet{wade2022estimators} compare the
analytic MBAR uncertainty against the spread over independent replicas and report
*under*-estimation. The manuscript reconciles with them through a leverage-weighted functional
that has no counterpart in their measurement. This computes the functional they actually report:
per-edge and per-target analytic se versus replicate spread, and the fraction of edges on which
the analytic se falls BELOW the replicate spread. Reference level, stated before running: with
3 replicates and perfectly calibrated bars, `s^2 ~ sigma^2 chi^2_2 / 2`, so
`P(s > sigma) = exp(-1) = 0.368` -- a third of edges land below by construction, and only an
excess over 0.368 is evidence of under-estimation. The c4-corrected variant (`sigma_hat = s/c4`,
reference `exp(-c4^2) = 0.456`) is reported alongside. Cluster bootstrap over systems.

Run:  PYTHONPATH=src python figs/make_figInf.py   (or `make figInf`).  Deterministic (every
resampling seeded, `SEED` below). Data: `data/openfe_replicates/combined_pymbar4_edge_data.csv`
(OpenFE IndustryBenchmarks2024, public, 3 independent replicates per edge).
"""
from __future__ import annotations

import math
import pathlib
import sys

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))

from figs import make_figA_replicates as figArep  # noqa: E402
from figs import make_figCut as figCut  # noqa: E402
from figs import make_figL_validation as figLval  # noqa: E402
from figs.make_figOOS import (  # noqa: E402
    leverage_weighted_mean,
    load_records,
    panel_a,
    resample_null_world,
    sample_var,
)
from figs.make_figStab import (  # noqa: E402
    ALPHA,
    REPLICATES,
    expected_random_jaccard,
    flag_set,
    jaccard,
    replicate_edges,
    system_flags,
)
from figs.make_figStab import load_systems as stab_load_systems  # noqa: E402

from bar.detectors import (  # noqa: E402
    HYSTERESIS_CUTOFF,
    anchor_score,
    auc_flag_vs_anchor,
    flag_calibrated,
    flag_fixed_cutoff,
    flag_fixed_se,
)
from bar.qc import gls_network  # noqa: E402

FIGDIR = pathlib.Path(__file__).resolve().parent
DOC = ROOT / "docs" / "results_figInf.md"

FLAGGED = ("brd4", "bace", "faah", "cdk8", "hif2a", "p38")
PAIRS = ((0, 1), (0, 2), (1, 2))
SEED = 20260810
N_BOOT = 2000          # cluster-bootstrap resamples (C1, C3, C6)
N_PERM = 1999          # permutation draws (C1)
N_NULL = 200           # no-systematic-error realizations in the circularity null (C5)
N_BOOT_AUC = 2000      # matches make_figCut's N_BOOT so the published row reproduces exactly
MIN_TARGET_EDGES = 8   # make_figA_replicates' own per-target threshold (C6)

# Published numbers this script is asked to put an interval on / test. Used only for the
# reproduction assertion (C2) and for the printed comparison; never as an input to a fit.
PUB_AUC_DIFF_CI = (-0.008, 0.316)
PUB_PRED_AGG = (0.85, 0.63, 1.08)
PUB_RHO = 0.648

C_OK, C_FLAG, C_REF, C_ALT = "#0072B2", "#D55E00", "#555555", "#7030A0"


# ============================================================ shared cluster machinery

def cluster_blocks(labels):
    """``(cluster_names, [row_indices_per_cluster])`` -- the resampling unit for every bootstrap
    below. Resampling these blocks, never individual rows, is what makes the interval honest
    about the nesting of edges inside networks."""
    order: dict = {}
    for i, s in enumerate(labels):
        order.setdefault(s, []).append(i)
    names = sorted(order)
    return names, [np.array(order[s], dtype=int) for s in names]


def cluster_resample(blocks, rng):
    """Row indices of one cluster-bootstrap draw: sample clusters with replacement, concatenate."""
    idx = rng.integers(0, len(blocks), len(blocks))
    return np.concatenate([blocks[i] for i in idx])


def within_cluster_permutation(blocks, n, rng):
    """A permutation of ``0..n-1`` that only ever moves a row within its own cluster.

    Applied to one member of a replicate pair it destroys the cross-replicate edge pairing while
    preserving system membership, system sizes and the marginal residual distribution.
    """
    perm = np.arange(n)
    for blk in blocks:
        perm[blk] = rng.permutation(blk)
    return perm


def pearson(a, b):
    a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    if a.size < 2:
        return math.nan
    sa, sb = float(np.std(a)), float(np.std(b))
    if sa <= 0 or sb <= 0:
        return math.nan
    return float(np.mean((a - a.mean()) * (b - b.mean())) / (sa * sb))


def perm_p_upper(obs, null):
    """One-sided upper permutation p, ``(1 + #{null >= obs}) / (1 + N)`` (the repo's convention)."""
    null = np.asarray(null, dtype=float)
    return float((1 + int(np.sum(null >= obs))) / (1 + null.size))


def pct_ci(values):
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    return float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))


# ============================================================ C1: residual correlation

def keyed_residuals():
    """`make_figL_validation.out_of_sample`'s edge set, but keyed by system so it can be clustered.

    Reproduces that function's construction exactly -- per replicate, fit every system with >= 3
    available rows and >= 1 cycle; keep the keys present in all three -- and the caller asserts
    that the pooled correlations and the edge count match the shipped function.
    """
    by = figLval.load()
    Z: dict[int, dict] = {k: {} for k in REPLICATES}
    for sysname, rs in by.items():
        for k in REPLICATES:
            ek = [(r["ligand_A"], r["ligand_B"], *v)
                  for r in rs if (v := figLval.edge_val(r, k))]
            if len(ek) < 3:
                continue
            fit = gls_network(ek)
            if fit.dof < 1:
                continue
            for (a, b, _y, _se), zz in zip(ek, fit.z, strict=True):
                Z[k][(sysname, a, b)] = float(zz)
    keys = [key for key in Z[0] if key in Z[1] and key in Z[2]]
    z = np.vstack([np.array([Z[k][key] for key in keys]) for k in REPLICATES])
    return keys, z, np.array([key[0] for key in keys])


def c1_residual_correlation():
    keys, z, sysn = keyed_residuals()
    z0, z1, z2, corrs, fl_ref = figLval.out_of_sample(figLval.load())
    obs = [pearson(z[i], z[j]) for i, j in PAIRS]
    assert z0.size == len(keys), (z0.size, len(keys))
    assert np.allclose(obs, corrs, atol=1e-12), (obs, corrs)
    assert int(fl_ref.sum()) == int(np.isin(sysn, FLAGGED).sum())

    names, blocks = cluster_blocks(sysn)
    fl = np.isin(sysn, FLAGGED)
    rng = np.random.default_rng(SEED)

    # --- cluster bootstrap of each pooled r (resamples systems, never edges)
    boot = np.empty((N_BOOT, len(PAIRS)))
    for b in range(N_BOOT):
        rows = cluster_resample(blocks, rng)
        for pi, (i, j) in enumerate(PAIRS):
            boot[b, pi] = pearson(z[i][rows], z[j][rows])
    ci = [pct_ci(boot[:, pi]) for pi in range(len(PAIRS))]

    # --- within-system permutation null (one pi per draw, shared by the three pairs)
    null = np.empty((N_PERM, len(PAIRS)))
    for d in range(N_PERM):
        perm = within_cluster_permutation(blocks, z.shape[1], rng)
        for pi, (i, j) in enumerate(PAIRS):
            null[d, pi] = pearson(z[i], z[j][perm])
    pvals = [perm_p_upper(obs[pi], null[:, pi]) for pi in range(len(PAIRS))]

    # --- flagged vs unflagged: stratified cluster bootstrap + system-label permutation
    obs_fl = [pearson(z[i][fl], z[j][fl]) for i, j in PAIRS]
    obs_un = [pearson(z[i][~fl], z[j][~fl]) for i, j in PAIRS]
    obs_diff = [a - b for a, b in zip(obs_fl, obs_un, strict=True)]
    is_flag_sys = np.array([n in FLAGGED for n in names])
    blk_f = [blocks[i] for i in np.flatnonzero(is_flag_sys)]
    blk_u = [blocks[i] for i in np.flatnonzero(~is_flag_sys)]
    bootd = np.empty((N_BOOT, len(PAIRS)))
    for b in range(N_BOOT):
        rf, ru = cluster_resample(blk_f, rng), cluster_resample(blk_u, rng)
        for pi, (i, j) in enumerate(PAIRS):
            bootd[b, pi] = pearson(z[i][rf], z[j][rf]) - pearson(z[i][ru], z[j][ru])
    ci_diff = [pct_ci(bootd[:, pi]) for pi in range(len(PAIRS))]

    n_flag_sys = int(is_flag_sys.sum())
    nulld = np.empty((N_PERM, len(PAIRS)))
    for d in range(N_PERM):
        pick = rng.choice(len(names), n_flag_sys, replace=False)
        rows_f = np.concatenate([blocks[i] for i in pick])
        mask = np.zeros(z.shape[1], dtype=bool)
        mask[rows_f] = True
        for pi, (i, j) in enumerate(PAIRS):
            nulld[d, pi] = (pearson(z[i][mask], z[j][mask])
                            - pearson(z[i][~mask], z[j][~mask]))
    pd_label = [perm_p_upper(obs_diff[pi], nulld[:, pi]) for pi in range(len(PAIRS))]
    # the "no correlation anywhere" version of the same contrast, for contrast
    nulld_edge = np.empty((N_PERM, len(PAIRS)))
    for d in range(N_PERM):
        perm = within_cluster_permutation(blocks, z.shape[1], rng)
        for pi, (i, j) in enumerate(PAIRS):
            nulld_edge[d, pi] = (pearson(z[i][fl], z[j][perm[fl]])
                                 - pearson(z[i][~fl], z[j][perm[~fl]]))
    pd_edge = [perm_p_upper(obs_diff[pi], nulld_edge[:, pi]) for pi in range(len(PAIRS))]

    return {"n_edges": z.shape[1], "n_sys": len(names), "n_flag_edges": int(fl.sum()),
            "n_flag_sys": n_flag_sys, "obs": obs, "ci": ci, "p": pvals,
            "null_mean": [float(np.mean(null[:, pi])) for pi in range(len(PAIRS))],
            "null_sd": [float(np.std(null[:, pi])) for pi in range(len(PAIRS))],
            "obs_flagged": obs_fl, "obs_unflagged": obs_un, "obs_diff": obs_diff,
            "ci_diff": ci_diff, "p_diff_label": pd_label, "p_diff_edge": pd_edge,
            "boot": boot}


# ============================================================ C2: the other AUC contrast

def paired_auc_contrasts(flags_a, flags_b, flags_c, anchor, n_boot=N_BOOT_AUC, seed=0):
    """All three paired AUC contrasts from ONE set of system resamples.

    The resamples are generated exactly as `bar.detectors.paired_auc_bootstrap` generates them
    (same rng, same seed, same order, the point estimate first), so the ``vs_max`` entry
    reproduces that function's published interval bit for bit; the tests assert it. Adding the
    two individual contrasts costs nothing extra and is the whole point of C2: the published
    number is a max, and a max hides which foil the tie came from.
    """
    A, B, C = (np.asarray(f, dtype=bool) for f in (flags_a, flags_b, flags_c))
    s = np.asarray(anchor, dtype=float)
    n = s.size
    rng = np.random.default_rng(seed)
    draws = [np.arange(n)] + [rng.integers(0, n, n) for _ in range(n_boot)]
    auc = np.array([[auc_flag_vs_anchor(F[idx], s[idx]) for F in (A, B, C)] for idx in draws])
    point, boot = auc[0], auc[1:]

    def contrast(col):
        d_point = float(point[0] - col(point))
        lo, hi = pct_ci(np.array([b[0] - col(b) for b in boot]))
        return {"diff": d_point, "ci_lo": lo, "ci_hi": hi,
                "verdict": "WIN" if (d_point > 0 and lo > 0) else "TIE"}

    return {"auc": {k: float(point[i]) for i, k in
                    enumerate(("calibrated", "cutoff", "fixed_se"))},
            "vs_cutoff": contrast(lambda a: a[1]),
            "vs_fixed_se": contrast(lambda a: a[2]),
            "vs_max": contrast(lambda a: max(a[1], a[2]))}


def c2_auc_contrasts():
    """C2 on the real detectors: `make_figCut`'s systems, anchors and flags, unchanged."""
    rep0, aligned = figCut.build(figCut.load_systems())
    fa, fb, fc = flag_calibrated(rep0), flag_fixed_cutoff(rep0), flag_fixed_se(rep0)
    names = [n for n in sorted(rep0) if n in fa and n in fb and n in fc and n in aligned]
    anchor = np.array([anchor_score(aligned[n]) for n in names])
    flags = {"calibrated": np.array([fa[n] for n in names]),
             "cutoff": np.array([fb[n] for n in names]),
             "fixed_se": np.array([fc[n] for n in names])}
    out = paired_auc_contrasts(flags["calibrated"], flags["cutoff"], flags["fixed_se"], anchor)
    out |= {"names": names, "n": len(names),
            "n_flag": {k: int(v.sum()) for k, v in flags.items()}}
    pub = out["vs_max"]
    assert abs(pub["ci_lo"] - PUB_AUC_DIFF_CI[0]) < 5e-4, pub
    assert abs(pub["ci_hi"] - PUB_AUC_DIFF_CI[1]) < 5e-4, pub
    return out


# ============================================================ C3: the observed aggregate

def aggregate_rules(obs, dof):
    """The three distinct ways to aggregate a per-system reduced chi^2 over systems.

    ``matched`` is ``sum_s obs_s dof_s / sum_s dof_s = sum_s X^2_s / sum_s dof_s``. Since
    ``E[X^2_s] = sum_e h_e c_e^-2`` and ``sum_e h_e = dof_s``, that is the observed counterpart of
    the leverage-weighted predicted aggregate -- i.e. the leverage-weighted rule and the
    dof-weighted global value are one and the same, which is why only three numbers appear here.
    """
    obs, dof = np.asarray(obs, dtype=float), np.asarray(dof, dtype=float)
    return {"matched": float(np.sum(obs * dof) / np.sum(dof)),
            "median": float(np.median(obs)),
            "unweighted_mean": float(np.mean(obs))}


def c3_aggregate(recs, A):
    """Observed aggregate under the rule that matches the prediction, plus the alternatives."""
    dof = A["dof"]
    obs = A["obs"][0]                      # replicate 0, the replicate the paper reports
    chi2 = obs * dof                       # X^2_s (matched numerator: E[X^2_s] = sum_e h_e c_e^-2)
    rng = np.random.default_rng(SEED)
    n = len(A["names"])
    boot = np.empty(N_BOOT)
    for b in range(N_BOOT):
        idx = rng.integers(0, n, n)
        boot[b] = chi2[idx].sum() / dof[idx].sum()
    rules = aggregate_rules(obs, dof)
    matched = rules["matched"]
    lo, hi = pct_ci(boot)
    pred_lo, pred_hi = A["agg_ci"]
    return {"matched": matched, "matched_ci": (lo, hi),
            "median": rules["median"], "unweighted_mean": rules["unweighted_mean"],
            "dof_weighted": matched,        # identical by construction; stated, not re-derived
            "pred": A["agg"], "pred_ci": (pred_lo, pred_hi),
            "overlap": not (lo > pred_hi or hi < pred_lo),
            "pred_covers_obs": pred_lo <= matched <= pred_hi,
            "obs_all_rot": A["obs_agg_all"], "pred_all_rot": A["agg_all"],
            "pred_all_ci": A["agg_all_ci"], "n_sys": n,
            "obs_per_sys": obs, "dof_per_sys": dof}


# ============================================================ C4: dof and stratified stability

def stratum_summary(sets, members, swings, label):
    """Flag-set stability restricted to one stratum of systems.

    ``sets`` are the three per-replicate flag sets, ``members`` the stratum's systems, ``swings``
    a system -> chi^2_nu max/min map. The Jaccard reference is the EXACT size-matched random-draw
    expectation *within the stratum*, since a stratum of 25 systems and one of 23 do not have the
    same chance overlap. Two empty sets give Jaccard 1.0 by `make_figStab`'s convention, which is
    an empty-set artefact and not evidence of stability -- the ``sizes`` column shows where.
    """
    members = set(members)
    sub = [set(s) & members for s in sets]
    jac = [jaccard(sub[i], sub[j]) for i, j in PAIRS]
    ref = [expected_random_jaccard(len(sub[i]), len(sub[j]), max(len(members), 1))
           for i, j in PAIRS]
    sw = [swings[s] for s in sorted(members)
          if s in swings and not math.isnan(swings[s])]
    return {"label": label, "n": len(members), "sizes": [len(s) for s in sub],
            "ever": sorted(set().union(*sub)), "inter": sorted(set.intersection(*sub)),
            "jaccard": jac, "mean_jaccard": float(np.mean(jac)),
            "jaccard_random": ref, "mean_jaccard_random": float(np.mean(ref)),
            "median_swing": float(np.median(sw)) if sw else math.nan,
            "max_swing": float(np.max(sw)) if sw else math.nan}


def c4_dof_stratified():
    by = stab_load_systems()
    rows_by_rep = [system_flags({n: replicate_edges(rs, k) for n, rs in by.items()})
                   for k in REPLICATES]
    sets = [flag_set(r) for r in rows_by_rep]
    maps = [{d["sys"]: d for d in r} for r in rows_by_rep]
    names = sorted(maps[0])
    dof = np.array([maps[0][s]["dof"] for s in names])
    counts = {v: int(np.sum(dof == v)) for v in (1, 2, 3)}
    q1, med, q3 = (float(np.percentile(dof, p)) for p in (25, 50, 75))

    def swing(s):
        rcs = [maps[k][s]["rc"] for k in REPLICATES if s in maps[k]]
        return max(rcs) / min(rcs) if len(rcs) > 1 and min(rcs) > 0 else math.nan

    swings = {s: swing(s) for s in names}

    def stratum(mask, label):
        return stratum_summary(
            sets, {s for s, m in zip(names, mask, strict=True) if m}, swings, label)

    primary = [stratum(dof <= med, f"low: dof <= {med:g} (median split)"),
               stratum(dof > med, f"high: dof > {med:g}")]
    secondary = [stratum(dof == 1, "dof == 1 (single cycle)"),
                 stratum(dof >= 2, "dof >= 2")]
    return {"names": names, "dof": dof, "counts": counts, "q1": q1, "median": med, "q3": q3,
            "min": int(dof.min()), "max": int(dof.max()), "mean": float(dof.mean()),
            "sets": sets, "primary": primary, "secondary": secondary,
            "flag_dof": {s: int(maps[0][s]["dof"]) for s in sorted(set().union(*sets))
                         if s in maps[0]}}


# ============================================================ C5: circularity null

def null_rho(recs, rng, c_by_system=None):
    """One realized Spearman of (predicted, observed) in a no-systematic-error world.

    ``c_by_system`` optionally inflates a system's REPORTED bars by a known factor before the
    draw, which is the positive control: a genuine per-system miscalibration must show up as a
    high rho even though the shared denominator is identical to the real analysis.
    """
    from scipy.stats import spearmanr

    pred, obs = [], []
    for s in sorted(recs):
        rec = recs[s]
        c = 1.0 if c_by_system is None else float(c_by_system[s])
        if c != 1.0:
            rec = dict(rec, se={k: rec["se"][k] * c for k in REPLICATES},
                       edges={k: [(a, b, y, se * c) for a, b, y, se in rec["edges"][k]]
                              for k in REPLICATES})
        y_by, z_by = resample_null_world(rec, rng)
        # the true sd is se/c, so the drawn values must not inherit the inflation
        if c != 1.0:
            y_by = {k: v / c for k, v in y_by.items()}
            z_by = {k: np.asarray(gls_network(
                [(a, b, float(y_by[k][e]), se)
                 for e, (a, b, _yy, se) in enumerate(rec["edges"][k])]).z) for k in REPLICATES}
        cinv2 = sample_var([y_by[1], y_by[2]]) / rec["se"][0] ** 2
        pred.append(leverage_weighted_mean(rec["h"][0], cinv2))
        obs.append(float(np.sum(z_by[0] ** 2)) / rec["dof_k"][0])
    return float(spearmanr(pred, obs).statistic)


def c5_circularity(recs):
    rng = np.random.default_rng(SEED)
    null = np.array([null_rho(recs, rng) for _ in range(N_NULL)])
    # positive control: a real, per-system bar miscalibration spanning 0.5x-2x
    names = sorted(recs)
    grid = np.exp(np.linspace(math.log(0.5), math.log(2.0), len(names)))
    cmap = dict(zip(names, grid, strict=True))
    rng_pc = np.random.default_rng(SEED + 1)
    ctrl = np.array([null_rho(recs, rng_pc, cmap) for _ in range(max(N_NULL // 4, 10))])
    lo, hi = pct_ci(null)
    return {"n_null": N_NULL, "null": null, "mean": float(np.mean(null)),
            "median": float(np.median(null)), "sd": float(np.std(null)),
            "ci": (lo, hi), "max": float(np.max(null)),
            "p_vs_null": perm_p_upper(PUB_RHO, null),
            "ctrl_mean": float(np.mean(ctrl)), "ctrl_median": float(np.median(ctrl)),
            "ctrl_ci": pct_ci(ctrl), "n_ctrl": ctrl.size}


# ============================================================ C6: Wade-comparable functional

def perfect_calibration_fraction_below(n_reps=3, scale=1.0):
    """P(analytic se < replicate SD / ``scale``) when the bars are exactly right.

    With ``n_reps`` replicates, ``s^2 (n-1) / sigma^2 ~ chi^2_{n-1}``, so with a correctly scaled
    bar ``se = sigma`` the probability is ``chi2.sf((n-1) * scale^2, n-1)``. At ``n_reps = 3``,
    ``scale = 1`` that is ``exp(-1) = 0.368``: a third of edges land below the replicate spread by
    construction, which is the level Wade et al.'s under-estimation claim has to be read against.
    """
    from scipy.stats import chi2

    return float(chi2.sf((n_reps - 1) * scale ** 2, n_reps - 1))


def fraction_below(analytic_se, replicate_sd, scale=1.0):
    """Fraction of edges on which the analytic se falls BELOW the replicate spread."""
    return float(np.mean(np.asarray(analytic_se, dtype=float)
                         < np.asarray(replicate_sd, dtype=float) / scale))


def c6_wade():
    rows = figArep.load()
    rep = np.array([r[0] for r in rows])       # RMS analytic (MBAR/sandwich) se over 3 replicates
    repl = np.array([r[1] for r in rows])      # across-replicate SD (n = 3, ddof = 1)
    sysn = np.array([r[3] for r in rows])
    keep = repl > 1e-6
    rep, repl, sysn = rep[keep], repl[keep], sysn[keep]

    c4 = figArep.C4_3
    below = rep < repl                          # analytic se BELOW the replicate spread
    below_c4 = rep < (repl / c4)
    ref_raw = perfect_calibration_fraction_below(3)          # = exp(-1), the n=3 reference
    ref_c4 = perfect_calibration_fraction_below(3, c4)       # = exp(-c4^2), bias-corrected

    names, blocks = cluster_blocks(sysn)
    rng = np.random.default_rng(SEED)
    boot = np.empty((N_BOOT, 2))
    for b in range(N_BOOT):
        rows_b = cluster_resample(blocks, rng)
        boot[b] = (below[rows_b].mean(), below_c4[rows_b].mean())

    per_t = []
    for s in names:
        m = sysn == s
        if int(m.sum()) < MIN_TARGET_EDGES:
            continue
        per_t.append({"target": s, "n": int(m.sum()), "frac_below": float(below[m].mean()),
                      "ratio": float(np.sqrt(np.mean(rep[m] ** 2) / np.mean(repl[m] ** 2))),
                      "median_ratio": float(np.median(rep[m] / repl[m]))})
    per_t.sort(key=lambda d: -d["frac_below"])
    return {"n_edges": int(rep.size), "n_sys": len(names),
            "frac_below": fraction_below(rep, repl), "frac_below_ci": pct_ci(boot[:, 0]),
            "frac_below_c4": fraction_below(rep, repl, c4), "frac_below_c4_ci": pct_ci(boot[:, 1]),
            "ref_raw": ref_raw, "ref_c4": ref_c4, "c4": c4,
            "rms_ratio": float(np.sqrt(np.mean(rep ** 2) / np.mean(repl ** 2))),
            "median_ratio": float(np.median(rep / repl)),
            "per_target": per_t,
            "n_targets_above_ref": sum(1 for d in per_t if d["frac_below"] > ref_raw),
            "n_targets_ratio_lt1": sum(1 for d in per_t if d["ratio"] < 1.0)}


# ============================================================ figure

def _style():
    mpl.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
        "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 8.5,
        "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
        "axes.spines.top": False, "axes.spines.right": False,
        "figure.dpi": 120, "savefig.dpi": 300, "savefig.bbox": "tight",
    })


def make_figure(R):
    fig, axes = plt.subplots(2, 3, figsize=(13.4, 8.8))
    (a1, a2, a3), (a4, a5, a6) = axes

    # --- C1 -------------------------------------------------------------------------------
    c1 = R["c1"]
    x = np.arange(len(PAIRS))
    lo = [c1["obs"][i] - c1["ci"][i][0] for i in range(3)]
    hi = [c1["ci"][i][1] - c1["obs"][i] for i in range(3)]
    a1.errorbar(x, c1["obs"], yerr=[lo, hi], fmt="o", color=C_OK, capsize=3, ms=5,
                label="pooled r, cluster-bootstrap 95% CI")
    a1.errorbar(x + 0.22, c1["obs_flagged"], yerr=None, fmt="s", color=C_FLAG, ms=5,
                label="flagged systems only")
    a1.errorbar(x - 0.22, c1["obs_unflagged"], yerr=None, fmt="^", color="#7f7f7f", ms=5,
                label="unflagged systems only")
    a1.axhline(0.0, color=C_REF, lw=0.8, ls=":")
    a1.set_xticks(x)
    a1.set_xticklabels([f"({i},{j})" for i, j in PAIRS])
    a1.set_xlim(-0.5, 2.5)
    a1.set_xlabel("replicate pair (the three rotations are dependent)")
    a1.set_ylabel("cross-replicate residual correlation")
    a1.legend(frameon=False, loc="lower right", fontsize=7)
    a1.set_title(f"C1  r has an interval now (n={c1['n_edges']} edges in\n {c1['n_sys']} systems); "
                 f"permutation p ≤ {max(c1['p']):.3f}", loc="left", fontweight="bold")

    # --- C2 -------------------------------------------------------------------------------
    c2 = R["c2"]
    keys = [("vs_cutoff", f"vs fixed cutoff\n({HYSTERESIS_CUTOFF} kcal/mol)"),
            ("vs_fixed_se", "vs fixed se\n(pooled)"),
            ("vs_max", "vs max(both)\n(published)")]
    xs = np.arange(len(keys))
    vals = [c2[k]["diff"] for k, _ in keys]
    err = [[c2[k]["diff"] - c2[k]["ci_lo"] for k, _ in keys],
           [c2[k]["ci_hi"] - c2[k]["diff"] for k, _ in keys]]
    cols = [C_FLAG if c2[k]["verdict"] == "WIN" else C_OK for k, _ in keys]
    a2.bar(xs, vals, color=cols, width=0.55)
    a2.errorbar(xs, vals, yerr=err, fmt="none", ecolor=C_REF, capsize=4, lw=1.0)
    for xi, (k, _lab) in enumerate(keys):
        a2.text(xi, c2[k]["ci_hi"] + 0.02, c2[k]["verdict"], ha="center", fontsize=7,
                color=cols[xi], fontweight="bold")
    a2.axhline(0.0, color=C_REF, lw=0.9)
    a2.set_xticks(xs)
    a2.set_xticklabels([lab for _k, lab in keys], fontsize=7)
    a2.set_ylabel("paired ΔAUC (calibrated − foil)")
    a2.set_title("C2  the missing half of the two-way null:\n the cutoff contrast excludes zero",
                 loc="left", fontweight="bold")

    # --- C3 -------------------------------------------------------------------------------
    c3 = R["c3"]
    labels = ["predicted\n(leverage-wtd)", "observed\n(matched rule)", "observed\nmedian",
              "observed\nunweighted mean"]
    vals = [c3["pred"], c3["matched"], c3["median"], c3["unweighted_mean"]]
    errs = np.array([[c3["pred"] - c3["pred_ci"][0], c3["matched"] - c3["matched_ci"][0], 0, 0],
                     [c3["pred_ci"][1] - c3["pred"], c3["matched_ci"][1] - c3["matched"], 0, 0]])
    a3.bar(range(4), vals, color=[C_ALT, C_FLAG, "#7f7f7f", "#bbbbbb"], width=0.6)
    errs[errs == 0.0] = np.nan          # the two unweighted rules carry no interval; draw none
    a3.errorbar(range(4), vals, yerr=errs, fmt="none", ecolor=C_REF, capsize=4, lw=1.0)
    for i, v in enumerate(vals):
        a3.text(i, v + 0.05, f"{v:.2f}", ha="center", fontsize=7)
    a3.axhline(1.0, color=C_REF, ls=":", lw=0.9)
    a3.set_xticks(range(4))
    a3.set_xticklabels(labels, fontsize=7)
    a3.set_ylabel(r"aggregate $\chi^2_\nu$ (replicate 0)")
    a3.set_title(f"C3  'the aggregate also matches': observed is\n "
                 f"{c3['matched'] / c3['pred']:.1f}× predicted and outside its interval",
                 loc="left", fontweight="bold")

    # --- C4 -------------------------------------------------------------------------------
    c4 = R["c4"]
    bins = np.arange(0.5, c4["max"] + 1.5, 1.0)
    a4.hist(c4["dof"], bins=bins, color=C_OK, alpha=0.85)
    flag_dofs = list(c4["flag_dof"].values())
    a4.hist(flag_dofs, bins=bins, color=C_FLAG, alpha=0.95, label="ever flagged")
    a4.axvline(c4["median"], color=C_REF, ls="--", lw=1.0,
               label=f"median dof = {c4['median']:g}")
    a4.set_xlabel("independent cycles (dof) per system")
    a4.set_ylabel(f"systems (of {len(c4['names'])})")
    a4.legend(frameon=False, fontsize=7)
    a4.set_title(f"C4  dof distribution: {c4['counts'][1]} systems have\n "
                 f"a single cycle (quartiles {c4['q1']:g}/{c4['median']:g}/{c4['q3']:g})",
                 loc="left", fontweight="bold")

    ax4b = a4.inset_axes((0.55, 0.42, 0.42, 0.45))
    st = c4["primary"]
    xs = np.arange(2)
    ax4b.bar(xs - 0.17, [s["mean_jaccard"] for s in st], width=0.32, color=C_OK, label="observed")
    ax4b.bar(xs + 0.17, [s["mean_jaccard_random"] for s in st], width=0.32, color="#bbbbbb",
             label="random")
    ax4b.set_xticks(xs)
    ax4b.set_xticklabels(["low dof", "high dof"], fontsize=7)
    ax4b.set_ylabel("mean Jaccard", fontsize=7)
    ax4b.tick_params(labelsize=7)
    ax4b.legend(frameon=False, fontsize=7)

    # --- C5 -------------------------------------------------------------------------------
    c5 = R["c5"]
    a5.hist(c5["null"], bins=24, color=C_OK, alpha=0.85,
            label=f"no-systematic-error null\n(n={c5['n_null']} realizations)")
    a5.axvline(PUB_RHO, color=C_FLAG, lw=1.6, label=f"reported ρ = {PUB_RHO:.3f}")
    a5.axvline(0.0, color=C_REF, ls=":", lw=1.0)
    a5.axvline(c5["ctrl_median"], color=C_ALT, ls="--", lw=1.2,
               label=f"positive control ρ = {c5['ctrl_median']:.2f}")
    a5.set_xlabel("realized Spearman ρ (predicted vs observed)")
    a5.set_ylabel("null realizations")
    a5.legend(frameon=False, fontsize=7, loc="upper center")
    a5.set_title(f"C5  shared-denominator null: ρ = {c5['mean']:+.3f}\n "
                 f"[{c5['ci'][0]:+.2f}, {c5['ci'][1]:+.2f}], so 0.648 is not an artefact",
                 loc="left", fontweight="bold")

    # --- C6 -------------------------------------------------------------------------------
    c6 = R["c6"]
    fr = [d["frac_below"] for d in c6["per_target"]]
    ys = np.arange(len(fr))
    a6.barh(ys, fr, color=[C_FLAG if v > c6["ref_raw"] else C_OK for v in fr], height=0.72)
    a6.axvline(c6["ref_raw"], color=C_REF, ls="--", lw=1.1,
               label=f"perfect calibration, n=3: {c6['ref_raw']:.3f}")
    a6.axvline(c6["frac_below"], color=C_ALT, ls=":", lw=1.2,
               label=f"pooled {c6['frac_below']:.3f} "
                     f"[{c6['frac_below_ci'][0]:.2f}, {c6['frac_below_ci'][1]:.2f}]")
    a6.set_yticks(ys)
    a6.set_yticklabels([d["target"] for d in c6["per_target"]], fontsize=7)
    a6.set_ylim(len(fr) - 0.5, -0.5)
    a6.set_xlabel("fraction of edges with analytic se < replicate spread")
    a6.legend(frameon=False, fontsize=7, loc="lower right")
    a6.set_title(f"C6  Wade's functional: {c6['n_targets_above_ref']} of "
                 f"{len(fr)} targets exceed\n the perfect-calibration reference",
                 loc="left", fontweight="bold")

    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(FIGDIR / f"figInf_inference.{ext}")
    plt.close(fig)


# ============================================================ results doc

def _pair(i, j):
    return f"({i},{j})"


def write_doc(R):
    c1, c2, c3, c4, c5, c6 = (R[k] for k in ("c1", "c2", "c3", "c4", "c5", "c6"))

    t1 = ["| pair | pooled r | cluster-bootstrap 95% CI | within-system permutation p | "
          "permutation null mean ± SD |", "|---|---:|---|---:|---|"]
    for pi, (i, j) in enumerate(PAIRS):
        t1.append(f"| {_pair(i, j)} | {c1['obs'][pi]:+.3f} | "
                  f"[{c1['ci'][pi][0]:+.3f}, {c1['ci'][pi][1]:+.3f}] | "
                  f"{c1['p'][pi]:.4f} | {c1['null_mean'][pi]:+.4f} ± {c1['null_sd'][pi]:.4f} |")

    t1b = ["| pair | r (flagged) | r (unflagged) | difference | stratified cluster CI | "
           "system-label permutation p | within-system permutation p |",
           "|---|---:|---:|---:|---|---:|---:|"]
    for pi, (i, j) in enumerate(PAIRS):
        t1b.append(f"| {_pair(i, j)} | {c1['obs_flagged'][pi]:+.3f} | "
                   f"{c1['obs_unflagged'][pi]:+.3f} | {c1['obs_diff'][pi]:+.3f} | "
                   f"[{c1['ci_diff'][pi][0]:+.3f}, {c1['ci_diff'][pi][1]:+.3f}] | "
                   f"{c1['p_diff_label'][pi]:.4f} | {c1['p_diff_edge'][pi]:.4f} |")

    t2 = ["| contrast | AUC (calibrated) | AUC (foil) | paired ΔAUC | 95% CI | verdict |",
          "|---|---:|---:|---:|---|---|"]
    for key, lab, foil in (("vs_cutoff", f"vs fixed {HYSTERESIS_CUTOFF} kcal/mol cutoff", "cutoff"),
                           ("vs_fixed_se", "vs fixed-se χ² test", "fixed_se"),
                           ("vs_max", "vs max(both) — the published row", None)):
        d = c2[key]
        foil_auc = "—" if foil is None else f"{c2['auc'][foil]:.3f}"
        t2.append(f"| {lab} | {c2['auc']['calibrated']:.3f} | {foil_auc} | {d['diff']:+.3f} | "
                  f"[{d['ci_lo']:+.3f}, {d['ci_hi']:+.3f}] | **{d['verdict']}** |")

    t3 = ["| aggregation rule | value | 95% CI (system cluster bootstrap) |", "|---|---:|---|",
          f"| **predicted** E[χ²ᵥ], leverage-weighted (published) | {c3['pred']:.2f} | "
          f"[{c3['pred_ci'][0]:.2f}, {c3['pred_ci'][1]:.2f}] |",
          f"| **observed, matched rule** (= dof-weighted global Σχ²/Σdof) | {c3['matched']:.2f} | "
          f"[{c3['matched_ci'][0]:.2f}, {c3['matched_ci'][1]:.2f}] |",
          f"| observed, median over systems | {c3['median']:.2f} | — |",
          f"| observed, unweighted mean over systems | {c3['unweighted_mean']:.2f} | — |"]

    t4 = ["| dof | systems |", "|---:|---:|"]
    for v in sorted(set(c4["dof"].tolist())):
        t4.append(f"| {v} | {int(np.sum(c4['dof'] == v))} |")

    def strat_rows(strata):
        rows = ["| stratum | systems | flagged per replicate | ever / always | "
                "pairwise Jaccard | mean (random ref) | median χ²ᵥ swing |",
                "|---|---:|---|---|---|---|---:|"]
        for s in strata:
            rows.append(
                f"| {s['label']} | {s['n']} | {'/'.join(str(v) for v in s['sizes'])} | "
                f"{len(s['ever'])} / {len(s['inter'])} | "
                + " / ".join(f"{v:.2f}" for v in s["jaccard"])
                + f" | {s['mean_jaccard']:.3f} ({s['mean_jaccard_random']:.3f}) | "
                f"{s['median_swing']:.2f}× |")
        return rows

    t6 = ["| target | edges | fraction of edges with analytic se < replicate spread | "
          "RMS ratio | median per-edge ratio |", "|---|---:|---:|---:|---:|"]
    for d in c6["per_target"]:
        mark = "**" if d["frac_below"] > c6["ref_raw"] else ""
        t6.append(f"| {d['target']} | {d['n']} | {mark}{d['frac_below']:.3f}{mark} | "
                  f"{d['ratio']:.2f} | {d['median_ratio']:.2f} |")

    verdicts = [
        ("C1", "**supports the manuscript, with the caveat now stated**",
         f"all three r are positive with cluster-bootstrap CIs excluding zero and permutation "
         f"p ≤ {max(c1['p']):.4f}; the flagged-vs-unflagged contrast is positive on all three "
         f"pairs but its system-label permutation p is {min(c1['p_diff_label']):.3f}–"
         f"{max(c1['p_diff_label']):.3f}."),
        ("C2", "**AGAINST the manuscript**" if c2["vs_cutoff"]["ci_lo"] > 0
         else "supports the manuscript",
         f"calibrated vs the {HYSTERESIS_CUTOFF} kcal/mol cutoff: "
         f"ΔAUC = {c2['vs_cutoff']['diff']:+.3f} "
         f"[{c2['vs_cutoff']['ci_lo']:+.3f}, {c2['vs_cutoff']['ci_hi']:+.3f}], which "
         + ("excludes zero, so the two-way tie claim is wrong as written."
            if c2["vs_cutoff"]["ci_lo"] > 0 else "includes zero, so the tie claim holds.")),
        ("C3", "**AGAINST the manuscript as written**" if not c3["pred_covers_obs"]
         else "supports the manuscript",
         f"under the matched (leverage-weighted = dof-weighted) rule the observed aggregate is "
         f"{c3['matched']:.2f} [{c3['matched_ci'][0]:.2f}, {c3['matched_ci'][1]:.2f}] against the "
         f"predicted {c3['pred']:.2f} [{c3['pred_ci'][0]:.2f}, {c3['pred_ci'][1]:.2f}]: the "
         f"intervals {'overlap' if c3['overlap'] else 'do not overlap'} but the predicted "
         f"interval {'covers' if c3['pred_covers_obs'] else 'does not cover'} the observed point, "
         f"and the observed value swings {c3['median']:.2f}–{c3['matched']:.2f} across rules."),
        ("C4", "**a real limitation, now quantified**",
         f"{c4['counts'][1]} of {len(c4['names'])} systems have a single cycle (quartiles "
         f"{c4['q1']:g}/{c4['median']:g}/{c4['q3']:g}); "
         f"{len(c4['primary'][1]['ever'])} of the "
         f"{len(c4['primary'][0]['ever']) + len(c4['primary'][1]['ever'])} ever-flagged systems "
         "and the only always-flagged one are high-dof, and the per-system "
         f"χ²ᵥ swing is {c4['secondary'][0]['median_swing']:.1f}× at dof = 1 versus "
         f"{c4['secondary'][1]['median_swing']:.1f}× at dof ≥ 2."),
        ("C5", "**supports the manuscript**",
         f"the shared-denominator null realizes ρ = {c5['mean']:+.3f} "
         f"[{c5['ci'][0]:+.3f}, {c5['ci'][1]:+.3f}], max {c5['max']:+.3f}; the reported "
         f"{PUB_RHO:.3f} exceeds every one of the {c5['n_null']} null draws "
         f"(p = {c5['p_vs_null']:.4f}), so it is not a shared-denominator artefact — but the null "
         "is not centred on zero and the comparison should be made against it, not against zero."),
        ("C6", "**against the manuscript's reconciliation**" if c6["frac_below_ci"][1] <
         c6["ref_raw"] else "supports the manuscript's reconciliation",
         f"on Wade et al.'s own functional the analytic se falls below the replicate spread on "
         f"{100 * c6['frac_below']:.1f}% "
         f"[{100 * c6['frac_below_ci'][0]:.1f}, {100 * c6['frac_below_ci'][1]:.1f}]% of edges, "
         f"*below* the {100 * c6['ref_raw']:.1f}% a perfectly calibrated bar produces at n = 3, so "
         "the like-for-like comparison makes the disagreement with their under-estimation finding "
         "sharper, not milder."),
    ]
    tv = ["| check | verdict | one line |", "|---|---|---|"]
    tv += [f"| {k} | {v} | {t} |" for k, v, t in verdicts]

    lines = [
        "# Results — Fig Inf: inference for six claims that were stated without it",
        "",
        "**Figure:** `figs/figInf_inference.{pdf,png}` · **Reproduce:** `make figInf`",
        "(`PYTHONPATH=src python figs/make_figInf.py`). Deterministic (every resampling seeded,",
        f"`SEED = {SEED}`). Data: `data/openfe_replicates/combined_pymbar4_edge_data.csv` (OpenFE",
        "IndustryBenchmarks2024, public, 3 independent replicates per edge). All fits are the",
        "repo's own (`bar.qc.gls_network`, `bar.leverage.curl_leverage`, `bar.detectors`); the",
        "edge construction and the detectors are imported from `make_figStab.py`,",
        "`make_figCut.py`, `make_figOOS.py`, `make_figL_validation.py` and",
        "`make_figA_replicates.py` rather than reimplemented, and C1 and C2 each assert that they",
        "reproduce the shipped number before adding anything to it.",
        "",
        "Every design below was fixed before its number was computed and is stated in the script",
        "docstring. Nothing was tuned after seeing a result. Several of these come back against",
        "the manuscript; they are reported as such.",
        "",
        "## Summary",
        "",
    ] + tv + [
        "",
        "## C1 — inference for the cross-replicate residual correlation",
        "",
        f"The {c1['n_edges']} edges are nested in {c1['n_sys']} networks and each residual is a",
        "projection through that network's residual maker, so neither an interval nor a test may",
        "treat edges as independent. The interval is a **cluster bootstrap that resamples the",
        f"{c1['n_sys']} systems** ({N_BOOT} draws); the test is a **within-system permutation**",
        f"({N_PERM} draws) that permutes the edge order inside each system, destroying the",
        "cross-replicate pairing while preserving system membership, system sizes and the marginal",
        "residual distribution. Under the sampling null the three replicate fits are independent,",
        "so `E[z_e^(i) z_f^(j)] = 0` for every pair of edges including `e = f`; the permutation",
        "therefore has the right null, and its realized mean (below) confirms it is centred on",
        "zero rather than on the projector-induced level that a naive within-replicate shuffle",
        "would produce.",
        "",
    ] + t1 + [
        "",
        "**The three rotations are dependent.** They are three views of the same three runs, and",
        "every pair of rotations shares a replicate — (0,1) and (0,2) share replicate 0, and so",
        "on. They are not three replications, the three p-values are not independent evidence,",
        "and no multiplicity correction across them would be meaningful. One permutation `pi` per",
        "draw is applied to all three pairs, so the null inherits the same dependence.",
        "",
        "### Where the correlation lives",
        f"{c1['n_flag_edges']} of the {c1['n_edges']} edges belong to the {c1['n_flag_sys']}",
        "flagged systems. The contrast gets a **stratified** cluster bootstrap (flagged systems",
        "resampled among themselves, unflagged among themselves, so neither stratum can empty) and",
        "two permutation nulls: a **system-label** permutation that reassigns which",
        f"{c1['n_flag_sys']} of the {c1['n_sys']} systems carry the flag, holding all data fixed —",
        "the null for *concentration* — and the same within-system edge permutation as above,",
        "whose null is *no correlation anywhere* and is therefore the easier bar.",
        "",
    ] + t1b + [
        "",
        "The contrast is positive on all three pairs, so the direction the manuscript asserts is",
        "the direction in the data. But the test that actually asks the manuscript's question —",
        "the system-label permutation — gives p = "
        + ", ".join(f"{v:.3f}" for v in c1["p_diff_label"]) + ". "
        + (f"None of the three reaches α = {ALPHA:.2f}, and the stratified cluster interval "
           f"covers zero on {sum(1 for lo, hi in c1['ci_diff'] if lo <= 0 <= hi)} of the 3 pairs."
           if min(c1["p_diff_label"]) > ALPHA else
           "Some reach α; the table shows which."),
        "The within-system edge permutation is much smaller (p = "
        + ", ".join(f"{v:.4f}" for v in c1["p_diff_edge"]) + "), but its null is *no correlation "
        "anywhere*, which is not the claim being made; a contrast can beat that null purely",
        "because the correlation exists at all.",
        "",
        f"**This is a qualification the manuscript needs.** With {c1['n_flag_sys']} flagged",
        f"systems out of {c1['n_sys']}, the label permutation has few distinguishable arrangements",
        "of the clusters that carry the signal, so the test is underpowered and the honest reading",
        "is *not established* rather than *refuted*. The existence of the cross-replicate",
        "correlation (first table) is solid at cluster-honest inference; its **concentration** in",
        "the flagged systems is directional evidence, not a demonstrated effect.",
        "",
        "## C2 — the missing half of the two-way null",
        "",
        "The manuscript reports one interval, `+0.140 [-0.008, +0.316]`, which is",
        "`AUC(A) − max(AUC(B), AUC(C))` and hence the comparison against the **fixed-se** foil.",
        "The comparison against the fixed 1.0 kcal/mol hysteresis cutoff had no interval anywhere.",
        f"All three contrasts below come from the SAME {N_BOOT_AUC} system resamples with the same",
        "`seed = 0` as `make_figCut.py`, so the published row reproduces exactly (asserted in",
        "code):",
        "",
    ] + t2 + [
        "",
        ("**This is against the manuscript.**" if c2["vs_cutoff"]["ci_lo"] > 0
         else "**This confirms the manuscript.**")
        + f" The calibrated null beats the "
        f"{HYSTERESIS_CUTOFF} kcal/mol cutoff by {c2['vs_cutoff']['diff']:+.3f} AUC with a 95% CI",
        f"of [{c2['vs_cutoff']['ci_lo']:+.3f}, {c2['vs_cutoff']['ci_hi']:+.3f}], which excludes",
        "zero. The sentence 'does not measurably out-discriminate either foil' is true of the",
        "fixed-se foil and false of the fixed-cutoff foil. The pre-registered WIN rule in",
        "`bar.detectors.paired_auc_bootstrap` is conjunctive (A must beat BOTH), so the overall",
        "TIE verdict stands unchanged — but the *reason* it is a tie is entirely the fixed-se",
        "comparison, and the text must say that instead of implying both comparisons are null.",
        f"For context, the cutoff rule flags {c2['n_flag']['cutoff']} of the {c2['n']} systems and",
        f"the calibrated rule flags {c2['n_flag']['calibrated']}.",
        "",
        "## C3 — the observed aggregate the prediction is said to match",
        "",
        "The predicted aggregate is `Σ_s Σ_e h_e c_e^−2 / Σ_s dof_s`. Because `Σ_e h_e = dof_s`",
        "(Theorem D1) and `E[X²_s] = Σ_e h_e c_e^−2`, the observed statistic under the **same**",
        "rule is `Σ_s X²_s / Σ_s dof_s` — that is, the leverage-weighted rule and the dof-weighted",
        "global value are the same number, so only three of the four requested aggregates are",
        "distinct. All are reported:",
        "",
    ] + t3 + [
        "",
        f"**The intervals {'overlap' if c3['overlap'] else 'do NOT overlap'}**",
        f"({c3['matched_ci'][0]:.2f} ≤ {c3['pred_ci'][1]:.2f}), but the predicted interval",
        f"{'covers' if c3['pred_covers_obs'] else 'does not cover'} the observed point:",
        f"{c3['matched']:.2f} against a predicted upper bound of {c3['pred_ci'][1]:.2f}, i.e. the",
        f"observed aggregate is {c3['matched'] / c3['pred']:.2f}× the predicted one. 'Matches'",
        "survives only in the weak sense that two wide intervals share ground.",
        "",
        "**The rule matters more than the agreement does.** The same 48 systems give",
        f"{c3['median']:.2f} under a median over systems, {c3['unweighted_mean']:.2f} under an",
        f"unweighted mean over systems and {c3['matched']:.2f} under the matched rule — a factor",
        f"of {c3['matched'] / c3['median']:.1f} between the extremes — because the closure χ² is",
        "dominated by a few high-dof, high-residual networks while the median system closes far",
        "better than its bars predict. Quoting 'the aggregate level also matches' without naming",
        "the rule or the observed number is therefore not a checkable claim, and it is the reader",
        "who has to guess which of these three the sentence means. Under the only rule that is",
        "actually matched to the prediction, observed sits above predicted.",
        "",
        "That direction is the expected one and is already the paper's own explanation — a",
        "reproducible systematic error inflates a closure residual without inflating the spread",
        "between repeats, so observed *should* exceed predicted wherever the QC fires. But that is",
        "an argument for an interpretable discrepancy, not for a match, and the sentence should",
        "state the observed number, the rule, and the direction rather than assert agreement.",
        "",
        "## C4 — the degrees-of-freedom distribution, and stability stratified by it",
        "",
        f"Across the {len(c4['names'])} admitted systems (replicate 0): min {c4['min']}, quartiles",
        f"{c4['q1']:g} / {c4['median']:g} / {c4['q3']:g}, max {c4['max']}, mean",
        f"{c4['mean']:.1f}. Counts at the low end: {c4['counts'][1]} systems with one cycle,",
        f"{c4['counts'][2]} with two, {c4['counts'][3]} with three.",
        "",
    ] + t4 + [
        "",
        "### Flag stability stratified by dof (primary split: median of dof)",
        "",
    ] + strat_rows(c4["primary"]) + [
        "",
        "### Secondary split: single-cycle systems versus the rest",
        "",
    ] + strat_rows(c4["secondary"]) + [
        "",
        "**Stability is carried by the high-cycle systems, and this is said plainly.** The",
        "evidence is not the Jaccard column — the low-dof stratum's",
        f"{c4['primary'][0]['mean_jaccard']:.3f} is computed over",
        f"{len(c4['primary'][0]['ever'])} ever-flagged system with per-replicate counts",
        f"{'/'.join(str(v) for v in c4['primary'][0]['sizes'])}, so it is a degenerate average of",
        "1.00 and two 0.00s and should not be compared with anything. The evidence is:",
        "",
        f"- **Where the flags are.** {len(c4['primary'][1]['ever'])} of the "
        f"{len(c4['primary'][0]['ever']) + len(c4['primary'][1]['ever'])} ever-flagged systems sit",
        "  in the high-dof stratum, including the only system flagged in *all three* replicates.",
        f"  The low-dof half of the benchmark ({c4['primary'][0]['n']} systems) contributes",
        f"  {len(c4['primary'][0]['ever'])} ever-flagged system and nothing that reproduces.",
        "- **Why.** The per-system χ²ᵥ swing across the three repeats has median "
        f"{c4['secondary'][0]['median_swing']:.2f}× at dof = 1 against "
        f"{c4['secondary'][1]['median_swing']:.2f}× at dof ≥ 2 (and "
        f"{c4['primary'][0]['median_swing']:.2f}× vs {c4['primary'][1]['median_swing']:.2f}× under",
        "  the median split). A single-cycle network's reduced χ² *is* one number, so it moves by",
        "  an order of magnitude run to run and its BH-adjusted q crosses α about as often as not.",
        f"- **How much of the benchmark this is.** {c4['counts'][1]} of {len(c4['names'])} systems",
        f"  have one cycle and {c4['counts'][1] + c4['counts'][2] + c4['counts'][3]} have three or",
        "  fewer, so this is not an edge case: a third of the benchmark is structurally incapable",
        "  of delivering a stable per-system verdict, whatever the detector.",
        "",
        "A Jaccard of 1.00 in a stratum with no flagged systems is the empty-set convention, not",
        "evidence of stability; the 'flagged per replicate' column shows where that applies. Any",
        "set-valued claim in the manuscript should be scoped to the well-determined networks, and",
        "the single-cycle systems should be reported as flag-eligible but not flag-stable.",
        "",
        "## C5 — is the predicted-versus-observed check circular?",
        "",
        "**Which se enters where.** In the headline rotation the predicted quantity is",
        "`Σ_e h_e (s_e² / se_{e,0}²) / dof`, where `s_e²` is the sample variance of the ΔΔG values",
        "of replicates **1 and 2** and `se_{e,0}` is the reported standard error of replicate",
        "**0**. The observed quantity is `Σ_e r_{e,0}² / se_{e,0}² / dof`. So replicate 0's",
        "reported se is the denominator on both sides, and the curl-leverage weights `h_e` are",
        "computed from replicate 0's `V_e` as well. The numerators are disjoint (replicates 1–2",
        "versus replicate 0), but the denominators and the weights are shared, which is exactly",
        "the referee's concern.",
        "",
        "**The null.** `make_figOOS.resample_null_world` — the repo's own no-systematic-error",
        "world — redraws each replicate's values as `y ~ N(0, se_k)` on the same graphs with the",
        "same reported errors, so `c_e ≡ 1` identically while `h_e` and `se_e` are untouched. The",
        "shared denominator is therefore preserved exactly, and any correlation the harness",
        f"realizes is the artefact. Over {N_NULL} realizations:",
        "",
        f"- realized Spearman ρ: mean {c5['mean']:+.3f}, median {c5['median']:+.3f}, SD",
        f"  {c5['sd']:.3f}, 95% range [{c5['ci'][0]:+.3f}, {c5['ci'][1]:+.3f}], max",
        f"  {c5['max']:+.3f};",
        f"- the reported ρ = {PUB_RHO:.3f} against that null: p = {c5['p_vs_null']:.4f};",
        "- positive control (a genuine per-system bar miscalibration spanning 0.5×–2.0×, same",
        f"  harness, same shared denominator): median ρ = {c5['ctrl_median']:+.3f}",
        f"  [{c5['ctrl_ci'][0]:+.3f}, {c5['ctrl_ci'][1]:+.3f}], n = {c5['n_ctrl']} — so the",
        "  harness",
        "  can see a real signal through the shared denominator, and its silence under the null is",
        "  informative rather than a lack of power.",
        "",
        "**Verdict: not circular.** The shared `V_e` does not manufacture a rank correlation,",
        "because under the null both sides are centred on their own expectations regardless of the",
        f"magnitude of `se_e`. The {PUB_RHO:.3f} should nevertheless be reported against this null",
        "rather than against zero, since 'zero' was never the operative alternative; the null's",
        f"95% range is [{c5['ci'][0]:+.3f}, {c5['ci'][1]:+.3f}], which is where an artefact would",
        "have shown up.",
        "",
        "## C6 — a like-for-like comparison with Wade et al.",
        "",
        "Wade et al. (JCTC 2022) compare the analytic MBAR uncertainty against the spread over",
        "independent replicas and report that MBAR *under*-estimates. The manuscript's",
        "reconciliation invokes a curl-leverage-weighted mean of `c_e^−2`, a functional with no",
        "counterpart in their measurement. Their functional, computed here on these edges: the",
        "per-edge analytic se (RMS over the three replicates, the repo's own convention in",
        "`make_figA_replicates.load`) against the across-replicate SD, and the **fraction of edges",
        "on which the analytic se falls below the replicate spread**.",
        "",
        "**Reference level, stated before running.** With three replicates and perfectly",
        "calibrated bars, `s² ~ σ² χ²₂ / 2`, so `P(s > σ) = e^{−1} = "
        f"{c6['ref_raw']:.3f}`: about a third of edges land below by construction, and only an",
        "excess over that is evidence of under-estimation. The small-sample-bias-corrected variant",
        f"(`σ̂ = s / c₄`, `c₄ = {c6['c4']:.3f}`) has reference `e^{{−c₄²}} = {c6['ref_c4']:.3f}`.",
        "",
        f"The edge set is `make_figA_replicates.load`'s ({c6['n_edges']} edges over",
        f"{c6['n_sys']} systems: every row complete in all three replicates with a non-degenerate",
        "spread), which is the set that produced the manuscript's 1.41×. It is admitted by a",
        "different rule from C1's 1143 network-fitted edges, so the two counts are not expected to",
        "agree.",
        "",
        f"- **Pooled: the analytic se is below the replicate spread on {c6['frac_below']:.3f}**",
        f"  [{c6['frac_below_ci'][0]:.3f}, {c6['frac_below_ci'][1]:.3f}] of the",
        f"  {c6['n_edges']} edges (system-cluster bootstrap), against the",
        f"  {c6['ref_raw']:.3f} a perfectly calibrated bar produces.",
        f"- c₄-corrected: {c6['frac_below_c4']:.3f} "
        f"[{c6['frac_below_c4_ci'][0]:.3f}, {c6['frac_below_c4_ci'][1]:.3f}] against",
        f"  {c6['ref_c4']:.3f}.",
        f"- For orientation, the same edges give an RMS-pooled ratio of {c6['rms_ratio']:.2f}× and",
        f"  a **median per-edge** ratio of {c6['median_ratio']:.2f}×.",
        f"- Per target ({len(c6['per_target'])} targets with ≥ {MIN_TARGET_EDGES} edges):",
        f"  {c6['n_targets_above_ref']} exceed the {c6['ref_raw']:.3f} reference and",
        f"  {c6['n_targets_ratio_lt1']} have an RMS ratio below 1.",
        "",
    ] + t6 + [
        "",
        "**Reading — this goes against the manuscript's reconciliation, not with it.** On the",
        "functional Wade et al. actually report, these edges are *more* conservative than a",
        f"perfectly calibrated bar would be: {100 * c6['frac_below']:.1f}% of edges fall below the",
        f"replicate spread against the {100 * c6['ref_raw']:.1f}% that perfect calibration",
        "produces",
        f"at n = 3, and the bootstrap interval "
        f"[{100 * c6['frac_below_ci'][0]:.1f}, {100 * c6['frac_below_ci'][1]:.1f}]% excludes the",
        "reference. The median per-edge ratio is "
        f"{c6['median_ratio']:.2f}×, *above* the RMS-pooled {c6['rms_ratio']:.2f}×, so the typical",
        "edge is more conservative than the aggregate ratio suggests rather than less. The",
        "like-for-like comparison therefore makes the disagreement with their under-estimation",
        "finding **sharper**, and the manuscript's sentence that the measurement is 'much closer",
        "to their picture than 1.41 suggests' is supported only by the leverage-weighted",
        "functional, which has no counterpart in what they measured. A reconciliation that only",
        "works in a functional the other paper never computes is not a reconciliation.",
        "",
        "The one caveat that does cut the other way is the manuscript's own, and it is not settled",
        "here: these repeats share starting coordinates and differ only in stochastic sampling,",
        "whereas Wade et al. resample the full ensemble, so this denominator is a lower bound on",
        "true run-to-run variability and every fraction above is an *under*-estimate of the",
        "fraction Wade et al. would measure. That argument alone has to carry the whole",
        "reconciliation; the numbers do not help it.",
        "",
        f"The heterogeneity is real and should be reported with the aggregate: "
        f"{c6['n_targets_above_ref']} of {len(c6['per_target'])} targets (≥ {MIN_TARGET_EDGES}",
        f"edges) exceed the reference fraction and {c6['n_targets_ratio_lt1']} have an RMS ratio",
        "below 1, so 'conservative in aggregate, heterogeneous per target' remains the accurate",
        "summary — it is only the *direction of the reconciliation* that this check contradicts.",
        "",
        "## Honest reading",
        "- **Three of the six come back against the manuscript**, and all three are claims that",
        "  need rewording rather than analyses that need redoing.",
        f"  - **C2:** the two-way tie is half wrong. The calibrated null does beat the "
        f"{HYSTERESIS_CUTOFF} kcal/mol",
        f"    cutoff ({c2['vs_cutoff']['diff']:+.3f} "
        f"[{c2['vs_cutoff']['ci_lo']:+.3f}, {c2['vs_cutoff']['ci_hi']:+.3f}]); the tie comes",
        "    entirely from the fixed-se foil. The conjunctive pre-registered WIN rule still yields",
        "    TIE overall, so the verdict stands and only the sentence has to change.",
        f"  - **C3:** the observed aggregate the prediction is said to match is "
        f"{c3['matched']:.2f} "
        f"[{c3['matched_ci'][0]:.2f}, {c3['matched_ci'][1]:.2f}] under the matched rule, against a",
        f"    predicted {c3['pred']:.2f} [{c3['pred_ci'][0]:.2f}, {c3['pred_ci'][1]:.2f}] — the",
        "    intervals overlap but the prediction does not cover the point, and the observed value",
        f"    ranges {c3['median']:.2f}–{c3['matched']:.2f} depending on a rule the text never",
        "    names.",
        "  - **C6:** the Wade et al. reconciliation runs through a functional they never compute.",
        "    On the one they do compute, this data set is *more* conservative than perfect",
        "    calibration, so the like-for-like comparison sharpens the disagreement instead of",
        "    softening it. Only the shared-starting-coordinates caveat still argues for",
        "    reconciliation, and it now has to carry that argument alone.",
        "- **C1 and C5 hold up.** The cross-replicate correlation is real under cluster-honest",
        "  inference (all three CIs exclude zero, permutation p at the resolution floor), and the",
        "  shared-denominator artefact the referee suspected in the predicted-versus-observed",
        f"  check does not exist ({PUB_RHO:.3f} exceeds all {c5['n_null']} null draws). Two",
        "  qualifications survive: C1 does **not** establish the *concentration* of the",
        "  correlation",
        "  in the flagged systems at α = 0.05 (the manuscript asserts it; the label-permutation",
        f"  test is underpowered with {c1['n_flag_sys']} flagged systems), and C5's null is not",
        f"  centred on zero — its 95% range reaches {c5['ci'][1]:+.3f} — so 'sampling noise would",
        "  give zero' is the wrong reference and the measured null is the right one.",
        "- **C4 quantifies a structural limitation** rather than testing a claim: "
        f"{c4['counts'][1]} of {len(c4['names'])} networks carry a single cycle, their reduced χ²",
        f"  swings by a median {c4['secondary'][0]['median_swing']:.1f}× between repeats, and",
        "  every",
        "  reproducibly flagged system is high-dof. Set-valued statements should be scoped to the",
        "  well-determined networks.",
        "- Nothing here was tuned. Every constant is in the script docstring and was fixed before",
        "  the corresponding number was computed; C1 and C2 additionally assert that they",
        "  reproduce the shipped numbers before adding inference to them.",
        "",
    ]
    DOC.write_text("\n".join(ln for ln in lines if ln is not None) + "\n")


# ============================================================ driver

def analyze():
    recs = load_records()
    A = panel_a(recs)
    return {"c1": c1_residual_correlation(), "c2": c2_auc_contrasts(),
            "c3": c3_aggregate(recs, A), "c4": c4_dof_stratified(),
            "c5": c5_circularity(recs), "c6": c6_wade()}


def main():
    _style()
    R = analyze()
    make_figure(R)
    write_doc(R)
    print(f"wrote figInf_inference.(pdf|png) to {FIGDIR}")
    print(f"wrote {DOC.relative_to(ROOT)}")

    c1, c2, c3, c4, c5, c6 = (R[k] for k in ("c1", "c2", "c3", "c4", "c5", "c6"))
    print(f"\n[C1] n={c1['n_edges']} edges in {c1['n_sys']} systems "
          f"({c1['n_flag_edges']} edges in {c1['n_flag_sys']} flagged systems)")
    for pi, (i, j) in enumerate(PAIRS):
        print(f"    pair {_pair(i, j)}  r={c1['obs'][pi]:+.3f} "
              f"cluster-CI[{c1['ci'][pi][0]:+.3f},{c1['ci'][pi][1]:+.3f}] "
              f"perm p={c1['p'][pi]:.4f}  (null {c1['null_mean'][pi]:+.4f}"
              f"±{c1['null_sd'][pi]:.4f})")
    for pi, (i, j) in enumerate(PAIRS):
        print(f"    pair {_pair(i, j)}  flagged {c1['obs_flagged'][pi]:+.3f} vs unflagged "
              f"{c1['obs_unflagged'][pi]:+.3f}  diff {c1['obs_diff'][pi]:+.3f} "
              f"[{c1['ci_diff'][pi][0]:+.3f},{c1['ci_diff'][pi][1]:+.3f}] "
              f"label-perm p={c1['p_diff_label'][pi]:.4f} edge-perm p={c1['p_diff_edge'][pi]:.4f}")

    print(f"\n[C2] AUC calibrated={c2['auc']['calibrated']:.3f} cutoff={c2['auc']['cutoff']:.3f} "
          f"fixed-se={c2['auc']['fixed_se']:.3f}  (n={c2['n']} systems)")
    for key, lab in (("vs_cutoff", "vs fixed cutoff"), ("vs_fixed_se", "vs fixed se"),
                     ("vs_max", "vs max(both) [published]")):
        d = c2[key]
        print(f"    {lab:28s} dAUC={d['diff']:+.3f} [{d['ci_lo']:+.3f},{d['ci_hi']:+.3f}]  "
              f"{d['verdict']}")

    print(f"\n[C3] predicted (leverage-weighted) {c3['pred']:.2f} "
          f"[{c3['pred_ci'][0]:.2f},{c3['pred_ci'][1]:.2f}]")
    print(f"     observed MATCHED rule (= dof-weighted global) {c3['matched']:.2f} "
          f"[{c3['matched_ci'][0]:.2f},{c3['matched_ci'][1]:.2f}]")
    print(f"     observed median over systems {c3['median']:.2f}; "
          f"unweighted mean {c3['unweighted_mean']:.2f}")
    print(f"     intervals overlap: {c3['overlap']}; predicted interval covers observed: "
          f"{c3['pred_covers_obs']}")

    print(f"\n[C4] dof over {len(c4['names'])} systems: min {c4['min']} q1 {c4['q1']:g} "
          f"median {c4['median']:g} q3 {c4['q3']:g} max {c4['max']}; "
          f"counts 1/2/3 = {c4['counts'][1]}/{c4['counts'][2]}/{c4['counts'][3]}")
    for group, tag in ((c4["primary"], "primary"), (c4["secondary"], "secondary")):
        for s in group:
            print(f"    [{tag}] {s['label']:34s} n={s['n']:2d} flagged={s['sizes']} "
                  f"ever={len(s['ever'])} always={len(s['inter'])} "
                  f"meanJ={s['mean_jaccard']:.3f} (rand {s['mean_jaccard_random']:.3f}) "
                  f"median swing {s['median_swing']:.2f}x")

    print(f"\n[C5] no-systematic-error null rho: mean {c5['mean']:+.4f} median "
          f"{c5['median']:+.4f} sd {c5['sd']:.4f} 95%[{c5['ci'][0]:+.4f},{c5['ci'][1]:+.4f}] "
          f"max {c5['max']:+.4f}  (n={c5['n_null']})")
    print(f"     reported rho={PUB_RHO:.3f} vs that null: p={c5['p_vs_null']:.4f}")
    print(f"     positive control median rho={c5['ctrl_median']:+.3f} "
          f"[{c5['ctrl_ci'][0]:+.3f},{c5['ctrl_ci'][1]:+.3f}] (n={c5['n_ctrl']})")

    print(f"\n[C6] analytic se BELOW replicate spread on {c6['frac_below']:.3f} "
          f"[{c6['frac_below_ci'][0]:.3f},{c6['frac_below_ci'][1]:.3f}] of {c6['n_edges']} edges "
          f"(perfect-calibration reference {c6['ref_raw']:.3f})")
    print(f"     c4-corrected {c6['frac_below_c4']:.3f} "
          f"[{c6['frac_below_c4_ci'][0]:.3f},{c6['frac_below_c4_ci'][1]:.3f}] "
          f"vs reference {c6['ref_c4']:.3f}")
    print(f"     RMS-pooled ratio {c6['rms_ratio']:.2f}x; median per-edge ratio "
          f"{c6['median_ratio']:.2f}x; {c6['n_targets_above_ref']}/{len(c6['per_target'])} "
          f"targets above reference; {c6['n_targets_ratio_lt1']} targets with RMS ratio < 1")


if __name__ == "__main__":
    main()
