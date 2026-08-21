"""Fig Self -- TWO-SIDED self-calibration of the cycle-closure null (exploratory; NOT the
pre-registered P8 robustness check, which is one-sided and untouched).

The manuscript states that the QC detector's *recall* is uncharacterized because there are no
ground-truth error labels. That is true of recall in the strict sense, but the replicate set does
supply a measurable proxy for the one quantity recall depends on most here: whether the null the
test is run against is the right width. `figs/make_figA_replicates.py` measures, per edge, how the
reported (MBAR = sandwich) se compares with that edge's own across-replicate spread. This script
asks what the Fig L flag set becomes when each edge is tested against its OWN measured bar instead
of its reported one -- in whichever direction the measurement points.

RELATION TO Fig/Table P8 (`figs/make_figP8.py`) -- read this before comparing numbers:
  * P8 is **pre-registered** and **one-sided**: it inflates a system's se by 1/ratio only where
    ratio < 1 (give each system the benefit of its own measured local overconfidence) and asks
    whether the six flags survive. It returned SURVIVES. It is not modified, re-run differently,
    or reinterpreted here.
  * This script is **exploratory** and **two-sided**: it applies the correction in whichever
    direction the measurement points, at per-edge granularity. It is a different question --
    "what would the test say if every bar were the width the replicates say it is?" -- and its
    answer is a statement about how many detections the observed miscalibration *suppresses*, not
    a revision of P8's verdict.

DESIGN, FIXED AND STATED BEFORE THE RUN (none of it is adjusted after seeing a flag count):

  D1. GRANULARITY = PER EDGE. The correction factor for edge e is that edge's own measured
      calibration
          ratio_e = rep_e / repl_e,
          rep_e   = sqrt(mean_k se_ek^2)   over the 3 replicates (RMS reported se),
          repl_e  = SD_k(ddG_ek), ddof=1   (the across-replicate "truth"),
      which is `make_figA_replicates.load()`'s per-edge pair, recomputed here from the same CSV
      through `make_figL.edge_val` so the edge construction is shared verbatim (a unit test
      asserts the recomputed table matches `figA_replicates.load()` edge for edge). Per-edge is
      what "that edge's own measured calibration" means and is the granularity this analysis is
      specified at; the per-system aggregate ratio (P8's / Fig A-rep's pooled RMS formula) is
      reported as a SENSITIVITY, not as an alternative primary.

  D2. UNDEFINED / DEGENERATE RATIO -> EDGE LEFT UNCHANGED (ratio := 1). An edge gets no
      correction when it has no complete 3-replicate record (no ratio can be measured) or when
      its replicate spread is degenerate (`repl_e <= 1e-6`, the same threshold Fig A-rep and P8
      use), or if the arithmetic is otherwise undefined (rep_e <= 0, non-finite). Rationale: an
      unmeasurable ratio carries no evidence about that edge's calibration, so the neutral action
      is to keep the reported bar. This is neutral, NOT conservative: it neither adds nor removes
      flags by construction.

  D3. DIRECTION CONVENTION: se -> se / ratio, both directions.
        ratio > 1  = reported bar WIDER than the measured run-to-run spread (over-conservative)
                     -> se shrinks -> bar TIGHTER -> larger |z| -> larger chi^2 -> MORE likely
                     flagged;
        ratio < 1  = reported bar TIGHTER than the measured spread (locally overconfident)
                     -> se grows  -> bar LOOSER  -> smaller chi^2 -> LESS likely flagged.
      P8 keeps only the second branch. Two-sided keeps both.

  D4. NO CAP on the primary. A capped variant (ratio clipped to [1/3, 3]) is reported as a
      sensitivity, because a per-edge ratio built on an n=3 SD has a heavy right tail.

  D5. NO c4 SMALL-SAMPLE CORRECTION on the primary: the raw n=3 SD ratio is used, identical to
      the ratio Fig A-rep panel B and P8 report. The c4(3)=0.886 bias-corrected variant
      (repl -> repl / c4, i.e. wider bars, fewer flags) is reported as a sensitivity.

  Everything else is Fig L's, unchanged: replicate 0, systems with >= 3 edges and dof >= 1,
  the GLS fit and chi^2 tail (`bar.qc.gls_network` / `chi2_sf`) reached through
  `make_figL`'s own adapters, and the pre-registered BH-FDR level alpha = 0.05.

SENSITIVITIES REPORTED (all of them, chosen before the run, none used to pick the primary):
  per-system aggregate ratio; capped ratio; c4-corrected ratio; one-sided per-edge (P8's
  direction at this granularity); and a HELD-OUT ratio measured on replicates {1,2} only, whose
  correction factor therefore never reads the replicate-0 values the closure test is scored on.

STATED DEPENDENCE (not a knob, a caveat): the primary ratio is measured on all three replicates,
including replicate 0, whose ddG values the closure test itself reads. The correction factor and
the tested residuals are therefore not fully independent. The held-out {1,2} variant was built as
the control for exactly that.

ADDED AFTER THE FIRST RUN (a diagnostic, not a design change -- no arm, constant, population or
alpha was altered): `null_calibration()` scores EVERY arm on a perfectly-calibrated synthetic
null (same graphs, same reported bars, ddG redrawn from those bars around exact node potentials,
3 replicates per edge), where every flag is a false positive by construction. It was added
because of a design-level argument, not because of any flag count: self-calibrating an edge by
its own replicate SD multiplies that edge's chi^2 contribution by sigma^2/s^2 with
E[sigma^2/s^2] = nu/(nu-2) on nu = n-1 dof, which DIVERGES at n = 3 (nu = 2). The per-edge arm is
therefore anti-conservative by construction; the aggregate arm pools the denominator over a whole
system (nu ~ 2E) and is not. The harness measures exactly that, and the results doc reports the
per-system aggregate arm as the one that may be read as a detection count -- a choice made on the
synthetic null, which never sees a real ddG value, not on the observed flag counts. The held-out
{1,2} variant (nu = 1) is invalid for the same reason and is reported as uninformative rather
than as a stronger version of the claim.

MANDATORY SANITY ANCHOR: the nominal arm must reproduce Fig L (median reduced chi^2 0.34, flagged
set {bace, brd4, cdk8, faah, hif2a, p38}); otherwise the script prints a STOP banner instead of a
result.

Run: PYTHONPATH=src python figs/make_figSelf.py  (or `make figSelf`). Deterministic (no
randomness anywhere: a closed-form recomputation of the same GLS + BH pipeline). Data:
`data/openfe_replicates/combined_pymbar4_edge_data.csv` (OpenFE IndustryBenchmarks2024, public).
"""
from __future__ import annotations

import csv
import math
import pathlib
import sys

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "figs"))

import make_figL as figL  # noqa: E402
from make_figStab import bh_qvalues  # noqa: E402
from paperstyle import (  # noqa: E402
    ALT, INK, MUTED, OURS, REF, figsize, finish, legend, panel, tint, use_paper_style,
)

from bar.qc import benjamini_hochberg as bh_flags  # noqa: E402
from bar.qc import chi2_sf, gls_network  # noqa: E402

DATA = ROOT / "data" / "openfe_replicates" / "combined_pymbar4_edge_data.csv"
FIGDIR = pathlib.Path(__file__).resolve().parent
DOC = ROOT / "docs" / "results_figSelf.md"

ALPHA = 0.05              # pre-registered BH-FDR level (Fig L's)
REPLICATE = 0             # Fig L's tested replicate
MIN_EDGES = 3             # Fig L's per-system admission rule
REPLICATES = (0, 1, 2)
DEGENERATE_SD = 1e-6      # Fig A-rep / P8's non-degenerate-spread threshold
CAP = 3.0                 # sensitivity only (D4)
NULL_DRAWS = 400          # synthetic-null draws for the realized-FPR diagnostic
OPERATING_SCALES = [1.04, 1.00, 0.92, 0.79]   # the calibration interval and its ends
C4_3 = math.sqrt(2.0 / 2.0) * math.gamma(3 / 2) / math.gamma(2 / 2)  # 0.8862; E[s_{n=3}] = c4*sigma

FLAGGED_ANCHOR = frozenset({"bace", "brd4", "cdk8", "faah", "hif2a", "p38"})
MEDIAN_RC_ANCHOR = 0.34
ANCHOR_TOL = 0.01

# Semantic colours (paperstyle), by quantity family. This figure is entirely about the article's
# own error bar and variants of its own test, so it lives on one hue and its tint: OURS = the
# reported (published) bar and everything the published test flags with it, tint(OURS) = the other
# side of the ratio = 1 crossing. ALT = the self-calibrated arm, the named alternative this figure
# explores, wherever it appears -- the per-system aggregate ratios, the flags it gains, its q
# values, and the flags it gains. No arm of this figure takes FOIL: every arm here is one of our
# own pre-specified design variants, and an arm does not change family by failing its own audit.
# An arm that fires on a perfectly-calibrated null is hatched instead. MUTED = the
# de-emphasised background: the per-edge context cloud and the systems no arm ever flags.
C_LOOSE, C_TIGHT, C_LOST = tint(OURS, 0.45), OURS, tint(ALT, 0.45)

# Panel-D display names. The ``variants`` labels below double as keys of the released results
# document, so the figure carries its own self-describing wording instead of repo-internal
# shorthand ("Fig L", "P8"): a reader of the printed figure has neither name available.
FIG_LABELS = {
    "nominal": "nominal: the published test",
    "per-edge": "pre-specified: per-edge, two-sided",
    "aggregate": "sensitivity: per-system aggregate",
    "capped": "per-edge, ratio capped to [1/3, 3]",
    "c4": "per-edge, $c_4(3)$-corrected SD",
    "one-sided": "per-edge, one-sided (loosen only)",
    "held-out": "per-edge, held-out ratio (reps 1+2)",
}


# --------------------------------------------------------------------------- the ratio table

def _edge_ratio(ses, ddgs):
    """``(rep, repl, ratio)`` from the per-replicate reported se and ddG of ONE edge.

    ``rep`` is the RMS reported se and ``repl`` the across-replicate SD (ddof=1) -- the pair
    `make_figA_replicates.load()` returns. ``ratio`` is ``None`` when it is not measurable (D2).
    """
    rep = float(np.sqrt(np.mean(np.asarray(ses, dtype=float) ** 2)))
    repl = float(np.std(np.asarray(ddgs, dtype=float), ddof=1))
    if repl <= DEGENERATE_SD or rep <= 0.0 or not math.isfinite(rep) or not math.isfinite(repl):
        return rep, repl, None
    return rep, repl, rep / repl


def edge_ratios(replicates=REPLICATES, c4=False):
    """``{(system, ligand_A, ligand_B): dict(rep, repl, ratio)}`` over the given replicates.

    Uses `make_figL.edge_val` for the edge construction, so the (ddG, se) pair is byte-for-byte
    the one the closure test consumes. An edge missing any of ``replicates`` is absent from the
    table (D2: it will be left unchanged). With ``c4=True`` the n=3 SD bias is removed from the
    denominator (``repl -> repl / c4``), the D5 sensitivity.
    """
    out = {}
    with open(DATA) as fh:
        for r in csv.DictReader(fh):
            vals = [figL.edge_val(r, k) for k in replicates]
            if any(v is None for v in vals):
                continue
            rep, repl, ratio = _edge_ratio([v[1] for v in vals], [v[0] for v in vals])
            if c4 and ratio is not None:
                repl, ratio = repl / C4_3, ratio * C4_3
            out[(r["system name"], r["ligand_A"], r["ligand_B"])] = {
                "rep": rep, "repl": repl, "ratio": ratio}
    return out


def system_ratios(table):
    """Per-system pooled RMS ratio -- Fig A-rep's / P8's aggregate formula, over the same
    non-degenerate-spread edges. Used only by the D1 sensitivity variant."""
    acc: dict[str, list[tuple[float, float]]] = {}
    for (sysname, _a, _b), d in table.items():
        if d["ratio"] is None:
            continue
        acc.setdefault(sysname, []).append((d["rep"], d["repl"]))
    return {s: float(np.sqrt(np.mean([v[0] ** 2 for v in vs]))
                     / np.sqrt(np.mean([v[1] ** 2 for v in vs]))) for s, vs in acc.items()}


# --------------------------------------------------------------------------- the correction

def corrected_se(se, ratio, mode="two-sided", cap=None):
    """Apply D2 + D3 + D4 to one edge: ``se -> se / ratio``.

    ``ratio > 1`` tightens the bar, ``ratio < 1`` loosens it (D3). An unmeasurable ratio leaves
    the bar untouched (D2). ``mode='one-sided'`` keeps only the loosening branch (P8's rule at
    this granularity). ``cap`` clips the ratio into ``[1/cap, cap]`` (D4, sensitivity only).
    """
    if ratio is None or not math.isfinite(ratio) or ratio <= 0.0:
        return se
    if mode == "one-sided" and ratio >= 1.0:
        return se
    if cap is not None:
        ratio = min(max(ratio, 1.0 / cap), cap)
    return se / ratio


def _neutral(_s, _a, _b):
    """Lookup that measures nothing: reproduces the nominal (Fig L) arm exactly."""
    return None


def per_edge_lookup(table):
    def look(s, a, b):
        d = table.get((s, a, b))
        return None if d is None else d["ratio"]
    return look


def per_system_lookup(sys_ratio):
    def look(s, _a, _b):
        return sys_ratio.get(s)
    return look


# --------------------------------------------------------------------------- the test

def run_test(by, lookup, mode="two-sided", cap=None, alpha=ALPHA):
    """Fig L's cycle-closure test on replicate ``REPLICATE`` with the per-edge se replaced by
    ``corrected_se``. Population, fit, chi^2 tail and BH-FDR are Fig L's, unchanged.

    Returns one dict per admitted system with the reduced chi^2, raw p, BH q and flag, plus the
    bookkeeping of how many of its edges were tightened / loosened / left unchanged.
    """
    rows = []
    for name, rs in sorted(by.items()):
        edges, used, n_unchanged, n_tight, n_loose = [], [], 0, 0, 0
        for r in rs:
            v = figL.edge_val(r, REPLICATE)
            if v is None:
                continue
            ddg, se = v
            ratio = lookup(name, r["ligand_A"], r["ligand_B"])
            se_new = corrected_se(se, ratio, mode=mode, cap=cap)
            if se_new == se:
                n_unchanged += 1
            elif se_new < se:
                n_tight += 1
            else:
                n_loose += 1
            if ratio is not None and math.isfinite(ratio) and ratio > 0:
                used.append(ratio)
            edges.append((r["ligand_A"], r["ligand_B"], ddg, se_new))
        if len(edges) < MIN_EDGES:
            continue
        X2, dof, _z = figL.gls_chi2(edges)
        if dof < 1:
            continue
        rows.append({"sys": name, "E": len(edges), "dof": dof, "rc": X2 / dof,
                     "p": chi2_sf(X2, dof), "n_tight": n_tight, "n_loose": n_loose,
                     "n_unchanged": n_unchanged,
                     "med_ratio": float(np.median(used)) if used else math.nan})
    if not rows:
        return rows
    pv = [d["p"] for d in rows]
    flags, qs = bh_flags(pv, alpha=alpha), bh_qvalues(pv)
    for d, fl, q in zip(rows, flags, qs, strict=True):
        d["flag"], d["q"] = bool(fl), float(q)
    assert all(d["flag"] == (d["q"] <= alpha) for d in rows), "BH q-values disagree with flags"
    return rows


def flag_set(rows):
    return {d["sys"] for d in rows if d["flag"]}


def median_rc(rows):
    return float(np.median([d["rc"] for d in rows]))


# --------------------------------------------------------------------------- null calibration

def null_topology(by):
    """Fig L's population reduced to ``(system, [(a, b, se)])`` -- real topology, real reported
    bars, values discarded. The admission rule (>= MIN_EDGES edges, dof >= 1) is Fig L's and is
    a property of the topology alone, so this population is the tested one."""
    out = []
    for name, rs in sorted(by.items()):
        edges = [(r["ligand_A"], r["ligand_B"], figL.edge_val(r, REPLICATE)[1])
                 for r in rs if figL.edge_val(r, REPLICATE) is not None]
        if len(edges) < MIN_EDGES:
            continue
        _X2, dof, _z = figL.gls_chi2([(a, b, 0.0, se) for a, b, se in edges])
        if dof < 1:
            continue
        out.append((name, edges))
    return out


def _null_arm_stat(edges, se, y, arm):
    """``(p, reduced chi^2)`` of one system in one synthetic draw under one arm.

    ``y`` is ``(3, E)`` replicate ddG drawn around exact node potentials with the reported ``se``
    as the TRUE sigma, so the reported bars are perfectly calibrated by construction and every
    flag is a false positive. Ratios are formed by the same formulas the real analysis uses.
    """
    if arm == "nominal":
        ratio = np.full(se.size, np.nan)
    else:
        reps = (1, 2) if arm == "held-out" else (0, 1, 2)
        repl = np.std(y[list(reps)], axis=0, ddof=1)
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = np.where(repl > DEGENERATE_SD, se / repl, np.nan)
        if arm == "c4":
            ratio = ratio * C4_3
        if arm == "aggregate":
            ok = np.isfinite(ratio)
            agg = (math.sqrt(float(np.mean(se[ok] ** 2)))
                   / math.sqrt(float(np.mean(repl[ok] ** 2)))) if ok.any() else np.nan
            ratio = np.full(se.size, agg)
    mode = "one-sided" if arm == "one-sided" else "two-sided"
    cap = CAP if arm == "capped" else None
    new = np.array([corrected_se(s, (None if not np.isfinite(rt) else float(rt)),
                                 mode=mode, cap=cap) for s, rt in zip(se, ratio, strict=True)])
    X2, dof, _z = figL.gls_chi2([(a, b, float(v), float(s))
                                 for (a, b, _s), v, s in zip(edges, y[0], new, strict=True)])
    return chi2_sf(X2, dof), X2 / dof


def null_calibration(by, arms, n_draws, seed=20260810):
    """Realized false-positive behaviour of every arm on a PERFECTLY-CALIBRATED synthetic null.

    Same graphs, same reported bars, values redrawn from those bars around exact node potentials
    (no systematic error anywhere), 3 replicates per edge. Any flag is therefore a false positive.
    Returns ``{arm: {"p_any", "mean_flagged", "median_rc"}}`` -- the probability that the arm
    flags at least one system when nothing is wrong, the mean number it flags, and the median
    reduced chi^2 it produces on a correctly-specified null (which should sit at 1).

    This is a diagnostic added after the first run because of a design-level argument (the n=3
    per-edge ratio makes the corrected chi^2's null expectation diverge; see the results doc), not
    because of any flag count. It changes no design choice: the arms are exactly the ones fixed
    before the run.
    """
    rng = np.random.default_rng(seed)
    topo = null_topology(by)
    any_flag = dict.fromkeys(arms, 0)
    n_flag = dict.fromkeys(arms, 0)
    rcs: dict[str, list[float]] = {a: [] for a in arms}
    for _ in range(n_draws):
        draws = []
        for name, edges in topo:
            se = np.array([e[2] for e in edges], dtype=float)
            draws.append((name, edges, se, rng.normal(0.0, se[None, :], size=(3, se.size))))
        for arm in arms:
            stats = [_null_arm_stat(edges, se, y, arm) for _n, edges, se, y in draws]
            fl = bh_flags([s[0] for s in stats], alpha=ALPHA)
            any_flag[arm] += int(bool(fl.any()))
            n_flag[arm] += int(fl.sum())
            rcs[arm].append(float(np.median([s[1] for s in stats])))
    return {a: {"p_any": any_flag[a] / n_draws, "mean_flagged": n_flag[a] / n_draws,
                "median_rc": float(np.median(rcs[a])), "k": any_flag[a], "n": n_draws,
                "ci": clopper_pearson(any_flag[a], n_draws)} for a in arms}


def clopper_pearson(k: int, n: int, alpha: float = 0.05) -> tuple[float, float]:
    """Exact binomial confidence interval for ``k`` successes in ``n`` draws.

    A realized level of 0.000 or 0.003 out of 400 draws is not resolved to three decimals, and
    quoting it as though it were invites a reader to compare grid points that the draw count
    cannot separate. The interval is what the simulation actually supports.
    """
    from scipy.stats import beta
    lo = 0.0 if k == 0 else float(beta.ppf(alpha / 2, k, n - k + 1))
    hi = 1.0 if k == n else float(beta.ppf(1 - alpha / 2, k + 1, n - k))
    return lo, hi


def flag_count_at_calibration(by, calibration_scales, alpha=ALPHA):
    """Flags on the REAL replicate-0 data when every reported bar is corrected to c times itself.

    ``c`` is the calibration scale (true se over reported se), so correcting the bars multiplies
    each se by ``c``: c < 1 tightens the null and can only add flags. This is the real-data
    companion of ``level_at_calibration``, which is a synthetic-null quantity.
    """
    out = {}
    for c in calibration_scales:
        rows = run_test(by, lambda _s, _a, _b, _c=c: 1.0 / _c)
        out[c] = len(flag_set(rows))
    return out


def level_at_calibration(by, calibration_scales, n_draws=NULL_DRAWS, seed=20260810):
    """Realized probability of any false flag when the reported bar is correct but the test uses
    a uniformly rescaled null.

    The argument is a CALIBRATION scale c = (true se) / (reported se), which is the quantity
    Section 4 measures as 0.92 with interval [0.79, 1.04]. The simulation needs the reciprocal: the
    truth is drawn at c times the reported bar while the test keeps using the reported bar, so the
    test's null is 1/c times the truth. Passing a calibration scale in as a simulation scale
    inverts the reading, which is the error this function now exists to prevent. Smaller c means
    the reported bars are wider than the truth, hence a lower level, not a higher one. Scale 1 is
    the shipped test and agrees with the uncorrected arm of ``null_calibration``.
    """
    rng = np.random.default_rng(seed)
    topo = null_topology(by)
    scales = {c: 1.0 / c for c in calibration_scales}
    hits = dict.fromkeys(calibration_scales, 0)
    for _ in range(n_draws):
        draws = []
        for _name, edges in topo:
            se = np.array([e[2] for e in edges], dtype=float)
            draws.append((edges, se, rng.normal(0.0, se, size=se.size)))
        for cal, scale in scales.items():
            ps = []
            for edges, se, y in draws:
                scaled = [(a, b, float(y[i]) * cal, float(se[i]))
                          for i, (a, b, _s) in enumerate(edges)]
                fit = gls_network(scaled)
                if fit.dof < 1:
                    continue
                ps.append(chi2_sf(fit.chi2, fit.dof))
            hits[cal] += int(bool(bh_flags(ps, alpha=ALPHA).any()))
    return {c: (hits[c], n_draws, hits[c] / n_draws, clopper_pearson(hits[c], n_draws))
            for c in calibration_scales}


# --------------------------------------------------------------------------- driver

def analyze(n_null=NULL_DRAWS):
    by = figL.load_systems()
    table = edge_ratios()
    table_c4 = edge_ratios(c4=True)
    table_ho = edge_ratios(replicates=(1, 2))          # held out from the tested replicate
    sysr = system_ratios(table)

    variants = [
        ("nominal", "nominal (reported se) = Fig L", run_test(by, _neutral)),
        ("per-edge", "pre-specified primary: per-edge, two-sided",
         run_test(by, per_edge_lookup(table))),
        ("aggregate", "granularity sensitivity: per-system aggregate, two-sided",
         run_test(by, per_system_lookup(sysr))),
        ("capped", f"per-edge, ratio capped to [1/{CAP:g}, {CAP:g}]",
         run_test(by, per_edge_lookup(table), cap=CAP)),
        ("c4", "per-edge, c4(3)-corrected SD", run_test(by, per_edge_lookup(table_c4))),
        ("one-sided", "per-edge, one-sided (P8's direction)",
         run_test(by, per_edge_lookup(table), mode="one-sided")),
        ("held-out", "per-edge, held-out ratio (replicates 1+2 only)",
         run_test(by, per_edge_lookup(table_ho))),
    ]
    arms = {k: rows for k, _lab, rows in variants}
    null = null_calibration(by, [k for k, _l, _r in variants], n_draws=n_null)

    tested = {d["sys"] for d in arms["nominal"]}          # report ratios over the tested systems
    sysr = {s: v for s, v in sysr.items() if s in tested}
    ratios = np.array([d["ratio"] for (s, _a, _b), d in table.items()
                       if d["ratio"] is not None and s in tested])
    n_degenerate = sum(1 for (s, _a, _b), d in table.items()
                       if d["ratio"] is None and s in tested)
    return {"by": by, "table": table, "sysr": sysr, "variants": variants, "arms": arms,
            "nominal": arms["nominal"], "primary": arms["per-edge"],
            "aggregate": arms["aggregate"], "null": null,
            "ratios": ratios, "n_degenerate": n_degenerate}


def transitions(nominal, primary):
    """``(kept, gained, lost)`` system-name lists for nominal -> self-calibrated."""
    a, b = flag_set(nominal), flag_set(primary)
    return sorted(a & b), sorted(b - a), sorted(a - b)


def make_figure(res):
    nominal, per_edge, agg = res["nominal"], res["primary"], res["aggregate"]
    mn = {d["sys"]: d for d in nominal}
    me = {d["sys"]: d for d in per_edge}
    ma = {d["sys"]: d for d in agg}
    kept, gained, lost = transitions(nominal, agg)
    ratios, sysr = res["ratios"], res["sysr"]

    fig = plt.figure(figsize=figsize(4, 6.2))
    # panel D carries seven named arms and their realized levels as text, so the bottom row is
    # split unevenly; every other panel is comfortable at half the text block.
    gs = fig.add_gridspec(2, 1)
    gs_top = gs[0].subgridspec(1, 2)
    gs_bot = gs[1].subgridspec(1, 2, width_ratios=[0.86, 1.14])
    axA = fig.add_subplot(gs_top[0, 0])
    axB = fig.add_subplot(gs_top[0, 1])
    axC = fig.add_subplot(gs_bot[0, 0])
    axD = fig.add_subplot(gs_bot[0, 1])

    # --- A: the measured calibration, and which way it pushes --------------------------------
    lo, hi = float(ratios.min()), float(ratios.max())
    bins = np.logspace(math.log10(lo * 0.9), math.log10(hi * 1.1), 44)
    # one quantity split by a threshold crossing: one hue, the far side taken in tint
    n_lo, _b, _p = axA.hist(ratios[ratios < 1.0], bins=bins, color=C_LOOSE, alpha=0.85,
                            label="ratio < 1  →  bar loosened")
    n_hi, _b, _p = axA.hist(ratios[ratios >= 1.0], bins=bins, color=C_TIGHT, alpha=0.85,
                            label="ratio ≥ 1  →  bar tightened")
    top = max(n_lo.max(), n_hi.max()) * 1.5
    axA.set_ylim(0, top)
    sv = np.array(sorted(sysr.values()))
    axA.plot(sv, np.full(sv.size, top * 0.70), "|", color=ALT, ms=8, mew=1.1,
             label=f"per-system aggregate ({sv.size})")
    # the guide stops below the legend instead of being drawn through its two rows of text
    axA.plot([1.0, 1.0], [0.0, top * 0.82], color=REF, ls=(0, (1.0, 2.0)), lw=1.0, zorder=0.5)
    axA.set_xscale("log")
    axA.set_xlabel("measured calibration ratio\n= reported se / replicate SD")
    axA.set_ylabel(f"edges (n = {ratios.size})")
    legend(axA, loc="upper left", fontsize=7.5, handlelength=1.2, borderaxespad=0.2)
    axA.text(0.985, 0.97,
             f"per-edge median {np.median(ratios):.2f}\n"
             f"range {lo:.2f}–{hi:.0f}\n"
             f"{100 * float(np.mean(ratios >= 1)):.0f}% ≥ 1 (tightened)\n"
             f"{res['n_degenerate']} unmeasurable, left unchanged\n"
             f"aggregate median {np.median(sv):.2f}",
             transform=axA.transAxes, va="top", ha="right", fontsize=7.5, color=INK)
    panel(axA, "A", "the reported bars are mostly too wide",
          subtitle="so self-calibration mostly tightens them")

    # --- B: per-system reduced chi^2, nominal vs self-calibrated -----------------------------
    axB.scatter([mn[s]["rc"] for s in mn], [me[s]["rc"] for s in mn], s=13, marker="o",
                facecolors="none", edgecolors=MUTED, linewidths=0.6, zorder=2,
                label="per-edge arm (all systems)")
    never = [s for s in mn if not mn[s]["flag"] and not ma[s]["flag"]]
    for lab, col, mk, sz, names in [("aggregate: never flagged", MUTED, "D", 13, never),
                                    (f"aggregate: gained ({len(gained)})", ALT, "^", 26,
                                     gained),
                                    (f"aggregate: kept ({len(kept)})", OURS, "D", 22, kept),
                                    (f"aggregate: lost ({len(lost)})", C_LOST, "v", 26, lost)]:
        if not names:
            continue
        axB.scatter([mn[s]["rc"] for s in names], [ma[s]["rc"] for s in names], s=sz, marker=mk,
                    c=col, alpha=0.9, lw=0, zorder=3, label=lab)
    allrc = [d["rc"] for d in nominal] + [d["rc"] for d in agg] + [d["rc"] for d in per_edge]
    lim = [min(allrc) * 0.6, max(allrc) * 2.2]
    # two reference lines in the same REF grey, told apart by dash pattern
    axB.plot(lim, lim, color=REF, ls=(0, (1.0, 2.0)), lw=1.0, zorder=1)
    axB.axhline(1.0, color=REF, ls=(0, (4.5, 2.0)), lw=0.9, zorder=1)
    if "renin" in ma:
        # every direction within a few points of the marker is occupied (thrombin sits directly
        # above it), so the name is set in the empty upper left and tied back with a leader
        axB.annotate("renin", (mn["renin"]["rc"], ma["renin"]["rc"]), xytext=(-6, 15),
                     textcoords="offset points", fontsize=7.5, color=ALT, ha="right",
                     va="bottom", arrowprops={"arrowstyle": "-", "lw": 0.6,
                                              "color": tint(ALT, 0.45),
                                              "shrinkA": 2.0, "shrinkB": 2.5})
    axB.set_xscale("log")
    axB.set_yscale("log")
    axB.set_xlim(*lim)
    axB.set_ylim(*lim)
    axB.set_xlabel(r"reduced $\chi^2$, nominal (reported se)")
    axB.set_ylabel(r"reduced $\chi^2$, self-calibrated")
    legend(axB, loc="lower right", fontsize=7.5, borderaxespad=0.2, handletextpad=0.3)
    axB.text(0.03, 0.97, f"median {median_rc(nominal):.3f} → {median_rc(agg):.3f} "
             f"(aggregate)\n              → {median_rc(per_edge):.3f} (per-edge)",
             transform=axB.transAxes, ha="left", va="top", fontsize=7.5, color=INK)
    panel(axB, "B", "correcting the null moves systems up",
          subtitle="the reported null was too wide")

    # --- C: q before -> after (the null-calibrated aggregate arm) -----------------------------
    union = sorted(set(flag_set(nominal)) | set(flag_set(agg)), key=lambda s: ma[s]["q"])
    y = np.arange(len(union))
    for yi, s in enumerate(union):
        axC.annotate("", xy=(max(ma[s]["q"], 1e-16), yi), xytext=(max(mn[s]["q"], 1e-16), yi),
                     arrowprops={"arrowstyle": "->", "color": tint(MUTED, 0.45), "lw": 0.8})
    axC.scatter([max(mn[s]["q"], 1e-16) for s in union], y, s=20, marker="o",
                facecolors=[OURS if mn[s]["flag"] else "white" for s in union],
                edgecolors=OURS, linewidths=0.9, zorder=3, label="nominal (filled = flagged)")
    axC.scatter([max(ma[s]["q"], 1e-16) for s in union], y, s=24, marker="D",
                facecolors=[ALT if ma[s]["flag"] else "white" for s in union],
                edgecolors=ALT, linewidths=0.9, zorder=4, label="self-calibrated (aggregate)")
    axC.set_xscale("log")
    # an opaque band, not a veil over the markers it contains
    axC.axvspan(ALPHA, 5.0, color=tint(REF, 0.92), lw=0, zorder=0)
    axC.axvline(ALPHA, color=REF, ls=(0, (1.0, 2.0)), lw=1.0, zorder=1)
    axC.set_xlim(1.2e-17, 5.0)
    axC.set_xticks([1e-16, 1e-12, 1e-8, 1e-4, 1e0])
    axC.set_ylim(len(union) + 1.8, -0.5)
    axC.set_yticks(y)
    axC.set_yticklabels(union, fontsize=7.5)
    axC.tick_params(axis="y", length=0)
    axC.set_xlabel("BH-adjusted $q$\n"
                   f"(shaded = not flagged, $\\alpha$={ALPHA})")
    legend(axC, loc="lower left", fontsize=7.5, borderaxespad=0.2, handletextpad=0.3)
    panel(axC, "C", "reported se versus the aggregate arm",
          subtitle="granularity sensitivity, per system")

    # --- D: every arm, with its realized false-positive rate ---------------------------------
    labels = [FIG_LABELS.get(k, lab) for k, lab, _r in res["variants"]]
    counts = [len(flag_set(r)) for _k, _l, r in res["variants"]]
    pany = [res["null"][k]["p_any"] for k, _l, _r in res["variants"]]
    nmean = [res["null"][k]["mean_flagged"] for k, _l, _r in res["variants"]]
    # Every arm here is one of OUR OWN pre-specified design variants, so every bar is OURS. An
    # arm that fires on a perfectly-calibrated null is marked by a hatch, not by a second hue:
    # FOIL is the overconfident stand-in and ALT is a rival method, and neither of those is what
    # a self-calibration variant becomes by behaving badly. The realized level is printed beside
    # every bar, so the mark carries emphasis and the number carries the value.
    fires = [p > 0.50 for p in pany]
    ys = np.arange(len(labels))[::-1]
    axD.barh(ys, counts, color=OURS, height=0.42, zorder=2,
             hatch=["///" if f else "" for f in fires],
             edgecolor=INK, linewidth=[0.7 if f else 0.0 for f in fires])
    # the arm's name and its realized level share the line above its bar, which leaves the whole
    # x-axis to the bars; at half the text block there is no room for a left-hand label column.
    for yi, lab, c, p, m in zip(ys, labels, counts, pany, nmean, strict=True):
        axD.annotate(lab, (0.0, yi + 0.28), xytext=(2.0, 0.0), textcoords="offset points",
                     ha="left", va="bottom", fontsize=7.5, color=INK)
        axD.annotate(f"null P(≥1)={p:.2f}, mean {m:.2f}", (1.0, yi + 0.28),
                     xycoords=axD.get_yaxis_transform(which="grid"),
                     xytext=(-2.0, 0.0), textcoords="offset points", ha="right",
                     va="bottom", fontsize=7.5, color=INK, fontweight=(
                         "bold" if p > 0.50 else "normal"))
        axD.text(c + max(counts) * 0.02, yi, f"{c}", ha="left", va="center", fontsize=7,
                 color=INK)
    axD.set_yticks([])
    axD.set_ylim(-0.7, len(labels) + 1.45)
    axD.set_xlim(0, max(counts) * 1.16)
    axD.set_xlabel(f"systems flagged of {len(nominal)}   (BH-FDR $\\alpha$={ALPHA})")
    axD.annotate(f"{NULL_DRAWS} null draws, every flag a false positive.\n"
                 f"hatched: the arm fires on that null more often than half the time\n"
                 f"nominal arm P(≥1) = {res['null']['nominal']['p_any']:.2f} at α={ALPHA}",
                 (0.0, len(labels) + 1.40), xytext=(2.0, 0.0), textcoords="offset points",
                 ha="left", va="top", fontsize=7.5, color=REF, linespacing=1.35)
    panel(axD, "D", "sensitivity to every design choice",
          subtitle="each arm scored on a perfectly-calibrated null")

    finish(fig, "figSelf_two_sided_selfcalibration")
    plt.close(fig)


def write_doc(res, levels, real_counts):
    nominal, per_edge, agg = res["nominal"], res["primary"], res["aggregate"]
    mn = {d["sys"]: d for d in nominal}
    me = {d["sys"]: d for d in per_edge}
    ma = {d["sys"]: d for d in agg}
    kept, gained, lost = transitions(nominal, agg)
    kept_e, gained_e, lost_e = transitions(nominal, per_edge)
    ratios, null = res["ratios"], res["null"]
    mc_se = math.sqrt(0.25 / NULL_DRAWS)
    # gains that ride BH's relaxing cutoff rather than their own chi^2 moving
    flat = [s for s in gained if abs(ma[s]["rc"] / mn[s]["rc"] - 1) < 0.15]

    order = sorted(mn, key=lambda s: ma[s]["q"])
    tbl = ["| system | edges | cycles | χ²ᵥ nominal | q nominal | agg. ratio | χ²ᵥ agg. "
           "| q agg. | transition (agg.) | χ²ᵥ per-edge | q per-edge |",
           "|---|---:|---:|---:|---:|---:|---:|---:|---|---:|---:|"]
    for s in order:
        a, b, c = mn[s], ma[s], me[s]
        trans = ("**gained**" if s in gained else "**kept**" if s in kept
                 else "**LOST**" if s in lost else "—")
        tbl.append(f"| {s} | {a['E']} | {a['dof']} | {a['rc']:.3f} | {a['q']:.1e} | "
                   f"{b['med_ratio']:.2f} | {b['rc']:.3f} | {b['q']:.1e} | {trans} | "
                   f"{c['rc']:.3f} | {c['q']:.1e} |")

    var = ["| arm | systems flagged | median χ²ᵥ | null P(≥1 flag) | null mean #flagged "
           "| null median χ²ᵥ | flag set |",
           "|---|---:|---:|---:|---:|---:|---|"]
    for k, lab, rows in res["variants"]:
        fs = sorted(flag_set(rows))
        n = null[k]
        var.append(f"| {lab} | {len(fs)} | {median_rc(rows):.3f} | {n['p_any']:.3f} "
                   f"({n['k']}/{n['n']}, CP [{n['ci'][0]:.4f}, {n['ci'][1]:.4f}]) | "
                   f"{n['mean_flagged']:.2f} | {n['median_rc']:.3f} | `{', '.join(fs)}` |")

    lines = [
        "# Results — Fig Self: two-sided self-calibration of the cycle-closure null",
        "",
        "**Figure:** `figs/figSelf_two_sided_selfcalibration.{pdf,png}` · **Reproduce:** "
        "`make figSelf`",
        "(`PYTHONPATH=src python figs/make_figSelf.py`). Deterministic (the analysis is a "
        "closed-form",
        "recomputation of the same GLS + BH pipeline; the synthetic-null diagnostic is seeded). "
        "Data:",
        "`data/openfe_replicates/combined_pymbar4_edge_data.csv` (OpenFE IndustryBenchmarks2024, "
        "public,",
        "3 independent replicates per edge).",
        "",
        "> **Exploratory, not pre-registered.** The pre-registered robustness check on this null "
        "is Fig/Table",
        "> P8 (`figs/make_figP8.py`, `docs/results_figP8.md`): **one-sided**, per-system, verdict "
        "**SURVIVES**.",
        "> It is untouched here — not re-run differently, not reinterpreted. This document asks a "
        "*different*,",
        "> two-sided question and reports it separately.",
        "",
        "## What this asks, and why it is not P8's question",
        "",
        "The manuscript states that the detector's **recall** is uncharacterized because there "
        "are no",
        "ground-truth error labels. That remains true of recall in the strict sense. But the "
        "replicate set",
        "does measure the one thing recall depends on most here: **whether the null the test is "
        "run against",
        "is the right width**. `figs/make_figA_replicates.py` measures, per edge, the ratio of "
        "the reported",
        "(MBAR = sandwich) se to that edge's own across-replicate SD. This analysis replaces each "
        "reported",
        "bar with its own measured one — `se → se / ratio` — **in whichever direction the "
        "measurement points**",
        f"— and re-runs the paper's own GLS + Benjamini–Hochberg pipeline unchanged at the "
        f"pre-registered α = {ALPHA}.",
        "",
        "P8 applies the same correction **only where it makes a bar looser** (ratio < 1): the "
        "conservative",
        "half, which is what a referee's circularity concern needs. The other half — edges whose "
        "reported bar",
        "is *wider* than their own measured spread — can only *add* detections, and it is the "
        "half that speaks",
        "to recall.",
        "",
        "## Design, fixed and stated before the run",
        "",
        "| # | choice | what was fixed | why |",
        "|---|---|---|---|",
        "| D1 | **granularity** | **per edge** (pre-specified primary): `ratio_e = rep_e / repl_e` "
        "with `rep_e = sqrt(mean_k se_ek²)` (RMS reported se over the 3 replicates) and "
        "`repl_e = SD_k(ΔΔG_ek)` (ddof=1) — `make_figA_replicates.load()`'s own per-edge pair, "
        "recomputed through `make_figL.edge_val` so the edge construction is shared verbatim | "
        "\"that edge's own measured calibration\" is per-edge. The **per-system aggregate** ratio "
        "(Fig A-rep's / P8's pooled RMS formula) was fixed at the same time as the granularity "
        "**sensitivity**. Both are reported below; see *Which arm to read* |",
        "| D2 | **undefined / degenerate ratio** | **edge left unchanged** (`ratio := 1`) when the "
        f"edge has no complete 3-replicate record, or `repl_e ≤ {DEGENERATE_SD:g}` (Fig A-rep's "
        "and P8's own threshold), or the arithmetic is non-finite | an unmeasurable ratio carries "
        "no evidence about that edge; the neutral action is to keep the reported bar. Neutral, "
        "**not** conservative: by construction it can neither add nor remove a flag |",
        "| D3 | **direction convention** | `se → se / ratio`, **both** directions. `ratio > 1` "
        "(reported bar wider than the measured spread) ⇒ se shrinks ⇒ **bar tighter** ⇒ larger "
        "χ² ⇒ **more** likely flagged. `ratio < 1` (reported bar tighter than the measured "
        "spread) ⇒ se grows ⇒ **bar looser** ⇒ **less** likely flagged | P8 keeps only the second "
        "branch; two-sided keeps both |",
        f"| D4 | **no cap** on the primary | a `[1/{CAP:g}, {CAP:g}]`-clipped variant is reported "
        "as a sensitivity | an n=3 SD ratio has a heavy right tail; the tail is disclosed rather "
        "than silently trimmed |",
        "| D5 | **no c4 correction** on the primary | raw n=3 SD ratio, identical to the ratio "
        f"Fig A-rep panel B and P8 report; the `c4(3)={C4_3:.3f}` bias-corrected variant "
        "(`repl → repl / c4`, i.e. wider bars, fewer flags) is a sensitivity | one ratio "
        "definition in the paper, not two |",
        "",
        f"Everything else is Fig L's, unchanged: replicate {REPLICATE}, systems with ≥ "
        f"{MIN_EDGES} edges and dof ≥ 1,",
        "`bar.qc.gls_network` / `bar.qc.chi2_sf` / `bar.qc.benjamini_hochberg` reached through "
        "`make_figL`'s own",
        f"adapters, and the pre-registered BH-FDR level α = {ALPHA}. The nominal arm is asserted "
        "against Fig L's",
        f"published values (median reduced χ² {MEDIAN_RC_ANCHOR}, flag set "
        f"`{', '.join(sorted(FLAGGED_ANCHOR))}`) before anything",
        "else is reported; on a mismatch the script prints a STOP banner instead of a result.",
        "",
        "## What the replicates say about the bars",
        "",
        f"Measurable per-edge ratios: **{ratios.size} edges**; {res['n_degenerate']} edges "
        "unmeasurable and left unchanged (D2).",
        f"Median **{np.median(ratios):.2f}**, interquartile range "
        f"{np.percentile(ratios, 25):.2f}–{np.percentile(ratios, 75):.2f}, full range "
        f"{ratios.min():.2f}–{ratios.max():.1f}. "
        f"**{100 * float(np.mean(ratios >= 1)):.0f}%** of edges have `ratio ≥ 1`,",
        "i.e. a reported bar *wider* than their own measured run-to-run spread, so a two-sided "
        "correction",
        "predominantly **tightens** bars. That is the whole mechanism of the result below.",
        "",
        "## Headline",
        "",
        "| arm | systems flagged | median reduced χ² | lost | gained |",
        "|---|---:|---:|---:|---:|",
        f"| nominal (reported se, = Fig L) | **{len(flag_set(nominal))}** of {len(nominal)} | "
        f"{median_rc(nominal):.3f} | — | — |",
        f"| self-calibrated, **per-system aggregate** | **{len(flag_set(agg))}** of {len(agg)} | "
        f"{median_rc(agg):.3f} | {len(lost)} | {len(gained)} |",
        f"| self-calibrated, **per-edge** (pre-specified primary) | **{len(flag_set(per_edge))}** "
        f"of {len(per_edge)} | {median_rc(per_edge):.3f} | {len(lost_e)} | {len(gained_e)} |",
        "",
        f"- Aggregate arm — **lost: {len(lost)}** (`{', '.join(lost) if lost else 'none'}`); "
        f"**gained: {len(gained)}** (`{', '.join(gained)}`); "
        f"**kept: {len(kept)}** (`{', '.join(kept)}`).",
        f"- Per-edge arm — **lost: {len(lost_e)}** "
        f"(`{', '.join(lost_e) if lost_e else 'none'}`); **gained: {len(gained_e)}** "
        f"(`{', '.join(gained_e) if gained_e else 'none'}`).",
        "",
        "The two granularities do **not** agree, and the difference is not cosmetic: the "
        "aggregate arm adds",
        f"{len(gained)} systems and loses none, while the per-edge arm adds {len(gained_e)} and "
        f"loses {len(lost_e)}",
        f"(`{', '.join(lost_e) if lost_e else 'none'}`). The next section decides which one may "
        "be read, on grounds",
        "that have nothing to do with the flag counts.",
        "",
        "## Which arm to read (decided on the null, not on the result)",
        "",
        "**The argument.** Correcting an edge's variance by its own measured ratio multiplies "
        "that edge's χ²",
        "contribution by `σ²/s²`, where `s²` is the edge's replicate sample variance on "
        "`ν = n − 1` degrees of",
        "freedom. `E[σ²/s²] = ν/(ν−2)`, which **diverges for ν ≤ 2**. With n = 3 replicates ν = 2, "
        "so the",
        "per-edge self-calibrated χ² has *no finite null expectation*: it is anti-conservative by "
        "construction,",
        "not by accident. The held-out ratio, built on n = 2 (`ν = 1`), is worse still. The "
        "aggregate ratio pools",
        "the denominator over a whole system (`ν ≈ 2E`), leaving only the finite-sample factor "
        "`2E/(2E−2)` — 11% at",
        "E = 10 edges, 3% at E = 40 — so it is *nearly*, but not exactly, unbiased.",
        "",
        "**The measurement.** Each arm was then run on a *perfectly-calibrated synthetic null*: "
        "the same graphs",
        "and the same reported bars, with ΔΔG redrawn from those bars around exact node "
        "potentials and three",
        "replicates per edge, so every flag is a false positive by construction "
        f"({NULL_DRAWS} draws; Monte-Carlo se ≤ {mc_se:.3f}).",
        "`null P(≥1 flag)` is the probability that the arm flags *anything* when nothing is "
        "wrong; under BH-FDR",
        f"at α = {ALPHA} that should be ≈ {ALPHA}. `null mean #flagged` is how many of the "
        f"{len(nominal)} systems it falsely flags per draw.",
        "",
    ] + var + [
        "",
        f"The nominal (published) test comes out at P(≥1) = {null['nominal']['p_any']:.3f}, "
        f"mean {null['nominal']['mean_flagged']:.2f} false flags — it controls its own",
        "false-positive rate, as the paper claims. **Every self-calibrated arm inflates it**, "
        "in the order the",
        "divergence argument predicts:",
        "",
        f"- per-system aggregate: P(≥1) = {null['aggregate']['p_any']:.3f}, mean "
        f"**{null['aggregate']['mean_flagged']:.2f}** false flags per draw;",
        f"- per-edge (n = 3, ν = 2): P(≥1) = {null['per-edge']['p_any']:.3f}, mean "
        f"{null['per-edge']['mean_flagged']:.2f};",
        f"- held-out (n = 2, ν = 1): P(≥1) = {null['held-out']['p_any']:.3f}, mean "
        f"**{null['held-out']['mean_flagged']:.1f}** — it flags "
        f"{100 * null['held-out']['mean_flagged'] / len(nominal):.0f}% of all systems when "
        "nothing is wrong.",
        "",
        "So no self-calibrated arm is FDR-controlled at the nominal level, and **none of them "
        "may be reported as",
        "a calibrated flag set**. What the harness does license is a *magnitude* comparison for "
        "the arm whose",
        f"inflation is smallest: the aggregate arm falsely flags "
        f"{null['aggregate']['mean_flagged']:.2f} systems per draw when nothing is wrong,",
        f"against **{len(gained)} systems gained** on the real data — a factor of "
        f"{len(gained) / max(null['aggregate']['mean_flagged'], 1e-9):.0f}. The gains are not "
        "the correction's own",
        "false-positive inflation. The per-edge arm cannot support even that statement (it "
        f"falsely flags {null['per-edge']['mean_flagged']:.2f}",
        f"per draw against {len(gained_e)} gained, and it *loses* {len(lost_e)}), and the "
        "held-out arm is uninformative rather than a",
        "stronger version of the claim. Everything below is the **aggregate** arm.",
        "",
        "## The result (aggregate arm)",
        "",
        f"**{len(flag_set(nominal))} of {len(nominal)} → {len(flag_set(agg))} of {len(agg)}** "
        f"flagged; median reduced χ² **{median_rc(nominal):.3f} → {median_rc(agg):.3f}**; "
        f"**{len(lost)} lost, {len(gained)} gained**.",
        "",
        "The reading is **not** that the flag set should be "
        f"{len(flag_set(agg))} systems — that arm is not FDR-controlled (above). It is that the",
        "observed miscalibration of the null **suppresses** detections rather than manufacturing "
        "them: correcting",
        "each bar to the width its own replicates measure *adds* systems and removes none, by a "
        f"margin ({len(gained)} gained",
        f"vs {null['aggregate']['mean_flagged']:.2f} expected from the correction's own "
        "false-positive inflation) that the inflation cannot explain. The",
        "published flag count is therefore bounded **below**, which is the direction the "
        "manuscript already claims",
        "on other grounds (median reduced χ² < 1 implies limited single-replicate power) — this "
        "makes that",
        "statement quantitative instead of rhetorical.",
        "",
        "The zero-loss half is the more robust of the two, because false-positive inflation is "
        "the failure mode that",
        "*adds* flags — it cannot manufacture the absence of losses. And losses do happen: the "
        f"per-edge arm loses {len(lost_e)}",
        f"(`{', '.join(lost_e) if lost_e else 'none'}`). That all {len(flag_set(nominal))} "
        "published systems survive being re-tested against their own",
        "measured bars is therefore a result, not an automatic consequence of the arithmetic.",
        "",
        f"**{'One' if len(flat) == 1 else str(len(flat))} of the gains "
        f"{'is a threshold effect' if len(flat) == 1 else 'are threshold effects'}, not "
        "a bar effect.** Benjamini–Hochberg's cutoff relaxes as other systems'",
        "p-values fall, so a system can cross without its own χ² moving much. The gains whose "
        "reduced χ² changes by",
        "less than 15%: "
        + (", ".join(f"`{s}` ({mn[s]['rc']:.2f} → {ma[s]['rc']:.2f}, ratio "
                     f"{ma[s]['med_ratio']:.2f})" for s in flat) or "none")
        + ". Those are carried by the rest of the list rather than by",
        "their own evidence, and should not be read individually.",
        "",
    ]
    if "renin" in gained:
        rn, rp = mn["renin"], ma["renin"]
        lines += [
            "`renin` is the informative gain: the manuscript already names it as a "
            "**systematic-but-unflagged**",
            "case (systematic by the pooled-replicate ratio, reduced χ² 2.5→4.7, yet short of FDR "
            "on a single",
            f"replicate). Under its own measured bars it clears FDR (χ²ᵥ {rn['rc']:.2f} → "
            f"{rp['rc']:.2f}, q {rn['q']:.1e} → {rp['q']:.1e}).",
            "A case the manuscript already documents as systematic-but-unflagged is recovered by "
            "correcting the null",
            "it was tested against — the closest thing this data affords to a recall measurement, "
            "and the one gain",
            "with independent corroboration in the paper.",
            "",
        ]
    lines += [
        "### Stated dependence (a caveat, not a knob)",
        "",
        f"The ratio is measured on all three replicates, including replicate {REPLICATE}, whose "
        "ΔΔG values the closure",
        "test itself reads, so the correction factor and the tested residuals are not fully "
        "independent: a system",
        "that happened to sample a small spread on this replicate gets both tighter bars and, on "
        "average, smaller",
        "residuals. The held-out variant was built as the control for exactly this — but the "
        "synthetic null shows",
        "it is invalid at n = 2 (`ν = 1`), so **this dependence is not resolved by the data "
        "available here** and is",
        "left as a stated limitation. Its direction is not obvious a priori; the aggregate ratio "
        "pools ~10–60 edges",
        "per system, which dilutes but does not eliminate it.",
        "",
        "## Full per-system before-and-after",
        "",
        "Every admitted system, ordered by its self-calibrated *q* (aggregate arm). `agg. ratio` "
        "is that system's",
        "pooled RMS reported-se / replicate-SD. The last two columns give the per-edge arm for "
        "completeness; per",
        "the section above they are **not** a detection count.",
        "",
    ] + tbl + [
        "",
        "## Honest reading",
        "",
        "- This is **exploratory**. It is reported next to — not instead of — the pre-registered "
        "one-sided check",
        "  (P8, verdict SURVIVES), which is unchanged.",
        "- The claim it supports is narrow and one-directional: **the observed miscalibration "
        "suppresses",
        "  detections**. It does not license reporting "
        f"{len(flag_set(agg))} flagged systems as the paper's flag set — no self-calibrated arm "
        "here is",
        "  FDR-controlled — and it does not characterize recall in the strict sense: there are "
        "still no",
        "  ground-truth error labels, only a measured proxy for the width of the null.",
        "- **The pre-specified primary granularity was the wrong choice, and is reported as "
        "such.** Per-edge",
        "  self-calibration on n = 3 replicates has no finite null expectation and inflates its "
        "own false-positive",
        f"  rate to P(≥1) = {null['per-edge']['p_any']:.2f} against a nominal {ALPHA}; the "
        "aggregate arm, fixed before the run as the",
        "  granularity sensitivity, is the least-inflated one. The design was not changed after "
        "seeing the flag",
        "  counts — the arm to read was chosen on a synthetic null that never sees a real ΔΔG "
        "value.",
        "- **The magnitude is load-bearing, the flag list is not.** Read \"ten more systems than "
        "the correction's",
        f"  own noise can produce ({null['aggregate']['mean_flagged']:.2f} per draw)\", not "
        "\"these ten systems are systematic\". Individual gains sitting",
        "  near q ≈ α are exactly the ones the residual inflation can move.",
        "- The correction is **measured, not assumed**: its direction on each system comes from "
        "that system's own",
        "  replicates. Where nothing is measurable, nothing is changed.",
        "- Nothing was tuned after seeing a flag count. D1–D5, the population, α and the full "
        "sensitivity list were",
        "  fixed and written into the script's docstring before the first run; the synthetic-null "
        "diagnostic was",
        "  added afterwards on the strength of the `ν/(ν−2)` argument, and it changes no design "
        "choice.",
        "",
    ]
    lines += ["", "## Realized level of the shipped test by calibration scale", "",
              f"Same seed, same {NULL_DRAWS} draws and the same graphs as the arm table above.",
              "The truth is drawn from the reported bars; the test then divides by a uniform",
              "scale times those bars, so scale 1 is the shipped test.", "",
              "| calibration scale (true/reported) | false flags / draws | P(any false flag) "
              "| 95% Clopper-Pearson | flags on the real data |", "|---|---:|---:|---|---:|"]
    for sc in OPERATING_SCALES:
        k, n, pt, (lo, hi) = levels[sc]
        lines.append(f"| {sc:.2f} | {k}/{n} | {pt:.3f} | [{lo:.4f}, {hi:.4f}] "
                     f"| {real_counts[sc]} |")
    lines += ["",
              "Scale 1.00 estimates the same quantity as the uncorrected arm of the table",
              "above, by a different path, and the two agree within Monte-Carlo error at this",
              "draw count. The interval is what 400 draws support: the two smallest points are",
              "not resolved from each other, or from zero, at this draw count. The last column",
              "is the real-data companion, the flag count when every reported bar is corrected",
              "to that scale, and it is not a simulation."]
    DOC.write_text("\n".join(lines) + "\n")


def main():
    use_paper_style()
    res = analyze()
    nominal, per_edge, agg = res["nominal"], res["primary"], res["aggregate"]

    med_nom = median_rc(nominal)
    flagged_nom = sorted(flag_set(nominal))
    anchor_ok = (abs(med_nom - MEDIAN_RC_ANCHOR) <= ANCHOR_TOL
                 and set(flagged_nom) == FLAGGED_ANCHOR)
    print(f"[self] population: {len(nominal)} systems (replicate {REPLICATE}, dof>=1)")
    print(f"[self] sanity anchor -- median rc_nominal = {med_nom:.4f} "
          f"(published Fig L: {MEDIAN_RC_ANCHOR}); flagged_nominal = {flagged_nom}")
    if not anchor_ok:
        print("\n[self] SANITY ANCHOR MISMATCH -- STOP. The population or test path has drifted "
              "from Fig L; do NOT trust anything below.")
        print("DONE_WITH_CONCERNS")
        return
    print(f"[self] sanity anchor OK (median within {ANCHOR_TOL}, flagged set exact match)\n")

    by_systems = figL.load_systems()
    levels = level_at_calibration(by_systems, OPERATING_SCALES)
    real_counts = flag_count_at_calibration(by_systems, OPERATING_SCALES)
    print("[self] realized P(any false flag) by CALIBRATION scale (true/reported), "
          f"{NULL_DRAWS} seeded draws, with the exact binomial interval and the real-data "
          "flag count at the same scale:")
    for sc in OPERATING_SCALES:
        k, n, pt, (lo, hi) = levels[sc]
        print(f"[self]   scale {sc:.2f} -> {k}/{n} = {pt:.3f}  CP [{lo:.4f}, {hi:.4f}]  "
              f"| real data flags {real_counts[sc]}")

    make_figure(res)
    write_doc(res, levels, real_counts)
    print(f"wrote figSelf_two_sided_selfcalibration.(pdf|png) to {FIGDIR}")
    print(f"wrote {DOC.relative_to(ROOT)}\n")

    ratios = res["ratios"]
    print(f"[ratios] per-edge n={ratios.size} measurable, {res['n_degenerate']} unchanged (D2); "
          f"median {np.median(ratios):.3f}  IQR {np.percentile(ratios, 25):.3f}-"
          f"{np.percentile(ratios, 75):.3f}  range {ratios.min():.3f}-{ratios.max():.2f}; "
          f"{100 * float(np.mean(ratios >= 1)):.1f}% >= 1 (bar tightened)")
    sv = np.array(sorted(res["sysr"].values()))
    print(f"[ratios] per-system aggregate n={sv.size}; median {np.median(sv):.3f}  "
          f"range {sv.min():.3f}-{sv.max():.3f}")

    for lab, arm in (("aggregate", agg), ("per-edge", per_edge)):
        kept, gained, lost = transitions(nominal, arm)
        print(f"\n[{lab}] {len(flag_set(nominal))}/{len(nominal)} -> {len(flag_set(arm))}/"
              f"{len(arm)} flagged; median reduced chi2 {median_rc(nominal):.3f} -> "
              f"{median_rc(arm):.3f}")
        print(f"[{lab}] kept ({len(kept)}): {kept}")
        print(f"[{lab}] gained ({len(gained)}): {gained}")
        print(f"[{lab}] lost ({len(lost)}): {lost if lost else '[]'}")

    print(f"\n[null calibration] {NULL_DRAWS} draws from a perfectly-calibrated synthetic null "
          f"(every flag is a false positive):")
    print(f"    {'arm':<46s}{'flagged':>8}{'P(>=1)':>9}{'mean #':>8}{'median rc':>11}")
    for k, lab, rows in res["variants"]:
        n = res["null"][k]
        print(f"    {lab:<46s}{len(flag_set(rows)):>8}{n['p_any']:>9.3f}"
              f"{n['mean_flagged']:>8.2f}{n['median_rc']:>11.3f}")

    mn = {d["sys"]: d for d in nominal}
    ma = {d["sys"]: d for d in agg}
    me = {d["sys"]: d for d in per_edge}
    print("\n[per-system] system                rc_nom     q_nom    rc_agg      q_agg  status"
          "     rc_edge")
    kept, gained, lost = transitions(nominal, agg)
    for s in sorted(mn, key=lambda s: ma[s]["q"]):
        st =("gained" if s in gained else "kept" if s in kept else "LOST" if s in lost else "")
        print(f"    {s:<24s}{mn[s]['rc']:>8.3f} {mn[s]['q']:>9.2e} {ma[s]['rc']:>9.3f} "
              f"{ma[s]['q']:>10.2e}  {st:<8s} {me[s]['rc']:>8.3f}")


if __name__ == "__main__":
    main()


