"""Fig P8 (SI table, no figure file) -- per-system flag robustness to each system's
OWN measured calibration error (peer-review response, ml-uq R1-4 / P8: the "mild
circularity" concern).

The Fig L cycle-closure QC test flags 6 of 48 systems using the calibrated sandwich
null. But the per-system calibration of that null is itself heterogeneous
(Fig A-rep / docs/results_figA_replicates.md): four of the six flagged systems are
themselves LOCALLY overconfident (their own reported se runs a bit tight against their
own across-replicate SD), and the one markedly overconfident protonation-variant target
in that check (bace_p3_arg368_in) sits inside the flagged `bace` system. A referee
called this a mild circularity: a system could be flagged partly because its own error
bars are too tight, not because its physics is inconsistent. This script asks the
direct question -- give each flagged (and every other) system the benefit of its own
measured local miscalibration and see whether the flags survive.

REUSED RATIO DEFINITION (do not re-derive): the per-system ratio here is computed by
calling the EXACT same formula as `figs/make_figA_replicates.py`'s per-target pooled
ratio (its `main()`, the `tratio` list feeding Fig A-rep panel B):

    ratio_s = sqrt(mean(rep[sys==s]**2)) / sqrt(mean(repl[sys==s]**2))

where, for every edge with a non-degenerate replicate spread (`repl > 1e-6`):
  - `rep`  = that edge's reported se = sqrt(complex_dDG**2 + solvent_dDG**2), RMS'd
             over the 3 independent OpenFE replicates via `make_figA_replicates.load()`;
  - `repl` = the across-replicate empirical SD of that edge's binding ddG (the 3
             independent-replicate "truth").
This is an RMS ratio (root-mean-square of reported se over RMS of replicate SD), NOT a
naive ratio of per-edge means -- the two differ modestly (see docs/results_figP8.md for
the numeric contrast); the paper's Fig A-rep / this table report the RMS version. We
call `make_figA_replicates.load()` directly and reduce it with `figA_replicates`'s own
per-target formula (unmodified), then look up each Fig-L system's ratio by the shared
`system name` key -- the same CSV field both scripts group by, so the systems line up
exactly (verified: `bace` here is 49 edges, matching Fig L's `bace` row).

The cycle-closure test itself is `bar.qc.gls_network` / `bar.qc.chi2_sf` /
`bar.qc.benjamini_hochberg`, called through `figs/make_figL.py`'s own thin adapters
(`load_systems`, `edge_val`, `gls_chi2`) so the nominal numbers are byte-for-byte Fig
L's numbers, not a re-derivation.

FROZEN PRE-REGISTRATION (recorded before this run; not tuned after seeing numbers):
  - Population: replicate 0, systems with dof >= 1 (Fig L's own population), BH-FDR
    alpha = 0.05.
  - Self-calibration inflates a system's per-edge se by 1/ratio ONLY where ratio < 1
    (give the system the benefit of its own measured local miscalibration); systems
    with ratio >= 1 are already conservative and are left unchanged.
  - SURVIVES iff every system flagged under the nominal null is still flagged under the
    self-calibrated null.
  - PARTIAL iff at least one but not all flagged systems drop out.
  - FAILS iff no flagged system survives.
  - All three outcomes are reportable; the inflation rule, population, and alpha are not
    adjusted in response to any number produced by this script.

MANDATORY SANITY ANCHOR: rc_nominal / flag_nominal must reproduce Fig L's published
values -- median reduced chi^2 0.34, flagged set {bace, brd4, cdk8, faah, hif2a, p38}.
If not, this script prints a STOP banner instead of a verdict.

Run: PYTHONPATH=src python figs/make_figP8.py   (or `make figP8`). Deterministic.
No new MD; no figure file (this deliverable is a table).
"""
from __future__ import annotations

import pathlib
import sys

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "figs"))

from bar.qc import benjamini_hochberg as bh_flags, chi2_sf  # noqa: E402

import make_figA_replicates as figArep  # noqa: E402
import make_figL as figL  # noqa: E402

FLAGGED_ANCHOR = frozenset({"bace", "brd4", "cdk8", "faah", "hif2a", "p38"})
MEDIAN_RC_ANCHOR = 0.34
ANCHOR_TOL = 0.01  # median reduced-chi2 tolerance


def compute_ratios() -> dict[str, float]:
    """Per-system RMS ratio, IDENTICAL formula to make_figA_replicates.py's per-target
    pooled ratio (its `main()`, `tratio` list) -- reused here without that function's
    display-only >=8-edge filter (we want every system Fig L's own population needs,
    including brd4 at 8 edges which clears it anyway)."""
    rows = figArep.load()
    rep = np.array([r[0] for r in rows])
    repl = np.array([r[1] for r in rows])
    sysn = np.array([r[3] for r in rows])
    keep = repl > 1e-6  # non-degenerate replicate spread, same filter as figArep.main()
    rep, repl, sysn = rep[keep], repl[keep], sysn[keep]
    ratios = {}
    for s in sorted(set(sysn)):
        m = sysn == s
        if m.sum() == 0:
            continue
        ratios[s] = float(np.sqrt(np.mean(rep[m] ** 2)) / np.sqrt(np.mean(repl[m] ** 2)))
    return ratios


def analyze():
    by = figL.load_systems()
    ratios = compute_ratios()
    rows = []
    missing_ratio = []
    for sys_name, rs in sorted(by.items()):
        e0 = [(r["ligand_A"], r["ligand_B"], *figL.edge_val(r, 0)) for r in rs if figL.edge_val(r, 0)]
        if len(e0) < 3:
            continue
        X2, dof, _ = figL.gls_chi2(e0)
        if dof < 1:
            continue

        ratio = ratios.get(sys_name)
        rc_nominal = X2 / dof
        p_nominal = chi2_sf(X2, dof)

        if ratio is None:
            missing_ratio.append(sys_name)
            inflate = 1.0
            note = "no ratio (no non-degenerate-spread edges) -- left unchanged"
        elif ratio < 1.0:
            inflate = 1.0 / ratio
            note = f"self-cal: se x {inflate:.3f} (1/ratio, ratio<1)"
        else:
            inflate = 1.0
            note = "ratio>=1 (already conservative) -- unchanged"

        e0_sc = [(a, b, ddg, se * inflate) for (a, b, ddg, se) in e0]
        X2_sc, dof_sc, _ = figL.gls_chi2(e0_sc)
        rc_selfcal = X2_sc / dof_sc
        p_selfcal = chi2_sf(X2_sc, dof_sc)

        rows.append(dict(
            sys=sys_name, E=len(e0), dof=dof, ratio=ratio, inflate=inflate, note=note,
            rc_nominal=rc_nominal, p_nominal=p_nominal,
            rc_selfcal=rc_selfcal, p_selfcal=p_selfcal,
        ))

    flags_nom = bh_flags([r["p_nominal"] for r in rows])
    flags_sc = bh_flags([r["p_selfcal"] for r in rows])
    for r, fn, fs in zip(rows, flags_nom, flags_sc):
        r["flag_nominal"] = bool(fn)
        r["flag_selfcal"] = bool(fs)
    return rows, missing_ratio


def main():
    rows, missing_ratio = analyze()

    # ---- mandatory sanity anchor ----
    med_nom = float(np.median([r["rc_nominal"] for r in rows]))
    flagged_nom = sorted(r["sys"] for r in rows if r["flag_nominal"])
    anchor_ok = (abs(med_nom - MEDIAN_RC_ANCHOR) <= ANCHOR_TOL
                 and set(flagged_nom) == FLAGGED_ANCHOR)

    print(f"[P8] population: {len(rows)} systems (replicate 0, dof>=1)")
    if missing_ratio:
        print(f"[P8] systems with no ratio (no non-degenerate-spread edges), left unchanged: "
              f"{missing_ratio}")
    print(f"[P8] sanity anchor -- median rc_nominal = {med_nom:.4f} "
          f"(published Fig L: {MEDIAN_RC_ANCHOR}); "
          f"flagged_nominal = {flagged_nom}")
    print(f"[P8] sanity anchor -- published Fig L flagged set = {sorted(FLAGGED_ANCHOR)}")
    if not anchor_ok:
        print("\n[P8] SANITY ANCHOR MISMATCH -- STOP. The population or test path has "
              "drifted from Fig L; do NOT trust the table below.")
        print("DONE_WITH_CONCERNS")
        return

    print(f"[P8] sanity anchor OK (median matches within {ANCHOR_TOL}, flagged set exact match)\n")

    # ---- full per-system table ----
    header = (f"{'system':<10}{'E':>4}{'dof':>5}{'ratio':>8}{'rc_nom':>9}{'flag_nom':>10}"
              f"{'rc_selfcal':>12}{'flag_selfcal':>13}  note")
    print(header)
    print("-" * len(header))
    for r in sorted(rows, key=lambda d: (not d["flag_nominal"], d["rc_nominal"])):
        ratio_s = f"{r['ratio']:.3f}" if r["ratio"] is not None else "n/a"
        print(f"{r['sys']:<10}{r['E']:>4}{r['dof']:>5}{ratio_s:>8}{r['rc_nominal']:>9.3f}"
              f"{str(r['flag_nominal']):>10}{r['rc_selfcal']:>12.3f}"
              f"{str(r['flag_selfcal']):>13}  {r['note']}")

    # ---- pre-registered readout: flagged systems only ----
    flagged_rows = [r for r in rows if r["flag_nominal"]]
    print("\n[P8] pre-registered readout -- systems flagged under the NOMINAL null:")
    survivors, dropped = [], []
    for r in sorted(flagged_rows, key=lambda d: d["sys"]):
        status = "STILL FLAGGED" if r["flag_selfcal"] else "DROPPED"
        (survivors if r["flag_selfcal"] else dropped).append(r["sys"])
        ratio_s = f"{r['ratio']:.3f}" if r["ratio"] is not None else "n/a"
        print(f"    {r['sys']:<8} ratio={ratio_s:<7} rc_nominal={r['rc_nominal']:<8.3f}"
              f" rc_selfcal={r['rc_selfcal']:<8.3f} -> {status}")

    n_flag = len(flagged_rows)
    if len(survivors) == n_flag:
        verdict = "SURVIVES"
    elif len(survivors) == 0:
        verdict = "FAILS"
    else:
        verdict = "PARTIAL"

    print(f"\n[P8] VERDICT: {verdict}  ({len(survivors)}/{n_flag} nominally-flagged systems "
          f"still flagged after self-calibration)")
    if verdict == "PARTIAL":
        print(f"[P8] dropped systems (circularity conceded for these): {sorted(dropped)}")
        print(f"[P8] surviving systems: {sorted(survivors)}")
    print("[P8] Note: the two distribution-free backstops (causal edge removal, Fig Lval; "
          "cross-replicate residual correlation, Fig Lval) do not depend on this null at all "
          "and are unaffected by this concern in either direction.")


if __name__ == "__main__":
    main()
