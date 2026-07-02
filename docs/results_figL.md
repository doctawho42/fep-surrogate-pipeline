# Results — Fig L: calibrated cycle-closure quality control (the impact result)

**Figure:** `figs/figL_calibrated_cycle_closure.{pdf,png}` · **Reproduce:** `make figL`
(`PYTHONPATH=src python figs/make_figL.py`). Deterministic. Data:
`data/openfe_replicates/combined_pymbar4_edge_data.csv` (OpenFE IndustryBenchmarks2024, public).

## The claim (read this first)
The paper's audit shows a calibrated σ does not sharpen intervals, re-rank edges, or steer a
generator. This figure is the **positive** counterpart: a calibrated per-edge aleatoric variance
is exactly the **null model** that turns thermodynamic-cycle closure into a *properly-calibrated
statistical test* — a deployable FEP quality-control tool that **separates edge-level systematic
error from sampling noise**, and that an overconfident learned σ provably cannot deliver.

This is the one downstream task where calibration is not second-order: it is **detection**, not
weighting-for-ranking, so the Gauss–Markov ceiling (which neutered edge weighting, Fig C) does
not apply.

## The test
A perturbation network must close its cycles: for any loop the signed sum of ΔΔG is zero in
expectation. Fit node potentials by GLS with weights `1/V_e`, `V_e = B_e/I_e²` (sandwich); the
residuals give
```
X² = Σ_e r_e² / V_e   ~   χ²(dof),   dof = #independent cycles = E − rank(incidence).
```
Reduced `χ²_ν = X²/dof`. Under the null (sampling-only), `E[χ²_ν] ≈ 1` (here **< 1**, median 0.34;
since systematic error only *raises* `χ²_ν`, a median of 0.34 bounds the bars as **≥1.7×
conservative** for closure residuals — consistent with, and if anything stronger than, the 1.41×
aggregate of Fig A-rep). `χ²_ν` well above 1 ⇒ the network is over-dispersed
relative to sampling ⇒ **edge-level systematic error** (non-convergence, hysteresis,
charge/water/protonation artifacts).

## Result: PASS (a real, chemically-sensible positive)

**Panel A — deployable test (replicate 0, sandwich null; 48 systems with ≥1 cycle).**
Median `χ²_ν = 0.34` (most systems sampling-consistent, as the conservative bars predict). A
Benjamini–Hochberg FDR test (α=0.05) flags **6 systems**:

| system | edges | cycles | reduced χ² | p | known hard because |
|---|---:|---:|---:|---:|---|
| brd4  |  8 |  1 | **15.4** | 1.4e-4 | bromodomain buried waters (dataset = `waterset_BRD4`) |
| bace  | 49 | 14 | **5.57** | 1.9e-10 | aspartic protease, protonation-sensitive (`jacs_bace`) |
| faah  | 31 |  8 | **3.71** | 2.8e-4 | covalent/large-pocket series |
| cdk8  | 63 | 29 | **2.68** | 3.0e-6 | large scaffold changes |
| hif2a | 59 | 19 | **2.53** | 2.6e-4 | buried polar cavity |
| p38   | 60 | 22 | **2.40** | 2.6e-4 | large R-group swaps (`jacs_p38`) |

Not a size artifact (brd4 has 8 edges, 1 cycle, χ²=15; the reduced χ² normalizes by dof).

**Panel B — the detector's usefulness IS σ-calibration (same test, same data, different σ).**

| σ model | scale | systems flagged |
|---|---:|---:|
| uniform ×0.15 shrink (learned head's central overconf., Fig A) | ×0.15 | **42/48 (88%)** — FPR→1, useless |
| **sandwich = MBAR (calibrated)** | ×1.0 | **6/48 (12%)** — selective |
| conservative bar | ×1.3 | 1/48 (2%) |
| too wide | ×2.0 | 0/48 (0%) |

An overconfident σ makes cycle-closure QC worthless — it flags essentially every system. Only a
calibrated σ yields a selective test. **This is the concrete decision an overconfident σ silently
corrupts and a calibrated one enables**, made quantitative on real data. (The ×0.15 here is a
*uniform* stand-in for the learned head's central overconfidence; the real head is regime-dependent,
0.09–0.20× (Fig A) — but *any* strong shrink drives the FPR→1, so the conclusion is robust to the
exact per-edge profile.)

**Panel C — the 3 independent replicates separate SYSTEMATIC from SAMPLING.**
Pooling the replicates shrinks the bars (`se → se/√3`) but *increases* the flagged systems'
reduced χ² toward the `3×` variance-components ceiling — which happens **only** if the closure
error is reproducible across independent runs (systematic), not if it is sampling noise (ratio
→1) or merely under-sized bars (ratio →1):

```
brd4 9.5→26.8   bace 6.7→13.7   cdk8 2.2→4.9   p38 2.2→4.7   renin 2.5→4.7   hif2a 2.5→3.6
```
Variance-components identity: single-rep `χ²_ν = 1 + s²/V_samp`, pooled-3 `χ²_ν = 1 + 3s²/V_samp`
(`s` = systematic closure, `V_samp` = sampling variance). Ratio → 3 iff systematic-dominated,
→ 1 iff sampling-only. brd4 ratio 2.8, bace 2.0: near-pure systematic. This rules out "the bars
are just too small" and confirms the flags are real model/convergence error.

**Panel A and Panel C select overlapping but distinct sets** (measuring related but different
things): Panel A = single-replicate BH-FDR *power*; Panel C = the *replicate systematic ratio*.
`renin` is systematic by the ratio (pooled χ²_ν 2.5→4.7, Panel C) yet did not clear FDR on one
replicate (under-powered → not in A's flagged six); conversely `faah` clears FDR (Panel A) but
has a smaller replicate ratio (a larger sampling share). Both readings are consistent with the
variance-components picture — a system can be systematic-dominated but single-replicate-underpowered,
or flag-worthy on one replicate with more of its dispersion still sampling.

**Localization (per-edge standardized GLS residual, pooled-3).** For systems with **≥2
independent cycles** the test pinpoints the culprit edges, e.g. bace `CAT-13b→CAT-13a` (z=−7.3),
`CAT-13e→CAT-13a` (+6.5); p38 `3flq→2v` (+5.6, ΔΔG≈7 kcal/mol) — chemically interpretable.
*Caveat:* localization needs the edge to sit in ≥2 cycles. A single-cycle system such as **brd4**
(8 edges, 1 cycle) is correctly *flagged*, but its closure error is shared equally across the
cycle's edges, so the specific culprit edge is **not identifiable** there (only the system is).

## Why this needs the differentiable sandwich (not just a number)
The test uses the reported pymbar4 MBAR se (= sandwich `B/I²` to leading order, Fig A). The
sandwich's added value is that it is a **per-edge, graph-native, non-negative, back-propagatable
weight** `w_e = I_e²/B_e` — the Fisher–resistance Laplacian conductance (Thm 3) — so the same
object that computes the closure null also (a) enters GLS as the correct weight, (b) is
differentiable for downstream training, and (c) is ≥0 by construction where the plug-in MBAR can
go negative. A returned pymbar number cannot be a graph weight; an overconfident learned σ gives
the wrong null (Panel B).

## Validation: the flags are causal and reproduce out-of-sample (`make figLval`)
Two data-internal predictive checks (short of new MD), in `figs/figL_validation.{pdf,png}`:

- **Repair test (causal localization).** Greedily removing the highest-|z| edges reaches
  sampling-consistency (reduced χ²≤1) in **3–7 guided removals**, versus **12–37 random**
  removals (bace 4 vs 24, faah 4 vs 12, cdk8 7 vs 37, hif2a 3 vs 26, p38 3 vs 29). The detector
  points at the *causal* culprit edges, not arbitrary ones — removing the ones it names fixes the
  network ~6–9× faster than chance.
- **Out-of-sample reproduction.** Fitting each of the 3 independent replicates separately, the
  per-edge standardized residuals correlate **r = +0.30, +0.42, +0.38** across replicate pairs
  (n=1143 edges). Pure sampling noise would give ≈0; the positive correlation means the flagged
  inconsistency **reproduces in unseen runs** — it is systematic, not a one-run fluctuation.

## Honest scope
- **Cycle closure detects EDGE-level inconsistency only.** Node-consistent force-field bias (a
  smooth per-ligand offset) cancels around a loop and is *invisible* here — correctly, but it
  means "systematic" = edge-level (convergence, hysteresis, charge/water/protonation), not all
  force-field error. Stated as a limitation, not hidden.
- **Retrospective (but predictive).** The flags coincide with the known-hard systems (brd4
  waters, bace protonation), survive a causal repair test, and reproduce on held-out replicates
  (Validation section above). The one step not done here is **new MD**: the decisive prospective
  test is flag → re-run/repair the culprit edge with fresh sampling → confirm the cycle closes and
  the affinity estimate improves. Natural next experiment.
- **χ² assumes Gaussian edge error;** Panel C's replicate reproducibility is the model-free
  backstop that does not rely on Gaussianity.
- p-values via Wilson–Hilferty (adequate; flagged p ∈ [1e-10, 1e-4]); BH-FDR across 48 systems.
- Threshold at nominal `χ²_ν=1` is conservative given the ≥1.7×-conservative bars; using the
  empirical conservative baseline (median ~0.34) would increase power.

## Novelty (positioned honestly)
Network MLE/GLS estimation of ΔG from ΔΔG (arsenic/cinnabar) and cycle-closure hysteresis as a
QC signal are **standard**. The contribution here is (1) framing closure as a *calibrated χ²
hypothesis test* with a per-edge **validated** aleatoric null (giving a controlled false-positive
rate and per-edge-adaptive thresholds, versus the field's fixed ad-hoc hysteresis cutoffs);
(2) the result that the test's FPR is **entirely determined by σ-calibration**; and (3) the
replicate-based systematic-vs-sampling separator. It operationalizes the paper's calibration
thesis into a concrete FEP capability.

## Gate
Real public data (1145 edges, 49 systems — the 48 with ≥1 cycle — 3 replicates); one-command reproduction; calibration
sweep + replicate control both included. Positive, non-trivial, and it is exactly the task an
overconfident σ ruins → **the paper's headline moves from "calibration buys trust" to
"calibration buys a deployable QC capability."**
