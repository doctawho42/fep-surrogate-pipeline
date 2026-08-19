# Results — Fig Stab: replicate stability of the cycle-closure flag set

**Figure:** `figs/figStab_replicate_stability.{pdf,png}` · **Reproduce:** `make figStab`
(`PYTHONPATH=src python figs/make_figStab.py`). Deterministic (closed-form random-draw
reference; no Monte-Carlo, no seed). Data:
`data/openfe_replicates/combined_pymbar4_edge_data.csv` (OpenFE IndustryBenchmarks2024,
public, 3 independent replicates per edge).

## What this asks
The Fig L headline flag set is computed on **replicate 0 only**. This figure re-runs the
*identical* test — same edge construction, same `bar.qc.gls_network` / `chi2_sf` /
`benjamini_hochberg`, same pre-registered α = 0.05 — on replicates 1 and 2, and on the
pooled network. Nothing is tuned: every constant is the one already in `make_figL.py`.

## Headline: the flag set is NOT stable across replicates

| replicate | systems tested | flagged |
|---|---:|---|
| 0 (the paper) | 48 | 6: `bace, brd4, cdk8, faah, hif2a, p38` |
| 1 | 48 | 6: `bace, brd4, ciordia_retro, hif2a, p38, renin` |
| 2 | 48 | 3: `bace, cdk8, renin` |
| pooled (≥2 reps) | 48 | 7: `bace, brd4, cdk8, ciordia_retro, hif2a, p38, renin` |

### Under Benjamini--Yekutieli (arbitrary dependence)

| replicate | flagged | vs BH |
|---|---|---|
| 0 | 6: `bace, brd4, cdk8, faah, hif2a, p38` | identical to BH |
| 1 | 5: `bace, ciordia_retro, hif2a, p38, renin` | drops `brd4` |
| 2 | 3: `bace, cdk8, renin` | identical to BH |

- **Intersection (1 system):** `bace`.
- **Union (8 systems):** `bace, brd4, cdk8, ciordia_retro, faah, hif2a, p38, renin`.
- **Pairwise Jaccard:** 0↔1 0.500, 0↔2 0.286, 1↔2 0.286 (mean 0.357).
- **Random-draw reference:** the exact expected Jaccard of two *independent uniform* draws of the observed sizes from the same 48 systems is 0.072, 0.048, 0.048 (mean 0.056); for two 6-of-48 draws it is 0.072.

So the observed overlap is well above chance — the test is not returning noise — but it
is far below the reproducibility a set-valued claim needs: only
`bace` survives all three runs, and 7 of the 8 ever-flagged systems appear in some runs and not others.

## Where the flips come from
Panel B: the flips are **threshold crossings**, not sign changes. A system's BH-adjusted
*q* swings by orders of magnitude between independent runs of the same protocol, so any
system whose *q* lives near α will cross it about half the time. The clearest case is
`faah`, one of the paper's six: reduced χ² = 3.71, 0.45, 1.47 across the three replicates — a swing of 8.3× that crosses the nominal χ²ᵥ = 1 threshold in both directions.

## The pooled network (the paper's own pooling rule)
Fig L panel C defines the pooling rule and it is reused here **verbatim**: per edge, the
mean ΔΔG over the available replicates and se = √(mean_k se_k² / n) — the rule that gives
panel C's `se → se/√3` — keeping edges with ≥ 2 available replicates.

**Stated assumption.** Panel C does *not* define a flagging rule on the pooled network —
it only compares reduced χ² single-vs-pooled. Running BH-FDR on the pooled p-values is
therefore this figure's extension; it uses the same pre-registered α = 0.05.
The one free knob in the pooling rule is how many replicates an edge must have. Both
settings are reported rather than chosen:

- `min_reps = 2` (panel C's own setting): 7 flagged — `bace, brd4, cdk8, ciordia_retro, hif2a, p38, renin`
- `min_reps = 3` (every edge complete in all three runs): 7 flagged — `bace, brd4, cdk8, ciordia_retro, hif2a, p38, renin`

The two agree exactly, so the pooled result does not depend on that choice.

The pooled network flags **7 of the 8 ever-flagged systems** — everything in the union except `faah`. Pooling triples the effective information per edge, so this is the highest-powered
single read available from this data; it is *not* an independent confirmation of any
single-replicate set (it re-uses all three runs).

## Which flags reproduce, model-free
The distribution-free counterpart of the χ² test: per system, the mean pairwise Pearson
correlation of the per-edge **signed** standardized residuals across the three
replicates, on edges complete in all three. Sampling noise → 0; a reproducible
edge-level systematic error → positive.

| system | flagged in reps | cross-replicate signed-z reproduction | aligned edges |
|---|---|---:|---:|
| brd4 | 0, 1 | +0.995 | 8 |
| cdk8 | 0, 2 | +0.623 | 63 |
| p38 | 0, 1 | +0.555 | 60 |
| bace | 0, 1, 2 | +0.553 | 49 |
| ciordia_retro | 1 | +0.518 | 45 |
| renin | 1, 2 | +0.467 | 42 |
| hif2a | 0, 1 | +0.206 | 59 |
| faah | 0 | -0.085 | 31 |

`faah` is the outlier: its closure error is essentially **not** the same error run to
run (-0.085), which is the mechanism behind its χ² swing above, and it is the
only union member the pooled test does not flag.

## What does NOT move

**1. The per-system ranking.** Spearman ρ of the reduced χ² across all 48 systems: (0,1) ρ = 0.65, (0,2) ρ = 0.64, (1,2) ρ = 0.72. The same networks are the worst-closing ones every run; it is the *dichotomisation*
at α, not the underlying signal, that is unstable.

**2. The selectivity-vs-σ-calibration curve** (Fig L panel B, per replicate):

| σ scale | replicate 0 | replicate 1 | replicate 2 |
|---|---:|---:|---:|
| ×0.15 | 42/48 (88%) | 41/48 (85%) | 41/48 (85%) |
| ×1 | 6/48 (12%) | 6/48 (12%) | 3/48 (6%) |
| ×1.3 | 1/48 (2%) | 1/48 (2%) | 1/48 (2%) |
| ×2 | 0/48 (0%) | 0/48 (0%) | 0/48 (0%) |

The panel-B conclusion — an overconfident σ drives the false-positive rate toward 1, a
calibrated σ is selective, an over-wide σ loses all power — reproduces on every
replicate. That claim is about the *detector*, and it does not depend on which systems
come out.

## Robustness to the FDR level
The instability is not an artefact of α = 0.05:

| α | set sizes (rep 0/1/2) | intersection | union | pairwise Jaccard |
|---:|---|---:|---:|---|
| 0.01 | 6/5/3 | 1 | 8 | 0.375 / 0.286 / 0.333 |
| 0.05 | 6/6/3 | 1 | 8 | 0.500 / 0.286 / 0.286 |
| 0.1 | 7/7/4 | 2 | 9 | 0.556 / 0.375 / 0.375 |
| 0.2 | 7/7/8 | 4 | 11 | 0.556 / 0.500 / 0.500 |

At every level tested the three replicates disagree on a majority of the union.

## Full per-system table
Reduced χ² / BH-adjusted *q* per replicate and pooled; **bold** = flagged at α = 0.05.
Systems are ordered by their strongest (smallest) *q* anywhere; edges and cycles are
replicate-0 counts (they vary by ≤ a few edges between replicates where a leg failed).

| system | edges | cycles | rep 0 χ²ᵥ / q | rep 1 χ²ᵥ / q | rep 2 χ²ᵥ / q | pooled χ²ᵥ / q | flagged in (reps) |
|---|---:|---:|---:|---:|---:|---:|---:|
| **bace** | 49 | 14 | **5.57 / 3.2e-09** | **7.33 / 7.0e-14** | **7.10 / 2.9e-13** | **13.68 / 1.4e-31** | 3/3 |
| **cdk8** | 63 | 29 | **2.68 / 6.1e-05** | 1.45 / 2.8e-01 | **2.62 / 1.1e-04** | **4.90 / 1.7e-15** | 2/3 |
| **p38** | 60 | 22 | **2.40 / 2.0e-03** | **2.48 / 1.7e-03** | 1.81 / 1.1e-01 | **4.70 / 2.7e-11** | 2/3 |
| **renin** | 42 | 14 | 1.29 / 8.1e-01 | **2.81 / 3.1e-03** | **3.32 / 3.8e-04** | **4.67 / 1.5e-07** | 2/3 |
| **hif2a** | 59 | 19 | **2.53 / 2.0e-03** | **2.96 / 2.4e-04** | 2.05 / 5.4e-02 | **3.61 / 1.6e-06** | 2/3 |
| **brd4** | 8 | 1 | **15.40 / 1.4e-03** | **8.05 / 3.6e-02** | 4.94 / 1.7e-01 | **26.81 / 1.8e-06** | 2/3 |
| **ciordia_retro** | 45 | 14 | 0.89 / 1.0e+00 | **3.86 / 3.1e-05** | 0.40 / 1.0e+00 | **3.12 / 4.5e-04** | 1/3 |
| **faah** | 31 | 8 | **3.71 / 2.0e-03** | 0.45 / 1.0e+00 | 1.47 / 7.1e-01 | 2.17 / 1.2e-01 | 1/3 |
| mcl1 | 76 | 24 | 1.82 / 5.5e-02 | 1.79 / 6.8e-02 | 1.04 / 1.0e+00 | 1.38 / 3.8e-01 | 0/3 |
| liga | 13 | 3 | 2.10 / 4.7e-01 | 0.76 / 1.0e+00 | 1.80 / 6.9e-01 | 3.41 / 9.4e-02 | 0/3 |
| tnks2 | 35 | 9 | 0.63 / 1.0e+00 | 0.43 / 1.0e+00 | 2.07 / 1.7e-01 | 2.21 / 9.4e-02 | 0/3 |
| tyk2 | 29 | 11 | 0.79 / 1.0e+00 | 1.06 / 1.0e+00 | 0.78 / 1.0e+00 | 2.06 / 9.4e-02 | 0/3 |
| bace_p3_arg368_in | 28 | 8 | 1.88 / 3.5e-01 | 1.88 / 2.8e-01 | 2.32 / 1.4e-01 | 1.03 / 9.7e-01 | 0/3 |
| mup1 | 7 | 2 | 0.55 / 1.0e+00 | 3.05 / 2.8e-01 | 0.29 / 1.0e+00 | 0.85 / 9.7e-01 | 0/3 |
| scyt_dehyd | 7 | 1 | 0.54 / 1.0e+00 | 1.34 / 8.5e-01 | 1.23 / 9.3e-01 | 3.02 / 3.3e-01 | 0/3 |
| taf12 | 8 | 1 | 3.21 / 3.9e-01 | 0.15 / 1.0e+00 | 0.33 / 1.0e+00 | 1.19 / 7.4e-01 | 0/3 |
| cdk2 | 27 | 10 | 0.48 / 1.0e+00 | 1.33 / 7.7e-01 | 1.22 / 9.3e-01 | 1.55 / 4.0e-01 | 0/3 |
| btk | 8 | 3 | 0.19 / 1.0e+00 | 2.15 / 4.0e-01 | 1.19 / 9.3e-01 | 0.83 / 1.0e+00 | 0/3 |
| t4_lysozyme | 14 | 3 | 0.25 / 1.0e+00 | 0.57 / 1.0e+00 | 2.14 / 5.0e-01 | 1.55 / 6.0e-01 | 0/3 |
| thrombin | 54 | 19 | 1.30 / 7.4e-01 | 0.86 / 1.0e+00 | 0.61 / 1.0e+00 | 1.31 / 5.2e-01 | 0/3 |
| bace1 | 3 | 1 | 0.05 / 1.0e+00 | 1.09 / 8.9e-01 | 0.43 / 1.0e+00 | 1.20 / 7.4e-01 | 0/3 |
| ephx2 | 4 | 1 | 0.77 / 1.0e+00 | 1.67 / 7.7e-01 | 1.33 / 9.3e-01 | 0.35 / 1.0e+00 | 0/3 |
| galectin | 36 | 11 | 0.37 / 1.0e+00 | 0.22 / 1.0e+00 | 0.88 / 1.0e+00 | 1.16 / 7.8e-01 | 0/3 |
| shp2 | 37 | 12 | 0.74 / 1.0e+00 | 1.21 / 8.7e-01 | 0.64 / 1.0e+00 | 0.84 / 1.0e+00 | 0/3 |
| jnk1 | 29 | 7 | 0.12 / 1.0e+00 | 0.12 / 1.0e+00 | 1.19 / 9.3e-01 | 0.45 / 1.0e+00 | 0/3 |
| hsp90_kung | 13 | 3 | 0.19 / 1.0e+00 | 0.18 / 1.0e+00 | 0.83 / 1.0e+00 | 0.89 / 9.7e-01 | 0/3 |
| bace_ciordia_prospective | 11 | 3 | 0.01 / 1.0e+00 | 0.07 / 1.0e+00 | 0.04 / 1.0e+00 | 0.05 / 1.0e+00 | 0/3 |
| chk1 | 15 | 3 | 0.28 / 1.0e+00 | 0.27 / 1.0e+00 | 0.69 / 1.0e+00 | 0.24 / 1.0e+00 | 0/3 |
| cmet | 39 | 16 | 0.17 / 1.0e+00 | 0.36 / 1.0e+00 | 0.42 / 1.0e+00 | 0.30 / 1.0e+00 | 0/3 |
| dlk | 6 | 2 | 0.02 / 1.0e+00 | 0.06 / 1.0e+00 | 0.05 / 1.0e+00 | 0.11 / 1.0e+00 | 0/3 |
| eg5 | 43 | 16 | 0.48 / 1.0e+00 | 0.40 / 1.0e+00 | 0.30 / 1.0e+00 | 0.64 / 1.0e+00 | 0/3 |
| egfr | 7 | 3 | 0.07 / 1.0e+00 | 0.01 / 1.0e+00 | 0.02 / 1.0e+00 | 0.01 / 1.0e+00 | 0/3 |
| factor_xa | 3 | 1 | 0.31 / 1.0e+00 | 0.01 / 1.0e+00 | 0.02 / 1.0e+00 | 0.22 / 1.0e+00 | 0/3 |
| hiv1_protease | 19 | 7 | 0.28 / 1.0e+00 | 0.43 / 1.0e+00 | 0.94 / 1.0e+00 | 0.65 / 1.0e+00 | 0/3 |
| hne | 23 | 7 | 0.13 / 1.0e+00 | 0.28 / 1.0e+00 | 0.07 / 1.0e+00 | 0.29 / 1.0e+00 | 0/3 |
| hsp90_2rings | 7 | 2 | 0.34 / 1.0e+00 | 0.29 / 1.0e+00 | 0.12 / 1.0e+00 | 0.65 / 1.0e+00 | 0/3 |
| hsp90_single_ring | 8 | 2 | 0.16 / 1.0e+00 | 0.02 / 1.0e+00 | 0.24 / 1.0e+00 | 0.02 / 1.0e+00 | 0/3 |
| hsp90_woodhead | 4 | 1 | 0.01 / 1.0e+00 | 0.07 / 1.0e+00 | 0.15 / 1.0e+00 | 0.19 / 1.0e+00 | 0/3 |
| irak4_s2 | 7 | 3 | 0.25 / 1.0e+00 | 0.15 / 1.0e+00 | 0.04 / 1.0e+00 | 0.28 / 1.0e+00 | 0/3 |
| irak4_s3 | 4 | 1 | 0.07 / 1.0e+00 | 0.07 / 1.0e+00 | 0.05 / 1.0e+00 | 0.02 / 1.0e+00 | 0/3 |
| itk | 5 | 2 | 0.68 / 1.0e+00 | 0.08 / 1.0e+00 | 0.10 / 1.0e+00 | 0.03 / 1.0e+00 | 0/3 |
| jak1 | 7 | 2 | 0.02 / 1.0e+00 | 0.32 / 1.0e+00 | 0.10 / 1.0e+00 | 0.12 / 1.0e+00 | 0/3 |
| jak2_set1 | 14 | 5 | 0.35 / 1.0e+00 | 0.83 / 1.0e+00 | 0.15 / 1.0e+00 | 0.69 / 1.0e+00 | 0/3 |
| jak2_set2 | 12 | 5 | 0.06 / 1.0e+00 | 1.02 / 1.0e+00 | 0.49 / 1.0e+00 | 0.36 / 1.0e+00 | 0/3 |
| keranen_p2 | 17 | 6 | 0.29 / 1.0e+00 | 0.19 / 1.0e+00 | 0.04 / 1.0e+00 | 0.11 / 1.0e+00 | 0/3 |
| ptp1b | 36 | 12 | 0.21 / 1.0e+00 | 0.20 / 1.0e+00 | 0.58 / 1.0e+00 | 0.49 / 1.0e+00 | 0/3 |
| syk | 67 | 22 | 0.11 / 1.0e+00 | 0.20 / 1.0e+00 | 0.25 / 1.0e+00 | 0.17 / 1.0e+00 | 0/3 |
| urokinase | 4 | 1 | 0.13 / 1.0e+00 | 0.01 / 1.0e+00 | 0.25 / 1.0e+00 | 0.00 / 1.0e+00 | 0/3 |

## Honest reading
- The paper's six flagged systems are a **replicate-0 realisation** of a test whose
  set-valued output is unstable, not a reproducible list. Reported as a set, the honest
  statement is: `bace` reproduces in every run; the union of 8 systems is what a single run can produce; and the pooled, highest-powered read gives 7.
- The **underlying quantity is stable** (ρ = 0.64–0.72); the instability is created by
  thresholding a continuous, low-power statistic. This is consistent with the median
  reduced χ² ≈ 0.34 already reported in `results_figL.md`, which implies limited
  single-replicate power, and with that document's existing note that the flag set is a
  *lower bound* on the systems carrying systematic error.
- **The detector-level claims are unaffected.** Panel B of Fig L (selectivity is a
  function of σ-calibration) reproduces on all three replicates. What does not survive
  is the membership of any particular six-system list.
- Nothing here was tuned after seeing a result: α, the σ scales, the edge construction
  and the pooling rule are all taken unchanged from `make_figL.py`.

