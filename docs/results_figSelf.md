# Results — Fig Self: two-sided self-calibration of the cycle-closure null

**Figure:** `figs/figSelf_two_sided_selfcalibration.{pdf,png}` · **Reproduce:** `make figSelf`
(`PYTHONPATH=src python figs/make_figSelf.py`). Deterministic (the analysis is a closed-form
recomputation of the same GLS + BH pipeline; the synthetic-null diagnostic is seeded). Data:
`data/openfe_replicates/combined_pymbar4_edge_data.csv` (OpenFE IndustryBenchmarks2024, public,
3 independent replicates per edge).

> **Exploratory, not pre-registered.** The pre-registered robustness check on this null is Fig/Table
> P8 (`figs/make_figP8.py`, `docs/results_figP8.md`): **one-sided**, per-system, verdict **SURVIVES**.
> It is untouched here — not re-run differently, not reinterpreted. This document asks a *different*,
> two-sided question and reports it separately.

## What this asks, and why it is not P8's question

The manuscript states that the detector's **recall** is uncharacterized because there are no
ground-truth error labels. That remains true of recall in the strict sense. But the replicate set
does measure the one thing recall depends on most here: **whether the null the test is run against
is the right width**. `figs/make_figA_replicates.py` measures, per edge, the ratio of the reported
(MBAR = sandwich) se to that edge's own across-replicate SD. This analysis replaces each reported
bar with its own measured one — `se → se / ratio` — **in whichever direction the measurement points**
— and re-runs the paper's own GLS + Benjamini–Hochberg pipeline unchanged at the pre-registered α = 0.05.

P8 applies the same correction **only where it makes a bar looser** (ratio < 1): the conservative
half, which is what a referee's circularity concern needs. The other half — edges whose reported bar
is *wider* than their own measured spread — can only *add* detections, and it is the half that speaks
to recall.

## Design, fixed and stated before the run

| # | choice | what was fixed | why |
|---|---|---|---|
| D1 | **granularity** | **per edge** (pre-specified primary): `ratio_e = rep_e / repl_e` with `rep_e = sqrt(mean_k se_ek²)` (RMS reported se over the 3 replicates) and `repl_e = SD_k(ΔΔG_ek)` (ddof=1) — `make_figA_replicates.load()`'s own per-edge pair, recomputed through `make_figL.edge_val` so the edge construction is shared verbatim | "that edge's own measured calibration" is per-edge. The **per-system aggregate** ratio (Fig A-rep's / P8's pooled RMS formula) was fixed at the same time as the granularity **sensitivity**. Both are reported below; see *Which arm to read* |
| D2 | **undefined / degenerate ratio** | **edge left unchanged** (`ratio := 1`) when the edge has no complete 3-replicate record, or `repl_e ≤ 1e-06` (Fig A-rep's and P8's own threshold), or the arithmetic is non-finite | an unmeasurable ratio carries no evidence about that edge; the neutral action is to keep the reported bar. Neutral, **not** conservative: by construction it can neither add nor remove a flag |
| D3 | **direction convention** | `se → se / ratio`, **both** directions. `ratio > 1` (reported bar wider than the measured spread) ⇒ se shrinks ⇒ **bar tighter** ⇒ larger χ² ⇒ **more** likely flagged. `ratio < 1` (reported bar tighter than the measured spread) ⇒ se grows ⇒ **bar looser** ⇒ **less** likely flagged | P8 keeps only the second branch; two-sided keeps both |
| D4 | **no cap** on the primary | a `[1/3, 3]`-clipped variant is reported as a sensitivity | an n=3 SD ratio has a heavy right tail; the tail is disclosed rather than silently trimmed |
| D5 | **no c4 correction** on the primary | raw n=3 SD ratio, identical to the ratio Fig A-rep panel B and P8 report; the `c4(3)=0.886` bias-corrected variant (`repl → repl / c4`, i.e. wider bars, fewer flags) is a sensitivity | one ratio definition in the paper, not two |

Everything else is Fig L's, unchanged: replicate 0, systems with ≥ 3 edges and dof ≥ 1,
`bar.qc.gls_network` / `bar.qc.chi2_sf` / `bar.qc.benjamini_hochberg` reached through `make_figL`'s own
adapters, and the pre-registered BH-FDR level α = 0.05. The nominal arm is asserted against Fig L's
published values (median reduced χ² 0.34, flag set `bace, brd4, cdk8, faah, hif2a, p38`) before anything
else is reported; on a mismatch the script prints a STOP banner instead of a result.

## What the replicates say about the bars

Measurable per-edge ratios: **1143 edges**; 0 edges unmeasurable and left unchanged (D2).
Median **1.81**, interquartile range 0.91–3.60, full range 0.12–170.0. **72%** of edges have `ratio ≥ 1`,
i.e. a reported bar *wider* than their own measured run-to-run spread, so a two-sided correction
predominantly **tightens** bars. That is the whole mechanism of the result below.

## Headline

| arm | systems flagged | median reduced χ² | lost | gained |
|---|---:|---:|---:|---:|
| nominal (reported se, = Fig L) | **6** of 48 | 0.343 | — | — |
| self-calibrated, **per-system aggregate** | **16** of 48 | 1.572 | 0 | 10 |
| self-calibrated, **per-edge** (pre-specified primary) | **7** of 48 | 1.313 | 2 | 3 |

- Aggregate arm — **lost: 0** (`none`); **gained: 10** (`ciordia_retro, cmet, galectin, irak4_s2, keranen_p2, mcl1, renin, thrombin, tnks2, tyk2`); **kept: 6** (`bace, brd4, cdk8, faah, hif2a, p38`).
- Per-edge arm — **lost: 2** (`faah, hif2a`); **gained: 3** (`irak4_s2, thrombin, tyk2`).

The two granularities do **not** agree, and the difference is not cosmetic: the aggregate arm adds
10 systems and loses none, while the per-edge arm adds 3 and loses 2
(`faah, hif2a`). The next section decides which one may be read, on grounds
that have nothing to do with the flag counts.

## Which arm to read (decided on the null, not on the result)

**The argument.** Correcting an edge's variance by its own measured ratio multiplies that edge's χ²
contribution by `σ²/s²`, where `s²` is the edge's replicate sample variance on `ν = n − 1` degrees of
freedom. `E[σ²/s²] = ν/(ν−2)`, which **diverges for ν ≤ 2**. With n = 3 replicates ν = 2, so the
per-edge self-calibrated χ² has *no finite null expectation*: it is anti-conservative by construction,
not by accident. The held-out ratio, built on n = 2 (`ν = 1`), is worse still. The aggregate ratio pools
the denominator over a whole system (`ν ≈ 2E`), leaving only the finite-sample factor `2E/(2E−2)` — 11% at
E = 10 edges, 3% at E = 40 — so it is *nearly*, but not exactly, unbiased.

**The measurement.** Each arm was then run on a *perfectly-calibrated synthetic null*: the same graphs
and the same reported bars, with ΔΔG redrawn from those bars around exact node potentials and three
replicates per edge, so every flag is a false positive by construction (400 draws; Monte-Carlo se ≤ 0.025).
`null P(≥1 flag)` is the probability that the arm flags *anything* when nothing is wrong; under BH-FDR
at α = 0.05 that should be ≈ 0.05. `null mean #flagged` is how many of the 48 systems it falsely flags per draw.

| arm | systems flagged | median χ²ᵥ | null P(≥1 flag) | null mean #flagged | null median χ²ᵥ | flag set |
|---|---:|---:|---:|---:|---:|---|
| nominal (reported se) = Fig L | 6 | 0.343 | 0.040 | 0.04 | 0.860 | `bace, brd4, cdk8, faah, hif2a, p38` |
| pre-specified primary: per-edge, two-sided | 7 | 1.313 | 0.588 | 1.12 | 1.053 | `bace, brd4, cdk8, irak4_s2, p38, thrombin, tyk2` |
| granularity sensitivity: per-system aggregate, two-sided | 16 | 1.572 | 0.240 | 0.33 | 0.916 | `bace, brd4, cdk8, ciordia_retro, cmet, faah, galectin, hif2a, irak4_s2, keranen_p2, mcl1, p38, renin, thrombin, tnks2, tyk2` |
| per-edge, ratio capped to [1/3, 3] | 6 | 0.868 | 0.380 | 0.54 | 1.019 | `bace, brd4, cdk8, p38, thrombin, tyk2` |
| per-edge, c4(3)-corrected SD | 6 | 1.031 | 0.230 | 0.28 | 0.827 | `bace, brd4, cdk8, p38, thrombin, tyk2` |
| per-edge, one-sided (P8's direction) | 2 | 0.308 | 0.003 | 0.00 | 0.635 | `bace, cdk8` |
| per-edge, held-out ratio (replicates 1+2 only) | 21 | 1.923 | 1.000 | 18.50 | 1.807 | `bace, brd4, cdk2, cdk8, eg5, egfr, faah, factor_xa, hif2a, irak4_s2, itk, keranen_p2, liga, mcl1, p38, ptp1b, renin, shp2, syk, thrombin, tyk2` |

The nominal (published) test comes out at P(≥1) = 0.040, mean 0.04 false flags — it controls its own
false-positive rate, as the paper claims. **Every self-calibrated arm inflates it**, in the order the
divergence argument predicts:

- per-system aggregate: P(≥1) = 0.240, mean **0.33** false flags per draw;
- per-edge (n = 3, ν = 2): P(≥1) = 0.588, mean 1.12;
- held-out (n = 2, ν = 1): P(≥1) = 1.000, mean **18.5** — it flags 39% of all systems when nothing is wrong.

So no self-calibrated arm is FDR-controlled at the nominal level, and **none of them may be reported as
a calibrated flag set**. What the harness does license is a *magnitude* comparison for the arm whose
inflation is smallest: the aggregate arm falsely flags 0.33 systems per draw when nothing is wrong,
against **10 systems gained** on the real data — a factor of 30. The gains are not the correction's own
false-positive inflation. The per-edge arm cannot support even that statement (it falsely flags 1.12
per draw against 3 gained, and it *loses* 2), and the held-out arm is uninformative rather than a
stronger version of the claim. Everything below is the **aggregate** arm.

## The result (aggregate arm)

**6 of 48 → 16 of 48** flagged; median reduced χ² **0.343 → 1.572**; **0 lost, 10 gained**.

The reading is **not** that the flag set should be 16 systems — that arm is not FDR-controlled (above). It is that the
observed miscalibration of the null **suppresses** detections rather than manufacturing them: correcting
each bar to the width its own replicates measure *adds* systems and removes none, by a margin (10 gained
vs 0.33 expected from the correction's own false-positive inflation) that the inflation cannot explain. The
published flag count is therefore bounded **below**, which is the direction the manuscript already claims
on other grounds (median reduced χ² < 1 implies limited single-replicate power) — this makes that
statement quantitative instead of rhetorical.

The zero-loss half is the more robust of the two, because false-positive inflation is the failure mode that
*adds* flags — it cannot manufacture the absence of losses. And losses do happen: the per-edge arm loses 2
(`faah, hif2a`). That all 6 published systems survive being re-tested against their own
measured bars is therefore a result, not an automatic consequence of the arithmetic.

**One of the gains is a threshold effect, not a bar effect.** Benjamini–Hochberg's cutoff relaxes as other systems'
p-values fall, so a system can cross without its own χ² moving much. The gains whose reduced χ² changes by
less than 15%: `mcl1` (1.82 → 1.87, ratio 1.01). Those are carried by the rest of the list rather than by
their own evidence, and should not be read individually.

`renin` is the informative gain: the manuscript already names it as a **systematic-but-unflagged**
case (systematic by the pooled-replicate ratio, reduced χ² 2.5→4.7, yet short of FDR on a single
replicate). Under its own measured bars it clears FDR (χ²ᵥ 1.29 → 2.53, q 8.1e-01 → 6.6e-03).
A case the manuscript already documents as systematic-but-unflagged is recovered by correcting the null
it was tested against — the closest thing this data affords to a recall measurement, and the one gain
with independent corroboration in the paper.

### Stated dependence (a caveat, not a knob)

The ratio is measured on all three replicates, including replicate 0, whose ΔΔG values the closure
test itself reads, so the correction factor and the tested residuals are not fully independent: a system
that happened to sample a small spread on this replicate gets both tighter bars and, on average, smaller
residuals. The held-out variant was built as the control for exactly this — but the synthetic null shows
it is invalid at n = 2 (`ν = 1`), so **this dependence is not resolved by the data available here** and is
left as a stated limitation. Its direction is not obvious a priori; the aggregate ratio pools ~10–60 edges
per system, which dilutes but does not eliminate it.

## Full per-system before-and-after

Every admitted system, ordered by its self-calibrated *q* (aggregate arm). `agg. ratio` is that system's
pooled RMS reported-se / replicate-SD. The last two columns give the per-edge arm for completeness; per
the section above they are **not** a detection count.

| system | edges | cycles | χ²ᵥ nominal | q nominal | agg. ratio | χ²ᵥ agg. | q agg. | transition (agg.) | χ²ᵥ per-edge | q per-edge |
|---|---:|---:|---:|---:|---:|---:|---:|---|---:|---:|
| cdk8 | 63 | 29 | 2.678 | 6.1e-05 | 1.40 | 5.286 | 3.6e-17 | **kept** | 9.409 | 8.5e-40 |
| p38 | 60 | 22 | 2.396 | 2.0e-03 | 1.43 | 4.880 | 7.8e-12 | **kept** | 3.346 | 2.3e-06 |
| bace | 49 | 14 | 5.572 | 3.2e-09 | 0.96 | 5.134 | 1.4e-08 | **kept** | 10.922 | 4.5e-24 |
| thrombin | 54 | 19 | 1.301 | 7.4e-01 | 1.76 | 4.044 | 7.8e-08 | **gained** | 3.581 | 2.3e-06 |
| galectin | 36 | 11 | 0.370 | 1.0e+00 | 3.67 | 4.997 | 7.5e-07 | **gained** | 0.980 | 6.5e-01 |
| tyk2 | 29 | 11 | 0.788 | 1.0e+00 | 2.38 | 4.461 | 7.4e-06 | **gained** | 4.392 | 1.2e-05 |
| ciordia_retro | 45 | 14 | 0.895 | 1.0e+00 | 1.94 | 3.361 | 1.3e-04 | **gained** | 1.685 | 2.0e-01 |
| hif2a | 59 | 19 | 2.534 | 2.0e-03 | 0.96 | 2.358 | 4.4e-03 | **kept** | 1.807 | 9.0e-02 |
| renin | 42 | 14 | 1.293 | 8.1e-01 | 1.40 | 2.535 | 6.6e-03 | **gained** | 1.582 | 2.3e-01 |
| brd4 | 8 | 1 | 15.397 | 1.4e-03 | 0.76 | 8.968 | 1.1e-02 | **kept** | 9.711 | 1.5e-02 |
| faah | 31 | 8 | 3.705 | 2.0e-03 | 0.90 | 2.993 | 1.1e-02 | **kept** | 2.414 | 8.0e-02 |
| irak4_s2 | 7 | 3 | 0.248 | 1.0e+00 | 4.37 | 4.725 | 1.1e-02 | **gained** | 4.759 | 1.7e-02 |
| tnks2 | 35 | 9 | 0.628 | 1.0e+00 | 2.06 | 2.654 | 1.7e-02 | **gained** | 1.855 | 2.0e-01 |
| mcl1 | 76 | 24 | 1.825 | 5.5e-02 | 1.01 | 1.873 | 2.0e-02 | **gained** | 1.043 | 6.5e-01 |
| keranen_p2 | 17 | 6 | 0.288 | 1.0e+00 | 3.21 | 2.970 | 2.1e-02 | **gained** | 1.658 | 3.0e-01 |
| cmet | 39 | 16 | 0.174 | 1.0e+00 | 3.43 | 2.049 | 2.4e-02 | **gained** | 1.525 | 2.3e-01 |
| hne | 23 | 7 | 0.132 | 1.0e+00 | 4.18 | 2.311 | 6.6e-02 | — | 1.496 | 3.4e-01 |
| hsp90_2rings | 7 | 2 | 0.336 | 1.0e+00 | 3.24 | 3.518 | 7.9e-02 | — | 1.812 | 3.4e-01 |
| liga | 13 | 3 | 2.103 | 4.7e-01 | 1.13 | 2.691 | 1.1e-01 | — | 2.372 | 2.2e-01 |
| itk | 5 | 2 | 0.680 | 1.0e+00 | 2.02 | 2.785 | 1.5e-01 | — | 1.921 | 3.3e-01 |
| taf12 | 8 | 1 | 3.208 | 3.9e-01 | 1.00 | 3.213 | 1.7e-01 | — | 2.935 | 2.3e-01 |
| cdk2 | 27 | 10 | 0.477 | 1.0e+00 | 1.88 | 1.689 | 1.7e-01 | — | 2.074 | 1.1e-01 |
| dlk | 6 | 2 | 0.015 | 1.0e+00 | 12.47 | 2.359 | 2.0e-01 | — | 1.127 | 5.8e-01 |
| eg5 | 43 | 16 | 0.480 | 1.0e+00 | 1.74 | 1.452 | 2.2e-01 | — | 1.490 | 2.4e-01 |
| hiv1_protease | 19 | 7 | 0.278 | 1.0e+00 | 2.39 | 1.584 | 2.6e-01 | — | 0.499 | 9.1e-01 |
| egfr | 7 | 3 | 0.069 | 1.0e+00 | 4.74 | 1.560 | 3.6e-01 | — | 2.507 | 2.0e-01 |
| jak2_set1 | 14 | 5 | 0.351 | 1.0e+00 | 2.03 | 1.450 | 3.6e-01 | — | 0.663 | 8.0e-01 |
| factor_xa | 3 | 1 | 0.306 | 1.0e+00 | 2.14 | 1.400 | 4.1e-01 | — | 1.380 | 4.8e-01 |
| shp2 | 37 | 12 | 0.742 | 1.0e+00 | 1.27 | 1.191 | 4.7e-01 | — | 0.893 | 7.4e-01 |
| scyt_dehyd | 7 | 1 | 0.543 | 1.0e+00 | 1.35 | 0.982 | 5.1e-01 | — | 0.924 | 5.8e-01 |
| ptp1b | 36 | 12 | 0.206 | 1.0e+00 | 2.27 | 1.059 | 6.1e-01 | — | 1.933 | 1.1e-01 |
| hsp90_single_ring | 8 | 2 | 0.157 | 1.0e+00 | 2.27 | 0.808 | 6.6e-01 | — | 0.854 | 6.5e-01 |
| urokinase | 4 | 1 | 0.127 | 1.0e+00 | 2.10 | 0.560 | 6.6e-01 | — | 0.575 | 6.5e-01 |
| ephx2 | 4 | 1 | 0.766 | 1.0e+00 | 0.77 | 0.452 | 7.1e-01 | — | 0.410 | 7.2e-01 |
| hsp90_kung | 13 | 3 | 0.191 | 1.0e+00 | 1.96 | 0.738 | 7.3e-01 | — | 1.247 | 5.4e-01 |
| irak4_s3 | 4 | 1 | 0.069 | 1.0e+00 | 1.80 | 0.225 | 8.4e-01 | — | 0.875 | 5.8e-01 |
| mup1 | 7 | 2 | 0.554 | 1.0e+00 | 0.89 | 0.438 | 8.4e-01 | — | 0.488 | 7.8e-01 |
| btk | 8 | 3 | 0.190 | 1.0e+00 | 1.57 | 0.470 | 8.9e-01 | — | 0.614 | 7.8e-01 |
| chk1 | 15 | 3 | 0.280 | 1.0e+00 | 1.19 | 0.399 | 9.1e-01 | — | 0.919 | 6.5e-01 |
| jnk1 | 29 | 7 | 0.121 | 1.0e+00 | 2.23 | 0.600 | 9.1e-01 | — | 0.389 | 9.5e-01 |
| bace1 | 3 | 1 | 0.052 | 1.0e+00 | 1.17 | 0.072 | 9.2e-01 | — | 0.072 | 8.9e-01 |
| hsp90_woodhead | 4 | 1 | 0.010 | 1.0e+00 | 1.39 | 0.020 | 9.9e-01 | — | 0.096 | 8.9e-01 |
| t4_lysozyme | 14 | 3 | 0.247 | 1.0e+00 | 0.99 | 0.241 | 9.9e-01 | — | 0.231 | 9.3e-01 |
| bace_ciordia_prospective | 11 | 3 | 0.012 | 1.0e+00 | 3.06 | 0.113 | 1.0e+00 | — | 0.104 | 9.7e-01 |
| bace_p3_arg368_in | 28 | 8 | 1.876 | 3.5e-01 | 0.41 | 0.321 | 1.0e+00 | — | 0.577 | 8.9e-01 |
| jak1 | 7 | 2 | 0.017 | 1.0e+00 | 1.97 | 0.065 | 1.0e+00 | — | 0.032 | 9.7e-01 |
| jak2_set2 | 12 | 5 | 0.060 | 1.0e+00 | 1.12 | 0.075 | 1.0e+00 | — | 0.491 | 8.9e-01 |
| syk | 67 | 22 | 0.109 | 1.0e+00 | 1.97 | 0.423 | 1.0e+00 | — | 1.169 | 5.1e-01 |

## Honest reading

- This is **exploratory**. It is reported next to — not instead of — the pre-registered one-sided check
  (P8, verdict SURVIVES), which is unchanged.
- The claim it supports is narrow and one-directional: **the observed miscalibration suppresses
  detections**. It does not license reporting 16 flagged systems as the paper's flag set — no self-calibrated arm here is
  FDR-controlled — and it does not characterize recall in the strict sense: there are still no
  ground-truth error labels, only a measured proxy for the width of the null.
- **The pre-specified primary granularity was the wrong choice, and is reported as such.** Per-edge
  self-calibration on n = 3 replicates has no finite null expectation and inflates its own false-positive
  rate to P(≥1) = 0.59 against a nominal 0.05; the aggregate arm, fixed before the run as the
  granularity sensitivity, is the least-inflated one. The design was not changed after seeing the flag
  counts — the arm to read was chosen on a synthetic null that never sees a real ΔΔG value.
- **The magnitude is load-bearing, the flag list is not.** Read "ten more systems than the correction's
  own noise can produce (0.33 per draw)", not "these ten systems are systematic". Individual gains sitting
  near q ≈ α are exactly the ones the residual inflation can move.
- The correction is **measured, not assumed**: its direction on each system comes from that system's own
  replicates. Where nothing is measurable, nothing is changed.
- Nothing was tuned after seeing a flag count. D1–D5, the population, α and the full sensitivity list were
  fixed and written into the script's docstring before the first run; the synthetic-null diagnostic was
  added afterwards on the strength of the `ν/(ν−2)` argument, and it changes no design choice.


## Realized level of the shipped test by calibration scale

Same seed, same 400 draws and the same graphs as the arm table above.
The truth is drawn from the reported bars; the test then divides by a uniform
scale times those bars, so scale 1 is the shipped test.

| calibration scale (true/reported) | P(any false flag) |
|---|---|
| 1.04 | 0.100 |
| 1.00 | 0.045 |
| 0.92 | 0.003 |
| 0.79 | 0.000 |

Scale 1.00 estimates the same quantity as the uncorrected arm of the table
above, by a different path, and the two agree within Monte-Carlo error at this
draw count.
