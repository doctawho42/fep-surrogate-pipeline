# Results — Fig Inf: inference for six claims that were stated without it

**Figure:** `figs/figInf_inference.{pdf,png}` · **Reproduce:** `make figInf`
(`PYTHONPATH=src python figs/make_figInf.py`). Deterministic (every resampling seeded,
`SEED = 20260810`). Data: `data/openfe_replicates/combined_pymbar4_edge_data.csv` (OpenFE
IndustryBenchmarks2024, public, 3 independent replicates per edge). All fits are the
repo's own (`bar.qc.gls_network`, `bar.leverage.curl_leverage`, `bar.detectors`); the
edge construction and the detectors are imported from `make_figStab.py`,
`make_figCut.py`, `make_figOOS.py`, `make_figL_validation.py` and
`make_figA_replicates.py` rather than reimplemented, and C1 and C2 each assert that they
reproduce the shipped number before adding anything to it.

Every design below was fixed before its number was computed and is stated in the script
docstring. Nothing was tuned after seeing a result. Several of these come back against
the manuscript; they are reported as such.

## Summary

| check | verdict | one line |
|---|---|---|
| C1 | **supports the manuscript, with the caveat now stated** | all three r are positive with cluster-bootstrap CIs excluding zero and permutation p ≤ 0.0005; the flagged-vs-unflagged contrast is positive on all three pairs but its system-label permutation p is 0.098–0.170. |
| C2 | **AGAINST the manuscript** | calibrated vs the 1.0 kcal/mol cutoff: ΔAUC = +0.290 [+0.094, +0.490], which excludes zero, so the two-way tie claim is wrong as written. |
| C3 | **AGAINST the manuscript as written** | under the matched (leverage-weighted = dof-weighted) rule the observed aggregate is 1.32 [0.81, 1.80] against the predicted 0.85 [0.63, 1.08]: the intervals overlap but the predicted interval does not cover the observed point, and the observed value swings 0.34–1.32 across rules. |
| C4 | **a real limitation, now quantified** | 9 of 48 systems have a single cycle (quartiles 2/5/12); 7 of the 8 ever-flagged systems and the only always-flagged one are high-dof, and the per-system χ²ᵥ swing is 14.2× at dof = 1 versus 3.9× at dof ≥ 2. |
| C5 | **supports the manuscript** | the shared-denominator null realizes ρ = +0.059 [-0.204, +0.364], max +0.489; the reported 0.648 exceeds every one of the 200 null draws (p = 0.0050), so it is not a shared-denominator artefact — but the null is not centred on zero and the comparison should be made against it, not against zero. |
| C6 | **against the manuscript's reconciliation** | on Wade et al.'s own functional the analytic se falls below the replicate spread on 27.9% [22.4, 33.2]% of edges, *below* the 36.8% a perfectly calibrated bar produces at n = 3, so the like-for-like comparison makes the disagreement with their under-estimation finding sharper, not milder. |

## C1 — inference for the cross-replicate residual correlation

The 1143 edges are nested in 48 networks and each residual is a
projection through that network's residual maker, so neither an interval nor a test may
treat edges as independent. The interval is a **cluster bootstrap that resamples the
48 systems** (2000 draws); the test is a **within-system permutation**
(1999 draws) that permutes the edge order inside each system, destroying the
cross-replicate pairing while preserving system membership, system sizes and the marginal
residual distribution. Under the sampling null the three replicate fits are independent,
so `E[z_e^(i) z_f^(j)] = 0` for every pair of edges including `e = f`; the permutation
therefore has the right null, and its realized mean (below) confirms it is centred on
zero rather than on the projector-induced level that a naive within-replicate shuffle
would produce.

| pair | pooled r | cluster-bootstrap 95% CI | within-system permutation p | permutation null mean ± SD |
|---|---:|---|---:|---|
| (0,1) | +0.302 | [+0.073, +0.474] | 0.0005 | +0.0041 ± 0.0387 |
| (0,2) | +0.420 | [+0.171, +0.613] | 0.0005 | +0.0078 ± 0.0414 |
| (1,2) | +0.379 | [+0.238, +0.475] | 0.0005 | +0.0043 ± 0.0400 |

**The three rotations are dependent.** They are three views of the same three runs, and
every pair of rotations shares a replicate — (0,1) and (0,2) share replicate 0, and so
on. They are not three replications, the three p-values are not independent evidence,
and no multiplicity correction across them would be meaningful. One permutation `pi` per
draw is applied to all three pairs, so the null inherits the same dependence.

### Where the correlation lives
270 of the 1143 edges belong to the 6
flagged systems. The contrast gets a **stratified** cluster bootstrap (flagged systems
resampled among themselves, unflagged among themselves, so neither stratum can empty) and
two permutation nulls: a **system-label** permutation that reassigns which
6 of the 48 systems carry the flag, holding all data fixed —
the null for *concentration* — and the same within-system edge permutation as above,
whose null is *no correlation anywhere* and is therefore the easier bar.

| pair | r (flagged) | r (unflagged) | difference | stratified cluster CI | system-label permutation p | within-system permutation p |
|---|---:|---:|---:|---|---:|---:|
| (0,1) | +0.434 | +0.142 | +0.292 | [-0.079, +0.656] | 0.1550 | 0.0005 |
| (0,2) | +0.585 | +0.199 | +0.386 | [+0.013, +0.659] | 0.0975 | 0.0005 |
| (1,2) | +0.454 | +0.302 | +0.152 | [-0.086, +0.369] | 0.1695 | 0.0210 |

The contrast is positive on all three pairs, so the direction the manuscript asserts is
the direction in the data. But the test that actually asks the manuscript's question —
the system-label permutation — gives p = 0.155, 0.098, 0.170. None of the three reaches α = 0.05, and the stratified cluster interval covers zero on 2 of the 3 pairs.
The within-system edge permutation is much smaller (p = 0.0005, 0.0005, 0.0210), but its null is *no correlation anywhere*, which is not the claim being made; a contrast can beat that null purely
because the correlation exists at all.

**This is a qualification the manuscript needs.** With 6 flagged
systems out of 48, the label permutation has few distinguishable arrangements
of the clusters that carry the signal, so the test is underpowered and the honest reading
is *not established* rather than *refuted*. The existence of the cross-replicate
correlation (first table) is solid at cluster-honest inference; its **concentration** in
the flagged systems is directional evidence, not a demonstrated effect.

## C2 — the missing half of the two-way null

The manuscript reports one interval, `+0.140 [-0.008, +0.316]`, which is
`AUC(A) − max(AUC(B), AUC(C))` and hence the comparison against the **fixed-se** foil.
The comparison against the fixed 1.0 kcal/mol hysteresis cutoff had no interval anywhere.
All three contrasts below come from the SAME 2000 system resamples with the same
`seed = 0` as `make_figCut.py`, so the published row reproduces exactly (asserted in
code):

| contrast | AUC (calibrated) | AUC (foil) | paired ΔAUC | 95% CI | verdict |
|---|---:|---:|---:|---|---|
| vs fixed 1.0 kcal/mol cutoff | 0.865 | 0.575 | +0.290 | [+0.094, +0.490] | **WIN** |
| vs fixed-se χ² test | 0.865 | 0.725 | +0.140 | [-0.006, +0.322] | **TIE** |
| vs max(both) — the published row | 0.865 | — | +0.140 | [-0.008, +0.316] | **TIE** |

**This is against the manuscript.** The calibrated null beats the 1.0 kcal/mol cutoff by +0.290 AUC with a 95% CI
of [+0.094, +0.490], which excludes
zero. The sentence 'does not measurably out-discriminate either foil' is true of the
fixed-se foil and false of the fixed-cutoff foil. The pre-registered WIN rule in
`bar.detectors.paired_auc_bootstrap` is conjunctive (A must beat BOTH), so the overall
TIE verdict stands unchanged — but the *reason* it is a tie is entirely the fixed-se
comparison, and the text must say that instead of implying both comparisons are null.
For context, the cutoff rule flags 26 of the 48 systems and
the calibrated rule flags 6.

## C3 — the observed aggregate the prediction is said to match

The predicted aggregate is `Σ_s Σ_e h_e c_e^−2 / Σ_s dof_s`. Because `Σ_e h_e = dof_s`
(Theorem D1) and `E[X²_s] = Σ_e h_e c_e^−2`, the observed statistic under the **same**
rule is `Σ_s X²_s / Σ_s dof_s` — that is, the leverage-weighted rule and the dof-weighted
global value are the same number, so only three of the four requested aggregates are
distinct. All are reported:

| aggregation rule | value | 95% CI (system cluster bootstrap) |
|---|---:|---|
| **predicted** E[χ²ᵥ], leverage-weighted (published) | 0.85 | [0.63, 1.08] |
| **observed, matched rule** (= dof-weighted global Σχ²/Σdof) | 1.32 | [0.81, 1.80] |
| observed, median over systems | 0.34 | — |
| observed, unweighted mean over systems | 1.14 | — |

**The intervals overlap**
(0.81 ≤ 1.08), but the predicted interval
does not cover the observed point:
1.32 against a predicted upper bound of 1.08, i.e. the
observed aggregate is 1.55× the predicted one. 'Matches'
survives only in the weak sense that two wide intervals share ground.

**The rule matters more than the agreement does.** The same 48 systems give
0.34 under a median over systems, 1.14 under an
unweighted mean over systems and 1.32 under the matched rule — a factor
of 3.8 between the extremes — because the closure χ² is
dominated by a few high-dof, high-residual networks while the median system closes far
better than its bars predict. Quoting 'the aggregate level also matches' without naming
the rule or the observed number is therefore not a checkable claim, and it is the reader
who has to guess which of these three the sentence means. Under the only rule that is
actually matched to the prediction, observed sits above predicted.

That direction is the expected one and is already the paper's own explanation — a
reproducible systematic error inflates a closure residual without inflating the spread
between repeats, so observed *should* exceed predicted wherever the QC fires. But that is
an argument for an interpretable discrepancy, not for a match, and the sentence should
state the observed number, the rule, and the direction rather than assert agreement.

## C4 — the degrees-of-freedom distribution, and stability stratified by it

Across the 48 admitted systems (replicate 0): min 1, quartiles
2 / 5 / 12, max 29, mean
7.8. Counts at the low end: 9 systems with one cycle,
6 with two, 8 with three.

| dof | systems |
|---:|---:|
| 1 | 9 |
| 2 | 6 |
| 3 | 8 |
| 5 | 2 |
| 6 | 1 |
| 7 | 3 |
| 8 | 2 |
| 9 | 1 |
| 10 | 1 |
| 11 | 2 |
| 12 | 2 |
| 14 | 3 |
| 16 | 2 |
| 19 | 2 |
| 22 | 2 |
| 24 | 1 |
| 29 | 1 |

### Flag stability stratified by dof (primary split: median of dof)

| stratum | systems | flagged per replicate | ever / always | pairwise Jaccard | mean (random ref) | median χ²ᵥ swing |
|---|---:|---|---|---|---|---:|
| low: dof <= 5 (median split) | 25 | 1/1/0 | 1 / 0 | 1.00 / 0.00 / 0.00 | 0.333 (0.013) | 6.45× |
| high: dof > 5 | 23 | 5/5/3 | 7 / 1 | 0.43 / 0.33 / 0.33 | 0.365 (0.110) | 2.39× |

### Secondary split: single-cycle systems versus the rest

| stratum | systems | flagged per replicate | ever / always | pairwise Jaccard | mean (random ref) | median χ²ᵥ swing |
|---|---:|---|---|---|---|---:|
| dof == 1 (single cycle) | 9 | 1/1/0 | 1 / 0 | 1.00 / 0.00 / 0.00 | 0.333 (0.037) | 14.16× |
| dof >= 2 | 39 | 5/5/3 | 7 / 1 | 0.43 / 0.33 / 0.33 | 0.365 (0.063) | 3.94× |

**Stability is carried by the high-cycle systems, and this is said plainly.** The
evidence is not the Jaccard column — the low-dof stratum's
0.333 is computed over
1 ever-flagged system with per-replicate counts
1/1/0, so it is a degenerate average of
1.00 and two 0.00s and should not be compared with anything. The evidence is:

- **Where the flags are.** 7 of the 8 ever-flagged systems sit
  in the high-dof stratum, including the only system flagged in *all three* replicates.
  The low-dof half of the benchmark (25 systems) contributes
  1 ever-flagged system and nothing that reproduces.
- **Why.** The per-system χ²ᵥ swing across the three repeats has median 14.16× at dof = 1 against 3.94× at dof ≥ 2 (and 6.45× vs 2.39× under
  the median split). A single-cycle network's reduced χ² *is* one number, so it moves by
  an order of magnitude run to run and its BH-adjusted q crosses α about as often as not.
- **How much of the benchmark this is.** 9 of 48 systems
  have one cycle and 23 have three or
  fewer, so this is not an edge case: a third of the benchmark is structurally incapable
  of delivering a stable per-system verdict, whatever the detector.

A Jaccard of 1.00 in a stratum with no flagged systems is the empty-set convention, not
evidence of stability; the 'flagged per replicate' column shows where that applies. Any
set-valued claim in the manuscript should be scoped to the well-determined networks, and
the single-cycle systems should be reported as flag-eligible but not flag-stable.

## C5 — is the predicted-versus-observed check circular?

**Which se enters where.** In the headline rotation the predicted quantity is
`Σ_e h_e (s_e² / se_{e,0}²) / dof`, where `s_e²` is the sample variance of the ΔΔG values
of replicates **1 and 2** and `se_{e,0}` is the reported standard error of replicate
**0**. The observed quantity is `Σ_e r_{e,0}² / se_{e,0}² / dof`. So replicate 0's
reported se is the denominator on both sides, and the curl-leverage weights `h_e` are
computed from replicate 0's `V_e` as well. The numerators are disjoint (replicates 1–2
versus replicate 0), but the denominators and the weights are shared, which is exactly
the referee's concern.

**The null.** `make_figOOS.resample_null_world` — the repo's own no-systematic-error
world — redraws each replicate's values as `y ~ N(0, se_k)` on the same graphs with the
same reported errors, so `c_e ≡ 1` identically while `h_e` and `se_e` are untouched. The
shared denominator is therefore preserved exactly, and any correlation the harness
realizes is the artefact. Over 200 realizations:

- realized Spearman ρ: mean +0.059, median +0.047, SD
  0.147, 95% range [-0.204, +0.364], max
  +0.489;
- the reported ρ = 0.648 against that null: p = 0.0050;
- positive control (a genuine per-system bar miscalibration spanning 0.5×–2.0×, same
  harness, same shared denominator): median ρ = +0.601
  [+0.363, +0.719], n = 50 — so the
  harness
  can see a real signal through the shared denominator, and its silence under the null is
  informative rather than a lack of power.

**Verdict: not circular.** The shared `V_e` does not manufacture a rank correlation,
because under the null both sides are centred on their own expectations regardless of the
magnitude of `se_e`. The 0.648 should nevertheless be reported against this null
rather than against zero, since 'zero' was never the operative alternative; the null's
95% range is [-0.204, +0.364], which is where an artefact would
have shown up.

## C6 — a like-for-like comparison with Wade et al.

Wade et al. (JCTC 2022) compare the analytic MBAR uncertainty against the spread over
independent replicas and report that MBAR *under*-estimates. The manuscript's
reconciliation invokes a curl-leverage-weighted mean of `c_e^−2`, a functional with no
counterpart in their measurement. Their functional, computed here on these edges: the
per-edge analytic se (RMS over the three replicates, the repo's own convention in
`make_figA_replicates.load`) against the across-replicate SD, and the **fraction of edges
on which the analytic se falls below the replicate spread**.

**Reference level, stated before running.** With three replicates and perfectly
calibrated bars, `s² ~ σ² χ²₂ / 2`, so `P(s > σ) = e^{−1} = 0.368`: about a third of edges land below by construction, and only an
excess over that is evidence of under-estimation. The small-sample-bias-corrected variant
(`σ̂ = s / c₄`, `c₄ = 0.886`) has reference `e^{−c₄²} = 0.456`.

The edge set is `make_figA_replicates.load`'s (1145 edges over
49 systems: every row complete in all three replicates with a non-degenerate
spread), which is the set that produced the manuscript's 1.41×. It is admitted by a
different rule from C1's 1143 network-fitted edges, so the two counts are not expected to
agree.

- **Pooled: the analytic se is below the replicate spread on 0.279**
  [0.224, 0.332] of the
  1145 edges (system-cluster bootstrap), against the
  0.368 a perfectly calibrated bar produces.
- c₄-corrected: 0.322 [0.264, 0.379] against
  0.456.
- For orientation, the same edges give an RMS-pooled ratio of 1.41× and
  a **median per-edge** ratio of 1.81×.
- Per target (34 targets with ≥ 8 edges):
  9 exceed the 0.368 reference and
  6 have an RMS ratio below 1.

| target | edges | fraction of edges with analytic se < replicate spread | RMS ratio | median per-edge ratio |
|---|---:|---:|---:|---:|
| brd4 | 8 | **0.625** | 0.76 | 0.55 |
| bace_p3_arg368_in | 28 | **0.536** | 0.41 | 0.88 |
| renin | 42 | **0.500** | 1.40 | 1.19 |
| t4_lysozyme | 14 | **0.500** | 0.99 | 1.00 |
| mcl1 | 76 | **0.487** | 1.01 | 1.09 |
| hif2a | 59 | **0.458** | 0.96 | 1.07 |
| faah | 31 | **0.452** | 0.90 | 1.10 |
| bace | 49 | **0.429** | 0.96 | 1.11 |
| btk | 8 | **0.375** | 1.57 | 1.76 |
| chk1 | 15 | 0.333 | 1.19 | 1.30 |
| thrombin | 54 | 0.333 | 1.76 | 1.55 |
| shp2 | 37 | 0.324 | 1.27 | 1.36 |
| hiv1_protease | 19 | 0.316 | 2.39 | 1.88 |
| liga | 13 | 0.308 | 1.13 | 1.55 |
| jak2_set1 | 14 | 0.286 | 2.03 | 1.62 |
| cdk2 | 27 | 0.259 | 1.88 | 1.73 |
| tnks2 | 35 | 0.257 | 2.06 | 2.41 |
| cdk8 | 63 | 0.254 | 1.40 | 1.36 |
| galectin | 36 | 0.250 | 3.67 | 2.46 |
| hsp90_single_ring | 8 | 0.250 | 2.27 | 2.45 |
| taf12 | 8 | 0.250 | 1.00 | 2.03 |
| hsp90_kung | 13 | 0.231 | 1.96 | 2.02 |
| ciordia_retro | 45 | 0.222 | 1.94 | 2.13 |
| p38 | 60 | 0.217 | 1.43 | 1.96 |
| eg5 | 43 | 0.209 | 1.74 | 1.72 |
| jak2_set2 | 12 | 0.167 | 1.12 | 1.94 |
| cmet | 39 | 0.154 | 3.43 | 2.98 |
| syk | 67 | 0.104 | 1.97 | 2.48 |
| jnk1 | 29 | 0.103 | 2.23 | 1.94 |
| bace_ciordia_prospective | 11 | 0.091 | 3.06 | 3.37 |
| tyk2 | 29 | 0.069 | 2.38 | 2.64 |
| keranen_p2 | 17 | 0.059 | 3.21 | 3.21 |
| ptp1b | 36 | 0.028 | 2.27 | 3.29 |
| hne | 23 | 0.000 | 4.18 | 4.00 |

**Reading — this goes against the manuscript's reconciliation, not with it.** On the
functional Wade et al. actually report, these edges are *more* conservative than a
perfectly calibrated bar would be: 27.9% of edges fall below the
replicate spread against the 36.8% that perfect calibration
produces
at n = 3, and the bootstrap interval [22.4, 33.2]% excludes the
reference. The median per-edge ratio is 1.81×, *above* the RMS-pooled 1.41×, so the typical
edge is more conservative than the aggregate ratio suggests rather than less. The
like-for-like comparison therefore makes the disagreement with their under-estimation
finding **sharper**, and the manuscript's sentence that the measurement is 'much closer
to their picture than 1.41 suggests' is supported only by the leverage-weighted
functional, which has no counterpart in what they measured. A reconciliation that only
works in a functional the other paper never computes is not a reconciliation.

The one caveat that does cut the other way is the manuscript's own, and it is not settled
here: these repeats share starting coordinates and differ only in stochastic sampling,
whereas Wade et al. resample the full ensemble, so this denominator is a lower bound on
true run-to-run variability and every fraction above is an *under*-estimate of the
fraction Wade et al. would measure. That argument alone has to carry the whole
reconciliation; the numbers do not help it.

The heterogeneity is real and should be reported with the aggregate: 9 of 34 targets (≥ 8
edges) exceed the reference fraction and 6 have an RMS ratio
below 1, so 'conservative in aggregate, heterogeneous per target' remains the accurate
summary — it is only the *direction of the reconciliation* that this check contradicts.

## Honest reading
- **Three of the six come back against the manuscript**, and all three are claims that
  need rewording rather than analyses that need redoing.
  - **C2:** the two-way tie is half wrong. The calibrated null does beat the 1.0 kcal/mol
    cutoff (+0.290 [+0.094, +0.490]); the tie comes
    entirely from the fixed-se foil. The conjunctive pre-registered WIN rule still yields
    TIE overall, so the verdict stands and only the sentence has to change.
  - **C3:** the observed aggregate the prediction is said to match is 1.32 [0.81, 1.80] under the matched rule, against a
    predicted 0.85 [0.63, 1.08] — the
    intervals overlap but the prediction does not cover the point, and the observed value
    ranges 0.34–1.32 depending on a rule the text never
    names.
  - **C6:** the Wade et al. reconciliation runs through a functional they never compute.
    On the one they do compute, this data set is *more* conservative than perfect
    calibration, so the like-for-like comparison sharpens the disagreement instead of
    softening it. Only the shared-starting-coordinates caveat still argues for
    reconciliation, and it now has to carry that argument alone.
- **C1 and C5 hold up.** The cross-replicate correlation is real under cluster-honest
  inference (all three CIs exclude zero, permutation p at the resolution floor), and the
  shared-denominator artefact the referee suspected in the predicted-versus-observed
  check does not exist (0.648 exceeds all 200 null draws). Two
  qualifications survive: C1 does **not** establish the *concentration* of the
  correlation
  in the flagged systems at α = 0.05 (the manuscript asserts it; the label-permutation
  test is underpowered with 6 flagged systems), and C5's null is not
  centred on zero — its 95% range reaches +0.364 — so 'sampling noise would
  give zero' is the wrong reference and the measured null is the right one.
- **C4 quantifies a structural limitation** rather than testing a claim: 9 of 48 networks carry a single cycle, their reduced χ²
  swings by a median 14.2× between repeats, and
  every
  reproducibly flagged system is high-dof. Set-valued statements should be scoped to the
  well-determined networks.
- Nothing here was tuned. Every constant is in the script docstring and was fixed before
  the corresponding number was computed; C1 and C2 additionally assert that they
  reproduce the shipped numbers before adding inference to them.

