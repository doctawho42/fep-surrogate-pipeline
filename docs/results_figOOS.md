# Results — Fig OOS: what survives of the QC localization claim, out-of-sample

**Figure:** `figs/figOOS_out_of_sample.{pdf,png}` · **Reproduce:** `make figOOS`
(`PYTHONPATH=src python figs/make_figOOS.py`). Deterministic (every resampling seeded,
`SEED = 20260810`). Data: `data/openfe_replicates/combined_pymbar4_edge_data.csv`
(OpenFE IndustryBenchmarks2024, public, 3 independent replicates per edge).

## Why this figure exists
The edge-removal *causal repair* panel was retracted: the statistic it drives down
(`sum_e z_e²`) is the quantity the removal order maximises (descending `z_e²`), so the
guided-vs-random contrast cannot fail, and a synthetic world with no localized systematic
error reproduces it. Both results below are constructed so that the quantity being
predicted is measured on data the selection never saw, which is the property the
retracted panel lacked. Nothing here is tuned: every constant is stated in the script
docstring and was fixed before any outcome was computed.

## Panel A — the replicate spread predicts the closure χ² it never saw

Write the reported per-edge variance as `V_e^rep = c_e² V_e^true`. Whitening the GLS fit
by the *reported* errors gives `E[X²] = tr(M D)` with `M = I − H` the residual maker and
`D = diag(c_e^−2)`, so with `h_e = M_ee` the curl-leverage of
`bar.leverage.curl_leverage` and `Σ_e h_e = dof` (Theorem D1),

> `E[χ²ᵥ] = Σ_e h_e c_e^−2 / dof` — the curl-leverage-weighted **mean** of `c_e^−2`.

The replicate set supplies an estimate of `c_e^−2` that the closure test never uses:
`s_e² / se_e²`, with `s_e²` the unbiased sample variance of the replicate ΔΔG values
(`E[s²] = σ²` at any n ≥ 2). **Independence is enforced by rotation**: replicate 0's
observed χ²ᵥ is predicted from `s_e` computed on replicates 1 and 2 only, and the same
construction is rotated over the other two choices. That is the headline. The
all-three-replicate variants are reported for comparison but are **in-sample**, since the
predicted replicate then contributes to its own predictor.

| variant | predictor uses | independent? | ρ | p | n | LOO ρ range |
|---|---|---|---:|---:|---:|---|
| **headline**: predict replicate 0 | replicates 1, 2 | yes | 0.648 | 6.3e-07 | 48 | [0.627, 0.687] |
| rotation: predict replicate 1 | replicates 0, 2 | yes | 0.617 | 3.0e-06 | 48 | [0.598, 0.665] |
| rotation: predict replicate 2 | replicates 0, 1 | yes | 0.703 | 2.5e-08 | 48 | [0.684, 0.741] |
| rotation-averaged | the held-out pair | yes | 0.891 | 2.0e-17 | 48 | [0.885, 0.907] |
| predict replicate 0 | all three | **no** | 0.754 | 6.0e-10 | 48 | [0.739, 0.785] |
| predict replicate 1 | all three | **no** | 0.825 | 5.3e-13 | 48 | [0.816, 0.853] |
| predict replicate 2 | all three | **no** | 0.790 | 2.4e-11 | 48 | [0.777, 0.810] |
| rotation-averaged | all three | **no** | 0.892 | 1.8e-17 | 48 | [0.885, 0.908] |

The in-sample variants are uniformly higher, as they must be: the shared replicate
correlates the two axes directly. Only the headline row is quotable.

**Diagnostics.** Reported on both predictors, because the difference is large enough to
matter: a diagnostic computed on the in-sample predictor is not a robustness check on the
headline.

| diagnostic | headline (independent) | same on the in-sample predictor |
|---|---:|---:|
| ρ, all systems | 0.648 (p = 6.3e-07, n = 48) | 0.754 (p = 6.0e-10) |
| ρ, excluding the 6 flagged systems (n = 42) | 0.545 | 0.672 |
| partial ρ given dof | 0.600 | 0.722 |
| ρ, systems with dof ≥ 5 (n = 25) | 0.702 | 0.790 |
| leave-one-out ρ range | [0.627, 0.687] | [0.739, 0.785] |

**Aggregate.** On the headline (independent, replicate 0) predictor the pooled
`E[χ²ᵥ] = 0.85`, system-cluster bootstrap 95% CI
[0.63, 1.08] (2000 resamples of the 48 systems);
pooling all three rotations gives 0.87 [0.62, 1.12].
Read as a bar width, the reported se is **1.09×** the replicate
standard deviation [0.96, 1.26]:
on this functional the bars are close to calibrated, not conservative by a large factor.
That is *not* in conflict with any replicate-validated ratio reported elsewhere. As the
SI records, a leverage-weighted mean of `c_e^−2` and a ratio of typical magnitudes are
different functionals of a heavy-tailed per-edge distribution. The observed side of the
same functional (dof-weighted mean χ²ᵥ) is 1.32 on replicate 0
(1.37 over all three), above the predicted level; the excess is what
the flagged systems contribute.

**Flagged systems sit above the identity line**, which is the expected direction: a
reproducible systematic error inflates the closure residual but not the replicate spread,
so observed exceeds predicted exactly where the QC fires.

| system | predicted E[χ²ᵥ] (independent) | predicted (in-sample) | observed χ²ᵥ (rep 0) | obs / pred (independent) | obs / pred (in-sample) |
|---|---:|---:|---:|---:|---:|
| bace | 2.48 | 2.18 | 5.57 | **2.2** | 2.6 |
| brd4 | 1.92 | 1.69 | 15.40 | **8.0** | 9.1 |
| cdk8 | 0.70 | 0.73 | 2.68 | **3.8** | 3.7 |
| faah | 1.28 | 1.56 | 3.71 | **2.9** | 2.4 |
| hif2a | 1.69 | 1.39 | 2.53 | **1.5** | 1.8 |
| p38 | 0.93 | 1.16 | 2.40 | **2.6** | 2.1 |

## Panel B — out-of-sample localization

Select edges on one replicate, evaluate on the two held out, rotate over all three
choices of selection replicate. Arms, all of the same size *K*:

- `static` — top-*K* by |z| on the selection replicate (the one-shot ordering);
- `guided` — `bar.qc.repair_order`'s greedy removal set (the paper's actual repair rule);
- `leverage` — top-*K* by curl-leverage *h* (topology only, blind to the values);
- `matched` — the **same |z| ranks** taken in an unflagged donor system (10 donors drawn from the 34 eligible, each scored against its
  own random null): is this just what high-|z| edges do anywhere?
- `random` — uniform *K*-subsets, which also supply the null and the p-values.

`K = max(1, #removals repair_order needs on the selection replicate)`; testable systems need dof ≥ 2.

### The statistics were calibrated before use, and two failed
Under the sampling null `Var(z_e) = h_e`, so selecting on |z| preferentially selects
high-leverage edges whose held-out `z²` is larger **by topology alone**. The harness runs
the entire pipeline in a world with no systematic error (same graphs, same reported se,
resampled values; 120 realizations × 3 rotations × 5 systems, 1800 tests) and measures each statistic's realized
false-positive rate at nominal α = 0.05. A statistic is used only if its realized rate is
within a factor of two of nominal, i.e. |α̂ − 0.05| ≤ 0.05:

| statistic | realized α (per test) | realized α (Stouffer) | verdict |
|---|---:|---:|---|
| held-out consistency gain  $\Delta\chi^2_\nu$ | 0.054 | 0.042 | USED |
| leverage-normalised held-out $z^2/h$ | 0.067 | 0.031 | USED |
| raw held-out $z^2$ | 0.198 | 0.489 | **DISCARDED** (invalid) |
| signed cross-replicate product $z_s z_h$ | 0.231 | 0.103 | **DISCARDED** (invalid) |

**2 of the 4 statistics are discarded**, and they are shown
here rather than silently dropped. Raw held-out `z²` fires at 0.198 because of the leverage confound above; dividing
by *h* restores the nominal level (0.067). The signed
cross-replicate product fires at 0.231 for a different
reason: selecting on |z| makes the multiplier `z_s` large, so the selected set's
statistic has a much wider null distribution than the equal-size random draws it is
compared with, and it exceeds their upper tail far more often than 5% of the time even
though both are centred on zero. That is measured here, not assumed: an exploratory pass
before this script reported that statistic as passing at α̂ = 0.077, which this run does
not reproduce. The `leverage` arm in the tables below carries the same confound as an
explicit control arm.

### Result
Effect sizes are means over the testable systems; p is the Stouffer combination over
systems within a rotation. `rand` is **each arm's own** random reference: for `matched`
that is the donor systems' random draws, since the matched statistic is measured in the
donors, not in the flagged system.

| statistic | arm | rotation 0 | rotation 1 | rotation 2 |
|---|---|---|---|---|
| held-out consistency gain  $\Delta\chi^2_\nu$ | `static` | +1.21 (rand -0.02), p=1.3e-03 | +1.42 (rand -0.02), p=1.3e-03 | +0.89 (rand +0.01), p=4.4e-03 |
| held-out consistency gain  $\Delta\chi^2_\nu$ | `guided` | +1.62 (rand -0.02), p=5.9e-05 | +1.53 (rand -0.02), p=1.2e-04 | +1.38 (rand +0.01), p=2.4e-04 |
| held-out consistency gain  $\Delta\chi^2_\nu$ | `leverage` | -0.60 (rand -0.02), p=9.1e-01 | -0.57 (rand -0.02), p=6.9e-01 | +0.00 (rand +0.01), p=2.3e-01 |
| held-out consistency gain  $\Delta\chi^2_\nu$ | `matched` | +0.29 (rand +0.22), p=1.0e+00 | +0.19 (rand +0.10), p=9.8e-01 | +0.19 (rand +0.09), p=2.6e-01 |
| leverage-normalised held-out $z^2/h$ | `static` | +9.97 (rand +2.88), p=1.4e-04 | +9.22 (rand +3.31), p=7.3e-04 | +10.11 (rand +3.26), p=2.5e-04 |
| leverage-normalised held-out $z^2/h$ | `guided` | +7.39 (rand +2.88), p=8.6e-04 | +6.60 (rand +3.31), p=1.4e-02 | +9.53 (rand +3.26), p=6.0e-04 |
| leverage-normalised held-out $z^2/h$ | `leverage` | +1.76 (rand +2.88), p=8.3e-01 | +3.18 (rand +3.31), p=6.0e-01 | +3.34 (rand +3.26), p=2.1e-01 |
| leverage-normalised held-out $z^2/h$ | `matched` | +1.81 (rand +0.99), p=9.7e-02 | +1.36 (rand +0.67), p=1.3e-03 | +1.49 (rand +0.81), p=7.2e-05 |

Per (system, rotation), the `static` arm's held-out Δχ²ᵥ ranges -0.39 to +5.39 against -0.34 to +0.17 for equal-size random.

### Per-system detail

| system | dof | K (rot 0/1/2) | held-out Δχ²ᵥ, `static` (rot 0/1/2) | same, `random` | per-system p, Δχ²ᵥ (rot 0/1/2) | signed product on the selected edges | whole-network cross-replicate z corr |
|---|---:|---|---|---|---|---:|---:|
| bace | 14 | 4/6/3 | +4.61/+5.39/+2.76 | -0.29/-0.34/-0.15 | 0.024/0.005/0.030 | +10.07 | +0.55 |
| faah | 8 | 4/1/1 | -0.39/-0.36/+0.34 | +0.17/+0.10/+0.07 | 0.791/0.893/0.339 | +0.19 | -0.09 |
| cdk8 | 29 | 7/2/5 | +1.03/+0.31/+0.61 | -0.01/+0.00/+0.02 | 0.002/0.088/0.018 | +5.25 | +0.62 |
| hif2a | 19 | 3/5/4 | +0.35/+0.72/-0.07 | -0.01/+0.11/+0.13 | 0.076/0.056/0.556 | +1.62 | +0.21 |
| p38 | 22 | 3/3/2 | +0.46/+1.03/+0.80 | +0.03/+0.03/+0.00 | 0.117/0.007/0.054 | +4.83 | +0.56 |

### Riders — these are not optional
- **It is an aggregate result.** Only 2 of the 5 testable flagged
  systems reach individual significance (p ≤ 0.05 on Δχ²ᵥ in at least two of the three
  rotations). The claim is about the set, not about any single system.
- **`faah` runs the wrong way** (negative mean held-out Δχ²ᵥ). Read it with the last two columns: a system whose signed residuals do not reproduce
  across runs has nothing for a one-replicate selection to transfer.
- **`brd4` is untestable**, having one independent cycle: removing a cycle edge leaves no closure to measure, so it is excluded rather than scored.
- **One-shot |z| ordering vs `repair_order`:** the greedy refit (`guided`) beats the one-shot |z| ordering (`static`) on **3 of the 6** (valid statistic × rotation) cells. The two arms are the same size by
  construction, so this is a like-for-like comparison of the ordering against the greedy
  refit; the tables above carry both, and neither is asserted in advance.
- **The `matched` control separates the two statistics.** On leverage-normalised held-out $z^2/h$ the same |z| ranks beat their own random null inside **unflagged** donor systems too, in at least two rotations, so that statistic measures a generic property of high-|z| edges rather than anything specific to the flagged set. On held-out consistency gain  $\Delta\chi^2_\nu$ the matched control does not fire, so that is where the flag-specific claim lives.
- The three rotations are not independent of one another (they reuse the same three
  runs), so the three p-values are three views of one dataset, not three replications.

## Honest reading
- Panel A is a genuine out-of-sample prediction: the replicate spread and the closure χ²
  are computed from disjoint data, related only through a theorem, and they agree in rank
  (ρ = 0.648) and in aggregate level. That is the strongest statement this
  benchmark supports about the calibrated null being the right null.
- Panel B's flag-specific content is narrower than it first looks. The `matched` control
  says the leverage-normalised held-out `z²/h` result is a generic property of high-|z|
  edges, reproducible in unflagged systems too; it is the held-out consistency gain
  Δχ²ᵥ, where the matched control stays silent, that is specific to the flagged set.
- Panel B recovers a weaker version of the retracted claim: the QC's per-edge ordering
  carries information about *unseen* runs of the same protocol. It does not show that
  removing those edges improves agreement with experiment (`make figLcausal` measured
  that separately, and found a null), and the selection is still evaluated on further
  runs of the same force field and protocol, so it validates reproducibility of the
  flagged error, not its physical diagnosis.
- Neither panel is a prospective test. The decisive experiment, re-running a flagged edge
  with fresh sampling and confirming the cycle closes, remains the natural next step.

