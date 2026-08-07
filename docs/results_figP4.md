# Results — Fig P4: QC calibration sweep under the REAL heterogeneous learned-sigma profile (peer-review item P4)

**Figure:** `figs/figP4_sigma_profile.{pdf,png}` · **Reproduce:** `make figP4`
(or `PYTHONPATH=src python figs/make_figP4.py`). Deterministic (fixed shuffle
`seed=20260808`, no bootstrap). Reuses `src/bar/sigma_profile.py` (Task 1:
`PROFILE_POINTS`, `rank_transfer`, `shuffled`) and `src/bar/qc.py`
(`gls_network`, `chi2_sf`, `benjamini_hochberg`) — no new modules, no new MD.

**Data provenance:** `data/openfe_replicates/combined_pymbar4_edge_data.csv`,
loaded exactly as `figs/make_figL.py` does — the public OpenFE replicate set.
Replicate 0 only. Each edge's `(a, b, ΔΔG, se)` comes from
`complex_repeat_0_DG/dDG` minus `solvent_repeat_0_DG/dDG` (se combined in
quadrature); each edge's overlap is the **minimum** of its
`complex_repeat_0_smallest_overlap` and `solvent_repeat_0_smallest_overlap`
(frozen: the worse leg limits the edge). Systems with fewer than 3 valid
replicate-0 edges, or with `gls_network(edges).dof < 1` (no independent
cycle), are dropped. **48 systems, 1143 edges** feed the sweep — the same
system set Fig L and Fig Cut use.

![Fig P4](../figs/figP4_sigma_profile.png)

## Motivation

Fig L panel B shrinks every edge's se by a single uniform `×0.15` to stand in
for an overconfident learned sigma, and reports ~88% of systems flagged.
Referees objected that a uniform shrink is near-mechanically forced — it
inflates every edge's χ² contribution by `1/0.15² ≈ 44×` regardless of that
edge's actual overlap, so the result says little about a real learned head,
whose miscalibration is heterogeneous and overlap-dependent (Fig A's swept
reported/true-se ratio ranges `0.09`–`0.20×` across overlap `0.26`–`0.78`).
This task pushes that **measured** profile through the identical
GLS + BH-FDR test, **per edge**, and reports what actually happens.

## Pre-registration (restated verbatim, frozen before any run)

- Replicate **0** only.
- Systems with `dof >= 1` only (independent-cycle networks).
- BH-FDR at `alpha = 0.05` (`bar.qc.benjamini_hochberg`).
- Shuffle seed **`20260808`** (`bar.sigma_profile.shuffled`).
- The profile curve is `bar.sigma_profile.PROFILE_POINTS` — verbatim from the
  committed `docs/results_figA.md` Panel-A table:
  `[(0.26, 0.20), (0.46, 0.11), (0.53, 0.09), (0.78, 0.15)]`. Never re-fit or
  re-tuned.
- Ranking is **global** (pool every edge across all 48 systems, rank-transfer
  once, then split back per system) — **not** per-system.
- Each edge's overlap = the **minimum** of its complex-leg and solvent-leg
  `smallest_overlap` at replicate 0.
- Verdict threshold: **`2×`**.
- **`DEGRADES`** iff the real-profile arm flags at least `2×` the calibrated
  arm's flagged count; **otherwise `COMPARABLE`**.
- `COMPARABLE` is a **fully acceptable, pre-registered outcome**. It would
  mean the "an overconfident learned sigma destroys the QC" claim is **not**
  supported by a realistic heterogeneous profile, and the manuscript must
  then soften it. Neither branch was to be treated as a failure to engineer
  around; the threshold, seed, arms, overlap definition, and profile were not
  to be adjusted in response to any observed number.

Four arms, all on the identical 48-system, 1143-edge network:

1. **calibrated sandwich (×1)** — no scaling; the existing Fig L / Fig Cut
   baseline.
2. **real learned profile (PRIMARY)** — each edge's se scaled by the profile
   ratio at that edge's global percentile rank of overlap.
3. **shuffled profile (control)** — the identical multiset of per-edge ratios
   from arm 2, randomly permuted (seed `20260808`) so the *association* with
   overlap is destroyed but the marginal distribution of ratios is unchanged.
4. **uniform ×0.15 (stress test)** — the original Fig L panel-B factor,
   retained for continuity but explicitly **not** the representative
   learned-head comparator.

## Complete verbatim stdout of the run

```
$ PYTHONPATH=src /Users/nikitapolomosnov/anaconda3/envs/fluor_screening/bin/python figs/make_figP4.py
wrote figP4_sigma_profile.(pdf|png) to /Users/nikitapolomosnov/PycharmProjects/fluor_screening/figs

[P4] profile points (overlap -> ratio): [(0.26, 0.2), (0.46, 0.11), (0.53, 0.09), (0.78, 0.15)]
[P4] edges=1143  systems=48  ratio range 0.090-0.200
[P4] calibrated sandwich (x1)               flagged  6/48 (12%)
[P4] real learned profile (PRIMARY)         flagged 42/48 (88%)
[P4] shuffled profile (control)             flagged 42/48 (88%)
[P4] uniform x0.15 (stress test)            flagged 42/48 (88%)
[P4] real-vs-calibrated ratio: 7.00x (DEGRADES needs >= 2x)
P4 VERDICT: DEGRADES
[P4] flagged by real profile but not calibrated: ['bace1', 'bace_p3_arg368_in', 'btk', 'cdk2', 'chk1', 'ciordia_retro', 'cmet', 'eg5', 'ephx2', 'factor_xa', 'galectin', 'hiv1_protease', 'hne', 'hsp90_2rings', 'hsp90_kung', 'hsp90_single_ring', 'irak4_s2', 'itk', 'jak2_set1', 'jak2_set2', 'jnk1', 'keranen_p2', 'liga', 'mcl1', 'mup1', 'ptp1b', 'renin', 'scyt_dehyd', 'shp2', 'syk', 't4_lysozyme', 'taf12', 'thrombin', 'tnks2', 'tyk2', 'urokinase']
[P4] flagged by calibrated but not real profile: []
```

Re-run via `make figP4` reproduces byte-for-byte identical printed numbers
(verified).

## Results

**Ratio range** (per-edge se-scaling factor after global rank-transfer):
**0.090–0.200×** — i.e. every edge's se is shrunk somewhere between 5× and
11× (variance inflated 25×–121×) relative to the calibrated sandwich, never
inflated. This is a narrow dynamic range: the whole profile sits well below
`×1`.

**Four arms' flagged counts and percentages** (n = 48 systems each):

| arm | flagged | percent |
|---|---|---|
| calibrated sandwich (×1) | 6/48 | 12% |
| **real learned profile (PRIMARY)** | **42/48** | **88%** |
| shuffled profile (control) | 42/48 | 88% |
| uniform ×0.15 (stress test) | 42/48 | 88% |

**Real-vs-calibrated ratio: 7.00×** (real 42, calibrated 6; `DEGRADES` needed
`>= 2×`).

**P4 VERDICT: DEGRADES.**

## Discordant systems

- **Flagged by the real profile but not by the calibrated sandwich (36
  systems):** `bace1, bace_p3_arg368_in, btk, cdk2, chk1, ciordia_retro,
  cmet, eg5, ephx2, factor_xa, galectin, hiv1_protease, hne, hsp90_2rings,
  hsp90_kung, hsp90_single_ring, irak4_s2, itk, jak2_set1, jak2_set2, jnk1,
  keranen_p2, liga, mcl1, mup1, ptp1b, renin, scyt_dehyd, shp2, syk,
  t4_lysozyme, taf12, thrombin, tnks2, tyk2, urokinase`.
- **Flagged by the calibrated sandwich but not by the real profile:**
  **none** (empty set) — the 6 calibrated-sandwich flags (`bace, brd4, cdk8,
  faah, hif2a, p38` per Fig L / Fig Cut) are a strict subset of the 42
  real-profile flags.

A supplementary check (not part of the pre-registered output, run to rule
out a bug given how closely the three shrinking arms agree — see "Honest
reading" below) shows the three shrinking arms are **not** identically the
same 42 systems: real vs. shuffled differ by exactly one system each way
(`real - shuffled = {bace1}`, `shuffled - real = {egfr}`); real vs. uniform
differ the same way (`real - uniform = {bace1}`, `uniform - real =
{egfr}`). The per-edge ratio assignment genuinely differs across arms
(1143/1143 unique ratio values in the real arm); the near-identical flagged
*counts* are a real feature of this profile's narrow dynamic range (see
below), not a computational artifact.

**Post-hoc tie-handling correction (not a re-run):** a review found
`rank_transfer`'s partial-tie handling assigned strictly increasing ranks via
`np.argsort` rather than sharing a percentile, which would matter only if the
input contained repeated overlap values. Since all 1143 real edges have
distinct overlaps (previous paragraph), this was a no-op here; `bar.sigma_profile.rank_transfer`
was fixed to use average-rank tie handling (`scipy.stats.rankdata`) and this
exact run was reproduced byte-for-byte identical (verbatim stdout above,
`P4 VERDICT: DEGRADES` and all four arms' flagged counts unchanged) after the
fix, confirming no reported number in this document moved.

## The rank-transfer caveat, in full

Fig A's controlled overlap sweep and OpenFE's `pymbar` `smallest_overlap` are
**different measurements on different scales**, and the profile transfer
must never be read as more than it is:

- **Fig A's normalized BAR overlap** (the profile's domain) spans
  **`0.26`–`0.78`** across its 4 swept points — a controlled Gaussian
  work-distribution Monte Carlo with `n_f = n_r = 20`.
- **OpenFE's `pymbar` `smallest_overlap`** (this task's per-edge overlap,
  minimum of the complex and solvent legs at replicate 0) spans
  **`0.0001`–`0.233`** across the 1143 real edges — a much lower, much wider
  range, measured by a different overlap statistic on real alchemical legs.

Because the two quantities are not interchangeable in raw value, the
transfer in `bar.sigma_profile.rank_transfer` is by **percentile rank, not
raw value**: an edge at the p-th percentile of the real (OpenFE)
overlap distribution receives the ratio the learned head showed at the p-th
percentile of the Fig A sweep. This preserves the profile's **ordering and
spread** — which is what the referees' objection is actually about (is the
miscalibration heterogeneous, does its structure matter) — without
pretending the two overlap scales are commensurable.

**This means the per-edge ratio assigned here is an approximation of how a
learned head's overconfidence would vary across these specific real edges,
not a direct prediction of the learned head evaluated on these edges.** No
learned head was ever run on the OpenFE benchmark; the profile is entirely
inherited from the independent Fig A controlled sweep and re-expressed by
rank. Any manuscript text drawing on this figure must state this transfer
explicitly and must not describe the per-edge ratios as "the learned head's
output on the OpenFE edges."

## Honest reading

- **The real, shuffled, and uniform arms land within one system of each
  other (42/48, 42/48, 42/48) despite genuinely different per-edge
  assignments.** This is not a coincidence of a narrow test: the profile's
  entire dynamic range (`0.090`–`0.200×`) inflates every edge's variance by
  at least 25× and at most ~121×, and — echoing the repo's Fig C finding that
  "weight choice is 2nd-order for GLS aggregate significance" — once a
  system's edges are *all* shrunk by a factor in this range, whether the
  network as a whole crosses the BH-FDR threshold is overwhelmingly
  determined by how much shrink there is on average, not by which specific
  edges get more or less of it within this range. Heterogeneity in *this*
  measured profile does not materially change which systems get flagged
  relative to a uniform stand-in of comparable central tendency — the
  discordant-system check above (`{bace1}` vs `{egfr}`) shows the real,
  shuffled, and uniform arms swap exactly one system each, not zero, so the
  per-edge structure is not entirely inert, but its effect on the flagged
  *count* is second-order at this sample size.
- **The pre-registered verdict is decided by the real-vs-calibrated
  comparison only** (42/48 vs 6/48, 7.00×), which is unambiguous and far
  past the `2×` threshold regardless of the shuffled/uniform arms' exact
  agreement with the real arm.
- **The 6 calibrated-sandwich flags are a strict subset of the 42
  real-profile flags** — the real profile does not disagree with the
  calibrated detector about which systems are wrong, it additionally flags
  36 more that the calibrated sandwich judges consistent. Combined with the
  narrow ratio range noted above, this is consistent with the real profile
  behaving, for the purposes of this specific test, close to (but not
  exactly) a uniform shrink of comparable magnitude to Fig L's `×0.15`
  stress test — which is itself the honest, pre-registered reading, not an
  engineered-around one.
- No numbers, arms, the seed, the overlap definition, or the threshold were
  adjusted after inspecting this result.

## What this changes in the manuscript

**The `DEGRADES` branch is selected.** The referees' objection — that a
uniform `×0.15` shrink is near-mechanically forced and says little about a
real learned head — is **not upheld** by this test: pushing the actual
measured, overlap-dependent profile through the identical GLS + BH-FDR
pipeline, per edge, still flags 42/48 (88%) of systems, 7.00× the calibrated
sandwich's 6/48 (12%), comfortably past the pre-registered `2×` bar. The
manuscript's existing claim — "an overconfident learned sigma destroys the
cycle-closure QC test" — **stands**, and Fig L's uniform `×0.15` panel can
continue to be cited as representative in direction and magnitude, not just
as a mechanically-forced stress test. What should change in the text:

- Cite this task's real-profile result (42/48, 7.00× calibrated) alongside
  or in place of the uniform `×0.15` result (also 42/48) as the primary
  evidence, since it answers the referees' specific objection and is no
  longer open to the "mechanically forced" critique — the shrink is
  per-edge and overlap-dependent, not uniform.
- Add the rank-transfer caveat (above) wherever this result is cited: the
  transfer is by percentile rank across two different overlap scales, so it
  approximates how the learned head's known overconfidence pattern would
  play out on these edges — it is not the learned head's own output on the
  OpenFE benchmark.
- Note, as an honest qualifier and not a walk-back, that the real, shuffled,
  and uniform arms are numerically close (all 42/48) because this
  particular profile's dynamic range is narrow enough that average shrink
  dominates the flagged count at this sample size; the manuscript should not
  claim the heterogeneous *structure* of the profile (as opposed to its
  magnitude) was shown to matter for this particular downstream test —
  only that the measured, non-uniform profile — not an artificially uniform
  one — still destroys the QC's selectivity.

Any Task 3 text drawing on this figure is contingent on exactly the numbers
and verdict above; nothing here was adjusted after the run.
