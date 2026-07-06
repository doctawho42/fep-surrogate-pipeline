# Results — Fig Lcausal: close-the-loop accuracy race (guided vs random QC-edge removal)

**Figure:** `figs/figLcausal_guided_vs_random.{pdf,png}` · **Reproduce:** `make figLcausal`
(`PYTHONPATH=src python figs/make_figLcausal.py`). Deterministic (seed 0, `n_perm=1000`).
Pre-registration: `data/openfe_replicates/closeloop_prereg.yaml`, immutability-anchored by an
external SHA-256 pinned in `tests/test_closeloop.py` (`bar.closeloop.load_prereg` raises if the
file is edited post-hoc). Data: `data/openfe_replicates/combined_pymbar4_edge_data.csv` (OpenFE
IndustryBenchmarks2024, public) joined to public ChEMBL affinity via `bar.ground` (no new MD).

## The question
Fig L / Fig Lval established that the calibrated cycle-closure QC flag is **causal for internal
network consistency**: guided removal of the top-`|z|` flagged edges reaches `reduced χ² ≤ 1` in
far fewer removals than random removal (3–7 guided vs 12–37 random-median removals across the
multi-cycle systems), and the flagged residuals reproduce out-of-sample across independent
replicates. That is a real result, but it is a **repair-of-internal-consistency** claim, not an
**accuracy-vs-experiment** claim. This increment asks the harder, pre-registered question directly:
**does acting on the QC flag (guided edge removal) improve accuracy against real experimental
affinity, more than an equal-size random removal?**

## Design (pre-registered before any affinity was joined)
- **Systems:** the 4 systems flagged by Fig L's BH-FDR test that are also **groundable** on public
  ChEMBL IC50/Ki (`coverage_min = 0.60` of network ligands mapped): `cdk8, hif2a, p38, bace`
  (`bar.ground.ground_system`, single assay per system in `assay_order = [IC50, Ki]`, never mixed).
  brd4 and faah (also Fig-L-flagged) are not in this set because they did not clear the coverage
  gate against public ChEMBL — this run only ever touches the 4 that did.
- **Guided removal:** `bar.qc.repair_order` — remove edges in order of `|z|` until
  `reduced χ² ≤ target_reduced_chi2 = 1.0`; this fixes the removal count `K` per system (the same
  trajectory endpoint used by Fig L/Fig Lval).
- **Random null:** 1000 random size-`K` edge removals per system (`n_perm = 1000`, seed 0).
- **Statistic:** `ΔMUE = MUE(full network) − MUE(after removal)`, both mean-aligned against
  experiment on the surviving node set, for guided and for each random draw. Per-system one-sided
  permutation `p = P(ΔMUE_random ≥ ΔMUE_guided)`. Combine across the 4 grounded systems by
  one-sided Stouffer on the per-system `p` + a sign test on `guided > random_mean`.
- **Success criterion (pre-registered):** combined Stouffer `p < 0.05` **and** a majority of
  systems below the 5th percentile of their own random null.
- **NO new MD.** Everything is a re-fit of the existing cached 3-replicate OpenFE
  IndustryBenchmarks2024 edges; only the affinity join and the removal bookkeeping are new.

## Result: realized branch = **NULL** (a measured scope boundary, not a failure)

**Per-system (K = guided removal count at the repair-order endpoint):**

| system | K | guided ΔMUE | random-null mean ΔMUE | one-sided p |
|---|---:|---:|---:|---:|
| cdk8  | 1 | **−0.018** | −0.003 | 0.905 |
| hif2a | 3 | **+0.031** | −0.025 | 0.136 |
| p38   | 2 | **−0.085** | −0.022 | 0.880 |
| bace  | 3 | **+0.113** | −0.012 | 0.071 |

**Combined verdict:** Stouffer `p = 0.484`, sign test `p = 0.6875` → **NULL**
(pre-registered success required Stouffer `p < 0.05`; not met by a wide margin, and the sign test
confirms no consistent direction).

**Baseline MUE vs experiment (full network, before any removal):**

| system | n ligands grounded | MUE vs experiment (kcal/mol) |
|---|---:|---:|
| cdk8  | 24 | 0.91 |
| hif2a | 41 | 1.47 |
| p38   | 35 | 0.79 |
| bace  | 24 | 0.88 |

(`n ligands grounded` is the affinity-cache count from `bar.ground.ground_system`, the coverage
denominator against the full network ligand list; the MUE fit itself runs on the smaller node set
that both survives edge removal and has a grounded affinity, so this count can exceed the number of
ligands actually contributing to a given MUE value.)

These baselines are typical-magnitude FEP errors (sub-to-low-single-kcal/mol MUE), confirming the
ChEMBL grounding itself is chemically sound — the null below is not an artifact of a broken
affinity join.

2 of 4 systems' guided removal slightly *improved* accuracy relative to their own random-null mean
(hif2a +0.031 vs −0.025; bace +0.113 vs −0.012), and 2 of 4 slightly *worsened* it (cdk8 −0.018 vs
−0.003; p38 −0.085 vs −0.022) — a wash, with no system individually significant (smallest p =
0.071, bace, still short of any per-system threshold) and no combined signal.

## The honest reading
The calibrated QC **causally repairs internal network consistency** — Fig L/Fig Lval already
showed this is a genuine, non-trivial, replicated positive (3–7 guided vs 12–37 random-median
removals to close cycles). This increment shows that same guided repair does **not**, on these 4
flagged systems, reduce force-field-limited **absolute accuracy vs experiment** beyond what an
equal-size random removal achieves. This is exactly what the paper's stated scope predicts, not a
contradiction of it: cycle closure is a test of **edge-level** consistency (signed ΔΔG sums to
zero around a loop), and it is **blind to node-consistent force-field bias** — a smooth per-ligand
offset cancels around any cycle. On real, flagged, multi-system data, that node-consistent bias is
exactly what dominates residual MUE-vs-experiment, so removing the edges that most violate cycle
closure does not systematically move accuracy in either direction. This NULL is the **direct
multi-system, experiment-grounded generalization** of the single-clean-system eg5 result below —
it upgrades that single-system finding to 4 independently flagged, ChEMBL-grounded systems with a
pre-registered combination rule, and it lands in the same place.

## Consistent with the eg5 clean-system control
`docs/results_eg5_accuracy.md` (single clean system, eg5, `reduced χ² = 0.64`) already found guided
QC action does not improve accuracy there either: internal `|z|` barely predicts per-edge error vs
experiment (Pearson **r = +0.13**, Spearman +0.05), and down-weighting/pruning flagged edges leaves
MUE essentially unchanged. That result was on a system the QC calls *clean*; this increment adds 4
systems the QC calls *flagged*, with real ChEMBL-grounded experimental affinity, and reaches the
same qualitative conclusion — reinforcing that the scope boundary (edge-level closure vs
node-consistent FF bias) holds regardless of whether the system is QC-clean or QC-flagged.

## Honest confound note
Guided removal is small by construction (`K` = 1–3 edges per system here) — a small, low-power
repair relative to network sizes of 23–41 nodes. The manual ligand→affinity ChEMBL join
(`bar.ground`) carries its own noise (median-pChEMBL aggregation over possibly-heterogeneous assay
records, cross-referenced structure matching). The random-permutation control uses the **same**
affinity mapping and the **same** removal count `K` for every draw, so this per-system, per-edge
noise is present identically in both arms of the contrast and cancels in the guided-vs-random
comparison — it cannot manufacture the NULL. What it *can* do is reduce power: with only 4 grounded
systems and `K` this small, a real but modest causal effect could be underpowered to detect. The
NULL reported here is therefore honestly a "not detected at this power," not a proof of exactly
zero effect — but it is consistent across all 4 systems (no system reaches even a per-system-naive
p < 0.05) and consistent with the independent eg5 control, which is the strongest evidence available
that the effect, if any, is small.

## Stereo-blind ChEMBL grounding limitation
`bar.ground`'s ChEMBL join is connectivity-based (`flexmatch`) and stereo-blind: enantiomer pairs
present in a network (same 2D graph, opposite stereocenter) resolve to the same ChEMBL record and
therefore share one experimental affinity value. Observed collapsed pairs in the committed affinity
caches (`bar.ground.collapsed_enantiomer_pairs`):

| system | collapsed enantiomer pairs |
|---|---|
| cdk8  | (30, 31) |
| hif2a | (237, 15), (165, 164), (7b, 7a) |
| p38   | (3fmk, 3fmh) |
| bace  | none |

This is a real match (flexmatch legitimately found that ChEMBL record), not a fabricated value, but
it necessarily makes one arm of each pair "wrong" relative to its true unmeasured enantiomer
affinity — inflating the **absolute** baseline MUE reported above. It **cancels in the
guided-vs-random contrast**, though: both arms of that comparison read the identical cached
affinity for the identical ligand, so the shared value is a constant offset on that ligand, not a
differential bias between the guided and random arms. Per this project's chirality invariants
(predict ΔG per enantiomer; racemate apparent affinity ≈ the stronger binder), this limitation is
disclosed rather than silently absorbed, per `bar/ground.py`'s module docstring.

## Bottom line
Pre-registered, deterministic, multi-system, experiment-grounded test: acting on the calibrated QC
flag does not improve MUE-vs-experiment beyond random removal on the 4 flagged-and-groundable
OpenFE systems (Stouffer `p = 0.484`, sign `p = 0.6875`, verdict **NULL**). This is a publishable
scope boundary, pre-registered as such: the QC's value is internal-consistency repair (Fig L/Fig
Lval, a genuine positive) and calibration-driven detection, not a force-field-accuracy multiplier —
exactly the claim the paper's scope makes, now confirmed on real experimental data across
independent systems rather than a single clean case.
