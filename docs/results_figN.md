# Results — Fig N: Paper-2 Increment-2 Step-0 collapse-stratum validity gate — **TERMINAL C**

**This is the load-bearing go/no-go for Paper-2 Increment 2.** Honest verdict: even after
aggregating the largest public breadth attempt available (LIT-PCBA + ChEMBL-diverse + BindingDB,
44 targets / 27,547 ligands), the pre-registered near-zero-similarity (`s < 0.15`) **collapse
stratum** contains only **6 query ligands across 3 fold-disjoint pocket clusters** — far below the
pre-registered power floor (≥30 queries / ≥8 fold clusters). This is **not** an ambiguous or
borderline result: the shape-null recovery gradient is cleanly monotone (high 0.941 → mid 0.588 →
orphan 0.370 → collapse 0.000) and the shape-null genuinely **does collapse** to random *within*
the collapse stratum (recovery@1 = 0.000, AUROC = 0.591, bootstrap CI = (0.0, 0.0) — all pass the
amended P2 collapse criterion). The gate still fails, because **P1 (power) fails outright**: 6 ≪ 30
and 3 ≪ 8. Three escalating breadth attempts (LIT-PCBA alone → +ChEMBL-diverse → +BindingDB) show
the collapse stratum barely moves (12 → 3 → 6), confirming the mechanism is structural, not a
data-gathering shortfall: **retrospective public bioactivity data cannot populate the orphan
regime at scale.** No Step-1 structure test is run (there is no valid, well-powered collapse
stratum to test it on).

**Code:** `src/screen/{bench_sources,stratify,fold_cluster,validity_gate}.py` +
`src/screen/sources/{litpcba,chembl_diverse,bindingdb}.py` (+ tests in
`tests/test_paper2_bench.py`) · `figs/make_figN.py` · `make figN`.

## Implementation note (read first — a real bug caught during this run)
The first foreground run of `make_figN.py` printed **n=782** for the collapse stratum, not the
expected ~6. Diagnosis: `_collapse_label` produces four strata named `"high"`, `"mid"`,
`"orphan"`, `"collapse"` — but `screen.validity_gate.verdict()` only ever reads whatever rows are
labeled `"orphan"` (it was written for Increment 1's single orphan cut). The original
Step-0-orchestrator code relabeled `collapse → "orphan"` to reuse `verdict()`, but this silently
**merged** the true collapse stratum (`s < 0.15`, n=6) with the pre-existing, larger `"orphan"`
stratum (`0.15 ≤ s < 0.35`, n≈776) under the same label, because both strata happen to share the
name `"orphan"`. Fixed by relabeling the pre-existing `"orphan"` rows out of the way
(`"orphan" → "_mid_orphan"`) *before* remapping `collapse → "orphan"`, so `verdict()`'s `"orphan"`
mask picks up only the true collapse rows. Independently re-derived `s < 0.15 = 6` /
`s < 0.10 = 4` directly from `stratify()` output (bypassing `verdict()` entirely) to confirm the
fix — this matches exactly what the prior commit (`a7c3b8d`, the BindingDB collapse-count
re-test) had reported from an ad-hoc script. **The corrected, load-bearing numbers in this
document are all post-fix**; the recovery **gradient** itself (computed directly from
per-stratum masks, not via `verdict()`) was never affected by this bug and is reported unchanged.

## What was tested
The fullest available aggregate: **44 targets, 27,547 unique ligands**
(`data/paper2_bench/triples_aggregate_bdb.parquet`, built by
`screen.bench_sources.build_triples(include_bindingdb=True)`), combining:
- **LIT-PCBA** (AVE-debiased, Increment 1's 15 targets / 7955 ligands),
- **ChEMBL-diverse** (14 additional fold-diverse targets spanning distinct Pfam families —
  kinases, proteases, nuclear receptors, ion-channel-adjacent enzymes — each with a holo PDB and
  ≥50 actives),
- **BindingDB** (15 further chemotype/fold-diverse targets, scraped from the live
  `ByUniProtids.jsp` search page — the only tractable per-target BindingDB path found in this
  environment; the documented `/rest/` JSON API times out and the bulk `BindingDB_All` dump is
  multi-GB).

Every ligand is a query with a known true target. Each query is stratified by `s` = max Tanimoto
(achiral ECFP4) to the union of **all 44 pockets'** active pools, **excluding the query's own
fingerprint by value** (leave-one-out, reusing `screen.stratify.stratify`):

| stratum | cut | n queries | shape-null recovery@1 |
|---|---|---|---|
| high | s ≥ 0.50 | — (majority of the 27,547) | 0.941 |
| mid | 0.35 ≤ s < 0.50 | — | 0.588 |
| orphan | 0.15 ≤ s < 0.35 | 776 | 0.370 |
| **collapse** | **s < 0.15** | **6** | **0.000** |

(`make_figN.py` uses its own `_collapse_label` cuts — high ≥0.50, mid ≥0.35, orphan ≥0.15,
collapse <0.15 — a stricter collapse cut than Increment 1's `deep_orphan` at 0.20, per spec §2.
The orphan-stratum n=776 is derived by exact arithmetic — the pre-fix run's merged
`verdict()`-reported n=782 for the relabeled `"orphan"` group, minus the independently-confirmed
collapse n=6 — not directly re-measured under `_collapse_label`'s cuts; the high/mid split is not
separately reported here since neither gate depends on it, only the monotone recovery@1 sequence
does. All four **recovery@1** values ARE the direct, unaffected `make_figN.py` output.)

## The three-way breadth-escalation history (the falsifier for "just add more data")
The pre-registered plan (spec §3) was: if the accessible aggregate's collapse stratum is thin,
escalate to BindingDB breadth. All three points were actually run:

| aggregate | targets | ligands | s<0.15 (or deep_orphan s<0.20 for Increment 1) | fold clusters |
|---|---|---|---|---|
| **Increment 1: LIT-PCBA alone** | 15 | 7,955 | 12 (`deep_orphan`, s<0.20 — Fig M) | not gated (P1 already failed on the wider `orphan` cut in Fig M) |
| **+ ChEMBL-diverse** | 29 | 26,952 | **3** (s<0.15; s<0.10 = 1) | **2** (targets GBA, HIVPR) |
| **+ BindingDB (fullest aggregate)** | 44 | 27,547 | **6** (s<0.15; s<0.10 = 4) | **3** |

(All three data points independently reproduced in this session directly from `stratify()`
output, not taken on faith from prior commit messages — see the Implementation Note above for
why that independent re-derivation mattered here.)

All three points sit far below the ≥30 queries / ≥8 folds power floor. Going from 15→29→44
targets and 7,955→26,952→27,547 ligands moved the collapse count only 12→3→6 — **not** a
monotone climb toward the power floor, and the ChEMBL-diverse step actually *shrank* it before
BindingDB partially recovered it. Both directions are visible: **more targets and more actives
densify chemical space rather than sparsify it.** Every new active pool added is itself a
drug-like, medicinal-chemistry-curated set, so it tends to sit in already-populated regions of
chemical space; each additional target makes it *less* likely that some random library ligand is
simultaneously far (`s < 0.15`) from *every* one of the (growing) set of pools. This directly
refutes the plan's working hypothesis in spec §3 ("more diverse islands ⇒ more near-zero query
ligands") — the mechanism runs the other way.

## Result: Step-0 gate — P1 (power) fails; P2 (collapse) actually passes

| stratum | n | recovery@1 | AUROC | 95% CI (recovery@1) |
|---|---|---|---|---|
| **collapse** | **6** | **0.000** | **0.591** | **(0.0, 0.0)** |

- **N = 44 pockets**, random baseline = 1/N = 0.023.
- **n_fold_clusters among collapse-query true targets = 3** (`screen.fold_cluster.n_disjoint_clusters`
  on the RCSB-Pfam-derived `fold` column).
- **Monotone-gradient sanity check: PASS.** recovery@1 strictly decreases high (0.941) → mid
  (0.588) → orphan (0.370) → collapse (0.000); the stratification is doing what it should, and the
  shape-null genuinely gets weaker as absolute similarity drops — consistent with, and extending,
  Fig M's high→orphan gradient.
- **P1 (power): FAIL.** 6 < 30 required collapse queries; 3 < 8 required fold-disjoint clusters.
  This is the load-bearing failure — the gate cannot even be run at adequate power.
- **P2 (amended, sample-size-adaptive collapse test): PASSES** on all three sub-conditions —
  recovery@1 = 0.000 ≤ 3/44 = 0.068; AUROC = 0.591 ≤ 0.60; 95% CI upper bound = 0.0 ≤ 1/44 = 0.023.
  Unlike Increment 1 (where the orphan-stratum shape-null stayed far above random), **the shape
  null genuinely disarms once similarity is pushed low enough** — but the regime where that
  happens is populated by only 6 ligands total, across only 3 distinguishable folds. The one
  useful signal (collapse is real) and the one fatal problem (there is almost no data there) are
  simultaneously true.
- **Verdict: FAIL / TERMINAL C** (P1 failure is sufficient by itself; §2 of the design spec treats
  power failure post-escalation as terminal, not iterable).

## Why (the honest mechanism)
1. **Aggregating more public actives densifies drug-like chemical space rather than sparsifying
   it.** Every source added (LIT-PCBA, ChEMBL-diverse, BindingDB) is itself a curated,
   medicinal-chemistry-relevant actives pool. Adding more such pools means a randomly drawn query
   ligand is *more* likely to be within Tanimoto 0.15 of *some* pool, not less — because drug-like
   space is finite and densely re-visited across unrelated programs (aromatic rings, common
   scaffolds, common substituents recur constantly). The three-point escalation (12 → 3 → 6)
   confirms this is a real, reproducible mechanism, not sampling noise: breadth barely moves the
   needle, and does not move it monotonically.
2. **The regime that would actually test structure — genuinely orphan, near-zero similarity to
   every known active — is a vanishingly small, saturating fraction of chemical space as sampled
   by any public bioactivity corpus.** ChEMBL, BindingDB, and LIT-PCBA are all built from
   published, patentable, medicinal-chemistry campaigns; by construction they cluster around
   tractable, drug-like chemotypes. The cage's actual regime (~0 Tanimoto similarity to any known
   active, an unusual difluoronaphthalenone+Michael-donor scaffold) is a genuine outlier relative
   to that distribution — which is exactly why it is scientifically interesting, and exactly why
   no public corpus, however large, is likely to contain enough neighbors-of-nobody to build a
   well-powered benchmark stratum around it.
3. **This deepens Fig H and the Increment-1 K1 (`docs/results_figM.md`), rather than merely
   repeating them.** Fig M showed retrospective public data can't disarm the ligand-similarity
   shortcut at a *moderate* similarity cut (0.20–0.35) because relative distinguishability
   survives. Fig N shows that even when you push the cut low enough that the shortcut *does*
   disarm (collapse stratum recovery@1 = 0.000), you can no longer find enough such ligands to
   test anything at power — the same underlying fact (public actives corpora are dense, curated,
   and drug-like) produces both failures from opposite ends.

## Sources reachable (scope of this verdict)
All three pre-registered sources were pulled and are committed: **LIT-PCBA** (Increment 1, reused),
**ChEMBL-diverse** (14 new fold-diverse targets via the EBI REST API + RCSB Pfam), and
**BindingDB** (15 new targets via the live `ByUniProtids.jsp` HTML search — the bulk TSV dump and
the documented `/rest/` JSON API were both impractical in this environment; see the docstring in
`src/screen/sources/bindingdb.py`). **The BindingDB breadth fallback specified in spec §3 is not
"still open" — it was tried, and it is the source of the strongest (though still failing) data
point (6 queries / 3 folds, up from 3 queries / 2 folds pre-BindingDB).** There is no further
pre-registered escalation path left in this plan; a materially different aggregation strategy
(e.g., deliberately hunting for
scaffold-disjoint or synthetically-inaccessible chemotypes rather than more of the same
medicinal-chemistry-curated actives) is a different, unregistered research design, not a rerun of
this one.

## Decision: the structure question is retrospectively undecidable — Paper-2's cage arm requires prospective validation, not a bigger public benchmark
- **The Step-1 structure gate (the free smina relative-scored test, `structure_gate.py` /
  `dock.py`) is NOT built and NOT run.** There is no valid, well-powered collapse stratum to run it
  on — running it on 6 queries across 3 folds would produce a number with no statistical meaning,
  and Plan 2 in the design spec explicitly conditions Step-1 on Step-0 passing.
- **This is a negative on the retrospective-public-benchmark *approach*, not a claim that
  structure-based target-ID is impossible.** The honest scope of this result: no public
  bioactivity aggregation strategy tried here (or plausibly triable at similar effort) reaches the
  near-zero-similarity regime at the power this gate pre-registered. A prospective panel (assay a
  curated, deliberately dissimilar target panel against the cage directly) or a fundamentally
  different construction (e.g., synthetically-designed hard-negative benchmarks built from
  scaffold-hopping rather than from existing bioactivity databases) remain open and are not ruled
  out by this result.
- **Consistent with, and reinforcing, both Fig H and Fig M.** Paper 1 ships as the methods+theory
  spine, unaffected. The cage/target-finding arm remains a **Paper 2** research question, and this
  increment closes the retrospective-benchmarking branch of it: the cage's near-zero-similarity
  regime is a genuine outlier that public corpora, however aggregated, do not populate at scale.
- **NIOCH operational debt stays separate and deliverable** (`docs/results_figH.md`): the cage can
  still be screened and returned as a ranked hypothesis list + recommended assays ("pending
  assays"), as an operational report, not a paper claim.

## Gate
`make check` green (ruff + mypy + full pytest incl. the Paper-2 benchmark tests, 128 passed).
Fullest-aggregate (44 targets / 27,547 ligands) collapse-stratum (`s<0.15`) shape-null: n=6,
recovery@1=0.000, AUROC=0.591, CI=(0.0, 0.0) — 3 fold-disjoint clusters. Monotone gradient
confirmed (0.941 → 0.588 → 0.370 → 0.000). P2 (collapse) passes; **P1 (power) fails** (6≪30,
3≪8) → **TERMINAL C → no Step-1 structure test; Paper-2's cage arm requires prospective
validation, not retrospective public benchmarking.**
