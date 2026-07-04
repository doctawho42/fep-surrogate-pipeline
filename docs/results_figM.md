# Results — Fig M: Paper-2 Increment-1 orphan-benchmark validity gate — **VALIDITY_KILL (K1)**

**This is the load-bearing go/no-go for Paper-2 Increment 1.** Honest verdict: on the LIT-PCBA
(AVE-debiased) benchmark, the orphan-stratum ligand-shape null does **not** collapse toward random
even under correct leave-one-out evaluation — it recovers the true target at recovery@1 = 0.474
(AUROC 0.863), well above random (1/N = 0.067) and far above the amended P2 collapse thresholds.
The orphan cut **weakens but does not disarm** the ligand-similarity shortcut, so the benchmark
cannot retrospectively test whether structure adds anything beyond ligand shape. Well powered (698
orphan queries, 12 fold-disjoint clusters), so this is a **definitive K1**, not a thinness result.

**Code:** `src/screen/{bench_sources,stratify,fold_cluster,validity_gate}.py` (+ tests in
`tests/test_paper2_bench.py`) · `figs/make_figM.py` · `make figM`.

## Correction notice (read first — honesty/reproducibility)
The **first** run of this gate reported orphan recovery@1 = 0.964 and diagnosed the null as winning
"by construction." That number was a **self-match artifact**: `shape_score_matrix` scored each query
against pools that still contained the query's own fingerprint (every LIT-PCBA ligand is itself an
active of its target), so each query trivially self-matched at Tanimoto 1.0. Fixed in commit
`c329c33` — `shape_score_matrix` now value-excludes the query from every pocket's pool (leave-one-out),
with a regression test (`test_shape_matrix_self_excludes_query_leave_one_out`). The corrected
recovery is **0.474**, not 0.964. The **verdict is unchanged (still VALIDITY_KILL, still K1)**, but
the mechanism is subtler than the artifact suggested — see the diagnosis below. All numbers in this
document are from the corrected (leave-one-out) run.

## What was tested
The full LIT-PCBA (AVE-debiased) triple table: **15 targets, 7955 unique ligands**
(`data/paper2_bench/triples.parquet`, built by `screen.bench_sources.build_triples`). Every ligand is
a query with a known true target. Each query is stratified by `s` = max Tanimoto (achiral ECFP4) to
the union of **all 15 pockets'** active pools, **excluding the query's own fingerprint by value**:

| stratum | cut | n queries |
|---|---|---|
| high | s ≥ 0.50 | 4631 |
| mid | 0.35 ≤ s < 0.50 | 2614 |
| **orphan** | **0.20 ≤ s < 0.35** | **698** |
| deep_orphan | s < 0.20 | 12 |

The Step-0 gate is the **shape null**: no structure, no docking, no learned model — just
max-Tanimoto-to-pocket-actives (leave-one-out) as the per-pocket score. If this trivial baseline
already recovers the target well above random in the orphan stratum, the stratum is not actually
orphan with respect to the benchmark's own actives pools, and no downstream structure-based scoring
can be retrospectively validated as adding value there.

## Result: the shape-null does NOT collapse (P2 fails)

| stratum | n | recovery@1 | AUROC | 95% CI (recovery@1) |
|---|---|---|---|---|
| **orphan** | **698** | **0.474** | **0.863** | (0.437, 0.511) |
| high | 4631 | 0.790 | — | — |

- **N = 15 pockets**, random baseline = 1/N = 0.067.
- **n_fold_clusters among orphan-query true targets = 12** (of 14 distinct orphan targets;
  `ADRB2`/`OPRK1` share Pfam `PF00001`, `ESR1_ant`/`PPARG` share InterPro `IPR001628`), via
  `screen.fold_cluster.n_disjoint_clusters` on the RCSB-derived `fold` column.
- **P1 (power): PASS.** 698 ≥ 30 orphan queries and 12 ≥ 8 fold-disjoint clusters — well populated,
  not thin. Rules out an "inconclusive / not enough data" verdict.
- **P2 (amended, sample-size-adaptive collapse test): FAIL** on all three sub-conditions:
  - recovery@1 = 0.474 vs required ≤ 3/N = 0.200 (should be near random; is **7.1× random**),
  - AUROC = 0.863 vs required ≤ 0.60 (should be near-chance; is far above),
  - 95% CI lower bound = 0.437 vs required ≤ 1/N = 0.067 (the entire bootstrap distribution sits
    far above random — not a borderline call).
- **Verdict: VALIDITY_KILL.**
- The recovery **gradient** high 0.790 → orphan 0.474 confirms the stratification is real (orphan
  queries ARE harder for the shape null than high-similarity ones), but 0.474 is still far from the
  ~0.067 collapse target: weakened, not disarmed.

## Why (the honest diagnosis)
1. **The orphan cut is about absolute Tanimoto to the whole library, not target-distinguishability.**
   `s` measures a query's nearest neighbour across all 15 pockets' pools combined (self-excluded). A
   query can have low absolute similarity to everything (s in 0.20–0.35) while still being **relatively**
   closer to its own target's active cloud than to any other target's, because target-specific chemical
   series occupy distinguishable regions of chemical space even at modest absolute Tanimoto. Reverse
   screening needs only that relative signal, so it survives the absolute-similarity cut.
2. **The regime where similarity truly is noise (near-zero similarity) is too thin here.** That regime
   is the deep_orphan stratum (s < 0.20) — exactly the cage's actual regime (~0 similarity to any known
   active) — but LIT-PCBA yields only **12** such queries, far below the ≥30 needed to gate on. So the
   one stratum where the null might collapse cannot be tested for lack of data.
3. **Class-imbalanced target sizes make the null easier.** LIT-PCBA per-target active counts are wildly
   unbalanced (ALDH1 = 5363 … ADRB2 = 17), so a few large, chemically coherent pools are trivially
   separable from the rest regardless of any single query's absolute novelty.
4. **This mirrors and deepens Fig H** (`docs/results_figH.md`): a public, retrospective, curated actives
   corpus keeps every query close enough to *some* labelled cluster — via relative similarity — that the
   true orphan regime is never sampled at scale. LIT-PCBA's AVE-debiasing removes analogue bias *within*
   a target's train/validation split, not *between* targets; cross-target distinguishability was never
   what it was built to remove.

## Sources reachable (scope of this verdict)
Only **LIT-PCBA** was pulled (the ChEMBL supplement and BindingDB breadth fallback in spec §3 were
deliberately deferred, gated on exactly this kind of finding). That deferral does **not** rescue the
verdict: LIT-PCBA alone gives a well-powered orphan stratum (698 queries, 12 folds), so this is a
genuine P2 failure, not thinness. It remains relevant to scope: **this is a K1 specific to
LIT-PCBA-style, already-curated, per-target-coherent actives pools.** A benchmark built from
BindingDB/ChEMBL breadth with genuinely cross-target-confusable, scaffold-disjoint ligands *might*
behave differently — but that is future work, not a basis for softening today's verdict.

## Decision: the structure question is retrospectively undecidable on this benchmark — the cage stays a prospective case study, not a benchmarked claim
- **Increment 2 (the free smina orphan gate) does NOT get a green light.** Testing whether a
  structure-based score beats ligand-shape in this "orphan" stratum would run against a null that is
  not blind (shape recovers at 0.474 there): any structure "win" could not be cleanly attributed to
  structure vs residual similarity.
- **Consistent with, and reinforcing, Fig H.** Paper 1 ships as the methods+theory spine unaffected.
  The cage/target-finding arm is **not** rescued by this increment; it remains a **Paper 2** research
  question requiring either (a) a genuinely cross-target-confusable orphan benchmark (not yet built —
  the BindingDB breadth attempt, or a relative-ambiguity stratification, are the candidate next steps),
  or (b) prospective validation instead of retrospective benchmarking.
- **NIOCH operational debt stays separate and deliverable** (`docs/results_figH.md`): the cage can be
  screened and returned as a ranked hypothesis list + recommended assays ("pending assays"), as a
  report, not a paper claim.

## Gate
`make check` green (ruff + mypy + full pytest incl. the Paper-2 benchmark tests). Corrected
leave-one-out orphan-stratum shape-null recovery@1 = 0.474 / AUROC 0.863 (well-powered: 698 queries /
12 fold clusters), far above the amended P2 collapse thresholds → **VALIDITY_KILL (K1) → Increment 2
does not proceed on this benchmark; the cage stays a prospective NIOCH case study, no benchmarked claim.**
