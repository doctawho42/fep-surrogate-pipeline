# Results — Fig M: Paper-2 Increment-1 orphan-benchmark validity gate — **VALIDITY_KILL (K1)**

**This is the load-bearing gate for Paper-2 Increment 1** (`docs/target_finding_plan.md`
extended per the amended spec). Honest verdict: on the LIT-PCBA (AVE-debiased) benchmark,
the orphan-stratum ligand-shape null does **not** collapse toward random — it recovers the
true target almost perfectly (recovery@1 = 0.964, AUROC = 0.997). The benchmark, as
constructed, cannot retrospectively test whether structure adds anything beyond ligand
shape: the shape null is never actually blind here. This is a **definitive K1
Validity-Kill**, not a thinness/inconclusive result — the orphan stratum is well powered
(698 queries, 12 fold-disjoint clusters), so there is no "not enough data" escape hatch.

**Code:** `src/screen/{bench_sources,stratify,fold_cluster,validity_gate}.py` (+ tests in
`tests/test_paper2_bench.py`) · `figs/make_figM.py` · `make figM`.

## What was tested
The full LIT-PCBA (AVE-debiased) triple table: **15 targets, 7955 unique ligands**
(`data/paper2_bench/triples.parquet`, built by `screen.bench_sources.build_triples`).
Every ligand is a query with a known true target (its LIT-PCBA folder). Each query is
stratified by `s` = max Tanimoto (achiral ECFP4) to the union of **all 15 pockets'**
active-ligand pools, **excluding the query's own fingerprint by value** (so a query is
never trivially "most similar to itself"):

| stratum | cut | n queries |
|---|---|---|
| high | s ≥ 0.50 | 4631 |
| mid | 0.35 ≤ s < 0.50 | 2614 |
| **orphan** | **0.20 ≤ s < 0.35** | **698** |
| deep_orphan | s < 0.20 | 12 |

The Step-0 "gate" is the **shape null**: no structure, no docking, no learned model —
just max-Tanimoto-to-pocket-actives used directly as the per-pocket score
(`screen.validity_gate.shape_score_matrix`). If this trivial baseline already recovers
the true target near-perfectly in the orphan stratum, the stratum is not actually orphan
with respect to the benchmark's own actives pools, and no downstream structure-based
scoring can be retrospectively validated as adding value there — the null is never blind.

## Result: the shape-null does NOT collapse (P2 fails decisively)

| stratum | n | recovery@1 | AUROC | 95% CI (recovery@1) |
|---|---|---|---|---|
| **orphan** | **698** | **0.964** | **0.997** | (0.950, 0.977) |
| high | 4631 | 0.971 | 0.998 | (0.967, 0.976) |

- **n_pockets (N) = 15**, random baseline = 1/N = 0.067.
- **n_fold_clusters among orphan-query true targets = 12** (of 14 distinct orphan
  targets; `ADRB2`/`OPRK1` share Pfam `PF00001` and `ESR1_ant`/`PPARG` share InterPro
  `IPR001628`, collapsing 14 targets to 12 disjoint fold clusters via
  `screen.fold_cluster.n_disjoint_clusters` on the `fold` column `build_triples` attaches
  via the RCSB Pfam/InterPro lookup).
- **P1 (statistical power): PASS.** 698 ≥ 30 orphan queries and 12 ≥ 8 fold-disjoint
  clusters — the orphan stratum is well populated, not thin. This rules out an
  "inconclusive, not enough data" verdict.
- **P2 (shape-null collapse, amended criterion): FAIL.** All three sub-conditions
  fail simultaneously and by a wide margin:
  - recovery@1 = 0.964 vs required ≤ 3/N = 0.200 (should be **near random**, is instead
    **14.4× random**)
  - AUROC = 0.997 vs required ≤ 0.60 (should be near-chance, is instead near-perfect)
  - 95% CI lower bound = 0.950 vs required ≤ 1/N = 0.067 (the *entire* bootstrap
    distribution sits far above random — not a borderline call)
- **Verdict: VALIDITY_KILL.**

## Why (the honest diagnosis)
1. **The orphan cut is about absolute Tanimoto to the whole library, not about
   target-distinguishability.** `s` measures the query's nearest neighbor across *all*
   15 pockets' actives pools combined (self excluded), including other actives of its
   *own* target. A query can have low absolute similarity to everything (s in
   0.20–0.35, "orphan" by the stratification cut) while still being closer to its own
   target's active-ligand cloud than to any other target's — because target-specific
   chemical series occupy distinguishable regions of chemical space even at modest
   absolute Tanimoto. Checked directly: the 371 ALDH1 orphan queries have mean s = 0.31
   (comfortably inside the orphan band) yet the shape null still recovers ALDH1 as their
   top hit essentially every time.
2. **Class-imbalanced target sizes make the shape null easier, not harder, in this
   corpus.** LIT-PCBA's per-target active counts are wildly unbalanced (ALDH1 = 5363,
   VDR = 655, ... ADRB2 = 17), so a handful of large, chemically coherent actives pools
   dominate the nearest-neighbor signal and are trivially separable from the other 14
   targets, regardless of how "orphan" (low-absolute-similarity) any individual query is.
3. **This mirrors the exact failure mode Fig H already diagnosed for the ChEMBL
   retrospective benchmark** (`docs/results_figH.md`): a public, retrospective,
   already-curated actives corpus keeps every query close enough to *some* known,
   labeled cluster that the "true orphan" regime (the cage's actual regime — ~0 similarity
   to any known actives, structurally novel scaffold) is never actually sampled. LIT-PCBA
   was chosen because it is the AVE-debiased, "orphan-honest" anchor per
   `bench_sources.py`'s docstring — and even there, the debiasing operates on
   train/validation splits *within* a target, not on making one target's actives
   indistinguishable from another's. Cross-target distinguishability was never the
   thing LIT-PCBA was built to remove.

## Sources reachable (scope of this verdict)
Per the amended Task-2 scope, **only LIT-PCBA** was pulled for this run — the ChEMBL
supplement and BindingDB breadth fallback mentioned in the original spec (§3) were
**deliberately deferred**, gated on exactly this kind of orphan-thinness finding. That
deferral is now moot for the purposes of a P1 (power) argument: LIT-PCBA alone gives
698 well-clustered orphan queries, so there is no thinness case to make. It remains
relevant to the *K1* diagnosis, however: **this is a K1 Validity-Kill specific to
LIT-PCBA-style, already-curated, per-target-coherent actives pools.** A benchmark built
from BindingDB/ChEMBL breadth with genuinely scaffold-disjoint, cross-target-confusable
ligands might behave differently, but that is future work, not a basis for softening
today's verdict. **Do not read this as "the benchmark was too thin" — it was not; P1
passed comfortably. The kill is P2, and P2 failed because the null is well-powered and
still doesn't collapse.**

## Decision (ratified): the structure question is retrospectively undecidable on this
## benchmark — the cage stays a prospective case study, not a benchmarked claim
- **Increment 2 (the free smina orphan gate) does NOT get a green light.** Building a
  structure-based scorer and testing whether it beats ligand-shape in the "orphan"
  stratum of this benchmark would be testing against a null that is not actually blind:
  ligand-shape already wins there by construction of the actives corpus, so any
  structure-based score would be evaluated on a benchmark that cannot distinguish
  "structure adds signal" from "structure inherits the same shape-driven separability".
- This is **consistent with, and reinforces, the Fig H finding**: Paper 1 ships as the
  methods+theory spine (Figs A–G + `paper_draft.tex`), unaffected. The cage/target-finding
  arm is **not** rescued by this increment either — it remains a **Paper 2** research
  question requiring either (a) a genuinely cross-target-confusable benchmark (not yet
  identified/built) or (b) abandoning retrospective benchmarking in favor of prospective
  validation.
- **NIOCH operational debt stays separate and deliverable**, as already decided in
  `docs/results_figH.md`: the cage can still be screened and handed back as a ranked
  hypothesis list + recommended assays ("pending assays"), independent of where any
  publication lands — but as a **report, not a paper claim**.

## Gate
`make check` green (ruff + mypy + full pytest suite, incl. the Paper-2 benchmark tests
in `tests/test_paper2_bench.py`). The orphan-stratum shape-null does **not** collapse on
LIT-PCBA (recovery@1 0.964 / AUROC 0.997, both far above the amended P2 thresholds,
computed on a well-powered orphan stratum of 698 queries / 12 fold clusters) →
**VALIDITY_KILL (K1) → Increment 2 (smina orphan gate) does not proceed on this
benchmark; the cage stays a prospective NIOCH case study, no benchmarked claim.**
