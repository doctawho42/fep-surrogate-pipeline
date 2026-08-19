# Pre-registration provenance (repository record, not manuscript content)

The manuscript prints **what** was pre-registered and what the frozen criteria say. It does not
print the anchoring bookkeeping: the digest, the immutability semantics, and the defects in the
frozen file's own wording. That material is recorded here, where it belongs.

## The file

`data/openfe_replicates/hodge_prereg.yaml` — the design file for the decomposition increment
(H1 visible fraction, H2 influence-ranked repair race, H3 auditability map). Committed before the
commit that computes any quantity it governs.

## The anchor, and exactly what it guarantees

SHA-256:

```
16279af73d90e1e9b39f33046ba6ee432b80ad750a6e5687a0e05b3e54b01585
```

Asserted on every run by `tests/test_hodge.py`, so any edit to the file fails the suite.

**What this does and does not establish.** The digest lives in the same repository as the file it
anchors, so it detects accidental drift and nothing more. There is no third-party timestamp, no
external registry entry, and no notarisation of any kind. It is not evidence against backdating,
and no claim in the manuscript rests on it being such evidence.

## Three errors in the frozen file, left standing

The file is left exactly as frozen, because correcting it would change the digest that anchors it.
Its own text is wrong in three places. Where the manuscript describes the same objects, the
manuscript governs.

| in the file | correct |
|---|---|
| calls `\|z_e\| * sqrt((1-h_e)/h_e)` "the DFFITS form" | it is the classical influence form scaled by `sqrt(h_e)`, not DFFITS |
| header says immutability is "anchored by an **external** SHA-256", and calls the file "Immutable" | the anchor is in-repository; the only guarantee is drift detection |
| cites "Theorem 5" for the conserved budget, and a part "5(i')" | the conservation law is the estimation--detection theorem and the decomposition is a different one; there is no part (i'), the intended reference is part (i). The file's numbering also predates the move of the backward-pass theorem into the Supporting Information, so every theorem number it quotes is one higher than the article's |

## Criteria that live only in code

Three pre-specified criteria are recorded in the analysis code rather than in any YAML, and were
audited against the descriptions the article prints. The audit and its outcome are in the
Supporting Information, since a mismatch between a frozen criterion and its printed description is
a finding a reader needs, not bookkeeping:

- the curl-leverage falsifier rule (`figs/make_figLev.py`);
- the fixed-cutoff head-to-head criterion (`src/bar/detectors.py`, `paired_auc_bootstrap`);
- the `|z|`-rule conjunctive criterion (`src/bar/closeloop.py`, `system_effect`, and
  `figs/make_figHodge.py`, `repair_race`).

One discrepancy was found there and is reported in print: `closeloop_prereg.yaml` describes the
conjunct as "below the 5th percentile" while `delta_mue` is `MUE(full) - MUE(after)`, so an
improvement is positive and the code correctly tests the **upper** tail, above the 95th percentile.
The returned key `below_5pct` is named after the file rather than after the statistic. No result
depends on the discrepancy.
