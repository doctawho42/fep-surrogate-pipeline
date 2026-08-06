# Results — Fig A multi-seed foil spread (peer-review item P6c)

**Reproduce:** `make figAseeds` (prints the table below; nothing is written to disk).

## What this answers
A referee noted that the Fig A learned-foil calibration ratios (plain Gaussian-NLL,
large-budget oracle, β-NLL, 5-member ensemble — all `reported se / true se`, where 1.0 is
perfect calibration and lower is more overconfident) come from a **single training seed**
and carry no uncertainty. This run retrains all four foils at 5 frozen seeds and reports the
across-seed spread, isolating training-seed variability from Monte-Carlo re-sampling noise.

**Post-commit correction (this revision):** a review of the first version of this run found
that the `ensemble` foil's 5 rows were **not** 5 independent retrains. `_train_mve_ensemble`
builds its 5 members from seeds `sd, sd+1, ..., sd+4`; the outer sweep passed the *same*
`sd = seed0 + i` used for the single-net foils, so consecutive rows' member sets were a
sliding window sharing 4 of 5 networks (`{629-633}, {630-634}, {631-635}, {632-636},
{633-637}`), falsely narrowing the reported ensemble spread. **Only the ensemble call was
changed**: row `i` now trains its 5 members from a disjoint block, base seed
`seed0 + 5*i` (`{629-633}, {634-638}, {639-643}, {644-648}, {649-653}` — no shared networks
across rows, 25 distinct member networks total). `plain`, `oracle`, and `betanll` use the
same single-net-per-row scheme as before (`sd = seed0 + i`) and were **not touched**; re-running
confirms all 15 of their values reproduce bit-identically (see below). Only `ensemble`
changed, and it changed substantially — see "Honest reading" below for what the corrected
spread shows.

## Frozen pre-registration
- **5 seeds:** `RNG_SEED .. RNG_SEED+4` = `20260629, 20260630, 20260631, 20260632, 20260633`
  (row index for `plain`/`oracle`/`betanll`; `ensemble` maps row `i` to member-seed block
  `seed0 + 5*i .. seed0 + 5*i + 4`, disjoint per row — see correction note above).
- **Separation sweep:** `FOIL_SEPS = np.linspace(0.8, 3.2, 11)` — the same sweep Fig A itself uses.
- **`reps = 1500`** Monte-Carlo replicates per separation.
- **Evaluation RNG stream held FIXED across training seeds** (`eval_seed = RNG_SEED` for every
  training seed), so the across-seed spread reflects only training-seed variability in the
  learned foils, not re-sampling noise in the Monte-Carlo truth.
- Report whatever the sweep produces; do not adjust seeds or reps in response to the result.

## Per-seed table (verbatim `make figAseeds` output, re-run after the ensemble-seed fix)
```
      seed      plain     oracle    betanll   ensemble
  20260629     0.1541     0.1967     0.1126     0.1984
  20260630     0.1319     0.1994     0.1024     0.2053
  20260631     0.1229     0.2571     0.1110     0.2096
  20260632     0.1495     0.1941     0.1218     0.2655
  20260633     0.2103     0.2035     0.1256     1.2200
```

## Per-foil mean, min, max, sd (across the 5 seeds)
```
     plain: mean 0.1537  min 0.1229  max 0.2103  sd 0.0341
    oracle: mean 0.2101  min 0.1941  max 0.2571  sd 0.0265
   betanll: mean 0.1147  min 0.1024  max 0.1256  sd 0.0092
  ensemble: mean 0.4198  min 0.1984  max 1.2200  sd 0.4481
```

## `plain` / `oracle` / `betanll` reproduced bit-identically
All 15 values (3 foils × 5 seeds) are byte-for-byte identical to the pre-fix run: `plain`
`[0.1541, 0.1319, 0.1229, 0.1495, 0.2103]`, `oracle` `[0.1967, 0.1994, 0.2571, 0.1941,
0.2035]`, `betanll` `[0.1126, 0.1024, 0.1110, 0.1218, 0.1256]` — unchanged to 4 decimal
places at every seed. This is expected: their training call (`_train_mve(seed=sd, ...)`)
was not touched, the row seeds and the fixed evaluation stream are unchanged, and each is a
single net with no cross-row seed reuse. Row `i=0`'s `ensemble` value is also unchanged
(`0.1984`), because member block `{629-633}` is identical under both the old and the
disjoint scheme when `i=0`.

## Honest reading — the ensemble spread is now wide, and one value crosses 1.0
The corrected ensemble spread is **not** a tight cluster near 0.20-0.22 as the pre-fix
(overlapping-seed) run suggested. With genuinely independent member sets it ranges from
`0.1984` (seed `20260629`) up to `1.2200` (seed `20260633`, member block `20260649-20260653`)
— mean `0.42`, sd `0.45`, more than double the mean. **The seed-`20260633` value crosses
1.0** (nominally *overestimating* true se by 22%), which we flag explicitly per the
pre-registration's instruction to report whatever the sweep produces rather than adjust it.

Diagnosis (post-hoc, does not change the reported numbers): the block `20260649-20260653`
was never sampled by the pre-fix sliding-window scheme (whose member seeds only ever reached
`637`), so this instability was invisible before the fix. Retraining that specific block in
isolation shows one member (seed `20260649`) predicts a wildly inflated variance
(`sigma ≈ 8.3`, vs `~0.15-0.35` for its 4 co-members) at the sweep's highest separation
(`s=3.2`), which sits at the edge of the members' training-separation range
(`s ~ U(0.6, 3.0)`). A single such outlier member dominates the ensemble's `mean(sigma_m^2)`
aleatoric term at that separation and pulls the row's sweep-averaged ratio far above the
other rows. This is a real, reproducible property of Gaussian-NLL MVE training at a
realistic (`n_train=200`) budget — occasional catastrophic single-member variance
mis-calibration near the edge of the training distribution — not a bug in the fix or the
sweep.

**Manuscript-facing implication:** the prior claim "no foil's interval reaches or crosses
1.0" (and the prior "ensemble ≈ oracle" / "ensemble ~0.20-0.22×, narrow spread" framing) no
longer holds for the `ensemble` foil once measured with truly independent member sets. The
`plain`, `oracle`, and `betanll` columns are unaffected and still support the original
overconfidence claim exactly as before (none reaches 1.0; ranges as listed above). The
`ensemble` foil's 5-seed spread is wide and includes a value at/above calibration; text
citing a narrow "ensemble reaches only the large-budget floor" spread should be corrected to
reflect this — the ensemble's per-seed reliability is itself unreliable at this training
budget, driven by occasional single-member variance blow-up rather than the ensemble as a
whole being well-calibrated. The ordering "β-NLL < plain < ensemble ≈ oracle" claimed for
every seed does **not** hold at seed `20260631`, where oracle (`0.2571`) sits well above
ensemble (`0.2096`, a ~23% relative gap in the corrected run — it was `0.2007` and a ~28% gap
pre-fix), while ensemble slightly exceeds oracle at seeds `20260629` and `20260630`, and
greatly exceeds it at `20260632` and `20260633` (the outlier row). The accurate statement is:
β-NLL is the lowest (most overconfident) foil at every seed; plain and oracle interleave
seed-by-seed; ensemble is the least predictable foil, ranging from tied-with-plain to
crossing 1.0.
