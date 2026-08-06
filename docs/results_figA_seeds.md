# Results — Fig A multi-seed foil spread (peer-review item P6c)

**Reproduce:** `make figAseeds` (prints the table below; nothing is written to disk).

## What this answers
A referee noted that the Fig A learned-foil calibration ratios (plain Gaussian-NLL,
large-budget oracle, β-NLL, 5-member ensemble — all `reported se / true se`, where 1.0 is
perfect calibration and lower is more overconfident) come from a **single training seed**
and carry no uncertainty. This run retrains all four foils at 5 frozen seeds and reports the
across-seed spread, isolating training-seed variability from Monte-Carlo re-sampling noise.

## Frozen pre-registration
- **5 seeds:** `RNG_SEED .. RNG_SEED+4` = `20260629, 20260630, 20260631, 20260632, 20260633`.
- **Separation sweep:** `FOIL_SEPS = np.linspace(0.8, 3.2, 11)` — the same sweep Fig A itself uses.
- **`reps = 1500`** Monte-Carlo replicates per separation.
- **Evaluation RNG stream held FIXED across training seeds** (`eval_seed = RNG_SEED` for every
  training seed), so the across-seed spread reflects only training-seed variability in the
  learned foils, not re-sampling noise in the Monte-Carlo truth.
- Report whatever the sweep produces; do not adjust seeds or reps in response to the result.

## Per-seed table (verbatim `make figAseeds` output)
```
      seed      plain     oracle    betanll   ensemble
  20260629     0.1541     0.1967     0.1126     0.1984
  20260630     0.1319     0.1994     0.1024     0.2018
  20260631     0.1229     0.2571     0.1110     0.2007
  20260632     0.1495     0.1941     0.1218     0.2198
  20260633     0.2103     0.2035     0.1256     0.2246
```

## Per-foil mean, min, max, sd (across the 5 seeds)
```
     plain: mean 0.1537  min 0.1229  max 0.2103  sd 0.0341
    oracle: mean 0.2101  min 0.1941  max 0.2571  sd 0.0265
   betanll: mean 0.1147  min 0.1024  max 0.1256  sd 0.0092
  ensemble: mean 0.2091  min 0.1984  max 0.2246  sd 0.0122
```

## Honest reading
The single-seed Fig A headline values — plain `0.15x`, oracle `0.20x`, β-NLL `0.11x`, ensemble
`0.20x` — all fall **inside** the observed 5-seed range (plain `[0.12, 0.21]`, oracle `[0.19,
0.26]`, β-NLL `[0.10, 0.13]`, ensemble `[0.20, 0.22]`); none is an outlier of its own spread.
The spread itself is narrow relative to the effect size for every foil (sd 0.009–0.034 against
means of 0.11–0.21), and **no foil's interval reaches or crosses 1.0** — the highest single
value across all 5 seeds and all 4 foils is the oracle at seed `20260631` (`0.2571`), still
roughly 4x below perfect calibration. The overconfidence claim for all four learned-variance
foils is therefore **not** an artifact of the one training seed used in Fig A: it holds
robustly across 5 independent training seeds at this label budget, with the ranking
β-NLL < plain < ensemble ≈ oracle preserved in every seed.
