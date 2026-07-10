# Fig AC (g-sweep) — autocorrelation robustness of the sandwich calibration

**Referee ask #1 [MUST]** (`docs/anticipated_referee_responses.md`): report the sandwich se at raw
`n` vs `n_eff = n/g`, and the panel-B calibration ratio under both. **No new MD** — reuses the same
alchemtest BACE1 RBFE complex-leg works as Fig A panel B.

Reproduce: `make gsweep` (→ `figs/figAC_gsweep.pdf`, `figAC_gsweep.png`). Deterministic
(`RNG_SEED = 20260710`; 400-sample bootstrap).

## What it tests

The sandwich `B = n_f Var_f[p] + n_r Var_r[p]` (Theorem "sandwich") assumes **independent** samples.
Real MD works are autocorrelated, so feeding raw samples over-counts the information and
under-estimates the variance (the estimator docstring, `src/bar/estimator.py`, is explicit about
this). The paper's protocol subsamples each leg to its statistical inefficiency `g`
(`paper_body.tex:172`). This panel quantifies the effect on real BACE1 data.

## Numbers (19 BACE1 binding edges, run 2026-07-10)

- Statistical inefficiency `g`: geometric-mean **1.48**, range [1.15, 2.16]; median count `n = 500 → 349`.
- **sandwich(n_eff)/truth = 1.00 [0.93, 1.07]** — calibrated; reproduces Fig A-B (which reported
  1.00 [0.96, 1.07]).
- **sandwich(raw)/truth = 0.83 [0.67, 0.97]** — overconfident, and equal to the
  `1/sqrt(g_bar) = 0.82` autocorrelation prediction.
- **sandwich(raw)/raw-boot = 1.00 [0.92, 1.05]** — the sandwich tracks its *own* (naive) bootstrap
  under raw `n` too: the closed form is internally consistent; the only quantity that must be right
  is the effective sample count.

Truth = correlation-aware reference = decorrelated (n_eff) iid bootstrap SD, the same reference used
in Fig A panel B.

## Message

The calibration **verdict is invariant to autocorrelation once `n_eff` is used**, as it is in every
real-data panel of the paper. Raw `n` would be overconfident by `~sqrt(g_bar)` (here 0.83×), but the
paper never uses raw `n`. This converts the single most-likely fresh-FEP-referee ask (a
demonstration that the i.i.d. sandwich survives real MD autocorrelation) from a revision trigger
into a pre-empted, in-manuscript SI panel — at zero new-MD cost and with no null risk.

SI text + figure: `docs/paper_si.tex` (`\label{fig:gsweep}`).
