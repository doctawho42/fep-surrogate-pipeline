# Results — Fig B: uncertainty decomposition & OOD calibration (the differentiator)

**Figure:** `figs/figB_ood_decomposition.{pdf,png}` · **Reproduce:** `make figB`.
Deterministic (seed 7). 12-member deep ensemble; ~6 s on CPU.

![Fig B](../figs/figB_ood_decomposition.png)

## Claim tested (plan Fig B)
Total Belief variance = aleatoric (sandwich `B/I²`) + epistemic (ensemble) + Δ
(force-field). The **decomposed per-edge** uncertainty holds **conditional** coverage
across the OOD spectrum, while a marginal split-conformal band (flat width) does not —
even though both are conformal-calibrated to the same ~90% **marginal** coverage.
**Falsifier:** split-conformal matches per-edge conditional coverage.

## Fair comparison
The ensemble is trained on the in-distribution support `[-2,2]` only (so OOD
extrapolation is real). Calibration and test are **exchangeable** (both span the full
range), so split-conformal earns its ~90% marginal guarantee fairly. We then probe
conditional coverage across trunk-space domain distance. All three are conformal-style:
- **ours** `μ ± q·σ_total` (normalized/Mondrian conformal; `σ_total=√(σ²_a+σ²_ens+σ²_Δ)`)
- **CQR** `μ ± 1.645·σ_a` + flat conformal adjustment `Q` (the stronger, *aleatoric*-
  adaptive conformal baseline)
- **split-conformal** `μ ± q` (flat width)

## Result: PASS (falsifier refuted, including vs CQR)

| method | marginal coverage | mean interval width | max \|bin cov − 0.90\| |
|---|---|---|---|
| **ours (decomposed)** | 0.897 | **3.76** | **0.293** |
| CQR (aleatoric-adaptive) | 0.888 | 7.64 | 0.712 |
| split-conformal (flat) | 0.890 | 8.38 | 0.734 |

Two headline findings, both honest:
1. **Same marginal coverage, ~2× sharper than even CQR.** Ours, CQR, and split-conformal
   all hit ~0.90 marginal, but ours does it with intervals **2.0× narrower than CQR**
   (3.76 vs 7.64) and **2.2× narrower than split-conformal** — it *adapts* per edge
   (tight where easy, wide where hard) instead of a flat (or aleatoric-only) band.
2. **Much better conditional calibration.** Max deviation of per-bin coverage from 0.90
   is **0.293 (ours)** vs **0.712 (CQR)** vs **0.734 (split-conformal)** — the flat band
   over-covers in-distribution and under-covers OOD; **CQR adapts to aleatoric
   heteroscedasticity but its flat conformal adjustment cannot grow with the
   epistemic/Δ OOD blow-up, so it fails conditionally too**; ours stays close to 0.90.

**Panel B** shows *why*: aleatoric (sandwich) is roughly flat, while **epistemic
(ensemble) and Δ (force-field) grow with domain distance** and dominate OOD. The
decomposition is interpretable — you know whether uncertainty is sampling noise,
model uncertainty, or force-field error — which a conformal band cannot tell you.

## Honest scope (plan §9 — this is the test, not a guarantee)
- **Far-OOD all degrade.** At domain distance > 1.5, coverage is ours 0.769 /
  CQR 0.697 / split-conformal 0.701. Extreme extrapolation is genuinely hard; the
  ensemble + Δ do not perfectly capture it. Ours is best everywhere and degrades the
  most gracefully, but does **not** hold a perfect 0.90 in the farthest bins. Reported,
  not hidden.
- **The sandwich is aleatoric-only.** OOD coverage is carried by the ensemble
  (epistemic) and the Δ term (force-field) — exactly the division of labour in the
  proofs sheet's "honest scope" remark. Fig B confirms the decomposition delivers
  conditional coverage in a controlled setting; it is not a worst-case guarantee.
- **CQR included (the tougher baseline) and it confirms the argument:** CQR is adaptive
  to *aleatoric* heteroscedasticity (slightly sharper than split: 7.64 vs 8.38) but its
  conformal adjustment is a flat scalar, so it cannot model the *epistemic/Δ* OOD blow-up
  — its conditional calibration (0.712) and far-OOD coverage (0.697) are no better than
  split-conformal, and it is 2× wider than ours. The physical decomposition is what wins.
- Controlled 1-D regression standing in for the surrogate; the marginal-vs-conditional
  coverage gap of split-conformal is a real, well-established property, not an artifact.

## Gate
`make check` green (41 tests) **and** Fig B regenerable by one command (`make figB`).
Decomposed per-edge uncertainty beats **both split-conformal and CQR** on sharpness
(2.0–2.2×) and conditional calibration (2.4×); neither conformal baseline matches
per-edge adaptivity OOD → **PASS**.
