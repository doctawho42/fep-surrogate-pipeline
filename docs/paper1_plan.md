# Paper 1 — Experimental Plan
### Self-calibrating differentiable free-energy surrogate for sample-efficient active learning

*Companion to the proofs sheet (`bar_proofs.tex`, v2). All calibration is measured against the **sandwich** variance `B/I²`, not `1/I`.*

---

## 0. Claim and falsification

**Minimal falsifiable claim.** On held-out congeneric series from public relative-binding-free-energy (RBFE) benchmarks, a surrogate trained through the differentiable BAR bottleneck — predicting the **sandwich** variance — achieves
(i) calibration (ECE / interval coverage) at least as good as a deep-ensemble + conformal surrogate, and
(ii) as the model in gauge-aware, cost-aware Knowledge-Gradient active learning, reaches top-*k* binders in **fewer FEP evaluations** than (a) ensemble + conformal under standard BO and (b) random / standard acquisition.

**Kill criterion (pre-registered).** If the method is *not* better on calibration **and** *not* fewer FEP-calls-to-top-10 on **≥2** benchmarks → the central methodological claim is dead; we do **not** rescue it by tweaking a module. Fallback = Theorems 1+3 (the O(1) differentiable BAR layer + the Fisher–resistance gauge-aware active learner) as a standalone methods+theory paper.

---

## 1. Datasets

**RBFE with precomputed FEP labels** (no FEP generated from scratch — the surrogate learns from existing public FEP; active learning reveals held-out precomputed labels as the "oracle"):
- **Wang et al. 2015 (Schrödinger FEP+)** — 8 targets, congeneric series. The field-standard.
- **Schindler et al. 2020 (Merck)** — 8 targets, larger / harder chemistry.
- **ToxBench** (arXiv 2507.08966) — AB-FEP labels, ERα; large single-target set for a within-target deep dive + absolute-FE check.
- **OpenFE / OpenFF** open-source benchmarks — reproducibility anchor.
- Construction / curation per **best-practices** (arXiv 2105.06222): congeneric definition, no analog leakage.

**Retrospective target identification** (separate, riskier figure): known ligand→target pairs from ChEMBL/BindingDB; pocket library from PDBbind apo/holo + AlphaFold structures. Hide the target, attempt recovery.

---

## 2. Splits and leakage control
- **Scaffold-disjoint** (Bemis–Murcko) → within-target generalization.
- **Target-disjoint** → cross-target transfer (for the unification/transfer claim).
- **Edge-level holdout** → the active-learning simulation reveals FEP edges sequentially.
- ≥ **5 random seeds**; report mean ± bootstrap CI. Pre-register the primary metrics (ECE + FEP-calls-to-top-10).

---

## 3. Methods and baselines

**Ours.** Frozen equivariant trunk (SO(3) + parity-odd `0o` channel) → BAR-bottleneck affinity head predicting `(μ_k, log σ_k)` per λ-window → fixed BAR layer → ΔΔG + **sandwich** variance `B/I²`; epistemic via shallow-head ensemble; learned Δ-residual head for force-field correction (own σ); decision layer = gauge-aware, cost-aware KG over the FEP-edge graph (weights = inverse sandwich variances).

**Baselines.**
- *Accuracy:* direct GNN / MACE surrogate predicting ΔΔG (no BAR layer); the public FEP reference itself (oracle ceiling).
- *Calibration:* deep ensemble + split-conformal / CQR; MC-dropout; last-layer Laplace; the AstraZeneca calibration recipe (arXiv 2407.14185).
- *Active learning:* random; standard BO (EI / UCB / qEI) on the same surrogate; GP-BoTorch qKG; **MFBind** multi-fidelity (arXiv 2402.10387); batch-BO foundation surrogate (arXiv 2511.10590).
- *Δ-learning positioning:* QuantumBind-RBFE (2501.01811), Exscientia ML/MM end-state corrections (2410.16818) — prior art our Δ-head builds on.

**Novelty vs closest prior art (state explicitly in the paper).** MFBind does multi-fidelity but **not** BAR-bottleneck calibration or gauge-aware identifiability. The AstraZeneca study does calibration but **not** physics-grounded / differentiable / per-edge. QuantumBind & Exscientia do Δ-corrections but **not** a calibrated differentiable surrogate bottleneck driving active learning. Contribution = the **integration** + two theorems (O(1) backward; Fisher–resistance gauge-aware AL) + sandwich self-calibration.

---

## 4. Metrics (reported variance = sandwich)
- **Accuracy:** MAE, RMSE, Kendall-τ, Spearman on ΔΔG vs FEP labels.
- **Calibration:** ECE, NLL, prediction-interval coverage at 90 %, sharpness (mean width), reliability diagrams. *Model variance = sandwich `B/I²`.*
- **Active learning:** FEP-calls-to-top-*k* (k = 1, 5, 10), simple-regret-vs-budget curves, cost-weighted regret, area-under-regret-curve.
- **Target-ID:** top-*k* recovery, AUROC of recovery, calibration of P(binds).
- **Transfer:** Δ(accuracy, calibration) on target-disjoint split with vs without the shared trunk.
- **Statistics:** paired comparisons across seeds/targets (Wilcoxon signed-rank), bootstrap CIs, effect sizes.

---

## 5. Headline experiments / figures

| Fig | What it shows | Supports | Falsifier |
|----|----|----|----|
| **A** "target the sandwich" | reported-se/true-se for `1/I` vs sandwich across **overlap** regimes (controlled sim + real FEP edges binned by overlap). Sandwich ≈ 1 everywhere; `1/I` over-conservative by a *varying* factor (≈1.2–2.3× in se) → no constant rescale fixes it → per-edge sandwich required. | Thm 2 (corrected) | sandwich not ≈ 1, or `1/I` error constant |
| **B** decomposition & OOD *(the differentiator)* | On test edges sorted by trunk-space domain distance: coverage held by aleatoric(sandwich)+epistemic(ensemble)+Δ; epistemic grows & dominates OOD; vs conformal (flat marginal band, no per-edge adaptation) and a single-number heuristic. | calibration-beats-conformal | conformal matches per-edge adaptivity |
| **C** active-learning efficiency | FEP-calls-to-top-*k*: ours vs ensemble+conformal+standard BO vs random vs MFBind; regret-vs-budget. | Thm 3 in action; main efficiency claim | no fewer FEP calls |
| **D** gauge-aware identifiability | KG = 0 on gauge / redundant cross-class directions; budget concentrates on decision-relevant effective-resistance drops; ablate gauge-awareness → wasted budget. | Thm 3(iv) | gauge-awareness gives no budget saving |
| **E** chirality completeness | designed enantiomer-pair test set: even (O(3)) readout identical (collapse); odd (`0o`) readout separates; downstream ΔΔG error on chiral pairs with/without `0o`. | chirality completeness theorem | even readout already separates |
| **F** retrospective target-ID *(stretch)* | top-*k* recovery + P(binds) calibration on hidden-target benchmark; honest caveats. | reverse-screening capability | recovery ≤ shape baseline |

---

## 6. Ablations (each tied to a theorem)
- Remove BAR layer (predict ΔΔG directly) → calibration degrades — **Thm 2**.
- `1/I` vs sandwich variance (Fig A) — **Thm 2 correction**.
- Remove `0o` pseudoscalar channel → enantiomer collapse (Fig E) — **chirality theorem**.
- Remove cycle-closure / Fisher–resistance prior → variance & accuracy on networks — **Thm 3**.
- σ-floor + overlap-regularizer on/off → training stability in low-overlap windows — **Thm 1 remark**.
- Δ-residual head on/off → OOD accuracy — **Thm 2 scope**.
- Gauge-awareness on/off → wasted budget (Fig D) — **Thm 3(iv)**.

---

## 7. Compute & feasibility (solo PhD)
- **No new FEP** for the core — precomputed public labels are the oracle.
- Trunk = a public equivariant model (MACE / NequIP), pretrained then frozen — no large from-scratch training.
- Active-learning experiments = cheap simulation over revealed precomputed labels.
- Only Fig F (target-ID) needs pocket-library processing — heavier, optional.
- Estimate: Fig A–E on a single GPU + modest CPU over a few months; Fig F is the stretch.

---

## 8. Risk ladder & fallback (what survives each failure)
- **Fig B fails (calibration ≤ conformal):** keep Thm 1+3, Fig A (sandwich correctness), Fig C/D (efficiency) → methods+theory paper on the differentiable layer + gauge-aware AL.
- **Fig C fails (efficiency null):** keep calibration + sandwich + theory.
- **Fig F weak:** drop to perspective / case study.
- **Transfer null:** report shared trunk as engineering, not a transfer claim.
- **Floor:** the four proven theorems + Fig A + Fig E are publishable as a focused methods+theory paper regardless of the application.

---

## 9. Honest residual risks
- The sandwich is the **aleatoric** (sampling) variance under a *correct* predicted work distribution; OOD epistemic error is the ensemble's job (Fig B is the test, not a guarantee).
- Thm 3 assumes Gaussian edge noise; strongly non-Gaussian BAR residuals make the Laplacian weights approximate.
- Public FEP labels carry their own force-field error — the "oracle" is imperfect; frame accuracy *relative to FEP*, never to experiment.
- `1/I`-vs-sandwich (Fig A) is about getting the **aleatoric** variance right, demonstrated across overlap regimes — it is **not** the OOD story (that is Fig B). Do not conflate them.
