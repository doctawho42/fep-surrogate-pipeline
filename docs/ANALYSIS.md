# ANALYSIS — Phase 0 (think, pressure-test, verify, plan)

*Companion to `docs/paper1_plan.md` and `docs/bar_proofs.tex`. Produced before any
implementation code, per the Master Brief §4. Every numerical claim below was
reproduced independently (script: `scratchpad/verify_proofs.py`, **26/26 checks pass**).*

---

## 1. Reconstructed design (in my own words)

The pipeline is a **self-calibrating differentiable surrogate for relative binding
free energy (RBFE / ΔΔG)**, built around one structural idea: put a *fixed,
physically-exact BAR estimator* at the bottleneck and let a neural trunk feed it
per-λ-window work statistics, rather than regressing ΔΔG directly.

```
 ligand pair (3D, per enantiomer)
        │
   ┌────▼─────────────────────────┐
   │ frozen equivariant trunk     │  SO(3) features + parity-odd 0o channel
   │ (MACE / NequIP, odd irreps)  │  ← invariant 1: 0o needed for chirality
   └────┬─────────────────────────┘
        │  per-λ-window readout
   ┌────▼─────────────────────────┐
   │ BAR-bottleneck affinity head │  predicts (μ_k, log σ_k) per window k
   └────┬─────────────────────────┘
        │  μ_k, σ_k  (work distribution params)
   ┌────▼─────────────────────────┐
   │ FIXED BAR layer (Thm 1+2)    │  root-find ΔF; sandwich Var = B/I²
   │  O(1) backward via IFT       │  dΔF/dθ = (1/I)·∂_θ S
   └────┬─────────────────────────┘
        │  ΔΔG ± sandwich se   (+ ensemble σ_ens + Δ-head σ_Δ)
   ┌────▼─────────────────────────┐
   │ decision layer (Thm 3)       │  gauge-aware, cost-aware KG over FEP-edge graph
   │  Laplacian w_e = I_e²/B_e    │  effective-resistance acquisition
   └──────────────────────────────┘
```

**Why a BAR bottleneck at all?** Because it makes the surrogate's *uncertainty*
a physical quantity (sampling variance of an MLE) instead of a learned guess.
Calibration becomes "free" — read off the sandwich — *provided the predicted work
distributions are correct* (Assumption A1). That proviso is the whole risk surface
(see §4) and is exactly why the paper layers an ensemble (epistemic) and a Δ-head
(force-field bias) on top.

### The four theorems, restated

- **T1 — O(1) differentiable backward (information share).** With `ΔF̂` defined
  implicitly by `S(ΔF̂,θ)=0` and curvature `∂_d S = −I = −Σ p(1−p)`, the implicit
  function theorem gives `dΔF̂/dθ = (1/I)·∂_θ S`. For window means this is the
  **information share**: `dΔF̂/dμ_f = I_f/I`, `dΔF̂/dμ_r = −I_r/I`, and
  `|share_f|+|share_r| = 1`. The backward pass is *one division by the curvature* —
  no unrolling of the root-find. (Invariant #3.)

- **T2 — self-calibration via the sandwich, NOT 1/I.** `ΔF̂` is an M-estimator;
  its asymptotic variance is the **sandwich** `A⁻¹BA⁻¹ = B/I²` with `A=I` and
  `B = n_f·Var_f[p] + n_r·Var_r[p]`. The information equality `A=B` (which would
  give `Var=1/I`) holds for *prospective* logistic regression but **fails** under
  BAR's stratified two-sample design. `1/I` over-estimates the variance by a
  *regime-dependent* factor; it is correct only in the good-overlap limit.
  This is the corrected version of the originally-claimed result. (Invariants #1, #2.)

- **T3 — Fisher–resistance correspondence.** Treat measured edges
  `y_e = (φ_j−φ_i)+ε_e`, `ε_e ∼ N(0, V_e)` with `V_e = B_e/I_e²` the per-edge
  sandwich variance. The Gaussian posterior precision of the node potentials φ is the
  **weighted graph Laplacian** `L = Σ_e w_e (e_i−e_j)(e_i−e_j)ᵀ` with conductances
  `w_e = 1/V_e = I_e²/B_e` (NOT raw overlaps `I_e`). Then (i) contrast variance =
  effective resistance `Ω_ab = (e_b−e_a)ᵀL⁺(e_b−e_a)`; (ii) a single edge has `Ω=V_e`;
  (iii) Sherman–Morrison gives O(1) updates `Ω' = Ω/(1+gΩ)`, and KG = decision-weighted
  Ω-reduction; (iv) `ker L ⊇ span{1}`, so gauge/redundant directions have zero KG —
  acquisition is automatically gauge-invariant. (Invariant #4.)

- **T4 — chirality completeness.** The pseudoscalar `χ(M)=det[p₂−p₁,p₃−p₁,p₄−p₁]`
  is the parity-odd `0o` irrep (triple product). For any linear `A`,
  `χ(A·M)=det(A)·χ(M)`, so χ is SO(3)-invariant but flips sign under reflection.
  O(3)-invariant readouts (functions of `‖v‖`, `v_i·v_j`, pairwise distances)
  **collapse enantiomers** (`f(M')=f(M)`); only an *odd* contraction separates them.
  Pairwise distances + `sgn χ` is a complete SO(3)-invariant: O(3) loses exactly the
  chirality bit, and a `0o` channel restores it. (Invariants #5, #6.)

---

## 2. Verification results (independently reproduced)

Script `scratchpad/verify_proofs.py`, env `fluor_screening` (Python 3.13;
numpy 2.5, scipy 1.18, **pymbar 4.0.3**, networkx 3.6). **26/26 checks pass.**

### T1 — backward / information share
| quantity | reproduced | proof sheet |
|---|---|---|
| `ΔF̂` on `x_f={0,1}, x_r={−0.5}, M=0` | −0.595708 | −0.59571 |
| `I_f/I` (closed form) | 0.596825 | 0.596825 |
| `dΔF̂/dμ_f` (finite diff) ↔ `I_f/I` | 0.596825 | machine precision |
| `dΔF̂/d(x_r mean)` = `+I_r/I` | +0.403175 | (unified-coord) |
| `dΔF̂/dμ_r` = `−I_r/I` | −0.403175 | theorem sign |
| `‖share_f‖+‖share_r‖` | 1.0 | 1 |
| `∂_d S` = `−I` | −0.618660 | `−I` |

> **Sign convention nailed down (matters for Phase 1).** The head predicts unified-
> coordinate means `μ_k`. The IFT gradient w.r.t. *any unified-coordinate mean* is
> `+I_share/I`. The theorem's `−I_r/I` for `μ_r` is purely because the reverse
> coordinate enters *negated* (`x_j = −(μ_r+η_j)`). The custom backward must apply
> `+I_f/I`, `+I_r/I` to the per-window means in unified coordinates and let the
> coordinate map carry the reverse sign — **do not hard-code a minus sign in the head.**

### T2 — sandwich vs 1/I, and the pymbar disambiguation
Gaussian work model (forward `x∼N(m_f,s²)`, reverse `x∼N(m_r,s²)`, `m_f−m_r=s²`,
`ΔF=(m_f+m_r)/2`); separation in σ = `s`. `n_f=n_r=20`, 2000 replicates.

Headline (sep = 1σ): empirical `SD(ΔF̂)=0.162`, sandwich se `0.156` (ratio **0.96**),
naive `1/√I = 0.354` (ratio **2.19**). Matches the sheet (0.160 / 0.159 / 0.354).

**Overlap sweep — the Fig A result, confirmed in simulation:**

| separation | emp SD | **naive/emp** | **sandwich/emp** | pymbar-MBAR/emp | pymbar-BAR/emp |
|---|---|---|---|---|---|
| 1.0 σ | 0.159 | **2.23** | 0.98 | 1.00 | 0.99 |
| 1.7 σ | 0.283 | **1.50** | 0.98 | 1.01 | 0.99 |
| 2.4 σ | 0.445 | **1.23** | 0.96 | 1.01 | 0.98 |

The naive error is a **non-constant** factor (2.23 → 1.23): *no constant rescale of
`1/I` can calibrate it across overlap regimes* — the per-edge sandwich is required.
This is the falsifiable content of Fig A.

**pymbar (CLAUDE.md ambiguity resolved).** pymbar 4 removed the top-level `BAR`
symbol; the estimator is `pymbar.other_estimators.bar(w_F, w_R,
uncertainty_method=...)` with two methods:
- `'MBAR'` (default) — agrees with MBAR for two states exactly; **tracks empirical SD
  best across all overlaps (≈1.00)** → this is the canonical production sandwich.
- `'BAR'` — Bennett's original; slightly smaller, still sandwich-like (≈0.98).

**Crucially, *neither* pymbar method equals the naive `1/I`** (both are ~sandwich).
So "never report `1/I`" is fully consistent with "match pymbar". Mapping:
unified `x_i = W^f_i`, `x_j = −W^r_j` ⟹ call `bar(w_F=x_f, w_R=−x_r)`.

**Decoding pymbar's exact formulas** (read from source, `other_estimators.py`
L373–525). In our notation pymbar's `fF = 1−p(x_i)` and `fR = p(x_j)`, so its
`vartemp = n_f⟨p(1−p)⟩_f + n_r⟨p(1−p)⟩_r = I`. Hence, exactly:
- `Var_MBAR = 1/I − (1/n_f + 1/n_r)`  (default; = pymbar `'MBAR'`)
- `Var_BAR  = ⟨(1−p)²⟩_f/(n_f⟨1−p⟩_f²) + ⟨p²⟩_r/(n_r⟨p⟩_r²) − nrat`  (Bennett 1976; pymbar `'BAR'`)
- `Var_sandwich = B/I²`,  `B = n_f Var_f[p] + n_r Var_r[p]`  (Theorem 2)

**Are they the same?** Algebraically *no* (`B/I² = 1/I − nrat` would require
`B = I − I²·nrat`, false). But they are **asymptotically equal**: a controlled MC
study (`scratchpad/which_sandwich.py`, N∈{50,200,1000}, 4000 reps) shows all three
match the MC-truth variance to ~1% across separations 1.0–2.5σ. The ~4% gap I first
saw was a pure tiny-N (n=20) + poor-overlap artifact, where the estimator itself is
ill-behaved. So the proof sheet's "coincides with Bennett/Shirts–Chodera" is correct
*asymptotically*. **No contradiction.**

> **Production decision.** The layer's headline variance is the theorem's
> **`B/I²`** — it is what Paper 1 is named after and is always ≥ 0 (robust;
> `1/I−nrat` can go negative at extreme poor overlap and pymbar then returns nan).
> Invariant #1 ("match pymbar numerically") is discharged by *also* exposing the two
> pymbar closed forms, which equal `pymbar.bar(...)` to **machine precision** (proven
> by re-deriving the same expressions), plus an asymptotic-agreement test of `B/I²`
> vs MBAR. Mapping for the cross-check: `bar(w_F=x_f, w_R=−x_r)`.

### T3 — Fisher–resistance (all machine precision)
Triangle, conductances `w(1,2)=1, w(1,3)=0.5, w(2,3)=2`: `Ω₁₂ = 0.714286 = 1/1.4`
(series-parallel). `ker L = span{1}` (smallest eig ≈ 1e-16, null vec ∝ 1). Adding
edge `(1,2)` with `g=0.7`: `Ω' = 0.476190 = Ω/(1+gΩ)` (Sherman–Morrison) and the
direct Laplacian recompute agree to 0. Single edge: `Ω = 1/w = V_e`.

### T4 — chirality (all machine precision)
`χ = triple product`; `χ(A·M)=det(A)·χ(M)`; reflection (det −1) flips sign; rotation
(det +1) invariant; pairwise distances identical for `M` and its mirror (max diff 0)
while `sgn χ` flips. O(3) readouts are provably blind; `0o` separates.

---

## 3. Pressure-test: do the proofs hold? are the invariants consistent?

- **T1** is a clean IFT argument; the only subtlety is the envelope remark — if the
  downstream loss *is* the BAR log-likelihood, `∂_d ℓ=0` at the root kills the
  implicit term and the IFT gradient is unnecessary. Our ΔF̂ feeds a *downstream*
  prediction loss, so the IFT gradient is needed. Internally consistent.
- **T2** is the load-bearing correction. The original plan reportedly claimed
  `Var=1/I`; the sandwich corrects it. The proof (one-step expansion + delta method,
  `Var(S)=B` by two-sample independence and `Var[1−p]=Var[p]`) is standard
  M-estimation. Verified numerically to ≈1% — solid. The key conceptual point: BAR is
  a *stratified* (retrospective, fixed `n_f,n_r`) two-sample design, which breaks the
  information equality that holds for prospective logistic regression.
- **T3** depends on **Gaussian edge noise** (Assumption, stated in plan §9). For
  strongly non-Gaussian BAR residuals the Laplacian weights are approximate. Honest.
  The (i)–(iv) algebra is exact given the Gaussian posterior.
- **T4** is exact linear algebra; no scope caveats.

**Consistency of invariants:** #1/#2 (sandwich) feed #4 (Laplacian weights
`I_e²/B_e`), which is internally consistent — the *same* `B_e/I_e²` is the per-edge
variance in both. #5/#6 (chirality) are orthogonal to the BAR machinery. No conflicts.

### Fig A finding — direction of the `1/I` error (corrects a Corollary wording)
Confirmed three ways (controlled MC, population `I/B`, real benzene edges): `1/I`
over-estimates the BAR variance **most at HIGH overlap** (`I/B` up to ~17× as `p→½`)
and converges to correct at low overlap. This matches the Theorem-2 remark and
CLAUDE.md ("2.2× at high overlap, shrinking toward low overlap"), but the proofs
sheet **Corollary** says the opposite ("low-overlap regimes… naive worst") — its
directional wording is backwards (the math is fine). Implication (a *strengthening*):
the sandwich correction matters most for the well-overlapped, reliable edges an
active learner exploits; `1/I` would distrust the best edges (≤13× on real data) and
mis-weight the Laplacian. See `docs/results_figA.md`. **Not a kill** — the central
Fig A claim (sandwich calibrated; `1/I` a non-constant factor) holds robustly.

### Named weak points (carry into the risk ladder)
1. **Sandwich is aleatoric-only.** It is the sampling variance *under correct
   physics* (A1). It says nothing about OOD epistemic error — that is the *ensemble's*
   job (Fig B), and force-field bias is the *Δ-head's* job. Total belief variance
   `= B/I² + σ²_ens + σ²_Δ`. **Do not conflate Fig A (aleatoric, across overlap) with
   Fig B (epistemic, across domain distance).**
2. **A1 (correct specification)** is false in practice (force-field error). The frozen
   BAR layer cannot self-correct bias; this is *by design* delegated to the Δ-head.
   The sandwich's calibration guarantee is conditional on A1.
3. **Gaussian-edge-noise (T3).** Non-Gaussian residuals → approximate Laplacian.
4. **"Oracle" is imperfect.** Public FEP labels carry their own force-field error;
   all accuracy is framed *relative to FEP*, never to experiment.
5. **pymbar plug-in vs analytic sandwich** finite-sample gap (~3–4% at low overlap) —
   resolved by matching pymbar-MBAR as the reference (§2).

---

## 4. Build plan — falsifiable figures, in order (deliverable · acceptance · kill)

Discipline (CLAUDE.md): **ship Fig A end-to-end before scaffolding anything else.**
No trunk / heads / solver / generator modules until a shipped figure needs them.

### Phase 1 — BAR layer + **Fig A**  ← current target
- **Deliverable.** `src/bar/`: BAR root-find (Brent/Newton), sandwich variance
  `B/I²` matching pymbar-MBAR, overlap `I`, and an O(1) `torch.autograd.Function`
  with custom backward (`+I_f/I`, `+I_r/I` on unified means). Graph util: weighted
  Laplacian, effective resistance, Sherman–Morrison O(1) update, `ker L`.
- **Acceptance tests** (pytest, machine precision where stated):
  gradient = information share (the −0.59571 / 0.596825 instance); shares sum to 1;
  curvature `=−I`; sandwich se ≈ MC SD (ratio ∈ [0.95,1.05], n=20, ≥2000 reps);
  naive `1/I` ≈ 2.2× at high overlap + the overlap-sweep ratios; sandwich matches
  pymbar-MBAR within tolerance; effective resistance series-parallel; Sherman–Morrison
  O(1); `ker L` check.
- **Fig A.** Pull a handful of congeneric edges from a public FEP+ set; per edge
  compute overlap, `1/I` interval, sandwich interval; plot reported-se/true-se (or
  coverage) vs overlap for both. **Acceptance:** sandwich ≈ calibrated across overlap;
  `1/I` off by a *non-constant* factor.
- **Kill / gate.** Full suite green + Fig A regenerable by one command → proceed to
  Fig E. **If the sandwich is NOT calibrated → STOP and report with evidence.**

### Fig E — chirality completeness
- Minimal equivariant readout with a toggleable `0o` (triple-product) channel; designed
  enantiomer-pair test set. **Show:** even (O(3)) readout gives identical outputs
  (collapse); `0o` separates; ΔΔG error on chiral pairs with/without `0o`.
- **Kill:** if the even readout already separates → the theorem's premise is violated;
  investigate (do not proceed).

### Fig C — active-learning efficiency
- BAR-bottleneck head on a frozen public trunk + gauge-aware cost-aware KG over the
  FEP-edge graph. Baselines: random, standard BO (EI/UCB/qEI), GP-BoTorch qKG,
  MFBind-style multi-fidelity, ensemble+conformal. **Metric:** FEP-calls-to-top-k +
  regret-vs-budget. **Kill:** no fewer FEP calls than baselines on ≥2 sets.

### Fig D — gauge-aware identifiability
- KG=0 on gauge/redundant cross-class directions; budget concentrates on
  decision-relevant effective-resistance drops; ablate gauge-awareness → wasted budget.

### Fig B — decomposition & OOD *(the differentiator)*
- aleatoric(sandwich)+epistemic(ensemble)+Δ; coverage vs trunk-space domain distance;
  beats conformal (marginal-only). **Kill:** conformal matches per-edge adaptivity.

### Fig F — retrospective target-ID *(stretch, optional)*
- hidden-target recovery + P(binds) calibration; honest caveats.

**Central kill criterion (pre-registered, plan §0):** if NOT better on calibration
**and** NOT fewer FEP-calls-to-top-10 on ≥2 benchmarks → central claim dead; fall
back to the T1+T3 methods+theory paper. Do not rescue by tweaking a module.

---

## 5. Task breakdown (Phase 0 → Phase 1)

- [x] Read CLAUDE.md, paper1_plan.md, bar_proofs.tex; reconstruct design + theorems.
- [x] Re-derive & verify all four theorems numerically (26/26).
- [x] Cross-check sandwich vs pymbar; resolve which uncertainty method = sandwich.
- [x] Repo skeleton: `src/bar/`, `tests/`, `figs/`, pyproject (pinned), pytest, ruff,
      mypy, Makefile, CI workflow.
- [x] BAR layer: root-find + sandwich + mbar/bennett (= pymbar to machine precision) + overlap.
- [x] O(1) autograd.Function with custom backward (information-share gradient; gradcheck).
- [x] Graph util: Laplacian / effective resistance / Sherman–Morrison / kerL.
- [x] Port verification into `tests/` as theorem-numbered unit tests; `make check` green (26).
- [x] Source public FEP edges (benzene solvation, alchemtest) for Fig A real panel.
- [x] Produce Fig A + `docs/results_figA.md`; commit. **GATE → PASS.**

### Next
- [x] **Fig E (chirality):** `src/bar/chiral.py` (9 tests) + `figs/make_figE.py`.
      Even readout collapse = 0.0e+00 (exact); `0o` cuts enantiomer ΔΔG MAE 11.7×
      (1.667→0.142). Kill criterion not triggered. See `docs/results_figE.md`.
- [~] **Fig C (active learning) — minimal-first done; efficiency NOT demonstrated.**
      `src/bar/active.py` (gauge-aware cost-aware Sherman-Morrison KG, 6 tests) +
      `figs/make_figC.py`. Honest finding: decision-focused KG < A-optimal (boundary
      focus starves the global estimate); sandwich-vs-naive weighting ≈ null for
      *ranking* (its value is *calibration*, Fig A — unbiased GLS mean for any weights).
      Kill NOT triggered (conjunction; calibration passed); risk-ladder §8 fallback
      intact. Definitive test (integrated qKG + real networks + full baselines)
      **deferred pending decision**. See `docs/results_figC.md`. (Efficiency fallback ACCEPTED.)
- [x] **Fig D (gauge identifiability):** `figs/make_figD.py` (reuses `active.py`).
      Gauge (all-ones) KG = 1e-21 (exact 0, Thm 3(iv)); gauge-aware vs gauge-unaware
      regret@budget 0.010 vs 0.268 — gauge-awareness saves budget. PASS. See `docs/results_figD.md`.
- [x] **Fig B (OOD decomposition):** `figs/make_figB.py`. Decomposed per-edge
      uncertainty (sandwich+ensemble+Δ) vs split-conformal at equal ~0.90 marginal
      coverage: ours 2.2× sharper (3.76 vs 8.38), conditional max-bin-err 0.293 vs 0.734.
      Falsifier refuted; far-OOD all degrade (honest). PASS. See `docs/results_figB.md`.
- [x] **Fig F (target-ID, STRETCH → case study):** `figs/make_figF.py`. Real ChEMBL
      reverse screening, 2 regimes (8 diverse + 8 within-family aminergic GPCRs),
      scaffold-disjoint. Recovery 0.84–0.97 vs random 0.12; P(binds) ECE 0.008; but
      fingerprint model ≈ shape baseline in both regimes (falsifier holds) — ligand
      similarity is the signal, beating it needs structure/physics. Motivating negative;
      risk-ladder §8 case study. See `docs/results_figF.md`.
- [ ] *(optional strengthen Fig A)* real **binding** edge panel:
      `alchemtest.amber.load_bace_example` (BACE1 RBFE) and/or OpenFE
      IndustryBenchmarks2024 Zenodo `u_ln.txt` archives (3 repeats → true-replicate se).

---

## 6. Stack & environment (decided)
- **Project env:** conda env `fluor_screening` (Python 3.13). Phase-0/1 deps:
  numpy 2.5, scipy 1.18, pymbar 4.0.3, networkx 3.6, matplotlib 3.11, pytest 9.1,
  torch 2.12, ruff, mypy. **Lazy:** rdkit, e3nn/MACE, gpytorch/botorch, numpyro,
  netcal/torchcp are NOT installed until Fig E/C/B need them.
- **Git:** repo root is the catch-all `~/PycharmProjects`; `fluor_screening/` is
  untracked within it. Keep all commits scoped to `fluor_screening/`.
