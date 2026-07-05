# Cage prospective-validation hand-off protocol

**Status:** pre-registered prospective target-ID protocol for the difluoronaphthalenone "cage."
**Audience:** the executing partner lab (НИОХ) and anyone auditing the result afterward.
**Companion artifacts:** `data/cage/prospective_prereg.yaml` (the frozen forecast + decision
constants — the single source of truth for every number in this document), `src/screen/prospective.py`
(the decision scaffold the partner's summary statistics feed), `docs/superpowers/specs/2026-07-05-cage-prospective-loop-design.md`
(the design rationale), `docs/cage_assay_request.md` (the earlier, less formal assay request that
this protocol supersedes), and `docs/results_boltz_cage.md` (the in-silico record behind the
forecast below).

## 1. Header and scope

This is a **pre-registered prospective target-identification protocol**, not a validated target
claim. Nothing in Sections 2–9 asserts that the cage binds any named protein; the whole point of
the exercise is to find out, honestly, whether it binds anything at all, and if so what — with the
decision rule fixed in writing before any wet data exists, so that no one (including the authors)
can quietly move the goalposts once results come in.

The forecast that this protocol tests is frozen in `data/cage/prospective_prereg.yaml`. That file's
git commit **is** the timestamp of the forecast: whatever revision of that YAML is checked into the
repository at the moment the wet lab receives compound is the forecast being tested, and it is
never edited after that point. If a genuine correction is ever needed, it happens as a new,
separately committed file or a clearly marked addendum — never as a silent edit to the frozen
version.

The **un-blinding rule** is simple and absolute: predictions are frozen first (the prereg commit),
data is collected second, `score_forecast` (in `src/screen/prospective.py`) is run third against
the observations the partner returns, and the scorecard it produces is the result. No forecast
edits are permitted after the prereg commit, for any reason, including "the forecast turned out to
be trivially wrong in an embarrassing way" or "we thought of a better threshold after seeing the
data." If the decision constants in the prereg turn out to be miscalibrated in hindsight, that is
itself a finding to report, not grounds for revision.

## 2. Compounds and forms

Four species are tested, corresponding to the two ring-fusion stereocentres (one diastereomer plus
its mirror image) crossed with the as-synthesized O-acetate and its hydrolysis product, the free
alcohol (deacetyl). All four SMILES below are copied verbatim from `data/cage/prospective_prereg.yaml`
(the `species` block); do not resynthesize or re-derive them from any other document.

| species ID | enantiomer | form | SMILES |
|---|---|---|---|
| RR-OAc | R,R | acetate | `CC(=O)O[C@]12C[C@H](c3ccccc3C1(F)F)c1c([nH]c(=O)[nH]c1=O)O2` |
| SS-OAc | S,S | acetate | `CC(=O)O[C@@]12C[C@@H](c3ccccc3C1(F)F)c1c([nH]c(=O)[nH]c1=O)O2` |
| RR-OH | R,R | deacetyl | `O=c1[nH]c2c(c(=O)[nH]1)[C@@H]1C[C@](O)(O2)C(F)(F)c2ccccc21` |
| SS-OH | S,S | deacetyl | `O=c1[nH]c2c(c(=O)[nH]1)[C@H]1C[C@@](O)(O2)C(F)(F)c2ccccc21` |

Every enantiomer is tested **separately** throughout this protocol — never pooled, never run only
as the racemate. This is not a stylistic preference; it is load-bearing for the co-primary
chirality hypothesis in Section 5, and it is how the racemate's apparent affinity is understood
once (and if) real numbers exist: a racemate's apparent affinity is expected to approximate the
affinity of the stronger-binding single enantiomer, so racemate-only data would silently discard
the signal this protocol is built to detect.

The O-acetate is the as-synthesized form; the free alcohol (deacetyl, the "OH" species above) is
the form expected to be biologically active if the acetate is hydrolyzed intracellularly or in
media. Because that hydrolysis (esterase-mediated OAc→OH interconversion) can happen inside live
cells or in serum-containing media during the assay time-course, any claim that an observed effect
is attributable cleanly to one form rather than the other is gated on the stability check described
in Section 3 — it is not to be assumed from the nominal dosed species.

## 3. Mandatory pre-assay QC

Two QC steps are mandatory and precede any data collection or interpretation. Neither is optional,
and neither is a formality: skipping either turns any subsequent result into an artifact that
cannot be distinguished from a real one.

**Chiral-purity QC.** Each of the four dosed stocks (RR-OAc, SS-OAc, RR-OH, SS-OH) must be assayed
for enantiomeric excess (ee) by chiral HPLC or SFC before use, and the measured ee reported
alongside every downstream result for that stock. This is not a box-ticking exercise: the
co-primary chirality statistic in Section 5 is only interpretable on enantiopure material. A stock
contaminated with its mirror image will show artificially reduced R,R-vs-S,S discordance simply
because the "R,R" well also contains some S,S (and vice versa) — and a resulting 0-covering
(concordant) enantiopreference confidence interval from a **low-ee stock is uninformative, not
evidence for aggregation**. Concretely: if either RR or SS stock ee falls below a level the partner
lab and analysis team agree is adequate for the discordance test to have power, that arm's
enantiopreference result must be flagged as uninterpretable in the returned summary rather than
folded into the F5 (promiscuity) call.

**OAc→OH stability QC.** Before any result is attributed cleanly to either the acetate or the
deacetyl form, an LC-MS stability check must confirm the identity of at least one representative
species over the assay's full time-course, in the actual assay matrix (media or lysate, as
relevant to the assay in question). If esterase activity in the matrix converts a meaningful
fraction of dosed OAc to OH (or vice versa) within that window, any claim that isolates "OAc
biology" from "OH biology" is invalid, and the result must instead be reported as reflecting
whichever mixture of species was actually present — not the nominal dosed form.

## 4. Stage 1 — un-blinder (primary: Cell Painting)

Stage 1 is target-agnostic by design: because in-silico is exhausted (Section 2 of the design spec;
see `docs/results_boltz_cage.md`) and gave no confident target, the wet lab breadth-first
un-blinder — not a targeted assay — is the primary experiment.

**Primary readout: Cell Painting.** Run the deacetyl (OH) form, both enantiomers separately, plus
vehicle. Match the resulting morphological profile against a reference-compound library to infer
mechanism of action. Before any morphological signature is scored, a **cytotoxicity/viability
range-finder** determines a sub-lethal dose: an "effect" observed only at an overtly cytotoxic
concentration is a generic stress phenotype, not a specific mechanism-of-action signature, and must
be reported separately from a genuine hit, never merged with it.

A **hit** at Stage 1 requires all three of the following:
1. a reproducible morphological signature whose distance from vehicle exceeds τ_CP = 0.35
   (morphological distance to vehicle, assay-normalized) by the LCB rule (Section 8) — i.e.
   `effect − z·σ_assay ≥ τ_CP` with `z = 1.645`;
2. at a dose confirmed sub-lethal by the range-finder in (a) above;
3. survival of the aggregation counter-screen described immediately below.

**Modality-split aggregation counter-screen.** The barbiturate/cyclic-ureide head of the cage is a
known promiscuous H-bonder with colloidal-aggregation risk, so any apparent hit must be checked
against non-specific aggregation before it is trusted. Critically, the classic in-well detergent
spike used to rule out aggregation (Triton at sub-solubilizing concentration) is **not valid for a
live-cell primary readout**: it lyses live cells, which would destroy the Cell Painting signal
outright rather than distinguish specific engagement from aggregation. The counter-screen is
therefore split by modality:
- **For Cell Painting:** a **cell-free** counter-screen on the same stock — dynamic light
  scattering (DLS) and/or nephelometric turbidity, and/or a parallel biochemical enzyme-inhibition
  assay run with and without an in-well Triton spike — together with the Cell Painting panel's own
  viability channel as a built-in sanity check.
- **For DSF and reporter assays (Sections 4 fallback and 5):** the classic in-well detergent spike,
  at roughly 0.01–0.03% Triton, is appropriate because these are cell-free or tolerant-of-detergent
  formats.

The underlying `aggregation_guard` logic (Section 8) is the same regardless of modality; only the
source of the detergent-arm signal differs.

**Fallbacks (documented, not primary).** If Cell Painting access or reference-library matching is
unavailable, two fallbacks are pre-registered:
- **DSF/nanoDSF thermal-shift panel** of available purified proteins. A hit is ΔTm > 2 °C (τ_DSF =
  2.0, in degrees Celsius) by the LCB rule.
- **Affinity chemoproteomics pulldown**, gated on the partner lab's competence to synthesize a
  suitable affinity probe from the cage scaffold — this is explicitly a НИОХ-competence-dependent
  fallback, not a default expectation.

## 5. Co-primary — chirality specificity

Independent of whatever target (if any) Stage 1 or Stage 2 turns up, the R,R-vs-S,S enantiomer
comparison is run as a **co-primary** endpoint in every single stage of this protocol, not merely
as a secondary control. This is the first prospective test, on real wet-lab ground truth, of
Paper-1 invariants #5 and #6: that chirality is real, physically discriminable signal (a specific
binding pocket can produce reproducible R,R-vs-S,S discordance), and that a purely colloidal or
promiscuous-aggregation artifact should show no such discordance because aggregation is not
expected to be stereospecific.

The statistic is the confidence interval on the difference (R,R − S,S), computed per the
`enantiopreference` function in `src/screen/prospective.py`: given per-arm effect estimates `rr` and
`ss` with their own measured standard errors `σ_rr` and `σ_ss`,

```
delta = rr - ss
sigma_delta = sqrt(sigma_rr**2 + sigma_ss**2)
CI = delta ± z * sigma_delta          (z = 1.645)
```

If this CI **excludes** zero, the result is **specificity-consistent** (real discordance between
enantiomers, unlikely for a non-stereospecific aggregator). If the CI **covers** zero, the result is
**aggregation-consistent** (no detectable discordance). This CI is reported at every stage
regardless of whatever the primary target-engagement verdict is — it is co-primary, not
conditional, and (per Section 3) is only interpretable when both stocks in the comparison have
adequate chiral purity.

To be explicit about how this statistic interacts with escalation decisions: the enantiopreference
CI is a **specificity qualifier**, not itself a trigger for escalation to Stage 3. A discordant
result does not by itself commit resources to further work; it only distinguishes, among results
that already pass the LCB commit rule, which look more consistent with specific engagement versus
which look more consistent with promiscuous binding.

## 6. Stage 2 — cheap bound: GR→AR reporter (both modes)

The one mechanistically grounded in-silico hypothesis worth a dedicated bounding experiment is the
steroid-receptor signal from docking (GR/AR ranked highest), even though Boltz-2 co-folding
independently pushed back on it (gap vs. agonist anchor: GR −0.524, AR −0.402, ER −0.700 — see
`docs/results_boltz_cage.md`). Because the two in-silico methods disagree, Stage 2 treats this as a
cheap bound to check, not the primary search.

Run GR then AR transactivation/reporter assays, all four species (both enantiomers, both forms),
with dexamethasone (GR) or DHT (AR) as positive control and vehicle as negative control. Critically,
run each receptor in **both agonist mode** (cage alone → does it activate the receptor?) **and
antagonist mode** (cage plus a sub-maximal dose of dexamethasone or DHT → does it suppress the
agonist response?). This dual-mode requirement exists because Boltz's binding_confidence score does
not distinguish an agonist-competent pose from an antagonist-competent one — a compound could
plausibly occupy the pocket and block the natural agonist without itself activating the receptor,
and the in-silico record has no way to rule that out.

Rank the priority order GR ≥ AR > ER, consistent with the in-silico record (GR and AR are the
targets docking and Boltz both engaged with, even in disagreement; ER was the target Boltz rejected
most strongly, so it is lowest priority and not part of the primary reporter panel).

"Steroid closed" (i.e., this hypothesis is set aside, consistent with the Boltz-2 negative) requires
**both** the agonist-mode **and** the antagonist-mode result to be negative by the LCB rule
(`effect − z·σ_assay < τ_reporter`, with `τ_reporter = 0.30`, fraction of agonist-max response, and
`z = 1.645`) — a negative in only one mode is not sufficient to close the hypothesis, precisely
because the two modes probe different pharmacology.

## 7. Stage 3 — escalation of confirmed hits only

Only hits that survive the full decision logic in Section 8 (LCB commit, post-multiplicity-correction,
aggregation-guard survival) proceed to Stage 3. Escalation follows the standard confirmatory ladder:
surface-plasmon resonance or isothermal titration calorimetry (Kd) first, then a functional assay
appropriate to the confirmed target, then co-crystallography if the earlier steps hold up. At every
step of this ladder, the R,R-vs-S,S enantiopreference comparison (Section 5) is re-run and
re-confirmed — a target confirmed at the reporter stage that loses its enantiodiscordance under
tighter biophysical measurement is a meaningful downgrade in confidence, not a detail to drop.

## 8. Decision logic

All decisions in this protocol are made mechanically from measured quantities using four functions
implemented in `src/screen/prospective.py`, with every constant taken verbatim from
`data/cage/prospective_prereg.yaml`. None of these thresholds are adjustable after the prereg
commit.

**Commit rule (LCB).** A result commits to escalation only if its lower confidence bound clears the
relevant threshold: `LCB = effect − z·σ_assay ≥ τ`, where `z = 1.645` throughout, and τ is
modality-specific: `τ_CP = 0.35` (Cell Painting, morphological distance to vehicle, assay-normalized),
`τ_DSF = 2.0` (thermal shift, degrees Celsius), `τ_reporter = 0.30` (nuclear-receptor reporter,
fraction of agonist-max response). `σ_assay` is always a **measured** replicate or Hill-fit standard
error returned by the partner lab — never an ensemble-model or otherwise learned uncertainty. This
distinction matters: Paper 1's own results (Fig A) show that learned/ensemble variance heads can be
substantially overconfident relative to measured variance, and this protocol exists specifically to
avoid repeating that mistake one level up, in a wet-lab setting where the consequences of an
overconfident false commit are real resources spent.

**Multiplicity correction (Holm–Bonferroni).** Because Stage 1 (the un-blinder) and Stage 2 (GR and
AR, each in agonist and antagonist mode) together generate several independent LCB decisions per
active enantiomer form, no single LCB pass is permitted to trigger Stage-3 escalation on its own
without correcting for the resulting multiple-comparisons problem. The escalation-triggering family
is `K = {un-blinder} ∪ ({GR, AR} × {agonist, antagonist})`, evaluated per active enantiomer form.
Holm–Bonferroni is applied across this family before any member of it is allowed to trigger
escalation. The enantiopreference CI (Section 5) sits outside this family: it is a specificity
qualifier reported alongside each comparison, not itself a member of the escalation-triggering set,
and so is not itself subject to this correction as a trigger.

**Stop rule (calibrated stopping, measured σ only).** The decision to stop collecting further
replicates and commit to the current best call is governed by `stop_rule`, which computes an
internally calibrated confidence `conf = Φ((effect − τ)/σ_assay)` — monotone increasing in
`(effect − τ)/σ_assay` — and stops once `conf ≥ bound`, with `bound = 0.90` fixed in the prereg.
`stop_rule` structurally takes only a measured `σ_assay`, never a pre-computed or model-derived
confidence, so this stopping decision cannot be gamed by substituting an optimistic learned
uncertainty for the real measured spread across replicates.

**Aggregation guard (σ-aware).** A signal is judged to **survive** the aggregation counter-screen
(and thus not be a colloidal/systematic artifact) only if its detergent-arm signal, discounted by
its own measured uncertainty, still clears a fraction of the no-detergent signal:
`survives ⟺ (detergent_signal − z·σ_detergent) ≥ frac · signal`, with `frac = 0.5` (a signal
retaining at least half its magnitude under detergent challenge, at 95% one-sided confidence,
counts as surviving) and `z = 1.645`. A signal that fails this test is flagged as a systematic
(aggregation) artifact regardless of how strongly it passed the LCB commit rule on its own.

**Minimum replicates.** All arms are run with a minimum of 3 technical replicates
(`min_replicates = 3` in the prereg), derived from a stated power calculation: with
`σ_assay = σ_rep/√n`, passing the LCB with margin `m = effect − τ` at `z = 1.645` requires
`n ≥ (z·σ_rep/m)²`. For the DSF anchor (`τ = 2 °C`, an assumed engaged-target shift of ~3 °C giving
`m = 1 °C`, and an assumed `σ_rep ≈ 1 °C`), this gives `n ≥ (1.645·1/1)² ≈ 2.7`, rounded up to 3 — a
replicate floor, not a guarantee of adequate power for every assay modality; modalities with
noisier per-replicate variance than the DSF anchor may need more than the floor to actually reach
the intended detection probability, and the partner lab should flag if measured `σ_rep` for a given
assay materially exceeds this assumption.

## 9. Exhaustive disposition table

Every wet-lab outcome that can arise from crossing (LCB pass/fail, post-Holm) × (enantiodiscordant
yes/no) × (survives the aggregation counter-screen yes/no) has a pre-registered disposition, so that
no real result is left without an assigned interpretation decided in advance of seeing it:

| LCB (post-Holm) | enantiodiscordant | survives detergent | disposition |
|---|---|---|---|
| pass | yes | yes | **HIT** — specific engagement of a target/pathway; escalate to Stage 3. Supersedes the F2 in-silico steroid-hypothesis framing if the hit is not steroid-family. |
| pass | no | yes | **ambiguous-engagement** — real engagement by the LCB rule, but fails the co-primary specificity bar; reported as-is, does **not** escalate. |
| pass | any | no | **aggregation-artifact** — LCB-significant but abolished by detergent challenge; consistent with F5, not treated as a hit. |
| fail | yes | yes | **inconclusive-but-suggestive** — held for one pre-registered confirmatory replicate; never silently folded into a "promiscuous binder" conclusion. |
| fail | no | any | **F5-confirmed-promiscuous** — a weak, promiscuous binder; a real, publishable negative that bounds the chemotype. |
| fail | yes | no | **aggregation-artifact** — discordance without detergent survival is a systematic-artifact flag, not evidence of specificity. |

Reporter-specific dispositions layer on top of this table for Stage 2:
- **Steroid rescue:** a GR or AR reporter positive (dose-responsive, enantioselective,
  agonist-comparable, detergent-surviving, LCB-passing post-Holm) in **either** agonist or
  antagonist mode rescues the steroid hypothesis despite the Boltz-2 negative.
- **Steroid closed:** requires **both** agonist-mode and antagonist-mode negatives by LCB (Section
  6) — consistent with the Boltz-2 read.
- **Non-enantioselective reporter positive:** a dose-responsive, agonist-comparable,
  detergent-surviving but non-enantioselective reporter result is **ambiguous-engagement**, as in
  the main table above — not a steroid rescue.

**Explicit non-kills.** A negative result across the board (the cage reads as an F5-confirmed
promiscuous binder, or Stage 2 reads steroid-closed) is **not** a failure of this protocol or of the
underlying research program — it is the false-lead guard working exactly as designed, and a
legitimate, publishable negative that bounds the chemotype. Similarly, any CPU-nondeterministic
drift observed between repeated in-silico scoring runs (noted in the design spec as a known
property of the Boltz pipeline) is a reproducibility note about the computational method, not a
kill condition for this wet-lab protocol.

## 10. Reporting

Results are appended to **this file**, in a new section below this line, only **after** wet data has
been collected and scored — the frozen pre-registration in `data/cage/prospective_prereg.yaml` is
never edited to reflect results, and no number in Sections 1–9 above is altered post-hoc.

The partner lab does not need to, and should not, transmit raw images, spectra, or plate-level raw
data back for this protocol's decision layer to function. What is needed, per (species, arm)
combination, is:
- the effect estimate and its **measured** standard error (σ_assay, from replicate or Hill-fit
  variance — never a model-derived or subjective confidence);
- the corresponding aggregation counter-screen result (the detergent-arm signal and its measured σ,
  in the modality-appropriate form described in Section 4);
- the chiral-purity (ee) QC result for the stock(s) used in that comparison;
- for Stage 1, the viability/range-finder outcome establishing that the scored dose was sub-lethal;
- for any OAc-vs-OH form attribution, the LC-MS stability QC result.

These summary statistics are exactly the inputs that `src/screen/prospective.py` consumes:
`decision_lcb`, `enantiopreference`, `aggregation_guard`, and `stop_rule` each take measured
effect/σ pairs, and `score_forecast` takes a dictionary of these observations to produce the final
`FalsificationReport` scoring the F1/F2/F4/F5 forecasts from `data/cage/prospective_prereg.yaml`
(F3 is explicitly not a scored forecast — it is recorded context, not a falsifiable claim). The
appended results section should report that scorecard verbatim, alongside the disposition (Section
9) that the observed (LCB, enantiodiscordance, survival) triple maps to.

---

## Results

*(To be appended after wet-lab data is collected and `score_forecast` has been run. Nothing below
this line exists yet; this section is a placeholder marking where results belong, not a claim that
any exist.)*
