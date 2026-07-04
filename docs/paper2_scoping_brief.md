# Paper 2 — target-finding, scoping brief (prepared 2026-07-04)

## 1. The question (one falsifiable sentence)

In the **orphan regime** — where a query ligand has near-zero maximum Tanimoto similarity to every reference active of every candidate pocket, so ligand-similarity is provably uninformative by construction — does a **calibrated structure/physics score** (relative-scored Boltz-2 co-folding, with a genuine ABFE sandwich variance on the shortlist) recover the true target **above both** the ligand-shape null and a raw-structure-score null, with trustworthy calibration?

## 2. Why attempt #1 failed (the honest Fig H post-mortem — not to be softened)

Attempt #1 ran a structure-based reverse-docking gate and it **failed cleanly**. Do not rescue this.

- **Docking lost to ligand-shape outright.** Top-1 recovery: raw smina **0.16** vs ligand-shape **0.95**; AUROC **0.62** vs **0.99** (`docs/results_figH.md`). Per-pocket z-normalisation made docking *worse* (top-1 0.08), so this was not merely an uncalibrated score — raw cross-pocket docking carries little genuine target-discriminating signal. It picks the **greasiest pocket**: AChE recovered 4/5 of its own actives and attracted everyone else, while EGFR/CDK2/CA-II/FXa got 0/5.
- **The benchmark never reached the orphan regime.** Median test-ligand→train Tanimoto was **0.64**; only ~5% of queries fell below 0.4. On such a benchmark the ligand-shape null already wins, so a structure "win" would be **unattributable** — it could be residual similarity leakage. The gate as run *could not have passed for a real reason.* This conflates two separable failures: (i) score-comparability (raw docking is uncomparable across pockets) and (ii) benchmark-regime (the null was never disarmed).
- **The gate was mis-specified.** It was an *absolute* contest (structure vs shape on the same easy queries) instead of an *attribution* contest (structure vs shape **in the stratum where shape is provably near-random**). The single validity test — does the shape null collapse to random in the orphan stratum? — was never applied.
- **The steroid hypothesis was rejected by Boltz-2.** On the NIOCH cage, docking proposed steroid receptors (GR 80–89% / AR 75–89% percentile-among-actives), but Boltz-2 co-folding **rejected** it: cage best-form binding_confidence GR 0.45 / AR 0.54 / ER 0.28 vs agonists 0.98/0.94/0.98 (gaps −0.52/−0.40/−0.70; `docs/results_boltz_cage.md`). Boltz **also** rejected the AChE greasy-pocket artefact docking had falsely scored ~100% (cage 0.27 vs donepezil 0.58) and passed its GR/AR/ER positive controls — which is what makes its steroid negative *credible*. Net in-silico read: **no confident target; a weak, promiscuous binder.**

The one genuinely useful signal from attempt #1: **Boltz-2, read as a gap-vs-anchor, discriminated where raw docking could not.** That, not docking, is the seed for attempt #2.

## 3. Design ingredients

### (a) Calibrated Boltz-ABFE binding estimate with UQ (reuse Paper-1 machinery)
Transplant the Fig-B decomposition (epistemic + aleatoric + conformal) from BAR edges onto co-folded complexes; `src/bar/calib.py` is feature-dimension-agnostic and reused near-verbatim.
- **Features** per (ligand, candidate-pocket): the Boltz vector `[binding_confidence, optimization_score, structure_confidence, iptm, ptm, complex_plddt, complex_iplddt]` + pocket-restricted PAE summary stats (mean/max interface PAE from `pae.npz` over pocket residues) + the **chirality-odd 0o pseudoscalar** from `src/bar/trunk.py:chir_pseudoscalar` (invariants #5/#6: ECFP is enantiomer-blind; the signed-volume channel is the sole chirality carrier, and the cage is chiral with 2 stereocentres).
- **Stage 1 — epistemic σ:** bootstrap deep ensemble (`calib.py:train_ensemble`, n_members=8) regressing features → experimental pIC50/pKd on a labelled fold; **member spread inflates on orphan/OOD complexes** — this *is* the regime discriminator.
- **Stage 2 — aleatoric σ:** run-to-run Boltz score noise from seed/MSA-subsampling replicates on a subset. **Honest label:** this is the *analogue* of the sandwich decomposition, **not** the literal `B/I²` (which requires alchemical work samples). Do not call it "the sandwich" — that belongs to the ABFE tier only. Mislabeling it would recreate the reviewer-MBAR problem one tier up.
- **Stage 3 — guarantee:** normalized/Mondrian split-conformal (`calib.py:conformal_q`, σ_total = √(σ_epi² + σ_ale²)) fit on a **target-disjoint** calibration fold with hard decoy pockets → coverage-guaranteed interval + a calibrated P(binds) via isotonic score→probability.

### (b) Target-disjoint, low-similarity ORPHAN benchmark with a ligand-similarity null
The load-bearing deliverable. Must be **stratified**, **fold-disjoint**, and **relative-scored**.
- **Pocket library:** ~30–50 holo pockets spanning diverse folds — PDBbind refined + LIT-PCBA (15 targets, unbiased actives) + DUD-E/DEKOIS for property-matched decoy pockets and hard negatives. **Cluster and hold out whole CATH/SCOP (or Pfam) fold-clusters**, not ligands, so a held-out true target never shares a fold with any training pocket (kills the pocket-similarity shortcut).
- **Queries:** each has a single known true target (drop promiscuous ligands, per Fig F). For each, compute `s = max Tanimoto (Morgan r=2)` to *any* training active of *any* pocket. **Stratify**: high (s≥0.5), mid (0.35–0.5), **ORPHAN (s<0.35, ideally <0.2)** — the orphan stratum decides the gate. Enrich it with LIT-PCBA/BindingDB novel-chemotype ligand→target pairs.
- **Scorer — relative, never raw cross-pocket:** score each (ligand, pocket) as a **percentile/gap against that pocket's own reference actives** docked identically (the structural analogue of calibrated `B/I²` vs raw `1/I`). Tier it: **(t1)** free smina percentile-among-actives; **(t2)** Boltz-2 gap-vs-anchor binding_confidence. Report the gate on the best tier; *if t1 ties the null but t2 beats it, that is the honest finding (raw docking dead, cofolding alive).*
- **Mandatory nulls, measured per stratum:** (i) ligand-shape (max Tanimoto to a pocket's actives); (ii) random (1/N_pockets); (iii) **pocket-size/hydrophobicity null** (rank by pocket SASA/lipophilicity) — to prove structure is not re-learning pocket volume (the Fig H disease).

### (c) The wet-lab un-blinder loop (real ground truth)
Because every targeted in-silico lead is non-robust, the wet-lab entry point must be **target-agnostic**, run **breadth-first**, on the deacetyl (OH) form, **both enantiomers**, with a **mandatory detergent counter-screen** (the barbiturate/cyclic-ureide head is a flagged promiscuous H-bonder + colloidal-aggregation risk).
- **Stage 1 (un-blinder, run first):** cheapest-available of — DSF/nanoDSF thermal-shift panel of purified proteins (cheapest engagement read); Cell Painting morphological profiling vs a reference-compound library (best orphan-MoA read, one plate); chemoproteomics pulldown (most direct, highest cost, needs a synthesizable affinity probe).
- **Stage 2 (cheap bound, not the lead):** GR then AR reporter/transactivation, both enantiomers, both forms, dexamethasone/DHT positive + vehicle negative controls. A clean **negative bounds** the one mechanistically-grounded hypothesis; a positive would **rescue** the steroid hypothesis despite Boltz. Rank GR ≥ AR > ER (only GR is cross-screen-reproducible; ER is the weakest, not worth a first slot).
- **The enantiomer pair is the built-in specificity control** (exercises Theorem 4 / invariants #5–6 on real ground truth): a promiscuous aggregator should show **no** reproducible R,R-vs-S,S preference; true pocket binding can. It is also the in-silico tie-breaker docking and Boltz disagreed on.
- **Decision logic (ported from Fig I/G/L, NOT the surrogate σ):** commit a target to expensive escalation only if `LCB = effect − z·σ_assay ≥ τ`, where **σ_assay is measured replicate/Hill-fit uncertainty**, τ a pre-registered biologically-meaningful effect (e.g. >2 °C Tm shift). Stop rule mirrors Fig G; the detergent/enantiomer-discordance check is the Fig-L systematic-vs-sampling detector. The surrogate donates its *decision calculus*, never a validated per-target probability.

### (d) The pocket library, if needed
The (b) pocket library doubles as the deployment library for the cage. For orphan queries there is **no known binder to seed the pocket**, so pocket detection is itself uncertain — a confound the calibration must **not** launder into false confidence. Require a **minimum reference-active count per pocket** (reversedock flagged thin 100% rows as artefacts); exclude or widen-CI thin pockets.

## 4. The falsifiable gate + honest KILL conditions

**Pre-register before any scoring. The decision is made in the ORPHAN stratum (s<0.35) only.**

**STEP 0 — validity gate (must pass first, else result void):** the ligand-shape null must **collapse to near-random** in the orphan stratum — top-1 within CI of 1/N_pockets (e.g. ≤0.10 for N=30 vs random 0.033), shape-null AUROC ≤ 0.60. If the shape null still beats random, the stratum is **not truly orphan** → re-stratify (tighten s). Never report a structure win against a still-informative null. *(This is the single test Fig H never applied.)*

**STEP 1 — the gate:** with the shape null neutralised, the best structure tier (relative-scored) must beat **both** nulls by a pre-registered, bootstrap-significant margin:
- orphan top-1 recovery ≥ 2× random **AND** ≥ shape-null + 0.15 absolute, with 5-seed bootstrap 95% CI on (structure − shape) **strictly > 0**;
- orphan recovery-AUROC ≥ 0.70 (vs shape ~0.50–0.60);
- **must also beat the pocket-size null** by the same CI test (proves complementarity, not greasiness).
- Powered: ≥30–50 orphan queries across ≥8 fold-disjoint pocket clusters, ≥5 seeds.
- Sharpest single number: **calibrated recovery AUROC − ligand-shape recovery AUROC > 0, CI excluding 0.**

**STEP 2 — calibration (secondary):** orphan recovery-ECE ≤ 0.05 under deployment imbalance, so the top-k confidence that drives the VoI/escalation layer is trustworthy.

**SANITY CHECK:** in the **high-similarity** stratum the shape null should still win/tie (reproducing Fig H) — confirming the orphan effect is a genuine stratum interaction, not a scorer artefact.

**HONEST KILL CONDITIONS — the arm is DEAD (cage reverts to a "pending assays" report, no methods claim) if ANY of:**
1. **Primary kill — structure ties the null in the orphan regime.** With Step-0 validity confirmed, the best structure tier's top-1 CI overlaps **both** the shape null and random. If structure cannot separate targets *where similarity is provably uninformative*, there is no attributable structure signal and the premise is refuted (Fig H repeating one fidelity tier up).
2. **Validity kill — no honest orphan stratum exists.** After tightening s, the shape null never collapses to random while retaining ≥50 orphan queries across ≥8 folds → the question is **retrospectively undecidable on public data**; frame only as a prospective NIOCH case study, never a benchmarked claim.
3. **Greasy-pocket kill.** Structure beats shape but **not** the pocket-size/hydrophobicity null → the signal is re-learned pocket volume, not recognition (the orphan cage would be misassigned to the greasiest pocket).
4. **Calibration kill.** Structure wins recovery but orphan recovery-ECE > 0.10–0.15, or conformal coverage misses target after target-disjoint calibration → the UQ is untrustworthy and the "calibrated decision" thesis collapses (downgrade to ranking-only case study).
5. **Unlearnable-score kill.** Boltz binding_confidence shows no monotone, target-disjoint relationship to experimental affinity on a labelled fold → nothing to calibrate.
6. **No-distance-signal kill.** Epistemic ensemble σ does **not** rise on orphan/OOD complexes → the method cannot know when it is extrapolating.
7. **Not-load-bearing kill.** A cheaper baseline you did not need Boltz for (pocket-shape descriptor + isotonic) matches the Boltz ranker → co-folding is not load-bearing.

**Explicit NON-kills (do not over-claim):** the cage reading **negative** is not a kill — that is the false-lead guard working. Run-to-run point-estimate drift (torch CPU non-determinism, a Paper-1 caveat) is a reproducibility note, not a kill.

## 5. Sequencing & cost (cheapest-decisive-experiment first)

Sequence it **Fig-H-style: curate first (cheap to fail), then gate on the free tier, escalate to paid only if warranted.**

1. **Curate the orphan benchmark FIRST** — fold-cluster the pocket library, stratify queries by similarity, build the orphan stratum from LIT-PCBA/BindingDB novel chemotypes, run the **leakage audit** (held-out vs train pocket TM-score + site-motif overlap; drop near-duplicates). This is the true blocker — if <~30 genuine orphan queries across ≥8 folds can't be assembled, hit **Validity Kill (2)** now and stop. *~2–4 weeks data engineering, in-silico, near-zero $.* Reuses `src/screen/recovery.py` metrics unchanged.
2. **Free smina-tier orphan gate + Step-0 validity check** — go/no-go on zero paid compute. Reuses `src/screen/dock.py` as-is (hours–days). If smina passes, great; if smina ties the null **while Step-0 validity holds**, escalate.
3. **Paid Boltz-2 tier — orphan stratum only.** Gate-scale ~40 queries × ~30 pockets × ~3 seed replicates × $0.025 ≈ **$90–150**, tens of min/job (rate-limit/wall-clock bound). Always `estimate-cost` / `start=false` first. Score relative to each pocket's own reference actives (percentile-among-actives is the main cost driver — budget it).
4. **Calibration wrapper + gate** — Boltz-feature extractor (parse `index.jsonl` + `pae.npz` pocket summaries) + calibrated-probability wrapper; reuse `calib.py` (ensemble + conformal + isotonic) and `recovery.py`. *~1 week.*
5. **STRETCH — multi-fidelity ABFE ladder.** Only if the gate passes: escalate the shortlist to a real alchemical ABFE (OpenFE/OpenMM) that returns the genuine sandwich variance `B/I²` (invariant #1, the ACTUAL BAR bottleneck re-entering as a first-class citizen), governed by the existing gauge-aware cost-aware KG (`src/bar/active.py`). *Weeks–months, GPU-hours, shortlist only.*
6. **Wet-lab un-blinder (§3c)** — dry design/protocol port is LOW (~2–4 days, reuses `active.py`, `data/cage/`, `docs/cage_assay_request.md`). Wet execution is MEDIUM and gated on an external partner (Stage 1 one Cell Painting plate or DSF panel; Stage 2 commercial GR/AR reporter kits ~1–2 weeks; Stage 3 SPR/ITC/co-crystal on confirmed hits only). Deliver as a **pre-registered protocol** so any partner lab can execute decisively. Compute for the in-silico half is already sunk.

**Do not build the full trunk+heads+solver up front — that is the project's stated failure mode.** Curate → free gate → paid gate → ladder, killing cheaply at each step.

## 6. First step

**Implementation begins with a `superpowers:brainstorming` session + a written plan — NOT code.** This brief is the *input* to that brainstorm, not a build order. The brainstorm must first stress-test the load-bearing risk (can a genuine orphan stratum of ≥30 queries across ≥8 fold-disjoint clusters even be assembled from public data?), because that single question decides whether Paper 2 is testable at all — and if it cannot, the honest outcome is **Validity Kill (2)**: the cage stays a pre-registered prospective case study and there is no benchmarked methods claim. Only after the brainstorm converges on a written, pre-registered plan (gate thresholds, strata, nulls, leakage audit fixed *before* any scoring) does any code get written.