# Phase-0 plan — the target-finding arm (the cage merge), before any build

*Mirrors `docs/ANALYSIS.md`: think + design + define the falsifiable gate BEFORE building.
This is the SECOND arm. The first arm (RBFE calibration: Figs A–E, the four theorems,
gauge AL) is validated. This arm is NOT — the ligand-only version (Fig F) failed, the
structure version is unbuilt. "Merging the screen" honestly = build + validate this arm.*

## 0. Two things that must not be glued together
- **Operational debt (NIOCH):** run the architecture over the cage → a ranked list of
  target hypotheses with calibrated uncertainty + a VoI-prioritised experiment plan.
  This is a **report, not a paper claim**; deliver it regardless of where publication
  lands; it need not be publishable. **Blocked on one input: the cage structure** (SMILES
  or SDF) — not in the repo.
- **Publication claim:** cell screening is **not prospectively verifiable** (target
  unknown, no assays). It must **not** be a load-bearing result, or the MBAR-reviewer
  problem returns in target-ID form. It enters the paper only as **motivation** (intro)
  and as a **deployment case study** (generated campaign, "pending assays"), never as
  "we found the target".

## 1. The reframe — the cage justifies the architecture
Standard ligand-based target-ID does not work on this molecule: ~0 ChEMBL
similarity hits, and Fig F already showed ligand similarity saturates at the shape
baseline. The molecule therefore **dictates** a structure-based, calibrated,
chirality-aware, VoI-driven method: no ligand shortcut, experiments are expensive, start
from scratch, two stereocentres. "Forced by the problem" at paper level — the cage
motivates the architecture in §1, it does not prop it up in results.

**Architecture flag (carry into calibration):** the **barbiturate ring (NH–CO–NH–CO)**
will dominate the complementarity signal, and barbiturates are promiscuous (GABA-A and
beyond). Without calibration + retrospective validation the method will "find" the
obvious-but-promiscuous pharmacophore and fool us. Calibration here is **protection
against a false lead on an orphan scaffold** — exactly why the machine exists.

## 2. The falsifiable gate (decides one paper vs two)
> **Does retrospective, target-disjoint, structure-based target-ID beat the baselines
> (the ligand-only Fig F shape baseline, and a docking-score-only ranker)?**

- **Beats →** one paper: method (Figs A–G) + cage-as-motivation + a validated
  structure-based target-ID figure + deployment case study. Stronger, more compelling.
- **Does not →** two papers: ship the methods+theory spine (current Figs A–G) now; the
  cage campaign becomes Paper 2, maturing with data. **Do not rescue the merge if the
  second arm does not stand.**

## 3. Three architectural questions — proposed answers (RATIFY before building)

### Q1 — an honest split (target-disjoint, not ligand-disjoint)
Holding out ligands but keeping their targets lets the model recover by ligand-similarity
(the Fig F trap). Holding out whole *targets* is necessary but not sufficient: if a held-out
target shares a **fold/family** with a training target, recovery is a pocket-similarity
shortcut. So: **fold/family-disjoint at the pocket level** — cluster the pocket library by
structural similarity (CATH/SCOP family or pocket-shape clustering) and hold out whole
clusters. Add **hard decoy pockets** (druggable apo/AF/unrelated holo sites) so high
recovery is non-trivial. The ligand-only Fig F method must *not* improve on this split
(it has no structure) — that is the foil to beat.

### Q2 — calibrating recovery under extreme class imbalance
Per query the positive is 1 true target among N decoys → per-pocket "P(binds)" is
degenerate under imbalance. **Calibrate the decision-relevant quantity instead:** the
probability that the true target is in the top-k of the ranking (recovery calibration),
mirroring Fig F's P(binds) ECE and Fig G's MC top-k confidence — imbalance-robust by
construction (it is about the *ranking*, not a per-pocket binary). The raw
score→probability map (docking/ABFE score → P) is fit on a **target-disjoint calibration
fold** with isotonic/beta calibration against **hard negatives** (ligand vs decoy pockets
of matched druggability), and reported as recovery-ECE under the deployment imbalance.

### Q3 — VoI: which target to test first, over docking → ABFE → assay
This is the decision layer, and it is the **gauge-aware cost-aware KG of `active.py`**
lifted to the target-ranking decision over a **multi-fidelity ladder** (MFBind-style):
- belief = calibrated distribution over which pocket is the true target;
- actions = run tier-t on candidate i: docking (cheap, noisy) → Boltz-ABFE (medium) →
  wet assay (expensive, near-definitive), each with cost `c_t` and information yield `η_t`;
- **VoI = expected reduction in decision loss (P(commit to the wrong target)) per cost.**
Calibration (Q2) is what makes the VoI honest — Fig G's lesson: a miscalibrated belief
mis-estimates the information yield and mis-orders the experiments.

## 4. Bounded build sequence (NOT the full trunk+heads+solver — anti-scope-creep)
1. **Curate the retrospective benchmark** (public): ligand→target pairs with holo
   structures (PDBbind / BindingDB+PDB subset) + decoy pockets (apo/AF/unrelated);
   fold-disjoint split (Q1). *(cheap, ratification-independent — can start once design ok)*
2. **Score** (ligand × candidate-pocket) with ONE orchestrated structure tool —
   **decision needed**: Boltz-ABFE (paid MCP; accurate; per-job cost) vs install
   docking (vina/smina; free; noisier). Do not rebuild — orchestrate.
3. **Calibrate + rank** (Q2) → target-disjoint top-k recovery + recovery-ECE vs the Fig F
   ligand baseline + a docking-only baseline. **← THE GATE.**
4. **If gate passes:** VoI decision layer (Q3, reuse `active.py`) → "which target first" plan.
5. **Deployment:** run the cage through the validated pipeline → ranked target hypotheses
   + VoI experiment plan ("pending assays"). *(needs the cage structure)*

## 5. What I need to proceed (decisions / inputs)
1. **Cage structure** (SMILES/SDF) — hard blocker for steps 5 and the NIOCH debt.
2. **Structure tool for step 2:** Boltz-ABFE (paid, per-job confirmation + cost) vs
   install docking (free, noisier). Cost/quality call.
3. **Ratify §3** (split / calibration / VoI) or correct it.
4. **Scope confirmation:** target one paper (method + deployment), two-paper fallback,
   gated on §2. Build the second arm only if it stands.

The screening run and the validation harness are mine to build (per the brief); the four
above are the user's calls. I do not build until §3 is ratified and the tool/cost gate is
resolved — building before that is the infrastructure-expansion trap.
