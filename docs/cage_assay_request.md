# Assay request — difluoronaphthalenone "cage" fragment (target ID)

**Status:** in-silico hypothesis-generation complete; this is a request for the decisive wet-lab
experiments. Not a validated target claim. Companion detail: `docs/nioch_cage_report.md`
(docking), `docs/results_boltz_cage.md` (Boltz-2 cross-check).

## The compound
A rigid fused difluoronaphthalenone + 1,3-dicarbonyl (barbituric-acid) Michael adduct,
O-acetylated. Two ring-fusion stereocentres → **one diastereomer + its mirror** (test as two
single enantiomers). The O-acetate is the as-synthesised form; the free alcohol (deacetyl) is the
likely active/hydrolysed species. **Test all four species, enantiomers separately** (the molecule
is chiral; a chirality preference is real signal, and racemate apparent affinity ≈ the stronger
binder).

| species | enantiomer | form | MW / formula | SMILES |
|---|---|---|---|---|
| RR-OAc | R,R (given) | acetate | 350.3 / C16H12F2N2O5 | `CC(=O)O[C@]12C[C@H](c3ccccc3C1(F)F)c1c([nH]c(=O)[nH]c1=O)O2` |
| SS-OAc | S,S (mirror) | acetate | 350.3 / C16H12F2N2O5 | `CC(=O)O[C@@]12C[C@@H](c3ccccc3C1(F)F)c1c([nH]c(=O)[nH]c1=O)O2` |
| RR-OH | R,R | deacetyl (active) | 308.2 / C14H10F2N2O4 | `O=c1[nH]c2c(c(=O)[nH]1)[C@@H]1C[C@](O)(O2)C(F)(F)c2ccccc21` |
| SS-OH | S,S | deacetyl (active) | 308.2 / C14H10F2N2O4 | `O=c1[nH]c2c(c(=O)[nH]1)[C@H]1C[C@@](O)(O2)C(F)(F)c2ccccc21` |

**Flags:** orphan scaffold (no close ChEMBL analogue); the cyclic-ureide head is a strong,
**promiscuous H-bonder** — guard against non-specific / aggregation artefacts.

## The question
What protein does it bind? (de-novo target ID for an orphan fragment.)

## What the computation says (honestly)
No confident target. Docking raised one coherent hypothesis — **steroid receptors (GR, AR)** —
but a stronger orthogonal model (**Boltz-2 co-folding**) did **not** corroborate it: the cage
scores far below the known agonists in their own pockets (best-form binding_confidence vs agonist:
GR 0.45 vs 0.98, AR 0.54 vs 0.94, ER 0.28 vs 0.98), and Boltz correctly rejected a docking
greasy-pocket artefact (AChE), so its negative read is credible. The two methods disagree → the
steroid signal is **not robust**, and the cage looks like a weak, **non-specific** binder
in-silico. In-silico is exhausted; the wet lab is the un-blinder.

## Recommended experiments (value-of-information order)

**Tier 1 — the un-blinder (do this first): broad, target-agnostic profiling.**
Because the target is unknown and in-silico is non-robust, breadth beats depth. Any one of:
- **Cellular thermal-shift / CETSA or a thermal-shift (DSF/nanoDSF) panel** of available purified
  proteins — cheap, detects engagement by Tm shift.
- **Phenotypic profiling (Cell Painting)** — morphological fingerprint; match to reference
  compounds to infer mechanism, ideal for an orphan.
- **Affinity-based chemoproteomics** (if an affinity probe can be made from the scaffold) —
  pulls down and identifies binding proteins directly; the most direct un-blinder.
Test the deacetyl (active) form, both enantiomers; include a vehicle and, given the promiscuity
flag, a **detergent / counter-screen control** to exclude colloidal aggregation.

**Tier 2 — a cheap *bound*, not the lead: targeted nuclear-receptor reporter.**
GR and AR transactivation/reporter assays (commercial kits), **both enantiomers, both forms**
(acetate as-given + deacetyl; the acetate scored marginally higher in Boltz). Include the known
agonist (dexamethasone for GR, DHT for AR) as positive control and a vehicle negative. A clean
negative *bounds* the orphan claim; the two in-silico methods no longer agree it will be positive,
so this is confirmatory/bounding, not the primary search.

**Tier 3 — escalate confirmed hits only:** binding (SPR/ITC, Kd) → functional → co-crystal, and
confirm the enantiomer preference at each step.

## Decision criteria
- **Tier 1 hit** (a reproducible engaged protein / a clean MoA fingerprint) → that protein/pathway
  becomes the lead; escalate (Tier 3). This supersedes the steroid hypothesis.
- **GR/AR reporter positive** (dose-responsive, enantioselective, agonist-comparable, survives the
  detergent counter-screen) → the steroid hypothesis is rescued despite Boltz; escalate.
- **GR/AR reporter negative** → consistent with Boltz; the steroid hypothesis is closed; rely on
  Tier 1.
- **Everything negative / non-specific only** → the orphan scaffold is a weak promiscuous binder;
  report as such (a real, publishable negative bounding the chemotype).

## Honest caveats
In-silico here is hypothesis-generating, not predictive (orphan + similarity-bound + method
disagreement). The promiscuity flag means any single moderate readout needs a counter-screen
before it is called a hit. Test enantiomers separately throughout — the chirality difference is
the cleanest internal control that a signal is specific.
