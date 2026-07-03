# Results — structural context of the cycle-closure QC flags (Fig L supplement)

**Figure:** `figs/figQC_structures.{pdf,png}` · **Reproduce:** `make qcstruct`
(`PYTHONPATH=src python figs/analyze_qc_structures.py`). Deterministic. Public PDBs cached in
`data/pdb/` (`4DJW`, `3MXF`), fetched from RCSB.

## The point (and the honest limit)
Fig L flags six systems (brd4, bace, faah, cdk8, hif2a, p38) as edge-level over-dispersed. The
paper calls them "chemically sensible." This note substantiates that with **real co-crystal
structures**: the flagged systems are the ones with an *independently-documented* hard structural
feature in the binding site. It is **retrospective structural context, not a per-edge causal
claim** — cycle closure is edge-level and blind to node-consistent bias (Fig L scope), so this
cannot say a particular water sits on a particular flagged edge. The claim is only that the
flagged set coincides with the features a domain expert would flag a priori.

## Concrete structural facts (data-driven, from the PDBs)

| system | PDB (ligand) | documented hard feature | measured from the structure |
|---|---|---|---|
| **BACE1** | 4DJW (0KP) | aspartic-protease catalytic dyad; protonation-state dependent | Asp93 **2.67 Å** and Asp289 **2.72 Å** from the ligand — the dyad H-bonds the inhibitor |
| **BRD4 BD1** | 3MXF (JQ1) | conserved acetyl-lysine-pocket water network | **8** waters ≤3.5 Å, **16** ≤4.5 Å of the ligand; recognition Asn140 at **3.15 Å** |
| faah | — | covalent serine hydrolase, large acyl channel | (literature; not measured here) |
| cdk8 | — | large scaffold changes | (Wang 2015 set) |
| hif2a | — | buried polar cavity | (literature) |
| p38 | — | large R-group swaps | (Wang 2015 set) |

BACE1 (aspartic dyad in direct contact → protonation) and BRD4 (buried KAc waters) are the two the
paper names; both are read straight from the crystal structure. Citations: BRD4 JQ1 complex
`filippakopoulos2010bet` (Nature 2010, DOI 10.1038/nature09504); the BRD4 conserved-water free-energy
challenge `aldeghi2016binding` (Chem. Sci. 2016, DOI 10.1039/c5sc02678d) — both verified against
Crossref. BACE1 catalytic-dyad protonation is the textbook aspartic-protease difficulty; the
structural fact (dyad H-bonding the ligand at ~2.7 Å) is read from PDB 4DJW.

## Honest scope
- **Retrospective + interpretive.** These are known-hard features, and the QC recovers exactly the
  known-hard systems; that is the sense of "chemically sensible." It is not proof that the specific
  systematic error the QC detects is *caused* by these features on the flagged edges. The decisive
  test is still prospective (re-run/repair a flagged edge with fresh MD; main text).
- **Two of six measured.** Only BACE1 and BRD4 (the systems the paper names) are read from
  structures here; the other four rely on the benchmark's documented chemistry.
- Distances are heavy-atom minimum distances to the co-crystallised ligand in a single public
  structure; they establish contact/proximity, not energetics.
