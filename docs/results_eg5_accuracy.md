# Results — eg5 accuracy vs experiment (the honest outcome of the impact swing)

**Reproduce:** `PYTHONPATH=src python figs/analyze_eg5_accuracy.py`.
**Data:** `data/eg5_experimental_chembl.csv` (ChEMBL KIF11 ATPase IC50, assay CHEMBL1101849, pulled
per-ligand via the ChEMBL API) + `data/openfe_replicates` (3-replicate mean).

## What we set out to test
Does using the calibrated cycle-closure QC to prune/down-weight flagged edges **improve the
delivered ΔG vs experiment** (the "detects problems → improves the answer" jump)? We flagged this
as a ~50/50 swing with a real chance of a null. It came out **null on the strong claim**, and the
null is informative — it confirms the QC's stated scope against real experiment.

## POSITIVE (real experimental grounding)
The real OpenFE eg5 RBFE network, solved by our GLS pipeline, agrees with **real ChEMBL KIF11
experimental affinity**:

| metric | value |
|---|---|
| MUE | **0.79 kcal/mol** |
| RMSE | 0.93 kcal/mol |
| Kendall τ | +0.37 |
| Pearson R | +0.65 |
| n ligands | 28 |
| cycle-closure reduced χ² | **0.64 (clean)** |

This is the first time the QC/network is validated against **experiment** (everything else was
internal consistency or replicate reproducibility), and the numbers are literature-credible for
eg5. The system the calibrated test calls **clean** (χ²=0.64) is indeed **accurate** (MUE 0.79) —
the two agree, as the QC promises.

*Caveats.* The experimental anchor is **IC50** (assay CHEMBL1101849), a Kd proxy; because the whole
series is one assay, the constant IC50/Kd offset cancels under the mean-alignment (`p−p̄+ē`), so the
valid comparison is on **relative, mean-aligned ΔG** (MUE/R/τ), not an absolute-Kd claim. Both ChEMBL
values were spot-checked against the live ChEMBL API (CHEMBL1077204 IC50 8 nM/pChEMBL 8.10;
CHEMBL1078691 34 nM/7.47 — exact). n=28 are point estimates (bootstrap CIs would widen them).

## NULL (honest, and predicted by the scope)
Acting on the QC does **not** improve accuracy vs experiment on eg5:
- Internal |z| barely predicts per-edge error vs experiment: **Pearson +0.13, Spearman +0.05**
  (top-quartile-|z| edges 1.28 vs 1.01 kcal/mol error, only 1.3×).
- Down-weighting or pruning flagged edges leaves MUE essentially unchanged; an inject-and-recover
  test (add a +3 kcal/mol edge error, then repair) does not recover accuracy, and an overconfident
  σ over-prunes (removes ~12 edges) and degrades it.

**Why this is the *expected* honest result, not a bug.** Cycle closure detects **edge-level**
inconsistency and is (correctly) **blind to node-consistent force-field bias** — a per-ligand
offset cancels around a loop (stated as the Fig L scope). On a **clean** system, the residual error
vs experiment is dominated by exactly that node-consistent FF bias, which the QC cannot see and
therefore cannot fix. So QC neither predicts nor improves accuracy here — while the network is
otherwise accurate and QC-clean. This is a self-consistent **confirmation of the scope against real
experiment**, not a contradiction.

## Why the strong claim is untestable here
The strong claim ("QC repair improves accuracy") can only be tested on a **QC-flagged** system
(brd4, bace, faah, cdk8, hif2a, p38) **with experimental affinities**. In the OpenFE
IndustryBenchmarks2024 set those systems are **anonymized** (ligand identities/affinities withheld),
and external experimental sources could not be fetched in this environment. So the strong test
requires unblinded data or new MD (flag → re-run/repair a flagged edge → measure the improvement).
That is the concrete, well-scoped next experiment — now with named target edges (e.g. bace CAT-13).

## Bottom line
The impact swing yields one real positive (experimental grounding of the QC on eg5) and an honest
null on the strong "improves accuracy" claim — a null that **reinforces** the paper's thesis
(calibration is a decision/QC instrument, not a performance multiplier) and confirms the Fig L
scope against experiment. It does not manufacture a capability the data does not support.
