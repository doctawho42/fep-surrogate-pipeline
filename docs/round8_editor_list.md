# Round 8 editor: the 13 required changes (verbatim)

Recovered from the panel journal. Untracked on purpose: this is the review record, not a
manuscript source. Status is maintained here as items land.

## 1

Rewrite the operating-point paragraph under one declared convention, stating the mapping (simulation scale = 1/calibration scale). Print the realized family-wise level as 0.000, 0.003, 0.045 and 0.100 at calibration scales 0.79, 0.92, 1.00 and 1.04, and label 0.215 and 0.912 as the cost of deliberately tightening a correct bar by eight and twenty-one per cent, not as points of the measured interval. Do NOT lead on the 0.003 point estimate: the anti-conservative end of your own measured interval realizes 0.100, twice nominal, and that is the number a reader needs. Fix OPERATING_SCALES and its comment, re-emit the results record, and print the draw count with a binomial interval at each grid point — 400 draws cannot resolve 0.000 or 0.003, and the same applies to Table S5's 0.003, which the body quotes as 'the only self-calibrated arm holding its level'. Separately, 'Across it the flag count runs from six to nine' is the one number in that paragraph with no released generator and it sits among numbers computed under the inverted convention; regenerate it and show it.

## 2

Rewrite Section 5 around a null rather than footnoting one. Relabel the 0.15x arm a counterfactual stand-in in Section 5, the Figure S16 caption, Table S7 and contributions bullet 5, using the words Section 4 now uses; print the assumed-sigma sweep; add Table S7 rows for the temperature-scaled and split-conformal baselines with their intervals; delete Table S7's claim of a clear benefit for the two decision tasks; and rewrite the Synthesis and contribution bullets 4 and 5 to state that the physics sigma matches rather than beats post-hoc recalibration on commit trust, and that the stopping separation is absent at the six per cent a same-budget pooled estimator produces.

## 3

Rebuild the abstract. JCTC requires three to four concise sentences; yours is 275 words in ten. In rebuilding it: replace 'only removes a fifth of its flags' with the counts ('removes 20 of its 26 flags'); name the replicate those counts come from, since the calibrated set is 6, 6 and 3 across the three and its membership churns; name the metric behind 'under one per cent' and quote the weakest reading (a factor of 1.36 below chance on bace under single-edge deletion, isotropic in kcal/mol) rather than the pooled point estimate; report both tails of the removal-rule result or neither; and carry the 1.5x/2x widening collapse. Move the Figure 1 float below the abstract so it stops splitting it across two pages.

## 4

Put the pose-and-preparation ceiling in the abstract and the Conclusion in your own words, and retire the two sentences that overclaim on it: 'a false-positive rate controlled to the extent the per-edge null is calibrated (validated against independent replicates, Figure 3)' and 'a false-positive rate tied to a replicate-validated null'. Both attach a controlled error rate to a null whose denominator excludes exactly the variance component the alternative hypothesis is built from, and nothing on this benchmark flags at a 1.5x widening. Replace 'replicate-validated' with the scoped phrase wherever it stands alone, and state the exclusion at each point.

## 5

Restrict the nesting claim by replicate. Qualify the Supporting Information's unqualified 'strictly nested' to replicate 0; name the replicate-1 break and the two systems responsible (ciordia_retro, mup1); restrict the abstract to the containment that holds on all three (the calibrated set inside the fixed cutoff's). State both detector definitions in Computational details — the aggregation rule for the 1.0 kcal/mol cutoff and the pooled-se chi-squared test appear nowhere in eighty-one pages, and three flag counts plus a headline claim rest on them. Report the three counts under the 56-system release grouping, which has not been run for the two baseline detectors.

## 6

Replace the TOC graphic. Its entire text is 'Which networks does the error bar flag? / calibrated / overconfident / 6/48 / 42/48 / 48 benchmark systems, cycle-closure test', with no stand-in label, in the one element of an ACS article that travels without its caption and whose comparison the body now calls a counterfactual. Replace it with the auditable-versus-unauditable split and its measurement against chance. Figure 1's artwork already reads 'overconfident stand-in' and needs no change.

## 7

Correct the multiple-testing justification. Delete 'the independence Benjamini-Hochberg needs'; let Benjamini-Yekutieli carry the guarantee, since it assumes no dependence structure and returns the identical six flags, and offer positive regression dependence as a remark with shared preparation and force field named as the mechanism. Fix the capital-after-semicolon splice the same commit left ('...dependent by construction; The nominal level...').

## 8

Print the visible-fraction bookkeeping in full: one table, four grounded systems by three conventions (whitened, isotropic in kcal/mol, unweighted), each with its matching chance level and its leave-one-edge-out range. State in the Supporting Information that the unweighted reading recomputes the projector in the unweighted metric, and give those numbers a released generator — the 5.45% and the per-system unweighted values currently have none. Carry the cdk8 exception into the Section 4.1 framing sentence and the Figure 6 caption, printing the pooled fraction with and without cdk8 at the point of use.

## 9

Repair the pre-registration apparatus and then audit it. The digest prints 32 of its 64 characters in the submitted PDF, so nothing can be checked against what a reader sees — break the string so it renders, and verify in the rendered output. Footnote the reproduced file's theorem-label errors and reconcile its 'external' and 'Immutable' wording with the prose, stating plainly that the anchor is in-repository, that no third-party timestamp exists, and that the claim is drift detection only. Print the pooled adjusted R-squared (0.980) beside the four per-system values and name the pooling convention. Correct the released adjudicator, which tested the unadjusted pooled criterion rather than the dof-adjusted rule the pre-registration froze, and say so in print. Then audit the remaining frozen criteria that live only in code — the |z| conjunctive criterion, the curl-leverage falsifier rule, the head-to-head criterion — against their printed text and report the outcome. One mismatch found by a referee is a defect; leaving the rest unchecked after that is a policy.

## 10

Fix the specific errors: eg5 at 0.79 kcal/mol beats three of the four flagged systems, not two, and ties p38 — correct the count and print the full ordering. Table S5's footnote 'it is the only arm that loses any' is false against the row two lines above it, where the per-edge one-sided arm retains two of the published six. Name misc_cdk8 in the cdk8 paragraph with its size and its three per-replicate p-values, including the 0.091 that is the smallest of the six. Complete references 31 and 34, which print the literal string 'others'. Remove the stray ';.' from the Supporting Information Available list.

## 11

Cite the leverage and effective-resistance identity to its own sources. Klein and Randic's resistance-distance paper does not contain a hat matrix or a leverage correspondence, and no graph-sparsification or hat-matrix source appears anywhere in the bibliography. Add a graph-sparsification source (Spielman and Srivastava on sparsification by effective resistances) and a hat-matrix source (Hoaglin and Welsch on the hat matrix in regression), keeping Klein and Randic for resistance distance itself and for the Theorem 3 identification, which is a fair one-line reading. I have not verified the volume, pages or year of either new entry and neither should you until you have checked them against the publisher record.

## 12

Report the outcome of the extended priority search into the transitivity/consistency and bias-adjustment strands of network meta-analysis, as found or not found, and add the paragraph stating what Theorem 5(iv) adds over Xu 2019's Fisher-information treatment and Giese and York 2021's constrained solution — or withdraw the claim if it is contained in either. Move Theorem 1 and its contributions bullet to the Supporting Information as you have agreed.

## 13

Proofread the rendered PDF, not the source. The truncated digest, the stray ';.', the two 'others' author lists and the semicolon splice all survived because the revision was checked against the LaTeX. I will check the resubmission against the compiled output.
