# Impact upgrade: from "calibration buys trust" to "calibration buys FEP quality control"

## The problem this fixes
The paper is honest but reads as *an audit of a mostly-null method*: the calibrated σ does not
sharpen intervals, re-rank edges, or steer a generator; the only positive was a vague "trust in
stop/commit." Reviewers (and the M1 comment) flagged this as a weak contribution. The fix is not
more spin — it is to find the one **downstream task where calibration is genuinely load-bearing**
and prove it on real data. That task is **detection**, not weighting.

## The result (new Fig L)
A calibrated per-edge aleatoric variance is exactly the **null model** that turns thermodynamic
**cycle closure** into a properly-calibrated statistical test, separating **edge-level systematic
error** (non-convergence, hysteresis, charge/water/protonation artifacts) from **sampling noise**.
On the public OpenFE IndustryBenchmarks2024 (1145 edges, 34 systems, 3 replicates):

1. **Deployable.** GLS cycle-closure χ² with the sandwich null: median reduced χ²=0.34 (most
   systems sampling-consistent), BH-FDR flags 6 chemically-sensible systems (brd4, bace, faah,
   cdk8, hif2a, p38).
2. **The money result — the detector's usefulness IS σ-calibration.** Same test, same data: an
   overconfident learned σ (×0.15, the paper's Fig-A foil) flags 88% of systems (FPR→1, useless);
   the calibrated sandwich flags 12% (selective); too-wide bars flag 0% (no power). This is the
   concrete decision an overconfident σ silently corrupts and a calibrated one enables.
3. **Replicates prove the flags are systematic.** Pooling the 3 independent replicates *increases*
   the flagged systems' χ² toward the 3× variance-components ceiling (brd4 9.5→26.8) — impossible
   for sampling noise or merely under-sized bars; only reproducible (systematic) error does this.

**Why it escapes the paper's own ceiling.** Fig C showed edge *weighting* is second-order for
ranking (Gauss–Markov: the contrast mean is unbiased for any positive weights). Cycle-closure QC
is **inference/detection**, not weighting a fixed estimate — Gauss–Markov does not apply, and the
calibrated *value* of V_e is decisive (Panel B).

## How this reframes the paper (minimal rewiring)
- **Title/abstract:** keep the tool-first framing; add one clause — the calibrated variance is not
  just "trustworthy," it *enables a quality-control test that separates systematic from sampling
  error, which an overconfident σ cannot.*
- **Contributions list:** add a 5th — "a calibrated cycle-closure QC test whose false-positive rate
  is controlled only when σ is calibrated; validated on 1145 real edges."
- **Structure:** the audit section stays (it is the honest map). Fig L becomes the **positive
  payoff** the audit was building toward: the section order becomes *audit (where it doesn't help)
  → **Fig L (where it decisively does)** → conclusion.* The narrative arc flips from "mostly null"
  to "null for the obvious uses, decisive for the right one."
- **Conclusion:** "calibration is a decision instrument, not a performance multiplier" → "calibration
  is a decision instrument — and its sharpest use is quality control: it is the null model that makes
  FEP cycle-closure a controlled-FPR test, a capability an overconfident learned σ destroys."

## Ready-to-paste LaTeX section

```latex
\section{Where the calibrated variance pays off: cycle-closure quality control}\label{sec:qc}
The audit (Section~\ref{sec:audit}) shows the calibrated $\sigma$ is second-order for sharpness,
ranking, and generation. Here is the task where it is not: separating \emph{systematic} error
(non-convergence, hysteresis, charge/water/protonation artifacts) from \emph{sampling} noise in a
perturbation network. This is detection, not weighting, so the Gauss--Markov ceiling of
Section~\ref{sec:audit} does not apply, and the calibrated \emph{value} of $V_e$ is decisive.

\paragraph{A calibrated cycle-closure test.} A perturbation network closes its cycles: the signed
sum of $\ddG$ around any loop is zero in expectation. Fitting node potentials by GLS with the
sandwich weights $w_e=I_e^2/B_e$ (Theorem~\ref{thm:fr}) gives residuals $r_e$ and the
goodness-of-fit statistic $X^2=\sum_e r_e^2/V_e\sim\chi^2(\mathrm{dof})$, with
$\mathrm{dof}=$ the number of independent cycles. Under correctly-sampled physics
$\mathbb E[X^2/\mathrm{dof}]\!\approx\!1$ (here $<1$, as the bars are ${\sim}1.4\times$
conservative, Fig.~\ref{fig:Arep}); a reduced $\chi^2$ well above one signals edge-level
systematic error. Cycle closure is (correctly) blind to node-consistent force-field bias, which
cancels around a loop.

\paragraph{Figure~\ref{fig:L}.} On the OpenFE IndustryBenchmarks2024 ($1145$ edges, $34$ systems,
$3$ replicates), the median reduced $\chi^2$ is $0.34$ and a Benjamini--Hochberg test flags six
chemically-sensible systems (\texttt{brd4}, \texttt{bace}, \texttt{faah}, \texttt{cdk8},
\texttt{hif2a}, \texttt{p38}; Fig.~\ref{fig:L}A). The test's value is entirely a function of
calibration (Fig.~\ref{fig:L}B): the overconfident learned $\sigma$ of Fig.~\ref{fig:A}
($\times0.15$) flags $88\%$ of systems (false-positive rate $\to1$, useless), the calibrated
sandwich flags $12\%$ (selective), and over-wide bars flag none. Finally, the three independent
replicates confirm the flags are systematic rather than sampling: pooling the replicates shrinks
the bars ($\mathrm{se}\to\mathrm{se}/\sqrt3$) yet \emph{increases} the flagged systems' reduced
$\chi^2$ toward the $3\times$ variance-components ceiling (Fig.~\ref{fig:L}C), which occurs only
for reproducible error. The calibrated bottleneck thus supplies, for free, the null model that
makes FEP cycle-closure a controlled-false-positive-rate test.
```
(Add a `\begin{figure}` for `figL_calibrated_cycle_closure.pdf` with the three-panel caption from
`results_figL.md`, and one sentence in Data\&Software: "Fig.~L reproduces via `make figL`.")

## Honest scope (carry into the paper, do not hide)
- Cycle closure detects **edge-level** inconsistency only; node-consistent force-field bias is
  invisible (cancels around loops). "Systematic" here = convergence/hysteresis/edge artifacts.
- **Retrospective + correlational.** Flags coincide with the known-hard systems — supportive, but
  the decisive experiment is **prospective**: flag → repair/re-run the culprit edge → confirm the
  cycle closes and the estimate improves. State it as the natural next test (and a falsifier: if
  repairing flagged edges does *not* improve closure, the flags are noise).
- Uses the reported MBAR se (= sandwich to leading order, Fig A); the differentiable sandwich's
  role is being the per-edge, graph-native, back-propagatable, non-negative weight.

## Novelty (position honestly)
Network MLE/GLS (arsenic/cinnabar) and cycle-closure hysteresis as a QC signal are standard. New
here: (1) closure as a *calibrated* $\chi^2$ test with a **validated** per-edge null → controlled
FPR and per-edge-adaptive thresholds (vs fixed ad-hoc hysteresis cutoffs); (2) the FPR is
determined by $\sigma$-calibration; (3) the replicate-based systematic-vs-sampling separator.

## Falsifiable claims (for the paper's honest-claims ledger)
- If per-edge σ is scaled away from calibrated, the cycle-closure FPR degrades monotonically
  (Panel B is the measured curve). **Falsifier:** an overconfident σ that still yields a selective
  test.
- Flagged systems carry reproducible (systematic) closure. **Falsifier:** flagged systems whose
  closure vanishes on independent replicates (Panel C would sit on the diagonal).
- Prospective: repairing a flagged edge improves closure. **Falsifier:** repair leaves closure
  unchanged.

## Validation done (short of new MD) — `make figLval`, `figs/figL_validation.{pdf,png}`
- **Causal repair test:** removing the flagged high-|z| edges reaches reduced χ²≤1 in 3–7 guided
  removals vs 12–37 random (bace 4 vs 24, cdk8 7 vs 37, p38 3 vs 29) — the flags are the causal
  culprits, not arbitrary edges.
- **Out-of-sample:** per-edge residuals correlate r=+0.30/+0.42/+0.38 across the 3 independent
  replicates (n=1143) — the flagged inconsistency reproduces on unseen runs (sampling noise → 0).

## Packaged as a tested API
- `src/bar/qc.py` — `gls_network`, `cycle_closure_test` (→ `NetworkFit` with `reduced_chi2`,
  `p_value`, per-edge `z`), `benjamini_hochberg`, `repair_order`. Pure NumPy, no SciPy/Torch.
- `tests/test_qc.py` — 6 unit tests (clean network ⇒ χ²≈1; injected systematic edge ⇒ detected +
  localized + repaired-first; calibration controls rejection; gauge-invariance; repair stops when
  acyclic). All green.

## Experimental grounding (done for eg5; honest null on the strong claim)
`figs/analyze_eg5_accuracy.py` + `docs/results_eg5_accuracy.md`, using real ChEMBL KIF11
affinities (`data/eg5_experimental_chembl.csv`): the eg5 network agrees with experiment at
**MUE 0.79 / R 0.65**, and its clean cycle-closure χ²=0.64 is consistent (QC-clean ⇒ accurate).
But acting on the QC does **not** improve accuracy on this clean system (internal |z| vs
error-vs-exp: Pearson +0.13) — exactly the Fig L scope (cycle closure is blind to the
node-consistent FF bias that dominates a clean system's error). The strong "QC-repair improves
accuracy" claim needs a QC-flagged system **with** experimental data; those are anonymized in the
benchmark, so it is scoped as the prospective next test. This null **reinforces** the paper's
thesis (calibration = decision/QC instrument, not a performance multiplier).

## Files
- `src/bar/qc.py`, `tests/test_qc.py`
- `figs/analyze_eg5_accuracy.py`, `data/eg5_experimental_chembl.csv`, `docs/results_eg5_accuracy.md`
- `figs/make_figL.py` (+ `make figL`), `figs/figL_calibrated_cycle_closure.{pdf,png}`
- `figs/make_figL_validation.py` (+ `make figLval`), `figs/figL_validation.{pdf,png}`
- `docs/results_figL.md` (full numbers, validation, scope, gate)
- Data already in-repo: `data/openfe_replicates/combined_pymbar4_edge_data.csv`
