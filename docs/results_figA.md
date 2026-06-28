# Results — Fig A: "target the sandwich"

**Figure:** `figs/figA_target_the_sandwich.{pdf,png}` · **Reproduce:** `make figA`
(or `python figs/make_figA.py`). Deterministic (seed 20260629).

![Fig A](../figs/figA_target_the_sandwich.png)

## Claim tested (plan Fig A)
The BAR uncertainty is the **sandwich** `Var = B/I²`, not the naive `1/I`. Across
overlap regimes the sandwich is calibrated (reported-se ≈ true-se), while `1/I` is
over-conservative by a factor that **varies** with overlap — so no constant rescale
can fix it; a per-edge sandwich is required.
**Falsifier:** sandwich not ≈ 1, or the `1/I` error is a *constant* factor.

## Result: PASS

**Panel A — controlled (MC truth).** Gaussian work model, `n_f=n_r=20`, 3000
replicates/point; true se = MC SD across replicates; bands = 95% bootstrap CI.

| overlap `4⟨p(1−p)⟩` | sandwich/true | naive `1/I` /true |
|---:|---:|---:|
| 0.17 (low overlap)  | 0.92 | **1.10** |
| 0.46                | 1.01 | **1.38** |
| 0.61                | 0.98 | **1.57** |
| 0.79                | 1.02 | **2.17** |
| 0.86 (high overlap) | 0.99 | **2.65** |

Sandwich/true ∈ [0.92, 1.02] across the whole range (calibrated). `1/I` ranges from
1.10 to 2.65 — a **non-constant** factor.

**Panel B — real FEP edges (benzene hydration, alchemtest GROMACS).** 19 adjacent-λ
BAR edges (Coulomb 4 + VDW 15); works decorrelated by the statistical inefficiency;
true se = autocorrelation-aware bootstrap on the decorrelated works (n ≈ 3.5–4k/edge).

- **sandwich/bootstrap = 1.007**, range [0.945, 1.076] — calibrated on real data.
- **naive/bootstrap = 4.0 mean**, rising monotonically **1.76 → 13.5** as overlap
  goes 0.69 → 0.99.

Real FEP edges live in the **high-overlap** regime (FEP is designed that way), which
is exactly where `1/I` is worst — so the correction is large (2–13×) on real edges.

## Direction of the effect — and a flag on the proofs sheet Corollary

The terminology-free fact, confirmed **three independent ways** (controlled MC,
population-level `I/B`, real-edge bootstrap):

> `1/I` over-estimates the BAR variance **most when the two work distributions
> overlap most** (small separation, `p≈½`), and converges to the correct variance
> when they overlap least.

Population-level `naive/sandwich = I/B` (single huge samples, no plug-in noise):

| separation | overlap coef | `I/B` (= naive/sandwich) |
|---:|---:|---:|
| 0.5 σ | 0.80 | **16.9** |
| 1.0 σ | 0.62 | 4.9 |
| 2.0 σ | 0.32 | 1.8 |
| 3.0 σ | 0.13 | 1.25 |

This **matches** the proofs sheet's Theorem-2 numerical remark (sep = 1σ → naive
2.21×) and `CLAUDE.md` invariant #1 ("2.2× in se at **high** overlap, shrinking
toward low overlap").

**⚠ Flag (not papered over):** the proofs sheet *Corollary* ("the unified objects")
states the opposite direction — "in OOD / **low-overlap** regimes `B_e` and `I_e`
diverge most, so the naive `1/I` is worst exactly where precision matters." The
evidence shows `B` and `I` diverge most at **high** overlap (`I/B`: 16.9 at high →
1.25 at low). The Corollary's *math* (V = B/I² coincides with 1/I only at info-
equality) is correct; only its **directional wording** is backwards and conflicts
with the same sheet's Theorem-2 remark and with CLAUDE.md. Recommended correction:
replace "low-overlap" with "high-overlap" in that Corollary sentence.

**Why this strengthens, not weakens, the paper.** The sandwich correction matters
*most* for the well-overlapped, low-variance edges — precisely the reliable edges an
active learner wants to exploit. Naive `1/I` would systematically *distrust the best
edges* (inflating their tiny variance up to ~13×), and mis-weight the Fisher–
resistance Laplacian (`w_e = I_e²/B_e`) on exactly those edges. The central Fig A
claim (sandwich calibrated; `1/I` off by a non-constant factor) holds robustly.

## Honest scope & next steps
- Panel B edges are **solvation** (hydration) FEP, not protein–ligand binding. The
  calibration result is a property of BAR statistics and is domain-agnostic, but a
  *binding* panel is stronger for the paper. Confirmed extensions (loaders in hand):
  `alchemtest.amber.load_bace_example` (a real BACE1 RBFE edge, u_nk extractable) and
  the **OpenFE IndustryBenchmarks2024** Zenodo archives, which ship raw `u_ln.txt`
  per edge **with 3 repeats** — giving a true independent-replicate se reference
  (vs the within-window bootstrap used here). Logged as the Fig A strengthening task.
- True se on real edges here is bootstrap-based (no repeats for benzene); the
  controlled panel supplies the independent MC-truth validation.

## Gate
Full test suite green (`make check`: 26 passed) **and** Fig A regenerable by one
command (`make figA`). Sandwich is calibrated → **proceed to Fig E (chirality)**.
