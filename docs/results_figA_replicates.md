# Fig A (replicates) — independent-replicate validation on real binding edges

The Fig A panel-B check bootstraps the SAME work samples (internal, not independent). Here we
use the OpenFE IndustryBenchmarks2024 **3 independent replicates per edge** to test whether the
reported per-replicate MBAR/sandwich se predicts the actual run-to-run spread. Binding ΔΔG per
replicate = complex − solvent leg; reported se = sqrt(complex_dDG² + solvent_dDG²); the truth is
the across-replicate empirical SD. `make figArep`.

## Result (1145 edges, 3 replicates each)
RMS reported se 0.837 kcal/mol vs RMS independent-replicate SD 0.592 kcal/mol -> **reported/replicate = 1.41**, bootstrap CI [1.26, 1.62].
After correcting the n=3 small-sample SD bias (E[s]=c4·σ, c4=0.886) the reported se still
over-predicts the true run-to-run SD by ~1.25×. The reported se exceeds the replicate
SD on 72% of edges; per system, 28 of 34 are conservative
(ratio >= 1) and 6 are below 1 — five of them borderline (0.76–0.99) and one
protonation-variant outlier (bace\_p3\_arg368\_in, 0.41). The calibration is therefore
conservative in aggregate but heterogeneous per system.

## Honest reading
On real protein–ligand binding edges the sandwich/MBAR uncertainty is **calibrated-to-conservative
in aggregate** against independent-replicate truth: pooled, it OVER-predicts run-to-run
reproducibility by ~1.4× and is not \emph{systematically} overconfident — the dangerous
failure mode. This is the opposite of the learned MVE head (≈7× *over*confident at realistic budget / ≈5× at large budget,
Fig A) and refutes the worry that the sampling sandwich would systematically under-state real
reproducibility. It is not uniformly conservative (6/34 systems dip below 1,
one markedly), so per-system calibration varies; but the aggregate and the vast majority of
systems are safe to act on, and no learned head matches even that.

## Scope
Reported se is OpenFE's pymbar4 MBAR uncertainty, which Fig A panel A establishes equals the
sandwich B/I² to leading order; we validate that reported quantity against replicate truth (we do
not recompute B/I² from the raw works, which are not in the released per-edge table). Robust to
overlap filtering (ratio ≈ 1.4 at smallest-overlap ≥ 0.10).
