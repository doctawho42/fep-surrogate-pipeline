# fep-surrogate-pipeline

Code, data and figures for *What Cycle Closure Can and Cannot See in a Relative Binding Free
Energy Network*.

A relative binding free energy calculation reports a free energy for each edge of a perturbation
network together with a standard error. Because the network must close its thermodynamic cycles,
the failure to close is routinely used as a consistency check. This work uses the reported
per-edge standard errors as the null of that check, which turns a descriptive diagnostic into a
test with a stated null, and then asks what such a test can see at all.

The answer is an algebraic bound. Per-edge error splits into a part carried by the network's cycles
and a part that is a difference of per-ligand offsets. The two are orthogonal. Only the first is
identifiable from the network's own edge values, and it is the second that biases the estimated
affinities. On the public OpenFE replicate benchmark the auditable subspace spans about a third of
the error space and carries under one per cent of the measured error against experiment. Resolving
the invisible part requires measurements from outside the network.

## Installation

Python 3.11 or newer. Dependencies are pinned to the versions the results were produced with.

```bash
python -m venv .venv && source .venv/bin/activate
pip install -e ".[dev,fig]"
```

The `fig` extra adds the alchemical and cheminformatics packages that the figure scripts read data
through (`alchemtest`, `alchemlyb`, `pandas`, `rdkit`, `scikit-learn`). The core library needs none
of them.

Every command below runs on the interpreter named by `PY`, which defaults to `python`. Point it at
the environment explicitly if you prefer:

```bash
make PY=/path/to/.venv/bin/python check
```

## Reproducing the figures

Each figure is produced by one command, which also writes a record of its numbers to `docs/` so
that a claim in the article can be traced to the run that produced it.

| Article figure | Command | Output |
|---|---|---|
| 1, overview | `make graphical` | `figs/graphical_abstract.pdf` |
| 2, sandwich variance against learned foils | `make figA` | `figs/figA_target_the_sandwich.pdf` |
| 3, replicate validation of the reported error | `make figArep` | `figs/figA_replicate_validation.pdf` |
| 4, the calibrated cycle-closure detector | `make figL` | `figs/figL_calibrated_cycle_closure.pdf` |
| 5, edge removal and cross-replicate residuals | `make figLval` | `figs/figL_validation.pdf` |
| 6, what the audit covers and where the error lives | `make figHodge` | `figs/figHodge_where_the_error_lives.pdf` |

The sixteen Supporting Information figures are produced by `make gsweep`, `figC`, `figD`,
`qcstruct`, `figLev`, `figLcausal`, `figCut`, `figStab`, `figOOS`, `figSelf`, `figP4`, `figP4b`,
`figK`, `figBreal`, `figJ` and `figG`. `make help` lists every target with a one-line description.

Two of these need more than the repository. `make qcstruct` fetches PDB entries 4DJW and 3MXF at
run time, and `make figH` runs docking and expects a `smina` binary, found through the `SMINA`
environment variable or on `PATH`.

## Building the manuscript

```bash
make jctc     # ACS submission build  -> docs/paper_jctc.pdf + docs/paper_jctc_si.pdf
make paper    # generic build         -> docs/paper_draft.pdf
```

Both compile the same body (`docs/paper_body.tex`) and Supporting Information
(`docs/paper_si.tex`) through different wrappers. ACS requires the Supporting Information as a
separate file, so `make jctc` produces two: the manuscript, and an S-numbered Supporting
Information with its own title page. The two reference each other's figures and sections through
`xr`, which reads the other document's `.aux`, so the target builds them in alternating passes
until both settle; a single pass leaves the cross-document numbers stale rather than undefined.

## Layout

```
src/bar/         the estimator and the quality-control instrument
  estimator.py     BAR root-find, sandwich variance B/I^2, information shares
  torch_layer.py   the same estimator as an autograd primitive, O(1) backward
  qc.py            generalized-least-squares network fit, closure chi^2, FDR control
  hodge.py         the gradient-cycle split, gradient R^2, influence ranking
  leverage.py      curl-leverage, the estimation-detection conservation law
  detectors.py     the calibrated null and the two baseline detectors
  closeloop.py     grounding against measured affinity, guided-versus-random removal
src/gen/         the trajectory-balance generator used in one downstream audit
src/screen/      the target-identification arm, reported as a negative result
figs/            one script per figure, plus the shared style module
docs/            the manuscript, the theorem notes, and one results record per figure
data/            the benchmark tables, the affinity joins, and the pre-registrations
tests/           the correctness invariants, as unit tests
```

`figs/paperstyle.py` holds the figure typography, the palette and the canonical widths, so a
figure's appearance is set in one place rather than per script.

## Data

Everything is public and nothing was recomputed. The alchemtest BACE1 and benzene alchemical sets
supply the work samples behind Figure 2. The OpenFE IndustryBenchmarks2024 release supplies 1145
binding edges over 49 systems with three independent protocol repeats each. The OpenFF
protein-ligand benchmark, which aggregates the FEP+ and Merck data sets, supplies the downstream
audit tasks. ChEMBL supplies the experimental affinities used for grounding. No free energy
simulation was run for this work.

The pre-registration files under `data/openfe_replicates/` fix the hypotheses, statistics and
success criteria of the decomposition and close-the-loop experiments. `docs/prereg_provenance.md`
records their digests and states exactly what those digests do and do not establish.

## Tests

```bash
make check     # ruff, mypy, and the test suite
```

288 tests over 32 files. They assert the correctness invariants rather than smoke-testing: that the
variance is the sandwich rather than the information-equality plug-in, that the autograd backward
matches finite differences, that curl-leverage sums to the degrees of freedom, that bridge edges
carry zero leverage, and that the Benjamini-Yekutieli flag set is contained in the
Benjamini-Hochberg one.

Two of the 288 resolve receptor sequences over the network and skip when it is unavailable, so
an offline run reports 286 passed and 2 skipped.

`python tests/verify_proofs.py` checks the theorem statements numerically, independently of the
suite.

## Citation

The article is under review. Until it appears, cite this repository at the tagged version.

```
N. L. Polomoshnov, What Cycle Closure Can and Cannot See in a Relative Binding
Free Energy Network, 2026. Zenodo. https://doi.org/10.5281/zenodo.22046692
```

The archive is https://doi.org/10.5281/zenodo.22046692, a concept identifier that resolves to the
most recently deposited version. Tag `v1.1.0` is the state every figure in the article was
generated from; `v1.0.0` is the state of the first submission.

## Licence

MIT. See `LICENSE`.
