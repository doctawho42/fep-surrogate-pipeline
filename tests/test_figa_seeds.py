"""Structural guards on the frozen multi-seed foil sweep (peer-review item P6c).

The training run itself is exercised by `make figAseeds`, not by the unit test: retraining four
foils at five seeds is far too slow for the suite. These tests pin the frozen constants and the
output contract so a later edit cannot silently reselect seeds or change the reported foils.
"""
from __future__ import annotations

import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "figs"))

import make_figA as M  # noqa: E402


def test_seed_set_is_frozen():
    assert M.N_FOIL_SEEDS == 5
    assert M.RNG_SEED == 20260629


def test_separation_sweep_matches_the_figure_sweep():
    # the spread must be measured over the same overlap sweep the figure uses
    assert len(M.FOIL_SEPS) == 11
    assert float(M.FOIL_SEPS[0]) == 0.8
    assert float(M.FOIL_SEPS[-1]) == 3.2


def test_foil_seed_spread_reports_one_row_per_seed_with_all_four_foils():
    """Runs the real function at a deliberately tiny budget purely to pin the output CONTRACT
    (row count, seed values, key set). The ratios themselves are not asserted here; the frozen
    full-budget numbers live in docs/results_figA_seeds.md."""
    rows = M.foil_seed_spread(seps=[1.5], n=8, reps=6, n_seeds=2, seed0=M.RNG_SEED)
    assert [r["seed"] for r in rows] == [M.RNG_SEED, M.RNG_SEED + 1]
    for r in rows:
        assert set(r) == {"seed", "plain", "oracle", "betanll", "ensemble"}
        for k in ("plain", "oracle", "betanll", "ensemble"):
            assert r[k] > 0.0
