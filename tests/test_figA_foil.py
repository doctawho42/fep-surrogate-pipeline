# tests/test_figA_foil.py
import numpy as np
import pytest

torch = pytest.importorskip("torch")
from figs.make_figA import _train_mve, controlled_panel, _gauss_edge  # noqa: E402


def test_foil_consumes_overlap_feature():
    """The honest foil must receive the overlap scalar I (6-dim feature vector)."""
    predict = _train_mve(n_train=40, seed=123, use_overlap=True)
    rng = np.random.default_rng(0)
    xf, xr = _gauss_edge(1.5, 20, rng)
    se = predict(xf, xr)
    assert np.isfinite(se) and se > 0


def test_controlled_panel_has_fair_and_oracle_columns():
    rows = controlled_panel(seps=np.linspace(1.0, 3.0, 3), reps=120, n_boot=40, seed=7)
    assert rows, "no rows"
    for d in rows:
        for k in ("sand", "mbar", "naive", "mve", "mve_oracle"):
            assert k in d and np.isfinite(d[k]), f"missing/non-finite {k}"
        assert d["sand"] > 0


def test_controlled_panel_deterministic():
    a = controlled_panel(seps=np.array([1.5, 2.5]), reps=80, n_boot=20, seed=42)
    b = controlled_panel(seps=np.array([1.5, 2.5]), reps=80, n_boot=20, seed=42)
    assert [round(d["mve"], 6) for d in a] == [round(d["mve"], 6) for d in b]
