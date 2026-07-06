# tests/test_figA_foil.py
import numpy as np
import pytest

torch = pytest.importorskip("torch")
from figs.make_figA import (  # noqa: E402
    _gauss_edge,
    _train_mve,
    _train_mve_ensemble,
    controlled_panel,
)


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


def test_controlled_panel_has_corrected_foil_columns():
    """The two corrected foils (beta-NLL, ensemble-of-MVEs) must be reported per the
    2026-07-05 mve-corrected-foil-experiment spec — closes the 'Gaussian-NLL strawman' hedge."""
    rows = controlled_panel(seps=np.linspace(1.0, 3.0, 3), reps=120, n_boot=40, seed=7)
    assert rows, "no rows"
    for d in rows:
        for k in ("mve_betanll", "mve_ensemble"):
            assert k in d and np.isfinite(d[k]) and d[k] > 0, f"missing/non-finite/non-positive {k}"


def test_beta_zero_reproduces_plain_nll():
    """beta=0.0 MUST reproduce the original (pre-beta) Gaussian-NLL behavior exactly."""
    predict_plain = _train_mve(n_train=40, seed=123, use_overlap=True)
    predict_beta0 = _train_mve(n_train=40, seed=123, use_overlap=True, beta=0.0)
    rng = np.random.default_rng(0)
    xf, xr = _gauss_edge(1.5, 20, rng)
    assert predict_plain(xf, xr) == predict_beta0(xf, xr)


def test_ensemble_se_dominates_single_member_aleatoric_se():
    """Epistemic term is non-negative: the ensemble mixture se must be >= a single member's
    aleatoric-only se on at least one evaluation point."""
    ens_predict = _train_mve_ensemble(n_train=200, seed=11, use_overlap=True, n_members=5)
    single_predict = _train_mve(n_train=200, seed=11, use_overlap=True, beta=0.0)
    rng = np.random.default_rng(3)
    found_ge = False
    for s in (0.8, 1.5, 2.2, 3.0):
        xf, xr = _gauss_edge(s, 20, rng)
        se_ens = ens_predict(xf, xr)
        se_single = single_predict(xf, xr)
        assert np.isfinite(se_ens) and se_ens > 0
        if se_ens >= se_single:
            found_ge = True
    assert found_ge, "ensemble se never exceeded a single member's aleatoric-only se"
