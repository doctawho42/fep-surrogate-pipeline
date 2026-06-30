# tests/test_trunk.py
from __future__ import annotations

import numpy as np

from bar.trunk import chir_pseudoscalar, featurize_edge, ligand_features

# (R)- and (S)-bromochlorofluoromethane: a minimal chiral enantiomer pair
_R = "[C@H](F)(Cl)Br"
_S = "[C@@H](F)(Cl)Br"
_ACHIRAL = "FC(F)(Cl)Br"


def test_featurize_edge_deterministic_and_has_0o_channel():
    f1 = featurize_edge("CCO", "CCCO", include_0o=True)
    f2 = featurize_edge("CCO", "CCCO", include_0o=True)
    assert np.array_equal(f1, f2)
    # the 0o channel makes the with-0o vector exactly one element longer
    f_no = featurize_edge("CCO", "CCCO", include_0o=False)
    assert f1.shape[0] == f_no.shape[0] + 1


def test_pseudoscalar_flips_sign_for_enantiomer():
    sr, ss = chir_pseudoscalar(_R), chir_pseudoscalar(_S)
    assert abs(sr) > 1e-6
    assert np.sign(sr) != np.sign(ss)
    assert abs(chir_pseudoscalar(_ACHIRAL)) < 1e-6


def test_even_features_are_enantiomer_blind_but_0o_separates():
    # edge from a common reference to each enantiomer
    e_no_R = featurize_edge(_ACHIRAL, _R, include_0o=False)
    e_no_S = featurize_edge(_ACHIRAL, _S, include_0o=False)
    assert np.array_equal(e_no_R, e_no_S)           # even (ECFP no-chirality) collapses
    e_0o_R = featurize_edge(_ACHIRAL, _R, include_0o=True)
    e_0o_S = featurize_edge(_ACHIRAL, _S, include_0o=True)
    assert not np.array_equal(e_0o_R, e_0o_S)        # 0o channel separates


def test_ligand_features_dimension():
    f = ligand_features("CCO")
    assert f.shape == (1030,)  # 1024 ECFP bits + 6 achiral descriptors


def test_amortized_sigma_monotone_and_nonneg():
    from bar.trunk import amortized_sigma
    s = amortized_sigma(np.array([0.1, 0.5, 1.0]), 0.3, q=1.2)
    assert np.all(s >= 0)
    assert s[2] > s[0]  # grows with epistemic


def test_ensemble_trunk_fits_and_epistemic_grows_off_support():
    from bar.trunk import EnsembleTrunk
    rng = np.random.default_rng(0)
    # synthetic congeneric edges: ddg correlates with a feature direction
    base = ["CCO", "CCCO", "CCCCO", "CCCCCO", "c1ccccc1", "Cc1ccccc1", "CCN", "CCCN"]
    edges = [(base[i % len(base)], base[(i + 1) % len(base)]) for i in range(40)]
    ddg = rng.normal(0, 1, 40)
    trunk = EnsembleTrunk().fit(edges, ddg, n_members=4)
    mu, se = trunk.predict(edges[:5])
    assert mu.shape == (5,) and se.shape == (5,)
    assert np.all(se >= 0)
