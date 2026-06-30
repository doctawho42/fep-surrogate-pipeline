# tests/test_trunk.py
from __future__ import annotations

import numpy as np

from bar.trunk import chir_pseudoscalar, featurize_edge, ligand_features  # noqa: F401

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
