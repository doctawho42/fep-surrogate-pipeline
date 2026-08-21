import pathlib
import sys
import urllib.error

import pytest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "scripts"))
import boltz_cage_inputs as bci  # noqa: E402
from rdkit import Chem  # noqa: E402
from rdkit.Chem import Descriptors  # noqa: E402


def test_cage_four_forms_valid_and_small():
    s = bci.cage_smiles()
    assert set(s) == {"RR_OAc", "SS_OAc", "RR_OH", "SS_OH"}
    for k, smi in s.items():
        m = Chem.MolFromSmiles(smi)
        assert m is not None, f"{k} invalid SMILES"
        assert m.GetNumHeavyAtoms() < 50, f"{k} too big for Boltz binding"


def test_deacetyl_loses_one_acetyl():
    s = bci.cage_smiles()
    for ena in ("RR", "SS"):
        oac = Chem.MolFromSmiles(s[f"{ena}_OAc"])
        oh = Chem.MolFromSmiles(s[f"{ena}_OH"])
        d = Descriptors.MolWt(oac) - Descriptors.MolWt(oh)
        assert abs(d - 42.0106) < 0.5, f"{ena} deacetyl mass delta {d:.2f} != ~42"


def test_enantiomers_are_mirror_not_identical():
    s = bci.cage_smiles()
    assert s["RR_OAc"] != s["SS_OAc"]


_TARGETS = {"GR", "AR", "ER", "DHODH", "AChE"}


def test_anchors_resolve_and_valid():
    a = bci.anchor_smiles()
    assert set(a) == _TARGETS
    for k, smi in a.items():
        assert Chem.MolFromSmiles(smi) is not None, f"{k} anchor invalid"


def _needs_network(fn):
    """Skip rather than fail when the sequence source is unreachable.

    Two tests below resolve receptor sequences over the network. A machine without
    connectivity, or a source that rate-limits, must not turn the suite red: the
    assertion is about the data, not about the network.
    """
    def wrapper():
        try:
            fn()
        except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError, OSError) as exc:
            pytest.skip(f"sequence source unreachable: {exc}")
    wrapper.__name__ = fn.__name__
    wrapper.__doc__ = fn.__doc__
    return wrapper


@_needs_network
def test_lbd_sequences_plausible():
    seqs = bci.lbd_sequences()
    assert set(seqs) == _TARGETS
    for k, s in seqs.items():
        assert 150 <= len(s) <= 600, f"{k} domain length {len(s)} implausible"
        assert set(s) <= set("ACDEFGHIKLMNPQRSTVWY"), f"{k} non-AA chars"


@_needs_network
def test_manifest_has_25_complexes():
    from collections import Counter
    m = bci.build_manifest()
    assert len(m["complexes"]) == 25  # 5 targets x (1 anchor + 4 cage forms)
    per = Counter(c["target"] for c in m["complexes"])
    assert all(v == 5 for v in per.values()) and set(per) == _TARGETS
    assert len([c for c in m["complexes"] if c["is_anchor"]]) == 5
    assert len({c["name"] for c in m["complexes"]}) == 25  # unique run_names
