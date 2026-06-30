import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "scripts"))
from rdkit import Chem  # noqa: E402
from rdkit.Chem import Descriptors  # noqa: E402

import boltz_cage_inputs as bci  # noqa: E402


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
