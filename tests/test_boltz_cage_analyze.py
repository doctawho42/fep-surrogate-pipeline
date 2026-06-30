import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "figs"))
import make_boltz_cage as mbc  # noqa: E402


def test_anchor_gap_is_cage_minus_anchor():
    conf = {"anchor": 0.80, "RR_OH": 0.62, "SS_OH": 0.50}
    g = mbc.anchor_gap(conf, "anchor")
    assert abs(g["RR_OH"] - (-0.18)) < 1e-9
    assert abs(g["SS_OH"] - (-0.30)) < 1e-9
    assert "anchor" not in g


def test_verdict_corroborates_when_GR_near_anchor_and_dhodh_rejected():
    gaps = {"GR": {"RR_OH": -0.05}, "AR": {"RR_OH": -0.20},
            "ER": {"RR_OH": -0.35}, "DHODH": {"RR_OH": -0.55}}
    assert "corroborat" in mbc.verdict(gaps).lower()


def test_verdict_challenges_when_GR_far_or_dhodh_not_rejected():
    gaps = {"GR": {"RR_OH": -0.40}, "DHODH": {"RR_OH": -0.10}}
    assert "corroborat" not in mbc.verdict(gaps).lower()
