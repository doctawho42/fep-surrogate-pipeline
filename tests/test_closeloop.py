"""Unit tests for the close-the-loop accuracy race (src/bar/closeloop.py)."""
from __future__ import annotations

import pathlib

import pytest

from bar import closeloop as C

# external witness; independent literal, never edited
PINNED_SHA256 = "35cb3b5a81c6878fddaa2f8cfc0107596fbdd9f1c6de3ab132524632c97ff86e"


def test_prereg_loads_and_matches_pinned_anchor():
    assert C.CLOSELOOP_PREREG_SHA256 == PINNED_SHA256
    assert C.sha256_of(C.CLOSELOOP_PREREG_PATH) == PINNED_SHA256
    pr = C.load_prereg()
    assert pr.systems == ["cdk8", "hif2a", "p38", "bace"]
    assert pr.coverage_min == 0.60
    assert pr.n_perm == 1000
    assert pr.min_grounded == 2
    assert pr.assay_order == ["IC50", "Ki"]


def test_prereg_immutable_against_external_anchor(tmp_path):
    tampered = tmp_path / "t.yaml"
    tampered.write_bytes(pathlib.Path(C.CLOSELOOP_PREREG_PATH).read_bytes() + b"\n# edit\n")
    with pytest.raises(ValueError):
        C.load_prereg(str(tampered), expected_sha256=PINNED_SHA256)
