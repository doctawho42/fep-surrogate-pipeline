"""Invariants of the curated-label join and the visible fraction measured on it.

Everything here is offline: the one network access in ``bar.curated`` is injectable, and the SDF
blobs the tests parse are built in memory by RDKit.
"""
from __future__ import annotations

import csv

import numpy as np
import pytest
from rdkit import Chem

from bar.curated import (
    Match,
    NetworkEdge,
    canonical_smiles,
    curated_edges,
    match_system,
    name_key,
    network_edges,
    pool,
    read_ligand_cache,
    read_system,
    resolve_ligands,
    sdf_smiles,
    visible_fraction,
    write_ligand_cache,
)

ETHANOL = "CCO"
PROPANOL = "CCCO"
BUTANOL = "CCCCO"
PENTANOL = "CCCCCO"


def _sdf(entries: dict[str, str]) -> str:
    """An SDF blob with one record per ``name -> SMILES`` entry."""
    blocks = []
    for name, smiles in entries.items():
        mol = Chem.MolFromSmiles(smiles)
        mol.SetProp("_Name", name)
        blocks.append(Chem.MolToMolBlock(mol) + "$$$$\n")
    return "".join(blocks)


# --------------------------------------------------------------------------------------
# The join keys.
# --------------------------------------------------------------------------------------
def test_name_key_is_separator_insensitive():
    assert name_key("30 flip") == name_key("30-flip") == name_key("30_flip")
    assert name_key("CHEMBL3402756_2.7 redocked") == name_key("CHEMBL3402756_2.7_redocked")
    assert name_key("  Example  7 ") == name_key("Example-7")


def test_name_key_does_not_merge_distinct_names():
    assert name_key("lig1") != name_key("lig2")
    assert name_key("30 flip") != name_key("30")


def test_canonical_smiles_is_a_structure_key():
    assert canonical_smiles("OCC") == canonical_smiles("CCO") == canonical_smiles("[H]OCC")
    assert canonical_smiles("not a molecule") is None


def test_canonical_smiles_keeps_stereochemistry_but_neutralize_drops_charge():
    left, right = canonical_smiles("N[C@H](C)C(=O)O"), canonical_smiles("N[C@@H](C)C(=O)O")
    assert left != right, "the exact key must separate enantiomers"
    charged = canonical_smiles("[NH3+][C@H](C)C(=O)O", neutralize=True)
    assert charged == canonical_smiles("N[C@H](C)C(=O)O", neutralize=True)
    assert canonical_smiles("[NH3+][C@H](C)C(=O)O") != canonical_smiles("N[C@H](C)C(=O)O")


def test_sdf_smiles_keys_by_normalized_name():
    parsed = sdf_smiles(_sdf({"30 flip": ETHANOL, "lig 2": PROPANOL}))
    assert parsed[name_key("30-flip")] == canonical_smiles(ETHANOL)
    assert parsed[name_key("lig_2")] == canonical_smiles(PROPANOL)


# --------------------------------------------------------------------------------------
# The two edge sources.
# --------------------------------------------------------------------------------------
def test_curated_edges_drops_unparseable_and_self_pairs(tmp_path):
    path = tmp_path / "toy.csv"
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["smiles_a", "smiles_b", "ddg"])
        writer.writerow([ETHANOL, PROPANOL, 1.5])
        writer.writerow([ETHANOL, "OCC", 9.9])          # the same structure twice
        writer.writerow([ETHANOL, "not a molecule", 3.0])
    edges = curated_edges("toy", tmp_path)
    assert edges == [(canonical_smiles(ETHANOL), canonical_smiles(PROPANOL), 1.5)]


def _combined_row(system: str, group: str, a: str, b: str, ddg: float, se: float,
                  replicates: int = 3) -> dict:
    row = {"ligand_A": a, "ligand_B": b, "system name": system, "system group": group}
    for k in range(3):
        present = k < replicates
        row[f"complex_repeat_{k}_DG (kcal/mol)"] = str(ddg) if present else ""
        row[f"complex_repeat_{k}_dDG (kcal/mol)"] = str(se) if present else ""
        row[f"solvent_repeat_{k}_DG (kcal/mol)"] = "0.0" if present else ""
        row[f"solvent_repeat_{k}_dDG (kcal/mol)"] = "0.0" if present else ""
    return row


def _write_combined(path, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def test_network_edges_averages_replicates_and_needs_two(tmp_path):
    path = tmp_path / "combined.csv"
    _write_combined(path, [
        _combined_row("sysA", "grp", "L1", "L2", 2.0, 0.3),
        _combined_row("sysA", "grp", "L2", "L3", 1.0, 0.3, replicates=1),
        _combined_row("sysB", "grp", "L1", "L2", 5.0, 0.3),
    ])
    edges = network_edges("sysA", path)
    assert len(edges) == 1, "an edge with one usable replicate carries no averaged se"
    assert edges[0].ddg == pytest.approx(2.0)
    # se = sqrt(mean se_k^2 / n) with se_k = sqrt(0.3^2 + 0^2) and n = 3
    assert edges[0].se == pytest.approx(0.3 / np.sqrt(3))


# --------------------------------------------------------------------------------------
# The match.
# --------------------------------------------------------------------------------------
def _net(pairs, group="grp"):
    return [NetworkEdge(a=a, b=b, ddg=ddg, se=se, group=group) for a, b, ddg, se in pairs]


def test_match_orients_the_curated_value_to_the_network_edge():
    resolved = {"grp": {name_key("L1"): canonical_smiles(ETHANOL),
                        name_key("L2"): canonical_smiles(PROPANOL)}}
    curated = [(canonical_smiles(PROPANOL), canonical_smiles(ETHANOL), 1.0)]  # L2 -> L1
    match = match_system(_net([("L1", "L2", 0.25, 0.1)]), resolved, curated)
    assert len(match.edges) == 1
    # the network edge runs L1 -> L2, so the curated +1.0 for L2 -> L1 enters as -1.0
    assert match.eps[0] == pytest.approx(0.25 - (-1.0))


def test_match_resolves_each_edge_inside_its_own_release():
    """The same ligand name in two releases may be two different molecules; cdk8 does this."""
    shared = canonical_smiles(BUTANOL)
    resolved = {
        "one": {name_key("L1"): canonical_smiles(ETHANOL), name_key("X"): shared},
        "two": {name_key("L1"): canonical_smiles(PROPANOL), name_key("X"): shared},
    }
    curated = [(canonical_smiles(ETHANOL), canonical_smiles(BUTANOL), 2.0)]
    net = _net([("L1", "X", 0.0, 0.1)], group="one") + _net([("L1", "X", 0.0, 0.1)], group="two")
    match = match_system(net, resolved, curated)
    assert len(match.edges) == 1, "only the release whose L1 is ethanol carries the curated edge"
    assert match.eps[0] == pytest.approx(-2.0)


def test_node_key_merges_one_molecule_named_twice_only_under_the_structure_key():
    """Two releases may give one molecule different names; cdk2, ptp1b and thrombin all do.

    The node key does not change which edges match, only how many nodes they span. Keying by
    name leaves the shared molecule as two nodes and the union in two components; keying by
    structure makes it one node and knits them together, which is what can add a cycle.
    """
    resolved = {
        "one": {name_key("A"): canonical_smiles(ETHANOL),
                name_key("B"): canonical_smiles(PROPANOL)},
        "two": {name_key("C"): canonical_smiles(PROPANOL),
                name_key("D"): canonical_smiles(BUTANOL)},
    }
    curated = [(canonical_smiles(ETHANOL), canonical_smiles(PROPANOL), 1.0),
               (canonical_smiles(PROPANOL), canonical_smiles(BUTANOL), 1.0)]
    net = _net([("A", "B", 0.0, 0.1)], group="one") + _net([("C", "D", 0.0, 0.1)], group="two")

    by_name = match_system(net, resolved, curated, by_structure=False)
    by_structure = match_system(net, resolved, curated, by_structure=True)

    assert len(by_name.edges) == len(by_structure.edges) == 2, "the edge match is unaffected"
    assert len({n for e in by_name.edges for n in e[:2]}) == 4, "B and C stay separate nodes"
    assert len({n for e in by_structure.edges for n in e[:2]}) == 3, "B and C are one node"
    assert np.allclose(by_name.eps, by_structure.eps), "the residuals are the same residuals"


def test_recovery_counts_ligand_resolution_not_edge_presence():
    """A curated edge counts as recovered when both its ligands are in the network, even if the
    network has no edge between them; the matched count is the smaller, separate quantity."""
    resolved = {"grp": {name_key("L1"): canonical_smiles(ETHANOL),
                        name_key("L2"): canonical_smiles(PROPANOL),
                        name_key("L3"): canonical_smiles(BUTANOL)}}
    curated = [
        (canonical_smiles(ETHANOL), canonical_smiles(PROPANOL), 1.0),   # present as an edge
        (canonical_smiles(ETHANOL), canonical_smiles(BUTANOL), 1.0),    # ligands present, edge not
        (canonical_smiles(ETHANOL), canonical_smiles(PENTANOL), 1.0),   # ligand absent
    ]
    net = _net([("L1", "L2", 0.0, 0.1), ("L2", "L3", 0.0, 0.1)])
    match = match_system(net, resolved, curated)
    assert match.n_recovered == 2
    assert len(match.edges) == 1
    assert match.n_unresolved_names == 0


def test_unresolved_names_are_counted_and_never_given_a_structure():
    resolved = {"grp": {name_key("L1"): canonical_smiles(ETHANOL)}}
    curated = [(canonical_smiles(ETHANOL), canonical_smiles(PROPANOL), 1.0)]
    match = match_system(_net([("L1", "L2", 0.0, 0.1)]), resolved, curated)
    assert match.n_unresolved_names == 1
    assert match.edges == []
    assert match.n_recovered == 0


# --------------------------------------------------------------------------------------
# The measurement.
# --------------------------------------------------------------------------------------
TRIANGLE = [("A", "B", 0.0, 1.0), ("B", "C", 0.0, 1.0), ("A", "C", 0.0, 1.0)]


def test_a_gradient_field_is_invisible_and_a_cycle_field_is_wholly_visible():
    potentials = {"A": 0.0, "B": 1.7, "C": -0.4}
    gradient = np.array([potentials[b] - potentials[a] for a, b, _, _ in TRIANGLE])
    f_grad, chance, _iso, dof = visible_fraction(TRIANGLE, gradient)
    assert f_grad == pytest.approx(0.0, abs=1e-12)
    assert dof == 1
    assert chance == pytest.approx(1 / 3)
    cycle = np.array([1.0, 1.0, -1.0])   # A->B->C, closed by the reversed A->C edge
    f_cycle, _c, _i, _d = visible_fraction(TRIANGLE, cycle)
    assert f_cycle == pytest.approx(1.0)


def test_the_two_chance_levels_agree_when_every_edge_has_the_same_error():
    _f, chance, iso, _dof = visible_fraction(TRIANGLE, np.array([1.0, 0.0, 0.0]))
    assert iso == pytest.approx(chance), "equal weights make tr(Pi W)/tr(W) equal dof/E"


def test_unequal_weights_separate_the_two_chance_levels():
    edges = [("A", "B", 0.0, 0.1), ("B", "C", 0.0, 1.0), ("A", "C", 0.0, 1.0)]
    _f, chance, iso, _dof = visible_fraction(edges, np.array([1.0, 0.0, 0.0]))
    assert iso != pytest.approx(chance)
    assert 0.0 < iso < 1.0


def _reading(edges, eps, system="sys", flagged=False):
    return read_system(system, Match(edges=list(edges), eps=np.asarray(eps, dtype=float),
                                     n_curated=len(edges), n_recovered=len(edges),
                                     n_unresolved_names=0), flagged)


def test_reading_reports_the_sub_network_and_flags_the_acyclic_case():
    triangle = _reading(TRIANGLE, [1.0, 1.0, -1.0])
    assert (triangle.n_nodes, triangle.n_edges, triangle.components, triangle.dof) == (3, 3, 1, 1)
    assert triangle.measurable
    path = _reading([("A", "B", 0.0, 1.0), ("B", "C", 0.0, 1.0)], [1.0, -1.0])
    assert path.dof == 0 and not path.measurable
    assert path.components == 1


def test_deletion_range_and_margin_are_taken_on_the_same_deleted_sub_network():
    """The projector moves when an edge goes, so the fraction and its chance level must be read
    off the same deletion; the worst margin is the minimum of those paired ratios."""
    edges = [("A", "B", 0.0, 1.0), ("B", "C", 0.0, 1.0), ("A", "C", 0.0, 1.0),
             ("C", "D", 0.0, 1.0), ("A", "D", 0.0, 1.0)]
    eps = np.array([0.4, -0.2, 0.1, 0.3, -0.5])
    reading = _reading(edges, eps)
    assert reading.dof == 2
    paired = []
    for k in range(len(edges)):
        sub = [e for i, e in enumerate(edges) if i != k]
        f_k, chance_k, iso_k, dof_k = visible_fraction(sub, np.delete(eps, k))
        if dof_k >= 1:
            paired.append((f_k, chance_k, iso_k))
    assert len(paired) == len(edges), "every deletion here leaves a cycle"
    assert reading.loo_visible == pytest.approx((min(p[0] for p in paired),
                                                 max(p[0] for p in paired)))
    assert reading.loo_worst_margin == pytest.approx(min(p[1] / p[0] for p in paired))
    assert reading.loo_worst_margin_iso == pytest.approx(min(p[2] / p[0] for p in paired))
    assert 0.0 <= reading.loo_visible[0] <= reading.loo_visible[1] <= 1.0


def test_pool_excludes_the_systems_that_carry_no_cycle():
    triangle = _reading(TRIANGLE, [1.0, 1.0, -1.0], system="cyclic")
    path = _reading([("A", "B", 0.0, 1.0), ("B", "C", 0.0, 1.0)], [1.0, -1.0], system="acyclic")
    pooled = pool([triangle, path])
    assert pooled["systems"] == 1
    assert pooled["E"] == 3
    assert pooled["f"] == pytest.approx(triangle.visible)
    assert np.isnan(pool([path])["f"])


def test_pool_is_a_ratio_of_sums_not_a_mean_of_ratios():
    one = _reading(TRIANGLE, [1.0, 1.0, -1.0], system="a")           # f = 1
    two = _reading(TRIANGLE, [1.0, -1.0, 2.0], system="b")           # a mixed field
    pooled = pool([one, two])
    expected = (one.cycle_sq + two.cycle_sq) / (one.total_sq + two.total_sq)
    assert pooled["f"] == pytest.approx(expected)
    assert pooled["f"] != pytest.approx(0.5 * (one.visible + two.visible))


# --------------------------------------------------------------------------------------
# Fetching and caching.
# --------------------------------------------------------------------------------------
def test_resolve_ligands_uses_the_injected_fetch_and_tolerates_a_missing_release():
    blobs = {"one": _sdf({"L1": ETHANOL}), "two": None}

    def fetch(url: str):
        return blobs["one"] if url.endswith("one/sys/ligands.sdf") else None

    resolved = resolve_ligands("sys", ["one", "two"], fetch)
    assert resolved["one"][name_key("L1")] == canonical_smiles(ETHANOL)
    assert resolved["two"] == {}, "a release that cannot be fetched contributes nothing"


def test_ligand_cache_round_trips(tmp_path):
    resolved = {"sys": {"one": {"l1": canonical_smiles(ETHANOL)}}}
    path = tmp_path / "cache.csv"
    write_ligand_cache(resolved, path)
    assert read_ligand_cache(path) == resolved
    assert read_ligand_cache(tmp_path / "absent.csv") == {}
