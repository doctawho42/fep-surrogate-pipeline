"""Curated per-edge experimental labels for the OpenFE replicate networks.

The visible-fraction measurement of Theorem 5 needs, per edge, a calculated ``ddG`` with a
reported standard error and an experimental ``ddG`` to compare it against. The article's headline
takes the experimental side from a stereo-blind, target-wide ChEMBL search, which is available for
four systems. The OpenFF protein-ligand benchmark (``data/fep_edges/*.csv``) carries a *curated*
per-edge experimental ``ddG`` for fifteen targets, fourteen of which name a system of the OpenFE
replicate benchmark verbatim. This module joins the two so the same measurement can be run on all
fourteen.

The join rule is the one the Supporting Information already used for its three-system cross-check,
recovered from that check's reported numbers and reproduced here to the digits quoted
(``hif2a`` 98 %, ``p38`` 91 %; see ``docs/results_figGround.md`` for the one deviation). It is:

1. a network ligand name is resolved to a structure through the OpenFE IndustryBenchmarks2024
   ``prepared_structures/<system group>/<system name>/ligands.sdf`` of the *release the edge
   belongs to*, so two releases that reuse a ligand name for different molecules -- which
   ``cdk8`` does, seven times -- cannot be mixed;
2. names are compared separator-insensitively: runs of whitespace, ``-`` and ``_`` collapse to one
   space, which is what turns the SDF's ``30 flip``, ``Example 7`` and ``CHEMBL3402756_2.7
   redocked`` into the benchmark's ``30-flip``, ``Example-7`` and ``CHEMBL3402756_2.7_redocked``.
   Nothing else about a name is rewritten and no unresolved name is ever given a structure;
3. structures are compared as RDKit isomeric canonical SMILES on both sides;
4. a curated edge matches a network edge when the two unordered canonical-SMILES pairs are equal;
   the curated ``ddG`` is then oriented to the network edge's direction.

Two counts follow from the rule and mean different things. The *recovery rate* is the share of a
system's curated edges whose two endpoints both resolve to a structure in the network, which is
what the Supporting Information quotes; the *matched* count is how many of those edges the network
actually carries, and is the size of the sub-network the measurement runs on. The second is always
the smaller, and it is the one that decides whether a system carries a measurement at all.

Pure NumPy plus RDKit for the structure keys. The one network access -- fetching the benchmark's
input-structure SDFs -- is injectable and cached to ``data/curated/openfe_ligand_smiles.csv``, so a
run after the first is offline and deterministic.
"""
from __future__ import annotations

import csv
import math
import re
import urllib.error
import urllib.request
from collections.abc import Callable, Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from numpy.typing import NDArray
from rdkit import Chem, RDLogger
from rdkit.Chem.MolStandardize import rdMolStandardize

from bar.hodge import hodge_split
from bar.qc import Edge, _incidence

RDLogger.DisableLog("rdApp.*")

COMBINED_CSV = Path("data/openfe_replicates/combined_pymbar4_edge_data.csv")
CURATED_DIR = Path("data/fep_edges")
LIGAND_CACHE = Path("data/curated/openfe_ligand_smiles.csv")

_PREPARED = (
    "https://raw.githubusercontent.com/OpenFreeEnergy/IndustryBenchmarks2024/main/"
    "industry_benchmarks/input_structures/prepared_structures"
)
_SEPARATORS = re.compile(r"[\s_\-]+")


# --------------------------------------------------------------------------------------
# 1. The join keys.
# --------------------------------------------------------------------------------------
def name_key(name: str) -> str:
    """Separator-insensitive ligand-name key: runs of whitespace, ``-`` and ``_`` become one space.

    The benchmark and its own input SDFs disagree only on separators (``30 flip`` against
    ``30-flip``, ``CHEMBL3402756_2.7 redocked`` against ``CHEMBL3402756_2.7_redocked``). Nothing
    else about a name is rewritten, so a name that still fails to resolve stays unresolved.
    """
    return _SEPARATORS.sub(" ", name).strip()


def canonical_smiles(smiles: str, neutralize: bool = False) -> str | None:
    """RDKit isomeric canonical SMILES, or ``None`` if the string does not parse.

    ``neutralize`` additionally strips formal charges through RDKit's standardizer before
    canonicalizing. That is NOT the frozen join rule: it is the key of the labelled sensitivity
    analysis, which exists because the benchmark's prepared input structures are protonated where
    its curated table is neutral -- ``thrombin`` is entirely protonated, and under the exact rule
    not one of its curated edges joins. Both readings are reported; neither replaces the other.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    if neutralize:
        mol = rdMolStandardize.Uncharger().uncharge(mol)
    return str(Chem.MolToSmiles(mol))


def sdf_smiles(text: str, neutralize: bool = False) -> dict[str, str]:
    """Parse an SDF blob into ``name_key -> canonical SMILES``; unparseable records are dropped."""
    out: dict[str, str] = {}
    supplier = Chem.SDMolSupplier()
    supplier.SetData(text, removeHs=False)
    for mol in supplier:
        if mol is None or not mol.HasProp("_Name"):
            continue
        key = name_key(str(mol.GetProp("_Name")))
        canonical = canonical_smiles(str(Chem.MolToSmiles(mol)), neutralize=neutralize)
        if key and canonical is not None:
            out[key] = canonical
    return out


# --------------------------------------------------------------------------------------
# 2. The two edge sources.
# --------------------------------------------------------------------------------------
@dataclass(frozen=True)
class NetworkEdge:
    """One replicate-averaged calculated edge of a system's network, with its release."""

    a: str
    b: str
    ddg: float
    se: float
    group: str


def _float(value: str | None) -> float:
    try:
        return float(value)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return math.nan


def network_edges(system: str, csv_path: Path = COMBINED_CSV) -> list[NetworkEdge]:
    """A system's edges, averaged over the replicates exactly as ``closeloop.load_system_edges``.

    The mean of the usable replicates' ``ddG`` and ``se = sqrt(mean se_k^2 / n)``, keeping edges
    with at least two usable replicates. Carries the release each edge came from, which the name
    resolution needs.
    """
    edges: list[NetworkEdge] = []
    with csv_path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            if row["system name"] != system:
                continue
            values, errors = [], []
            for k in (0, 1, 2):
                complex_dg = _float(row[f"complex_repeat_{k}_DG (kcal/mol)"])
                complex_se = _float(row[f"complex_repeat_{k}_dDG (kcal/mol)"])
                solvent_dg = _float(row[f"solvent_repeat_{k}_DG (kcal/mol)"])
                solvent_se = _float(row[f"solvent_repeat_{k}_dDG (kcal/mol)"])
                quad = (complex_dg, complex_se, solvent_dg, solvent_se)
                if not any(math.isnan(v) for v in quad):
                    values.append(complex_dg - solvent_dg)
                    errors.append(math.sqrt(complex_se ** 2 + solvent_se ** 2))
            if len(values) >= 2:
                edges.append(NetworkEdge(
                    a=row["ligand_A"], b=row["ligand_B"], ddg=float(np.mean(values)),
                    se=math.sqrt(float(np.mean(np.array(errors) ** 2)) / len(values)),
                    group=row["system group"],
                ))
    return edges


def system_groups(system: str, csv_path: Path = COMBINED_CSV) -> list[str]:
    """The releases (``system group``) a system's edges come from, in sorted order."""
    return sorted({e.group for e in network_edges(system, csv_path)})


def curated_edges(system: str, directory: Path = CURATED_DIR,
                  neutralize: bool = False) -> list[tuple[str, str, float]]:
    """The benchmark's curated ``(canonical A, canonical B, ddG)`` edges for a target.

    Rows whose SMILES do not parse, and rows whose two endpoints canonicalize to the same
    structure, are dropped: neither can be joined to a network edge.
    """
    out: list[tuple[str, str, float]] = []
    with (directory / f"{system}.csv").open(newline="") as handle:
        for row in csv.DictReader(handle):
            first = canonical_smiles(row["smiles_a"], neutralize=neutralize)
            second = canonical_smiles(row["smiles_b"], neutralize=neutralize)
            if first is None or second is None or first == second:
                continue
            out.append((first, second, float(row["ddg"])))
    return out


# --------------------------------------------------------------------------------------
# 3. Resolving network ligand names to structures.
# --------------------------------------------------------------------------------------
def _fetch(url: str) -> str | None:
    try:
        with urllib.request.urlopen(url, timeout=60) as response:
            return str(response.read().decode("utf-8", errors="replace"))
    except (urllib.error.URLError, TimeoutError, OSError):
        return None


def resolve_ligands(
    system: str,
    groups: Sequence[str],
    fetch: Callable[[str], str | None] = _fetch,
    neutralize: bool = False,
) -> dict[str, dict[str, str]]:
    """``release -> {name_key -> canonical SMILES}`` from the benchmark's prepared input SDFs.

    One SDF per release, at the path the benchmark's own directory layout dictates. A release
    whose SDF cannot be fetched contributes nothing rather than a guess.
    """
    out: dict[str, dict[str, str]] = {}
    for group in groups:
        text = fetch(f"{_PREPARED}/{group}/{system}/ligands.sdf")
        out[group] = sdf_smiles(text, neutralize=neutralize) if text is not None else {}
    return out


def write_ligand_cache(
    resolved: Mapping[str, Mapping[str, Mapping[str, str]]],
    path: Path = LIGAND_CACHE,
) -> None:
    """Write ``{system: {release: {name_key: SMILES}}}`` to a committed CSV."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["system", "group", "ligand_key", "smiles"])
        for system in sorted(resolved):
            for group in sorted(resolved[system]):
                for key in sorted(resolved[system][group]):
                    writer.writerow([system, group, key, resolved[system][group][key]])


def read_ligand_cache(path: Path = LIGAND_CACHE) -> dict[str, dict[str, dict[str, str]]]:
    """Read the committed ligand-structure cache; an empty mapping if it does not exist."""
    out: dict[str, dict[str, dict[str, str]]] = {}
    if not path.exists():
        return out
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            out.setdefault(row["system"], {}).setdefault(row["group"], {})
            out[row["system"]][row["group"]][row["ligand_key"]] = row["smiles"]
    return out


# --------------------------------------------------------------------------------------
# 4. The match.
# --------------------------------------------------------------------------------------
@dataclass(frozen=True)
class Match:
    """A system's curated labels joined onto its network."""

    edges: list[Edge]            # (a, b, calculated ddG, se) of the matched sub-network
    eps: NDArray                 # calculated minus curated, per matched edge, kcal/mol
    n_curated: int               # curated edges offered
    n_recovered: int             # ... whose two endpoints both resolve into the network
    n_unresolved_names: int      # network ligand names with no structure in their release


def match_system(
    net: Sequence[NetworkEdge],
    resolved: Mapping[str, Mapping[str, str]],
    curated: Sequence[tuple[str, str, float]],
    by_structure: bool = True,
) -> Match:
    """Join curated edges onto a network, by unordered canonical-SMILES pair.

    An edge of the network is kept when both its ligands resolve to a structure in their own
    release and the curated set contains that unordered structure pair; the curated value is then
    signed to the network edge's direction. Where the curated set offers the same structure pair
    twice, the first row wins, deterministically.

    ``by_structure`` decides how the matched sub-network's *nodes* are identified, which is a
    separate question from how its edges are matched and changes the topology rather than the
    membership. Under the default, two ligand names that resolve to the same molecule are one
    node, so releases that name a shared ligand differently are knitted into one component; under
    ``False`` they stay two nodes. This benchmark reuses one name for different molecules
    (``cdk8``, seven times) and gives one molecule different names in different releases
    (``cdk2``, ``ptp1b``, ``thrombin``), so a name is not an identifier in it and the structure
    key is the supported one. It is also the less favourable of the two here: merging nodes adds
    cycles, which can only raise the auditable share. Both readings are reported.
    """
    by_pair: dict[frozenset[str], tuple[str, str, float]] = {}
    for first, second, ddg in curated:
        by_pair.setdefault(frozenset((first, second)), (first, second, ddg))

    structures: set[str] = set()
    unresolved: set[tuple[str, str]] = set()
    for edge in net:
        for name in (edge.a, edge.b):
            smiles = resolved.get(edge.group, {}).get(name_key(name))
            if smiles is None:
                unresolved.add((edge.group, name))
            else:
                structures.add(smiles)
    n_recovered = sum(1 for a, b, _ in curated if a in structures and b in structures)

    edges: list[Edge] = []
    eps: list[float] = []
    for edge in net:
        head = resolved.get(edge.group, {}).get(name_key(edge.a))
        tail = resolved.get(edge.group, {}).get(name_key(edge.b))
        if head is None or tail is None or head == tail:
            continue
        hit = by_pair.get(frozenset((head, tail)))
        if hit is None:
            continue
        curated_head, _curated_tail, ddg = hit
        experimental = ddg if curated_head == head else -ddg
        key_a, key_b = (head, tail) if by_structure else (edge.a, edge.b)
        edges.append((key_a, key_b, edge.ddg, edge.se))
        eps.append(edge.ddg - experimental)
    return Match(edges=edges, eps=np.array(eps, dtype=float), n_curated=len(curated),
                 n_recovered=n_recovered, n_unresolved_names=len(unresolved))


# --------------------------------------------------------------------------------------
# 5. The measurement.
# --------------------------------------------------------------------------------------
def projector(edges: Sequence[Edge]) -> tuple[NDArray, NDArray]:
    """``(residual projector of the whitened incidence, weight matrix W = diag(1/se^2))``."""
    se = np.array([e[3] for e in edges], dtype=float)
    whitened = np.diag(1.0 / se) @ _incidence(list(edges))[1]
    proj = np.eye(len(edges)) - whitened @ np.linalg.pinv(whitened.T @ whitened) @ whitened.T
    return proj, np.diag(1.0 / se ** 2)


def visible_fraction(edges: Sequence[Edge], eps: NDArray) -> tuple[float, float, float, int]:
    """``(f, dof/E chance, isotropic-in-kcal/mol chance, dof)`` for one error field.

    ``f = ||Pi eps~||^2 / ||eps~||^2`` in the whitened metric the closure statistic lives in;
    the isotropic chance level is ``tr(Pi W)/tr(W)``, what an error isotropic in kcal/mol would
    show under the same projector.
    """
    edges = list(edges)
    se = np.array([e[3] for e in edges], dtype=float)
    proj, weights = projector(edges)
    tilde = np.asarray(eps, dtype=float) / se
    total = float(tilde @ tilde)
    split = hodge_split(edges, np.asarray(eps, dtype=float))
    fraction = float(tilde @ proj @ tilde) / total if total > 0 else 0.0
    chance_iso = float(np.trace(proj @ weights) / np.trace(weights))
    return fraction, split.dof / len(edges), chance_iso, split.dof


@dataclass(frozen=True)
class Reading:
    """Everything one system contributes: its match, its sub-network, and its measurement."""

    system: str
    n_curated: int
    n_recovered: int
    recovery: float
    n_unresolved_names: int
    n_nodes: int
    n_edges: int
    components: int
    dof: int
    visible: float
    chance: float
    chance_iso: float
    loo_visible: tuple[float, float]
    loo_chance: tuple[float, float]
    loo_chance_iso: tuple[float, float]
    loo_worst_margin: float
    loo_worst_margin_iso: float
    median_abs_eps: float
    median_se: float
    cycle_sq: float
    total_sq: float
    trace_proj_w: float
    trace_w: float
    flagged: bool

    @property
    def measurable(self) -> bool:
        """A sub-network with no independent cycle carries no visible fraction at all."""
        return self.dof >= 1


def read_system(system: str, match: Match, flagged: bool) -> Reading:
    """One system's match as a reading; a sub-network too small to fit is reported, not hidden."""
    edges, eps = match.edges, match.eps
    empty = (float("nan"), float("nan"))
    base = {
        "system": system, "n_curated": match.n_curated, "n_recovered": match.n_recovered,
        "recovery": match.n_recovered / match.n_curated if match.n_curated else float("nan"),
        "n_unresolved_names": match.n_unresolved_names,
    }
    if len(edges) < 2:
        return Reading(**base, n_nodes=len({n for e in edges for n in e[:2]}),  # type: ignore[arg-type]
                       n_edges=len(edges), components=0, dof=0, visible=float("nan"),
                       chance=float("nan"), chance_iso=float("nan"), loo_visible=empty,
                       loo_chance=empty, loo_chance_iso=empty,
                       loo_worst_margin=float("nan"), loo_worst_margin_iso=float("nan"),
                       median_abs_eps=float(np.median(np.abs(eps))) if len(eps) else float("nan"),
                       median_se=float("nan"), cycle_sq=0.0, total_sq=0.0,
                       trace_proj_w=0.0, trace_w=0.0, flagged=flagged)

    se = np.array([e[3] for e in edges], dtype=float)
    fraction, chance, chance_iso, dof = visible_fraction(edges, eps)
    proj, weights = projector(edges)
    tilde = eps / se
    nodes = _incidence(edges)[0]
    rank = int(np.linalg.matrix_rank(_incidence(edges)[1]))

    loo: list[tuple[float, float, float]] = []
    for k in range(len(edges)):
        sub = [e for i, e in enumerate(edges) if i != k]
        sub_eps = np.array([eps[i] for i in range(len(edges)) if i != k])
        if len(sub) < 2:
            continue
        f_k, chance_k, iso_k, dof_k = visible_fraction(sub, sub_eps)
        if dof_k >= 1:
            loo.append((f_k, chance_k, iso_k))
    ranges = tuple(
        (min(x[i] for x in loo), max(x[i] for x in loo)) if loo else empty for i in range(3)
    )
    # the margin is paired per deletion: deleting an edge moves the projector, so the fraction
    # and its chance level must be read off the same sub-network. Deletions that leave no visible
    # component at all carry no margin and are skipped rather than reported as infinite.
    positive = [x for x in loo if x[0] > 0]
    worst = tuple(min(x[i] / x[0] for x in positive) if positive else float("nan")
                  for i in (1, 2))
    return Reading(**base, n_nodes=len(nodes), n_edges=len(edges), components=len(nodes) - rank,  # type: ignore[arg-type]
                   dof=dof, visible=fraction, chance=chance, chance_iso=chance_iso,
                   loo_visible=ranges[0], loo_chance=ranges[1], loo_chance_iso=ranges[2],
                   loo_worst_margin=worst[0], loo_worst_margin_iso=worst[1],
                   median_abs_eps=float(np.median(np.abs(eps))), median_se=float(np.median(se)),
                   cycle_sq=float(tilde @ proj @ tilde), total_sq=float(tilde @ tilde),
                   trace_proj_w=float(np.trace(proj @ weights)), trace_w=float(np.trace(weights)),
                   flagged=flagged)


def pool(readings: Iterable[Reading]) -> dict[str, float]:
    """Pool readings the way the article pools: ratios of summed squared norms, not means.

    ``f`` is ``sum ||Pi eps~||^2 / sum ||eps~||^2``, the chance level is ``sum dof / sum E``, and
    the isotropic chance level is ``sum tr(Pi W) / sum tr(W)``. Systems whose matched sub-network
    carries no cycle contribute nothing to either side and are excluded.
    """
    usable = [r for r in readings if r.measurable]
    if not usable:
        return {"f": float("nan"), "chance": float("nan"), "chance_iso": float("nan"),
                "E": 0.0, "dof": 0.0, "systems": 0.0}
    total = sum(r.total_sq for r in usable)
    return {
        "f": sum(r.cycle_sq for r in usable) / total if total > 0 else float("nan"),
        "chance": sum(r.dof for r in usable) / sum(r.n_edges for r in usable),
        "chance_iso": sum(r.trace_proj_w for r in usable) / sum(r.trace_w for r in usable),
        "E": float(sum(r.n_edges for r in usable)),
        "dof": float(sum(r.dof for r in usable)),
        "systems": float(len(usable)),
    }
