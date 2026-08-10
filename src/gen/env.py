"""Evaluable generation environment for the calibrated-generation gate (Fig K).

A LigandTrie is a prefix-DAG over a fixed set of ligand SMILES: a GFlowNet builds a SMILES
char-by-char along trie edges, and terminals are exactly the input ligands (each END-terminated
so no SMILES is a prefix of another). Because the terminal set is the benchmark target's measured
ligands, every generation is evaluable. Spec:
docs/superpowers/specs/2026-06-30-gflownet-generator-design.md.
"""
from __future__ import annotations

import csv as _csv
import pathlib
from dataclasses import dataclass

import numpy as np
import torch
from numpy.typing import NDArray

from bar.calib import conformal_q
from bar.trunk import EnsembleTrunk, amortized_sigma, featurize_edge

_ROOT = pathlib.Path(__file__).resolve().parent.parent.parent
_ALL_EDGES = _ROOT / "data" / "fep_edges" / "all_edges.csv"
SIGMA_ALE = 0.4  # benchmark label noise ~0.4 kcal/mol (matches Fig J)
ALPHA = 0.10

END = "\n"


class LigandTrie:
    """Prefix-DAG over ligand SMILES. State = a prefix string (root = ""); actions append one
    character; terminals are the END-terminated input ligands."""

    def __init__(self, smiles_list: list[str]) -> None:
        self.smiles: list[str] = list(dict.fromkeys(smiles_list))
        self.children: dict[str, dict[str, str]] = {}
        self.terminal_smiles: dict[str, str] = {}
        for s in self.smiles:
            seq = s + END
            for i in range(len(seq)):
                self.children.setdefault(seq[:i], {})[seq[i]] = seq[: i + 1]
            self.children.setdefault(seq, {})
            self.terminal_smiles[seq] = s

    def actions(self, state: str) -> list[tuple[str, str]]:
        """List of (char, child_state) available from `state`."""
        return list(self.children.get(state, {}).items())

    def is_terminal(self, state: str) -> bool:
        return state in self.terminal_smiles


def _read_all_edges(csv_path: pathlib.Path) -> list[tuple[str, str, str, float]]:
    rows = []
    with open(csv_path) as f:
        for r in _csv.DictReader(f):
            rows.append((r["target"], r["smiles_a"], r["smiles_b"], float(r["ddg"])))
    return rows


def target_ligands(
    target: str, csv_path: pathlib.Path | None = None
) -> tuple[str, dict[str, float]]:
    """Reference SMILES + {ligand: true ΔΔG vs reference} for one target, from its edge graph.
    ΔΔG(b)-ΔΔG(a) = ddg on a directed edge a->b; recover per-ligand ΔΔG by a signed BFS from the
    highest-degree reference (drop ligands disconnected from it)."""
    import networkx as nx

    rows = [r for r in _read_all_edges(csv_path or _ALL_EDGES) if r[0] == target]
    g = nx.DiGraph()
    for _t, a, b, d in rows:
        g.add_edge(a, b, ddg=d)
    if g.number_of_nodes() == 0:
        raise ValueError(f"no edges for target {target}")
    ug = g.to_undirected()
    ref = max(ug.nodes, key=lambda n: ug.degree(n))
    true: dict[str, float] = {ref: 0.0}
    # sorted(): node_connected_component returns a set, and iterating a set of strings is
    # ordered by Python's randomized string hash. That order reaches the ligand list, the trie
    # and the sampler, and made this pipeline irreproducible between runs.
    for node in sorted(nx.node_connected_component(ug, ref)):
        if node == ref:
            continue
        path = nx.shortest_path(ug, ref, node)
        acc = 0.0
        for u, v in zip(path[:-1], path[1:], strict=True):
            acc += g[u][v]["ddg"] if g.has_edge(u, v) else -g[v][u]["ddg"]
        true[node] = acc
    return ref, true


@dataclass
class LeafRewards:
    true_ddg: dict[str, float]
    mu: dict[str, float]
    sigma_total: dict[str, float]
    sigma_mve: dict[str, float]


def _learned_sigma(train_edges: list[tuple[str, str]], train_ddg: NDArray,
                   eval_edges: list[tuple[str, str]], seed: int) -> NDArray:
    """Overconfident learned-σ foil: a 2-output (mean, logvar) Gaussian-NLL head on the same
    features (mirrors Fig J's MVE foil; src/gen keeps its own copy to avoid importing from figs)."""
    Xtr = np.stack([featurize_edge(a, b) for a, b in train_edges]).astype(np.float32)
    Xev = np.stack([featurize_edge(a, b) for a, b in eval_edges]).astype(np.float32)
    mu, sd = Xtr.mean(0), Xtr.std(0) + 1e-6
    torch.manual_seed(seed)
    net = torch.nn.Sequential(torch.nn.Linear(Xtr.shape[1], 64), torch.nn.SiLU(),
                              torch.nn.Linear(64, 2))
    opt = torch.optim.Adam(net.parameters(), lr=3e-3)
    xt = torch.tensor((Xtr - mu) / sd)
    yt = torch.tensor(np.asarray(train_ddg, dtype=np.float32))
    for _ in range(400):
        opt.zero_grad()
        o = net(xt)
        nll = (0.5 * o[:, 1] + 0.5 * (yt - o[:, 0]) ** 2 / torch.exp(o[:, 1])).mean()
        nll.backward()
        opt.step()
    with torch.no_grad():
        oe = net(torch.tensor((Xev - mu) / sd)).numpy()
    return np.sqrt(np.exp(oe[:, 1]))


def precompute_leaf_rewards(target: str, seed: int = 0, n_members: int = 8,
                            csv_path: pathlib.Path | None = None) -> tuple[str, LeafRewards]:
    """Train the trunk on all_edges MINUS `target` (target-disjoint -> the target's ligands are
    partially OOD), conformal-calibrate in-distribution, and fill per-ligand rewards."""
    path = csv_path or _ALL_EDGES
    ref, true = target_ligands(target, path)
    rng = np.random.default_rng(1000 + seed)
    other = [(a, b, d) for (t, a, b, d) in _read_all_edges(path) if t != target]
    rng.shuffle(other)
    cut = int(0.8 * len(other))
    fit_edges = [(a, b) for a, b, _ in other[:cut]]
    fit_ddg = np.array([d for _a, _b, d in other[:cut]])
    cal_edges = [(a, b) for a, b, _ in other[cut:]]
    cal_ddg = np.array([d for _a, _b, d in other[cut:]])

    trunk = EnsembleTrunk().fit(fit_edges, fit_ddg, n_members=n_members)
    mu_cal, se_cal = trunk.predict(cal_edges)
    st_cal = np.sqrt(se_cal ** 2 + SIGMA_ALE ** 2)
    q = conformal_q(np.abs(cal_ddg - mu_cal), st_cal, ALPHA)

    ligs = list(true)
    edges = [(ref, L) for L in ligs]
    mu_l, se_l = trunk.predict(edges)
    sig_total = amortized_sigma(se_l, SIGMA_ALE, q)
    sig_mve = _learned_sigma(fit_edges, fit_ddg, edges, seed)
    lr = LeafRewards(
        true_ddg=true,
        mu={L: float(mu_l[i]) for i, L in enumerate(ligs)},
        sigma_total={L: float(sig_total[i]) for i, L in enumerate(ligs)},
        sigma_mve={L: float(sig_mve[i]) for i, L in enumerate(ligs)},
    )
    return ref, lr
