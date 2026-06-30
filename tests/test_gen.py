from __future__ import annotations

import math
from collections import Counter

import numpy as np

from gen.env import END, LeafRewards, LigandTrie, precompute_leaf_rewards, target_ligands
from gen.gflownet import TabularTB


def test_trie_reaches_exactly_the_input_ligands():
    trie = LigandTrie(["CCO", "CCN", "CCO"])  # dup ignored
    assert trie.smiles == ["CCO", "CCN"]
    # terminals are exactly the END-terminated inputs, nothing else
    assert set(trie.terminal_smiles.values()) == {"CCO", "CCN"}
    assert all(t.endswith(END) for t in trie.terminal_smiles)
    # root branches into the shared first char 'C'
    assert [a[0] for a in trie.actions("")] == ["C"]
    # a full walk reaches a terminal
    state = ""
    while not trie.is_terminal(state):
        state = trie.actions(state)[0][1]
    assert trie.terminal_smiles[state] in {"CCO", "CCN"}


def test_no_smiles_is_a_prefix_terminal_of_another():
    # 'CCO' is a prefix of 'CCOC'; END-termination keeps them distinct terminals
    trie = LigandTrie(["CCO", "CCOC"])
    assert not trie.is_terminal("CCO")        # not a terminal (no END yet)
    assert trie.is_terminal("CCO" + END)
    assert trie.is_terminal("CCOC" + END)


# jnk1 is not present in data/fep_edges/all_edges.csv; p38 is the chosen substitute
# (29 nodes, 29-node largest connected component — well above the >=8 threshold).
_TARGET = "p38"


def test_target_ligands_recovers_ddg_graph():
    ref, true = target_ligands(_TARGET)
    assert ref in true and abs(true[ref]) < 1e-9   # reference has ΔΔG 0 vs itself
    assert len(true) >= 8                           # a real congeneric series
    assert all(np.isfinite(v) for v in true.values())


def test_precompute_leaf_rewards_fills_all_arms():
    ref, lr = precompute_leaf_rewards(_TARGET, seed=0, n_members=4)
    assert isinstance(lr, LeafRewards)
    ligs = set(lr.true_ddg)
    assert ligs == set(lr.mu) == set(lr.sigma_total) == set(lr.sigma_mve)
    assert all(s >= 0 for s in lr.sigma_total.values())
    assert all(s >= 0 for s in lr.sigma_mve.values())
    # the trunk's sigma_total varies across ligands (some OOD) — not a constant
    assert np.std(list(lr.sigma_total.values())) > 1e-6


def test_tb_samples_proportional_to_reward():
    # toy trie with known rewards: the trained sampler must match the Boltzmann distribution
    trie = LigandTrie(["AA", "AB", "BA"])
    log_reward = {"AA": math.log(1.0), "AB": math.log(3.0), "BA": math.log(2.0)}
    tb = TabularTB(trie, seed=0).train(log_reward, steps=1500, lr=0.1, batch=16, seed=0)
    c = Counter(tb.sample_terminals(3000, seed=1))
    freq = {k: c[k] / 3000 for k in log_reward}
    target = {"AA": 1 / 6, "AB": 3 / 6, "BA": 2 / 6}  # ∝ reward
    for k in target:
        assert abs(freq[k] - target[k]) < 0.06


def test_tb_sampler_is_diverse_not_collapsed():
    trie = LigandTrie(["AA", "AB", "BA", "BB"])
    log_reward = dict.fromkeys(["AA", "AB", "BA", "BB"], 0.0)  # uniform reward
    tb = TabularTB(trie, seed=0).train(log_reward, steps=800, lr=0.1, batch=16, seed=0)
    c = Counter(tb.sample_terminals(2000, seed=2))
    assert len(c) == 4  # covers all modes, no collapse
