from __future__ import annotations

from gen.env import END, LigandTrie


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
