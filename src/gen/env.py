"""Evaluable generation environment for the calibrated-generation gate (Fig K).

A LigandTrie is a prefix-DAG over a fixed set of ligand SMILES: a GFlowNet builds a SMILES
char-by-char along trie edges, and terminals are exactly the input ligands (each END-terminated
so no SMILES is a prefix of another). Because the terminal set is the benchmark target's measured
ligands, every generation is evaluable. Spec:
docs/superpowers/specs/2026-06-30-gflownet-generator-design.md.
"""
from __future__ import annotations

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
