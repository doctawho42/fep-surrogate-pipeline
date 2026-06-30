# src/gen/gflownet.py
"""Minimal tabular trajectory-balance GFlowNet over a LigandTrie. One learnable forward logit
per (state, action) edge + a scalar log-partition logZ; the trajectory-balance loss
(logZ + Σ log P_F − log R)² trains the forward policy so terminals are sampled ∝ reward
(Bengio et al. 2021). Spec: docs/superpowers/specs/2026-06-30-gflownet-generator-design.md.
"""
from __future__ import annotations

import torch

from gen.env import LigandTrie


class TabularTB:
    def __init__(self, trie: LigandTrie, seed: int = 0) -> None:
        self.trie = trie
        self.edges: list[tuple[str, str]] = [
            (st, ch) for st in trie.children for ch in trie.children[st]
        ]
        self.idx: dict[tuple[str, str], int] = {e: i for i, e in enumerate(self.edges)}
        torch.manual_seed(seed)
        self.logit = torch.zeros(max(len(self.edges), 1), requires_grad=True)
        self.logZ = torch.zeros(1, requires_grad=True)

    def _logp(self, state: str) -> tuple[list[tuple[str, str]], torch.Tensor]:
        acts = self.trie.actions(state)
        ids = [self.idx[(state, ch)] for ch, _ in acts]
        return acts, torch.log_softmax(self.logit[ids], dim=0)

    def _rollout(self, gen: torch.Generator) -> tuple[str, torch.Tensor]:
        state, logpf = "", torch.zeros(())
        while not self.trie.is_terminal(state):
            acts = self.trie.actions(state)
            if len(acts) == 1:
                state = acts[0][1]           # deterministic: forward logp = 0
                continue
            _, logp = self._logp(state)      # _logp recomputes acts internally; fine
            i = int(torch.multinomial(torch.exp(logp.detach()), 1, generator=gen).item())
            logpf = logpf + logp[i]
            state = acts[i][1]
        return state, logpf

    def train(self, log_reward: dict[str, float], steps: int = 2000, lr: float = 0.1,
              batch: int = 16, seed: int = 0) -> TabularTB:
        gen = torch.Generator().manual_seed(seed)
        opt = torch.optim.Adam([self.logit, self.logZ], lr=lr)
        for _ in range(steps):
            opt.zero_grad()
            loss = torch.zeros(())
            for _ in range(batch):
                term, logpf = self._rollout(gen)
                log_rew = log_reward[self.trie.terminal_smiles[term]]
                loss = loss + (self.logZ + logpf - log_rew) ** 2
            (loss / batch).backward()
            opt.step()
        return self

    def sample_terminals(self, k: int, seed: int = 0) -> list[str]:
        gen = torch.Generator().manual_seed(seed)
        with torch.no_grad():
            return [self.trie.terminal_smiles[self._rollout(gen)[0]] for _ in range(k)]
