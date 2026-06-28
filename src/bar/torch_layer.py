"""Differentiable BAR layer (Theorem 1, invariant #3).

A ``torch.autograd.Function`` whose forward solves ``S(dF)=0`` (delegated to the
numpy root-finder) and whose backward is **O(1)** via the implicit function theorem:
``dF/dx_k = (1/I) d_{x_k} S = p_k(1-p_k)/I`` for every unified-coordinate sample k.
The root-find is NOT unrolled — backward is one division by the curvature I.
"""
from __future__ import annotations

import numpy as np
import torch

from .estimator import default_log_count_ratio, solve_bar


class _BARDeltaF(torch.autograd.Function):
    @staticmethod
    def forward(ctx, x_f: torch.Tensor, x_r: torch.Tensor, M: float) -> torch.Tensor:
        xf = x_f.detach().cpu().numpy().astype(float)
        xr = x_r.detach().cpu().numpy().astype(float)
        d = solve_bar(xf, xr, M)
        # p(1-p) per sample at the root; saved for the O(1) backward.
        pf = 1.0 / (1.0 + np.exp(-(xf - d + M)))
        pr = 1.0 / (1.0 + np.exp(-(xr - d + M)))
        gf = pf * (1.0 - pf)
        gr = pr * (1.0 - pr)
        inv_I = 1.0 / (gf.sum() + gr.sum())
        ctx.save_for_backward(
            torch.as_tensor(gf * inv_I, dtype=x_f.dtype, device=x_f.device),
            torch.as_tensor(gr * inv_I, dtype=x_r.dtype, device=x_r.device),
        )
        return torch.as_tensor(d, dtype=x_f.dtype, device=x_f.device)

    @staticmethod
    def backward(ctx, grad_out: torch.Tensor):
        grad_f, grad_r = ctx.saved_tensors
        return grad_out * grad_f, grad_out * grad_r, None


def bar_delta_f(
    x_f: torch.Tensor, x_r: torch.Tensor, M: float | None = None
) -> torch.Tensor:
    """Differentiable BAR free-energy estimate dF_hat. ``M`` defaults to
    ln(n_f/n_r). Gradients w.r.t. x_f, x_r are the information-share contributions
    p(1-p)/I; summed over a window they give I_window/I (Theorem 1)."""
    m = default_log_count_ratio(x_f.shape[-1], x_r.shape[-1]) if M is None else float(M)
    return _BARDeltaF.apply(x_f, x_r, m)
