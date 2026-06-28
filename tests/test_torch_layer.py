"""Differentiable BAR layer (Theorem 1, invariant #3): O(1) custom backward via
the implicit function theorem, dF/dx_k = p_k(1-p_k)/I (information-share gradient).
The root-find is NOT unrolled."""
from __future__ import annotations

import numpy as np
import pytest

torch = pytest.importorskip("torch")

from bar.estimator import fisher_info, solve_bar  # noqa: E402
from bar.torch_layer import bar_delta_f  # noqa: E402

XF0 = [0.0, 1.0]
XR0 = [-0.5]


def test_forward_matches_numpy_solve():
    xf = torch.tensor(XF0, dtype=torch.float64)
    xr = torch.tensor(XR0, dtype=torch.float64)
    got = bar_delta_f(xf, xr, 0.0).item()
    assert got == pytest.approx(solve_bar(np.array(XF0), np.array(XR0), 0.0), abs=1e-10)


def test_backward_forward_window_share_is_If_over_I():
    # mu_f shifts ALL forward samples; grad accumulates to sum_i p_i(1-p_i)/I = I_f/I
    base_f = torch.tensor(XF0, dtype=torch.float64)
    mu_f = torch.zeros(1, dtype=torch.float64, requires_grad=True)
    xr = torch.tensor(XR0, dtype=torch.float64)
    bar_delta_f(base_f + mu_f, xr, 0.0).backward()
    assert mu_f.grad.item() == pytest.approx(0.596825, abs=1e-5)


def test_backward_reverse_window_share_is_Ir_over_I():
    # grad wrt the UNIFIED reverse-coordinate mean = +I_r/I (theorem's -I_r/I is the
    # mu_r convention where x_j = -(mu_r+eta); here we shift x_r directly)
    base_r = torch.tensor(XR0, dtype=torch.float64)
    mu_r = torch.zeros(1, dtype=torch.float64, requires_grad=True)
    xf = torch.tensor(XF0, dtype=torch.float64)
    bar_delta_f(xf, base_r + mu_r, 0.0).backward()
    af, ar = np.array(XF0), np.array(XR0)
    _, Ir, I = fisher_info(solve_bar(af, ar, 0.0), af, ar, 0.0)
    assert mu_r.grad.item() == pytest.approx(Ir / I, abs=1e-9)


def test_per_sample_gradient_is_p_times_one_minus_p_over_I():
    xf = torch.tensor(XF0, dtype=torch.float64, requires_grad=True)
    xr = torch.tensor(XR0, dtype=torch.float64, requires_grad=True)
    bar_delta_f(xf, xr, 0.0).backward()
    d = solve_bar(np.array(XF0), np.array(XR0), 0.0)
    from scipy.special import expit

    _, _, I = fisher_info(d, np.array(XF0), np.array(XR0), 0.0)
    pf = expit(np.array(XF0) - d)
    pr = expit(np.array(XR0) - d)
    np.testing.assert_allclose(xf.grad.numpy(), pf * (1 - pf) / I, rtol=1e-9)
    np.testing.assert_allclose(xr.grad.numpy(), pr * (1 - pr) / I, rtol=1e-9)


def test_information_shares_sum_to_one_via_autograd():
    xf = torch.tensor(XF0, dtype=torch.float64, requires_grad=True)
    xr = torch.tensor(XR0, dtype=torch.float64, requires_grad=True)
    bar_delta_f(xf, xr, 0.0).backward()
    assert xf.grad.sum().item() + xr.grad.sum().item() == pytest.approx(1.0, abs=1e-12)


def test_gradcheck_double_precision():
    rng = np.random.default_rng(3)
    xf = torch.tensor(rng.normal(0.3, 1, 6), dtype=torch.float64, requires_grad=True)
    xr = torch.tensor(rng.normal(-0.3, 1, 5), dtype=torch.float64, requires_grad=True)

    def f(a, b):
        return bar_delta_f(a, b, None)

    assert torch.autograd.gradcheck(f, (xf, xr), eps=1e-6, atol=1e-6)
