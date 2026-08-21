"""BAR information bottleneck — the differentiable layer (Theorems 1 & 2) and the
Fisher–resistance graph utilities (Theorem 3).

Public API is populated as Phase 1 lands. Correctness invariants:
  * variance = sandwich B/I^2, never 1/I (matches pymbar `bar(..., 'MBAR')`)
  * O(1) backward via the implicit function theorem (information-share gradient)
  * Laplacian edge weights = inverse sandwich variances I^2/B
"""
