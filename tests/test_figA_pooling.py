import numpy as np
from figs.make_figA import pooled_se_recovery


def test_pooled_recovery_shrinks_and_is_small_at_large_N():
    e200 = pooled_se_recovery(200, seed=7)
    e4000 = pooled_se_recovery(4000, seed=7)
    assert 0.0 < e4000 < 0.10          # recovers sandwich to <10% at 4000 single-label edges
    assert e4000 < e200                # error shrinks with budget (identifiable)


def test_pooled_recovery_deterministic():
    assert pooled_se_recovery(500, seed=3) == pooled_se_recovery(500, seed=3)
