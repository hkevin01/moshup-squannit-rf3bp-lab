import numpy as np

from rf3bp_lab.dynamics.models import cr3bp_rhs, rf3bp_pulsating_rhs
from rf3bp_lab.dynamics.params import SystemParams


def test_rhs_shapes() -> None:
    p = SystemParams()
    state = np.array([0.4, 0.1, 0.0, 0.0, 0.3, 0.0], dtype=float)
    a = cr3bp_rhs(0.0, state, p)
    b = rf3bp_pulsating_rhs(0.0, state, p)
    assert a.shape == (6,)
    assert b.shape == (6,)
    assert np.all(np.isfinite(a))
    assert np.all(np.isfinite(b))
