import numpy as np

from rf3bp_lab.dynamics.params import SystemParams
from rf3bp_lab.shooting.hierarchical import HierarchicalShooter


def test_hierarchical_shooter_runs() -> None:
    p = SystemParams()
    shooter = HierarchicalShooter(p)
    seed = np.array([0.55, 0.0, 0.0, 0.0, 0.42, 0.02], dtype=float)
    out = shooter.solve(seed)
    assert "period" in out
    assert out["period"] > 0.0
    assert np.isfinite(out["residual_norm"])
    assert "stages" in out
    assert len(out["stages"]) > 2
    assert out["stages"][0].label == "cr3bp"
