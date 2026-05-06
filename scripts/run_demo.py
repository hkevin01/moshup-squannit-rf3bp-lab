from __future__ import annotations

import os
import numpy as np

from rf3bp_lab.dynamics.models import cr3bp_rhs, propagate, rf3bp_pulsating_rhs
from rf3bp_lab.dynamics.params import SystemParams
from rf3bp_lab.shooting.hierarchical import HierarchicalShooter
from rf3bp_lab.utils.plotting import plot_trajectory


def main() -> None:
    p = SystemParams()
    seed = np.array([0.55, 0.0, 0.0, 0.0, 0.42, 0.02], dtype=float)

    shooter = HierarchicalShooter(p)
    result = shooter.solve(seed)

    s0 = result["state0"]
    period = float(result["period"])

    t_eval = np.linspace(0.0, period * 4.0, 3000)
    sol_cr3bp = propagate(cr3bp_rhs, s0, (0.0, t_eval[-1]), p, t_eval=t_eval)
    sol_rf3bp = propagate(rf3bp_pulsating_rhs, s0, (0.0, t_eval[-1]), p, t_eval=t_eval)

    print("Hierarchical shooting summary")
    print(f"period={period:.6f}")
    print(f"residual_norm={result['residual_norm']:.3e}")
    print(f"stage1_cost={result['stage1_cost']:.3e}, stage2_cost={result['stage2_cost']:.3e}")

    plot_trajectory(sol_cr3bp.y[:3], "CR3BP Seed Orbit")
    plot_trajectory(sol_rf3bp.y[:3], "RF3BP Pulsating + Perturbations")

    import matplotlib.pyplot as plt

    out_dir = os.environ.get("OUTPUT_DIR", ".")
    os.makedirs(out_dir, exist_ok=True)
    out_file = os.path.join(out_dir, "rf3bp_demo.png")

    if os.environ.get("SHOW_PLOTS", "0") == "1":
        plt.show()
    else:
        plt.savefig(out_file, dpi=150)
        print(f"saved plot: {out_file}")


if __name__ == "__main__":
    main()
