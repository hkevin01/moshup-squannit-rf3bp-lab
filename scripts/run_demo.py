from __future__ import annotations

import os

import matplotlib.pyplot as plt
import numpy as np

from rf3bp_lab.dynamics.models import FidelityWeights, cr3bp_rhs, propagate, rf3bp_breakdown, rf3bp_pulsating_rhs
from rf3bp_lab.dynamics.params import SystemParams
from rf3bp_lab.shooting.hierarchical import HierarchicalShooter
from rf3bp_lab.utils.plotting import plot_perturbation_norms, plot_stage_convergence, plot_trajectory


def _component_norm_history(t: np.ndarray, states: np.ndarray, p: SystemParams) -> dict[str, np.ndarray]:
    keys = ["nonspherical", "pulsation", "solar_gravity", "srp"]
    history = {k: np.zeros_like(t) for k in keys}
    for i in range(t.size):
        b_full = rf3bp_breakdown(float(t[i]), states[:, i], p)
        b_no_j2 = rf3bp_breakdown(
            float(t[i]),
            states[:, i],
            p,
            FidelityWeights(pulsation=1.0, nonspherical=0.0, solar_gravity=1.0, srp=1.0),
        )
        history["nonspherical"][i] = np.linalg.norm((b_full.grav_primary + b_full.grav_secondary) - (b_no_j2.grav_primary + b_no_j2.grav_secondary))
        history["pulsation"][i] = np.linalg.norm(b_full.pulsation)
        history["solar_gravity"][i] = np.linalg.norm(b_full.solar_gravity)
        history["srp"][i] = np.linalg.norm(b_full.srp)
    return history


def main() -> None:
    p = SystemParams()
    seed = np.array([0.55, 0.0, 0.0, 0.0, 0.42, 0.02], dtype=float)

    shooter = HierarchicalShooter(p)
    result = shooter.solve(seed)

    s0 = result["state0"]
    period = float(result["period"])
    stage_data = result["stages"]

    t_eval = np.linspace(0.0, period * 6.0, 4000)
    sol_cr3bp = propagate(cr3bp_rhs, s0, (0.0, t_eval[-1]), p, t_eval=t_eval)
    sol_rf3bp = propagate(rf3bp_pulsating_rhs, s0, (0.0, t_eval[-1]), p, t_eval=t_eval)

    print("Hierarchical shooting summary")
    print(f"period={period:.6f}")
    print(f"residual_norm={result['residual_norm']:.3e}")
    print(f"stage1_cost={result['stage1_cost']:.3e}, final_cost={result['stage2_cost']:.3e}")
    print(f"stages={len(stage_data)}")

    out_dir = os.environ.get("OUTPUT_DIR", ".")
    os.makedirs(out_dir, exist_ok=True)

    plot_trajectory(sol_cr3bp.y[:3], "CR3BP Reference Orbit")
    plt.savefig(os.path.join(out_dir, "trajectory_cr3bp.png"), dpi=160)
    plt.close()

    plot_trajectory(sol_rf3bp.y[:3], "RF3BP Orbit with Pulsation + Perturbations")
    plt.savefig(os.path.join(out_dir, "trajectory_rf3bp.png"), dpi=160)
    plt.close()

    perturbation = _component_norm_history(sol_rf3bp.t, sol_rf3bp.y, p)
    plot_perturbation_norms(sol_rf3bp.t, perturbation)
    plt.savefig(os.path.join(out_dir, "perturbation_norms.png"), dpi=160)
    plt.close()

    stage_labels = [s.label for s in stage_data]
    stage_residuals = np.array([s.residual_norm for s in stage_data])
    plot_stage_convergence(stage_labels, stage_residuals)
    plt.savefig(os.path.join(out_dir, "continuation_convergence.png"), dpi=160)
    plt.close()

    if os.environ.get("SHOW_PLOTS", "0") == "1":
        plt.show()

    print(f"saved plots in {out_dir}")


if __name__ == "__main__":
    main()