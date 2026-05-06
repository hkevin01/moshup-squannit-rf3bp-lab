from __future__ import annotations

import json
import os

import matplotlib.pyplot as plt
import numpy as np

from rf3bp_lab.dynamics.models import FidelityWeights, compare_cr3bp_rf3bp, cr3bp_rhs, propagate, rf3bp_breakdown, rf3bp_pulsating_rhs
from rf3bp_lab.dynamics.params import SystemParams
from rf3bp_lab.shooting.hierarchical import HierarchicalShooter
from rf3bp_lab.utils.plotting import plot_model_gap, plot_perturbation_norms, plot_result_dashboard, plot_stage_convergence, plot_trajectory


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


def _model_gap_history(t: np.ndarray, states: np.ndarray, p: SystemParams) -> tuple[np.ndarray, np.ndarray]:
    absolute_gap = np.zeros_like(t)
    relative_gap = np.zeros_like(t)
    for i in range(t.size):
        gap = compare_cr3bp_rf3bp(float(t[i]), states[:, i], p)
        absolute_gap[i] = gap.delta_norm
        relative_gap[i] = gap.relative_gap
    return absolute_gap, relative_gap


def _write_json(path: str, payload: dict) -> None:
    with open(path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2, sort_keys=True)


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

    gap_abs, gap_rel = _model_gap_history(sol_rf3bp.t, sol_rf3bp.y, p)
    plot_model_gap(sol_rf3bp.t, gap_abs, gap_rel)
    plt.savefig(os.path.join(out_dir, "model_gap_cr3bp_vs_rf3bp.png"), dpi=160)
    plt.close()

    perturbation_peaks = {k: float(np.max(v)) for k, v in perturbation.items()}
    summary_metrics = {
        "period": period,
        "residual_norm": float(result["residual_norm"]),
        "stage1_cost": float(result["stage1_cost"]),
        "final_cost": float(result["stage2_cost"]),
        "max_abs_gap": float(np.max(gap_abs)),
        "max_rel_gap": float(np.max(gap_rel)),
        "mean_abs_gap": float(np.mean(gap_abs)),
        "mean_rel_gap": float(np.mean(gap_rel)),
        "n_steps": int(sol_rf3bp.t.size),
    }
    plot_result_dashboard(
        summary_metrics,
        perturbation_peaks,
        stage_labels=stage_labels,
        stage_residuals=stage_residuals,
    )
    plt.savefig(os.path.join(out_dir, "result_snapshot_dashboard.png"), dpi=160)
    plt.close()

    results_dir = os.path.join(os.path.dirname(out_dir), "results") if os.path.basename(out_dir) == "figures" else os.path.join(out_dir, "results")
    os.makedirs(results_dir, exist_ok=True)
    _write_json(
        os.path.join(results_dir, "latest_demo_metrics.json"),
        {
            "summary_metrics": summary_metrics,
            "perturbation_peaks": perturbation_peaks,
            "stage_residuals": [float(s.residual_norm) for s in stage_data],
            "stage_labels": [s.label for s in stage_data],
        },
    )

    if os.environ.get("SHOW_PLOTS", "0") == "1":
        plt.show()

    print(f"saved plots in {out_dir}")


if __name__ == "__main__":
    main()