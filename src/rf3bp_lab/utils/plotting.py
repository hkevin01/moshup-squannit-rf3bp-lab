from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np


def plot_trajectory(states: np.ndarray, title: str = "Trajectory") -> None:
    fig = plt.figure(figsize=(7, 6))
    ax = fig.add_subplot(111, projection="3d")
    ax.plot(states[0], states[1], states[2], lw=1.5)
    ax.scatter(states[0, 0], states[1, 0], states[2, 0], s=40, marker="o", label="start")
    ax.scatter(states[0, -1], states[1, -1], states[2, -1], s=40, marker="x", label="end")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    ax.set_title(title)
    ax.legend()
    plt.tight_layout()


def plot_trajectory_planes(
    states_cr3bp: np.ndarray,
    states_rf3bp: np.ndarray,
    title: str = "Trajectory Projections: CR3BP vs RF3BP",
) -> None:
    """Four-panel view: 3-D plus three 2-D plane projections, both models overlaid.

    Parameters
    ----------
    states_cr3bp, states_rf3bp : ndarray, shape (3+, N)
        Position rows are [x, y, z, ...].  Only the first three rows are used.
    """
    C1, C2 = "#1d3557", "#e76f51"   # CR3BP blue, RF3BP orange
    fig = plt.figure(figsize=(12, 9))
    fig.suptitle(title, fontsize=13, fontweight="bold")

    # --- 3-D panel (top-left, larger) -----------------------------------------
    ax3d = fig.add_subplot(2, 2, 1, projection="3d")
    ax3d.plot(states_cr3bp[0], states_cr3bp[1], states_cr3bp[2],
              lw=1.4, color=C1, label="CR3BP", alpha=0.85)
    ax3d.plot(states_rf3bp[0], states_rf3bp[1], states_rf3bp[2],
              lw=1.4, color=C2, label="RF3BP", alpha=0.85, linestyle="--")
    ax3d.scatter(*states_cr3bp[:3, 0],  s=30, color=C1, marker="o", zorder=5)
    ax3d.scatter(*states_rf3bp[:3, 0],  s=30, color=C2, marker="o", zorder=5)
    ax3d.set_xlabel("x"); ax3d.set_ylabel("y"); ax3d.set_zlabel("z")
    ax3d.set_title("3-D view")
    ax3d.legend(fontsize=8)

    # --- XY projection (top-right) -------------------------------------------
    ax_xy = fig.add_subplot(2, 2, 2)
    ax_xy.plot(states_cr3bp[0], states_cr3bp[1], lw=1.3, color=C1, label="CR3BP")
    ax_xy.plot(states_rf3bp[0], states_rf3bp[1], lw=1.3, color=C2,
               label="RF3BP", linestyle="--", alpha=0.9)
    ax_xy.scatter(0, 0, s=60, color="#888", marker="+", zorder=5, label="CoM")
    ax_xy.set_xlabel("x"); ax_xy.set_ylabel("y")
    ax_xy.set_title("XY plane (equatorial)")
    ax_xy.set_aspect("equal", adjustable="datalim")
    ax_xy.grid(alpha=0.25); ax_xy.legend(fontsize=8)

    # --- XZ projection (bottom-left) -----------------------------------------
    ax_xz = fig.add_subplot(2, 2, 3)
    ax_xz.plot(states_cr3bp[0], states_cr3bp[2], lw=1.3, color=C1, label="CR3BP")
    ax_xz.plot(states_rf3bp[0], states_rf3bp[2], lw=1.3, color=C2,
               label="RF3BP", linestyle="--", alpha=0.9)
    ax_xz.set_xlabel("x"); ax_xz.set_ylabel("z")
    ax_xz.set_title("XZ plane (out-of-plane)")
    ax_xz.set_aspect("equal", adjustable="datalim")
    ax_xz.grid(alpha=0.25); ax_xz.legend(fontsize=8)

    # --- YZ projection (bottom-right) ----------------------------------------
    ax_yz = fig.add_subplot(2, 2, 4)
    ax_yz.plot(states_cr3bp[1], states_cr3bp[2], lw=1.3, color=C1, label="CR3BP")
    ax_yz.plot(states_rf3bp[1], states_rf3bp[2], lw=1.3, color=C2,
               label="RF3BP", linestyle="--", alpha=0.9)
    ax_yz.set_xlabel("y"); ax_yz.set_ylabel("z")
    ax_yz.set_title("YZ plane (meridional)")
    ax_yz.set_aspect("equal", adjustable="datalim")
    ax_yz.grid(alpha=0.25); ax_yz.legend(fontsize=8)

    plt.tight_layout(rect=[0, 0, 1, 0.96])


def plot_trajectory_velocity(
    time: np.ndarray,
    states_cr3bp: np.ndarray,
    states_rf3bp: np.ndarray,
    title: str = "Velocity Components Over Time",
) -> None:
    """Three-panel velocity time-history: vx, vy, vz for both models."""
    C1, C2 = "#1d3557", "#e76f51"
    labels = ["vx", "vy", "vz"]
    fig, axes = plt.subplots(3, 1, figsize=(10, 7), sharex=True)
    fig.suptitle(title, fontsize=12, fontweight="bold")
    for i, (ax, lbl) in enumerate(zip(axes, labels)):
        ax.plot(time, states_cr3bp[3 + i], lw=1.3, color=C1, label="CR3BP")
        ax.plot(time, states_rf3bp[3 + i], lw=1.3, color=C2,
                label="RF3BP", linestyle="--", alpha=0.9)
        ax.set_ylabel(lbl)
        ax.grid(alpha=0.25)
        ax.legend(fontsize=8, loc="upper right")
    axes[-1].set_xlabel("t (normalised)")
    plt.tight_layout(rect=[0, 0, 1, 0.96])


def plot_trajectory_deviation(
    time: np.ndarray,
    states_cr3bp: np.ndarray,
    states_rf3bp: np.ndarray,
    title: str = "Position Deviation: RF3BP - CR3BP",
) -> None:
    """Show the growing positional difference between the two models over time."""
    delta = states_rf3bp[:3] - states_cr3bp[:3]
    dist  = np.linalg.norm(delta, axis=0)

    C_X, C_Y, C_Z, C_D = "#e63946", "#2a9d8f", "#f4a261", "#1d3557"
    fig, (ax_comp, ax_norm) = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
    fig.suptitle(title, fontsize=12, fontweight="bold")

    ax_comp.plot(time, delta[0], lw=1.2, color=C_X, label="dx")
    ax_comp.plot(time, delta[1], lw=1.2, color=C_Y, label="dy")
    ax_comp.plot(time, delta[2], lw=1.2, color=C_Z, label="dz")
    ax_comp.set_ylabel("component deviation")
    ax_comp.grid(alpha=0.25)
    ax_comp.legend(ncol=3, fontsize=9)
    ax_comp.axhline(0, color="#aaa", lw=0.8, linestyle=":")

    ax_norm.plot(time, dist, lw=1.5, color=C_D)
    ax_norm.fill_between(time, 0, dist, alpha=0.15, color=C_D)
    ax_norm.set_ylabel(r"$\|\mathbf{r}_{RF3BP} - \mathbf{r}_{CR3BP}\|$")
    ax_norm.set_xlabel("t (normalised)")
    ax_norm.grid(alpha=0.25)

    plt.tight_layout(rect=[0, 0, 1, 0.96])


def plot_perturbation_norms(time: np.ndarray, component_norms: dict[str, np.ndarray], title: str = "Perturbation Norms") -> None:
    fig, ax = plt.subplots(figsize=(8.5, 4.8))
    for label, values in component_norms.items():
        ax.plot(time, values, lw=1.4, label=label)
    ax.set_xlabel("t")
    ax.set_ylabel("acceleration norm")
    ax.set_yscale("log")
    ax.set_title(title)
    ax.grid(alpha=0.25)
    ax.legend(ncol=2, fontsize=9)
    plt.tight_layout()


def plot_stage_convergence(stage_labels: list[str], residual_norms: np.ndarray, title: str = "Hierarchical Continuation Convergence") -> None:
    fig, ax = plt.subplots(figsize=(9.0, 4.6))
    x = np.arange(len(stage_labels))
    ax.plot(x, residual_norms, marker="o", lw=1.5)
    ax.set_xticks(x)
    ax.set_xticklabels(stage_labels, rotation=60, ha="right", fontsize=8)
    ax.set_yscale("log")
    ax.set_ylabel("residual norm")
    ax.set_title(title)
    ax.grid(alpha=0.3)
    plt.tight_layout()


def plot_model_gap(time: np.ndarray, absolute_gap: np.ndarray, relative_gap: np.ndarray, title: str = "RF3BP vs CR3BP Acceleration Gap") -> None:
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8.5, 6.4), sharex=True)

    ax1.plot(time, absolute_gap, color="#ae2012", lw=1.5)
    ax1.set_yscale("log")
    ax1.set_ylabel(r"$\|a_{RF3BP} - a_{CR3BP}\|$")
    ax1.grid(alpha=0.25)
    ax1.set_title(title)

    ax2.plot(time, relative_gap, color="#005f73", lw=1.5)
    ax2.set_yscale("log")
    ax2.set_xlabel("t")
    ax2.set_ylabel(r"$\|\Delta a\| / \|a_{CR3BP}\|$")
    ax2.grid(alpha=0.25)

    plt.tight_layout()


def plot_result_dashboard(
    summary_metrics: dict[str, float],
    perturbation_peaks: dict[str, float],
    stage_labels: list[str] | None = None,
    stage_residuals: np.ndarray | None = None,
    title: str = "RF3BP Demo Result Snapshot",
) -> None:
    fig = plt.figure(figsize=(13.0, 7.6), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, wspace=0.28, hspace=0.34)

    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, 0])
    ax4 = fig.add_subplot(gs[1, 1])

    labels = list(perturbation_peaks.keys())
    values = np.array([perturbation_peaks[k] for k in labels], dtype=float)
    colors = ["#b00020", "#e76f51", "#1d3557", "#2a9d8f"]
    bars = ax1.bar(labels, values, color=colors, alpha=0.92)
    ax1.set_yscale("log")
    ax1.set_ylabel("Peak Acceleration Norm")
    ax1.set_title("Perturbation Dominance", fontweight="bold")
    ax1.grid(alpha=0.2, axis="y")
    ax1.set_axisbelow(True)
    for bar, val in zip(bars, values):
        y = max(val * 1.2, 1e-20)
        ax1.text(bar.get_x() + bar.get_width() / 2.0, y, f"{val:.2e}", ha="center", va="bottom", fontsize=8)

    if stage_labels is not None and stage_residuals is not None and len(stage_labels) == len(stage_residuals):
        x = np.arange(len(stage_labels))
        ax2.plot(x, stage_residuals, marker="o", lw=2.0, color="#264653")
        ax2.fill_between(x, stage_residuals, np.max(stage_residuals), alpha=0.12, color="#264653")
        ax2.set_xticks(x)
        ax2.set_xticklabels(stage_labels, rotation=40, ha="right", fontsize=8)
        ax2.set_yscale("log")
        ax2.set_ylabel("Residual Norm")
        ax2.set_title("Continuation Progress", fontweight="bold")
        ax2.grid(alpha=0.2)
    else:
        ax2.axis("off")
        ax2.text(0.5, 0.5, "Continuation data unavailable", ha="center", va="center", fontsize=11)

    abs_gap = summary_metrics["max_abs_gap"]
    mean_gap = summary_metrics["mean_abs_gap"]
    rel_max = summary_metrics["max_rel_gap"]
    rel_mean = summary_metrics["mean_rel_gap"]
    gap_labels = ["Max |delta a|", "Mean |delta a|", "Max rho", "Mean rho"]
    gap_vals = [abs_gap, mean_gap, rel_max, rel_mean]
    gap_colors = ["#c1121f", "#f3722c", "#277da1", "#43aa8b"]
    y_idx = np.arange(len(gap_labels))
    ax3.barh(y_idx, gap_vals, color=gap_colors, alpha=0.9)
    ax3.set_xscale("log")
    ax3.set_yticks(y_idx)
    ax3.set_yticklabels(gap_labels)
    ax3.set_xlabel("Magnitude (log scale)")
    ax3.set_title("RF3BP vs CR3BP Gap Indicators", fontweight="bold")
    ax3.grid(alpha=0.2, axis="x")
    for idx, val in enumerate(gap_vals):
        ax3.text(val * 1.08, idx, f"{val:.2e}", va="center", fontsize=8)

    ax4.axis("off")
    ax4.set_title("Mission Snapshot", fontweight="bold")
    cards = [
        ("Period", f"{summary_metrics['period']:.6f}", "#d9ed92"),
        ("Residual", f"{summary_metrics['residual_norm']:.3e}", "#ffddd2"),
        ("Stage-1 Cost", f"{summary_metrics['stage1_cost']:.3e}", "#e9edc9"),
        ("Final Cost", f"{summary_metrics['final_cost']:.3e}", "#fde2e4"),
    ]
    x0, y0 = 0.03, 0.86
    dy = 0.2
    for i, (label, value, face) in enumerate(cards):
        y = y0 - i * dy
        ax4.text(
            x0,
            y,
            f"{label}\n{value}",
            transform=ax4.transAxes,
            fontsize=10,
            va="top",
            bbox={"boxstyle": "round,pad=0.35", "facecolor": face, "edgecolor": "#7a7a7a"},
        )

    dominant_name = labels[int(np.argmax(values))] if len(values) else "n/a"
    dominant_val = float(np.max(values)) if len(values) else 0.0
    ax4.text(
        0.52,
        0.36,
        "Visual Aids",
        transform=ax4.transAxes,
        fontsize=11,
        fontweight="bold",
    )
    ax4.text(
        0.52,
        0.25,
        f"Dominant perturbation: {dominant_name}",
        transform=ax4.transAxes,
        fontsize=9,
    )
    ax4.text(
        0.52,
        0.16,
        f"Peak value: {dominant_val:.2e}",
        transform=ax4.transAxes,
        fontsize=9,
    )
    ax4.text(
        0.52,
        0.07,
        "Lower residual trend indicates better stage-to-stage correction.",
        transform=ax4.transAxes,
        fontsize=8,
        color="#404040",
    )

    fig.suptitle(title, fontsize=15, fontweight="bold")