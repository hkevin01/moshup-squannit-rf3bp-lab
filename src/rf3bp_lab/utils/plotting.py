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