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
