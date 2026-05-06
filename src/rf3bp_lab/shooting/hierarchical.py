from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import least_squares

from rf3bp_lab.dynamics.models import cr3bp_rhs, propagate, rf3bp_pulsating_rhs
from rf3bp_lab.dynamics.params import SystemParams


@dataclass(frozen=True)
class ShootingConfig:
    period_guess: float = 6.0
    max_nfev: int = 120


class HierarchicalShooter:
    def __init__(self, params: SystemParams, config: ShootingConfig | None = None):
        self.params = params
        self.config = config or ShootingConfig()

    def _periodicity_residual(self, x: np.ndarray, rhs) -> np.ndarray:
        state0 = x[:6]
        tf = float(x[6])
        sol = propagate(rhs, state0, (0.0, tf), self.params)
        sf = sol.y[:, -1]
        # Simple planar periodic constraints around a Poincare-like return
        return np.array([
            sf[1] - state0[1],
            sf[2] - state0[2],
            sf[4] - state0[4],
            sf[5] - state0[5],
            sf[0] - state0[0],
            sf[3] - state0[3],
        ])

    def solve(self, seed_state: np.ndarray) -> dict[str, np.ndarray | float]:
        x0 = np.hstack([seed_state, self.config.period_guess])

        # Stage 1: CR3BP seed correction
        stage1 = least_squares(
            lambda x: self._periodicity_residual(x, cr3bp_rhs),
            x0,
            max_nfev=self.config.max_nfev,
        )

        # Stage 2: continuation into pulsating RF3BP + perturbations
        stage2 = least_squares(
            lambda x: self._periodicity_residual(x, rf3bp_pulsating_rhs),
            stage1.x,
            max_nfev=self.config.max_nfev,
        )

        return {
            "state0": stage2.x[:6],
            "period": float(stage2.x[6]),
            "residual_norm": float(np.linalg.norm(stage2.fun)),
            "stage1_cost": float(stage1.cost),
            "stage2_cost": float(stage2.cost),
        }
