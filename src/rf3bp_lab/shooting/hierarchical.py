from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import least_squares

from rf3bp_lab.dynamics.models import (
    FidelityWeights,
    cr3bp_rhs,
    propagate,
    rf3bp_pulsating_rhs_weighted,
)
from rf3bp_lab.dynamics.params import SystemParams


@dataclass(frozen=True)
class ShootingConfig:
    period_guess: float = 6.0
    max_nfev: int = 140
    continuation_steps: int = 5
    loss: str = "soft_l1"


@dataclass(frozen=True)
class StageResult:
    label: str
    fidelity: FidelityWeights
    cost: float
    residual_norm: float
    period: float


def _continuation_fidelities(steps: int) -> list[FidelityWeights]:
    values = np.linspace(0.0, 1.0, max(steps, 2))
    out: list[FidelityWeights] = [FidelityWeights(0.0, 0.0, 0.0, 0.0)]

    for v in values[1:]:
        out.append(FidelityWeights(pulsation=float(v), nonspherical=0.0, solar_gravity=0.0, srp=0.0))
    for v in values[1:]:
        out.append(FidelityWeights(pulsation=1.0, nonspherical=float(v), solar_gravity=0.0, srp=0.0))
    for v in values[1:]:
        out.append(FidelityWeights(pulsation=1.0, nonspherical=1.0, solar_gravity=float(v), srp=0.0))
    for v in values[1:]:
        out.append(FidelityWeights(pulsation=1.0, nonspherical=1.0, solar_gravity=1.0, srp=float(v)))

    return out


class HierarchicalShooter:
    def __init__(self, params: SystemParams, config: ShootingConfig | None = None):
        self.params = params
        self.config = config or ShootingConfig()

    def _periodicity_residual(self, x: np.ndarray, rhs) -> np.ndarray:
        state0 = x[:6]
        tf = float(x[6])
        sol = propagate(rhs, state0, (0.0, tf), self.params)
        sf = sol.y[:, -1]
        return np.array([
            sf[1] - state0[1],
            sf[2] - state0[2],
            sf[4] - state0[4],
            sf[5] - state0[5],
            sf[0] - state0[0],
            sf[3] - state0[3],
        ])

    def _solve_stage(self, x0: np.ndarray, rhs, label: str, fidelity: FidelityWeights) -> tuple[np.ndarray, StageResult]:
        fit = least_squares(
            lambda x: self._periodicity_residual(x, rhs),
            x0,
            max_nfev=self.config.max_nfev,
            loss=self.config.loss,
        )
        stage = StageResult(
            label=label,
            fidelity=fidelity,
            cost=float(fit.cost),
            residual_norm=float(np.linalg.norm(fit.fun)),
            period=float(fit.x[6]),
        )
        return fit.x, stage

    def solve(self, seed_state: np.ndarray) -> dict[str, np.ndarray | float | list[StageResult]]:
        x = np.hstack([seed_state, self.config.period_guess])

        x, stage_cr3bp = self._solve_stage(x, cr3bp_rhs, "cr3bp", FidelityWeights(0.0, 0.0, 0.0, 0.0))
        stages: list[StageResult] = [stage_cr3bp]

        for i, fidelity in enumerate(_continuation_fidelities(self.config.continuation_steps)):
            if i == 0:
                continue
            rhs = lambda t, y, p, fw=fidelity: rf3bp_pulsating_rhs_weighted(t, y, p, fw)
            x, stage = self._solve_stage(x, rhs, f"continuation_{i:02d}", fidelity)
            stages.append(stage)

        return {
            "state0": x[:6],
            "period": float(x[6]),
            "residual_norm": stages[-1].residual_norm,
            "stage1_cost": stages[0].cost,
            "stage2_cost": stages[-1].cost,
            "stages": stages,
        }