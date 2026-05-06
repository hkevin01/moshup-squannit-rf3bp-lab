from __future__ import annotations

from dataclasses import dataclass

import numpy as np

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
    max_nfev: int = 4
    continuation_steps: int = 2
    period_bounds: tuple[float, float] = (1.0, 10.0)
    state_bound: float = 2.0
    solve_tol: float = 5e-3
    jac_step: float = 1e-4
    damping: float = 0.7
    integrator_method: str = "DOP853"
    integrator_rtol: float = 1e-7
    integrator_atol: float = 1e-9
    integrator_max_step: float = 0.2


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
        lower, upper = self.config.period_bounds
        if (not np.isfinite(tf)) or tf < lower or tf > upper or np.any(~np.isfinite(state0)):
            return np.full(6, 1.0e3, dtype=float)
        sol = propagate(
            rhs,
            state0,
            (0.0, tf),
            self.params,
            method=self.config.integrator_method,
            rtol=self.config.integrator_rtol,
            atol=self.config.integrator_atol,
            max_step=self.config.integrator_max_step,
        )
        if not sol.success or sol.y.shape[1] == 0:
            return np.full(6, 1.0e3, dtype=float)
        sf = sol.y[:, -1]
        return np.array([
            sf[1] - state0[1],
            sf[2] - state0[2],
            sf[4] - state0[4],
            sf[5] - state0[5],
            sf[0] - state0[0],
            sf[3] - state0[3],
        ])

    def _clip_decision(self, x: np.ndarray) -> np.ndarray:
        state_lower = np.full(6, -self.config.state_bound, dtype=float)
        state_upper = np.full(6, self.config.state_bound, dtype=float)
        period_lower, period_upper = self.config.period_bounds
        lower = np.hstack([state_lower, period_lower])
        upper = np.hstack([state_upper, period_upper])
        return np.clip(x, lower, upper)

    def _finite_difference_jacobian(self, x: np.ndarray, rhs, residual: np.ndarray) -> np.ndarray:
        jac = np.zeros((residual.size, x.size), dtype=float)
        for j in range(x.size):
            step = self.config.jac_step * max(1.0, abs(float(x[j])))
            xp = x.copy()
            xp[j] += step
            xp = self._clip_decision(xp)
            rp = self._periodicity_residual(xp, rhs)
            jac[:, j] = (rp - residual) / step
        return jac

    def _solve_stage(self, x0: np.ndarray, rhs, label: str, fidelity: FidelityWeights) -> tuple[np.ndarray, StageResult]:
        x = self._clip_decision(x0)
        residual = self._periodicity_residual(x, rhs)
        best_x = x.copy()
        best_residual = residual.copy()

        for _ in range(self.config.max_nfev):
            residual_norm = float(np.linalg.norm(residual))
            if residual_norm <= self.config.solve_tol:
                break

            jac = self._finite_difference_jacobian(x, rhs, residual)
            step, *_ = np.linalg.lstsq(jac, -residual, rcond=None)

            candidate = self._clip_decision(x + self.config.damping * step)
            candidate_residual = self._periodicity_residual(candidate, rhs)

            if np.linalg.norm(candidate_residual) < np.linalg.norm(residual):
                x = candidate
                residual = candidate_residual
            else:
                reduced = self._clip_decision(x + 0.5 * self.config.damping * step)
                reduced_residual = self._periodicity_residual(reduced, rhs)
                if np.linalg.norm(reduced_residual) < np.linalg.norm(residual):
                    x = reduced
                    residual = reduced_residual
                else:
                    break

            if np.linalg.norm(residual) < np.linalg.norm(best_residual):
                best_x = x.copy()
                best_residual = residual.copy()

        x = best_x
        residual = best_residual
        stage = StageResult(
            label=label,
            fidelity=fidelity,
            cost=0.5 * float(best_residual @ best_residual),
            residual_norm=float(np.linalg.norm(best_residual)),
            period=float(best_x[6]),
        )
        return best_x, stage

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