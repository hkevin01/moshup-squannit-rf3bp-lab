from __future__ import annotations

"""ID: RF3BP-LAB-SHOOT-MULTI
Requirement: Provide multi-arc shooting continuation for complex trajectories.
Purpose: Solve periodic or constrained trajectories divided into N arcs.
Rationale: Improves convergence for long-duration orbits and enables midcourse constraints.
Inputs: Seed states, arc configuration, fidelity weights, and solver parameters.
Outputs: Corrected states, arc residuals, and convergence metrics.
Preconditions: Number of arcs >= 2, arc times are finite and ordered.
Postconditions: Continuity constraints are satisfied to solver tolerance.
Assumptions: Smooth dynamics, well-posed boundary conditions.
Side Effects: None.
Failure Modes: Poor arc partitioning can cause stiff residual mapping.
Error Handling: Damped Newton, step bounds, residual penalties.
Constraints: Designed for O(N) arcs; N << 100 recommended.
Verification: Unit tests in test_shooting.py.
References: Gonzalez-Camachoa et al., J. Guidance Control Dyn (2015).
"""

from dataclasses import dataclass
import numpy as np

from rf3bp_lab.dynamics.models import (
    FidelityWeights,
    rf3bp_pulsating_rhs_weighted,
    propagate,
)
from rf3bp_lab.dynamics.params import SystemParams


@dataclass(frozen=True)
class ArcDefinition:
    """ID: RF3BP-LAB-SHOOT-MULTI-ARCDEF
    Requirement: Define one shooting arc with time span and boundary treatment.
    Purpose: Specify arc geometry for multi-arc solver.
    Inputs: Start time, end time, arc index, is_periodic flag.
    Outputs: Immutable arc descriptor.
    """
    t_start: float
    t_end: float
    arc_index: int
    is_periodic: bool = False


@dataclass(frozen=True)
class MultiShootingConfig:
    """ID: RF3BP-LAB-SHOOT-MULTI-CONFIG
    Requirement: Hold configuration for multi-arc shooting solver.
    Purpose: Centralize tunables for reproducible shooting behavior.
    """
    num_arcs: int = 3
    max_nfev: int = 10
    solve_tol: float = 1e-3
    damping: float = 0.7
    arc_period_guess: float = 2.0
    jac_step: float = 1e-4
    state_bound: float = 2.0
    integrator_method: str = "DOP853"
    integrator_rtol: float = 1e-7
    integrator_atol: float = 1e-9


@dataclass(frozen=True)
class ArcResult:
    """ID: RF3BP-LAB-SHOOT-MULTI-ARCRESULT
    Requirement: Store outputs for one arc in multi-shooting solution.
    Purpose: Preserve arc-level diagnostics.
    """
    arc_index: int
    initial_state: np.ndarray
    final_state: np.ndarray
    residual_norm: float
    continuity_residual: np.ndarray | None = None


class MultiShootingBVP:
    """ID: RF3BP-LAB-SHOOT-MULTI-BVP
    Requirement: Solve multi-arc boundary-value problem via simultaneous shooting.
    Purpose: Find periodic or constrained orbits with arc continuity.
    Rationale: Better numerical stability than single-arc for long orbits.
    """

    def __init__(
        self,
        params: SystemParams,
        config: MultiShootingConfig | None = None,
    ):
        """ID: RF3BP-LAB-SHOOT-MULTI-BVP-INIT"""
        self.params = params
        self.config = config or MultiShootingConfig()
        self.arc_results: list[ArcResult] = []

    def solve(
        self,
        initial_states: np.ndarray,
        arc_defs: list[ArcDefinition],
        fidelity: FidelityWeights,
    ) -> dict:
        """ID: RF3BP-LAB-SHOOT-MULTI-BVP-SOLVE
        Requirement: Solve multi-arc BVP by simultaneous Newton correction.
        Inputs: Array of initial states (one per arc), arc definitions, fidelity weights.
        Outputs: Dictionary with corrected states, residuals, arc results.
        Preconditions: initial_states shape is (N_arcs, 6).
        Postconditions: Returns converged solution or partial results if iteration limit hit.
        Rationale: Each Newton iteration corrects all N_arcs states simultaneously.
        """

        if len(arc_defs) != len(initial_states):
            raise ValueError("Number of arcs must match number of initial states")

        n_arcs = len(arc_defs)
        state_dim = 6
        total_vars = n_arcs * state_dim

        # Copy initial states for iteration
        states = initial_states.copy()
        
        residual_history = []
        self.arc_results = []

        for nfev in range(self.config.max_nfev):
            # Propagate each arc and collect residuals
            arc_results = []
            all_residuals = np.zeros(total_vars)
            
            for i, arc_def in enumerate(arc_defs):
                t_span = (arc_def.t_start, arc_def.t_end)
                state0 = states[i].copy()
                
                # Propagate arc
                result = propagate(
                    lambda t, y, p: rf3bp_pulsating_rhs_weighted(t, y, p, fidelity),
                    state0,
                    t_span,
                    self.params,
                    method=self.config.integrator_method,
                    rtol=self.config.integrator_rtol,
                    atol=self.config.integrator_atol,
                )
                
                final_state = result.y[:, -1] if result.success else state0
                residual = final_state - states[(i + 1) % n_arcs]  # Continuity residual
                
                if arc_def.is_periodic and i == n_arcs - 1:
                    # Last arc closes the loop
                    residual = final_state - states[0]
                
                all_residuals[i * state_dim:(i + 1) * state_dim] = residual
                
                arc_results.append(ArcResult(
                    arc_index=i,
                    initial_state=state0,
                    final_state=final_state,
                    residual_norm=float(np.linalg.norm(residual)),
                    continuity_residual=residual,
                ))
            
            self.arc_results = arc_results
            residual_norm = float(np.linalg.norm(all_residuals))
            residual_history.append(residual_norm)
            
            if residual_norm < self.config.solve_tol:
                return {
                    'success': True,
                    'states': states,
                    'residual_norm': residual_norm,
                    'arc_results': arc_results,
                    'residual_history': residual_history,
                    'iterations': nfev + 1,
                }
            
            # Compute finite-difference Jacobian
            jacobian = np.zeros((total_vars, total_vars))
            for j in range(total_vars):
                states_pert = states.copy()
                states_pert.flat[j] += self.config.jac_step
                
                residuals_pert = np.zeros(total_vars)
                for i, arc_def in enumerate(arc_defs):
                    t_span = (arc_def.t_start, arc_def.t_end)
                    state0 = states_pert[i].copy()
                    
                    result = propagate(
                        lambda t, y, p: rf3bp_pulsating_rhs_weighted(t, y, p, fidelity),
                        state0,
                        t_span,
                        self.params,
                        method=self.config.integrator_method,
                        rtol=self.config.integrator_rtol,
                        atol=self.config.integrator_atol,
                    )
                    
                    final_state_pert = result.y[:, -1] if result.success else state0
                    residual_pert = final_state_pert - states_pert[(i + 1) % n_arcs]
                    
                    if arc_def.is_periodic and i == n_arcs - 1:
                        residual_pert = final_state_pert - states_pert[0]
                    
                    residuals_pert[i * state_dim:(i + 1) * state_dim] = residual_pert
                
                jacobian[:, j] = (residuals_pert - all_residuals) / self.config.jac_step
            
            # Solve Newton step with damping
            try:
                delta = np.linalg.solve(jacobian, -all_residuals)
            except np.linalg.LinAlgError:
                # Singular Jacobian - use pseudoinverse
                delta = -np.linalg.pinv(jacobian) @ all_residuals
            
            # Damping and bounds
            delta_clipped = np.clip(delta, -self.config.state_bound, self.config.state_bound)
            states += self.config.damping * delta_clipped.reshape(states.shape)
            
            # Clip state components
            states = np.clip(states, -self.config.state_bound, self.config.state_bound)
        
        return {
            'success': False,
            'states': states,
            'residual_norm': residual_history[-1] if residual_history else np.inf,
            'arc_results': self.arc_results,
            'residual_history': residual_history,
            'iterations': self.config.max_nfev,
        }
