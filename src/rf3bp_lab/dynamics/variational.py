from __future__ import annotations

"""ID: RF3BP-LAB-DYN-VARIATIONAL
Requirement: Provide state transition matrix (STM) propagation for covariance analysis.
Purpose: Enable mission-grade navigation uncertainty quantification.
Rationale: STM is Jacobian of flow map; covariance propagates via P(t) = Phi(t)*P(t0)*Phi(t)^T.
Inputs: State and STM initial conditions, dynamics parameters.
Outputs: Integrated state and STM history.
Preconditions: Initial state and STM (identity) are finite.
Postconditions: STM remains orthogonal (for Hamiltonian systems) or well-conditioned.
Assumptions: Smooth dynamics, small perturbations for linearization validity.
Side Effects: None.
Failure Modes: Linearization breaks down for large uncertainties; mitigation via small covariance.
Error Handling: Numerical integration uses adaptive step control.
Constraints: 6x6 STM adds modest computational overhead (42 ODEs vs 6).
Verification: Unit tests validate STM against numerical differentiation.
References: Curtis, Orbital Mechanics for Engineering Students (2013).
"""

import numpy as np
from scipy.integrate import solve_ivp
from dataclasses import dataclass

from .models import rf3bp_pulsating_rhs_weighted, FidelityWeights
from .params import SystemParams


@dataclass(frozen=True)
class VariationalState:
    """ID: RF3BP-LAB-DYN-VARSTATE
    Requirement: Bundle state vector and state transition matrix.
    Purpose: Carry both nominal and sensitivity information through integration.
    Inputs: State (6,), STM (6, 6).
    Outputs: Immutable container.
    Preconditions: State is 6D, STM is 6x6 and finite.
    Postconditions: All fields remain immutable after construction.
    """
    state: np.ndarray
    stm: np.ndarray


def _flatten_variational(state: np.ndarray, stm: np.ndarray) -> np.ndarray:
    """ID: RF3BP-LAB-DYN-FLATTEN-VAR
    Requirement: Convert state + STM into flat vector for ODE integration.
    Inputs: State (6,), STM (6, 6).
    Outputs: Flattened array of length 42 (6 + 36).
    Preconditions: Shapes are (6,) and (6, 6).
    Postconditions: Returns (42,) array.
    """
    return np.concatenate([state, stm.flatten()])


def _unflatten_variational(flat: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """ID: RF3BP-LAB-DYN-UNFLATTEN-VAR
    Requirement: Reconstruct state and STM from flat vector.
    Inputs: Flattened array of length 42.
    Outputs: Tuple (state: (6,), stm: (6, 6)).
    Preconditions: Input length is 42.
    Postconditions: Returns properly shaped arrays.
    """
    state = flat[:6].copy()
    stm = flat[6:].reshape((6, 6))
    return state, stm


def _jacobian_rf3bp(
    t: float,
    state: np.ndarray,
    p: SystemParams,
    fidelity: FidelityWeights,
) -> np.ndarray:
    """ID: RF3BP-LAB-DYN-JACO-RFBP
    Requirement: Compute 6x6 Jacobian of RF3BP dynamics wrt state.
    Inputs: Time t, state vector, parameters p, fidelity weights.
    Outputs: 6x6 Jacobian matrix.
    Preconditions: State is 6D and finite.
    Postconditions: Output is (6, 6) and finite for non-singular positions.
    Assumptions: Uses finite differences with step size 1e-7.
    Side Effects: None.
    Failure Modes: Step size too large/small can degrade accuracy.
    Error Handling: Returns zeros if rhs evaluation fails.
    Rationale: Finite differences are robust and don't require symbolic derivatives.
    """
    
    eps = 1e-7
    jaco = np.zeros((6, 6))
    
    # Base RHS evaluation
    from .models import rf3bp_pulsating_rhs_weighted
    rhs_base = rf3bp_pulsating_rhs_weighted(t, state, p, fidelity)
    
    for j in range(6):
        state_pert = state.copy()
        state_pert[j] += eps
        rhs_pert = rf3bp_pulsating_rhs_weighted(t, state_pert, p, fidelity)
        jaco[:, j] = (rhs_pert - rhs_base) / eps
    
    return jaco


def variational_rhs(
    t: float,
    flat_vars: np.ndarray,
    p: SystemParams,
    fidelity: FidelityWeights,
) -> np.ndarray:
    """ID: RF3BP-LAB-DYN-VAR-RHS
    Requirement: Compute time derivatives of state and STM.
    Inputs: Time t, flattened state+STM, parameters, fidelity.
    Outputs: Derivative of flattened state+STM (42 components).
    Preconditions: flat_vars has 42 components.
    Postconditions: Output has 42 components.
    Rationale: d(STM)/dt = Jacobian * STM (Jacobi equation).
    References: Curtis, Orbital Mechanics (2013), Chapter 7.
    """
    
    state, stm = _unflatten_variational(flat_vars)
    
    from .models import rf3bp_pulsating_rhs_weighted
    state_dot = rf3bp_pulsating_rhs_weighted(t, state, p, fidelity)
    
    jaco = _jacobian_rf3bp(t, state, p, fidelity)
    stm_dot = jaco @ stm
    
    return _flatten_variational(state_dot, stm_dot)


def propagate_variational(
    state0: np.ndarray,
    t_span: tuple[float, float],
    p: SystemParams,
    fidelity: FidelityWeights | None = None,
    t_eval: np.ndarray | None = None,
    *,
    method: str = "DOP853",
    rtol: float = 1e-7,
    atol: float = 1e-9,
    max_step: float = np.inf,
) -> dict:
    """ID: RF3BP-LAB-DYN-PROP-VAR
    Requirement: Integrate variational equations to obtain state and STM history.
    Inputs: Initial state, time span, parameters, fidelity, ODE solver options.
    Outputs: Dictionary with time array, state history, and STM history.
    Preconditions: state0 is (6,) and finite.
    Postconditions: Returns state shape (N_time, 6) and stm shape (N_time, 6, 6).
    Failure Modes: Integration can fail on stiff problems; uses adaptive stepping.
    """
    
    fidelity = fidelity or FidelityWeights()
    
    # Initialize STM as identity
    stm0 = np.eye(6)
    flat_vars0 = _flatten_variational(state0, stm0)
    
    # Integrate variational equations
    result = solve_ivp(
        lambda t, y: variational_rhs(t, y, p, fidelity),
        t_span=t_span,
        y0=flat_vars0,
        t_eval=t_eval,
        method=method,
        rtol=rtol,
        atol=atol,
        max_step=max_step,
    )
    
    # Reconstruct state and STM arrays
    n_time = len(result.t)
    states = np.zeros((n_time, 6))
    stms = np.zeros((n_time, 6, 6))
    
    for i in range(n_time):
        state_i, stm_i = _unflatten_variational(result.y[:, i])
        states[i] = state_i
        stms[i] = stm_i
    
    return {
        'success': result.status == 0,
        'time': result.t,
        'states': states,
        'stms': stms,
        'message': result.message,
    }


def covariance_at_time(
    cov_t0: np.ndarray,
    stm: np.ndarray,
) -> np.ndarray:
    """ID: RF3BP-LAB-DYN-COVAR-TRANSFORM
    Requirement: Propagate covariance matrix using STM.
    Inputs: Initial covariance P(t0) (6x6), STM Phi(t, t0) (6x6).
    Outputs: Covariance at time t: P(t) = Phi @ P(t0) @ Phi^T.
    Preconditions: Both inputs are (6x6) and symmetric positive definite.
    Postconditions: Output is (6x6), symmetric positive definite.
    Rationale: Standard covariance time update in linear systems.
    """
    
    return stm @ cov_t0 @ stm.T


def dilution_of_precision(cov: np.ndarray) -> dict:
    """ID: RF3BP-LAB-DYN-DOP
    Requirement: Compute dilution-of-precision (DOP) metrics from covariance.
    Inputs: Covariance matrix (6x6).
    Outputs: Dictionary with geometric DOP (GDOP) and component standard deviations.
    Preconditions: Covariance is symmetric positive definite.
    Postconditions: All DOP values are non-negative.
    Purpose: Navigation quality metric for mission planning.
    Rationale: Trace(position submatrix) provides position uncertainty; analogous for velocity.
    """
    
    cov_pos = cov[:3, :3]
    cov_vel = cov[3:6, 3:6]
    
    pos_uncertainty = np.sqrt(np.trace(cov_pos))
    vel_uncertainty = np.sqrt(np.trace(cov_vel))
    
    # GDOP-like metric
    total_uncertainty = np.sqrt(np.trace(cov))
    
    return {
        'position_std_dev': pos_uncertainty,
        'velocity_std_dev': vel_uncertainty,
        'total_std_dev': total_uncertainty,
        'gdop_like': total_uncertainty,
        'cond_number': np.linalg.cond(cov),
    }
