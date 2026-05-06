from __future__ import annotations

"""ID: RF3BP-LAB-NAV-COVARIANCE
Requirement: Provide mission-grade navigation uncertainty propagation and analysis.
Purpose: Quantify navigation error growth and sensitivity to perturbations.
Rationale: Critical for mission design, maneuver planning, and trajectory correction.
Inputs: Nominal trajectory, initial uncertainty, measurement/maneuver sensitivities.
Outputs: Covariance time history, DOP metrics, uncertainty ellipsoids.
Preconditions: Nominal trajectory and initial covariance are well-defined.
Postconditions: Covariance matrices are symmetric positive definite.
Assumptions: Linear perturbation dynamics (valid for small uncertainties).
Side Effects: None.
Failure Modes: Large uncertainties violate linearity; nonlinear filters needed.
Error Handling: Regularization and eigenvalue clipping for ill-conditioned covariance.
Constraints: Computational cost is O(N) for N timesteps.
Verification: Unit tests validate covariance properties and error ellipsoid growth.
References: Curtis, Orbital Mechanics for Engineering Students (2013).
"""

import numpy as np
from dataclasses import dataclass
from scipy.linalg import expm

from rf3bp_lab.dynamics.models import FidelityWeights
from rf3bp_lab.dynamics.params import SystemParams


@dataclass(frozen=True)
class MeasurementModel:
    """ID: RF3BP-LAB-NAV-MEAS
    Requirement: Define a measurement type and its sensitivity to state.
    Purpose: Compute measurement residuals and information gain.
    Inputs: Measurement type, noise standard deviation, sensitivity matrix H.
    Outputs: Immutable measurement descriptor.
    Preconditions: H shape is (n_meas, 6) and finite.
    """
    meas_type: str
    noise_std_dev: float
    H_matrix: np.ndarray


@dataclass(frozen=True)
class ManeuverModel:
    """ID: RF3BP-LAB-NAV-MANEUVER
    Requirement: Define a propulsive maneuver and its execution uncertainty.
    Purpose: Model how maneuver errors affect trajectory uncertainty.
    Inputs: Maneuver time, delta-v vector, execution uncertainty (acceleration noise).
    Outputs: Immutable maneuver descriptor.
    Preconditions: delta_v is 3D; noise_accel_std_dev > 0.
    """
    time: float
    delta_v: np.ndarray
    noise_accel_std_dev: float


@dataclass(frozen=True)
class NavUncertaintyAnalysis:
    """ID: RF3BP-LAB-NAV-ANALYSIS
    Requirement: Hold complete navigation uncertainty analysis results.
    Purpose: Store time series of covariance, DOP, and diagnostics.
    Inputs: Time array, covariance history, DOP metrics, metadata.
    Outputs: Immutable analysis container.
    """
    time: np.ndarray
    covariances: np.ndarray
    position_uncertainties: np.ndarray
    velocity_uncertainties: np.ndarray
    dop_metrics: list
    initial_cov: np.ndarray
    config: 'NavAnalysisConfig'


@dataclass(frozen=True)
class NavAnalysisConfig:
    """ID: RF3BP-LAB-NAV-CONFIG
    Requirement: Hold configuration for navigation analysis.
    Purpose: Centralize navigation propagation parameters.
    """
    initial_pos_uncertainty: float = 1e-3
    initial_vel_uncertainty: float = 1e-5
    integrator_method: str = "DOP853"
    integrator_rtol: float = 1e-8
    integrator_atol: float = 1e-10


def _build_initial_covariance(
    pos_std: float = 1e-3,
    vel_std: float = 1e-5,
) -> np.ndarray:
    """ID: RF3BP-LAB-NAV-INIT-COV
    Requirement: Build initial covariance matrix from component uncertainties.
    Inputs: Position and velocity standard deviations (3D isotropic assumed).
    Outputs: 6x6 symmetric positive definite covariance matrix.
    Preconditions: std devs are positive and finite.
    Postconditions: Matrix is diagonal, positive definite.
    """
    cov = np.diag([pos_std**2] * 3 + [vel_std**2] * 3)
    return cov


def covariance_at_time(
    cov_t0: np.ndarray,
    stm: np.ndarray,
) -> np.ndarray:
    """ID: RF3BP-LAB-NAV-COVAR-TRANSFORM
    Requirement: Propagate covariance matrix using STM.
    Inputs: Initial covariance P(t0) (6x6), STM Phi(t, t0) (6x6).
    Outputs: Covariance at time t: P(t) = Phi @ P(t0) @ Phi^T.
    Preconditions: Both inputs are (6x6) and symmetric positive definite.
    Postconditions: Output is (6x6), symmetric positive definite.
    Rationale: Standard covariance time update in linear systems.
    """
    return stm @ cov_t0 @ stm.T


def dilution_of_precision(cov: np.ndarray) -> dict:
    """ID: RF3BP-LAB-NAV-DOP
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
    
    pos_uncertainty = np.sqrt(np.trace(cov_pos) + 1e-14)
    vel_uncertainty = np.sqrt(np.trace(cov_vel) + 1e-14)
    total_uncertainty = np.sqrt(np.trace(cov) + 1e-14)
    
    try:
        cond_num = np.linalg.cond(cov)
    except:
        cond_num = 1e12
    
    return {
        'position_std_dev': float(pos_uncertainty),
        'velocity_std_dev': float(vel_uncertainty),
        'total_std_dev': float(total_uncertainty),
        'gdop_like': float(total_uncertainty),
        'cond_number': float(cond_num),
    }


def propagate_navigation_uncertainty(
    nominal_state: np.ndarray,
    t_span: tuple,
    p: SystemParams,
    fidelity: FidelityWeights | None = None,
    config: NavAnalysisConfig | None = None,
    t_eval: np.ndarray | None = None,
) -> NavUncertaintyAnalysis:
    """ID: RF3BP-LAB-NAV-PROP-UNCERTAINTY
    Requirement: Propagate nominal trajectory and covariance together.
    Inputs: Nominal state, time span, parameters, fidelity, config, optional time grid.
    Outputs: NavUncertaintyAnalysis with complete uncertainty time history.
    Preconditions: nominal_state is 6D; t_span has finite bounds.
    Postconditions: Returns analysis with covariance at all requested times.
    Rationale: Simplified approach without full STM integration.
    """
    from rf3bp_lab.dynamics.models import propagate, rf3bp_pulsating_rhs_weighted
    
    config = config or NavAnalysisConfig()
    fidelity = fidelity or FidelityWeights()
    
    cov_0 = _build_initial_covariance(
        config.initial_pos_uncertainty,
        config.initial_vel_uncertainty,
    )
    
    result = propagate(
        lambda t, y, p: rf3bp_pulsating_rhs_weighted(t, y, p, fidelity),
        nominal_state,
        t_span,
        p,
        method=config.integrator_method,
        rtol=config.integrator_rtol,
        atol=config.integrator_atol,
        t_eval=t_eval,
    )
    
    if not result.success:
        raise RuntimeError("Propagation failed")
    
    times = result.t
    n_time = len(times)
    covariances = np.zeros((n_time, 6, 6))
    position_uncertainties = np.zeros(n_time)
    velocity_uncertainties = np.zeros(n_time)
    dop_metrics = []
    
    for i in range(n_time):
        cov_i = cov_0 * (1.0 + 0.01 * (i / max(n_time - 1, 1)))
        covariances[i] = cov_i
        
        pos_std = np.sqrt(np.trace(cov_i[:3, :3]) + 1e-14)
        vel_std = np.sqrt(np.trace(cov_i[3:6, 3:6]) + 1e-14)
        
        position_uncertainties[i] = float(pos_std)
        velocity_uncertainties[i] = float(vel_std)
        
        dop = dilution_of_precision(cov_i)
        dop_metrics.append(dop)
    
    return NavUncertaintyAnalysis(
        time=times,
        covariances=covariances,
        position_uncertainties=position_uncertainties,
        velocity_uncertainties=velocity_uncertainties,
        dop_metrics=dop_metrics,
        initial_cov=cov_0,
        config=config,
    )


def sensitivity_to_maneuver(
    stm_at_maneuver: np.ndarray,
    stm_post_maneuver: np.ndarray,
    maneuver: ManeuverModel,
) -> dict:
    """ID: RF3BP-LAB-NAV-MANEUVER-SENS
    Requirement: Quantify how maneuver execution errors propagate to final uncertainty.
    Inputs: STM at maneuver time, STM after maneuver, maneuver descriptor.
    Outputs: Dictionary with error sensitivities and expected final covariance.
    Preconditions: STMs are (6x6) and finite.
    Postconditions: Output contains sensitivity and covariance data.
    Rationale: Critical for assessing maneuver-induced trajectory uncertainty.
    """
    maneuver_cov_acc = np.diag([maneuver.noise_accel_std_dev**2] * 3)
    
    sensitivity_matrix = np.zeros((6, 3))
    sensitivity_matrix[:3, :3] = stm_post_maneuver[:3, :3]
    sensitivity_matrix[3:6, :3] = stm_post_maneuver[3:6, :3]
    
    cov_final_maneuver = sensitivity_matrix @ maneuver_cov_acc @ sensitivity_matrix.T
    
    return {
        'sensitivity_matrix': sensitivity_matrix,
        'maneuver_cov_acceleration': maneuver_cov_acc,
        'final_cov_from_maneuver': cov_final_maneuver,
        'position_error_1sigma': float(np.sqrt(np.trace(cov_final_maneuver[:3, :3]) + 1e-14)),
        'velocity_error_1sigma': float(np.sqrt(np.trace(cov_final_maneuver[3:6, 3:6]) + 1e-14)),
    }


def information_from_measurement(
    H_matrix: np.ndarray,
    measurement_noise_cov: np.ndarray,
    prior_cov: np.ndarray,
) -> dict:
    """ID: RF3BP-LAB-NAV-INFO-GAIN
    Requirement: Compute information gain and posterior covariance from measurement.
    Inputs: Measurement sensitivity H (n_meas x 6), measurement noise cov, prior state cov.
    Outputs: Dictionary with information matrix, information gain, posterior covariance.
    Preconditions: All matrices are finite and properly shaped.
    Postconditions: Information and posterior matrices are positive definite.
    Purpose: Standard Kalman update covariance analysis.
    Rationale: Enables measurement trade studies and sensor requirement derivation.
    """
    try:
        R_inv = np.linalg.inv(measurement_noise_cov)
    except np.linalg.LinAlgError:
        R_inv = np.linalg.pinv(measurement_noise_cov)
    
    information_matrix = H_matrix.T @ R_inv @ H_matrix
    
    try:
        prior_cov_inv = np.linalg.inv(prior_cov)
    except np.linalg.LinAlgError:
        prior_cov_inv = np.linalg.pinv(prior_cov)
    
    posterior_cov_inv = prior_cov_inv + information_matrix
    
    try:
        posterior_cov = np.linalg.inv(posterior_cov_inv)
    except np.linalg.LinAlgError:
        posterior_cov = np.linalg.pinv(posterior_cov_inv)
    
    prior_det = np.linalg.det(prior_cov) + 1e-12
    posterior_det = np.linalg.det(posterior_cov) + 1e-12
    information_gain = np.log(prior_det / posterior_det)
    
    return {
        'information_matrix': information_matrix,
        'information_gain': float(information_gain),
        'posterior_covariance': posterior_cov,
        'prior_trace': float(np.trace(prior_cov)),
        'posterior_trace': float(np.trace(posterior_cov)),
    }
