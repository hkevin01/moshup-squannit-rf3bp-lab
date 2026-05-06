from __future__ import annotations

"""ID: RF3BP-LAB-DYN-MODELS
Requirement: Provide deterministic CR3BP and RF3BP-inspired dynamics building blocks.
Purpose: Support bounded-orbit experimentation near binary-asteroid systems.
Rationale: Keep physics terms decomposed for continuation and diagnostics.
Inputs: Normalized state/time and SystemParams.
Outputs: State derivatives, perturbation breakdowns, and propagation results.
Preconditions: Finite numeric inputs and physically meaningful parameter scales.
Postconditions: Returned arrays have fixed shapes and finite values for nominal inputs.
Assumptions: Rotating-frame normalization and simplified perturbation proxies.
Side Effects: None.
Failure Modes: Ill-conditioned geometry can amplify roundoff near singular distances.
Error Handling: Softening epsilons and bounds checks prevent non-finite operations.
Constraints: Lightweight runtime for iterative shooting workflows.
Verification: Unit tests in tests/test_dynamics.py and tests/test_shooting.py.
References: DOI 10.2514/1.G009686, DOI 10.1016/0771-050X(80)90013-3.
"""

from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp

from .params import SystemParams


@dataclass(frozen=True)
class FidelityWeights:
    """ID: RF3BP-LAB-DYN-FIDELITY
    Requirement: Scale perturbation families continuously from 0 to 1.
    Purpose: Enable hierarchical continuation between model fidelities.
    Inputs: Scalar weights for pulsation, nonspherical gravity, solar gravity, and SRP.
    Outputs: Immutable weight container.
    Preconditions: Weights are finite scalars.
    Postconditions: Values remain unchanged after construction.
    Failure Modes: Non-finite values can destabilize residual correction.
    """

    pulsation: float = 1.0
    nonspherical: float = 1.0
    solar_gravity: float = 1.0
    srp: float = 1.0


@dataclass(frozen=True)
class RelativeKinematics:
    """ID: RF3BP-LAB-DYN-RELKIN
    Requirement: Carry relative secondary-body position, velocity, acceleration, and jerk.
    Purpose: Expose pulsation-driven kinematic terms for diagnostics and dynamics.
    """

    r_vec: np.ndarray
    v_vec: np.ndarray
    a_vec: np.ndarray
    j_vec: np.ndarray


@dataclass(frozen=True)
class PerturbationBreakdown:
    """ID: RF3BP-LAB-DYN-BREAKDOWN
    Requirement: Store additive acceleration terms used by the weighted RF3BP rhs.
    Purpose: Support perturbation ranking and sensitivity analysis.
    """

    grav_primary: np.ndarray
    grav_secondary: np.ndarray
    coriolis: np.ndarray
    centrifugal: np.ndarray
    pulsation: np.ndarray
    solar_gravity: np.ndarray
    srp: np.ndarray


@dataclass(frozen=True)
class DynamicsGap:
    """ID: RF3BP-LAB-DYN-GAP
    Requirement: Capture vector/scalar acceleration differences between RF3BP and CR3BP.
    Purpose: Quantify fidelity gap for diagnostics, continuation tuning, and plotting.
    Inputs: Time t, state, and model evaluations.
    Outputs: CR3BP acceleration, RF3BP acceleration, delta vector, and relative norm.
    Preconditions: Acceleration vectors are finite 3-vectors.
    Postconditions: relative_gap is non-negative.
    Failure Modes: Near-zero CR3BP norm can destabilize relative metrics.
    Error Handling: Uses epsilon denominator floor.
    """

    cr3bp_acc: np.ndarray
    rf3bp_acc: np.ndarray
    delta_acc: np.ndarray
    delta_norm: float
    relative_gap: float


def _r12_scale(t: float, p: SystemParams) -> tuple[float, float, float]:
    """ID: RF3BP-LAB-DYN-R12SCALE
    Requirement: Compute normalized separation scale and its first two derivatives.
    Inputs: Time t, system parameters p.
    Outputs: (r12, r12_dot, r12_ddot).
    Preconditions: p.pulsation_nu and p.pulsation_e are finite.
    Postconditions: Returns finite scalars for nominal parameter ranges.
    """

    c = np.cos(p.pulsation_nu * t)
    s = np.sin(p.pulsation_nu * t)
    r = 1.0 + p.pulsation_e * c
    r_dot = -p.pulsation_e * p.pulsation_nu * s
    r_ddot = -p.pulsation_e * (p.pulsation_nu**2) * c
    return r, r_dot, r_ddot


def _j2_accel(rel: np.ndarray, mu_body: float, j2: float, r_body_scaled: float) -> np.ndarray:
    """ID: RF3BP-LAB-DYN-J2
    Requirement: Return point-mass plus J2-style acceleration for one gravitating body.
    Inputs: Relative vector rel, body mass fraction mu_body, J2 proxy, scaled radius.
    Outputs: 3-vector acceleration.
    Failure Modes: Near-zero radius can cause singular growth without regularization.
    Error Handling: Adds epsilon to squared distance.
    """

    x, y, z = rel
    r2 = np.dot(rel, rel) + 1e-12
    r = np.sqrt(r2)
    fac = -mu_body / (r2 * r)
    a_pm = fac * rel

    z2_r2 = (z * z) / r2
    k = 1.5 * j2 * mu_body * (r_body_scaled**2) / (r**5)
    ax = k * x * (5.0 * z2_r2 - 1.0)
    ay = k * y * (5.0 * z2_r2 - 1.0)
    az = k * z * (5.0 * z2_r2 - 3.0)
    return a_pm + np.array([ax, ay, az])


def _solar_terms(t: float, pos: np.ndarray, p: SystemParams) -> tuple[np.ndarray, np.ndarray]:
    """ID: RF3BP-LAB-DYN-SOLAR
    Requirement: Compute differential solar gravity and SRP proxies.
    Inputs: Time t, spacecraft position pos, system parameters p.
    Outputs: Tuple (a_sun, a_srp).
    Assumptions: Simplified planar Sun motion and constant SRP magnitude scaling.
    """

    sun_angle = 0.03 * t
    sun = p.sun_distance_scaled * np.array([np.cos(sun_angle), np.sin(sun_angle), 0.0])
    rel_sc = sun - pos
    rel_bary = sun

    d_sc = np.linalg.norm(rel_sc) + 1e-12
    d_b = np.linalg.norm(rel_bary) + 1e-12
    a_sun = p.sun_mu_scaled * (rel_sc / d_sc**3 - rel_bary / d_b**3)

    srp_dir = -rel_sc / d_sc
    a_srp = p.srp_accel_scaled * srp_dir
    return a_sun, a_srp


def _secondary_relative_kinematics_from_potential(t: float, p: SystemParams) -> RelativeKinematics:
    """ID: RF3BP-LAB-DYN-SECKIN
    Requirement: Derive relative acceleration and jerk from a two-body potential Hessian.
    Inputs: Time t and system parameters p.
    Outputs: RelativeKinematics dataclass.
    Rationale: Exposes pulsation-driven inertial terms with explicit derivatives.
    """

    r12, r12_dot, _ = _r12_scale(t, p)
    r_vec = np.array([r12, 0.0, 0.0], dtype=float)
    v_vec = np.array([r12_dot, 0.0, 0.0], dtype=float)

    # Two-body potential U = -mu / r with normalized mu = 1.
    r = np.linalg.norm(r_vec) + 1e-12
    rr = np.outer(r_vec, r_vec)
    identity = np.eye(3)

    a_vec = -(1.0 / r**3) * r_vec
    hess_u = (identity / r**3) - 3.0 * rr / r**5
    j_vec = -hess_u @ v_vec

    return RelativeKinematics(r_vec=r_vec, v_vec=v_vec, a_vec=a_vec, j_vec=j_vec)


def cr3bp_rhs(t: float, state: np.ndarray, p: SystemParams) -> np.ndarray:
    """ID: RF3BP-LAB-DYN-CR3BP-RHS
    Requirement: Evaluate the normalized rotating-frame CR3BP state derivative.
    Inputs: Time t, state [x,y,z,vx,vy,vz], parameters p.
    Outputs: 6-element derivative vector.
    Preconditions: state has 6 finite components.
    Postconditions: Output shape is (6,).
    """

    x, y, z, vx, vy, vz = state
    mu = p.mu

    r1_vec = np.array([x + mu, y, z])
    r2_vec = np.array([x - (1.0 - mu), y, z])
    r1 = np.linalg.norm(r1_vec) + 1e-12
    r2 = np.linalg.norm(r2_vec) + 1e-12

    ddx = 2.0 * vy + x - (1.0 - mu) * (x + mu) / r1**3 - mu * (x - (1.0 - mu)) / r2**3
    ddy = -2.0 * vx + y - (1.0 - mu) * y / r1**3 - mu * y / r2**3
    ddz = -(1.0 - mu) * z / r1**3 - mu * z / r2**3

    return np.array([vx, vy, vz, ddx, ddy, ddz], dtype=float)


def rf3bp_breakdown(t: float, state: np.ndarray, p: SystemParams, fidelity: FidelityWeights | None = None) -> PerturbationBreakdown:
    """ID: RF3BP-LAB-DYN-RF3BP-BREAKDOWN
    Requirement: Compute all additive acceleration components for the RF3BP-inspired model.
    Inputs: Time t, state vector, parameters p, optional fidelity weights.
    Outputs: PerturbationBreakdown with 7 acceleration components.
    Constraints: Designed for repeated calls inside continuation and plotting loops.
    """

    fidelity = fidelity or FidelityWeights()

    x, y, z, vx, vy, vz = state
    mu1 = 1.0 - p.mu
    mu2 = p.mu

    rel = _secondary_relative_kinematics_from_potential(t, p)
    r12 = np.linalg.norm(rel.r_vec)
    r12_dot = rel.v_vec[0]
    r12_acc = rel.a_vec[0]

    primary = np.array([-p.mu * r12, 0.0, 0.0])
    secondary = np.array([(1.0 - p.mu) * r12, 0.0, 0.0])
    pos = np.array([x, y, z])

    rel1 = pos - primary
    rel2 = pos - secondary

    rp = p.r_primary_m / p.r12_mean_m
    rs = p.r_secondary_m / p.r12_mean_m

    r1 = np.linalg.norm(rel1) + 1e-12
    r2 = np.linalg.norm(rel2) + 1e-12
    a_pm1 = -(mu1 / r1**3) * rel1
    a_pm2 = -(mu2 / r2**3) * rel2

    a_j2_1 = _j2_accel(rel1, mu1, p.j2_primary, rp) - a_pm1
    a_j2_2 = _j2_accel(rel2, mu2, p.j2_secondary, rs) - a_pm2

    a_g1 = a_pm1 + fidelity.nonspherical * a_j2_1
    a_g2 = a_pm2 + fidelity.nonspherical * a_j2_2

    coriolis = np.array([2.0 * p.omega0 * vy, -2.0 * p.omega0 * vx, 0.0])
    centrifugal = np.array([p.omega0**2 * x, p.omega0**2 * y, 0.0])

    pulsation_acc = -2.0 * (r12_dot / max(r12, 1e-12)) * np.array([vx, vy, vz]) - (r12_acc / max(r12, 1e-12)) * pos
    pulsation_acc *= fidelity.pulsation

    a_sun, a_srp = _solar_terms(t, pos, p)
    a_sun *= fidelity.solar_gravity
    a_srp *= fidelity.srp

    return PerturbationBreakdown(
        grav_primary=a_g1,
        grav_secondary=a_g2,
        coriolis=coriolis,
        centrifugal=centrifugal,
        pulsation=pulsation_acc,
        solar_gravity=a_sun,
        srp=a_srp,
    )


def rf3bp_pulsating_rhs_weighted(
    t: float,
    state: np.ndarray,
    p: SystemParams,
    fidelity: FidelityWeights,
) -> np.ndarray:
    """ID: RF3BP-LAB-DYN-RF3BP-WRHS
    Requirement: Build full state derivative by summing weighted perturbation terms.
    Inputs: Time t, state, parameters p, fidelity weights.
    Outputs: 6-element derivative vector.
    Return-Path: Single return statement with deterministic assembly.
    """

    x, y, z, vx, vy, vz = state
    b = rf3bp_breakdown(t, state, p, fidelity)
    acc = b.grav_primary + b.grav_secondary + b.coriolis + b.centrifugal + b.pulsation + b.solar_gravity + b.srp
    return np.array([vx, vy, vz, acc[0], acc[1], acc[2]], dtype=float)


def rf3bp_pulsating_rhs(t: float, state: np.ndarray, p: SystemParams) -> np.ndarray:
    """ID: RF3BP-LAB-DYN-RF3BP-RHS
    Requirement: Evaluate weighted RF3BP rhs using full-unity fidelity weights.
    """

    return rf3bp_pulsating_rhs_weighted(t, state, p, FidelityWeights())


def secondary_relative_kinematics(t: float, p: SystemParams) -> RelativeKinematics:
    """ID: RF3BP-LAB-DYN-SECKIN-PUBLIC
    Requirement: Expose secondary relative kinematics to callers and tests.
    """

    return _secondary_relative_kinematics_from_potential(t, p)


def cr3bp_effective_potential(state: np.ndarray, p: SystemParams) -> float:
    """ID: RF3BP-LAB-DYN-CR3BP-OMEGA
    Requirement: Compute rotating-frame CR3BP effective potential Omega(x,y,z).
    Purpose: Provide scalar diagnostics for equilibrium and energy analyses.
    Inputs: State vector with at least first three coordinates and parameters p.
    Outputs: Scalar potential value.
    Preconditions: State has finite position components.
    Postconditions: Finite scalar for non-collision positions.
    Failure Modes: Near-primary/secondary positions create singular growth.
    Error Handling: Distance softening epsilon protects against division-by-zero.
    References: Szebehely, Theory of Orbits (1967).
    """

    x, y, z = state[:3]
    mu = p.mu
    r1 = np.linalg.norm(np.array([x + mu, y, z])) + 1e-12
    r2 = np.linalg.norm(np.array([x - (1.0 - mu), y, z])) + 1e-12
    return 0.5 * (x * x + y * y) + (1.0 - mu) / r1 + mu / r2


def cr3bp_jacobi_constant(state: np.ndarray, p: SystemParams) -> float:
    """ID: RF3BP-LAB-DYN-CR3BP-JACOBI
    Requirement: Compute Jacobi constant C = 2*Omega - v^2 for CR3BP states.
    Purpose: Enable trajectory quality checks and phase-space diagnostics.
    Inputs: Full state vector [x,y,z,vx,vy,vz], parameters p.
    Outputs: Scalar Jacobi constant.
    Preconditions: State has 6 finite components.
    Postconditions: Finite scalar for non-collision states.
    Verification: Covered by unit tests that assert finite values.
    References: CR3BP standard integral from rotating-frame equations.
    """

    vx, vy, vz = state[3:6]
    return 2.0 * cr3bp_effective_potential(state, p) - (vx * vx + vy * vy + vz * vz)


def compare_cr3bp_rf3bp(
    t: float,
    state: np.ndarray,
    p: SystemParams,
    fidelity: FidelityWeights | None = None,
) -> DynamicsGap:
    """ID: RF3BP-LAB-DYN-COMPARE
    Requirement: Compute acceleration-level gap between CR3BP and weighted RF3BP dynamics.
    Purpose: Provide a direct mathematical answer to "how different" the two models are.
    Rationale: Orbit-design decisions depend on local acceleration mismatch magnitudes.
    Inputs: Time t, 6D state, parameters p, optional fidelity weights.
    Outputs: DynamicsGap dataclass with absolute and relative gap metrics.
    Preconditions: State has at least 6 finite elements.
    Postconditions: delta_norm >= 0 and relative_gap >= 0.
    Assumptions: Comparison is done in the same normalized rotating frame.
    Side Effects: None.
    Failure Modes: Singular configurations near primaries can inflate norms.
    Error Handling: Relative denominator uses max(norm, 1e-12).
    Constraints: Lightweight enough for per-step history plots.
    Verification: tests/test_dynamics.py coverage.
    References: DOI 10.2514/1.G009686, Szebehely (1967).
    """

    fw = fidelity or FidelityWeights()
    a_cr3bp = cr3bp_rhs(t, state, p)[3:6]
    a_rf3bp = rf3bp_pulsating_rhs_weighted(t, state, p, fw)[3:6]
    delta = a_rf3bp - a_cr3bp
    delta_norm = float(np.linalg.norm(delta))
    base_norm = max(float(np.linalg.norm(a_cr3bp)), 1e-12)
    relative_gap = delta_norm / base_norm
    return DynamicsGap(
        cr3bp_acc=a_cr3bp,
        rf3bp_acc=a_rf3bp,
        delta_acc=delta,
        delta_norm=delta_norm,
        relative_gap=relative_gap,
    )


def propagate(
    rhs,
    state0: np.ndarray,
    t_span: tuple[float, float],
    p: SystemParams,
    t_eval: np.ndarray | None = None,
    *,
    method: str = "DOP853",
    rtol: float = 1e-10,
    atol: float = 1e-12,
    max_step: float = np.inf,
):
    """ID: RF3BP-LAB-DYN-PROPAGATE
    Requirement: Integrate a parameterized rhs over a time span.
    Inputs: rhs(t, y, p), initial state, time span, params, solver controls.
    Outputs: scipy.integrate.OdeResult.
    Preconditions: rhs is callable and state dimensions are consistent.
    Postconditions: Returns solver status and sampled state history.
    Failure Modes: Integrator may fail on stiff/ill-scaled trajectories.
    Error Handling: Caller checks solve_ivp success flag.
    References: DOI 10.1016/0771-050X(80)90013-3.
    """

    return solve_ivp(
        lambda t, y: rhs(t, y, p),
        t_span=t_span,
        y0=state0,
        t_eval=t_eval,
        method=method,
        rtol=rtol,
        atol=atol,
        max_step=max_step,
    )