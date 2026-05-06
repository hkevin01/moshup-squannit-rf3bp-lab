from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp

from .params import SystemParams


@dataclass(frozen=True)
class FidelityWeights:
    pulsation: float = 1.0
    nonspherical: float = 1.0
    solar_gravity: float = 1.0
    srp: float = 1.0


@dataclass(frozen=True)
class RelativeKinematics:
    r_vec: np.ndarray
    v_vec: np.ndarray
    a_vec: np.ndarray
    j_vec: np.ndarray


@dataclass(frozen=True)
class PerturbationBreakdown:
    grav_primary: np.ndarray
    grav_secondary: np.ndarray
    coriolis: np.ndarray
    centrifugal: np.ndarray
    pulsation: np.ndarray
    solar_gravity: np.ndarray
    srp: np.ndarray


def _r12_scale(t: float, p: SystemParams) -> tuple[float, float, float]:
    c = np.cos(p.pulsation_nu * t)
    s = np.sin(p.pulsation_nu * t)
    r = 1.0 + p.pulsation_e * c
    r_dot = -p.pulsation_e * p.pulsation_nu * s
    r_ddot = -p.pulsation_e * (p.pulsation_nu**2) * c
    return r, r_dot, r_ddot


def _j2_accel(rel: np.ndarray, mu_body: float, j2: float, r_body_scaled: float) -> np.ndarray:
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
    x, y, z, vx, vy, vz = state
    b = rf3bp_breakdown(t, state, p, fidelity)
    acc = b.grav_primary + b.grav_secondary + b.coriolis + b.centrifugal + b.pulsation + b.solar_gravity + b.srp
    return np.array([vx, vy, vz, acc[0], acc[1], acc[2]], dtype=float)


def rf3bp_pulsating_rhs(t: float, state: np.ndarray, p: SystemParams) -> np.ndarray:
    return rf3bp_pulsating_rhs_weighted(t, state, p, FidelityWeights())


def secondary_relative_kinematics(t: float, p: SystemParams) -> RelativeKinematics:
    return _secondary_relative_kinematics_from_potential(t, p)


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