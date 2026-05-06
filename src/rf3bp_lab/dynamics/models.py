from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp

from .params import SystemParams


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


def rf3bp_pulsating_rhs(t: float, state: np.ndarray, p: SystemParams) -> np.ndarray:
    x, y, z, vx, vy, vz = state
    mu1 = 1.0 - p.mu
    mu2 = p.mu

    r12, r12_dot, r12_ddot = _r12_scale(t, p)

    primary = np.array([-p.mu * r12, 0.0, 0.0])
    secondary = np.array([(1.0 - p.mu) * r12, 0.0, 0.0])
    pos = np.array([x, y, z])

    rel1 = pos - primary
    rel2 = pos - secondary

    rp = p.r_primary_m / p.r12_mean_m
    rs = p.r_secondary_m / p.r12_mean_m

    a_g1 = _j2_accel(rel1, mu1, p.j2_primary, rp)
    a_g2 = _j2_accel(rel2, mu2, p.j2_secondary, rs)

    # Rotating-pulsating inertial terms
    coriolis = np.array([2.0 * p.omega0 * vy, -2.0 * p.omega0 * vx, 0.0])
    centrifugal = np.array([p.omega0**2 * x, p.omega0**2 * y, 0.0])

    # Nonuniform pulsation correction from scale acceleration + velocity coupling
    pulsation_acc = -2.0 * (r12_dot / max(r12, 1e-12)) * np.array([vx, vy, vz]) - (r12_ddot / max(r12, 1e-12)) * pos

    a_sun, a_srp = _solar_terms(t, pos, p)

    acc = a_g1 + a_g2 + coriolis + centrifugal + pulsation_acc + a_sun + a_srp
    return np.array([vx, vy, vz, acc[0], acc[1], acc[2]], dtype=float)


def propagate(rhs, state0: np.ndarray, t_span: tuple[float, float], p: SystemParams, t_eval: np.ndarray | None = None):
    return solve_ivp(lambda t, y: rhs(t, y, p), t_span=t_span, y0=state0, t_eval=t_eval, rtol=1e-10, atol=1e-12)
