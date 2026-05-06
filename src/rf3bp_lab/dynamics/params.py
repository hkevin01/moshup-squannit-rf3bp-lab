from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class SystemParams:
    # Normalized RF3BP parameters.
    mu: float = 0.02
    # Mean separation in meters (approximate for KW4/Moshup-Squannit)
    r12_mean_m: float = 2500.0
    # Pulsation model r12(t) = r0 * (1 + e*cos(nu*t))
    pulsation_e: float = 0.08
    pulsation_nu: float = 0.35
    # Nonspherical gravity via simple J2-style term
    j2_primary: float = 0.05
    j2_secondary: float = 0.02
    r_primary_m: float = 700.0
    r_secondary_m: float = 225.0
    # Third-body solar perturbation and SRP in normalized units
    sun_mu_scaled: float = 5.0e-4
    sun_distance_scaled: float = 2000.0
    srp_accel_scaled: float = 2.5e-6
    # Scale factors
    omega0: float = 1.0
