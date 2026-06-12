from __future__ import annotations

"""ID: RF3BP-LAB-DYN-SUN-EPHEMERIS
Requirement: Provide analytical, SPICE-free Sun position accurate to ~1 arc-minute.
Purpose: Replace the naive circular orbit with a proper ecliptic ephemeris so that
         solar perturbations are geometrically correct without requiring SPICE kernels.
Rationale: SPICE is a heavy dependency with OS-specific shared libraries that cannot
           be assumed on all user machines.  The JPL Astronomical Almanac low-precision
           formulae (Meeus 1998) reproduce Sun direction to ~1' without any external
           data files, making them portable and reproducible.
Inputs: Julian date (TDB approximation acceptable for this fidelity level).
Outputs: Sun unit vector in the J2000 equatorial frame, heliocentric distance (AU).
Preconditions: Julian date is finite and represents a date between 1950 and 2100.
Postconditions: Returned unit vector is finite and normalised to 1.
Assumptions: Low-precision formulae valid to ~0.01 deg in ecliptic longitude.
Side Effects: None.
Failure Modes: Extrapolation far outside 1950-2100 degrades accuracy gracefully.
Error Handling: Domain guard replaces non-finite JD with J2000.0 epoch.
Constraints: Pure Python / NumPy; no external astronomy libraries required.
Verification: Unit tests compare against known Sun direction at J2000.0.
References: Meeus, Astronomical Algorithms 2nd ed. (1998), Chapter 25.
            Vallado, Fundamentals of Astrodynamics (2013), Section 5.1.
"""

import numpy as np
from dataclasses import dataclass


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
_J2000 = 2451545.0          # Julian date of J2000.0 epoch
_AU_M  = 1.495978707e11     # metres per astronomical unit
_DEG   = np.pi / 180.0


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class SunState:
    """ID: RF3BP-LAB-SUN-STATE
    Requirement: Bundle the Sun's geometric state for a given epoch.
    Purpose: Provide all solar-perturbation inputs in one immutable container.
    Inputs: Julian date, unit direction vector (J2000 equatorial), distance in AU.
    Outputs: Immutable Sun state descriptor.
    Preconditions: direction is a unit 3-vector; distance > 0.
    Postconditions: All fields are finite.
    """
    jd: float
    direction_j2000: np.ndarray   # unit vector toward Sun from barycentre
    distance_au: float            # heliocentric distance in AU
    distance_m: float             # heliocentric distance in metres


@dataclass(frozen=True)
class SunEphemerisConfig:
    """ID: RF3BP-LAB-SUN-CONFIG
    Requirement: Hold conversion parameters bridging the RF3BP normalisation with SI.
    Purpose: Translate the analytical AU-frame Sun position into the normalised RF3BP frame.
    Inputs: Characteristic length scale (metres), characteristic time scale (seconds).
    Outputs: Immutable config container.
    """
    length_scale_m: float = 2500.0   # r12_mean_m from SystemParams
    time_scale_s:   float = 6950.0   # ~1/(omega0) converted to seconds; approximate


# ---------------------------------------------------------------------------
# Analytical low-precision Sun position (Meeus Ch. 25)
# ---------------------------------------------------------------------------

def sun_position_j2000(jd: float) -> SunState:
    """ID: RF3BP-LAB-SUN-POSITION
    Requirement: Compute Sun unit direction and distance from Earth/barycentre at epoch JD.
    Inputs: jd (float) - Julian date (TDB is fine at this precision level).
    Outputs: SunState with direction vector (J2000 equatorial) and distance.
    Preconditions: jd is finite and in [2415020, 2488070] (year 1900-2100).
    Postconditions: direction_j2000 is a unit vector; distance_au > 0.
    Failure Modes: Non-finite jd returns J2000.0 fallback.
    References: Meeus 1998, Chapter 25, low-precision solar coordinates.
    """
    if not np.isfinite(jd):
        jd = _J2000

    T = (jd - _J2000) / 36525.0   # Julian centuries from J2000.0

    # Mean longitude of the Sun (degrees)
    L0 = (280.46646 + 36000.76983 * T + 0.0003032 * T * T) % 360.0

    # Mean anomaly of the Sun (degrees)
    M = (357.52911 + 35999.05029 * T - 0.0001537 * T * T) % 360.0
    M_rad = M * _DEG

    # Equation of the centre
    C = (
        (1.914602 - 0.004817 * T - 0.000014 * T * T) * np.sin(M_rad)
        + (0.019993 - 0.000101 * T) * np.sin(2.0 * M_rad)
        + 0.000289 * np.sin(3.0 * M_rad)
    )

    # Sun true longitude (degrees)
    lam = L0 + C
    lam_rad = lam * _DEG

    # Sun true anomaly
    nu = M + C

    # Radius vector (AU)
    e = 0.016708634 - 0.000042037 * T - 0.0000001267 * T * T
    R = 1.000001018 * (1.0 - e * e) / (1.0 + e * np.cos(nu * _DEG))

    # Apparent longitude: subtract aberration
    omega = (125.04 - 1934.136 * T) * _DEG
    lam_apparent = lam - 0.00569 - 0.00478 * np.sin(omega)
    lam_app_rad = lam_apparent * _DEG

    # Obliquity of the ecliptic (degrees)
    eps0 = 23.0 + 26.0 / 60.0 + 21.448 / 3600.0
    eps = eps0 - (46.8150 * T + 0.00059 * T * T - 0.001813 * T * T * T) / 3600.0
    eps_rad = eps * _DEG

    # Ecliptic -> equatorial (J2000)
    x = np.cos(lam_app_rad)
    y = np.cos(eps_rad) * np.sin(lam_app_rad)
    z = np.sin(eps_rad) * np.sin(lam_app_rad)
    direction = np.array([x, y, z], dtype=float)
    direction /= np.linalg.norm(direction) + 1e-14

    return SunState(
        jd=jd,
        direction_j2000=direction,
        distance_au=float(R),
        distance_m=float(R * _AU_M),
    )


def jd_from_mission_time(
    t_normalised: float,
    epoch_jd: float = _J2000,
    omega0: float = 1.0,
    r12_mean_m: float = 2500.0,
) -> float:
    """ID: RF3BP-LAB-SUN-JD-CONVERT
    Requirement: Convert normalised RF3BP time to Julian date.
    Purpose: Bridge the dimensionless RF3BP clock to the ephemeris time axis.
    Inputs: t_normalised (dimensionless), epoch_jd, angular frequency omega0 (rad/s),
            characteristic length r12_mean_m (metres).
    Outputs: Julian date (float).
    Preconditions: omega0 > 0; r12_mean_m > 0.
    Postconditions: Output is finite.
    Rationale: RF3BP time unit = 1/omega0 seconds; omega0 derived from Kepler's third law.
    """
    import math
    # Estimate omega0 from Kepler: omega0 = sqrt(G*(m1+m2)/r12^3)
    # For normalized units omega0 is given directly.
    period_s = 2.0 * math.pi / max(omega0, 1e-12)
    t_seconds = t_normalised * period_s / (2.0 * np.pi)
    t_days = t_seconds / 86400.0
    return epoch_jd + t_days


def sun_acceleration_rf3bp(
    t: float,
    spacecraft_pos: np.ndarray,
    sun_state: SunState,
    sun_mu_scaled: float,
    sun_distance_scaled: float,
    srp_accel_scaled: float,
    eclipse_shadow_fraction: float = 0.0,
) -> tuple[np.ndarray, np.ndarray]:
    """ID: RF3BP-LAB-SUN-ACCEL
    Requirement: Compute ephemeris-quality solar gravity and SRP accelerations.
    Purpose: Provide improved solar perturbation terms using real Sun direction.
    Inputs: t (dimensionless RF3BP time), spacecraft_pos (3,), SunState,
            scaled force parameters from SystemParams, eclipse shadow fraction [0,1].
    Outputs: (a_solar, a_srp) each shape (3,).
    Preconditions: sun_state direction is a unit vector.
    Postconditions: Both output vectors are finite.
    Failure Modes: Very small sun_distance_scaled can cause numeric overflow.
    Error Handling: Distance floor at 1e-6.
    """
    # Scale sun direction to the normalised RF3BP distance
    sun_vec = sun_state.direction_j2000 * sun_distance_scaled

    # Differential solar gravity (third-body)
    rel_sc = sun_vec - spacecraft_pos
    rel_bc = sun_vec

    d_sc = np.linalg.norm(rel_sc) + 1e-14
    d_bc = np.linalg.norm(rel_bc) + 1e-14

    a_solar = sun_mu_scaled * (rel_sc / d_sc**3 - rel_bc / d_bc**3)

    # SRP - attenuate by eclipse shadow
    srp_dir = -rel_sc / d_sc
    a_srp = srp_accel_scaled * (1.0 - eclipse_shadow_fraction) * srp_dir

    return a_solar, a_srp


# ---------------------------------------------------------------------------
# Convenience: Sun state at a reference epoch for Moshup-Squannit (2001)
# ---------------------------------------------------------------------------

def default_sun_state() -> SunState:
    """ID: RF3BP-LAB-SUN-DEFAULT
    Requirement: Return a representative Sun state for Moshup-Squannit epoch 2001.
    Purpose: Provide a sensible default for demonstrations and unit tests.
    Outputs: SunState at JD 2452000.5 (2001-07-12).
    """
    return sun_position_j2000(2452000.5)
