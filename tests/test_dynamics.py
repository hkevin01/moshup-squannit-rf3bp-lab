import numpy as np

from rf3bp_lab.dynamics.models import (
    FidelityWeights,
    cr3bp_rhs,
    rf3bp_breakdown,
    rf3bp_pulsating_rhs,
    rf3bp_pulsating_rhs_weighted,
    secondary_relative_kinematics,
)
from rf3bp_lab.dynamics.params import SystemParams


def test_rhs_shapes() -> None:
    p = SystemParams()
    state = np.array([0.4, 0.1, 0.0, 0.0, 0.3, 0.0], dtype=float)
    a = cr3bp_rhs(0.0, state, p)
    b = rf3bp_pulsating_rhs(0.0, state, p)
    assert a.shape == (6,)
    assert b.shape == (6,)
    assert np.all(np.isfinite(a))
    assert np.all(np.isfinite(b))


def test_secondary_kinematics_is_finite() -> None:
    p = SystemParams()
    kin = secondary_relative_kinematics(1.2, p)
    assert kin.r_vec.shape == (3,)
    assert kin.v_vec.shape == (3,)
    assert kin.a_vec.shape == (3,)
    assert kin.j_vec.shape == (3,)
    assert np.all(np.isfinite(kin.r_vec))
    assert np.all(np.isfinite(kin.v_vec))
    assert np.all(np.isfinite(kin.a_vec))
    assert np.all(np.isfinite(kin.j_vec))


def test_weighted_rhs_removes_selected_perturbations() -> None:
    p = SystemParams()
    state = np.array([0.4, 0.1, 0.0, 0.0, 0.3, 0.0], dtype=float)

    full = rf3bp_pulsating_rhs_weighted(t=0.4, state=state, p=p, fidelity=FidelityWeights())
    no_solar = rf3bp_pulsating_rhs_weighted(
        t=0.4,
        state=state,
        p=p,
        fidelity=FidelityWeights(pulsation=1.0, nonspherical=1.0, solar_gravity=0.0, srp=0.0),
    )

    assert np.linalg.norm(full[3:] - no_solar[3:]) > 0.0


def test_breakdown_components_are_finite() -> None:
    p = SystemParams()
    state = np.array([0.4, -0.1, 0.01, -0.05, 0.2, -0.01], dtype=float)
    b = rf3bp_breakdown(0.5, state, p)

    assert np.all(np.isfinite(b.grav_primary))
    assert np.all(np.isfinite(b.grav_secondary))
    assert np.all(np.isfinite(b.coriolis))
    assert np.all(np.isfinite(b.centrifugal))
    assert np.all(np.isfinite(b.pulsation))
    assert np.all(np.isfinite(b.solar_gravity))
    assert np.all(np.isfinite(b.srp))
