import numpy as np

from rf3bp_lab.dynamics.models import (
    FidelityWeights,
    compare_cr3bp_rf3bp,
    cr3bp_effective_potential,
    cr3bp_jacobi_constant,
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


def test_cr3bp_scalar_diagnostics_are_finite() -> None:
    p = SystemParams()
    state = np.array([0.3, 0.2, 0.01, -0.05, 0.12, -0.01], dtype=float)

    omega = cr3bp_effective_potential(state, p)
    jacobi = cr3bp_jacobi_constant(state, p)

    assert np.isfinite(omega)
    assert np.isfinite(jacobi)


def test_cr3bp_rf3bp_gap_metrics_are_finite() -> None:
    p = SystemParams()
    state = np.array([0.35, 0.11, -0.02, 0.03, 0.18, 0.01], dtype=float)

    gap = compare_cr3bp_rf3bp(0.25, state, p)

    assert gap.cr3bp_acc.shape == (3,)
    assert gap.rf3bp_acc.shape == (3,)
    assert gap.delta_acc.shape == (3,)
    assert np.all(np.isfinite(gap.cr3bp_acc))
    assert np.all(np.isfinite(gap.rf3bp_acc))
    assert np.isfinite(gap.delta_norm)
    assert np.isfinite(gap.relative_gap)
    assert gap.delta_norm >= 0.0
    assert gap.relative_gap >= 0.0


def test_cr3bp_rf3bp_gap_changes_with_fidelity() -> None:
    p = SystemParams()
    state = np.array([0.35, 0.11, -0.02, 0.03, 0.18, 0.01], dtype=float)

    gap_low = compare_cr3bp_rf3bp(
        0.25,
        state,
        p,
        fidelity=FidelityWeights(pulsation=0.0, nonspherical=0.0, solar_gravity=0.0, srp=0.0),
    )
    gap_high = compare_cr3bp_rf3bp(0.25, state, p, fidelity=FidelityWeights())

    assert gap_high.delta_norm >= gap_low.delta_norm
