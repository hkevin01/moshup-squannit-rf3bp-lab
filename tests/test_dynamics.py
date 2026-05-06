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


# ================================================================================
# Advanced Features Tests
# ================================================================================

def test_harmonic_gravity_model_creation():
    """Test creation of spherical harmonic gravity model."""
    from rf3bp_lab.dynamics.advanced_gravity import create_default_harmonic_model
    
    model = create_default_harmonic_model(mu=1.0, radius=0.5, max_degree=4)
    
    assert model.mu == 1.0
    assert model.radius == 0.5
    assert len(model.coefficients) > 0
    assert model.use_normalized is True


def test_eclipse_detection_no_shadow():
    """Test eclipse detection when spacecraft is fully illuminated."""
    from rf3bp_lab.dynamics.advanced_gravity import detect_eclipse
    
    sc_pos = np.array([1.0, 0.0, 0.0])
    primary_pos = np.array([0.0, 0.0, 0.0])
    secondary_pos = np.array([-0.01, 0.0, 0.0])
    sun_direction = np.array([1.0, 0.0, 0.0]) / np.sqrt(1.0 + 1e-12)
    
    eclipse = detect_eclipse(
        sc_pos, primary_pos, 0.1, secondary_pos, 0.05, sun_direction
    )
    
    assert eclipse.shadow_fraction >= 0.0
    assert eclipse.shadow_fraction <= 1.0
    assert isinstance(eclipse.in_primary_shadow, bool)


def test_eclipse_aware_srp():
    """Test eclipse-aware SRP attenuation."""
    from rf3bp_lab.dynamics.advanced_gravity import detect_eclipse, eclipse_aware_srp
    
    srp_base = np.array([1e-6, 0.0, 0.0])
    
    sc_pos = np.array([1.0, 0.0, 0.0])
    primary_pos = np.array([0.0, 0.0, 0.0])
    secondary_pos = np.array([-0.01, 0.0, 0.0])
    sun_direction = np.array([1.0, 0.0, 0.0]) / np.sqrt(1.0 + 1e-12)
    
    eclipse = detect_eclipse(
        sc_pos, primary_pos, 0.1, secondary_pos, 0.05, sun_direction
    )
    
    srp_eclipse = eclipse_aware_srp(srp_base, eclipse)
    
    assert np.linalg.norm(srp_eclipse) <= np.linalg.norm(srp_base) * (1 + 1e-12)


def test_variational_state_initialization():
    """Test variational state initialization."""
    from rf3bp_lab.dynamics.variational import _flatten_variational, _unflatten_variational
    
    state = np.array([0.5, 0.0, 0.0, 0.0, 0.4, 0.02])
    stm = np.eye(6)
    
    flat = _flatten_variational(state, stm)
    assert flat.shape == (42,)
    
    state_rec, stm_rec = _unflatten_variational(flat)
    assert np.allclose(state, state_rec)
    assert np.allclose(stm, stm_rec)


def test_covariance_propagation_preserves_psd():
    """Test that covariance propagation preserves positive-definiteness."""
    from rf3bp_lab.dynamics.variational import covariance_at_time
    
    cov_0 = np.diag([1.0, 1.0, 1.0, 0.01, 0.01, 0.01])
    stm = np.eye(6)
    
    cov_t = covariance_at_time(cov_0, stm)
    
    evals = np.linalg.eigvalsh(cov_t)
    assert np.all(evals > -1e-10)


def test_dilution_of_precision_computes_finite():
    """Test DOP metric computation."""
    from rf3bp_lab.dynamics.variational import dilution_of_precision
    
    cov = np.diag([1.0, 1.0, 1.0, 0.1, 0.1, 0.1])
    
    dop = dilution_of_precision(cov)
    
    assert 'position_std_dev' in dop
    assert 'velocity_std_dev' in dop
    assert dop['position_std_dev'] > 0
    assert dop['velocity_std_dev'] > 0
