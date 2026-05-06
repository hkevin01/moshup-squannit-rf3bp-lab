import numpy as np

from rf3bp_lab.dynamics.params import SystemParams
from rf3bp_lab.shooting.hierarchical import HierarchicalShooter


def test_hierarchical_shooter_runs() -> None:
    p = SystemParams()
    shooter = HierarchicalShooter(p)
    seed = np.array([0.55, 0.0, 0.0, 0.0, 0.42, 0.02], dtype=float)
    out = shooter.solve(seed)
    assert "period" in out
    assert out["period"] > 0.0
    assert np.isfinite(out["residual_norm"])
    assert "stages" in out
    assert len(out["stages"]) > 2
    assert out["stages"][0].label == "cr3bp"


# ================================================================================
# Multi-Shooting and Orbit Family Tests
# ================================================================================

def test_multi_shooting_config_creation():
    """Test multi-shooting configuration creation."""
    from rf3bp_lab.shooting.multi_shooting import MultiShootingConfig
    
    config = MultiShootingConfig(num_arcs=3, max_nfev=5)
    
    assert config.num_arcs == 3
    assert config.max_nfev == 5
    assert config.solve_tol > 0


def test_arc_definition_creation():
    """Test arc definition creation."""
    from rf3bp_lab.shooting.multi_shooting import ArcDefinition
    
    arc = ArcDefinition(t_start=0.0, t_end=1.0, arc_index=0, is_periodic=True)
    
    assert arc.t_start == 0.0
    assert arc.t_end == 1.0
    assert arc.is_periodic is True


def test_orbit_family_point_creation():
    """Test orbit family point creation."""
    from rf3bp_lab.shooting.orbit_families import OrbitFamilyPoint
    
    state = np.array([0.5, 0.0, 0.0, 0.0, 0.4, 0.02])
    point = OrbitFamilyPoint(
        param_value=1.0,
        state=state,
        period=6.0,
        energy=-0.5,
        jacobi_constant=3.0,
        stability_index=1.1,
        objective_value=100.0,
    )
    
    assert point.param_value == 1.0
    assert point.period == 6.0
    assert np.allclose(point.state, state)


def test_family_config_creation():
    """Test family configuration creation."""
    from rf3bp_lab.shooting.orbit_families import FamilyConfig
    
    config = FamilyConfig(param_name="amplitude", num_points=15)
    
    assert config.param_name == "amplitude"
    assert config.num_points == 15
    assert config.param_min >= 0


def test_pareto_front_identification():
    """Test Pareto front identification for orbit families."""
    from rf3bp_lab.shooting.orbit_families import OrbitFamilyPoint, pareto_front
    
    state0 = np.array([0.5, 0.0, 0.0, 0.0, 0.4, 0.02])
    
    points = [
        OrbitFamilyPoint(1.0, state0, 5.0, -0.5, 3.0, 1.0),
        OrbitFamilyPoint(2.0, state0, 6.0, -0.4, 3.5, 1.1),
        OrbitFamilyPoint(3.0, state0, 7.0, -0.3, 4.0, 0.9),
    ]
    
    pareto_indices = pareto_front(points, ['energy', 'stability_index'])
    
    assert len(pareto_indices) > 0
    assert all(isinstance(i, (int, np.integer)) for i in pareto_indices)


def test_family_metrics_computation():
    """Test family metrics computation."""
    from rf3bp_lab.shooting.orbit_families import OrbitFamilyPoint, compute_family_metrics
    
    state0 = np.array([0.5, 0.0, 0.0, 0.0, 0.4, 0.02])
    
    points = [
        OrbitFamilyPoint(1.0, state0, 5.0, -0.5, 3.0, 1.0),
        OrbitFamilyPoint(2.0, state0, 6.0, -0.4, 3.5, 1.1),
        OrbitFamilyPoint(3.0, state0, 7.0, -0.3, 4.0, 0.9),
    ]
    
    metrics = compute_family_metrics(points)
    
    assert 'period_mean' in metrics
    assert 'energy_range' in metrics
    assert metrics['period_mean'] > 0
