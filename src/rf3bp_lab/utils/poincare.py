from __future__ import annotations

"""ID: RF3BP-LAB-UTILS-POINCARE
Requirement: Provide Poincare section and event-surface diagnostics for RF3BP orbits.
Purpose: Reveal the phase-space structure (islands, chaotic zones, KAM tori) of
         trajectories in the RF3BP rotating frame.
Rationale: Poincare sections are the standard tool for understanding periodic and
           quasi-periodic orbit families in conservative Hamiltonian systems.
           scipy.integrate.solve_ivp event detection provides exact crossing times.
Inputs: Initial state, propagation span, section plane definition, dynamics RHS.
Outputs: Arrays of (x, v) crossing coordinates, crossing times, stability metrics.
Preconditions: Initial state is 6D; dynamics RHS returns a 6D derivative.
Postconditions: Crossings are ordered by time; each crossing is a finite 6-vector.
Assumptions: Near-conservative dynamics; section crossings are well separated.
Side Effects: None.
Failure Modes: Trajectories that never cross the section return empty arrays.
Error Handling: Integration timeouts return partial results; degenerate sections are guarded.
Constraints: Accurate event location requires dense_output=True in solve_ivp.
Verification: Unit tests in tests/test_dynamics.py.
References: Koon et al., Nonlinear Dynamics in Aerospace Applications (2008).
            Szebehely, Theory of Orbits (1967), Chapter 4.
"""

import numpy as np
from dataclasses import dataclass
from typing import Callable
from scipy.integrate import solve_ivp


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class SectionPlane:
    """ID: RF3BP-LAB-POINCARE-PLANE
    Requirement: Define a Poincare section as a hyperplane n·(x - x0) = 0.
    Purpose: Specify where trajectories are sampled.
    Inputs: Normal vector n (3D), offset x0 (3D), crossing direction (+1 or -1).
    Outputs: Immutable section descriptor.
    Preconditions: n is non-zero.
    Postconditions: normal is stored as a unit vector.
    """
    normal: np.ndarray      # unit normal to the section plane (3D)
    offset: np.ndarray      # a point on the plane (3D)
    direction: int = 1      # +1: positive crossings only; -1: negative; 0: both

    def __post_init__(self) -> None:
        # Store as unit vector to make crossing check numerically stable
        n_norm = np.linalg.norm(self.normal)
        if n_norm < 1e-14:
            raise ValueError("SectionPlane normal must be non-zero.")
        object.__setattr__(self, "normal", self.normal / n_norm)

    def signed_distance(self, pos: np.ndarray) -> float:
        """Signed distance of pos from the section plane."""
        return float(np.dot(self.normal, pos - self.offset))


@dataclass(frozen=True)
class PoincareCrossing:
    """ID: RF3BP-LAB-POINCARE-CROSSING
    Requirement: Store one Poincare section crossing with full state and diagnostics.
    Purpose: Collect phase-space coordinates at each section crossing.
    Inputs: Time, full 6-state, in-plane (u,v) coordinates, crossing index.
    Outputs: Immutable crossing record.
    """
    t: float
    state: np.ndarray          # full 6-vector at crossing
    u: float                   # first in-plane position coordinate
    v: float                   # second in-plane position coordinate
    crossing_index: int


@dataclass(frozen=True)
class PoincareSection:
    """ID: RF3BP-LAB-POINCARE-SECTION
    Requirement: Hold the full result of a Poincare section computation.
    Purpose: Package crossings, velocity projections, and metadata.
    Inputs: List of crossings, section plane, propagation parameters.
    Outputs: Immutable section result container.
    """
    crossings: list[PoincareCrossing]
    plane: SectionPlane
    t_span: tuple[float, float]
    n_crossings: int

    @property
    def u_coords(self) -> np.ndarray:
        return np.array([c.u for c in self.crossings])

    @property
    def v_coords(self) -> np.ndarray:
        return np.array([c.v for c in self.crossings])

    @property
    def times(self) -> np.ndarray:
        return np.array([c.t for c in self.crossings])

    @property
    def states(self) -> np.ndarray:
        if not self.crossings:
            return np.empty((0, 6))
        return np.stack([c.state for c in self.crossings])


# ---------------------------------------------------------------------------
# Section plane constructors
# ---------------------------------------------------------------------------

def xz_plane_section(
    y_offset: float = 0.0,
    direction: int = 1,
) -> SectionPlane:
    """ID: RF3BP-LAB-POINCARE-XZ
    Requirement: Create a y = y_offset section plane (standard CR3BP section).
    Inputs: y_offset scalar, crossing direction.
    Outputs: SectionPlane with normal (0,1,0).
    """
    return SectionPlane(
        normal=np.array([0.0, 1.0, 0.0]),
        offset=np.array([0.0, y_offset, 0.0]),
        direction=direction,
    )


def xy_plane_section(
    z_offset: float = 0.0,
    direction: int = 1,
) -> SectionPlane:
    """ID: RF3BP-LAB-POINCARE-XY
    Requirement: Create a z = z_offset section plane.
    Inputs: z_offset scalar, crossing direction.
    Outputs: SectionPlane with normal (0,0,1).
    """
    return SectionPlane(
        normal=np.array([0.0, 0.0, 1.0]),
        offset=np.array([0.0, 0.0, z_offset]),
        direction=direction,
    )


# ---------------------------------------------------------------------------
# In-plane coordinate extraction
# ---------------------------------------------------------------------------

def _in_plane_coords(
    state: np.ndarray,
    plane: SectionPlane,
) -> tuple[float, float]:
    """ID: RF3BP-LAB-POINCARE-INPLANE
    Requirement: Project position onto the section plane's two in-plane axes.
    Inputs: 6D state, SectionPlane.
    Outputs: (u, v) tuple of in-plane coordinates.
    Rationale: Provides the canonical 2D Poincare map coordinates.
    """
    pos = state[:3]
    # Build two orthogonal axes in the section plane
    n = plane.normal
    # Choose a reference vector not parallel to n
    ref = np.array([1.0, 0.0, 0.0])
    if abs(np.dot(n, ref)) > 0.9:
        ref = np.array([0.0, 1.0, 0.0])
    e1 = ref - np.dot(ref, n) * n
    e1 /= np.linalg.norm(e1) + 1e-14
    e2 = np.cross(n, e1)
    e2 /= np.linalg.norm(e2) + 1e-14

    r_proj = pos - plane.offset
    return float(np.dot(r_proj, e1)), float(np.dot(r_proj, e2))


# ---------------------------------------------------------------------------
# Main section computation
# ---------------------------------------------------------------------------

def compute_poincare_section(
    state0: np.ndarray,
    t_span: tuple[float, float],
    rhs: Callable,
    plane: SectionPlane,
    *,
    max_crossings: int = 500,
    rtol: float = 1e-8,
    atol: float = 1e-10,
    max_step: float = 0.2,
) -> PoincareSection:
    """ID: RF3BP-LAB-POINCARE-COMPUTE
    Requirement: Integrate a trajectory and record all section crossings.
    Purpose: Produce a Poincare section map for phase-space structure discovery.
    Inputs: state0 (6,), t_span, dynamics rhs f(t,y), SectionPlane, solver options.
    Outputs: PoincareSection with all crossing records.
    Preconditions: rhs accepts (t, state) and returns (6,); state0 is finite.
    Postconditions: Crossings are time-ordered; states at crossings are finite.
    Failure Modes: No crossings if trajectory runs parallel to section.
    Error Handling: Returns empty section if integration fails.
    Rationale: solve_ivp dense_output + event detection gives sub-step accurate crossing times.
    References: scipy.integrate.solve_ivp documentation; Koon et al. 2008.
    """
    state0 = np.asarray(state0, dtype=float)
    crossings: list[PoincareCrossing] = []
    crossing_count = [0]  # mutable container for event closure

    def section_event(t: float, y: np.ndarray) -> float:
        return plane.signed_distance(y[:3])

    section_event.terminal = False   # type: ignore[attr-defined]
    section_event.direction = float(plane.direction)  # type: ignore[attr-defined]

    try:
        sol = solve_ivp(
            rhs,
            t_span,
            state0,
            method="DOP853",
            events=[section_event],
            dense_output=True,
            rtol=rtol,
            atol=atol,
            max_step=max_step,
        )
    except Exception:
        return PoincareSection(
            crossings=[], plane=plane, t_span=t_span, n_crossings=0
        )

    if sol.t_events is None or len(sol.t_events[0]) == 0:
        return PoincareSection(
            crossings=[], plane=plane, t_span=t_span, n_crossings=0
        )

    for ti in sol.t_events[0]:
        if crossing_count[0] >= max_crossings:
            break
        state_at_cross = sol.sol(ti)  # dense interpolation
        u, v = _in_plane_coords(state_at_cross, plane)
        crossings.append(PoincareCrossing(
            t=float(ti),
            state=np.array(state_at_cross, dtype=float),
            u=u,
            v=v,
            crossing_index=crossing_count[0],
        ))
        crossing_count[0] += 1

    return PoincareSection(
        crossings=crossings,
        plane=plane,
        t_span=t_span,
        n_crossings=len(crossings),
    )


# ---------------------------------------------------------------------------
# Stability index from crossing return-time spread
# ---------------------------------------------------------------------------

def poincare_stability_index(section: PoincareSection) -> float:
    """ID: RF3BP-LAB-POINCARE-STABILITY
    Requirement: Estimate stability from the spread in successive crossing times.
    Purpose: Distinguish periodic (constant spacing) from chaotic (variable spacing) orbits.
    Inputs: PoincareSection with at least 3 crossings.
    Outputs: Coefficient of variation of return times (0 = perfectly periodic).
    Preconditions: section.n_crossings >= 3.
    Postconditions: Output is non-negative scalar.
    Rationale: Periodic orbits have constant return times; chaotic ones do not.
    """
    if section.n_crossings < 3:
        return float("nan")
    times = section.times
    dt = np.diff(times)
    mean_dt = np.mean(dt)
    if mean_dt < 1e-14:
        return float("nan")
    return float(np.std(dt) / mean_dt)


# ---------------------------------------------------------------------------
# Plotting helper
# ---------------------------------------------------------------------------

def plot_poincare_section(
    section: PoincareSection,
    title: str = "Poincare Section",
    color: str = "#264653",
) -> None:
    """ID: RF3BP-LAB-POINCARE-PLOT
    Requirement: Produce a 2D Poincare section scatter plot.
    Purpose: Visualise phase-space islands, tori, and chaotic zones.
    Inputs: PoincareSection result, plot title string, point colour.
    Outputs: matplotlib figure (caller saves or shows).
    Preconditions: section has at least one crossing.
    Postconditions: Figure is created with labelled axes.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7.0, 6.5))

    if section.n_crossings == 0:
        ax.text(0.5, 0.5, "No section crossings recorded.",
                ha="center", va="center", transform=ax.transAxes)
    else:
        u = section.u_coords
        v = section.v_coords
        sc = ax.scatter(u, v, s=4, c=section.times, cmap="viridis", alpha=0.8)
        plt.colorbar(sc, ax=ax, label="time")
        ax.set_xlabel("u  (in-plane position)")
        ax.set_ylabel("v  (in-plane position)")

        stab = poincare_stability_index(section)
        if np.isfinite(stab):
            ax.set_title(f"{title}  |  {section.n_crossings} crossings  |  "
                         f"stability index = {stab:.3f}", fontweight="bold")
        else:
            ax.set_title(f"{title}  |  {section.n_crossings} crossings", fontweight="bold")

    ax.grid(alpha=0.2)
    plt.tight_layout()
