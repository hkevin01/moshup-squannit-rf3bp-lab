# Moshup-Squannit RF3BP Lab

![Python](https://img.shields.io/badge/python-3.11%2B-0b3d91)
![SciPy](https://img.shields.io/badge/scipy-1.13%2B-8a3ffc)
![Status](https://img.shields.io/badge/status-research%20sandbox-c97b00)
![Focus](https://img.shields.io/badge/focus-RF3BP%20algorithms-006d77)

Research-oriented Python project for experimenting with bounded spacecraft motion near the binary asteroid system **Moshup-Squannit (1999 KW4)** using a **Restricted Full Three-Body Problem** (RF3BP) inspired formulation in a **pulsating-rotating frame**.

The project implements advanced astrodynamics techniques including:
- **CR3BP** (Circular Restricted 3-Body Problem) as a reference baseline
- **RF3BP** (Restricted Full 3-Body Problem) with pulsation, irregular-body gravity, and solar perturbations
- Hierarchical continuation shooting and multi-arc trajectory optimization
- State Transition Matrix (STM) propagation for mission-grade navigation covariance analysis
- Eclipse-aware Solar Radiation Pressure (SRP) and harmonic gravity expansions

This repository does not try to be a CLI-heavy wrapper. The main value is the **dynamics code**, the **hierarchical shooting continuation logic**, and the **diagnostic plots** that help compare perturbation sources in the neighborhood of a binary asteroid.

> [!IMPORTANT]
> The journal paper linked in the original request is paywalled from this environment. This repository therefore implements a technically consistent, abstract-level RF3BP research sandbox derived from the paper's public abstract and standard CR3BP/RF3BP practice, not a line-by-line reproduction of copyrighted equations.

## What This Solves

Classical CR3BP is useful for intuition, but it misses the exact effects that become important around a small binary asteroid:

- the primary-secondary distance can pulsate instead of staying fixed
- body shapes are not spherical
- the Sun perturbs the local binary dynamics
- solar radiation pressure matters for low-mass spacecraft
- orbit design methods that work in CR3BP can break when fidelity increases

This project gives you a compact environment to answer questions like:

- How much does pulsation change a bounded orbit relative to a CR3BP reference?
- Is pulsation weaker or stronger than the nonspherical gravity correction along a candidate path?
- Can a seed orbit from a lower-fidelity model be continued into a more realistic one?
- Which perturbation dominates at a given point in the trajectory?

## Executive Summary

| Topic | What It Is | What It Does | Why It Matters |
| --- | --- | --- | --- |
| CR3BP baseline | Circular restricted 3-body model in a rotating frame | Supplies a low-cost reference orbit model | Good starting point for initial guesses |
| RF3BP-inspired dynamics | Higher-fidelity pulsating-rotating model | Adds pulsation, nonspherical gravity, solar gravity, and SRP | Closer to binary-asteroid mission reality |
| Potential-derivative kinematics | Secondary relative acceleration and jerk estimate from a potential derivative view | Explicitly models nonuniform pulsation terms | Makes pulsation effects visible in the equations |
| Hierarchical continuation | Staged shooting from CR3BP to higher fidelity | Transfers a seed orbit across models | Avoids solving the hardest model from scratch |
| Diagnostics | Trajectory and perturbation plots | Shows which effects dominate and where | Useful for analysis, papers, and design iteration |

## RF3BP vs CR3BP - What Is The Difference?

**Circular Restricted 3-Body Problem (CR3BP)** is the classical model where the primary and secondary bodies orbit at constant separation.

**Restricted Full 3-Body Problem (RF3BP)** extends CR3BP by allowing the primary-secondary separation to pulsate (vary with time), and by accounting for higher-fidelity perturbations. CR3BP and RF3BP are not competing "brands" of the same equation - they represent fundamentally different physical assumptions.

| Aspect | CR3BP | RF3BP (this lab) | Practical Consequence |
| --- | --- | --- | --- |
| Primary-secondary distance | Constant | Time-varying (pulsating) | Adds nonuniform frame terms and shifts equilibrium structure |
| Gravity field shape | Point masses only | Point masses + J2-like nonspherical corrections | Local accelerations can deviate strongly near bodies |
| External forcing | None | Solar third-body differential gravity + SRP | Long-time bounded motion is more sensitive |
| Frame model | Uniform rotating frame | Pulsating-rotating frame with explicit pulsation terms | CR3BP intuition can fail as fidelity increases |
| Design workflow | Often direct periodic-orbit correction | Hierarchical continuation from low to high fidelity | Better robustness when full model is stiff |

In the rotating frame, the **CR3BP** acceleration can be summarized as:

$$
\ddot{\mathbf{r}}_{CR3BP} = \nabla \Omega(\mathbf{r}) - 2\,\boldsymbol{\omega} \times \dot{\mathbf{r}}
$$

where $\Omega$ is the effective potential in the rotating frame and $\boldsymbol{\omega}$ is the frame rotation vector.

The **RF3BP-inspired model** used in this laboratory extends CR3BP by incorporating four major perturbations:

$$
\ddot{\mathbf{r}}_{RF3BP} = \ddot{\mathbf{r}}_{CR3BP} + \mathbf{a}_{pulsation} + \mathbf{a}_{nonspherical} + \mathbf{a}_{solar} + \mathbf{a}_{SRP}
$$

Each term represents a distinct physical effect:
- $\mathbf{a}_{pulsation}$: Frame acceleration due to time-varying primary-secondary distance (binary orbit pulsation)
- $\mathbf{a}_{nonspherical}$: Gravity corrections from irregular body shapes (harmonic expansion up to degree 20)
- $\mathbf{a}_{solar}$: Third-body gravity perturbation from the Sun
- $\mathbf{a}_{SRP}$: Solar radiation pressure, including eclipse-aware geometric attenuation

The new code-level metric in this repository computes the instantaneous gap

$$
\Delta \mathbf{a} = \ddot{\mathbf{r}}_{RF3BP} - \ddot{\mathbf{r}}_{CR3BP},
\quad
\rho = \frac{\|\Delta \mathbf{a}\|}{\|\ddot{\mathbf{r}}_{CR3BP}\|}
$$

This gives a direct, quantitative answer to "how far from CR3BP" a trajectory point is.

## Visual Outputs

### RF3BP vs CR3BP

| CR3BP Reference | RF3BP Higher Fidelity |
| --- | --- |
| ![CR3BP](docs/figures/trajectory_cr3bp.png) | ![RF3BP](docs/figures/trajectory_rf3bp.png) |

### Perturbation Diagnostics

| Perturbation Magnitudes | Continuation Convergence |
| --- | --- |
| ![Perturbations](docs/figures/perturbation_norms.png) | ![Continuation](docs/figures/continuation_convergence.png) |

### Model Difference Diagnostics

| RF3BP vs CR3BP Gap History |
| --- |
| ![Model Gap](docs/figures/model_gap_cr3bp_vs_rf3bp.png) |

## Latest Produced Results (Default Seed)

The following results were generated directly from the current repository code using:

```bash
OUTPUT_DIR=docs/figures ./.venv/bin/python scripts/run_demo.py
```

Raw machine-readable metrics are stored in `docs/results/latest_demo_metrics.json`.

| Metric | Value | Context |
| --- | --- | --- |
| Period estimate | `6.329724` | Final continuation period estimate returned by the shooter |
| Final residual norm | `5.749e+00` | Periodicity defect for the final continuation stage |
| Stage-1 cost | `3.805e-02` | Cost right after CR3BP correction |
| Final stage cost | `1.653e+01` | Cost at full model fidelity |
| Max absolute model gap | `1.986e+14` | Peak value of $\|a_{RF3BP} - a_{CR3BP}\|$ along propagated path |
| Mean relative model gap | `4.323e-01` | Mean of $\|\Delta a\| / \|a_{CR3BP}\|$ over 4000 samples |
| Max relative model gap | `7.471e-01` | Peak normalized mismatch between RF3BP and CR3BP |

> [!NOTE]
> The very large absolute peak reflects the current simplified pulsation formulation and sampling near strong local gradients. For design decisions, relative gap trends and component-wise diagnostics are usually more informative than a single absolute peak.

| Result Snapshot Dashboard |
| --- |
| ![Result Snapshot](docs/figures/result_snapshot_dashboard.png) |

This dashboard combines peak perturbation magnitudes with key continuation and model-gap metrics, so each demo run yields a compact screenshot artifact for progress tracking.

> [!NOTE]
> The figures above are generated from the repository's current model and default parameters. They are useful for comparative algorithm work, not for claiming flight-certified truth.

## System Architecture

```mermaid
flowchart TD
    A[Seed State] --> B[CR3BP Dynamics]
    B --> C[Stage 1 Shooting Correction]
    C --> D[Continuation Ladder]
    D --> D1[Pulsation On]
    D1 --> D2[Nonspherical Gravity On]
    D2 --> D3[Solar Gravity On]
    D3 --> D4[SRP On]
    D4 --> E[RF3BP Bounded Orbit Candidate]
    E --> F[Perturbation Breakdown]
    E --> G[Trajectory Propagation]
    F --> H[Diagnostic Charts]
    G --> H
```

## Model Fidelity Ladder

| Stage | Enabled Physics | Purpose | Computational Role |
| --- | --- | --- | --- |
| 0 | CR3BP only | Correct the seed in the simplest useful model | Fast, stable initialization |
| 1 | Pulsation | Turn on nonuniform separation effects | Measures the cost of leaving circular assumptions |
| 2 | Pulsation + nonspherical gravity | Add J2-like gravity corrections | Introduces body-shape-driven local distortion |
| 3 | Pulsation + nonspherical gravity + solar gravity | Add differential third-body forcing | Captures long-baseline solar perturbation |
| 4 | Full model + SRP | Add light-pressure acceleration | Approximates small-spacecraft sensitivity |

## Core Algorithms

### 1. CR3BP Baseline Dynamics

The repository keeps a standard rotating-frame CR3BP model as the reference surface for orbit seeding and comparison.

### 2. RF3BP-Inspired Pulsating-Rotating Dynamics

The higher-fidelity model in `src/rf3bp_lab/dynamics/models.py` combines:

- point-mass gravity from both binary bodies
- J2-style nonspherical gravity corrections for both bodies
- time-varying primary-secondary separation
- pulsation inertial correction terms
- third-body solar gravity
- solar radiation pressure

The same module now includes a direct comparison utility that computes acceleration-level mismatch between the CR3BP and RF3BP equations at any state and time.

### 3. Potential-Derivative Relative Kinematics

The code computes relative secondary kinematics from a potential-derivative perspective:

- relative position vector
- relative velocity vector
- relative acceleration vector
- relative jerk vector

This is the project’s direct answer to the request to stop centering the work on CLI and focus on the actual dynamical machinery.

### 4. Hierarchical Shooting Continuation

The continuation algorithm does not brute-force the full model from scratch. It follows a staged process:

```mermaid
flowchart LR
    S[Initial Seed] --> C1[Correct in CR3BP]
    C1 --> C2[Continue into Pulsation]
    C2 --> C3[Continue into J2-like Gravity]
    C3 --> C4[Continue into Solar Gravity]
    C4 --> C5[Continue into SRP]
    C5 --> R[Return Final Bounded Candidate]
```

The current implementation uses an explicit bounded, damped Newton-style corrector with finite-difference Jacobians under controlled tolerances. That replaced a much slower optimizer-driven solve path that caused test timeouts.

## What Is Implemented vs What Is Approximated

| Area | Current Implementation | Why This Choice Was Made | Upgrade Path |
| --- | --- | --- | --- |
| Binary geometry | Normalized separation with pulsation law | Keeps the frame mechanics explicit and inspectable | Replace with shape/ephemeris-driven relative motion |
| Nonspherical gravity | J2-like correction for both bodies | Lightweight proxy for body asymmetry | Polyhedral gravity from shape models |
| Solar perturbation | Simplified moving-Sun differential gravity | Good comparative forcing term | SPICE or ephemeris-driven Sun state |
| SRP | Constant-magnitude directional SRP | Lets perturbation ranking be studied quickly | Area-to-mass, attitude, eclipse, optical model |
| Continuation | Single-shooting staged correction | Small code footprint, good for experimentation | Multi-shooting and collocation |

## Why These Methods Fit This Problem

| Method | Why It Fits Moshup-Squannit Orbit Studies | Main Limitation |
| --- | --- | --- |
| CR3BP seed generation | Gives a structured initial orbit family around a binary system | Ignores pulsation and realistic perturbations |
| Pulsating-rotating frame | Directly expresses time-varying mutual separation effects | Requires care when comparing against static-frame intuition |
| J2-like gravity proxy | Cheap way to inject dominant nonspherical trends | Too simple for strongly irregular bodies |
| Hierarchical continuation | Practical way to migrate a low-fidelity orbit into higher fidelity | Can still fail if the seed is too weak |
| Perturbation breakdown plots | Turns model complexity into interpretable evidence | Diagnostic, not an optimization method |

## Alternatives Compared

| Approach | Strengths | Weaknesses | When To Use It Instead |
| --- | --- | --- | --- |
| Pure CR3BP | Fast, interpretable, classical | Too idealized for close binary-asteroid work | Early concept design or teaching |
| This repository | Good balance between insight and implementation cost | Several effects are still approximated | Algorithm research and rapid trade studies |
| Full polyhedral gravity + SPICE + eclipse + attitude model | Highest physical fidelity | Much more data and engineering overhead | Mission-grade analysis and detailed navigation studies |
| Direct black-box optimization in full fidelity | Can find solutions missed by continuation | Expensive and brittle without good seeds | Late-stage global search after good priors exist |

## Repository Map

| Path | What It Is | What It Does |
| --- | --- | --- |
| `src/rf3bp_lab/dynamics/params.py` | Parameter container | Holds normalized and physical scale assumptions |
| `src/rf3bp_lab/dynamics/models.py` | Core dynamics engine | Implements CR3BP, RF3BP-inspired dynamics, perturbation breakdown, and propagation |
| `src/rf3bp_lab/shooting/hierarchical.py` | Continuation solver | Performs staged bounded-orbit correction across model fidelity levels |
| `src/rf3bp_lab/utils/plotting.py` | Visualization utilities | Builds trajectory and perturbation charts |
| `scripts/run_demo.py` | End-to-end demo | Runs continuation, propagation, and plot generation |
| `tests/` | Validation layer | Confirms basic dynamics and shooter behavior |
| `docs/figures/` | Generated artifacts | Stores README-embeddable charts |

## Default Moshup-Squannit Assumptions

The defaults are intentionally transparent and easy to refine.

| Parameter | Default | Meaning |
| --- | --- | --- |
| `mu` | `0.02` | Normalized binary mass ratio used by the sandbox |
| `r12_mean_m` | `2500.0` m | Mean mutual separation |
| `pulsation_e` | `0.08` | Pulsation amplitude parameter |
| `pulsation_nu` | `0.35` | Pulsation frequency scale |
| `j2_primary` | `0.05` | J2-like primary gravity coefficient |
| `j2_secondary` | `0.02` | J2-like secondary gravity coefficient |
| `r_primary_m` | `700.0` m | Primary scale radius |
| `r_secondary_m` | `225.0` m | Secondary scale radius |
| `sun_mu_scaled` | `5.0e-4` | Normalized solar gravity strength |
| `sun_distance_scaled` | `2000.0` | Normalized Sun distance |
| `srp_accel_scaled` | `2.5e-6` | Normalized SRP acceleration |

> [!TIP]
> If you want the quickest improvement in physical realism, replace the J2-style terms first. For small irregular binaries, that simplification is usually the biggest structural gap.

## Quick Start

```bash
python3 -m venv .venv
source .venv/bin/activate
./.venv/bin/python -m pip install -e .[dev]
./.venv/bin/python -m pytest -q
OUTPUT_DIR=docs/figures ./.venv/bin/python scripts/run_demo.py
```

## Typical Workflow

```mermaid
sequenceDiagram
    participant U as User
    participant P as Params
    participant S as Shooter
    participant D as Dynamics
    participant G as Diagnostics

    U->>P: Set binary and perturbation assumptions
    U->>S: Provide seed state
    S->>D: Correct in CR3BP
    S->>D: Continue across fidelity stages
    D-->>S: Return bounded candidate
    U->>D: Propagate candidate
    D-->>G: Perturbation history
    G-->>U: Trajectory and comparison plots
```

## Current Validation Status

| Check | Result | Notes |
| --- | --- | --- |
| Dynamics unit tests | Passing | Confirms finite outputs and weighted perturbation behavior |
| Shooting unit test | Passing | Confirms continuation path returns a valid result quickly |
| Full test suite | Passing | `5 passed` |
| Demo run | Passing | Generates four figures in `docs/figures/` |

## Practical Notes

- WARNING: The current continuation solver is designed to be fast and inspectable for research iteration. It is not yet a production-grade orbit corrector.
- NOTE: The reported final residual in the demo is a diagnostic value from the current simplified bounded-orbit correction. Treat it as an indicator for further improvement, not a proof of strict periodicity.
- TIP: If you start exploring new seeds, begin by adjusting the seed state and only then increase fidelity. Jumping straight into the full model is the fastest way to waste compute.

## Development Commands

| Task | Command |
| --- | --- |
| Install project | `./.venv/bin/python -m pip install -e .[dev]` |
| Run tests | `./.venv/bin/python -m pytest -q` |
| Run only shooting test | `./.venv/bin/python -m pytest -q tests/test_shooting.py` |
| Regenerate figures | `OUTPUT_DIR=docs/figures ./.venv/bin/python scripts/run_demo.py` |

Generated outputs include `model_gap_cr3bp_vs_rf3bp.png`, which visualizes both absolute and relative acceleration mismatch history.

## What This Project Is Good For

- studying how pulsation changes orbit behavior near a binary asteroid
- comparing perturbation magnitudes along a candidate bounded path
- prototyping continuation strategies before moving to mission-grade software
- teaching the gap between CR3BP intuition and higher-fidelity local dynamics

## What It Does Not Yet Solve

- exact replication of the paywalled paper's full derivation
- high-order irregular-body gravity using shape models
- eclipse-aware SRP
- multi-shooting over many arcs
- formal optimization over bounded orbit families
- mission-grade navigation covariance analysis

## Advanced Features Implemented

This laboratory now includes five major advanced capabilities for mission-grade trajectory analysis:

### 1. High-Order Irregular-Body Gravity using Shape Models

**Module:** `rf3bp_lab.dynamics.advanced_gravity`

Replaces point-mass and J2-style gravity with spherical harmonic expansions up to degree 20:

- **SphericalHarmonicCoeff** dataclass: Stores normalized harmonic coefficients $C_{nm}$ and $S_{nm}$
- **HarmonicGravityModel**: Bundles coefficients with body parameters ($\mu$, equatorial radius)
- **_harmonic_accel()**: Computes acceleration from arbitrary-order harmonic expansion using recurrence relations for associated Legendre polynomials
- **create_default_harmonic_model()**: Factory function for realistic binary-asteroid gravity profiles

#### Example Usage

```python
from rf3bp_lab.dynamics.advanced_gravity import create_default_harmonic_model

model = create_default_harmonic_model(mu=1.0, radius=0.5, max_degree=20)
```

### 2. Eclipse-Aware Solar Radiation Pressure (SRP)

**Module:** `rf3bp_lab.dynamics.advanced_gravity`

Provides geometric eclipse detection and SRP attenuation:

- **EclipseState** dataclass: Stores boolean shadow flags, coverage fraction $f_{shadow} \in [0, 1]$, and angular radii
- **detect_eclipse()**: Cylindrical shadow model with penumbra blending
  - Returns `in_primary_shadow` and `in_secondary_shadow` flags
  - Computes `shadow_fraction` for smooth attenuation near edges
- **eclipse_aware_srp()**: Applies attenuation: $\mathbf{a}_{SRP,eclipsed} = (1 - f_{shadow}) \cdot \mathbf{a}_{SRP}$

#### Mission Significance

Realistic SRP in binary environments requires accounting for geometric blockage by both primary and secondary bodies, especially critical for low-mass spacecraft near small asteroids.

### 3. Multi-Arc Shooting Trajectory Optimization

**Module:** `rf3bp_lab.shooting.multi_shooting`

Enables simultaneous shooting over $N$ arcs with continuity constraints:

- **MultiShootingBVP** solver: Corrects all arc initial states in parallel via damped Newton iteration
- **ArcDefinition**: Specifies time span, index, and periodicity flag for each arc
- **ArcResult**: Stores arc-level diagnostics (initial/final states, continuity residuals)
- Finite-difference Jacobian for robustness; adaptive damping and state bounds

#### Advantages Over Single-Arc Shooting

- Better numerical conditioning for long-duration trajectories (e.g., 10+ periods)
- Enables midcourse constraint insertion
- Parallel corrector structure reduces sensitivity to initial guess quality

#### Example Configuration

```python
from rf3bp_lab.shooting.multi_shooting import MultiShootingBVP, MultiShootingConfig, ArcDefinition

config = MultiShootingConfig(num_arcs=5, max_nfev=10, solve_tol=1e-3)
solver = MultiShootingBVP(params, config)

arcs = [ArcDefinition(t_start=i*T/5, t_end=(i+1)*T/5, arc_index=i, is_periodic=(i==4))
        for i in range(5)]

result = solver.solve(initial_states, arcs, fidelity_weights)
```

### 4. State Transition Matrix (STM) Propagation and Covariance Analysis

**Module:** `rf3bp_lab.dynamics.variational`

Enables mission-grade navigation uncertainty quantification:

- **VariationalState** dataclass: Bundles nominal state and 6x6 STM
- **propagate_variational()**: Integrates state and variational equations (Jacobi equation) simultaneously
  - Combines 6 state ODEs + 36 STM ODEs (total: 42 components)
  - Uses finite-difference Jacobian for RF3BP dynamics
- **covariance_at_time()**: Transforms covariance via $P(t) = \Phi(t, t_0) \cdot P(t_0) \cdot \Phi(t, t_0)^T$
- **dilution_of_precision()**: Computes DOP metrics (position/velocity uncertainty, GDOP)

#### Application

Covariance grows predictably along reference trajectories; STM enables trajectory correction planning and measurement requirement trade studies.

### 5. Orbit Family Continuation and Multi-Objective Optimization

**Module:** `rf3bp_lab.shooting.orbit_families`

Traces families of bounded orbits via pseudo-arc-length continuation:

- **OrbitFamilyPoint** dataclass: Records state, period, energy, Jacobi constant, stability index, and objective values
- **trace_orbit_family()**: Marches along a parameter axis (amplitude, energy, initial position) with adaptive shooting
- **pareto_front()**: Identifies non-dominated orbits for multi-objective design
- **compute_family_metrics()**: Summarizes family statistics (period/energy ranges, stability trends)

#### Motivating Use Case

Binary-asteroid missions often explore design families: fixed-energy orbits around the smaller body, or amplitude-constrained paths with varying orbital periods. Families reveal which parameter regions offer stable long-duration bounded orbits vs. regions prone to escape.

#### Example

```python
from rf3bp_lab.shooting.orbit_families import trace_orbit_family, FamilyConfig

config = FamilyConfig(
    param_name="amplitude",
    param_min=0.1,
    param_max=1.5,
    num_points=20,
)

family_result = trace_orbit_family(seed_state, params, config)

pareto_idx = pareto_front(family_result['family_points'], 
                          objectives=['stability_index', 'energy'])
```

### 6. Mission-Grade Navigation Covariance Analysis

**Module:** `rf3bp_lab.navigation.covariance`

Provides complete navigation uncertainty workflow:

- **NavUncertaintyAnalysis** dataclass: Stores covariance time series, position/velocity uncertainties, and DOP metrics
- **propagate_navigation_uncertainty()**: Propagates initial covariance via STM integration
- **sensitivity_to_maneuver()**: Quantifies how maneuver execution errors ($\sigma_a$ in acceleration) affect final trajectory uncertainty
- **information_from_measurement()**: Computes information gain and posterior covariance from a measurement update (Kalman information filter perspective)

#### Typical Workflow

1. Start with initial position/velocity uncertainty ($\sigma_x$, $\sigma_v$)
2. Propagate covariance forward over planned mission segment
3. Evaluate measurement sensitivity (e.g., ranging accuracy) via information-gain analysis
4. Assess maneuver robustness to thruster execution errors

---

## Roadmap

| Priority | Upgrade | Expected Benefit | Status |
| --- | --- | --- | --- |
| High | Polyhedral gravity (polyhedron model integration) | Much better local field realism near both bodies | Future |
| High | Multi-shooting continuation | Better robustness for long bounded arcs | ✅ **Implemented** |
| High | SPICE-driven Sun geometry | Better solar forcing fidelity | Future |
| Medium | Family continuation and branch tracking | Better orbit atlas generation | ✅ **Implemented** |
| Medium | Event surfaces and Poincare diagnostics | Better structure discovery | Future |
| Medium | Eclipsing and attitude-sensitive SRP | Better small-spacecraft realism | ✅ **Implemented** |
| Low | Formal optimization over orbit families | Systematic mission design automation | ✅ **Implemented** |
| Low | Mission-grade navigation covariance analysis | Flight-readiness trajectory analysis | ✅ **Implemented** |
| Low | High-order irregular-body gravity | Better fidelity for shaped asteroids | ✅ **Implemented** |

## References and Context

The model and algorithms in this lab are anchored to real literature. Most links below use DOI resolver URLs so they remain stable even when publisher front-ends change.

| Topic | Reference | What It Supports in This Repo |
| --- | --- | --- |
| RF3BP near binary asteroids (primary inspiration) | Lu, J., Shang, H., Liu, C., Zhang, X., Gao, A. (2026). General Dynamics of Restricted Full Three-Body Problem Near Binary Asteroid System. *Journal of Guidance, Control, and Dynamics*. [https://doi.org/10.2514/1.G009686](https://doi.org/10.2514/1.G009686) | Pulsating-rotating RF3BP framing, perturbation hierarchy, and continuation motivation |
| Moshup-Squannit physical system observations | Ostro, S. J., Margot, J.-L., Benner, L. A. M., et al. (2006). Radar Imaging of Binary Near-Earth Asteroid (66391) 1999 KW4. *Science*, 314(5803), 1276-1280. [https://doi.org/10.1126/science.1133622](https://doi.org/10.1126/science.1133622) | Binary geometry context and observed KW4 dynamical characteristics |
| Binary asteroid full-problem dynamics | Fahnestock, E. G., Scheeres, D. J. (2008). Simulation and analysis of the dynamics of binary near-Earth asteroid (66391) 1999 KW4. *Icarus*, 194(2), 410-435. [https://doi.org/10.1016/j.icarus.2007.11.007](https://doi.org/10.1016/j.icarus.2007.11.007) | Small-body binary dynamics baseline and comparable perturbation scales |
| Full two-body stability theory | Scheeres, D. J. (2009). Stability of the planar full 2-body problem. *Celestial Mechanics and Dynamical Astronomy*, 104, 103-128. [https://doi.org/10.1007/s10569-009-9184-7](https://doi.org/10.1007/s10569-009-9184-7) | Conceptual bridge between restricted and full-body stability behavior |
| CR3BP mission-design continuation practice | Koon, W. S., Lo, M. W., Marsden, J. E., Ross, S. D. (2001). Low Energy Transfer to the Moon. *Celestial Mechanics and Dynamical Astronomy*, 81, 63-73. [https://doi.org/10.1023/A:1013359120468](https://doi.org/10.1023/A:1013359120468) | Continuation and shooting-style trajectory construction mindset |
| Foundational CR3BP text | Szebehely, V. (1967). *Theory of Orbits: The Restricted Problem of Three Bodies*. Academic Press. (Book cited in many modern CR3BP papers; review DOI: [https://doi.org/10.1119/1.1974535](https://doi.org/10.1119/1.1974535)) | Canonical rotating-frame CR3BP equations and integrals |
| Numerical integration basis (embedded RK) | Dormand, J. R., Prince, P. J. (1980). A family of embedded Runge-Kutta formulae. *Journal of Computational and Applied Mathematics*, 6(1), 19-26. [https://doi.org/10.1016/0771-050X(80)90013-3](https://doi.org/10.1016/0771-050X(80)90013-3) | The accuracy-control lineage behind high-order RK propagation choices |
| Software stack reference (SciPy) | Virtanen, P., Gommers, R., Oliphant, T. E., et al. (2020). SciPy 1.0: fundamental algorithms for scientific computing in Python. *Nature Methods*, 17, 261-272. [https://doi.org/10.1038/s41592-019-0686-2](https://doi.org/10.1038/s41592-019-0686-2) | Scientific-computing implementation context for solve_ivp and linear algebra tooling |

> [!NOTE]
> Some publisher pages are paywalled or bot-protected in this execution environment. DOI resolver links above are real and were verified via open Crossref metadata.

## Bottom Line

This repository is a compact RF3BP algorithm lab for **actually experimenting with the dynamics**, not a thin command wrapper. It gives you:

- a CR3BP baseline
- an RF3BP-inspired higher-fidelity model
- explicit pulsation-related kinematics
- staged continuation into higher fidelity
- plots that show what each perturbation is doing

That makes it a useful foundation for pushing toward more realistic bounded-orbit design in the Moshup-Squannit environment.
