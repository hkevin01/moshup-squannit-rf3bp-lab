# Moshup-Squannit RF3BP Lab

Python project to experiment with bounded spacecraft motion near the binary asteroid system Moshup-Squannit (1999 KW4), inspired by RF3BP pulsating-rotating formulations.

## What This Project Implements

- CR3BP baseline equations in a rotating frame
- A higher-fidelity RF3BP-inspired model with:
  - nonuniform pulsation through time-varying primary-secondary separation
  - nonspherical gravity approximation via J2-like terms for both bodies
  - third-body solar gravity perturbation
  - solar radiation pressure (SRP)
- Hierarchical shooting algorithm:
  - Stage 1: correct seed orbit in CR3BP
  - Stage 2: continue solution into pulsating + perturbation model

## Important Scope Note

The source journal article is access restricted in this environment. This project is a research sandbox that follows the abstract-level model ingredients and standard CR3BP/RF3BP conventions, not a line-by-line reproduction of proprietary equations.

## Moshup-Squannit Defaults

The default constants are approximate, based on public sources (JPL radar summary and open references):

- Primary diameter roughly 1.3 to 1.5 km
- Secondary scale roughly 0.4 to 0.6 km
- Mutual orbit radius roughly 2.5 km
- Mutual period roughly 17.4 h

These are configured in `SystemParams` and can be refined as better shape/gravity solutions are added.

## Quick Start

```bash
cd /home/kevin/Projects/moshup-squannit-rf3bp-lab
python3 -m venv .venv
source .venv/bin/activate
pip install -e .[dev]
pytest -q
python scripts/run_demo.py
```

## Structure

- `src/rf3bp_lab/dynamics/models.py`: CR3BP and RF3BP-inspired dynamics
- `src/rf3bp_lab/shooting/hierarchical.py`: staged shooting/continuation solver
- `scripts/run_demo.py`: end-to-end orbit generation and comparison
- `tests/`: smoke tests for dynamics and shooting

## Next Upgrades

- Replace J2 approximation with polyhedral gravity using shape models
- Add explicit primary-secondary relative acceleration/jerk terms from potential derivatives
- Add finite-time bounded-orbit search over multi-shooting segments
- Add ephemeris-driven Sun geometry and eclipsing for SRP
