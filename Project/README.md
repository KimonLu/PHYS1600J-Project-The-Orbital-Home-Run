# PHYS1600J Problem B: The Orbital Home Run

This directory contains the numerical models, derived data, figures, and an
archived paper build for the PHYS1600J group project. The current scientific
workflow uses LOLA LDEM64 terrain and GRAIL GRGM1200B gravity. The earlier
LDEM4 and degree-2 targeted-return demonstrations have been removed because
they are superseded by the high-fidelity workflow.

## Reproduction profiles

Install the Python dependencies from the project directory:

```text
python -m pip install -r scripts/requirements.txt
```

The lightweight analytic and degree-2 figures require no large downloads:

```text
python scripts/generate_all.py --profile core
```

Download the authoritative LDEM64 and GRGM1200B sources before running the
high-fidelity or web-data profiles:

```text
python scripts/download_science_data.py
python scripts/generate_all.py --profile science
python scripts/generate_all.py --profile web
```

`--profile full` runs every generator. The bounded height continuation can be
slow; it intentionally tests only 2 km and higher release heights.

To validate the packaged derived results without regenerating them:

```text
python scripts/generate_all.py --profile validate
```

The full LDEM64 extrema check is reported as skipped when the untracked NASA
source file is absent. All checks based on committed derived outputs still run.

## Script layout

Numbered scripts are ordered by model development:

1. `01_ideal_models.py` -- analytic two-body orbit family and speed scale.
2. `02_sensitivity.py` -- ideal launch-angle, height, and Monte Carlo tests.
3. `03_rotation_and_resonance.py` -- lunar rotation and return resonance.
4. `04_realistic_perturbations.py` -- degree-2 and fixed third-body hierarchy.
5. `05_terrain_envelope.py` -- conservative analytic LDEM64 elevation envelope.
6. `06_validation.py` -- central-model integration convergence.
7. `07_gravity_convergence.py` -- GRGM1200B truncation study.
8. `08_prepare_web_data.py` -- browser terrain and gravity products.
9. `09_validate_web_solver.py` -- browser-atlas trajectory cross-check.
10. `10_high_fidelity_case.py` -- direct degree-600 boundary-value correction.
11. `11_case_sensitivity.py` -- fixed-time and free-time Jacobian convergence.
12. `12_validate_web_terrain.py` -- lossless terrain-tile checks.
13. `13_terrain_visualizations.py` -- LDEM64 map and trajectory corridor.
14. `14_surface_feasibility.py` -- bounded great-circle terrain sampling.
15. `15_high_fidelity_validation.py` -- integrator and gravity-order checks.
16. `16_height_continuation.py` -- bounded continuation at 2 km and above.

Shared modules use descriptive names: `orbital_home_run.py`,
`lunar_gravity.py`, `lunar_terrain.py`, `general_solver.py`, and `plotting.py`.
`validate_results.py` checks the packaged cross-script results.

## Data layout

- `data/input/` contains small versioned constants and the external-data
  manifest.
- `data/external/` contains downloaded NASA products and is ignored by Git.
- `data/output/` contains versioned, machine-readable derived results.
- `Web/public/data/` contains browser-ready tiles and gravity-atlas files.

All high-fidelity case products start with `high_fidelity_`; bounded
near-surface products start with `near_surface_`; browser cross-check products
start with `web_`. Temporary logs, errors, caches, and smoke-test trajectories
are ignored or excluded from the repository.

The nominal return is a deterministic model result released 30.000 km above
the reference sphere and 19.243 km above the local bilinearly interpolated
LDEM64 surface. Millimetre-scale endpoint values are boundary-value solver
residuals, not physical trajectory accuracy.
