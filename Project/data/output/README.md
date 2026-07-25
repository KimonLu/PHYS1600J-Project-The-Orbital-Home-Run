# Generated output data

Scripts 01--07 generate the ideal, sensitivity, rotation, degree-2
perturbation, legacy LDEM4, validation, and controlled targeted-return files.

The high-fidelity workflow adds:

- `gravity_degree_convergence.csv`: identical-state GRGM1200B truncation runs
  at 12, 15, and 30 km.
- `gravity_model_selection.csv`: explicit degree-600 selection against a
  degree-1200 reference and a 10 m criterion.
- `gravity_convergence_run_manifest.json`: parameters and source provenance.
- `web_gravity_atlas_validation.json`: random-point acceleration comparison.
- `web_solver_crosscheck.json`: complete atlas-versus-direct trajectories.
- `web_terrain_tile_validation.json`: lossless LDEM64 decode and bilinear
  interpolation checks across tile and longitude seams.
- `high_fidelity_case_summary.json`: degree-600 targeted-return definition and
  terminal/clearance metrics.
- `high_fidelity_case_trajectory.csv`: complete inertial and body-fixed
  trajectory used by the paper.
- `high_fidelity_case_sensitivity.json`: projected local launch-state
  sensitivity.
- `validation_report.txt`: analytic, numerical, terrain, high-degree, and web
  consistency tests.

Raw LDEM64 and GRGM1200B are stored under `data/external/` and are not outputs
or Git artifacts. Run `scripts/download_science_data.py` to recover them from
the official sources with size and SHA-256 verification.
