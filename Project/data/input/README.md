# Input data

## `physical_constants.json`

Small, human-readable constants used by every default script. SI units are encoded in key names. Principal sources are JPL/DE440-compatible lunar parameters, NASA/IAU lunar rotation, LOLA/LROC terrain scales, official baseball dimensions, and the cited Statcast speed benchmark.

## `gravity_degree2.csv`

Representative **unnormalised** principal-axis degree-2 coefficients used in the transparent M3 model:

- `J2 = 2.033e-4`
- `C22 = 2.24e-5`

These are not presented as a substitute for GRGM900C/GRGM1200A. Their purpose is to show, with a small independently reproducible model, that Kepler closure is structurally broken by lunar non-sphericity.

## `external_data_manifest.csv`

Official NASA PDS URLs for optional LOLA LDEM4 and GRAIL GRGM1200A files. Run `scripts/download_optional_data.py` to download selected entries into this directory.

The paper's numbered figures do not require the optional files. If LDEM4 is present, `05_terrain_analysis.py` also produces a topographic map and a sampled body-fixed ground-track clearance table/figure.
