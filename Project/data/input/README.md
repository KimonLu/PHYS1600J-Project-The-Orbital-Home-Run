# Input data

## `physical_constants.json`

Small, human-readable constants used by every default script. SI units are encoded in key names. Principal sources are JPL/DE440-compatible lunar parameters, NASA/IAU lunar rotation, LOLA/LROC terrain scales, official baseball dimensions, and the cited Statcast speed benchmark.

## `gravity_degree2.csv`

**Unnormalised** principal-axis degree-2 coefficients used in the transparent M3 model:

- `J2 = 2.033560627118353e-4`
- `C22 = 2.239820214219053e-5`

The source JPL DE430 values use a 1738 km gravity reference radius. They are rescaled to the project's 1737.4 km computational radius so that each degree-2 product \(C_{2m}R^2\) remains invariant. The core propagator reads this CSV directly rather than duplicating the values in code.

These coefficients are not presented as a substitute for GRGM900C/GRGM1200A. Their purpose is to show, with a small independently reproducible model, that Kepler closure is structurally broken by lunar non-sphericity.

## Authoritative large products

The main realistic solver uses:

- **LOLA LDEM64**, a global, pixel-registered \(11520\times23040\) signed
  16-bit elevation grid at 64 pixels per degree (about 474 m per pixel at the
  equator), relative to the 1737.4 km reference sphere.
- **GRAIL GRGM1200B**, 4-pi geodesy-normalized lunar gravity coefficients
  through degree and order 1200, with formal coefficient uncertainties.

These source products are 531 MB and 83 MB respectively and are not committed
to Git. Run:

```text
python scripts/download_science_data.py
```

The downloader writes to `data/external/`, resumes partial transfers, checks
the exact byte counts, and verifies SHA-256. `external_data_manifest.csv`
records the official URLs, measured hashes, and scientific purpose.

LDEM4 remains only for reproducing the earlier low-resolution equatorial
demonstration. It is no longer the collision surface for the main solver.

LOLA also publishes LDEM128, LDEM256, and regional products at still finer
sampling. They are valuable for landing-site or local-hazard studies, but
LDEM64 is already finer than the defensible global trajectory uncertainty in
this project; increasing terrain sampling alone would not cure gravity,
orientation, ephemeris, rock-scale, or launch-state uncertainty.
