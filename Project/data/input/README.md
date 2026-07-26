# Input data

## Versioned inputs

`physical_constants.json` contains the small SI-unit constants used by the
analytic and numerical models. The global terrain-envelope height is the
directly scanned LDEM64 maximum of 10.757 km relative to the 1737.4 km sphere.

`gravity_degree2.csv` contains the unnormalised principal-axis degree-2
coefficients used only by the transparent perturbation hierarchy. The
high-fidelity solver does not substitute these coefficients for GRGM1200B.

`external_data_manifest.csv` records the official URLs, exact byte counts,
SHA-256 digests, local paths, purposes, and provenance of the large NASA
products.

## Untracked authoritative products

The high-fidelity workflow requires:

- LOLA LDEM64: a 11520 by 23040 signed 16-bit global elevation grid at
  64 pixels per degree and 0.5 m vertical quantization;
- GRAIL GRGM1200B: 4-pi geodesy-normalized gravity coefficients and formal
  uncertainties through degree and order 1200.

Download and verify both products with:

```text
python scripts/download_science_data.py
```

They are written under `data/external/`, which is intentionally ignored by
Git. Higher-resolution LOLA products exist, but terrain sampling alone does
not remove gravity, frame, ephemeris, launch-state, or rock-scale limitations.
