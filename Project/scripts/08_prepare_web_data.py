#!/usr/bin/env python3
"""Prepare browser-streamable scientific data without committing NASA sources."""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path

import numpy as np

from lunar_gravity import BRILLOUIN_RADIUS_M, GRAILGravity
from lunar_terrain import LDEM64, META


ROOT = Path(__file__).resolve().parents[1]
WEB_ROOT = next(
    (candidate for candidate in (ROOT / "Web", ROOT.parent / "Web") if candidate.exists()),
    ROOT / "Web",
)
WEB_PUBLIC = WEB_ROOT / "public"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def write_gzip(path: Path, payload: bytes, level: int = 9) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as raw:
        with gzip.GzipFile(
            filename="",
            mode="wb",
            compresslevel=level,
            fileobj=raw,
            mtime=0,
        ) as compressed:
            compressed.write(payload)


def prepare_terrain(tile_degrees: int = 10) -> dict[str, object]:
    if 180 % tile_degrees or 360 % tile_degrees:
        raise ValueError("tile_degrees must divide both 180 and 360")
    terrain = LDEM64()
    pixels = int(tile_degrees * META.pixels_per_degree)
    lat_tiles = 180 // tile_degrees
    lon_tiles = 360 // tile_degrees
    target = WEB_PUBLIC / "data" / "terrain"
    target.mkdir(parents=True, exist_ok=True)
    total_compressed = 0
    records = []
    for lat_index in range(lat_tiles):
        south = -90 + lat_index * tile_degrees
        north = south + tile_degrees
        row_start = int((90 - north) * META.pixels_per_degree)
        row_end = row_start + pixels
        for lon_index in range(lon_tiles):
            west = lon_index * tile_degrees
            east = west + tile_degrees
            col_start = int(west * META.pixels_per_degree)
            col_end = col_start + pixels
            # Preserve the source north-to-south, west-to-east raster order.
            tile = np.ascontiguousarray(
                terrain.grid[row_start:row_end, col_start:col_end], dtype="<i2"
            )
            # Horizontal delta coding is exactly reversible in int16 modular
            # arithmetic and greatly reduces the static-site payload.
            encoded = np.empty_like(tile)
            encoded[:, 0] = tile[:, 0]
            encoded[:, 1:] = np.diff(tile.astype(np.int32), axis=1).astype("<i2")
            filename = f"t_{lat_index:02d}_{lon_index:02d}.i16.gz"
            destination = target / filename
            write_gzip(destination, encoded.tobytes(order="C"))
            size = destination.stat().st_size
            total_compressed += size
            records.append(
                {
                    "lat_index": lat_index,
                    "lon_index": lon_index,
                    "south_deg": south,
                    "north_deg": north,
                    "west_deg": west,
                    "east_deg": east,
                    "file": filename,
                    "compressed_bytes": size,
                }
            )
        print(
            f"terrain latitude band {south:+03d} to {north:+03d} deg complete",
            flush=True,
        )
    manifest: dict[str, object] = {
        "dataset": "LOLA LDEM64",
        "source": (
            "https://pds-geosciences.wustl.edu/lro/lro-l-lola-3-rdr-v1/"
            "lrolol_1xxx/data/lola_gdr/cylindrical/img/ldem_64.img"
        ),
        "source_sha256": (
            "98f1824b1a999630bf7b1f59575fe4c1d56a24722a968d98bcce07e6d6d44d4b"
        ),
        "reference_radius_m": META.reference_radius_m,
        "pixels_per_degree": META.pixels_per_degree,
        "tile_degrees": tile_degrees,
        "tile_rows": pixels,
        "tile_columns": pixels,
        "data_type": "little-endian signed int16",
        "predictor": "first sample absolute, then horizontal int16 deltas",
        "row_order": "north_to_south",
        "column_order": "west_to_east",
        "scale_m_per_dn": META.scale_m_per_dn,
        "compressed_encoding": "gzip",
        "total_compressed_bytes": total_compressed,
        "tile_count": len(records),
        "tiles": records,
    }
    manifest_path = target / "manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    manifest["manifest_sha256"] = sha256(manifest_path)
    return manifest


def _resample_dh_grid(
    grid: np.ndarray, target_pixels_per_degree: int
) -> np.ndarray:
    """Bilinearly resample a DH2 grid to pixel-centred cylindrical samples."""
    nlat, nlon = grid.shape
    if nlon != 2 * nlat:
        raise ValueError(f"Expected a DH2 grid, received {grid.shape}")
    target_rows = 180 * target_pixels_per_degree
    target_cols = 360 * target_pixels_per_degree
    latitudes = 90.0 - (
        np.arange(target_rows, dtype=float) + 0.5
    ) / target_pixels_per_degree
    longitudes = (
        np.arange(target_cols, dtype=float) + 0.5
    ) / target_pixels_per_degree
    row = (90.0 - latitudes) * nlat / 180.0
    col = longitudes * nlat / 180.0
    r0 = np.floor(row).astype(np.int64)
    c0 = np.floor(col).astype(np.int64)
    r1 = np.minimum(r0 + 1, nlat - 1)
    c1 = (c0 + 1) % nlon
    fr = (row - r0)[:, None]
    fc = (col - c0)[None, :]
    return (
        (1.0 - fr) * (1.0 - fc) * grid[r0[:, None], c0[None, :]]
        + (1.0 - fr) * fc * grid[r0[:, None], c1[None, :]]
        + fr * (1.0 - fc) * grid[r1[:, None], c0[None, :]]
        + fr * fc * grid[r1[:, None], c1[None, :]]
    )


def _atlas_interpolate(
    shells: list[np.ndarray],
    altitude_shells_m: list[float],
    altitude_m: float,
    latitude_deg: float,
    longitude_deg: float,
    pixels_per_degree: int,
) -> np.ndarray:
    altitude = float(np.clip(altitude_m, altitude_shells_m[0], altitude_shells_m[-1]))
    high = int(np.searchsorted(altitude_shells_m, altitude, side="right"))
    high = min(max(1, high), len(altitude_shells_m) - 1)
    low = high - 1
    fraction_alt = (altitude - altitude_shells_m[low]) / (
        altitude_shells_m[high] - altitude_shells_m[low]
    )
    lat = float(np.clip(latitude_deg, -90 + 0.5 / pixels_per_degree, 90 - 0.5 / pixels_per_degree))
    lon = longitude_deg % 360.0
    row = (90.0 - lat) * pixels_per_degree - 0.5
    col = lon * pixels_per_degree - 0.5
    r0 = int(np.floor(row))
    c0 = int(np.floor(col)) % (360 * pixels_per_degree)
    r1 = min(r0 + 1, 180 * pixels_per_degree - 1)
    c1 = (c0 + 1) % (360 * pixels_per_degree)
    fr = row - np.floor(row)
    fc = col - np.floor(col)

    def horizontal(grid: np.ndarray) -> np.ndarray:
        return (
            (1 - fr) * (1 - fc) * grid[r0, c0]
            + (1 - fr) * fc * grid[r0, c1]
            + fr * (1 - fc) * grid[r1, c0]
            + fr * fc * grid[r1, c1]
        )

    a = horizontal(shells[low])
    b = horizontal(shells[high])
    return a + fraction_alt * (b - a)


def prepare_gravity(
    degree: int = 600,
    tile_degrees: int = 15,
    pixels_per_degree: int = 8,
) -> dict[str, object]:
    """Create a lazy-loadable degree-N acceleration atlas for the Web Worker."""
    if 180 % tile_degrees or 360 % tile_degrees:
        raise ValueError("tile_degrees must divide both 180 and 360")
    if tile_degrees * pixels_per_degree <= 1:
        raise ValueError("gravity tiles need multiple samples")
    from pyshtools.backends import shtools

    altitude_shells = [
        BRILLOUIN_RADIUS_M - META.reference_radius_m,
        12_000.0,
        15_000.0,
        20_000.0,
        30_000.0,
        50_000.0,
        80_000.0,
        120_000.0,
        200_000.0,
        320_000.0,
    ]
    gravity = GRAILGravity(maximum_degree=degree)
    target = WEB_PUBLIC / "data" / "gravity"
    target.mkdir(parents=True, exist_ok=True)
    rows = 180 * pixels_per_degree
    columns = 360 * pixels_per_degree
    tile_pixels = tile_degrees * pixels_per_degree
    lat_tiles = 180 // tile_degrees
    lon_tiles = 360 // tile_degrees
    shell_grids: list[np.ndarray] = []
    total_compressed = 0
    tile_count = 0

    for shell_index, altitude in enumerate(altitude_shells):
        radius = META.reference_radius_m + altitude
        radial, theta, phi, _, _ = shtools.MakeGravGridDH(
            gravity.coefficients.coeffs,
            gravity.coefficients.gm,
            gravity.coefficients.r0,
            a=radius,
            f=0.0,
            lmax=degree,
            sampling=2,
            lmax_calc=degree,
            omega=0.0,
            normal_gravity=0,
            extend=False,
        )
        components = np.empty((rows, columns, 3), dtype=np.float32)
        # Store only the non-central correction. The browser evaluates the
        # dominant -GM/r^2 term analytically; this makes interpolation between
        # altitude shells far more accurate and improves compression.
        components[:, :, 0] = (
            _resample_dh_grid(radial, pixels_per_degree)
            + gravity.coefficients.gm / radius**2
        ).astype(np.float32)
        components[:, :, 1] = _resample_dh_grid(
            theta, pixels_per_degree
        ).astype(np.float32)
        components[:, :, 2] = _resample_dh_grid(
            phi, pixels_per_degree
        ).astype(np.float32)
        shell_grids.append(components)
        del radial, theta, phi

        for lat_index in range(lat_tiles):
            # Browser tile index increases from south to north, while arrays
            # run from north to south.
            source_lat_index = lat_tiles - 1 - lat_index
            row_start = source_lat_index * tile_pixels
            row_end = row_start + tile_pixels
            for lon_index in range(lon_tiles):
                col_start = lon_index * tile_pixels
                col_end = col_start + tile_pixels
                tile = np.ascontiguousarray(
                    components[row_start:row_end, col_start:col_end, :],
                    dtype="<f4",
                )
                filename = (
                    f"g_{shell_index:02d}_{lat_index:02d}_{lon_index:02d}.f32.gz"
                )
                destination = target / filename
                write_gzip(destination, tile.tobytes(order="C"))
                total_compressed += destination.stat().st_size
                tile_count += 1
        print(
            f"gravity shell {shell_index + 1}/{len(altitude_shells)} "
            f"at {altitude / 1000:g} km complete",
            flush=True,
        )

    # Validate the browser atlas interpolation against direct degree-N values.
    rng = np.random.default_rng(1600)
    errors = []
    relative_errors = []
    for _ in range(120):
        altitude = float(rng.uniform(altitude_shells[0], altitude_shells[-1]))
        latitude = float(rng.uniform(-89.5, 89.5))
        longitude = float(rng.uniform(0.0, 360.0))
        radius = META.reference_radius_m + altitude
        lat_r = np.radians(latitude)
        lon_r = np.radians(longitude)
        body = radius * np.array(
            [
                np.cos(lat_r) * np.cos(lon_r),
                np.cos(lat_r) * np.sin(lon_r),
                np.sin(lat_r),
            ]
        )
        direct = gravity.body_acceleration(body, degree=degree)
        spherical = _atlas_interpolate(
            shell_grids,
            altitude_shells,
            altitude,
            latitude,
            longitude,
            pixels_per_degree,
        )
        spherical[0] -= gravity.coefficients.gm / radius**2
        radial_basis = body / radius
        east = np.array([-np.sin(lon_r), np.cos(lon_r), 0.0])
        north = np.array(
            [
                -np.sin(lat_r) * np.cos(lon_r),
                -np.sin(lat_r) * np.sin(lon_r),
                np.cos(lat_r),
            ]
        )
        atlas = (
            spherical[0] * radial_basis
            - spherical[1] * north
            + spherical[2] * east
        )
        error = float(np.linalg.norm(atlas - direct))
        errors.append(error)
        relative_errors.append(error / float(np.linalg.norm(direct)))

    validation = {
        "sample_count": len(errors),
        "random_seed": 1600,
        "rms_acceleration_error_m_s2": float(np.sqrt(np.mean(np.square(errors)))),
        "maximum_acceleration_error_m_s2": float(np.max(errors)),
        "rms_relative_acceleration_error": float(
            np.sqrt(np.mean(np.square(relative_errors)))
        ),
        "maximum_relative_acceleration_error": float(np.max(relative_errors)),
    }
    manifest: dict[str, object] = {
        "model": "GRAIL GRGM1200B",
        "degree": degree,
        "source": "https://pgda.gsfc.nasa.gov/data/MoonRM1/sha.grgm1200b_sigma",
        "source_sha256": (
            "f08a988b43f3eaa5a2089045a9b7e41e02f16542c7912b87ea34366fafa39bc5"
        ),
        "reference_radius_m": META.reference_radius_m,
        "gravity_coefficient_reference_radius_m": gravity.coefficients.r0,
        "gm_m3_s2": gravity.coefficients.gm,
        "stores_noncentral_correction": True,
        "brillouin_radius_m": BRILLOUIN_RADIUS_M,
        "pixels_per_degree": pixels_per_degree,
        "tile_degrees": tile_degrees,
        "tile_rows": tile_pixels,
        "tile_columns": tile_pixels,
        "altitude_shells_m": altitude_shells,
        "component_order": ["radial", "theta", "phi"],
        "data_type": "little-endian float32, components interleaved",
        "row_order": "north_to_south",
        "column_order": "west_to_east",
        "compressed_encoding": "gzip",
        "tile_count": tile_count,
        "total_compressed_bytes": total_compressed,
        "validation": validation,
    }
    manifest_path = target / "manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    (ROOT / "data" / "output" / "web_gravity_atlas_validation.json").write_text(
        json.dumps(validation, indent=2), encoding="utf-8"
    )
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--terrain", action="store_true")
    parser.add_argument("--gravity", action="store_true")
    parser.add_argument("--tile-degrees", type=int, default=10)
    parser.add_argument("--gravity-degree", type=int, default=600)
    args = parser.parse_args()
    if not args.terrain and not args.gravity:
        parser.error("choose --terrain and/or --gravity")
    if args.terrain:
        manifest = prepare_terrain(args.tile_degrees)
        print(
            f"Wrote {manifest['tile_count']} terrain tiles, "
            f"{manifest['total_compressed_bytes'] / 1e6:.1f} MB compressed"
        )
    if args.gravity:
        manifest = prepare_gravity(degree=args.gravity_degree)
        print(
            f"Wrote {manifest['tile_count']} gravity tiles, "
            f"{manifest['total_compressed_bytes'] / 1e6:.1f} MB compressed"
        )


if __name__ == "__main__":
    main()
