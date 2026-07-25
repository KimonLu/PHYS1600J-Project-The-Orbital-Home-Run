#!/usr/bin/env python3
"""Cross-check browser LDEM64 tiles, including tile and longitude seams."""
from __future__ import annotations

import gzip
import json
from pathlib import Path

import numpy as np

from lunar_terrain import LDEM64, META


ROOT = Path(__file__).resolve().parents[1]
WEB_ROOT = next(
    (candidate for candidate in (ROOT / "Web", ROOT.parent / "Web") if candidate.exists()),
    ROOT / "Web",
)
TERRAIN_DIR = WEB_ROOT / "public" / "data" / "terrain"
OUTPUT = ROOT / "data" / "output" / "web_terrain_tile_validation.json"


class TiledTerrain:
    def __init__(self) -> None:
        self.meta = json.loads(
            (TERRAIN_DIR / "manifest.json").read_text(encoding="utf-8")
        )
        self.cache: dict[tuple[int, int], np.ndarray] = {}

    def tile(self, latitude_index: int, longitude_index: int) -> np.ndarray:
        key = latitude_index, longitude_index
        if key not in self.cache:
            filename = f"t_{latitude_index:02d}_{longitude_index:02d}.i16.gz"
            encoded = np.frombuffer(
                gzip.decompress((TERRAIN_DIR / filename).read_bytes()),
                dtype="<i2",
            ).reshape(self.meta["tile_rows"], self.meta["tile_columns"])
            # Reproduce JavaScript int16 modular cumulative summation.
            decoded = np.cumsum(encoded.astype(np.int32), axis=1).astype("<i2")
            self.cache[key] = decoded
        return self.cache[key]

    def sample_dn(self, row: int, column: int) -> int:
        total_rows = int(180 * self.meta["pixels_per_degree"])
        total_columns = int(360 * self.meta["pixels_per_degree"])
        row = int(np.clip(row, 0, total_rows - 1))
        column %= total_columns
        band_from_north = row // int(self.meta["tile_rows"])
        latitude_index = (
            int(180 // self.meta["tile_degrees"]) - 1 - band_from_north
        )
        longitude_index = column // int(self.meta["tile_columns"])
        local_row = row % int(self.meta["tile_rows"])
        local_column = column % int(self.meta["tile_columns"])
        return int(self.tile(latitude_index, longitude_index)[local_row, local_column])

    def elevation_m(self, latitude_deg: float, longitude_deg: float) -> float:
        pixels_per_degree = self.meta["pixels_per_degree"]
        latitude = float(
            np.clip(
                latitude_deg,
                -90 + 0.5 / pixels_per_degree,
                90 - 0.5 / pixels_per_degree,
            )
        )
        longitude = longitude_deg % 360.0
        row = (90 - latitude) * pixels_per_degree - 0.5
        column = longitude * pixels_per_degree - 0.5
        r0 = int(np.floor(row))
        c0 = int(np.floor(column))
        r1 = min(r0 + 1, int(180 * pixels_per_degree - 1))
        c1 = (c0 + 1) % int(360 * pixels_per_degree)
        fr = row - r0
        fc = column - np.floor(column)
        z00 = self.sample_dn(r0, c0)
        z01 = self.sample_dn(r0, c1)
        z10 = self.sample_dn(r1, c0)
        z11 = self.sample_dn(r1, c1)
        return float(
            (
                (1 - fr) * (1 - fc) * z00
                + (1 - fr) * fc * z01
                + fr * (1 - fc) * z10
                + fr * fc * z11
            )
            * self.meta["scale_m_per_dn"]
        )


def main() -> None:
    direct = LDEM64()
    tiled = TiledTerrain()
    rng = np.random.default_rng(1600)
    points = [
        (float(rng.uniform(-89.99, 89.99)), float(rng.uniform(-720, 720)))
        for _ in range(400)
    ]
    # Exercise both sides and exact locations of every tile seam.
    for latitude in range(-80, 90, 10):
        for offset in (-1e-7, 0.0, 1e-7):
            points.append((latitude + offset, 123.456789))
    for longitude in range(0, 361, 10):
        for offset in (-1e-7, 0.0, 1e-7):
            points.append((12.345678, longitude + offset))

    differences = np.asarray(
        [
            tiled.elevation_m(latitude, longitude)
            - direct.elevation_m(latitude, longitude)
            for latitude, longitude in points
        ]
    )
    report = {
        "sample_count": len(points),
        "random_seed": 1600,
        "maximum_absolute_difference_m": float(np.max(np.abs(differences))),
        "rms_difference_m": float(np.sqrt(np.mean(differences**2))),
        "tile_count": tiled.meta["tile_count"],
        "interpretation": (
            "Lossless tile decode and bilinear interpolation reproduce direct "
            "LDEM64 queries, including 10-degree seams and 0/360 wrapping."
        ),
    }
    if report["maximum_absolute_difference_m"] > 1e-9:
        raise AssertionError(report)
    OUTPUT.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
