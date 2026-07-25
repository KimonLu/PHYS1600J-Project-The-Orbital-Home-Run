"""Vectorized access to the global LOLA LDEM64 elevation model."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np


ROOT = Path(__file__).resolve().parents[1]


@dataclass(frozen=True)
class TerrainMetadata:
    rows: int = 11_520
    columns: int = 23_040
    pixels_per_degree: float = 64.0
    scale_m_per_dn: float = 0.5
    reference_radius_m: float = 1_737_400.0
    horizontal_pixel_equator_m: float = 473.8
    frame: str = "Mean Earth / Polar Axis of DE421"


META = TerrainMetadata()


class LDEM64:
    """Memory-mapped, bilinearly interpolated LDEM64 terrain.

    The product is pixel-registered.  Row centres run from
    ``90 - 0.5/64`` to ``-90 + 0.5/64`` degrees; column centres run eastward
    from ``0.5/64`` to ``360 - 0.5/64`` degrees.  Interpolation is periodic in
    longitude and clamped at the polar pixel centres.
    """

    def __init__(self, path: Path | str | None = None) -> None:
        if path is None:
            path = ROOT / "data" / "external" / "ldem64" / "ldem_64.img"
        self.path = Path(path)
        expected = META.rows * META.columns * np.dtype("<i2").itemsize
        if not self.path.exists():
            raise FileNotFoundError(
                f"LDEM64 not found at {self.path}. Run download_science_data.py."
            )
        if self.path.stat().st_size != expected:
            raise ValueError(
                f"Invalid LDEM64 size: {self.path.stat().st_size:,}; "
                f"expected {expected:,} bytes"
            )
        self.grid = np.memmap(
            self.path,
            mode="r",
            dtype="<i2",
            shape=(META.rows, META.columns),
        )

    @staticmethod
    def _fractional_indices(
        latitude_deg: np.ndarray, longitude_deg: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray]:
        lat = np.clip(
            latitude_deg,
            -90.0 + 0.5 / META.pixels_per_degree,
            90.0 - 0.5 / META.pixels_per_degree,
        )
        lon = np.mod(longitude_deg, 360.0)
        row = (90.0 - lat) * META.pixels_per_degree - 0.5
        col = lon * META.pixels_per_degree - 0.5
        return row, col

    def elevation_m(
        self,
        latitude_deg: float | Iterable[float] | np.ndarray,
        longitude_deg: float | Iterable[float] | np.ndarray,
    ) -> float | np.ndarray:
        """Return bilinearly interpolated elevation above the 1737.4 km sphere."""
        lat, lon = np.broadcast_arrays(
            np.asarray(latitude_deg, dtype=float),
            np.asarray(longitude_deg, dtype=float),
        )
        row, col = self._fractional_indices(lat, lon)
        r0 = np.floor(row).astype(np.int64)
        c0 = np.floor(col).astype(np.int64)
        r1 = np.minimum(r0 + 1, META.rows - 1)
        c1 = np.mod(c0 + 1, META.columns)
        c0 = np.mod(c0, META.columns)
        fr = row - r0
        fc = col - np.floor(col)

        z00 = np.asarray(self.grid[r0, c0], dtype=float)
        z01 = np.asarray(self.grid[r0, c1], dtype=float)
        z10 = np.asarray(self.grid[r1, c0], dtype=float)
        z11 = np.asarray(self.grid[r1, c1], dtype=float)
        value_dn = (
            (1.0 - fr) * (1.0 - fc) * z00
            + (1.0 - fr) * fc * z01
            + fr * (1.0 - fc) * z10
            + fr * fc * z11
        )
        value = value_dn * META.scale_m_per_dn
        return float(value) if value.ndim == 0 else value

    def nearest_elevation_m(
        self, latitude_deg: float, longitude_deg: float
    ) -> float:
        row, col = self._fractional_indices(
            np.asarray(latitude_deg), np.asarray(longitude_deg)
        )
        rr = int(np.clip(np.rint(row), 0, META.rows - 1))
        cc = int(np.rint(col)) % META.columns
        return float(self.grid[rr, cc]) * META.scale_m_per_dn

    def extrema(self, rows_per_chunk: int = 256) -> tuple[float, float]:
        """Compute exact global extrema without loading the 531 MB grid at once."""
        minimum = np.iinfo(np.int16).max
        maximum = np.iinfo(np.int16).min
        for start in range(0, META.rows, rows_per_chunk):
            block = self.grid[start : start + rows_per_chunk]
            minimum = min(minimum, int(block.min()))
            maximum = max(maximum, int(block.max()))
        return minimum * META.scale_m_per_dn, maximum * META.scale_m_per_dn


def body_fixed_lat_lon_deg(position_body_m: np.ndarray) -> tuple[float, float]:
    radius = float(np.linalg.norm(position_body_m))
    if radius == 0.0:
        raise ValueError("Position cannot be the lunar centre")
    lat = np.degrees(np.arcsin(position_body_m[2] / radius))
    lon = np.degrees(np.arctan2(position_body_m[1], position_body_m[0])) % 360.0
    return float(lat), float(lon)
