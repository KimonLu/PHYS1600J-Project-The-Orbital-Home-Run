#!/usr/bin/env python3
"""Cross-check the browser gravity atlas against direct degree-600 propagation."""
from __future__ import annotations

import gzip
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp

from orbital_home_run import C, rot_z


ROOT = Path(__file__).resolve().parents[1]
WEB_ROOT = next(
    (candidate for candidate in (ROOT / "Web", ROOT.parent / "Web") if candidate.exists()),
    ROOT / "Web",
)
ATLAS = WEB_ROOT / "public" / "data" / "gravity"


class Atlas:
    def __init__(self) -> None:
        self.meta = json.loads((ATLAS / "manifest.json").read_text(encoding="utf-8"))
        self.cache: dict[tuple[int, int, int], np.ndarray] = {}

    def tile(self, shell: int, lat_index: int, lon_index: int) -> np.ndarray:
        key = (shell, lat_index, lon_index)
        if key not in self.cache:
            filename = f"g_{shell:02d}_{lat_index:02d}_{lon_index:02d}.f32.gz"
            with gzip.open(ATLAS / filename, "rb") as handle:
                data = np.frombuffer(handle.read(), dtype="<f4")
            size = self.meta["tile_rows"]
            self.cache[key] = data.reshape(size, size, 3)
        return self.cache[key]

    def sample_shell(
        self, shell: int, latitude: float, longitude: float
    ) -> np.ndarray:
        m = self.meta
        ppd = m["pixels_per_degree"]
        tile_deg = m["tile_degrees"]
        lat = float(np.clip(latitude, -90 + 0.5 / ppd, 90 - 0.5 / ppd))
        lon = longitude % 360.0
        lat_index = min(180 // tile_deg - 1, int((lat + 90) // tile_deg))
        lon_index = int(lon // tile_deg)
        north = -90 + (lat_index + 1) * tile_deg
        west = lon_index * tile_deg
        row = (north - lat) * ppd - 0.5
        col = (lon - west) * ppd - 0.5
        r0, c0 = math.floor(row), math.floor(col)
        fr, fc = row - r0, col - c0

        def one(rr: int, cc: int) -> np.ndarray:
            li, lo, r, c = lat_index, lon_index, rr, cc
            if r < 0:
                li = min(180 // tile_deg - 1, li + 1)
                r += m["tile_rows"]
            elif r >= m["tile_rows"]:
                li = max(0, li - 1)
                r -= m["tile_rows"]
            if c < 0:
                lo = (lo - 1) % (360 // tile_deg)
                c += m["tile_columns"]
            elif c >= m["tile_columns"]:
                lo = (lo + 1) % (360 // tile_deg)
                c -= m["tile_columns"]
            return self.tile(shell, li, lo)[r, c].astype(float)

        z00, z01, z10, z11 = (
            one(r0, c0),
            one(r0, c0 + 1),
            one(r0 + 1, c0),
            one(r0 + 1, c0 + 1),
        )
        return (
            (1 - fr) * (1 - fc) * z00
            + (1 - fr) * fc * z01
            + fr * (1 - fc) * z10
            + fr * fc * z11
        )

    def spherical(self, radius: float, latitude: float, longitude: float) -> np.ndarray:
        m = self.meta
        altitude = radius - m["reference_radius_m"]
        shells = m["altitude_shells_m"]
        altitude = float(np.clip(altitude, shells[0], shells[-1]))
        high = int(np.searchsorted(shells, altitude, side="right"))
        high = min(max(1, high), len(shells) - 1)
        low = high - 1
        fraction = (altitude - shells[low]) / (shells[high] - shells[low])
        value = self.sample_shell(low, latitude, longitude)
        if high != low:
            value = value + fraction * (
                self.sample_shell(high, latitude, longitude) - value
            )
        value[0] -= m["gm_m3_s2"] / radius**2
        return value

    def acceleration_body(self, position: np.ndarray) -> np.ndarray:
        radius = float(np.linalg.norm(position))
        latitude = math.degrees(math.asin(position[2] / radius))
        longitude = math.degrees(math.atan2(position[1], position[0])) % 360.0
        gr, gtheta, gphi = self.spherical(radius, latitude, longitude)
        lat, lon = math.radians(latitude), math.radians(longitude)
        radial = np.array(
            [math.cos(lat) * math.cos(lon), math.cos(lat) * math.sin(lon), math.sin(lat)]
        )
        north = np.array(
            [-math.sin(lat) * math.cos(lon), -math.sin(lat) * math.sin(lon), math.cos(lat)]
        )
        east = np.array([-math.sin(lon), math.cos(lon), 0.0])
        return gr * radial - gtheta * north + gphi * east


def main() -> None:
    table = pd.read_csv(ROOT / "data" / "output" / "gravity_degree_convergence.csv")
    results = []
    for altitude in sorted(table["altitude_m"].unique()):
        reference = table[
            (table["altitude_m"] == altitude) & (table["degree"] == 600)
        ].iloc[0]
        duration = float(reference["duration_s"])
        radius = C.radius + float(altitude)
        speed_relative = float(reference["surface_relative_speed_m_s"])
        state0 = np.array([radius, 0.0, 0.0, 0.0, speed_relative, 0.0])
        atlas = Atlas()
        omega = np.array([0.0, 0.0, C.omega])

        def rhs(_: float, state: np.ndarray) -> np.ndarray:
            position, velocity = state[:3], state[3:]
            acceleration = (
                atlas.acceleration_body(position)
                - 2.0 * np.cross(omega, velocity)
                - np.cross(omega, np.cross(omega, position))
            )
            return np.concatenate([velocity, acceleration])

        solution = solve_ivp(
            rhs,
            (0.0, duration),
            state0,
            method="DOP853",
            max_step=10.0,
            rtol=2e-10,
            atol=np.array([2e-4, 2e-4, 2e-4, 2e-7, 2e-7, 2e-7]),
        )
        if not solution.success:
            raise RuntimeError(solution.message)
        final_body = solution.y[:, -1]
        body_to_inertial = rot_z(C.omega * duration)
        position_i = body_to_inertial @ final_body[:3]
        velocity_i = body_to_inertial @ (
            final_body[3:] + np.cross(omega, final_body[:3])
        )
        reference_position = reference[
            ["final_x_m", "final_y_m", "final_z_m"]
        ].to_numpy(float)
        reference_velocity = reference[
            ["final_vx_m_s", "final_vy_m_s", "final_vz_m_s"]
        ].to_numpy(float)
        results.append(
            {
                "case_altitude_m": float(altitude),
                "degree": 600,
                "position_difference_atlas_vs_direct_m": float(
                    np.linalg.norm(position_i - reference_position)
                ),
                "velocity_difference_atlas_vs_direct_m_s": float(
                    np.linalg.norm(velocity_i - reference_velocity)
                ),
                "atlas_function_evaluations": int(solution.nfev),
                "gravity_tiles_loaded": len(atlas.cache),
            }
        )
    path = ROOT / "data" / "output" / "web_solver_crosscheck.json"
    path.write_text(json.dumps(results, indent=2), encoding="utf-8")
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
