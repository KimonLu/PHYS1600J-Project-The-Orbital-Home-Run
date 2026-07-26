#!/usr/bin/env python3
"""Quantify the terrain penalty for ideal near-surface great-circle orbits.

This is deliberately a terrain feasibility experiment, not a proof that a
high-fidelity return exists.  It samples the complete LDEM64 great circle at
1/64 degree and records the platform height required to keep a circular
central-gravity orbit above the bilinearly interpolated DEM.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from lunar_terrain import LDEM64
from plotting import save, setup


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "output"
FIGURES = ROOT / "figures"
RNG_SEED = 1600
N_RANDOM = 2_000
N_TRACK = 360 * 64
BALL_RADIUS_M = 0.0369
MAXIMUM_LAT_DEG = 5.4296875
MAXIMUM_LON_DEG = 201.3671875


def unit_vector(latitude_deg: float, longitude_deg: float) -> np.ndarray:
    lat = np.radians(latitude_deg)
    lon = np.radians(longitude_deg)
    return np.array([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)])


def tangent(latitude_deg: float, longitude_deg: float, azimuth_deg: float) -> np.ndarray:
    lat = np.radians(latitude_deg)
    lon = np.radians(longitude_deg)
    az = np.radians(azimuth_deg)
    north = np.array([-np.sin(lat) * np.cos(lon), -np.sin(lat) * np.sin(lon), np.cos(lat)])
    east = np.array([-np.sin(lon), np.cos(lon), 0.0])
    return np.cos(az) * north + np.sin(az) * east


def great_circle_requirement(
    terrain: LDEM64, latitude_deg: float, longitude_deg: float, azimuth_deg: float
) -> tuple[float, float, float]:
    r0 = unit_vector(latitude_deg, longitude_deg)
    t0 = tangent(latitude_deg, longitude_deg, azimuth_deg)
    phase = np.linspace(0.0, 2.0 * np.pi, N_TRACK, endpoint=False)
    points = np.cos(phase)[:, None] * r0 + np.sin(phase)[:, None] * t0
    latitude = np.degrees(np.arcsin(np.clip(points[:, 2], -1.0, 1.0)))
    longitude = np.degrees(np.arctan2(points[:, 1], points[:, 0])) % 360.0
    profile = np.asarray(terrain.elevation_m(latitude, longitude))
    launch_elevation = float(terrain.elevation_m(latitude_deg, longitude_deg))
    required_height = float(np.max(profile) - launch_elevation + BALL_RADIUS_M)
    return max(required_height, BALL_RADIUS_M), float(np.max(profile)), launch_elevation


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    terrain = LDEM64()
    rng = np.random.default_rng(RNG_SEED)
    rows: list[dict[str, float | str]] = []

    z = rng.uniform(-1.0, 1.0, N_RANDOM)
    latitudes = np.degrees(np.arcsin(z))
    longitudes = rng.uniform(0.0, 360.0, N_RANDOM)
    azimuths = rng.uniform(0.0, 360.0, N_RANDOM)
    for index, (lat, lon, az) in enumerate(zip(latitudes, longitudes, azimuths)):
        required, maximum, launch = great_circle_requirement(terrain, lat, lon, az)
        rows.append(
            {
                "sample": index,
                "family": "uniform_random_site_and_azimuth",
                "latitude_deg": lat,
                "longitude_deg_east": lon,
                "azimuth_deg": az,
                "launch_elevation_m": launch,
                "maximum_track_elevation_m": maximum,
                "minimum_platform_height_m": required,
            }
        )

    for azimuth in np.arange(0.0, 360.0, 1.0):
        required, maximum, launch = great_circle_requirement(
            terrain, MAXIMUM_LAT_DEG, MAXIMUM_LON_DEG, float(azimuth)
        )
        rows.append(
            {
                "sample": int(azimuth),
                "family": "ldem64_maximum_site_azimuth_scan",
                "latitude_deg": MAXIMUM_LAT_DEG,
                "longitude_deg_east": MAXIMUM_LON_DEG,
                "azimuth_deg": float(azimuth),
                "launch_elevation_m": launch,
                "maximum_track_elevation_m": maximum,
                "minimum_platform_height_m": required,
            }
        )

    frame = pd.DataFrame(rows)
    frame.to_csv(OUT / "near_surface_great_circle_samples.csv", index=False)
    random_values = frame.loc[
        frame["family"] == "uniform_random_site_and_azimuth", "minimum_platform_height_m"
    ].to_numpy()
    summit_values = frame.loc[
        frame["family"] == "ldem64_maximum_site_azimuth_scan", "minimum_platform_height_m"
    ].to_numpy()
    summary = {
        "interpretation": (
            "Central-gravity circular great-circle terrain requirement relative to the "
            "bilinearly interpolated LDEM64 surface; not a high-fidelity return search."
        ),
        "random_seed": RNG_SEED,
        "random_sample_count": N_RANDOM,
        "along_track_spacing_deg": 1.0 / 64.0,
        "random_required_height_m": {
            "minimum": float(np.min(random_values)),
            "median": float(np.median(random_values)),
            "p95": float(np.quantile(random_values, 0.95)),
            "maximum": float(np.max(random_values)),
            "fraction_at_or_below_1m": float(np.mean(random_values <= 1.0)),
            "fraction_at_or_below_100m": float(np.mean(random_values <= 100.0)),
        },
        "ldem64_maximum_site_scan": {
            "azimuth_count": 360,
            "minimum_required_height_m": float(np.min(summit_values)),
            "maximum_required_height_m": float(np.max(summit_values)),
        },
        "high_fidelity_height_continuation": {
            "status": "performed_separately_by_16_height_continuation.py",
            "output": "near_surface_height_continuation.csv",
            "scope": "2 km and above at the LDEM64 summit site",
        },
    }
    (OUT / "near_surface_great_circle_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )

    setup()
    fig, ax = plt.subplots(figsize=(8.4, 4.8))
    bins = np.linspace(0.0, max(20_000.0, float(np.max(random_values))), 55)
    ax.hist(random_values / 1000.0, bins=bins / 1000.0, color="#4477aa", alpha=0.82)
    ax.axvline(np.min(summit_values) / 1000.0, color="#cc6677", linewidth=2.0,
               label="Best azimuth from LDEM64 maximum")
    ax.set_xlabel("Required platform height above local LDEM64 terrain (km)")
    ax.set_ylabel("Random great-circle count")
    ax.set_title("Near-surface circular-orbit terrain requirement")
    ax.legend()
    save(fig, "fig22_surface_height_feasibility")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
