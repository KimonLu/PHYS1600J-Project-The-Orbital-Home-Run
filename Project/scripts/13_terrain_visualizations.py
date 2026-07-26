#!/usr/bin/env python3
"""Generate LDEM64 topography and high-fidelity corridor figures."""
from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from lunar_terrain import LDEM64, META
from plotting import OUT, save, setup


ROOT = Path(__file__).resolve().parents[1]
TRAJECTORY = OUT / "high_fidelity_case_trajectory.csv"
SUMMARY = OUT / "high_fidelity_case_summary.json"
BLOCK_SIZE = 16


def ldem64_block_mean_km(terrain: LDEM64) -> np.ndarray:
    """Return a display-only 16 x 16 block mean of the full LDEM64 grid."""
    if META.rows % BLOCK_SIZE or META.columns % BLOCK_SIZE:
        raise ValueError("LDEM64 dimensions are not divisible by the block size")
    output_rows = META.rows // BLOCK_SIZE
    output_columns = META.columns // BLOCK_SIZE
    result = np.empty((output_rows, output_columns), dtype=np.float32)
    for output_row in range(output_rows):
        start = output_row * BLOCK_SIZE
        block = np.asarray(
            terrain.grid[start : start + BLOCK_SIZE], dtype=np.float32
        )
        result[output_row] = block.reshape(
            BLOCK_SIZE, output_columns, BLOCK_SIZE
        ).mean(axis=(0, 2))
    result *= META.scale_m_per_dn / 1000.0
    return result


def make_global_topography(terrain: LDEM64) -> None:
    elevation_km = ldem64_block_mean_km(terrain)
    exact_minimum_m, exact_maximum_m = terrain.extrema()
    longitude_centres_deg = (
        np.arange(elevation_km.shape[1]) + 0.5
    ) * BLOCK_SIZE / META.pixels_per_degree
    latitude_centres_deg = 90.0 - (
        np.arange(elevation_km.shape[0]) + 0.5
    ) * BLOCK_SIZE / META.pixels_per_degree

    np.savez_compressed(
        OUT / "ldem64_global_heatmap_16x16.npz",
        elevation_km=elevation_km,
        longitude_centres_deg=longitude_centres_deg,
        latitude_centres_deg=latitude_centres_deg,
    )
    metadata = {
        "source_product": "LOLA LDEM64",
        "source_grid_rows": META.rows,
        "source_grid_columns": META.columns,
        "source_pixels_per_degree": META.pixels_per_degree,
        "source_scale_m_per_dn": META.scale_m_per_dn,
        "reference_radius_m": META.reference_radius_m,
        "display_aggregation": "16 x 16 arithmetic block mean",
        "display_grid_rows": int(elevation_km.shape[0]),
        "display_grid_columns": int(elevation_km.shape[1]),
        "exact_source_minimum_m": exact_minimum_m,
        "exact_source_maximum_m": exact_maximum_m,
        "display_minimum_km": float(elevation_km.min()),
        "display_maximum_km": float(elevation_km.max()),
        "collision_queries_use_full_resolution": True,
    }
    (OUT / "ldem64_global_heatmap_summary.json").write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
    )

    fig, ax = plt.subplots(figsize=(7.2, 3.75), constrained_layout=True)
    image = ax.imshow(
        elevation_km,
        extent=(0.0, 360.0, -90.0, 90.0),
        origin="upper",
        aspect="auto",
        cmap="viridis",
        vmin=exact_minimum_m / 1000.0,
        vmax=exact_maximum_m / 1000.0,
        interpolation="nearest",
        rasterized=True,
    )
    colourbar = fig.colorbar(image, ax=ax, pad=0.025)
    colourbar.set_label("Elevation above 1737.4 km sphere (km)")
    ax.set(
        title="Global LOLA LDEM64 topography",
        xlabel="Longitude (deg E)",
        ylabel="Latitude (deg)",
        xlim=(0.0, 360.0),
        ylim=(-90.0, 90.0),
    )
    ax.set_xticks(np.arange(0.0, 361.0, 60.0))
    ax.set_yticks(np.arange(-90.0, 91.0, 30.0))
    save(fig, "fig20_ldem64_topography")


def load_trajectory() -> dict[str, np.ndarray]:
    if not TRAJECTORY.exists():
        raise FileNotFoundError(
            f"{TRAJECTORY} is missing; run 10_high_fidelity_case.py first"
        )
    data = np.genfromtxt(TRAJECTORY, delimiter=",", names=True)
    return {name: np.asarray(data[name], dtype=float) for name in data.dtype.names}


def make_corridor_profile() -> None:
    trajectory = load_trajectory()
    summary = json.loads(SUMMARY.read_text(encoding="utf-8"))
    longitude = trajectory["longitude_deg_east"]
    unwrapped_longitude = np.degrees(np.unwrap(np.radians(longitude)))
    signed_progress = unwrapped_longitude - unwrapped_longitude[0]
    direction = 1.0 if signed_progress[-1] >= 0.0 else -1.0
    track_longitude = direction * signed_progress

    altitude_km = trajectory["altitude_above_reference_m"] / 1000.0
    terrain_km = trajectory["terrain_elevation_m"] / 1000.0
    clearance_km = trajectory["clearance_m"] / 1000.0
    minimum_index = int(np.argmin(clearance_km))

    profile_path = OUT / "high_fidelity_corridor_profile.csv"
    with profile_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "time_s",
                "track_longitude_from_launch_deg",
                "latitude_deg",
                "longitude_deg_east",
                "ball_center_elevation_above_reference_km",
                "terrain_elevation_above_reference_km",
                "clearance_km",
            ]
        )
        writer.writerows(
            zip(
                trajectory["time_s"],
                track_longitude,
                trajectory["latitude_deg"],
                longitude,
                altitude_km,
                terrain_km,
                clearance_km,
            )
        )

    lower_limit = min(-12.0, float(terrain_km.min()) - 1.0)
    upper_limit = float(altitude_km.max()) + 2.0
    fig, ax = plt.subplots(figsize=(7.2, 3.65), constrained_layout=True)
    ax.fill_between(
        track_longitude,
        lower_limit,
        terrain_km,
        color="0.82",
        label="LOLA LDEM64 terrain along ground track",
        zorder=1,
    )
    ax.plot(
        track_longitude,
        terrain_km,
        color="0.48",
        linewidth=0.8,
        zorder=2,
    )
    ax.plot(
        track_longitude,
        altitude_km,
        color="#1769aa",
        linewidth=1.8,
        label="Degree-600 ball-centre trajectory",
        zorder=3,
    )
    x_minimum = track_longitude[minimum_index]
    ax.plot(
        [x_minimum, x_minimum],
        [terrain_km[minimum_index], altitude_km[minimum_index]],
        color="#d95f02",
        linewidth=2.0,
        zorder=4,
    )
    ax.scatter(
        [x_minimum, x_minimum],
        [terrain_km[minimum_index], altitude_km[minimum_index]],
        color="#d95f02",
        s=16,
        zorder=5,
    )
    ax.annotate(
        f"minimum clearance = {clearance_km[minimum_index]:.3f} km",
        xy=(x_minimum, 0.5 * (terrain_km[minimum_index] + altitude_km[minimum_index])),
        xytext=(-8, 0),
        textcoords="offset points",
        va="center",
        ha="right",
        fontsize=8,
        color="#9c3d00",
    )
    ax.set(
        title="Terrain clearance along the degree-600 return corridor",
        xlabel="Body-fixed longitude travelled from launch (deg)",
        ylabel="Elevation above 1737.4 km sphere (km)",
        xlim=(float(track_longitude.min()), float(track_longitude.max())),
        ylim=(lower_limit, upper_limit),
    )
    ax.legend(loc="upper right")
    ax.text(
        0.01,
        0.02,
        (
            f"Ground track: {summary['latitude_deg']:.4f} deg launch latitude; "
            "latitude varies along the 3-D arc"
        ),
        transform=ax.transAxes,
        fontsize=7.5,
        color="0.30",
    )
    save(fig, "fig19_high_fidelity_corridor")


def main() -> None:
    setup()
    terrain = LDEM64()
    make_corridor_profile()
    make_global_topography(terrain)
    print("Generated fig19_high_fidelity_corridor and fig20_ldem64_topography")


if __name__ == "__main__":
    main()
