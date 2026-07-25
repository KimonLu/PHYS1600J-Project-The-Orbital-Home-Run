#!/usr/bin/env python3
from __future__ import annotations
import math
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from orbital_home_run import C, minimum_speed_for_clearance, maximum_safe_angle, initial_state, osculating_elements_from_state, propagate, acceleration_central, body_fixed_lat_lon
from plotting import setup, save, OUT, ROOT


def load_ldem4(path: Path) -> np.ndarray:
    """Load NASA PDS LDEM4 (720x1440 PC_REAL, elevation in km)."""
    data = np.fromfile(path, dtype="<f4")
    expected = 720*1440
    if data.size != expected:
        raise ValueError(f"Expected {expected} float32 values, found {data.size}")
    return data.reshape(720, 1440)


def main() -> None:
    setup()
    hmax = C.highest_elevation
    launch_altitudes = np.array([1, 10, 100, 1000, 5000, hmax+1, 20000, 50000], dtype=float)
    rows = []
    for h in launch_altitudes:
        r0 = C.radius+h
        for obstacle_name, obstacle in [("mean sphere", C.radius), ("global terrain envelope", C.radius+hmax)]:
            if r0 <= obstacle:
                vmin = math.nan
                amax = math.nan
            else:
                vmin = minimum_speed_for_clearance(r0, 0.0, obstacle)
                vchosen = 1.05*math.sqrt(C.mu/r0)
                amax = maximum_safe_angle(r0, vchosen, obstacle)
            rows.append([h, obstacle_name, obstacle-C.radius, vmin, math.degrees(amax) if math.isfinite(amax) else math.nan])
    pd.DataFrame(rows, columns=["launch_altitude_m", "clearance_model", "obstacle_elevation_m", "minimum_horizontal_speed_m_s", "max_safe_angle_at_1p05_local_circular_deg"]).to_csv(OUT / "terrain_envelope_requirements.csv", index=False)

    heights = np.linspace(hmax+0.1, 60000, 500)
    vmins = np.array([minimum_speed_for_clearance(C.radius+h, 0.0, C.radius+hmax) for h in heights])
    fig, ax = plt.subplots(figsize=(3.45, 2.8))
    ax.plot(heights/1e3, vmins/1e3)
    ax.set_xlabel("Launch altitude above mean radius (km)")
    ax.set_ylabel("Minimum horizontal speed (km/s)")
    ax.set_title("Conservative global-terrain clearance")
    ax.axvline(hmax/1e3, linewidth=0.8, linestyle="--")
    save(fig, "fig11_terrain_envelope_speed")

    ldem = ROOT/"data"/"input"/"ldem_4_float.img"
    if ldem.exists():
        status = OUT/"ldem4_status.txt"
        if status.exists():
            status.unlink()
        z = load_ldem4(ldem)
        stats = pd.DataFrame([
            ["minimum_km", float(np.nanmin(z))],
            ["maximum_km", float(np.nanmax(z))],
            ["mean_km", float(np.nanmean(z))],
            ["std_km", float(np.nanstd(z))],
        ], columns=["quantity", "value"])
        stats.to_csv(OUT/"ldem4_computed_statistics.csv", index=False)
        # Downsample only for a reproducible overview; collision analysis uses the full grid.
        fig, ax = plt.subplots(figsize=(6.8, 3.0))
        im = ax.imshow(z[::2, ::2], extent=[0,360,-90,90], origin="upper", aspect="auto")
        ax.set_xlabel("Longitude (deg E)")
        ax.set_ylabel("Latitude (deg)")
        ax.set_title("Optional LOLA LDEM4 topography")
        fig.colorbar(im, ax=ax, label="Elevation (km)")
        save(fig, "fig_optional_ldem4_map")

        # Demonstration ground track for the same 20 km, u=1.03 case used in M3--M4.
        r0 = C.radius + 20_000.0
        v_inertial = 1.03*math.sqrt(C.mu/r0)
        v_surface = v_inertial-C.omega*r0
        state0 = initial_state(r0, 0.0, 0.0, v_surface, 0.0, math.radians(90), True)
        period = osculating_elements_from_state(state0)["period"]
        tr = propagate(state0, period, acceleration_central, max_step=5.0, samples=2501)
        gt = []
        for t, state in zip(tr["t"], tr["state"]):
            lat, lon = body_fixed_lat_lon(state[:3], float(t))
            lat_deg = math.degrees(lat)
            lon_deg = math.degrees(lon) % 360.0
            row = min(719, max(0, int((90.0-lat_deg)*4.0)))
            col = min(1439, max(0, int(lon_deg*4.0)))
            terrain_m = float(z[row, col])*1000.0
            radial_alt_m = float(np.linalg.norm(state[:3])-C.radius)
            clearance_m = radial_alt_m-terrain_m-C.ball_radius
            gt.append([t, lat_deg, lon_deg, radial_alt_m, terrain_m, clearance_m])
        gtdf = pd.DataFrame(gt, columns=["time_s","latitude_deg","longitude_deg_east","radial_altitude_m","terrain_elevation_m","ball_clearance_m"])
        gtdf.to_csv(OUT/"ldem4_groundtrack_clearance.csv", index=False)
        fig, ax = plt.subplots(figsize=(3.45,2.8))
        ax.plot(gtdf["time_s"]/60.0, gtdf["ball_clearance_m"]/1000.0)
        ax.axhline(0.0, linewidth=0.8)
        ax.set_xlabel("Time (min)")
        ax.set_ylabel("Terrain clearance (km)")
        ax.set_title("Optional LDEM4 ground-track clearance")
        save(fig,"fig_optional_groundtrack_clearance")
    else:
        (OUT/"ldem4_status.txt").write_text(
            "LDEM4 was not bundled. Run scripts/download_optional_data.py LOLA_LDEM4 LOLA_LDEM4_LABEL, then rerun this script.\n",
            encoding="utf-8"
        )


if __name__ == "__main__":
    main()
