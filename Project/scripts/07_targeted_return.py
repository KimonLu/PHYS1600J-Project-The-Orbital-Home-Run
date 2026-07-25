#!/usr/bin/env python3
"""Terrain-aware, rotating-site return and differential correction.

The calculation deliberately answers a question that a one-inertial-period
closure test cannot answer: when does the ball meet the rotating ballpark?
For an equatorial prograde circular orbit the analytic central-gravity answer
is the synodic period 2*pi/(n-Omega).  The same case is then corrected in the
degree-2 plus fixed third-body model.
"""
from __future__ import annotations

import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Circle
from scipy.optimize import least_squares

from orbital_home_run import (
    C,
    body_fixed_lat_lon,
    initial_state,
    make_acceleration_model,
    propagate,
    rot_z,
)
from plotting import OUT, ROOT, save, setup


LDEM_ROWS = 720
LDEM_COLS = 1440
RETURN_TOLERANCE_M = 0.5
TERRAIN_MARGIN_M = 100.0


def load_ldem4(path: Path) -> np.ndarray:
    data = np.fromfile(path, dtype="<f4")
    if data.size != LDEM_ROWS * LDEM_COLS:
        raise ValueError(
            f"Expected {LDEM_ROWS * LDEM_COLS} LDEM4 samples, found {data.size}"
        )
    return data.reshape(LDEM_ROWS, LDEM_COLS)


def final_state(
    launch_radius: float,
    longitude: float,
    speed_surface: float,
    gamma: float,
    azimuth: float,
    duration: float,
    acceleration,
    max_step: float = 3.0,
) -> np.ndarray:
    state0 = initial_state(
        launch_radius,
        0.0,
        longitude,
        speed_surface,
        gamma,
        azimuth,
        True,
    )
    return propagate(
        state0,
        duration,
        acceleration,
        max_step=max_step,
        rtol=2e-11,
        atol=2e-6,
        samples=2,
    )["state"][-1]


def solve_scheduled_return(
    launch_radius: float,
    longitude: float,
    duration: float,
    speed_guess: float,
    acceleration,
) -> tuple[float, float]:
    """Solve eastward speed and elevation for a fixed-time equatorial return."""
    site_body = np.array(
        [
            launch_radius * math.cos(longitude),
            launch_radius * math.sin(longitude),
            0.0,
        ]
    )
    target = rot_z(C.omega * duration) @ site_body
    azimuth = 0.5 * math.pi

    def residual(x: np.ndarray) -> np.ndarray:
        final = final_state(
            launch_radius,
            longitude,
            float(x[0]),
            float(x[1]),
            azimuth,
            duration,
            acceleration,
            max_step=5.0,
        )
        # The selected force geometry is symmetric about the equatorial plane,
        # so z remains zero.  Scale planar residuals to kilometres.
        return (final[:2] - target[:2]) / 1000.0

    solution = least_squares(
        residual,
        np.array([speed_guess, 0.0]),
        x_scale=np.array([1.0, 1e-4]),
        diff_step=np.array([1e-5, 1e-4]),
        xtol=1e-11,
        ftol=1e-11,
        gtol=1e-11,
        max_nfev=60,
    )
    if not solution.success or np.linalg.norm(solution.fun) > 1e-6:
        raise RuntimeError(
            "Scheduled-return differential correction did not converge: "
            f"{solution.message}; residual={solution.fun}"
        )
    return float(solution.x[0]), float(solution.x[1])


def terrain_clearance(
    trajectory: dict[str, np.ndarray],
    equatorial_profile_m: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Nearest-cell LDEM4 clearance for an equatorial trajectory."""
    relative_longitude = []
    terrain = []
    altitude = []
    for t, state in zip(trajectory["t"], trajectory["state"]):
        latitude, longitude = body_fixed_lat_lon(state[:3], float(t))
        if abs(latitude) > 1e-9:
            raise ValueError("Targeted-return case left the equatorial symmetry plane")
        longitude_deg = math.degrees(longitude) % 360.0
        col = min(LDEM_COLS - 1, max(0, int(longitude_deg * 4.0)))
        relative_longitude.append(longitude_deg)
        terrain.append(float(equatorial_profile_m[col]))
        altitude.append(float(np.linalg.norm(state[:3]) - C.radius))
    altitude_a = np.asarray(altitude)
    terrain_a = np.asarray(terrain)
    clearance = altitude_a - terrain_a - C.ball_radius
    return np.asarray(relative_longitude), altitude_a, clearance


def projected_sensitivity(
    launch_radius: float,
    longitude: float,
    duration: float,
    parameters: np.ndarray,
    acceleration,
) -> pd.DataFrame:
    """Linearized closest-approach sensitivity with encounter time free.

    Fixed-time terminal derivatives are projected perpendicular to the
    ball-site relative velocity.  The removed component changes encounter
    time to first order rather than the spatial closest approach.
    """
    site_body = np.array(
        [
            launch_radius * math.cos(longitude),
            launch_radius * math.sin(longitude),
            0.0,
        ]
    )
    site = rot_z(C.omega * duration) @ site_body
    nominal = final_state(
        launch_radius,
        longitude,
        parameters[0],
        parameters[1],
        parameters[2],
        duration,
        acceleration,
        max_step=2.0,
    )
    relative_velocity = nominal[3:] - np.cross(
        np.array([0.0, 0.0, C.omega]), site
    )
    projector = np.eye(3) - np.outer(relative_velocity, relative_velocity) / (
        relative_velocity @ relative_velocity
    )

    steps = np.array([0.01, 1e-5, 1e-5])
    names = ["surface_speed", "elevation_angle", "azimuth"]
    units = ["m/s", "rad", "rad"]
    derivative = np.empty((3, 3))
    for k, step in enumerate(steps):
        plus = parameters.copy()
        minus = parameters.copy()
        plus[k] += step
        minus[k] -= step
        final_plus = final_state(
            launch_radius,
            longitude,
            plus[0],
            plus[1],
            plus[2],
            duration,
            acceleration,
            max_step=2.0,
        )
        final_minus = final_state(
            launch_radius,
            longitude,
            minus[0],
            minus[1],
            minus[2],
            duration,
            acceleration,
            max_step=2.0,
        )
        derivative[:, k] = (final_plus[:3] - final_minus[:3]) / (2.0 * step)

    projected = projector @ derivative
    gains = np.linalg.norm(projected, axis=0)
    tolerances = RETURN_TOLERANCE_M / gains
    return pd.DataFrame(
        {
            "parameter": names,
            "parameter_unit": units,
            "finite_difference_step": steps,
            "closest_miss_gain_m_per_unit": gains,
            "parameter_tolerance_for_0p5m_miss": tolerances,
            "tolerance_deg_if_angle": [
                math.nan,
                math.degrees(tolerances[1]),
                math.degrees(tolerances[2]),
            ],
        }
    )


def main() -> None:
    setup()
    ldem_path = ROOT / "data" / "input" / "ldem_4_float.img"
    if not ldem_path.exists():
        raise FileNotFoundError(
            "The targeted terrain-aware result requires NASA LOLA LDEM4. "
            "Run download_optional_data.py LOLA_LDEM4 first."
        )
    ldem = load_ldem4(ldem_path)

    # LDEM4 is pixel-registered at +/-0.125 deg around the equator.  Taking
    # the maximum of the two bracketing rows is a conservative equatorial
    # terrain profile at the product's 0.25-deg resolution.
    equatorial_profile_m = (
        np.maximum(ldem[359, :], ldem[360, :]).astype(float) * 1000.0
    )
    peak_col = int(np.argmax(equatorial_profile_m))
    peak_longitude_deg = (peak_col + 0.5) / 4.0
    peak_longitude = math.radians(peak_longitude_deg)
    equatorial_peak_m = float(equatorial_profile_m[peak_col])

    launch_altitude = equatorial_peak_m + TERRAIN_MARGIN_M + C.ball_radius
    launch_radius = C.radius + launch_altitude
    mean_motion = math.sqrt(C.mu / launch_radius**3)
    inertial_circular_speed = math.sqrt(C.mu / launch_radius)
    surface_circular_speed = inertial_circular_speed - C.omega * launch_radius
    inertial_period = 2.0 * math.pi / mean_motion
    synodic_period = 2.0 * math.pi / (mean_motion - C.omega)
    azimuth = 0.5 * math.pi

    central = make_acceleration_model()
    degree2 = make_acceleration_model(include_degree2=True)
    full = make_acceleration_model(
        include_degree2=True, include_earth=True, include_sun=True
    )

    speed_d2, gamma_d2 = solve_scheduled_return(
        launch_radius,
        peak_longitude,
        synodic_period,
        surface_circular_speed,
        degree2,
    )
    speed_full, gamma_full = solve_scheduled_return(
        launch_radius,
        peak_longitude,
        synodic_period,
        surface_circular_speed,
        full,
    )

    cases = [
        ("Central analytic", central, surface_circular_speed, 0.0),
        ("Full, uncorrected", full, surface_circular_speed, 0.0),
        ("Degree-2, targeted", degree2, speed_d2, gamma_d2),
        ("Full, targeted", full, speed_full, gamma_full),
    ]
    site_body = np.array(
        [
            launch_radius * math.cos(peak_longitude),
            launch_radius * math.sin(peak_longitude),
            0.0,
        ]
    )
    target = rot_z(C.omega * synodic_period) @ site_body
    trajectories: dict[str, dict[str, np.ndarray]] = {}
    summary_rows = []
    for name, acceleration, speed, gamma in cases:
        state0 = initial_state(
            launch_radius,
            0.0,
            peak_longitude,
            speed,
            gamma,
            azimuth,
            True,
        )
        trajectory = propagate(
            state0,
            synodic_period,
            acceleration,
            max_step=2.0,
            samples=4001,
        )
        trajectories[name] = trajectory
        longitude, altitude, clearance = terrain_clearance(
            trajectory, equatorial_profile_m
        )
        theta = np.unwrap(
            np.arctan2(
                trajectory["state"][:, 1],
                trajectory["state"][:, 0],
            )
        )
        inertial_revolutions = (theta[-1] - theta[0]) / (2.0 * math.pi)
        relative_revolutions = (
            theta[-1] - theta[0] - C.omega * synodic_period
        ) / (2.0 * math.pi)
        final_miss = float(
            np.linalg.norm(trajectory["state"][-1, :3] - target)
        )
        site_velocity = np.cross(
            np.array([0.0, 0.0, C.omega]),
            target,
        )
        terminal_relative_speed = float(
            np.linalg.norm(trajectory["state"][-1, 3:] - site_velocity)
        )
        summary_rows.append(
            [
                name,
                speed,
                math.degrees(gamma),
                synodic_period,
                final_miss,
                float(np.min(altitude)),
                float(np.max(altitude)),
                float(np.min(clearance)),
                inertial_revolutions,
                relative_revolutions,
                terminal_relative_speed,
                0.5 * C.ball_mass * terminal_relative_speed**2,
            ]
        )

    pd.DataFrame(
        summary_rows,
        columns=[
            "model",
            "surface_relative_speed_m_s",
            "elevation_angle_deg",
            "scheduled_return_time_s",
            "scheduled_site_miss_m",
            "minimum_altitude_above_mean_sphere_m",
            "maximum_altitude_above_mean_sphere_m",
            "minimum_ldem4_clearance_m",
            "inertial_revolutions",
            "body_fixed_revolutions",
            "terminal_relative_speed_m_s",
            "terminal_kinetic_energy_J",
        ],
    ).to_csv(OUT / "targeted_return_summary.csv", index=False)

    pd.DataFrame(
        [
            ["LDEM4_equatorial_peak_m", equatorial_peak_m],
            ["stadium_longitude_deg_east", peak_longitude_deg],
            ["terrain_clearance_margin_m", TERRAIN_MARGIN_M],
            ["ball_center_launch_altitude_m", launch_altitude],
            ["central_inertial_circular_speed_m_s", inertial_circular_speed],
            ["central_surface_relative_speed_m_s", surface_circular_speed],
            ["central_inertial_period_s", inertial_period],
            ["central_synodic_return_time_s", synodic_period],
            ["synodic_delay_beyond_inertial_period_s", synodic_period - inertial_period],
        ],
        columns=["quantity", "value"],
    ).to_csv(OUT / "targeted_return_case_definition.csv", index=False)

    sensitivity = projected_sensitivity(
        launch_radius,
        peak_longitude,
        synodic_period,
        np.array([speed_full, gamma_full, azimuth]),
        full,
    )
    sensitivity.to_csv(OUT / "targeted_return_sensitivity.csv", index=False)

    # True-scale trajectory plus the body-fixed longitude history makes both
    # the global orbit and the moving-site closure visually explicit.
    fig, axes = plt.subplots(1, 2, figsize=(7.1, 3.0))
    ax = axes[0]
    ax.add_patch(
        Circle(
            (0.0, 0.0),
            C.radius / 1e3,
            facecolor="0.90",
            edgecolor="0.25",
            linewidth=0.8,
            zorder=0,
        )
    )
    for name in ["Central analytic", "Full, targeted"]:
        state = trajectories[name]["state"]
        ax.plot(state[:, 0] / 1e3, state[:, 1] / 1e3, label=name)
    ax.scatter(
        [site_body[0] / 1e3],
        [site_body[1] / 1e3],
        marker="*",
        s=45,
        zorder=5,
        label="stadium at launch",
    )
    ax.set_aspect("equal")
    lim = (launch_radius + 45_000.0) / 1e3
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_xlabel("Moon-centred $x$ (km)")
    ax.set_ylabel("Moon-centred $y$ (km)")
    ax.set_title("A complete prograde lunar circuit")
    ax.legend(fontsize=6.2, loc="lower left")

    ax = axes[1]
    for name in ["Central analytic", "Full, uncorrected", "Full, targeted"]:
        trajectory = trajectories[name]
        theta = np.unwrap(
            np.arctan2(
                trajectory["state"][:, 1],
                trajectory["state"][:, 0],
            )
        )
        relative = np.degrees(
            theta - theta[0] - C.omega * trajectory["t"]
        )
        ax.plot(trajectory["t"] / 60.0, relative, label=name)
    ax.axhline(360.0, linewidth=0.8, linestyle="--", color="0.3")
    ax.set_xlabel("Time after launch (min)")
    ax.set_ylabel("Longitude travelled relative to field (deg)")
    ax.set_title("Return occurs after one relative revolution")
    ax.legend(fontsize=6.2)
    save(fig, "fig13_targeted_orbit_overview")

    fig, ax = plt.subplots(figsize=(7.1, 3.0))
    longitude_grid = (np.arange(LDEM_COLS) + 0.5) / 4.0
    ax.fill_between(
        longitude_grid,
        equatorial_profile_m / 1e3,
        -10.0,
        color="0.82",
        label="LOLA LDEM4 equatorial terrain",
    )
    for name in ["Central analytic", "Full, uncorrected", "Full, targeted"]:
        trajectory = trajectories[name]
        longitude, altitude, _ = terrain_clearance(
            trajectory, equatorial_profile_m
        )
        order = np.argsort(longitude)
        ax.plot(longitude[order], altitude[order] / 1e3, label=name)
    ax.axvline(
        peak_longitude_deg,
        linewidth=0.8,
        linestyle="--",
        color="0.3",
    )
    ax.set_xlim(0.0, 360.0)
    ax.set_ylim(float(np.min(equatorial_profile_m) / 1e3 - 0.5), launch_altitude / 1e3 + 1.5)
    ax.set_xlabel("Body-fixed longitude (deg E)")
    ax.set_ylabel("Elevation above mean radius (km)")
    ax.set_title("Terrain-aware equatorial orbital corridor")
    ax.legend(fontsize=6.2, ncol=2)
    save(fig, "fig14_targeted_clearance")

    fig, ax = plt.subplots(figsize=(3.45, 2.8))
    gains = sensitivity["closest_miss_gain_m_per_unit"].to_numpy()
    angle_gains = math.radians(1.0) * gains[1:]
    ax.bar(
        ["elevation", "azimuth"],
        angle_gains,
        color=["C1", "C2"],
    )
    ax.set_ylabel("Closest miss from a 1 deg error (m)")
    ax.set_title("Return targeting is angle-limited")
    save(fig, "fig15_targeted_sensitivity")


if __name__ == "__main__":
    main()
