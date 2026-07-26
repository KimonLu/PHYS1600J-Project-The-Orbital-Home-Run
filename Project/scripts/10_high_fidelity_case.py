#!/usr/bin/env python3
"""Differentially correct one terrain-safe, degree-600 scheduled return."""
from __future__ import annotations

import json
import math
from pathlib import Path
import time

import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares

from lunar_gravity import BRILLOUIN_RADIUS_M, GRAILGravity
from lunar_terrain import LDEM64, META as TERRAIN_META
from orbital_home_run import (
    C,
    EARTH_MOON_DISTANCE,
    MU_EARTH,
    MU_SUN,
    SUN_MOON_DISTANCE,
    acceleration_third_body,
    local_basis,
    rot_z,
)
from plotting import save, setup


ROOT = Path(__file__).resolve().parents[1]
DATA_OUT = ROOT / "data" / "output"
SCENARIO_LABEL = "fixed_short_arc_third_body_geometry"
LATITUDE_DEG = 5.4296875
LONGITUDE_DEG = 201.3671875
REFERENCE_ALTITUDE_M = 30_000.0
POSITION_RESIDUAL_TARGET_M = 0.01
JACOBIAN_COMPONENT_STEP_M_S = 0.02


def initial_state(velocity_local_neu: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    latitude = math.radians(LATITUDE_DEG)
    longitude = math.radians(LONGITUDE_DEG)
    radial, east, north = local_basis(latitude, longitude)
    position = (C.radius + REFERENCE_ALTITUDE_M) * radial
    velocity_relative = (
        velocity_local_neu[0] * north
        + velocity_local_neu[1] * east
        + velocity_local_neu[2] * radial
    )
    velocity_inertial = velocity_relative + np.cross(
        np.array([0.0, 0.0, C.omega]), position
    )
    return np.concatenate([position, velocity_inertial]), position


def acceleration(gravity: GRAILGravity):
    earth = np.array([EARTH_MOON_DISTANCE, 0.0, 0.0])
    sun = np.array([0.0, SUN_MOON_DISTANCE, 0.0])

    def model(time_s: float, position: np.ndarray) -> np.ndarray:
        return (
            gravity.inertial_acceleration(time_s, position)
            + acceleration_third_body(position, earth, MU_EARTH)
            + acceleration_third_body(position, sun, MU_SUN)
        )

    return model


def propagate(
    velocity_local_neu: np.ndarray,
    duration_s: float,
    model,
    samples: int = 2,
    *,
    method: str = "DOP853",
    max_step_s: float = 10.0,
    rtol: float = 2e-10,
    atol: np.ndarray | None = None,
    dense_output: bool = False,
) -> dict[str, np.ndarray]:
    state0, _ = initial_state(velocity_local_neu)

    def rhs(time_s: float, state: np.ndarray) -> np.ndarray:
        return np.concatenate([state[3:], model(time_s, state[:3])])

    if atol is None:
        atol = np.array([2e-4, 2e-4, 2e-4, 2e-7, 2e-7, 2e-7])
    solution = solve_ivp(
        rhs,
        (0.0, duration_s),
        state0,
        method=method,
        t_eval=np.linspace(0.0, duration_s, samples),
        max_step=max_step_s,
        rtol=rtol,
        atol=atol,
        dense_output=dense_output,
    )
    if not solution.success:
        raise RuntimeError(solution.message)
    result = {"time_s": solution.t, "state": solution.y.T, "nfev": solution.nfev}
    if dense_output:
        result["dense_solution"] = solution.sol
    return result


def terminal_position_jacobian(
    velocity_local_neu: np.ndarray,
    duration_s: float,
    model,
    step_m_s: float = JACOBIAN_COMPONENT_STEP_M_S,
) -> np.ndarray:
    """Central finite-difference terminal-position Jacobian in NEU components."""
    jacobian = np.empty((3, 3))
    for column in range(3):
        plus = velocity_local_neu.copy()
        minus = velocity_local_neu.copy()
        plus[column] += step_m_s
        minus[column] -= step_m_s
        plus_final = propagate(plus, duration_s, model, samples=2)["state"][-1, :3]
        minus_final = propagate(minus, duration_s, model, samples=2)["state"][-1, :3]
        jacobian[:, column] = (plus_final - minus_final) / (2.0 * step_m_s)
    return jacobian


def main() -> None:
    terrain = LDEM64()
    terrain_launch = terrain.elevation_m(LATITUDE_DEG, LONGITUDE_DEG)
    height_above_terrain = REFERENCE_ALTITUDE_M - terrain_launch
    radius = C.radius + REFERENCE_ALTITUDE_M
    mean_motion = math.sqrt(C.mu / radius**3)
    duration = 2.0 * math.pi / (mean_motion - C.omega)
    circular_surface_speed = math.sqrt(C.mu / radius) - C.omega * radius
    _, site_body = initial_state(np.array([0.0, circular_surface_speed, 0.0]))
    target_inertial = rot_z(C.omega * duration) @ site_body

    print("Loading degree-300 correction model", flush=True)
    gravity300 = GRAILGravity(maximum_degree=300)
    model300 = acceleration(gravity300)

    def residual300(velocity: np.ndarray) -> np.ndarray:
        final = propagate(velocity, duration, model300, samples=2)["state"][-1]
        value = (final[:3] - target_inertial) / 1000.0
        print(
            "trial",
            np.array2string(velocity, precision=7),
            "miss_km",
            f"{np.linalg.norm(value):.6g}",
            flush=True,
        )
        return value

    guess = np.array([0.0, circular_surface_speed, 0.0])
    correction = least_squares(
        residual300,
        guess,
        x_scale=np.ones(3),
        diff_step=2e-5,
        xtol=2e-10,
        ftol=2e-10,
        gtol=2e-10,
        max_nfev=30,
    )
    if not correction.success or np.linalg.norm(correction.fun) > 1e-5:
        raise RuntimeError(
            f"Degree-300 correction failed: {correction.message}; {correction.fun}"
        )
    velocity = correction.x.copy()
    degree300_residual_m = 1000.0 * np.asarray(correction.fun, dtype=float)

    print("Loading selected degree-600 model", flush=True)
    gravity600 = GRAILGravity(maximum_degree=600)
    model600 = acceleration(gravity600)
    started = time.perf_counter()
    correction_history: list[dict[str, object]] = []
    initial_degree600_miss = math.nan
    for iteration in range(4):
        direct = propagate(velocity, duration, model600, samples=2)
        residual600 = direct["state"][-1, :3] - target_inertial
        residual_norm = float(np.linalg.norm(residual600))
        if iteration == 0:
            initial_degree600_miss = residual_norm
        print(
            f"Degree-600 iteration {iteration}: residual={residual_norm:.9f} m",
            flush=True,
        )
        entry: dict[str, object] = {
            "iteration": iteration,
            "velocity_neu_m_s_before": velocity.tolist(),
            "residual_vector_m_before": residual600.tolist(),
            "residual_norm_m_before": residual_norm,
        }
        if residual_norm <= POSITION_RESIDUAL_TARGET_M:
            entry["converged_without_update"] = True
            correction_history.append(entry)
            break
        jacobian600 = terminal_position_jacobian(
            velocity,
            duration,
            model600,
            step_m_s=JACOBIAN_COMPONENT_STEP_M_S,
        )
        velocity_update = np.linalg.solve(jacobian600, -residual600)
        entry.update(
            {
                "jacobian_degree": 600,
                "jacobian_component_step_m_s": JACOBIAN_COMPONENT_STEP_M_S,
                "terminal_position_jacobian_s": jacobian600.tolist(),
                "velocity_update_neu_m_s": velocity_update.tolist(),
                "velocity_update_norm_m_s": float(np.linalg.norm(velocity_update)),
            }
        )
        correction_history.append(entry)
        velocity += velocity_update
    else:
        raise RuntimeError("Degree-600 Newton correction exceeded four iterations")

    final_trajectory = propagate(velocity, duration, model600, samples=2001)
    final_state = final_trajectory["state"][-1]
    final_residual_vector = final_state[:3] - target_inertial
    final_miss = float(np.linalg.norm(final_residual_vector))
    if final_miss > POSITION_RESIDUAL_TARGET_M:
        raise RuntimeError(
            f"Degree-600 residual {final_miss:.6g} m exceeds "
            f"{POSITION_RESIDUAL_TARGET_M:g} m target"
        )
    runtime = time.perf_counter() - started

    times = final_trajectory["time_s"]
    states = final_trajectory["state"]
    body_positions = np.array(
        [rot_z(-C.omega * time_s) @ state[:3] for time_s, state in zip(times, states)]
    )
    radii = np.linalg.norm(body_positions, axis=1)
    latitudes = np.degrees(np.arcsin(body_positions[:, 2] / radii))
    longitudes = (
        np.degrees(np.arctan2(body_positions[:, 1], body_positions[:, 0])) % 360.0
    )
    elevations = np.asarray(terrain.elevation_m(latitudes, longitudes))
    altitudes = radii - TERRAIN_META.reference_radius_m
    clearances = altitudes - elevations - C.ball_radius
    body_velocities = np.array(
        [
            rot_z(-C.omega * time_s)
            @ (
                state[3:]
                - np.cross(np.array([0.0, 0.0, C.omega]), state[:3])
            )
            for time_s, state in zip(times, states)
        ]
    )
    return_speed = float(np.linalg.norm(body_velocities[-1]))
    minimum_radius = float(np.min(radii))
    minimum_clearance = float(np.min(clearances))
    minimum_brillouin_clearance = minimum_radius - BRILLOUIN_RADIUS_M
    speed = float(np.linalg.norm(velocity))
    elevation_angle = math.degrees(math.asin(velocity[2] / speed))
    azimuth = math.degrees(math.atan2(velocity[1], velocity[0])) % 360.0

    table = pd.DataFrame(
        {
            "time_s": times,
            "x_body_m": body_positions[:, 0],
            "y_body_m": body_positions[:, 1],
            "z_body_m": body_positions[:, 2],
            "latitude_deg": latitudes,
            "longitude_deg_east": longitudes,
            "altitude_above_reference_m": altitudes,
            "terrain_elevation_m": elevations,
            "clearance_m": clearances,
        }
    )
    table.to_csv(DATA_OUT / "high_fidelity_case_trajectory.csv", index=False)
    diagnostics = {
        "scenario_label": SCENARIO_LABEL,
        "degree300_initialization": {
            "optimizer": "scipy.optimize.least_squares",
            "finite_difference_relative_step": 2e-5,
            "velocity_neu_m_s": correction.x.tolist(),
            "residual_vector_m": degree300_residual_m.tolist(),
            "residual_norm_m": float(np.linalg.norm(degree300_residual_m)),
            "function_evaluations": int(correction.nfev),
            "termination_message": correction.message,
        },
        "degree600_newton": {
            "residual_target_m": POSITION_RESIDUAL_TARGET_M,
            "maximum_iterations": 4,
            "history": correction_history,
            "final_velocity_neu_m_s": velocity.tolist(),
            "final_residual_vector_m": final_residual_vector.tolist(),
            "final_residual_norm_m": final_miss,
        },
    }
    (DATA_OUT / "high_fidelity_correction_diagnostics.json").write_text(
        json.dumps(diagnostics, indent=2) + "\n", encoding="utf-8"
    )

    summary = {
        "scenario_label": SCENARIO_LABEL,
        "frame_treatment": (
            "uniformly rotating aligned PA/ME frames; Earth and Sun fixed in "
            "the Moon-centred inertial frame over the short arc"
        ),
        "latitude_deg": LATITUDE_DEG,
        "longitude_deg_east": LONGITUDE_DEG,
        "launch_terrain_elevation_m": terrain_launch,
        "ball_center_height_above_terrain_m": height_above_terrain,
        "ball_center_altitude_above_reference_m": REFERENCE_ALTITUDE_M,
        "duration_s": duration,
        "gravity_degree": 600,
        "velocity_north_m_s": float(velocity[0]),
        "velocity_east_m_s": float(velocity[1]),
        "velocity_up_m_s": float(velocity[2]),
        "speed_m_s": speed,
        "elevation_deg": elevation_angle,
        "azimuth_deg_clockwise_from_north": azimuth,
        "scheduled_boundary_value_residual_m": final_miss,
        "scheduled_boundary_value_residual_vector_m": final_residual_vector.tolist(),
        "minimum_ldem64_clearance_m": minimum_clearance,
        "minimum_brillouin_clearance_m": minimum_brillouin_clearance,
        "return_relative_speed_m_s": return_speed,
        "kinetic_energy_at_return_j": 0.5 * C.ball_mass * return_speed**2,
        "initial_degree600_miss_before_update_m": initial_degree600_miss,
        "degree600_refinement_runtime_s": runtime,
        "status": (
            "NOMINAL_MODEL_RETURN"
            if final_miss <= 1.0
            and minimum_clearance > 0
            and minimum_brillouin_clearance >= 0
            else "INVALID"
        ),
    }
    (DATA_OUT / "high_fidelity_case_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    pd.DataFrame([summary]).to_csv(
        DATA_OUT / "high_fidelity_case_summary.csv", index=False
    )
    print(json.dumps(summary, indent=2), flush=True)

    setup()
    fig = plt.figure(figsize=(7.2, 5.8))
    ax3d = fig.add_subplot(2, 2, 1, projection="3d")
    unit = body_positions / C.radius
    sphere_u = np.linspace(0, 2 * np.pi, 70)
    sphere_v = np.linspace(0, np.pi, 35)
    ax3d.plot_surface(
        np.outer(np.cos(sphere_u), np.sin(sphere_v)),
        np.outer(np.sin(sphere_u), np.sin(sphere_v)),
        np.outer(np.ones_like(sphere_u), np.cos(sphere_v)),
        color="0.82",
        edgecolor="none",
        alpha=0.72,
    )
    ax3d.plot(unit[:, 0], unit[:, 1], unit[:, 2], color="C0", linewidth=1.2)
    ax3d.scatter(unit[0, 0], unit[0, 1], unit[0, 2], color="C2", s=18)
    ax3d.set_box_aspect((1, 1, 1))
    ax3d.set_title("Degree-600 body-fixed trajectory")
    ax3d.set_axis_off()

    ax_map = fig.add_subplot(2, 2, 2)
    texture = mpimg.imread(
        ROOT / "Web" / "public" / "assets" / "lroc_color_2k.jpg"
    )
    ax_map.imshow(texture, extent=(0, 360, -90, 90), origin="upper", alpha=0.65)
    discontinuity = np.abs(np.diff(longitudes)) > 180
    start = 0
    for stop in np.flatnonzero(discontinuity) + 1:
        ax_map.plot(longitudes[start:stop], latitudes[start:stop], color="C0")
        start = stop
    ax_map.plot(longitudes[start:], latitudes[start:], color="C0")
    ax_map.scatter([LONGITUDE_DEG], [LATITUDE_DEG], color="C2", s=16, zorder=3)
    ax_map.set_xlim(0, 360)
    ax_map.set_ylim(-90, 90)
    ax_map.set_xlabel("Longitude (deg E)")
    ax_map.set_ylabel("Latitude (deg)")
    ax_map.set_title("LDEM64 ground track")

    ax_clearance = fig.add_subplot(2, 1, 2)
    ax_clearance.plot(times / 60.0, clearances / 1000.0, color="C0")
    ax_clearance.axhline(0.0, color="C3", linestyle="--", linewidth=0.8)
    ax_clearance.set_xlabel("Time after launch (min)")
    ax_clearance.set_ylabel("Ball-to-LDEM64 clearance (km)")
    ax_clearance.set_title(
        f"Global collision test: minimum clearance {minimum_clearance / 1000:.2f} km"
    )
    save(fig, "fig17_high_fidelity_return")


if __name__ == "__main__":
    main()
