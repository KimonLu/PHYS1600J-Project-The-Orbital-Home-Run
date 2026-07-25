#!/usr/bin/env python3
"""Local closest-return sensitivity around the degree-600 case solution."""
from __future__ import annotations

import json
import math
from pathlib import Path
import runpy

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from lunar_gravity import GRAILGravity
from orbital_home_run import C, rot_z
from plotting import save, setup


ROOT = Path(__file__).resolve().parents[1]
DATA_OUT = ROOT / "data" / "output"


def main() -> None:
    case_module = runpy.run_path(
        str(ROOT / "scripts" / "11_high_fidelity_case.py"),
        run_name="high_fidelity_case_module",
    )
    propagate = case_module["propagate"]
    make_acceleration = case_module["acceleration"]
    summary = json.loads(
        (DATA_OUT / "high_fidelity_case_summary.json").read_text(encoding="utf-8")
    )
    velocity = np.array(
        [
            summary["velocity_north_m_s"],
            summary["velocity_east_m_s"],
            summary["velocity_up_m_s"],
        ]
    )
    duration = float(summary["duration_s"])
    gravity = GRAILGravity(maximum_degree=300)
    model = make_acceleration(gravity)
    delta = 0.02
    jacobian = np.empty((3, 3))
    for column in range(3):
        plus = velocity.copy()
        minus = velocity.copy()
        plus[column] += delta
        minus[column] -= delta
        plus_final = propagate(plus, duration, model, samples=2)["state"][-1, :3]
        minus_final = propagate(minus, duration, model, samples=2)["state"][-1, :3]
        jacobian[:, column] = (plus_final - minus_final) / (2.0 * delta)

    terminal = propagate(velocity, duration, model, samples=2)["state"][-1]
    site_body = case_module["initial_state"](velocity)[1]
    site_inertial = rot_z(C.omega * duration) @ site_body
    site_velocity = np.cross(np.array([0.0, 0.0, C.omega]), site_inertial)
    relative_velocity = terminal[3:] - site_velocity
    direction = relative_velocity / np.linalg.norm(relative_velocity)
    projector = np.eye(3) - np.outer(direction, direction)
    projected = projector @ jacobian

    speed = float(np.linalg.norm(velocity))
    elevation = math.radians(summary["elevation_deg"])
    azimuth = math.radians(summary["azimuth_deg_clockwise_from_north"])
    dv_dspeed = velocity / speed
    dv_delevation = speed * np.array(
        [
            -math.sin(elevation) * math.cos(azimuth),
            -math.sin(elevation) * math.sin(azimuth),
            math.cos(elevation),
        ]
    )
    dv_dazimuth = speed * math.cos(elevation) * np.array(
        [-math.sin(azimuth), math.cos(azimuth), 0.0]
    )
    speed_gain = float(np.linalg.norm(projected @ dv_dspeed))
    elevation_gain = float(np.linalg.norm(projected @ dv_delevation))
    azimuth_gain = float(np.linalg.norm(projected @ dv_dazimuth))
    result = {
        "gravity_degree_for_jacobian": 300,
        "base_solution_gravity_degree": 600,
        "finite_difference_velocity_step_m_s": delta,
        "closest_miss_gain_m_per_m_s_speed": speed_gain,
        "closest_miss_gain_m_per_rad_elevation": elevation_gain,
        "closest_miss_gain_m_per_rad_azimuth": azimuth_gain,
        "one_m_tolerance_speed_m_s": 1.0 / speed_gain,
        "one_m_tolerance_elevation_deg": math.degrees(1.0 / elevation_gain),
        "one_m_tolerance_azimuth_deg": math.degrees(1.0 / azimuth_gain),
        "projected_jacobian_singular_values_s": np.linalg.svd(
            projected, compute_uv=False
        ).tolist(),
    }
    (DATA_OUT / "high_fidelity_case_sensitivity.json").write_text(
        json.dumps(result, indent=2), encoding="utf-8"
    )
    pd.DataFrame([result]).to_csv(
        DATA_OUT / "high_fidelity_case_sensitivity.csv", index=False
    )
    print(json.dumps(result, indent=2))

    setup()
    fig, ax = plt.subplots(figsize=(3.45, 2.8))
    gains_per_degree = np.deg2rad(1.0) * np.array(
        [elevation_gain, azimuth_gain]
    )
    ax.bar(["elevation", "azimuth"], gains_per_degree, color=["C1", "C2"])
    ax.set_ylabel("Closest miss from a 1 deg error (m)")
    ax.set_title("Degree-600 case: angular sensitivity")
    save(fig, "fig18_high_fidelity_sensitivity")


if __name__ == "__main__":
    main()
