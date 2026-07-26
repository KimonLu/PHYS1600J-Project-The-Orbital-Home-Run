#!/usr/bin/env python3
"""Validate local return sensitivities at degrees 300 and 600.

The fixed-time Jacobian measures the scheduled endpoint.  Its projection
normal to the terminal site-relative velocity measures the first-order
closest miss when encounter time is free.  Neither quantity is a statistical
uncertainty or an independent simultaneous launcher tolerance.
"""
from __future__ import annotations

import json
import math
from pathlib import Path
import runpy
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from lunar_gravity import GRAILGravity
from orbital_home_run import C, rot_z
from plotting import save, setup


ROOT = Path(__file__).resolve().parents[1]
DATA_OUT = ROOT / "data" / "output"
DEGREES = (300, 600)
PARAMETERS = ("speed", "elevation", "azimuth")
STEPS = {
    "speed": (0.01, 0.02, 0.04),
    "elevation": (5e-6, 1e-5, 2e-5),
    "azimuth": (5e-6, 1e-5, 2e-5),
}


def plot_validation(frame: pd.DataFrame) -> None:
    """Plot order and finite-difference-step convergence from the machine table."""
    setup()
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.8), sharex=True)
    standard_steps = {
        "speed": 0.01,
        "elevation": math.radians(0.001),
        "azimuth": math.radians(0.001),
    }
    colors = {"speed": "C0", "elevation": "C1", "azimuth": "C2"}
    labels = {"speed": "speed", "elevation": "elevation", "azimuth": "azimuth"}
    for axis, column, title in [
        (axes[0], "fixed_time_gain_m_per_unit", "Fixed scheduled time"),
        (axes[1], "free_time_projected_gain_m_per_unit", "Free encounter time"),
    ]:
        for parameter in PARAMETERS:
            subset = frame.loc[frame["parameter"] == parameter].copy()
            middle = STEPS[parameter][1]
            for degree, linestyle in [(300, "--"), (600, "-")]:
                degree_rows = subset.loc[subset["gravity_degree"] == degree].sort_values(
                    "finite_difference_step"
                )
                x = degree_rows["finite_difference_step"].to_numpy(dtype=float) / middle
                y = (
                    degree_rows[column].to_numpy(dtype=float)
                    * standard_steps[parameter]
                )
                axis.plot(
                    x,
                    y,
                    marker="o",
                    linewidth=1.8,
                    linestyle=linestyle,
                    color=colors[parameter],
                    label=f"{labels[parameter]}, N={degree}",
                )
        axis.set_yscale("log")
        axis.set_xticks([0.5, 1.0, 2.0], ["half", "middle", "double"])
        axis.set_xlabel("Finite-difference step")
        axis.set_title(title)
        axis.grid(True, which="both", alpha=0.25)
    axes[0].set_ylabel("Response to 0.01 m/s or 0.001 deg (m)")
    handles, legend_labels = axes[1].get_legend_handles_labels()
    fig.legend(
        handles,
        legend_labels,
        loc="lower center",
        ncol=3,
        frameon=False,
        bbox_to_anchor=(0.5, -0.02),
    )
    fig.subplots_adjust(bottom=0.27, wspace=0.25)
    save(fig, "fig21_jacobian_validation")


def velocity_from_parameters(parameters: np.ndarray) -> np.ndarray:
    speed, elevation, azimuth = parameters
    return speed * np.array(
        [
            math.cos(elevation) * math.cos(azimuth),
            math.cos(elevation) * math.sin(azimuth),
            math.sin(elevation),
        ]
    )


def main() -> None:
    if "--plot-only" in sys.argv:
        plot_validation(pd.read_csv(DATA_OUT / "high_fidelity_jacobian_convergence.csv"))
        return
    case = runpy.run_path(
        str(ROOT / "scripts" / "10_high_fidelity_case.py"),
        run_name="high_fidelity_case_module",
    )
    propagate = case["propagate"]
    make_acceleration = case["acceleration"]
    initial_state = case["initial_state"]
    summary = json.loads(
        (DATA_OUT / "high_fidelity_case_summary.json").read_text(encoding="utf-8")
    )
    parameters = np.array(
        [
            summary["speed_m_s"],
            math.radians(summary["elevation_deg"]),
            math.radians(summary["azimuth_deg_clockwise_from_north"]),
        ]
    )
    velocity = velocity_from_parameters(parameters)
    duration = float(summary["duration_s"])
    site_body = initial_state(velocity)[1]
    site_inertial = rot_z(C.omega * duration) @ site_body
    site_velocity = np.cross(np.array([0.0, 0.0, C.omega]), site_inertial)

    models: dict[int, object] = {}
    terminals: dict[int, np.ndarray] = {}
    for degree in DEGREES:
        print(f"Loading degree-{degree} sensitivity model", flush=True)
        models[degree] = make_acceleration(GRAILGravity(maximum_degree=degree))
        terminals[degree] = propagate(
            velocity, duration, models[degree], samples=2
        )["state"][-1]

    relative_velocity_600 = terminals[600][3:] - site_velocity
    direction = relative_velocity_600 / np.linalg.norm(relative_velocity_600)
    projector = np.eye(3) - np.outer(direction, direction)

    rows: list[dict[str, float | int | str]] = []
    jacobians: dict[tuple[int, int], np.ndarray] = {}
    for degree in DEGREES:
        model = models[degree]
        for step_index in range(3):
            jacobian = np.empty((3, 3))
            for column, parameter in enumerate(PARAMETERS):
                step = STEPS[parameter][step_index]
                plus = parameters.copy()
                minus = parameters.copy()
                plus[column] += step
                minus[column] -= step
                plus_final = propagate(
                    velocity_from_parameters(plus), duration, model, samples=2
                )["state"][-1, :3]
                minus_final = propagate(
                    velocity_from_parameters(minus), duration, model, samples=2
                )["state"][-1, :3]
                derivative = (plus_final - minus_final) / (2.0 * step)
                jacobian[:, column] = derivative
                fixed_gain = float(np.linalg.norm(derivative))
                projected_derivative = projector @ derivative
                projected_gain = float(np.linalg.norm(projected_derivative))
                encounter_time_gain = float(
                    -np.dot(relative_velocity_600, derivative)
                    / np.dot(relative_velocity_600, relative_velocity_600)
                )
                rows.append(
                    {
                        "gravity_degree": degree,
                        "step_index": step_index,
                        "parameter": parameter,
                        "parameter_unit": "m/s" if parameter == "speed" else "rad",
                        "finite_difference_step": step,
                        "fixed_time_gain_m_per_unit": fixed_gain,
                        "free_time_projected_gain_m_per_unit": projected_gain,
                        "encounter_time_gain_s_per_unit": encounter_time_gain,
                        "one_m_fixed_time_tolerance": 1.0 / fixed_gain,
                        "one_m_free_time_tolerance": 1.0 / projected_gain,
                    }
                )
            jacobians[(degree, step_index)] = jacobian

    table = pd.DataFrame(rows)
    table.to_csv(DATA_OUT / "high_fidelity_jacobian_convergence.csv", index=False)
    chosen = table[table["step_index"] == 1].copy()
    chosen["one_m_fixed_time_tolerance_deg"] = np.where(
        chosen["parameter"] == "speed",
        np.nan,
        np.degrees(chosen["one_m_fixed_time_tolerance"]),
    )
    chosen["one_m_free_time_tolerance_deg"] = np.where(
        chosen["parameter"] == "speed",
        np.nan,
        np.degrees(chosen["one_m_free_time_tolerance"]),
    )
    chosen.to_csv(DATA_OUT / "high_fidelity_case_sensitivity.csv", index=False)

    comparison: dict[str, object] = {
        "common_projection_direction": (
            "degree-600 terminal velocity relative to the rotating site"
        ),
        "degree600_terminal_relative_speed_m_s": float(
            np.linalg.norm(relative_velocity_600)
        ),
        "one_m_sphere_one_sided_crossing_time_s": float(
            1.0 / np.linalg.norm(relative_velocity_600)
        ),
        "two_m_diameter_crossing_time_s": float(
            2.0 / np.linalg.norm(relative_velocity_600)
        ),
        "middle_step_results": {},
    }
    for parameter in PARAMETERS:
        p300 = chosen[
            (chosen["gravity_degree"] == 300)
            & (chosen["parameter"] == parameter)
        ].iloc[0]
        p600 = chosen[
            (chosen["gravity_degree"] == 600)
            & (chosen["parameter"] == parameter)
        ].iloc[0]
        small600 = table[
            (table["gravity_degree"] == 600)
            & (table["parameter"] == parameter)
            & (table["step_index"] == 0)
        ].iloc[0]
        comparison["middle_step_results"][parameter] = {
            "degree600_fixed_time_gain_m_per_unit": float(
                p600["fixed_time_gain_m_per_unit"]
            ),
            "degree600_free_time_gain_m_per_unit": float(
                p600["free_time_projected_gain_m_per_unit"]
            ),
            "degree300_vs_600_fixed_time_relative_difference": float(
                abs(p300["fixed_time_gain_m_per_unit"] - p600["fixed_time_gain_m_per_unit"])
                / p600["fixed_time_gain_m_per_unit"]
            ),
            "degree300_vs_600_free_time_relative_difference": float(
                abs(
                    p300["free_time_projected_gain_m_per_unit"]
                    - p600["free_time_projected_gain_m_per_unit"]
                )
                / p600["free_time_projected_gain_m_per_unit"]
            ),
            "degree600_middle_vs_small_fixed_time_relative_difference": float(
                abs(
                    p600["fixed_time_gain_m_per_unit"]
                    - small600["fixed_time_gain_m_per_unit"]
                )
                / p600["fixed_time_gain_m_per_unit"]
            ),
            "degree600_middle_vs_small_free_time_relative_difference": float(
                abs(
                    p600["free_time_projected_gain_m_per_unit"]
                    - small600["free_time_projected_gain_m_per_unit"]
                )
                / p600["free_time_projected_gain_m_per_unit"]
            ),
            "publish_degree300_proxy": bool(
                abs(
                    p300["free_time_projected_gain_m_per_unit"]
                    - p600["free_time_projected_gain_m_per_unit"]
                )
                / p600["free_time_projected_gain_m_per_unit"]
                <= 0.05
            ),
        }
    singular_values_fixed = np.linalg.svd(jacobians[(600, 1)], compute_uv=False)
    singular_values_projected = np.linalg.svd(
        projector @ jacobians[(600, 1)], compute_uv=False
    )
    comparison["degree600_middle_step_jacobian"] = {
        "fixed_time_matrix": jacobians[(600, 1)].tolist(),
        "projected_matrix": (projector @ jacobians[(600, 1)]).tolist(),
        "fixed_time_singular_values": singular_values_fixed.tolist(),
        "projected_singular_values": singular_values_projected.tolist(),
    }
    (DATA_OUT / "high_fidelity_case_sensitivity.json").write_text(
        json.dumps(comparison, indent=2) + "\n", encoding="utf-8"
    )

    plot_validation(frame)
    print(json.dumps(comparison, indent=2), flush=True)


if __name__ == "__main__":
    main()
