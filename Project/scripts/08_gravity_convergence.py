#!/usr/bin/env python3
"""Quantify trajectory displacement caused by GRGM1200B truncation degree.

The experiment uses the same initial state for every force model and compares
the Cartesian state after one central-model synodic circuit.  Two circular
reference altitudes are included: 12 km (terrain-critical) and 30 km
(operationally safer).  Degree 1200 is the reference, not an assertion of
ground truth; model covariance, reference frames, and unmodelled forces remain.
"""
from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
import math
from pathlib import Path
import time

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp

from lunar_gravity import GRAILGravity
from lunar_gravity import BRILLOUIN_RADIUS_M
from orbital_home_run import C, initial_state
from plotting import OUT, save, setup


ROOT = Path(__file__).resolve().parents[1]
DATA_OUT = ROOT / "data" / "output"
DEGREES = (2, 10, 20, 50, 100, 200, 300, 400, 600, 900, 1200)
ALTITUDES_M = (12_000.0, 15_000.0, 30_000.0)


@dataclass(frozen=True)
class Job:
    altitude_m: float
    degree: int
    max_step_s: float = 10.0


def propagate_job(job: Job) -> dict[str, float]:
    radius = C.radius + job.altitude_m
    mean_motion = math.sqrt(C.mu / radius**3)
    duration = 2.0 * math.pi / (mean_motion - C.omega)
    surface_relative_speed = math.sqrt(C.mu / radius) - C.omega * radius
    state0 = initial_state(
        radius,
        latitude=0.0,
        longitude=0.0,
        speed_surface=surface_relative_speed,
        gamma=0.0,
        azimuth=0.5 * math.pi,
        include_surface_rotation=True,
    )
    gravity = GRAILGravity(maximum_degree=job.degree)

    def rhs(time_s: float, state: np.ndarray) -> np.ndarray:
        return np.concatenate(
            [state[3:], gravity.inertial_acceleration(time_s, state[:3])]
        )

    started = time.perf_counter()
    solution = solve_ivp(
        rhs,
        (0.0, duration),
        state0,
        method="DOP853",
        t_eval=np.linspace(0.0, duration, 401),
        max_step=job.max_step_s,
        rtol=2e-10,
        atol=np.array([2e-4, 2e-4, 2e-4, 2e-7, 2e-7, 2e-7]),
    )
    runtime = time.perf_counter() - started
    if not solution.success:
        raise RuntimeError(solution.message)
    positions = solution.y[:3].T
    radii = np.linalg.norm(positions, axis=1)
    return {
        "altitude_m": job.altitude_m,
        "degree": job.degree,
        "duration_s": duration,
        "surface_relative_speed_m_s": surface_relative_speed,
        "final_x_m": solution.y[0, -1],
        "final_y_m": solution.y[1, -1],
        "final_z_m": solution.y[2, -1],
        "final_vx_m_s": solution.y[3, -1],
        "final_vy_m_s": solution.y[4, -1],
        "final_vz_m_s": solution.y[5, -1],
        "minimum_radius_m": float(np.min(radii)),
        "maximum_radius_m": float(np.max(radii)),
        "runtime_s": runtime,
        "function_evaluations": int(solution.nfev),
    }


def add_reference_differences(table: pd.DataFrame) -> pd.DataFrame:
    pieces = []
    for altitude, group in table.groupby("altitude_m"):
        group = group.copy()
        reference = group.loc[group["degree"].idxmax()]
        ref_position = reference[["final_x_m", "final_y_m", "final_z_m"]].to_numpy(float)
        ref_velocity = reference[
            ["final_vx_m_s", "final_vy_m_s", "final_vz_m_s"]
        ].to_numpy(float)
        positions = group[["final_x_m", "final_y_m", "final_z_m"]].to_numpy(float)
        velocities = group[
            ["final_vx_m_s", "final_vy_m_s", "final_vz_m_s"]
        ].to_numpy(float)
        group["position_difference_vs_degree1200_m"] = np.linalg.norm(
            positions - ref_position, axis=1
        )
        group["velocity_difference_vs_degree1200_m_s"] = np.linalg.norm(
            velocities - ref_velocity, axis=1
        )
        group["radial_envelope_difference_vs_degree1200_m"] = np.maximum(
            np.abs(group["minimum_radius_m"] - reference["minimum_radius_m"]),
            np.abs(group["maximum_radius_m"] - reference["maximum_radius_m"]),
        )
        group["minimum_brillouin_clearance_m"] = (
            group["minimum_radius_m"] - BRILLOUIN_RADIUS_M
        )
        pieces.append(group)
    return pd.concat(pieces, ignore_index=True).sort_values(
        ["altitude_m", "degree"]
    )


def choose_degree(table: pd.DataFrame, threshold_m: float = 10.0) -> int:
    """Smallest degree meeting the threshold on formally valid reference arcs."""
    valid_altitudes = table.loc[
        table["minimum_brillouin_clearance_m"] >= 0, "altitude_m"
    ].unique()
    for degree in DEGREES:
        rows = table[
            (table["degree"] == degree)
            & (table["altitude_m"].isin(valid_altitudes))
        ]
        if len(rows) == len(valid_altitudes) and (
            rows["position_difference_vs_degree1200_m"] <= threshold_m
        ).all():
            return int(degree)
    return 1200


def make_figure(table: pd.DataFrame, selected_degree: int) -> None:
    setup()
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0))
    for altitude, group in table.groupby("altitude_m"):
        label = f"{altitude / 1000:g} km reference altitude"
        axes[0].loglog(
            group["degree"],
            np.maximum(group["position_difference_vs_degree1200_m"], 1e-4),
            marker="o",
            markersize=3,
            label=label,
        )
        axes[1].loglog(
            group["degree"],
            group["runtime_s"],
            marker="o",
            markersize=3,
            label=label,
        )
    axes[0].axhline(10.0, color="0.35", linestyle="--", linewidth=0.8)
    axes[0].axvline(selected_degree, color="C3", linestyle=":", linewidth=1.0)
    axes[0].set_xlabel("Maximum spherical-harmonic degree")
    axes[0].set_ylabel("Final position difference from degree 1200 (m)")
    axes[0].set_title("Trajectory convergence")
    axes[0].legend(fontsize=6.5)
    axes[1].set_xlabel("Maximum spherical-harmonic degree")
    axes[1].set_ylabel("Propagation runtime (s)")
    axes[1].set_title("Accuracy--cost trade-off")
    axes[1].legend(fontsize=6.5)
    save(fig, "fig16_gravity_degree_convergence")


def main() -> None:
    jobs = [Job(altitude, degree) for altitude in ALTITUDES_M for degree in DEGREES]
    existing_path = DATA_OUT / "gravity_degree_convergence.csv"
    rows: list[dict[str, float]] = []
    completed: set[tuple[float, int]] = set()
    if existing_path.exists():
        existing = pd.read_csv(existing_path)
        raw_columns = [
            "altitude_m",
            "degree",
            "duration_s",
            "surface_relative_speed_m_s",
            "final_x_m",
            "final_y_m",
            "final_z_m",
            "final_vx_m_s",
            "final_vy_m_s",
            "final_vz_m_s",
            "minimum_radius_m",
            "maximum_radius_m",
            "runtime_s",
            "function_evaluations",
        ]
        for record in existing[raw_columns].to_dict(orient="records"):
            key = (float(record["altitude_m"]), int(record["degree"]))
            if key[0] in ALTITUDES_M and key[1] in DEGREES:
                rows.append(record)
                completed.add(key)
    jobs = [
        job for job in jobs if (job.altitude_m, job.degree) not in completed
    ]
    if completed:
        print(f"Reusing {len(completed)} completed propagations", flush=True)
    workers = min(4, len(jobs))
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(propagate_job, job): job for job in jobs}
        for future in as_completed(futures):
            job = futures[future]
            row = future.result()
            rows.append(row)
            print(
                f"alt={job.altitude_m / 1000:g} km degree={job.degree:4d} "
                f"runtime={row['runtime_s']:.2f} s",
                flush=True,
            )
    table = add_reference_differences(pd.DataFrame(rows))
    selected = choose_degree(table)
    table["selected_main_degree"] = selected
    DATA_OUT.mkdir(parents=True, exist_ok=True)
    table.to_csv(DATA_OUT / "gravity_degree_convergence.csv", index=False)
    pd.DataFrame(
        [
            {
                "reference_degree": 1200,
                "selection_threshold_m": 10.0,
                "selected_main_degree": selected,
                "terrain_model": "LOLA LDEM64",
                "interpretation": (
                    "Numerical truncation comparison only; not a total "
                    "physical uncertainty bound."
                ),
            }
        ]
    ).to_csv(DATA_OUT / "gravity_model_selection.csv", index=False)
    make_figure(table, selected)
    print(f"SELECTED DEGREE: {selected}")


if __name__ == "__main__":
    main()
