#!/usr/bin/env python3
"""Bounded near-surface continuation at the concrete launch location.

For each requested platform height, a degree-300 solve supplies an initial
velocity at the central-model synodic time (which lies in the declared
6200--7200 s window).  A degree-300 central-difference Jacobian then removes
one degree-600 model increment, after which the corrected state is propagated
directly at degree 600 and checked against LDEM64 and the Brillouin sphere.

This is a deliberately bounded, reproducible search: it tests one scheduled
time per height, not every possible velocity/time pair in the full window.
Failure therefore means only that this continuation did not find a feasible
solution in its declared domain.
"""
from __future__ import annotations

import json
import math
from pathlib import Path
import runpy

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import least_squares

from lunar_gravity import BRILLOUIN_RADIUS_M, GRAILGravity
from lunar_terrain import LDEM64, META as TERRAIN_META
from orbital_home_run import C, local_basis, rot_z
from plotting import save, setup


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "output"
LATITUDE_DEG = 5.4296875
LONGITUDE_DEG = 201.3671875
HEIGHTS_M = (2_000, 4_000, 8_000, 12_000, 16_000, 19_243)
RETURN_WINDOW_S = (6_200.0, 7_200.0)
JACOBIAN_STEP_M_S = 0.02
RETURN_RESIDUAL_LIMIT_M = 1.0


def main() -> None:
    case = runpy.run_path(
        str(ROOT / "scripts" / "10_high_fidelity_case.py"),
        run_name="high_fidelity_case_module",
    )
    make_acceleration = case["acceleration"]
    terrain = LDEM64()
    terrain_launch = float(terrain.elevation_m(LATITUDE_DEG, LONGITUDE_DEG))
    latitude = math.radians(LATITUDE_DEG)
    longitude = math.radians(LONGITUDE_DEG)
    radial, east, north = local_basis(latitude, longitude)
    omega_vector = np.array([0.0, 0.0, C.omega])

    gravity300 = GRAILGravity(maximum_degree=300)
    gravity600 = GRAILGravity(maximum_degree=600)
    model300 = make_acceleration(gravity300)
    model600 = make_acceleration(gravity600)

    def initial_state(height_m: float, velocity_neu: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        radius = TERRAIN_META.reference_radius_m + terrain_launch + height_m
        site = radius * radial
        velocity_relative = (
            velocity_neu[0] * north + velocity_neu[1] * east + velocity_neu[2] * radial
        )
        return np.concatenate([site, velocity_relative + np.cross(omega_vector, site)]), site

    def propagate(
        height_m: float,
        velocity_neu: np.ndarray,
        duration_s: float,
        model,
        samples: int = 2,
    ) -> np.ndarray:
        from scipy.integrate import solve_ivp

        state0, _ = initial_state(height_m, velocity_neu)

        def rhs(time_s: float, state: np.ndarray) -> np.ndarray:
            return np.concatenate([state[3:], model(time_s, state[:3])])

        solution = solve_ivp(
            rhs,
            (0.0, duration_s),
            state0,
            method="DOP853",
            t_eval=np.linspace(0.0, duration_s, samples),
            max_step=10.0,
            rtol=2e-10,
            atol=np.array([2e-4] * 3 + [2e-7] * 3),
        )
        if not solution.success:
            raise RuntimeError(solution.message)
        return solution.y.T

    checkpoint = OUT / "near_surface_height_continuation.csv"
    if checkpoint.exists():
        rows = (
            pd.read_csv(checkpoint)
            .loc[lambda frame: frame["height_above_local_ldem64_m"].isin(HEIGHTS_M)]
            .to_dict("records")
        )
    else:
        rows: list[dict[str, object]] = []
    completed = {float(row["height_above_local_ldem64_m"]) for row in rows}
    seed = np.array([0.0, 1660.0, 0.0])
    # Descending continuation starts from the already demonstrated platform.
    for height_m in reversed(HEIGHTS_M):
        if float(height_m) in completed:
            prior = next(
                row for row in rows
                if float(row["height_above_local_ldem64_m"]) == float(height_m)
            )
            if all(key in prior and pd.notna(prior[key]) for key in (
                "velocity_north_m_s", "velocity_east_m_s", "velocity_up_m_s"
            )):
                seed = np.array([
                    prior["velocity_north_m_s"],
                    prior["velocity_east_m_s"],
                    prior["velocity_up_m_s"],
                ], dtype=float)
            continue
        radius = TERRAIN_META.reference_radius_m + terrain_launch + height_m
        duration = 2.0 * math.pi / (math.sqrt(C.mu / radius**3) - C.omega)
        if not RETURN_WINDOW_S[0] <= duration <= RETURN_WINDOW_S[1]:
            rows.append(
                {
                    "height_above_local_ldem64_m": height_m,
                    "scheduled_time_s": duration,
                    "status": "OUTSIDE_DECLARED_RETURN_WINDOW",
                }
            )
            continue
        _, site = initial_state(height_m, seed)
        target = rot_z(C.omega * duration) @ site

        if height_m == HEIGHTS_M[-1]:
            nominal = json.loads(
                (OUT / "high_fidelity_case_summary.json").read_text(encoding="utf-8")
            )
            seed = np.array(
                [
                    nominal["velocity_north_m_s"],
                    nominal["velocity_east_m_s"],
                    nominal["velocity_up_m_s"],
                ]
            )

        def residual300(velocity_neu: np.ndarray) -> np.ndarray:
            return (propagate(height_m, velocity_neu, duration, model300)[-1, :3] - target) / 1000.0

        try:
            correction300 = least_squares(
                residual300,
                seed,
                diff_step=2e-5,
                xtol=2e-9,
                ftol=2e-9,
                gtol=2e-9,
                max_nfev=16,
            )
        except (ValueError, RuntimeError) as error:
            rows.append(
                {
                    "height_above_local_ldem64_m": height_m,
                    "reference_altitude_m": terrain_launch + height_m,
                    "scheduled_time_s": duration,
                    "status": "NO_RETURN_FOUND_IN_BOUNDED_SEARCH",
                    "failure_stage": "degree300_initialization",
                    "failure_message": str(error),
                }
            )
            pd.DataFrame(rows).sort_values("height_above_local_ldem64_m").to_csv(
                checkpoint, index=False
            )
            print(height_m, "degree300 initialization failed", error, flush=True)
            continue
        velocity300 = correction300.x
        seed = velocity300
        try:
            residual600_before = propagate(
                height_m, velocity300, duration, model600
            )[-1, :3] - target

            jacobian300 = np.empty((3, 3))
            for component in range(3):
                plus = velocity300.copy()
                minus = velocity300.copy()
                plus[component] += JACOBIAN_STEP_M_S
                minus[component] -= JACOBIAN_STEP_M_S
                jacobian300[:, component] = (
                    propagate(height_m, plus, duration, model300)[-1, :3]
                    - propagate(height_m, minus, duration, model300)[-1, :3]
                ) / (2.0 * JACOBIAN_STEP_M_S)
            velocity600 = velocity300 - np.linalg.solve(jacobian300, residual600_before)
            states600 = propagate(
                height_m, velocity600, duration, model600, samples=1441
            )
        except (ValueError, RuntimeError, np.linalg.LinAlgError) as error:
            rows.append(
                {
                    "height_above_local_ldem64_m": height_m,
                    "reference_altitude_m": terrain_launch + height_m,
                    "scheduled_time_s": duration,
                    "velocity_north_m_s": velocity300[0],
                    "velocity_east_m_s": velocity300[1],
                    "velocity_up_m_s": velocity300[2],
                    "degree300_initial_residual_m": float(
                        np.linalg.norm(correction300.fun) * 1000.0
                    ),
                    "status": "NO_RETURN_FOUND_IN_BOUNDED_SEARCH",
                    "failure_stage": "degree600_proxy_update_or_validation",
                    "failure_message": str(error),
                }
            )
            pd.DataFrame(rows).sort_values("height_above_local_ldem64_m").to_csv(
                checkpoint, index=False
            )
            print(height_m, "degree600 validation failed", error, flush=True)
            continue
        residual600 = states600[-1, :3] - target

        times = np.linspace(0.0, duration, len(states600))
        body_positions = np.array(
            [rot_z(-C.omega * time_s) @ state[:3] for time_s, state in zip(times, states600)]
        )
        radii = np.linalg.norm(body_positions, axis=1)
        latitudes = np.degrees(np.arcsin(body_positions[:, 2] / radii))
        longitudes = np.degrees(np.arctan2(body_positions[:, 1], body_positions[:, 0])) % 360.0
        elevations = np.asarray(terrain.elevation_m(latitudes, longitudes))
        clearance = radii - TERRAIN_META.reference_radius_m - elevations - C.ball_radius
        minimum_clearance = float(np.min(clearance))
        minimum_brillouin = float(np.min(radii) - BRILLOUIN_RADIUS_M)
        residual_norm = float(np.linalg.norm(residual600))
        if minimum_clearance <= 0.0:
            status = "IMPACT"
        elif minimum_brillouin <= 0.0:
            status = "BELOW_BRILLOUIN_DOMAIN"
        elif residual_norm <= RETURN_RESIDUAL_LIMIT_M:
            status = "NOMINAL_RETURN_FOUND_IN_BOUNDED_SEARCH"
        else:
            status = "NO_RETURN_FOUND_IN_BOUNDED_SEARCH"
        rows.append(
            {
                "height_above_local_ldem64_m": height_m,
                "reference_altitude_m": terrain_launch + height_m,
                "scheduled_time_s": duration,
                "velocity_north_m_s": velocity600[0],
                "velocity_east_m_s": velocity600[1],
                "velocity_up_m_s": velocity600[2],
                "degree300_initial_residual_m": float(np.linalg.norm(correction300.fun) * 1000.0),
                "degree600_residual_before_proxy_update_m": float(np.linalg.norm(residual600_before)),
                "degree600_residual_after_proxy_update_m": residual_norm,
                "minimum_ldem64_clearance_m": minimum_clearance,
                "minimum_brillouin_clearance_m": minimum_brillouin,
                "status": status,
            }
        )
        pd.DataFrame(rows).sort_values("height_above_local_ldem64_m").to_csv(
            checkpoint, index=False
        )
        print(height_m, status, residual_norm, minimum_clearance, minimum_brillouin, flush=True)

    frame = pd.DataFrame(rows).sort_values("height_above_local_ldem64_m")
    frame.to_csv(checkpoint, index=False)
    summary = {
        "search_domain": {
            "site_latitude_deg": LATITUDE_DEG,
            "site_longitude_deg_east": LONGITUDE_DEG,
            "site_ldem64_elevation_m": terrain_launch,
            "platform_heights_m": list(HEIGHTS_M),
            "return_window_s": list(RETURN_WINDOW_S),
            "time_candidates_per_height": 1,
            "time_candidate_definition": "central-gravity synodic time",
            "degree300_role": "velocity solve and proxy Jacobian",
            "degree600_role": "direct propagation before and after one proxy update",
            "return_residual_limit_m": RETURN_RESIDUAL_LIMIT_M,
        },
        "interpretation": (
            "Failure means only that this specified one-time-per-height continuation "
            "did not find a terrain-safe, Brillouin-valid degree-600 return."
        ),
        "status_counts": frame["status"].value_counts().to_dict(),
    }
    (OUT / "near_surface_height_continuation_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )

    setup()
    fig, (clearance_axis, residual_axis) = plt.subplots(
        2, 1, figsize=(8.8, 6.4), sharex=True,
        gridspec_kw={"height_ratios": [1.45, 1.0]},
    )
    x_km = frame["height_above_local_ldem64_m"].to_numpy(dtype=float) / 1000.0
    clearance_axis.scatter(
        x_km,
        frame["minimum_ldem64_clearance_m"].to_numpy(dtype=float) / 1000.0,
        s=70,
        marker="o",
        label="Minimum LDEM64 clearance",
    )
    clearance_axis.scatter(
        x_km,
        frame["minimum_brillouin_clearance_m"].to_numpy(dtype=float) / 1000.0,
        s=70,
        marker="s",
        label="Minimum Brillouin clearance",
    )
    clearance_axis.axhspan(-60.0, 0.0, color="#cc6677", alpha=0.08)
    clearance_axis.axhline(0.0, color="k", linewidth=0.9)
    clearance_axis.set_ylabel("Minimum clearance (km)")
    clearance_axis.set_title("Bounded degree-600 height continuation")
    clearance_axis.legend()

    residual = frame["degree600_residual_after_proxy_update_m"].to_numpy(dtype=float)
    residual_axis.scatter(x_km, residual, s=65, color="#4477aa")
    residual_axis.axhline(
        RETURN_RESIDUAL_LIMIT_M, color="k", linestyle="--", linewidth=1.0,
        label="1 m nominal-return criterion",
    )
    residual_axis.set_yscale("log")
    residual_axis.set_ylabel("Scheduled residual (m)")
    residual_axis.set_xlabel("Platform height above local LDEM64 terrain (km)")
    residual_axis.legend()
    short_status = {
        "NOMINAL_RETURN_FOUND_IN_BOUNDED_SEARCH": "return",
        "BELOW_BRILLOUIN_DOMAIN": "below Brillouin",
        "NO_RETURN_FOUND_IN_BOUNDED_SEARCH": "no return found",
        "IMPACT": "impact",
    }
    for x_value, status in zip(x_km, frame["status"]):
        residual_axis.annotate(
            short_status.get(status, status),
            (x_value, 0.04),
            xycoords=("data", "axes fraction"),
            rotation=55,
            ha="left",
            va="bottom",
            fontsize=8,
        )
    fig.tight_layout()
    save(fig, "fig23_height_continuation")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
