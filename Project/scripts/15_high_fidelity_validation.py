#!/usr/bin/env python3
"""Cross-check the nominal degree-600 solution numerically and spectrally."""
from __future__ import annotations

import json
from pathlib import Path
import runpy

import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar

from lunar_gravity import GRAILGravity, BRILLOUIN_RADIUS_M
from lunar_terrain import LDEM64, META as TERRAIN_META
from orbital_home_run import C, rot_z


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "output"


def main() -> None:
    case = runpy.run_path(
        str(ROOT / "scripts" / "10_high_fidelity_case.py"),
        run_name="high_fidelity_case_module",
    )
    propagate = case["propagate"]
    make_acceleration = case["acceleration"]
    initial_state = case["initial_state"]
    summary = json.loads(
        (OUT / "high_fidelity_case_summary.json").read_text(encoding="utf-8")
    )
    velocity = np.array(
        [
            summary["velocity_north_m_s"],
            summary["velocity_east_m_s"],
            summary["velocity_up_m_s"],
        ]
    )
    duration = float(summary["duration_s"])
    site_body = initial_state(velocity)[1]
    terrain = LDEM64()

    baseline_atol = np.array([2e-4, 2e-4, 2e-4, 2e-7, 2e-7, 2e-7])
    strict_atol = np.array([5e-5, 5e-5, 5e-5, 5e-8, 5e-8, 5e-8])
    configurations = [
        ("DOP853_step10", 600, "DOP853", 10.0, 2e-10, baseline_atol),
        ("DOP853_step5", 600, "DOP853", 5.0, 2e-10, baseline_atol),
        ("DOP853_step2p5", 600, "DOP853", 2.5, 2e-10, baseline_atol),
        ("RK45_strict_step2p5", 600, "RK45", 2.5, 5e-11, strict_atol),
        ("degree900", 900, "DOP853", 10.0, 2e-10, baseline_atol),
        ("degree1200", 1200, "DOP853", 10.0, 2e-10, baseline_atol),
    ]
    rows: list[dict[str, object]] = []
    final_positions: dict[str, np.ndarray] = {}
    for name, degree, method, max_step, rtol, atol in configurations:
        print(f"Running {name}", flush=True)
        model = make_acceleration(GRAILGravity(maximum_degree=degree))
        result = propagate(
            velocity,
            duration + 30.0,
            model,
            samples=2002,
            method=method,
            max_step_s=max_step,
            rtol=rtol,
            atol=atol,
            dense_output=True,
        )
        dense = result["dense_solution"]
        fixed_state = np.asarray(dense(duration), dtype=float)
        fixed_site = rot_z(C.omega * duration) @ site_body
        fixed_residual = fixed_state[:3] - fixed_site
        final_positions[name] = fixed_state[:3]

        def distance_to_site(time_s: float) -> float:
            state = np.asarray(dense(time_s), dtype=float)
            body_position = rot_z(-C.omega * time_s) @ state[:3]
            return float(np.linalg.norm(body_position - site_body))

        closest = minimize_scalar(
            distance_to_site,
            bounds=(duration - 60.0, duration + 30.0),
            method="bounded",
            options={"xatol": 1e-7},
        )
        sample_times = np.linspace(0.0, duration, 2001)
        states = np.asarray(dense(sample_times), dtype=float).T
        body_positions = np.array(
            [
                rot_z(-C.omega * time_s) @ state[:3]
                for time_s, state in zip(sample_times, states)
            ]
        )
        radii = np.linalg.norm(body_positions, axis=1)
        latitudes = np.degrees(np.arcsin(body_positions[:, 2] / radii))
        longitudes = (
            np.degrees(np.arctan2(body_positions[:, 1], body_positions[:, 0]))
            % 360.0
        )
        elevations = np.asarray(terrain.elevation_m(latitudes, longitudes))
        clearances = (
            radii
            - TERRAIN_META.reference_radius_m
            - elevations
            - C.ball_radius
        )
        rows.append(
            {
                "configuration": name,
                "gravity_degree": degree,
                "integrator": method,
                "maximum_step_s": max_step,
                "rtol": rtol,
                "position_atol_m": float(atol[0]),
                "velocity_atol_m_s": float(atol[3]),
                "function_evaluations": int(result["nfev"]),
                "fixed_time_residual_x_m": float(fixed_residual[0]),
                "fixed_time_residual_y_m": float(fixed_residual[1]),
                "fixed_time_residual_z_m": float(fixed_residual[2]),
                "fixed_time_residual_norm_m": float(np.linalg.norm(fixed_residual)),
                "closest_return_time_s": float(closest.x),
                "closest_return_distance_m": float(closest.fun),
                "minimum_ldem64_clearance_m": float(np.min(clearances)),
                "minimum_brillouin_clearance_m": float(
                    np.min(radii) - BRILLOUIN_RADIUS_M
                ),
            }
        )

    table = pd.DataFrame(rows)
    numerical_reference = final_positions["DOP853_step2p5"]
    gravity_reference = final_positions["degree1200"]
    table["endpoint_difference_vs_dop853_step2p5_m"] = [
        float(np.linalg.norm(final_positions[name] - numerical_reference))
        for name in table["configuration"]
    ]
    table["observed_endpoint_discrepancy_vs_degree1200_m"] = [
        float(np.linalg.norm(final_positions[name] - gravity_reference))
        for name in table["configuration"]
    ]
    numerical = table[table["gravity_degree"] == 600].copy()
    gravity = table[
        table["configuration"].isin(["DOP853_step10", "degree900", "degree1200"])
    ].copy()
    numerical.to_csv(OUT / "high_fidelity_integrator_crosscheck.csv", index=False)
    gravity.to_csv(OUT / "high_fidelity_gravity_order_crosscheck.csv", index=False)
    summary_out = {
        "interpretation": (
            "Numerical differences and truncation-model discrepancies are "
            "reported separately; neither is a physical trajectory-error bound."
        ),
        "maximum_degree600_numerical_endpoint_difference_m": float(
            numerical["endpoint_difference_vs_dop853_step2p5_m"].max()
        ),
        "degree600_vs_degree1200_endpoint_discrepancy_m": float(
            table.loc[
                table["configuration"] == "DOP853_step10",
                "observed_endpoint_discrepancy_vs_degree1200_m",
            ].iloc[0]
        ),
        "degree900_vs_degree1200_endpoint_discrepancy_m": float(
            table.loc[
                table["configuration"] == "degree900",
                "observed_endpoint_discrepancy_vs_degree1200_m",
            ].iloc[0]
        ),
    }
    (OUT / "high_fidelity_validation_summary.json").write_text(
        json.dumps(summary_out, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary_out, indent=2), flush=True)


if __name__ == "__main__":
    main()
