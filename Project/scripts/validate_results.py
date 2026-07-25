#!/usr/bin/env python3
"""Independent consistency checks for the packaged results."""
from __future__ import annotations
import json
import math
from pathlib import Path
import numpy as np
import pandas as pd

from orbital_home_run import (
    C, horizontal_surface_orbit, orbit_elements, maximum_safe_angle,
    exact_surface_shortfall, small_angle_shortfall, propagate, acceleration_central,
    acceleration_degree2, potential_degree2, make_acceleration_model, initial_state,
)

ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "data" / "output" / "validation_report.txt"


def check(condition: bool, message: str) -> str:
    if not condition:
        raise AssertionError(message)
    return f"PASS: {message}"


def main() -> None:
    lines: list[str] = []
    lines.append(check(abs(C.surface_gravity-1.62421887656) < 1e-9, "mu/R^2 surface gravity"))
    lines.append(check(abs(C.escape_speed/math.sqrt(2)/C.circular_speed-1) < 1e-14, "escape speed equals sqrt(2) circular speed"))

    for u in (1.001, 1.05, 1.2, 1.35):
        closed = horizontal_surface_orbit(u)
        numerical_elements = orbit_elements(C.radius, u*C.circular_speed, 0.0)
        lines.append(check(abs(closed["a"]-numerical_elements["a"]) < 1e-6, f"closed-form semimajor axis at u={u}"))
        lines.append(check(abs(closed["ra"]-numerical_elements["ra"]) < 1e-6, f"closed-form apoapsis at u={u}"))

    gamma = 1e-7
    exact = exact_surface_shortfall(1.2, gamma)["surface_shortfall"]
    asym = small_angle_shortfall(1.2, gamma)
    lines.append(check(abs(exact-asym)/exact < 1e-6, "small-angle surface-shortfall asymptotic relation"))

    r0 = C.radius+1.0
    v = 1.2*math.sqrt(C.mu/r0)
    gmax = maximum_safe_angle(r0, v)
    el_boundary = orbit_elements(r0, v, gmax)
    lines.append(check(abs(el_boundary["rp"]-C.radius) < 2e-3, "finite-height boundary gives periapsis at reference sphere"))

    u = 1.2
    closed = horizontal_surface_orbit(u)
    state0 = np.array([C.radius, 0, 0, 0, u*C.circular_speed, 0], dtype=float)
    tr = propagate(state0, closed["period"], acceleration_central, max_step=30, samples=1201)
    closure = np.linalg.norm(tr["state"][-1]-state0)
    lines.append(check(closure < 1e-5, "central numerical orbit closes to better than 10 micrometres in state norm"))

    # Independent finite-difference check of the degree-2 potential gradient.
    point = np.array([C.radius + 33_000.0, 211_000.0, -137_000.0])
    time = 1234.5
    step = 0.25
    gradient = np.empty(3)
    for k in range(3):
        delta = np.zeros(3)
        delta[k] = step
        gradient[k] = (
            potential_degree2(time, point + delta)
            - potential_degree2(time, point - delta)
        ) / (2.0 * step)
    analytic_acceleration = acceleration_degree2(time, point)
    lines.append(check(
        np.linalg.norm(gradient-analytic_acceleration)
        / np.linalg.norm(analytic_acceleration) < 2e-8,
        "degree-2 acceleration matches an independent potential-gradient check",
    ))

    # The uniformly rotating tesseral field conserves E-Omega*h_z.
    state_d2 = initial_state(
        C.radius + 20_000.0,
        0.0,
        0.0,
        1.03 * math.sqrt(C.mu / (C.radius + 20_000.0))
        - C.omega * (C.radius + 20_000.0),
        0.0,
        0.5 * math.pi,
        True,
    )
    period_d2 = orbit_elements(
        np.linalg.norm(state_d2[:3]),
        np.linalg.norm(state_d2[3:]),
        0.0,
    )["period"]
    tr_d2 = propagate(
        state_d2,
        period_d2,
        acceleration_degree2,
        max_step=5.0,
        samples=1201,
    )
    jacobi_like = []
    for t, state in zip(tr_d2["t"], tr_d2["state"]):
        energy = 0.5 * np.dot(state[3:], state[3:]) - potential_degree2(
            float(t), state[:3]
        )
        hz = np.cross(state[:3], state[3:])[2]
        jacobi_like.append(energy - C.omega * hz)
    jacobi_like = np.asarray(jacobi_like)
    lines.append(check(
        np.max(np.abs((jacobi_like-jacobi_like[0]) / jacobi_like[0])) < 2e-12,
        "rotating degree-2 Jacobi integral is conserved",
    ))

    targeted = pd.read_csv(ROOT / "data" / "output" / "targeted_return_summary.csv")
    full_targeted = targeted.loc[targeted["model"] == "Full, targeted"].iloc[0]
    targeted_case_df = pd.read_csv(
        ROOT / "data" / "output" / "targeted_return_case_definition.csv"
    )
    targeted_case = dict(
        zip(targeted_case_df["quantity"], targeted_case_df["value"])
    )
    lines.append(check(
        full_targeted["scheduled_site_miss_m"] < 1e-3,
        "differentially corrected full-model trajectory hits the scheduled rotating site within 1 mm",
    ))
    lines.append(check(
        full_targeted["minimum_ldem4_clearance_m"] > 99.0,
        "targeted return retains the designed 100 m LDEM4 terrain margin",
    ))
    lines.append(check(
        abs(full_targeted["body_fixed_revolutions"]-1.0) < 1e-10
        and full_targeted["inertial_revolutions"] > 1.0,
        "targeted return completes one body-fixed circuit and more than one inertial revolution",
    ))
    targeted_state0 = initial_state(
        C.radius + targeted_case["ball_center_launch_altitude_m"],
        0.0,
        math.radians(targeted_case["stadium_longitude_deg_east"]),
        full_targeted["surface_relative_speed_m_s"],
        math.radians(full_targeted["elevation_angle_deg"]),
        0.5 * math.pi,
        True,
    )
    full_model = make_acceleration_model(
        include_degree2=True,
        include_earth=True,
        include_sun=True,
    )
    targeted_coarse = propagate(
        targeted_state0,
        full_targeted["scheduled_return_time_s"],
        full_model,
        max_step=4.0,
        samples=2,
    )["state"][-1]
    targeted_fine = propagate(
        targeted_state0,
        full_targeted["scheduled_return_time_s"],
        full_model,
        max_step=2.0,
        samples=2,
    )["state"][-1]
    lines.append(check(
        np.linalg.norm(targeted_coarse[:3]-targeted_fine[:3]) < 1e-4,
        "targeted return is converged below 0.1 mm when maximum step is halved",
    ))

    # High-fidelity terrain and gravity products.
    ldem64 = np.memmap(
        ROOT / "data" / "external" / "ldem64" / "ldem_64.img",
        dtype="<i2",
        mode="r",
        shape=(11_520, 23_040),
    )
    lines.append(check(
        float(ldem64.min()) * 0.5 == -9_114.5
        and float(ldem64.max()) * 0.5 == 10_757.0,
        "full LDEM64 source extrema reproduce the reported elevation range",
    ))
    del ldem64

    model_selection = pd.read_csv(
        ROOT / "data" / "output" / "gravity_model_selection.csv"
    ).iloc[0]
    lines.append(check(
        int(model_selection["reference_degree"]) == 1200
        and int(model_selection["selected_main_degree"]) == 600
        and float(model_selection["selection_threshold_m"]) == 10.0,
        "degree 600 is the lowest tested truncation selected against degree 1200",
    ))
    convergence = pd.read_csv(
        ROOT / "data" / "output" / "gravity_degree_convergence.csv"
    )
    for altitude in (15_000.0, 30_000.0):
        row600 = convergence.loc[
            (convergence["altitude_m"] == altitude)
            & (convergence["degree"] == 600)
        ].iloc[0]
        lines.append(check(
            row600["position_difference_vs_degree1200_m"] < 10.0
            and row600["minimum_brillouin_clearance_m"] > 0.0,
            f"degree-600 comparison is within 10 m and formally exterior at {altitude/1000:g} km",
        ))
    warning12 = convergence.loc[
        (convergence["altitude_m"] == 12_000.0)
        & (convergence["degree"] == 600)
    ].iloc[0]
    lines.append(check(
        warning12["minimum_brillouin_clearance_m"] < 0.0,
        "12 km comparison is correctly flagged below the conservative Brillouin sphere",
    ))

    high_case = json.loads(
        (ROOT / "data" / "output" / "high_fidelity_case_summary.json").read_text(
            encoding="utf-8"
        )
    )
    lines.append(check(
        high_case["gravity_degree"] == 600
        and high_case["scheduled_return_miss_m"] < 0.01
        and high_case["minimum_ldem64_clearance_m"] > 10_000.0
        and high_case["minimum_brillouin_clearance_m"] > 0.0,
        "specific degree-600 solution closes within 1 cm and clears terrain and Brillouin sphere",
    ))
    lines.append(check(
        199_000.0 < high_case["kinetic_energy_at_return_j"] < 201_000.0,
        "specific return carries approximately 200 kJ of kinetic energy",
    ))

    web_crosscheck = json.loads(
        (ROOT / "data" / "output" / "web_solver_crosscheck.json").read_text(
            encoding="utf-8"
        )
    )
    web_by_altitude = {
        float(row["case_altitude_m"]): float(
            row["position_difference_atlas_vs_direct_m"]
        )
        for row in web_crosscheck
    }
    lines.append(check(
        web_by_altitude[15_000.0] < 10.0
        and web_by_altitude[30_000.0] < 10.0,
        "browser gravity atlas remains within 10 m of direct degree-600 trajectories at valid test altitudes",
    ))
    lines.append(check(
        web_by_altitude[12_000.0] > 10.0,
        "browser 12 km cross-check exposes rather than hides the near-surface error increase",
    ))
    terrain_tile_validation = json.loads(
        (
            ROOT
            / "data"
            / "output"
            / "web_terrain_tile_validation.json"
        ).read_text(encoding="utf-8")
    )
    lines.append(check(
        terrain_tile_validation["maximum_absolute_difference_m"] == 0.0,
        "browser LDEM64 tiles reproduce direct terrain queries across tile seams",
    ))

    sensitivity = json.loads(
        (ROOT / "data" / "output" / "high_fidelity_case_sensitivity.json").read_text(
            encoding="utf-8"
        )
    )
    lines.append(check(
        0.0 < sensitivity["one_m_tolerance_speed_m_s"] < 0.1
        and 0.0 < sensitivity["one_m_tolerance_elevation_deg"] < 0.01
        and 0.0 < sensitivity["one_m_tolerance_azimuth_deg"] < 0.01,
        "specific-case one-metre launch tolerances are finite and in the reported precision regime",
    ))

    general_summary_path = (
        ROOT
        / "data"
        / "output"
        / "general_degree600"
        / "general_solver_summary.json"
    )
    if general_summary_path.exists():
        general_case = json.loads(general_summary_path.read_text(encoding="utf-8"))
        lines.append(check(
            general_case["status"] == "RETURN"
            and general_case["closest_return_distance_m"] < 0.01
            and general_case["impact_time_s"] is None
            and general_case["outside_harmonic_brillouin_sphere"],
            "public general-solver entry point reproduces the degree-600 return without impact",
        ))

    REPORT.write_text("\n".join(lines)+"\n", encoding="utf-8")
    print(REPORT.read_text(encoding="utf-8"), end="")


if __name__ == "__main__":
    main()
