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

    # High-fidelity terrain and gravity products.
    ldem64_path = ROOT / "data" / "external" / "ldem64" / "ldem_64.img"
    if ldem64_path.exists():
        ldem64 = np.memmap(
            ldem64_path,
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
    else:
        lines.append(
            "SKIP: full LDEM64 source extrema (run download_science_data.py)"
        )

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
        and high_case["scheduled_boundary_value_residual_m"] < 0.01
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
    middle = sensitivity["middle_step_results"]
    lines.append(check(
        all(
            middle[name]["degree600_middle_vs_small_fixed_time_relative_difference"]
            < 0.01
            and middle[name]["degree600_middle_vs_small_free_time_relative_difference"]
            < 0.01
            and middle[name]["publish_degree300_proxy"]
            for name in ("speed", "elevation", "azimuth")
        ),
        "degree-600 Jacobian is step-stable and the degree-300 proxy passes the declared five-percent test",
    ))
    validation = json.loads(
        (ROOT / "data" / "output" / "high_fidelity_validation_summary.json").read_text(
            encoding="utf-8"
        )
    )
    lines.append(check(
        validation["maximum_degree600_numerical_endpoint_difference_m"] < 1e-4,
        "DOP853 step refinement and strict RK45 agree below 0.1 mm at the endpoint",
    ))
    continuation = pd.read_csv(
        ROOT / "data" / "output" / "near_surface_height_continuation.csv"
    )
    accepted = continuation.loc[
        continuation["status"] == "NOMINAL_RETURN_FOUND_IN_BOUNDED_SEARCH"
    ]
    lines.append(check(
        len(continuation) == 6
        and set(continuation["height_above_local_ldem64_m"])
        == {2_000, 4_000, 8_000, 12_000, 16_000, 19_243}
        and accepted["height_above_local_ldem64_m"].tolist() == [19_243],
        "bounded height continuation covers 2 km and above and accepts only the 19.243 km candidate",
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
