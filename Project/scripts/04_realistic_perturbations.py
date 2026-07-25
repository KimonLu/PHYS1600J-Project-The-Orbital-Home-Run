#!/usr/bin/env python3
from __future__ import annotations
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from orbital_home_run import (
    C, initial_state, osculating_elements_from_state, propagate,
    make_acceleration_model, rot_z, body_fixed_lat_lon
)
from plotting import setup, save, OUT


def main() -> None:
    setup()
    # A deliberately terrain-clear idealized stadium: 20 km above the mean sphere,
    # equatorial, eastward. We specify an inertial speed 1.03 times local circular,
    # then subtract the site's rotational velocity to obtain the launch speed
    # relative to the lunar surface.
    altitude = 20_000.0
    r0 = C.radius + altitude
    u_inertial = 1.03
    v_inertial = u_inertial*math.sqrt(C.mu/r0)
    v_surface = v_inertial-C.omega*r0
    state0 = initial_state(r0, 0.0, 0.0, v_surface, 0.0, math.radians(90), True)
    osc0 = osculating_elements_from_state(state0)
    T = osc0["period"]

    models = [
        ("M0 central", make_acceleration_model()),
        ("M1 degree-2", make_acceleration_model(include_degree2=True)),
        ("M2 degree-2 + Earth", make_acceleration_model(include_degree2=True, include_earth=True)),
        ("M3 degree-2 + Earth + Sun", make_acceleration_model(include_degree2=True, include_earth=True, include_sun=True)),
    ]
    trajectories = {}
    rows = []
    site_body = np.array([r0, 0.0, 0.0])
    for name, acc in models:
        tr = propagate(state0, T, acc, max_step=2.0, samples=3501)
        trajectories[name] = tr
        radii = np.linalg.norm(tr["state"][:, :3], axis=1)
        final = tr["state"][-1]
        state_pos_error = np.linalg.norm(final[:3]-state0[:3])
        state_vel_error = np.linalg.norm(final[3:]-state0[3:])
        site_i = rot_z(C.omega*T)@site_body
        site_miss = np.linalg.norm(final[:3]-site_i)
        lat, lon = body_fixed_lat_lon(final[:3], T)
        rows.append([
            name, T, (radii.min()-C.radius), (radii.max()-C.radius),
            state_pos_error, state_vel_error, site_miss,
            math.degrees(lat), math.degrees(lon),
        ])
    out = pd.DataFrame(rows, columns=[
        "model", "reference_period_s", "minimum_altitude_m", "maximum_altitude_m",
        "inertial_position_error_after_T_m", "inertial_velocity_error_after_T_m_s",
        "distance_to_rotating_launch_site_after_T_m", "final_body_fixed_latitude_deg",
        "final_body_fixed_longitude_deg"
    ])
    out.to_csv(OUT / "perturbation_model_comparison.csv", index=False)

    central_radius = np.linalg.norm(
        trajectories["M0 central"]["state"][:, :3], axis=1
    )
    fig, ax = plt.subplots(figsize=(3.45, 2.8))
    for name, tr in list(trajectories.items())[1:]:
        radial_residual = (
            np.linalg.norm(tr["state"][:, :3], axis=1) - central_radius
        )
        ax.plot(tr["t"]/60, radial_residual, label=name)
    ax.set_xlabel("Time (min)")
    ax.set_ylabel("Altitude residual from central model (m)")
    ax.set_title("Perturbations change the radial history")
    ax.legend(fontsize=6.5)
    save(fig, "fig09_perturbation_altitude")

    central = trajectories["M0 central"]["state"][:, :3]
    degree2_state = trajectories["M1 degree-2"]["state"][:, :3]
    earth_state = trajectories["M2 degree-2 + Earth"]["state"][:, :3]
    full_state = trajectories["M3 degree-2 + Earth + Sun"]["state"][:, :3]
    fig, ax = plt.subplots(figsize=(3.45, 2.8))
    incremental = [
        ("degree-2 minus central", degree2_state-central),
        ("Earth increment", earth_state-degree2_state),
        ("Sun increment", full_state-earth_state),
    ]
    for name, delta in incremental:
        sep = np.linalg.norm(delta, axis=1)
        ax.plot(trajectories["M0 central"]["t"]/60, sep, label=name)
    ax.set_yscale("log")
    ax.set_xlabel("Time (min)")
    ax.set_ylabel("Incremental position change (m)")
    ax.set_title("Each added force shown on its own scale")
    ax.legend(fontsize=6.5)
    save(fig, "fig10_perturbation_separation")

    # Order-of-magnitude accelerations at the launch point.
    from orbital_home_run import acceleration_central, acceleration_degree2, acceleration_third_body, MU_EARTH, MU_SUN, EARTH_MOON_DISTANCE, SUN_MOON_DISTANCE
    rvec = state0[:3]
    ac = acceleration_central(0, rvec)
    ad2 = acceleration_degree2(0, rvec)-ac
    ae = acceleration_third_body(rvec, np.array([EARTH_MOON_DISTANCE,0,0]), MU_EARTH)
    ass = acceleration_third_body(rvec, np.array([0,SUN_MOON_DISTANCE,0]), MU_SUN)
    pd.DataFrame([
        ["central lunar gravity", np.linalg.norm(ac), np.linalg.norm(ac)/np.linalg.norm(ac)],
        ["degree-2 correction", np.linalg.norm(ad2), np.linalg.norm(ad2)/np.linalg.norm(ac)],
        ["Earth tidal acceleration", np.linalg.norm(ae), np.linalg.norm(ae)/np.linalg.norm(ac)],
        ["Sun tidal acceleration", np.linalg.norm(ass), np.linalg.norm(ass)/np.linalg.norm(ac)],
    ], columns=["term", "acceleration_m_s2", "fraction_of_central"]).to_csv(OUT / "acceleration_scale_comparison.csv", index=False)

    pd.DataFrame([
        ["launch_altitude_m", altitude],
        ["inertial_speed_ratio", u_inertial],
        ["surface_relative_speed_m_s", v_surface],
        ["initial_osculating_period_s", T],
        ["initial_periapsis_altitude_m", osc0["rp"]-C.radius],
        ["initial_apolune_altitude_m", osc0["ra"]-C.radius],
    ], columns=["quantity", "value"]).to_csv(OUT / "perturbation_case_definition.csv", index=False)


if __name__ == "__main__":
    main()
