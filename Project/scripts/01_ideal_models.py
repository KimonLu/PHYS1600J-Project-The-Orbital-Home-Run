#!/usr/bin/env python3
from __future__ import annotations
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, Arc

from orbital_home_run import C, horizontal_surface_orbit
from plotting import setup, save, OUT


def main() -> None:
    setup()
    constants = pd.DataFrame([
        ["Lunar gravitational parameter", C.mu, "m^3 s^-2"],
        ["Mean lunar radius", C.radius, "m"],
        ["Surface gravity from mu/R^2", C.surface_gravity, "m s^-2"],
        ["Circular speed at mean radius", C.circular_speed, "m s^-1"],
        ["Escape speed at mean radius", C.escape_speed, "m s^-1"],
        ["Circular period at mean radius", C.circular_period, "s"],
        ["Equatorial rotation speed", C.omega*C.radius, "m s^-1"],
        ["Reference baseball mass", C.ball_mass, "kg"],
        ["Statcast record exit speed", C.statcast_speed, "m s^-1"],
    ], columns=["quantity", "value", "unit"])
    constants.to_csv(OUT / "constants_summary.csv", index=False)

    u = np.linspace(1.00005, math.sqrt(2)-0.003, 800)
    rows = []
    for x in u:
        o = horizontal_surface_orbit(float(x))
        rows.append([x, o["v"], (o["ra"]-C.radius)/1e3, o["period"]/60, o["e"], o["a"]/1e3])
    df = pd.DataFrame(rows, columns=["u", "speed_m_s", "apolune_altitude_km", "period_min", "eccentricity", "semimajor_axis_km"])
    df.to_csv(OUT / "ideal_orbit_family.csv", index=False)

    sample_u = [1.001, 1.01, 1.05, 1.10, 1.20, 1.30]
    sample = []
    for x in sample_u:
        o = horizontal_surface_orbit(x)
        sample.append([x, o["v"]/1e3, (o["ra"]-C.radius)/1e3, o["period"]/60, o["e"]])
    pd.DataFrame(sample, columns=["u", "speed_km_s", "apolune_altitude_km", "period_min", "eccentricity"]).to_csv(OUT / "ideal_reference_cases.csv", index=False)

    fig, ax = plt.subplots(figsize=(3.45, 2.7))
    ax.plot(df["u"], df["apolune_altitude_km"])
    ax.set_yscale("log")
    ax.set_xlabel(r"Dimensionless launch speed $u=v/v_c$")
    ax.set_ylabel("Apolune altitude (km)")
    ax.set_title("Ideal horizontal surface-orbit family")
    ax.axvline(1.0, linewidth=0.8, linestyle="--")
    ax.text(1.005, 12, "circular boundary", fontsize=7)
    save(fig, "fig01_orbit_family_altitude")

    fig, ax = plt.subplots(figsize=(3.45, 2.7))
    ax.plot(df["u"], df["period_min"])
    ax.set_yscale("log")
    ax.set_xlabel(r"Dimensionless launch speed $u=v/v_c$")
    ax.set_ylabel("Orbital period (min)")
    ax.set_title("Period diverges near escape")
    save(fig, "fig02_orbit_family_period")

    # Scale-free geometry schematic for the theorem and early re-impact.
    fig, ax = plt.subplots(figsize=(5.9, 3.2))
    moon = Circle((0, 0), 1.0, facecolor="0.92", edgecolor="0.15", linewidth=1.0)
    ax.add_patch(moon)
    # Ellipse with two intersections with reference sphere.
    th = np.linspace(-0.13, 2*np.pi-0.13, 700)
    a, e = 1.55, 0.36
    p = a*(1-e*e)
    rr = p/(1+e*np.cos(th))
    x, y = rr*np.cos(th), rr*np.sin(th)
    ax.plot(x, y, linewidth=1.5, label="Kepler ellipse")
    f0 = 0.42
    p1 = np.array([math.cos(f0), math.sin(f0)])
    p2 = np.array([math.cos(-f0), math.sin(-f0)])
    ax.scatter([p1[0], p2[0]], [p1[1], p2[1]], zorder=5)
    tangent = np.array([-math.sin(f0), math.cos(f0)])
    radial = p1
    direction = math.cos(0.18)*tangent + math.sin(0.18)*radial
    ax.add_patch(FancyArrowPatch(p1, p1+0.55*direction, arrowstyle="->", mutation_scale=10, linewidth=1.3))
    ax.text(p1[0]+0.2, p1[1]+0.25, r"launch: $\gamma>0$", fontsize=8)
    ax.text(p2[0]+0.04, p2[1]-0.17, "second surface intersection", fontsize=8)
    ax.add_patch(Arc((0,0), 0.55,0.55, theta1=-math.degrees(f0), theta2=math.degrees(f0), linewidth=1.0))
    ax.text(0.31, -0.02, r"$2f_0$", fontsize=9)
    ax.text(-0.35, 0.04, "Moon", fontsize=10)
    ax.set_aspect("equal")
    ax.set_xlim(-1.3, 2.05)
    ax.set_ylim(-1.25, 1.35)
    ax.axis("off")
    ax.set_title("Why a non-tangential surface launch re-impacts early")
    save(fig, "fig03_surface_intersection_geometry")

    speed = pd.DataFrame([
        ["Fastest recorded MLB exit speed", C.statcast_speed],
        ["Ideal lunar circular speed", C.circular_speed],
        ["Ideal lunar escape speed", C.escape_speed],
    ], columns=["case", "speed_m_s"])
    speed["kinetic_energy_J"] = 0.5*C.ball_mass*speed["speed_m_s"]**2
    speed.to_csv(OUT / "speed_energy_comparison.csv", index=False)


if __name__ == "__main__":
    main()
