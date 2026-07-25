#!/usr/bin/env python3
from __future__ import annotations
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from orbital_home_run import C, exact_surface_shortfall, small_angle_shortfall, orbit_elements, maximum_safe_angle
from plotting import setup, save, OUT


def main() -> None:
    setup()
    u_values = [1.01, 1.05, 1.2, 1.35]
    gamma = np.logspace(-9, -2.2, 500)
    rows = []
    fig, ax = plt.subplots(figsize=(3.45, 2.8))
    for u in u_values:
        exact = np.array([exact_surface_shortfall(u, float(g))["surface_shortfall"] for g in gamma])
        asym = np.array([small_angle_shortfall(u, float(g)) for g in gamma])
        for g, s, sa in zip(gamma, exact, asym):
            rows.append([u, g, math.degrees(g), s, sa, abs(sa-s)/s])
        ax.loglog(np.degrees(gamma), exact, label=fr"$u={u}$")
    ax.axhline(0.5, linewidth=0.8, linestyle="--")
    ax.axhline(100.0, linewidth=0.8, linestyle=":")
    ax.set_xlabel("Launch-angle error (deg)")
    ax.set_ylabel("Surface shortfall (m)")
    ax.set_title("A tiny upward error causes an early impact")
    ax.legend(ncol=2)
    save(fig, "fig04_angle_sensitivity")
    pd.DataFrame(rows, columns=["u", "gamma_rad", "gamma_deg", "exact_shortfall_m", "small_angle_shortfall_m", "relative_error"]).to_csv(OUT / "angle_sensitivity.csv", index=False)

    heights = np.logspace(-2, 5, 350)  # ball-centre height from 1 cm to 100 km
    fig, ax = plt.subplots(figsize=(3.45, 2.8))
    tol_rows = []
    for u in [1.01, 1.05, 1.2, 1.35]:
        angles = []
        for h in heights:
            r0 = C.radius + h
            v = u*math.sqrt(C.mu/r0)
            amax = maximum_safe_angle(r0, v, C.radius)
            angles.append(math.degrees(amax) if math.isfinite(amax) else math.nan)
            tol_rows.append([u, h, angles[-1]])
        ax.loglog(heights, angles, label=fr"$u_0={u}$")
    ax.axvline(C.ball_radius, linewidth=0.8, linestyle="--")
    ax.axvline(1.0, linewidth=0.8, linestyle=":")
    ax.set_xlabel("Ball-centre launch height above reference sphere (m)")
    ax.set_ylabel(r"Maximum safe $|\gamma|$ (deg)")
    ax.set_title("Finite height opens a narrow safe wedge")
    ax.legend(ncol=2)
    save(fig, "fig05_finite_height_tolerance")
    pd.DataFrame(tol_rows, columns=["u_local", "launch_height_m", "max_safe_angle_deg"]).to_csv(OUT / "finite_height_tolerance.csv", index=False)

    # Monte Carlo: practical launch dispersions around a 1 m high idealized tee.
    rng = np.random.default_rng(1600)
    n = 100000
    h = 1.0
    u_nom = 1.2
    sigma_u = 5e-4
    sigma_gamma = math.radians(0.05)
    us = rng.normal(u_nom, sigma_u, n)
    gammas = rng.normal(0.0, sigma_gamma, n)
    r0 = C.radius + h
    rp = np.empty(n)
    period = np.empty(n)
    for i, (u, g) in enumerate(zip(us, gammas)):
        v = u*math.sqrt(C.mu/r0)
        e = orbit_elements(r0, v, g)
        rp[i] = e["rp"]
        period[i] = e["period"]
    safe = rp >= C.radius
    mc = pd.DataFrame({"u_local": us, "gamma_deg": np.degrees(gammas), "periapsis_altitude_m": rp-C.radius, "period_s": period, "safe": safe})
    mc.to_csv(OUT / "monte_carlo_launches.csv", index=False)
    summary = pd.DataFrame([
        ["sample_size", n],
        ["nominal_u", u_nom],
        ["sigma_u", sigma_u],
        ["sigma_gamma_deg", math.degrees(sigma_gamma)],
        ["safe_fraction", float(np.mean(safe))],
        ["analytic_max_safe_angle_deg", math.degrees(maximum_safe_angle(r0, u_nom*math.sqrt(C.mu/r0), C.radius))],
    ], columns=["quantity", "value"])
    summary.to_csv(OUT / "monte_carlo_summary.csv", index=False)

    # Scatter a reproducible subset for legibility.
    idx = rng.choice(n, 5000, replace=False)
    fig, ax = plt.subplots(figsize=(3.45, 2.8))
    ax.scatter(mc.loc[idx, "gamma_deg"], mc.loc[idx, "periapsis_altitude_m"], s=3, alpha=0.25)
    ax.axhline(0, linewidth=1.0)
    ax.set_xlabel("Elevation error (deg)")
    ax.set_ylabel("Osculating periapsis altitude (m)")
    ax.set_title("Monte Carlo launch dispersion at 1 m height")
    ax.set_ylim(-3000, 10)
    save(fig, "fig06_monte_carlo_periapsis")


if __name__ == "__main__":
    main()
