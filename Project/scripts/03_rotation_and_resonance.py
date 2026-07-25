#!/usr/bin/env python3
from __future__ import annotations
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from orbital_home_run import C, horizontal_surface_orbit, rotation_arc_shift
from plotting import setup, save, OUT


def resonance_u(N: int) -> float:
    target_period = C.rotation_period/N
    q = (C.circular_period/target_period)**(2.0/3.0)
    u2 = 2.0-q
    return math.sqrt(u2) if 1.0 <= u2 < 2.0 else math.nan


def main() -> None:
    setup()
    lat = np.linspace(-90, 90, 361)
    cases = [1.001, 1.05, 1.2]
    rows = []
    fig, ax = plt.subplots(figsize=(3.45, 2.8))
    for u in cases:
        T = horizontal_surface_orbit(u)["period"]
        shift = np.array([rotation_arc_shift(T, math.radians(x))/1e3 for x in lat])
        for x, s in zip(lat, shift):
            rows.append([u, x, T/60, s])
        ax.plot(lat, shift, label=fr"$u={u}$")
    ax.set_xlabel("Launch latitude (deg)")
    ax.set_ylabel("One-orbit surface displacement (km)")
    ax.set_title("The rotating ballpark does not wait")
    ax.legend()
    save(fig, "fig07_rotation_miss_by_latitude")
    pd.DataFrame(rows, columns=["u", "latitude_deg", "period_min", "surface_arc_shift_km"]).to_csv(OUT / "rotation_miss_by_latitude.csv", index=False)

    max_n = int(math.floor(C.rotation_period/C.circular_period))
    resonance_rows = []
    for N in range(1, max_n+1):
        u = resonance_u(N)
        if not math.isfinite(u):
            continue
        o = horizontal_surface_orbit(u)
        resonance_rows.append([N, C.rotation_period/N/60, u, o["v"]/1e3, (o["ra"]-C.radius)/1e3, o["e"]])
    rdf = pd.DataFrame(resonance_rows, columns=["orbits_per_lunar_rotation_N", "period_min", "u", "speed_km_s", "apolune_altitude_km", "eccentricity"])
    rdf.to_csv(OUT / "rotation_resonances.csv", index=False)

    # Show high-order resonances that remain within a nominal low-lunar-orbit apolune.
    low = rdf[rdf["apolune_altitude_km"] < 5000].copy()
    fig, ax = plt.subplots(figsize=(3.45, 2.8))
    ax.plot(low["orbits_per_lunar_rotation_N"], low["apolune_altitude_km"])
    ax.set_yscale("log")
    ax.set_xlabel("Number of revolutions per lunar rotation, $N$")
    ax.set_ylabel("Apolune altitude (km)")
    ax.set_title("Exact multi-orbit return resonances")
    save(fig, "fig08_rotation_resonances")

    key = rdf.tail(12).copy()
    key.to_csv(OUT / "near_circular_resonances.csv", index=False)

    summary = pd.DataFrame([
        ["lunar_rotation_period_days", C.rotation_period/86400],
        ["equatorial_rotation_speed_m_s", C.omega*C.radius],
        ["near_surface_period_min", C.circular_period/60],
        ["near_surface_equatorial_shift_km", rotation_arc_shift(C.circular_period, 0.0)/1e3],
        ["maximum_integer_N", max_n],
    ], columns=["quantity", "value"])
    summary.to_csv(OUT / "rotation_summary.csv", index=False)


if __name__ == "__main__":
    main()
