#!/usr/bin/env python3
"""Conservative analytic clearance above the global LDEM64 elevation range."""

from __future__ import annotations
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from orbital_home_run import C, maximum_safe_angle, minimum_speed_for_clearance
from plotting import OUT, save, setup


def main() -> None:
    setup()
    hmax = C.highest_elevation
    launch_altitudes = np.array([1, 10, 100, 1000, 5000, hmax+1, 20000, 50000], dtype=float)
    rows = []
    for h in launch_altitudes:
        r0 = C.radius+h
        for obstacle_name, obstacle in [("mean sphere", C.radius), ("global terrain envelope", C.radius+hmax)]:
            if r0 <= obstacle:
                vmin = math.nan
                amax = math.nan
            else:
                vmin = minimum_speed_for_clearance(r0, 0.0, obstacle)
                vchosen = 1.05*math.sqrt(C.mu/r0)
                amax = maximum_safe_angle(r0, vchosen, obstacle)
            rows.append([h, obstacle_name, obstacle-C.radius, vmin, math.degrees(amax) if math.isfinite(amax) else math.nan])
    pd.DataFrame(rows, columns=["launch_altitude_m", "clearance_model", "obstacle_elevation_m", "minimum_horizontal_speed_m_s", "max_safe_angle_at_1p05_local_circular_deg"]).to_csv(OUT / "terrain_envelope_requirements.csv", index=False)

    heights = np.linspace(hmax+0.1, 60000, 500)
    vmins = np.array([minimum_speed_for_clearance(C.radius+h, 0.0, C.radius+hmax) for h in heights])
    fig, ax = plt.subplots(figsize=(3.45, 2.8))
    ax.plot(heights/1e3, vmins/1e3)
    ax.set_xlabel("Launch altitude above mean radius (km)")
    ax.set_ylabel("Minimum horizontal speed (km/s)")
    ax.set_title("Conservative global-terrain clearance")
    ax.axvline(hmax/1e3, linewidth=0.8, linestyle="--")
    save(fig, "fig11_terrain_envelope_speed")


if __name__ == "__main__":
    main()
