#!/usr/bin/env python3
"""Independent consistency checks for the packaged results."""
from __future__ import annotations
import math
from pathlib import Path
import numpy as np

from orbital_home_run import (
    C, horizontal_surface_orbit, orbit_elements, maximum_safe_angle,
    exact_surface_shortfall, small_angle_shortfall, propagate, acceleration_central,
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

    REPORT.write_text("\n".join(lines)+"\n", encoding="utf-8")
    print(REPORT.read_text(encoding="utf-8"), end="")


if __name__ == "__main__":
    main()
