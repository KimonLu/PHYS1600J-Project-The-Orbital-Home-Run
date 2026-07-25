"""Core mechanics for the PHYS1600J Orbital Home Run project.

All internal calculations use SI units. The module deliberately separates
closed-form two-body formulae from numerical propagation so that every
simulation can be verified against an analytic benchmark.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import csv
import json
import math
from typing import Callable

import numpy as np
from scipy.integrate import solve_ivp

ROOT = Path(__file__).resolve().parents[1]


@dataclass(frozen=True)
class Constants:
    mu: float
    radius: float
    rotation_period: float
    omega: float
    ball_mass: float
    ball_radius: float
    statcast_speed: float
    highest_elevation: float

    @property
    def circular_speed(self) -> float:
        return math.sqrt(self.mu / self.radius)

    @property
    def escape_speed(self) -> float:
        return math.sqrt(2.0 * self.mu / self.radius)

    @property
    def circular_period(self) -> float:
        return 2.0 * math.pi * math.sqrt(self.radius**3 / self.mu)

    @property
    def surface_gravity(self) -> float:
        return self.mu / self.radius**2


def load_constants(path: Path | None = None) -> Constants:
    if path is None:
        path = ROOT / "data" / "input" / "physical_constants.json"
    raw = json.loads(path.read_text(encoding="utf-8"))
    moon = raw["moon"]
    ball = raw["baseball"]
    period = float(moon["sidereal_rotation_period_s"])
    return Constants(
        mu=float(moon["mu_m3_s2"]),
        radius=float(moon["mean_radius_m"]),
        rotation_period=period,
        omega=2.0 * math.pi / period,
        ball_mass=float(ball["representative_mass_kg"]),
        ball_radius=float(ball["representative_radius_m"]),
        statcast_speed=float(ball["statcast_record_exit_speed_m_s"]),
        highest_elevation=float(moon["global_highest_elevation_m"]),
    )


C = load_constants()


def load_degree2_coefficients(
    path: Path | None = None,
) -> tuple[float, float]:
    """Load the unnormalised degree-2 coefficients used by the propagator."""
    if path is None:
        path = ROOT / "data" / "input" / "gravity_degree2.csv"
    values: dict[str, float] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            if row["normalization"].strip().lower() != "unnormalized":
                raise ValueError("Degree-2 propagator requires unnormalised coefficients")
            values[row["coefficient"].strip()] = float(row["value"])
    if "J2" not in values or "C22" not in values:
        raise ValueError("gravity_degree2.csv must define J2 and C22")
    return values["J2"], values["C22"]


J2_DEGREE2, C22_DEGREE2 = load_degree2_coefficients()


def orbit_elements(r0: float, v: float, gamma: float, mu: float = C.mu) -> dict[str, float]:
    """Return planar osculating elements for a launch at radius r0.

    gamma is the flight-path angle measured from local horizontal; positive is
    upward. The returned true anomaly is the initial true anomaly measured from
    periapsis in the direction of motion.
    """
    vr = v * math.sin(gamma)
    vt = v * math.cos(gamma)
    energy = 0.5 * v * v - mu / r0
    h = r0 * vt
    p = h * h / mu
    e_cos_f = h * vt / mu - 1.0
    e_sin_f = h * vr / mu
    e = math.hypot(e_cos_f, e_sin_f)
    f0 = math.atan2(e_sin_f, e_cos_f) if e > 1e-15 else 0.0
    if energy < 0.0:
        a = -mu / (2.0 * energy)
        period = 2.0 * math.pi * math.sqrt(a**3 / mu)
        rp = a * (1.0 - e)
        ra = a * (1.0 + e)
    else:
        a = math.inf if abs(energy) < 1e-15 else -mu / (2.0 * energy)
        period = math.inf
        rp = p / (1.0 + e)
        ra = math.inf
    return {
        "energy": energy,
        "h": h,
        "p": p,
        "a": a,
        "e": e,
        "f0": f0,
        "rp": rp,
        "ra": ra,
        "period": period,
        "vr": vr,
        "vt": vt,
    }


def horizontal_surface_orbit(u: float, radius: float = C.radius, mu: float = C.mu) -> dict[str, float]:
    """Closed-form orbit for horizontal launch from the reference sphere."""
    if not (0.0 < u < math.sqrt(2.0)):
        raise ValueError("u must be between 0 and sqrt(2) for a bound orbit")
    vc = math.sqrt(mu / radius)
    v = u * vc
    a = radius / (2.0 - u * u)
    e = abs(u * u - 1.0)
    if u >= 1.0:
        rp = radius
        ra = radius * u * u / (2.0 - u * u)
    else:
        ra = radius
        rp = radius * u * u / (2.0 - u * u)
    period = 2.0 * math.pi * math.sqrt(a**3 / mu)
    return {"u": u, "v": v, "a": a, "e": e, "rp": rp, "ra": ra, "period": period}


def exact_surface_shortfall(u: float, gamma: float, radius: float = C.radius) -> dict[str, float]:
    """Exact second-intersection shortfall for a launch from r=R.

    For a positive outbound gamma and u>1, the ball intersects the same sphere
    after sweeping 2*pi-2*f0, so it lands behind the launch point by 2*f0.
    """
    denom = u * u * math.cos(gamma) ** 2 - 1.0
    numer = u * u * math.sin(gamma) * math.cos(gamma)
    f0 = math.atan2(numer, denom)
    if f0 < 0.0:
        f0 += math.pi
    angular_shortfall = 2.0 * f0
    return {
        "f0": f0,
        "angular_shortfall": angular_shortfall,
        "surface_shortfall": radius * angular_shortfall,
    }


def small_angle_shortfall(u: float, gamma: float, radius: float = C.radius) -> float:
    return 2.0 * radius * u * u / (u * u - 1.0) * abs(gamma)


def minimum_speed_for_clearance(
    launch_radius: float,
    gamma: float,
    obstacle_radius: float = C.radius,
    mu: float = C.mu,
) -> float:
    """Minimum speed whose osculating periapsis is at obstacle_radius.

    A finite result also requires launch_radius*cos(gamma)>obstacle_radius.
    """
    denominator = launch_radius**2 * math.cos(gamma) ** 2 - obstacle_radius**2
    numerator = 2.0 * mu * obstacle_radius * (1.0 - obstacle_radius / launch_radius)
    if denominator <= 0.0 or numerator < 0.0:
        return math.nan
    return math.sqrt(numerator / denominator)


def maximum_safe_angle(
    launch_radius: float,
    v: float,
    obstacle_radius: float = C.radius,
    mu: float = C.mu,
) -> float:
    """Maximum |gamma| for a bound orbit with periapsis >= obstacle_radius."""
    rhs = (
        obstacle_radius**2
        + 2.0 * mu * obstacle_radius * (1.0 - obstacle_radius / launch_radius) / v**2
    ) / launch_radius**2
    if rhs > 1.0:
        return math.nan
    if rhs <= 0.0:
        return 0.5 * math.pi
    return math.acos(math.sqrt(rhs))


def rotation_miss_distance(period: float, latitude: float, radius: float = C.radius, omega: float = C.omega) -> float:
    """Chord distance between a site and its inertial position after period."""
    delta = omega * period
    return 2.0 * radius * abs(math.cos(latitude)) * abs(math.sin(0.5 * delta))


def rotation_arc_shift(period: float, latitude: float, radius: float = C.radius, omega: float = C.omega) -> float:
    return radius * abs(math.cos(latitude)) * abs(omega * period)


def local_basis(latitude: float, longitude: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    cl, sl = math.cos(latitude), math.sin(latitude)
    co, so = math.cos(longitude), math.sin(longitude)
    er = np.array([cl * co, cl * so, sl], dtype=float)
    east = np.array([-so, co, 0.0], dtype=float)
    north = np.array([-sl * co, -sl * so, cl], dtype=float)
    return er, east, north


def initial_state(
    launch_radius: float,
    latitude: float,
    longitude: float,
    speed_surface: float,
    gamma: float,
    azimuth: float,
    include_surface_rotation: bool = True,
    omega: float = C.omega,
) -> np.ndarray:
    """Build an inertial Cartesian state.

    Azimuth is measured clockwise from north. Longitude is body-fixed longitude
    at t=0, when the body-fixed and inertial axes coincide.
    """
    er, east, north = local_basis(latitude, longitude)
    horizontal = math.cos(azimuth) * north + math.sin(azimuth) * east
    r = launch_radius * er
    v_rel = speed_surface * (math.cos(gamma) * horizontal + math.sin(gamma) * er)
    if include_surface_rotation:
        v_rot = np.cross(np.array([0.0, 0.0, omega]), r)
    else:
        v_rot = np.zeros(3)
    return np.concatenate([r, v_rel + v_rot])


def rot_z(angle: float) -> np.ndarray:
    c, s = math.cos(angle), math.sin(angle)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


def acceleration_central(t: float, r_i: np.ndarray, mu: float = C.mu) -> np.ndarray:
    del t
    norm = np.linalg.norm(r_i)
    return -mu * r_i / norm**3


def acceleration_degree2(
    t: float,
    r_i: np.ndarray,
    mu: float = C.mu,
    radius: float = C.radius,
    omega: float = C.omega,
    j2: float = J2_DEGREE2,
    c22: float = C22_DEGREE2,
) -> np.ndarray:
    """Central + unnormalised J2/C22 acceleration in a rotating PA frame."""
    body_from_inertial = rot_z(-omega * t)
    r_b = body_from_inertial @ r_i
    x, y, z = r_b
    r2 = float(np.dot(r_b, r_b))
    r = math.sqrt(r2)
    ax_coeff = 0.5 * j2 + 3.0 * c22
    ay_coeff = 0.5 * j2 - 3.0 * c22
    az_coeff = -j2
    q = ax_coeff * x * x + ay_coeff * y * y + az_coeff * z * z
    coeffs = np.array([ax_coeff, ay_coeff, az_coeff])
    a_b = -mu * r_b / r**3 + mu * radius**2 * (
        2.0 * coeffs * r_b / r**5 - 5.0 * q * r_b / r**7
    )
    return body_from_inertial.T @ a_b


def potential_degree2(
    t: float,
    r_i: np.ndarray,
    mu: float = C.mu,
    radius: float = C.radius,
    omega: float = C.omega,
    j2: float = J2_DEGREE2,
    c22: float = C22_DEGREE2,
) -> float:
    """Positive gravitational potential U whose inertial gradient is acceleration.

    Mechanical specific energy is 0.5*v^2-U.  Since the tesseral field rotates,
    energy alone is not conserved, but E-Omega*h_z is conserved for this
    uniformly rotating degree-2 model.
    """
    r_b = rot_z(-omega * t) @ np.asarray(r_i, dtype=float)
    x, y, z = r_b
    r = float(np.linalg.norm(r_b))
    q = (
        (0.5 * j2 + 3.0 * c22) * x * x
        + (0.5 * j2 - 3.0 * c22) * y * y
        - j2 * z * z
    )
    return mu / r + mu * radius**2 * q / r**5


def propagate(
    state0: np.ndarray,
    duration: float,
    acceleration: Callable[[float, np.ndarray], np.ndarray] = acceleration_central,
    max_step: float = 5.0,
    rtol: float = 2e-11,
    atol: float = 2e-6,
    samples: int = 2001,
) -> dict[str, np.ndarray]:
    def rhs(t: float, y: np.ndarray) -> np.ndarray:
        return np.concatenate([y[3:], acceleration(t, y[:3])])

    t_eval = np.linspace(0.0, duration, samples)
    solution = solve_ivp(
        rhs,
        (0.0, duration),
        state0,
        method="DOP853",
        t_eval=t_eval,
        max_step=max_step,
        rtol=rtol,
        atol=atol,
    )
    if not solution.success:
        raise RuntimeError(solution.message)
    return {"t": solution.t, "state": solution.y.T}


def closest_return_metrics(
    trajectory: dict[str, np.ndarray],
    launch_site_body: np.ndarray,
    omega: float = C.omega,
    skip_fraction: float = 0.5,
) -> dict[str, float]:
    """Find the post-half-orbit closest approach to the rotating launch site."""
    t = trajectory["t"]
    state = trajectory["state"]
    start = int(skip_fraction * len(t))
    dists = np.empty(len(t) - start)
    for k, tk in enumerate(t[start:]):
        site_i = rot_z(omega * tk) @ launch_site_body
        dists[k] = np.linalg.norm(state[start + k, :3] - site_i)
    idx_local = int(np.argmin(dists))
    idx = start + idx_local
    return {
        "closest_time": float(t[idx]),
        "closest_distance": float(dists[idx_local]),
        "minimum_radius": float(np.min(np.linalg.norm(state[:, :3], axis=1))),
        "final_radius": float(np.linalg.norm(state[-1, :3])),
    }


def state_error(a: np.ndarray, b: np.ndarray) -> tuple[float, float]:
    return float(np.linalg.norm(a[:3] - b[:3])), float(np.linalg.norm(a[3:] - b[3:]))

MU_EARTH = 3.986004418e14
MU_SUN = 1.32712440018e20
EARTH_MOON_DISTANCE = 384400e3
SUN_MOON_DISTANCE = 149597870700.0


def acceleration_third_body(
    r_i: np.ndarray,
    body_position_i: np.ndarray,
    body_mu: float,
) -> np.ndarray:
    """Differential third-body acceleration in a Moon-centred inertial frame."""
    rel = body_position_i - r_i
    return body_mu * (rel / np.linalg.norm(rel) ** 3 - body_position_i / np.linalg.norm(body_position_i) ** 3)


def make_acceleration_model(
    include_degree2: bool = False,
    include_earth: bool = False,
    include_sun: bool = False,
    earth_position_i: np.ndarray | None = None,
    sun_position_i: np.ndarray | None = None,
) -> Callable[[float, np.ndarray], np.ndarray]:
    """Construct a reproducible hierarchy of lunar acceleration models.

    Earth and Sun are held fixed during a single low-lunar-orbit propagation.
    This is an explicit short-arc approximation used only for an order-of-
    magnitude perturbation comparison; the optional high-fidelity workflow is
    documented separately in the project README.
    """
    if earth_position_i is None:
        earth_position_i = np.array([EARTH_MOON_DISTANCE, 0.0, 0.0])
    if sun_position_i is None:
        sun_position_i = np.array([0.0, SUN_MOON_DISTANCE, 0.0])

    def model(t: float, r_i: np.ndarray) -> np.ndarray:
        if include_degree2:
            acc = acceleration_degree2(t, r_i)
        else:
            acc = acceleration_central(t, r_i)
        if include_earth:
            acc = acc + acceleration_third_body(r_i, earth_position_i, MU_EARTH)
        if include_sun:
            acc = acc + acceleration_third_body(r_i, sun_position_i, MU_SUN)
        return acc

    return model


def body_fixed_lat_lon(r_i: np.ndarray, t: float, omega: float = C.omega) -> tuple[float, float]:
    r_b = rot_z(-omega * t) @ r_i
    radius = np.linalg.norm(r_b)
    latitude = math.asin(float(r_b[2]) / radius)
    longitude = math.atan2(float(r_b[1]), float(r_b[0]))
    return latitude, longitude


def osculating_elements_from_state(state: np.ndarray, mu: float = C.mu) -> dict[str, float]:
    r_vec = np.asarray(state[:3], dtype=float)
    v_vec = np.asarray(state[3:], dtype=float)
    r = np.linalg.norm(r_vec)
    v = np.linalg.norm(v_vec)
    h_vec = np.cross(r_vec, v_vec)
    h = np.linalg.norm(h_vec)
    e_vec = np.cross(v_vec, h_vec) / mu - r_vec / r
    e = np.linalg.norm(e_vec)
    energy = 0.5 * v * v - mu / r
    a = -mu / (2.0 * energy) if energy < 0 else math.inf
    rp = a * (1.0 - e) if math.isfinite(a) else h * h / mu / (1.0 + e)
    ra = a * (1.0 + e) if math.isfinite(a) else math.inf
    period = 2.0 * math.pi * math.sqrt(a**3 / mu) if math.isfinite(a) else math.inf
    return {"a": float(a), "e": float(e), "rp": float(rp), "ra": float(ra), "period": float(period), "energy": float(energy), "h": float(h)}
