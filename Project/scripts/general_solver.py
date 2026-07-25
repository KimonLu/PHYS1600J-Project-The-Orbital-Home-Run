"""General terrain-aware lunar trajectory and return solver.

Inputs are a body-fixed launch site (latitude, east longitude, and ball-centre
height above local terrain) plus a surface-relative velocity (speed, elevation,
and azimuth clockwise from north).  Outputs include the inertial and body-fixed
3-D trajectories, the equirectangular ground track, first terrain collision,
and the closest post-launch approach to the rotating launch point.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
import math
from pathlib import Path
from typing import Callable, Literal

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from scipy.optimize import minimize_scalar

from lunar_gravity import BRILLOUIN_RADIUS_M, GRAILGravity
from lunar_terrain import LDEM64, META as TERRAIN_META
from orbital_home_run import (
    C,
    acceleration_central,
    acceleration_degree2,
    acceleration_third_body,
    initial_state,
    rot_z,
    EARTH_MOON_DISTANCE,
    MU_EARTH,
    MU_SUN,
    SUN_MOON_DISTANCE,
)


ROOT = Path(__file__).resolve().parents[1]
Status = Literal[
    "RETURN",
    "RETURN_UNCERTAIN",
    "IMPACT",
    "ESCAPE",
    "NO_RETURN_WITHIN_WINDOW",
    "NUMERICAL_FAILURE",
]


@dataclass(frozen=True)
class LaunchInput:
    latitude_deg: float
    longitude_deg_east: float
    height_above_terrain_m: float
    speed_m_s: float
    elevation_deg: float
    azimuth_deg_clockwise_from_north: float


@dataclass(frozen=True)
class ModelConfig:
    gravity_model: str = "grail"
    gravity_degree: int = 600
    include_earth_tide: bool = True
    include_solar_tide: bool = True
    duration_s: float = 8_000.0
    return_tolerance_m: float = 1.0
    model_uncertainty_m: float = 100.0
    minimum_return_time_s: float = 1_000.0
    max_step_s: float = 10.0
    output_samples: int = 2_001
    epoch_utc: str = "2026-01-01T00:00:00Z"


@dataclass
class TrajectoryResult:
    status: Status
    message: str
    launch: dict[str, float]
    model: dict[str, object]
    launch_terrain_elevation_m: float
    launch_radius_m: float
    closest_return_time_s: float | None
    closest_return_distance_m: float | None
    closest_return_relative_speed_m_s: float | None
    impact_time_s: float | None
    impact_latitude_deg: float | None
    impact_longitude_deg_east: float | None
    impact_speed_m_s: float | None
    minimum_clearance_m: float
    maximum_altitude_above_reference_m: float
    final_specific_energy_j_kg: float
    minimum_radius_m: float
    outside_harmonic_brillouin_sphere: bool
    trajectory: pd.DataFrame

    def summary_dict(self) -> dict[str, object]:
        result = asdict(self)
        result.pop("trajectory")
        return result

    def save(self, directory: Path | str, stem: str = "general_solver") -> None:
        directory = Path(directory)
        directory.mkdir(parents=True, exist_ok=True)
        self.trajectory.to_csv(directory / f"{stem}_trajectory.csv", index=False)
        (directory / f"{stem}_summary.json").write_text(
            json.dumps(self.summary_dict(), indent=2, ensure_ascii=False),
            encoding="utf-8",
        )


def _body_position(position_inertial_m: np.ndarray, time_s: float) -> np.ndarray:
    return rot_z(-C.omega * time_s) @ np.asarray(position_inertial_m)


def _lat_lon_radius_body(position_body_m: np.ndarray) -> tuple[float, float, float]:
    radius = float(np.linalg.norm(position_body_m))
    lat = math.degrees(math.asin(float(position_body_m[2]) / radius))
    lon = math.degrees(math.atan2(float(position_body_m[1]), float(position_body_m[0]))) % 360.0
    return lat, lon, radius


def build_acceleration(config: ModelConfig) -> Callable[[float, np.ndarray], np.ndarray]:
    gravity_name = config.gravity_model.strip().lower()
    grail: GRAILGravity | None = None
    if gravity_name == "grail":
        grail = GRAILGravity(maximum_degree=config.gravity_degree)
    elif gravity_name not in {"central", "degree2"}:
        raise ValueError("gravity_model must be central, degree2, or grail")

    earth_position = np.array([EARTH_MOON_DISTANCE, 0.0, 0.0])
    sun_position = np.array([0.0, SUN_MOON_DISTANCE, 0.0])

    def acceleration(time_s: float, position_i: np.ndarray) -> np.ndarray:
        if gravity_name == "central":
            value = acceleration_central(time_s, position_i)
        elif gravity_name == "degree2":
            value = acceleration_degree2(time_s, position_i)
        else:
            assert grail is not None
            value = grail.inertial_acceleration(
                time_s, position_i, degree=config.gravity_degree
            )
        if config.include_earth_tide:
            value = value + acceleration_third_body(
                position_i, earth_position, MU_EARTH
            )
        if config.include_solar_tide:
            value = value + acceleration_third_body(
                position_i, sun_position, MU_SUN
            )
        return value

    return acceleration


def solve_trajectory(
    launch: LaunchInput,
    config: ModelConfig,
    terrain: LDEM64 | None = None,
    acceleration: Callable[[float, np.ndarray], np.ndarray] | None = None,
) -> TrajectoryResult:
    """Propagate one launch and classify its first collision and closest return."""
    if terrain is None:
        terrain = LDEM64()
    if not -90.0 <= launch.latitude_deg <= 90.0:
        raise ValueError("latitude must lie in [-90, 90] degrees")
    if launch.height_above_terrain_m <= C.ball_radius:
        raise ValueError(
            "ball-centre height must exceed the baseball radius to start above terrain"
        )
    if launch.speed_m_s < 0.0:
        raise ValueError("speed must be non-negative")
    if config.duration_s <= 0.0 or config.output_samples < 2:
        raise ValueError("duration and output_samples must be positive")
    if config.minimum_return_time_s >= config.duration_s:
        raise ValueError("minimum_return_time_s must be below duration_s")

    longitude = launch.longitude_deg_east % 360.0
    terrain_at_launch = terrain.elevation_m(launch.latitude_deg, longitude)
    launch_radius = (
        TERRAIN_META.reference_radius_m
        + terrain_at_launch
        + launch.height_above_terrain_m
    )
    state0 = initial_state(
        launch_radius,
        math.radians(launch.latitude_deg),
        math.radians(longitude),
        launch.speed_m_s,
        math.radians(launch.elevation_deg),
        math.radians(launch.azimuth_deg_clockwise_from_north),
        include_surface_rotation=True,
    )
    launch_site_body = state0[:3].copy()
    if acceleration is None:
        acceleration = build_acceleration(config)

    def rhs(time_s: float, state: np.ndarray) -> np.ndarray:
        return np.concatenate([state[3:], acceleration(time_s, state[:3])])

    def clearance_event(time_s: float, state: np.ndarray) -> float:
        body = _body_position(state[:3], time_s)
        lat, lon, radius = _lat_lon_radius_body(body)
        surface_radius = TERRAIN_META.reference_radius_m + terrain.elevation_m(lat, lon)
        return radius - surface_radius - C.ball_radius

    clearance_event.terminal = True
    clearance_event.direction = -1.0

    try:
        solution = solve_ivp(
            rhs,
            (0.0, config.duration_s),
            state0,
            method="DOP853",
            max_step=config.max_step_s,
            rtol=2e-10,
            atol=np.array([2e-4, 2e-4, 2e-4, 2e-7, 2e-7, 2e-7]),
            events=clearance_event,
            dense_output=True,
        )
    except Exception as exc:
        empty = pd.DataFrame()
        return TrajectoryResult(
            "NUMERICAL_FAILURE",
            f"Propagation raised {type(exc).__name__}: {exc}",
            asdict(launch),
            asdict(config),
            terrain_at_launch,
            launch_radius,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            math.nan,
            math.nan,
            math.nan,
            math.nan,
            False,
            empty,
        )
    if not solution.success or solution.sol is None:
        empty = pd.DataFrame()
        return TrajectoryResult(
            "NUMERICAL_FAILURE",
            solution.message,
            asdict(launch),
            asdict(config),
            terrain_at_launch,
            launch_radius,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            math.nan,
            math.nan,
            math.nan,
            math.nan,
            False,
            empty,
        )

    end_time = float(solution.t[-1])
    times = np.linspace(0.0, end_time, config.output_samples)
    states = solution.sol(times).T

    body_positions = np.empty((len(times), 3))
    latitudes = np.empty(len(times))
    longitudes = np.empty(len(times))
    radii = np.linalg.norm(states[:, :3], axis=1)
    for index, (time_s, state) in enumerate(zip(times, states)):
        body = _body_position(state[:3], float(time_s))
        body_positions[index] = body
        latitudes[index], longitudes[index], _ = _lat_lon_radius_body(body)
    terrain_elevations = np.asarray(
        terrain.elevation_m(latitudes, longitudes), dtype=float
    )
    altitudes = radii - TERRAIN_META.reference_radius_m
    clearances = altitudes - terrain_elevations - C.ball_radius

    def return_distance(time_s: float) -> float:
        state = solution.sol(time_s)
        site_i = rot_z(C.omega * time_s) @ launch_site_body
        return float(np.linalg.norm(state[:3] - site_i))

    eligible = np.flatnonzero(times >= config.minimum_return_time_s)
    candidates: list[tuple[float, float]] = []
    if len(eligible):
        distances = np.array([return_distance(float(times[i])) for i in eligible])
        for j in range(1, len(eligible) - 1):
            if distances[j] <= distances[j - 1] and distances[j] <= distances[j + 1]:
                left = float(times[eligible[j - 1]])
                right = float(times[eligible[j + 1]])
                optimum = minimize_scalar(
                    return_distance,
                    bounds=(left, right),
                    method="bounded",
                    options={"xatol": 1e-8},
                )
                candidates.append((float(optimum.fun), float(optimum.x)))
        if not candidates:
            index = int(np.argmin(distances))
            candidates.append((float(distances[index]), float(times[eligible[index]])))
    if candidates:
        closest_distance, closest_time = min(candidates)
        closest_state = solution.sol(closest_time)
        closest_site_i = rot_z(C.omega * closest_time) @ launch_site_body
        site_velocity_i = np.cross(
            np.array([0.0, 0.0, C.omega]), closest_site_i
        )
        closest_relative_speed = float(
            np.linalg.norm(closest_state[3:] - site_velocity_i)
        )
    else:
        closest_time = closest_distance = closest_relative_speed = None

    impact_time = impact_lat = impact_lon = impact_speed = None
    impacted = bool(solution.t_events and len(solution.t_events[0]))
    if impacted:
        impact_time = float(solution.t_events[0][0])
        impact_state = solution.y_events[0][0]
        impact_body = _body_position(impact_state[:3], impact_time)
        impact_lat, impact_lon, _ = _lat_lon_radius_body(impact_body)
        surface_velocity_i = np.cross(
            np.array([0.0, 0.0, C.omega]), impact_state[:3]
        )
        impact_speed = float(np.linalg.norm(impact_state[3:] - surface_velocity_i))

    final_state = states[-1]
    final_radius = float(np.linalg.norm(final_state[:3]))
    final_energy = 0.5 * float(np.dot(final_state[3:], final_state[3:])) - C.mu / final_radius
    if (
        closest_distance is not None
        and closest_distance + config.model_uncertainty_m
        <= config.return_tolerance_m
        and (impact_time is None or closest_time <= impact_time)
    ):
        status: Status = "RETURN"
        message = (
            "The complete model-uncertainty interval lies inside the requested "
            f"{config.return_tolerance_m:g} m return tolerance."
        )
    elif (
        closest_distance is not None
        and closest_distance - config.model_uncertainty_m
        <= config.return_tolerance_m
        and (impact_time is None or closest_time <= impact_time)
    ):
        status = "RETURN_UNCERTAIN"
        message = (
            "Nominal miss is inside the tolerance plus model-uncertainty band; "
            "the physical return classification is not robust."
        )
    elif impacted:
        status = "IMPACT"
        message = "The ball intersects the LDEM64 surface before a qualified return."
    elif final_energy > 0.0 and final_radius > 2.0 * C.radius:
        status = "ESCAPE"
        message = "The trajectory is energetically unbound and leaves the lunar vicinity."
    else:
        status = "NO_RETURN_WITHIN_WINDOW"
        message = "No qualified return occurs before the propagation window ends."
    minimum_radius = float(np.min(radii))
    harmonic_domain_ok = (
        config.gravity_model.strip().lower() != "grail"
        or minimum_radius >= BRILLOUIN_RADIUS_M
    )
    if not harmonic_domain_ok:
        message += (
            " Part of the arc lies below the conservative Brillouin sphere; "
            "high-degree external harmonics are downward-continued there and "
            "the gravity result is not physically certified."
        )

    trajectory = pd.DataFrame(
        {
            "time_s": times,
            "x_inertial_m": states[:, 0],
            "y_inertial_m": states[:, 1],
            "z_inertial_m": states[:, 2],
            "vx_inertial_m_s": states[:, 3],
            "vy_inertial_m_s": states[:, 4],
            "vz_inertial_m_s": states[:, 5],
            "x_body_m": body_positions[:, 0],
            "y_body_m": body_positions[:, 1],
            "z_body_m": body_positions[:, 2],
            "latitude_deg": latitudes,
            "longitude_deg_east": longitudes,
            "altitude_above_reference_m": altitudes,
            "terrain_elevation_m": terrain_elevations,
            "clearance_m": clearances,
        }
    )
    return TrajectoryResult(
        status,
        message,
        asdict(launch),
        asdict(config),
        terrain_at_launch,
        launch_radius,
        closest_time,
        closest_distance,
        closest_relative_speed,
        impact_time,
        impact_lat,
        impact_lon,
        impact_speed,
        float(np.min(clearances)),
        float(np.max(altitudes)),
        float(final_energy),
        minimum_radius,
        harmonic_domain_ok,
        trajectory,
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--latitude", type=float, required=True)
    parser.add_argument("--longitude", type=float, required=True)
    parser.add_argument("--height", type=float, default=1.0)
    parser.add_argument("--speed", type=float, required=True)
    parser.add_argument("--elevation", type=float, default=0.0)
    parser.add_argument("--azimuth", type=float, default=90.0)
    parser.add_argument(
        "--gravity", choices=("central", "degree2", "grail"), default="grail"
    )
    parser.add_argument("--degree", type=int, default=600)
    parser.add_argument("--duration", type=float, default=8_000.0)
    parser.add_argument("--return-tolerance", type=float, default=1.0)
    parser.add_argument("--model-uncertainty", type=float, default=100.0)
    parser.add_argument("--output", type=Path, default=ROOT / "data" / "output")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    launch = LaunchInput(
        args.latitude,
        args.longitude,
        args.height,
        args.speed,
        args.elevation,
        args.azimuth,
    )
    config = ModelConfig(
        gravity_model=args.gravity,
        gravity_degree=args.degree,
        duration_s=args.duration,
        return_tolerance_m=args.return_tolerance,
        model_uncertainty_m=args.model_uncertainty,
    )
    result = solve_trajectory(launch, config)
    result.save(args.output)
    print(json.dumps(result.summary_dict(), indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
