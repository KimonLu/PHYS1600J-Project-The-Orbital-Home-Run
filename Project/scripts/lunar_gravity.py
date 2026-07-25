"""High-degree lunar gravity backed by NASA GRAIL GRGM1200B."""
from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path

import numpy as np

from orbital_home_run import C, rot_z


ROOT = Path(__file__).resolve().parents[1]


@dataclass(frozen=True)
class GravityMetadata:
    model: str = "GRGM1200B"
    maximum_degree: int = 1200
    gm_m3_s2: float = 4.9028001224452998e12
    reference_radius_m: float = 1_738_000.0
    normalization: str = "4pi geodesy normalized"
    frame: str = "lunar principal axes"


META = GravityMetadata()
# Conservative sphere enclosing the highest LDEM64 terrain plus a small
# frame/product margin. An external spherical-harmonic expansion is formally
# guaranteed to converge only outside a mass-enclosing sphere.
BRILLOUIN_RADIUS_M = 1_748_200.0


class GRAILGravity:
    """Point acceleration evaluator for a selectable harmonic truncation.

    PySHTOOLS returns spherical components ``(radial, colatitudinal, east)``.
    They are converted to Cartesian in the rotating body frame and then into
    the inertial frame used by the propagator.
    """

    def __init__(
        self,
        path: Path | str | None = None,
        maximum_degree: int = META.maximum_degree,
    ) -> None:
        try:
            import pyshtools as pysh
            from pyshtools.backends import shtools
        except ImportError as exc:
            raise ImportError(
                "High-degree gravity requires pyshtools. Install requirements.txt."
            ) from exc
        if not 0 <= maximum_degree <= META.maximum_degree:
            raise ValueError("maximum_degree must be between 0 and 1200")
        if path is None:
            path = ROOT / "data" / "external" / "grail" / "sha.grgm1200b_sigma"
        self.path = Path(path)
        if not self.path.exists():
            raise FileNotFoundError(
                f"GRGM1200B not found at {self.path}. Run download_science_data.py."
            )
        self.maximum_degree = int(maximum_degree)
        self._make_gravity_point = shtools.MakeGravGridPoint
        self.coefficients = pysh.SHGravCoeffs.from_file(
            self.path,
            format="shtools",
            header=True,
            r0_index=1,
            gm_index=0,
            lmax=self.maximum_degree,
            errors=True,
            normalization="4pi",
            csphase=1,
        )

    def body_acceleration(
        self, position_body_m: np.ndarray, degree: int | None = None
    ) -> np.ndarray:
        r = np.asarray(position_body_m, dtype=float)
        radius = float(np.linalg.norm(r))
        if radius <= META.reference_radius_m * 0.95:
            raise ValueError("GRGM point is unphysically far inside the Moon")
        lat = math.degrees(math.asin(float(r[2]) / radius))
        lon = math.degrees(math.atan2(float(r[1]), float(r[0]))) % 360.0
        lmax = self.maximum_degree if degree is None else int(degree)
        if not 0 <= lmax <= self.maximum_degree:
            raise ValueError("degree exceeds loaded coefficient limit")
        gr, gtheta, gphi = self._make_gravity_point(
            self.coefficients.coeffs,
            self.coefficients.gm,
            self.coefficients.r0,
            radius,
            lat,
            lon,
            lmax=lmax,
        )
        lat_r = math.radians(lat)
        lon_r = math.radians(lon)
        cl, sl = math.cos(lat_r), math.sin(lat_r)
        co, so = math.cos(lon_r), math.sin(lon_r)
        radial = np.array([cl * co, cl * so, sl])
        north = np.array([-sl * co, -sl * so, cl])
        east = np.array([-so, co, 0.0])
        return float(gr) * radial - float(gtheta) * north + float(gphi) * east

    @staticmethod
    def inside_brillouin_sphere(position_m: np.ndarray) -> bool:
        return float(np.linalg.norm(position_m)) < BRILLOUIN_RADIUS_M

    def inertial_acceleration(
        self, time_s: float, position_inertial_m: np.ndarray, degree: int | None = None
    ) -> np.ndarray:
        body_from_inertial = rot_z(-C.omega * time_s)
        position_body = body_from_inertial @ np.asarray(position_inertial_m)
        acceleration_body = self.body_acceleration(position_body, degree=degree)
        return body_from_inertial.T @ acceleration_body
