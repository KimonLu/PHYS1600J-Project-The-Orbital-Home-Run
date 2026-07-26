#!/usr/bin/env python3
"""Run documented reproduction profiles for the project."""

from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]
SCRIPT_DIR = ROOT / "scripts"

CORE_SCRIPTS = [
    "01_ideal_models.py",
    "02_sensitivity.py",
    "03_rotation_and_resonance.py",
    "04_realistic_perturbations.py",
    "05_terrain_envelope.py",
    "06_validation.py",
]

SCIENCE_SCRIPTS = [
    "07_gravity_convergence.py",
    "10_high_fidelity_case.py",
    "11_case_sensitivity.py",
    "14_surface_feasibility.py",
    "16_height_continuation.py",
    "15_high_fidelity_validation.py",
    "13_terrain_visualizations.py",
]

WEB_SCRIPTS = [
    "08_prepare_web_data.py",
    "09_validate_web_solver.py",
    "12_validate_web_terrain.py",
]


def run(script: str) -> None:
    print(f"==> {script}", flush=True)
    subprocess.run(
        [sys.executable, str(SCRIPT_DIR / script)],
        cwd=SCRIPT_DIR,
        check=True,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--profile",
        choices=("core", "science", "web", "full", "validate"),
        default="core",
        help=(
            "core is lightweight; science and web require downloaded NASA "
            "products; full runs every generator; validate checks packaged outputs"
        ),
    )
    args = parser.parse_args()

    if args.profile in {"core", "science", "full"}:
        for script in CORE_SCRIPTS:
            run(script)
    if args.profile in {"science", "full"}:
        for script in SCIENCE_SCRIPTS:
            run(script)
    if args.profile in {"web", "full"}:
        for script in WEB_SCRIPTS:
            run(script)
    if args.profile in {"science", "web", "full", "validate"}:
        run("validate_results.py")


if __name__ == "__main__":
    main()
