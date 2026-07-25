#!/usr/bin/env python3
from __future__ import annotations
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = [
    "01_ideal_models.py",
    "02_sensitivity.py",
    "03_rotation_and_resonance.py",
    "04_realistic_perturbations.py",
    "05_terrain_analysis.py",
    "06_validation.py",
    "validate_results.py",
]


def main() -> None:
    for script in SCRIPTS:
        print(f"==> {script}")
        subprocess.run([sys.executable, str(ROOT/"scripts"/script)], cwd=ROOT/"scripts", check=True)


if __name__ == "__main__":
    main()
