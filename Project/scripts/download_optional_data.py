#!/usr/bin/env python3
"""Download optional LOLA/GRAIL products listed in the project manifest.

The default paper figures do not require these large products. This script is
provided to extend the terrain and high-degree gravity analyses with official
NASA PDS data.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path
import urllib.request

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "data" / "input" / "external_data_manifest.csv"
DEST = ROOT / "data" / "input"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("names", nargs="*", help="Dataset names; omit for all")
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    selected = set(args.names)
    with MANIFEST.open(newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    for row in rows:
        if selected and row["dataset"] not in selected:
            continue
        target = DEST / row["local_filename"]
        if target.exists() and not args.overwrite:
            print(f"skip {target.name}: already exists")
            continue
        print(f"download {row['dataset']} -> {target}")
        urllib.request.urlretrieve(row["download_url"], target)


if __name__ == "__main__":
    main()
