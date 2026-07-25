#!/usr/bin/env python3
"""Download and verify the authoritative lunar data used by this project.

Large source products are deliberately kept outside Git.  The script supports
HTTP range-resume, checks the expected byte count, and can verify SHA-256
digests measured from the archived NASA products on 2026-07-25.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
import hashlib
from pathlib import Path
import shutil
import sys
import urllib.request


ROOT = Path(__file__).resolve().parents[1]
EXTERNAL = ROOT / "data" / "external"


@dataclass(frozen=True)
class Product:
    key: str
    url: str
    relative_path: Path
    size: int
    sha256: str
    description: str


PRODUCTS = {
    product.key: product
    for product in (
        Product(
            "ldem64",
            "https://pds-geosciences.wustl.edu/lro/lro-l-lola-3-rdr-v1/"
            "lrolol_1xxx/data/lola_gdr/cylindrical/img/ldem_64.img",
            Path("ldem64/ldem_64.img"),
            530_841_600,
            "98f1824b1a999630bf7b1f59575fe4c1d56a24722a968d98bcce07e6d6d44d4b",
            "LOLA LDEM64, 64 pixels per degree, signed 16-bit, 0.5 m/DN",
        ),
        Product(
            "ldem64-label",
            "https://pds-geosciences.wustl.edu/lro/lro-l-lola-3-rdr-v1/"
            "lrolol_1xxx/data/lola_gdr/cylindrical/img/ldem_64.lbl",
            Path("ldem64/ldem_64.lbl"),
            5_128,
            "4e7a94c76c716a147b21ccf07669d004246b0d1b017b7d63d50befaa1a34142d",
            "PDS3 label for LOLA LDEM64",
        ),
        Product(
            "grgm1200b",
            "https://pgda.gsfc.nasa.gov/data/MoonRM1/sha.grgm1200b_sigma",
            Path("grail/sha.grgm1200b_sigma"),
            83_006_841,
            "f08a988b43f3eaa5a2089045a9b7e41e02f16542c7912b87ea34366fafa39bc5",
            "GRAIL GRGM1200B gravity coefficients and formal uncertainties",
        ),
    )
}


def sha256(path: Path, chunk_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def download(product: Product, verify_hash: bool = True) -> Path:
    destination = EXTERNAL / product.relative_path
    destination.parent.mkdir(parents=True, exist_ok=True)
    partial = destination.with_suffix(destination.suffix + ".part")

    if destination.exists():
        if destination.stat().st_size == product.size:
            if not verify_hash or sha256(destination) == product.sha256:
                print(f"OK      {product.key}: {destination}")
                return destination
        raise RuntimeError(
            f"{destination} exists but does not match the published manifest; "
            "move it aside before downloading again"
        )

    offset = partial.stat().st_size if partial.exists() else 0
    if offset > product.size:
        raise RuntimeError(f"Partial file is larger than expected: {partial}")

    request = urllib.request.Request(product.url)
    if offset:
        request.add_header("Range", f"bytes={offset}-")
    print(f"GET     {product.key}: resuming at {offset:,} / {product.size:,} bytes")

    with urllib.request.urlopen(request, timeout=60) as response:
        status = getattr(response, "status", None)
        if offset and status != 206:
            raise RuntimeError(
                f"Server did not honour Range request for {product.key}; "
                f"remove {partial} and retry"
            )
        mode = "ab" if offset else "wb"
        with partial.open(mode) as output:
            shutil.copyfileobj(response, output, length=8 * 1024 * 1024)

    actual_size = partial.stat().st_size
    if actual_size != product.size:
        raise RuntimeError(
            f"Incomplete {product.key}: {actual_size:,} != {product.size:,} bytes"
        )
    if verify_hash:
        actual_hash = sha256(partial)
        if actual_hash != product.sha256:
            raise RuntimeError(
                f"SHA-256 mismatch for {product.key}: {actual_hash} != "
                f"{product.sha256}"
            )
    partial.replace(destination)
    print(f"SAVED   {product.key}: {destination}")
    return destination


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "products",
        nargs="*",
        choices=sorted(PRODUCTS),
        default=list(PRODUCTS),
        help="products to download (default: all)",
    )
    parser.add_argument(
        "--size-only",
        action="store_true",
        help="skip the SHA-256 pass after a matching-size file is found",
    )
    args = parser.parse_args()
    keys = args.products or list(PRODUCTS)
    for key in keys:
        download(PRODUCTS[key], verify_hash=not args.size_only)


if __name__ == "__main__":
    try:
        main()
    except (OSError, RuntimeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1) from exc
