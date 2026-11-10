# /// script
# requires-python = ">=3.11"
# dependencies = []
# ///
"""Download the BGS Great Britain Aeromagnetic Survey points as one CSV.

Source: https://ogcapi.bgs.ac.uk/collections/aeromag (OGC API Features, no login).
Output (gitignored, never committed while the licence is unconfirmed):

    data/raw/bgs_aeromag.csv     all points, ordered by feature_id
    data/raw/manifest.json       URL, date, numberMatched, row count, sha256

Pages are cached in data/raw/pages/ so an interrupted download resumes where it stopped.

    uv run scripts/fetch_bgs_aeromag.py
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import io
import json
import sys
import time
import urllib.request
from pathlib import Path

BASE_URL = "https://ogcapi.bgs.ac.uk/collections/aeromag/items"
PAGE_SIZE = 10_000
ROOT = Path(__file__).resolve().parents[1]


def get(url: str, retries: int = 5) -> bytes:
    """GET with exponential back-off; the API is labelled BETA and does drop requests."""
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=120) as resp:
                return resp.read()
        except OSError as err:  # URLError, HTTPError and timeouts are all OSError
            if attempt == retries - 1:
                raise
            wait = 2**attempt
            print(f"  {err}; retrying in {wait} s", file=sys.stderr)
            time.sleep(wait)
    raise AssertionError("unreachable")


def number_matched() -> int:
    body = json.loads(get(f"{BASE_URL}?f=json&limit=1"))
    return int(body["numberMatched"])


def fetch_page(offset: int, pages_dir: Path) -> Path:
    path = pages_dir / f"page_{offset:07d}.csv"
    if path.exists():
        return path
    url = f"{BASE_URL}?f=csv&limit={PAGE_SIZE}&offset={offset}&sortby=feature_id"
    data = get(url)
    # Write to a temporary name first so a crash never leaves a truncated page behind.
    tmp = path.with_suffix(".part")
    tmp.write_bytes(data)
    tmp.replace(path)
    return path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out-dir", type=Path, default=ROOT / "data" / "raw")
    args = parser.parse_args()

    pages_dir = args.out_dir / "pages"
    pages_dir.mkdir(parents=True, exist_ok=True)

    expected = number_matched()
    print(f"numberMatched = {expected}")
    offsets = range(0, expected, PAGE_SIZE)

    header: list[str] | None = None
    rows: list[list[str]] = []
    for i, offset in enumerate(offsets, start=1):
        page = fetch_page(offset, pages_dir)
        reader = csv.reader(io.StringIO(page.read_text(encoding="utf-8")))
        page_header = next(reader)
        if header is None:
            header = page_header
        elif page_header != header:
            raise SystemExit(f"column mismatch in {page.name}: {page_header}")
        rows.extend(reader)
        print(f"  page {i}/{len(offsets)}: {len(rows)} rows", end="\r")
    print()
    assert header is not None

    # Pass criteria from the plan: every point arrived exactly once.
    ids = [int(r[header.index("feature_id")]) for r in rows]
    if len(rows) != expected:
        raise SystemExit(f"FAIL: {len(rows)} rows, expected {expected}")
    if len(set(ids)) != len(ids):
        raise SystemExit(f"FAIL: {len(ids) - len(set(ids))} duplicate feature_id values")
    print(f"PASS: {len(rows)} rows == numberMatched, no duplicate feature_id")

    out_csv = args.out_dir / "bgs_aeromag.csv"
    with out_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(rows)

    manifest = {
        "source": BASE_URL,
        "downloaded_utc": dt.datetime.now(dt.UTC).isoformat(timespec="seconds"),
        "number_matched": expected,
        "rows": len(rows),
        "columns": header,
        "sha256": hashlib.sha256(out_csv.read_bytes()).hexdigest(),
        "licence": "unconfirmed (see docs/research-notes.md); do not redistribute",
    }
    (args.out_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"wrote {out_csv} and manifest.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
