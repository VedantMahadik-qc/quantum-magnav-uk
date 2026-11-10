# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy"]
# ///
"""Pick start points for `magnav demo --map-file`: the strongest- and weakest-gradient flights.

The demo flies 60 km east from its start. A candidate start (on a 5 km lattice) qualifies if
every grid node within the corridor [x-5 km, x+65 km] x [y-5 km, y+5 km] has data, so no
particle can leave the map. Among those, print the starts with the largest and smallest RMS
horizontal gradient (nT/km) over the corridor.

    uv run scripts/pick_demo_lines.py [data/grid/bgs_250m.npy]
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    npy = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "data" / "grid" / "bgs_250m.npy"
    values = np.load(npy)
    meta = json.loads(npy.with_suffix(".json").read_text())
    h = meta["dx"]
    gy, gx = np.gradient(values, h)  # nT/m; NaN next to gaps, which also disqualifies
    grad2 = (gx**2 + gy**2) * 1e6  # (nT/km)^2

    def idx(v: float) -> int:
        return int(round(v / h))

    results = []
    ny, nx = values.shape
    step = idx(5_000)
    for row in range(idx(5_000), ny - idx(5_000), step):
        for col in range(idx(5_000), nx - idx(65_000), step):
            win = grad2[row - idx(5_000):row + idx(5_000) + 1,
                        col - idx(5_000):col + idx(65_000) + 1]
            if np.isnan(win).any():
                continue
            x = meta["x0"] + col * h
            y = meta["y0"] + row * h
            results.append((float(np.sqrt(win.mean())), x, y))

    results.sort()
    print(f"{len(results)} fully covered corridors")
    for label, (g, x, y) in [("weakest", results[0]), ("strongest", results[-1])]:
        print(f"{label:9s}: --start-x {x:.0f} --start-y {y:.0f}   RMS gradient {g:.1f} nT/km")


if __name__ == "__main__":
    main()
