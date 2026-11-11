"""Shared helpers for the BGS data scripts (imported, not run directly)."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RAW_CSV = ROOT / "data" / "raw" / "bgs_aeromag.csv"


def load_points(path: Path = RAW_CSV) -> pd.DataFrame:
    if not path.exists():
        raise SystemExit(f"{path} not found: run `uv run scripts/fetch_bgs_aeromag.py` first")
    return pd.read_csv(path)


def repair_heights(d: pd.DataFrame) -> pd.DataFrame:
    """Set `aod = ground_elev + mtc` where `aod` was recorded as 0.

    449 HG64 points over the Highlands (ground 695-909 m) have aod = 0 although
    ground_elev + mtc is 1000-1214 m; every other point satisfies aod == ground_elev + mtc.
    Fitting them at sea level would place real data about 1 km too low.
    """
    bad = (d.aod == 0) & (d.ground_elev + d.mtc > 0)
    if bad.any():
        d = d.copy()
        d.loc[bad, "aod"] = d.ground_elev[bad] + d.mtc[bad]
        print(f"repaired aod = 0 at {int(bad.sum())} points (set to ground_elev + mtc)")
    return d


def line_key(d: pd.DataFrame) -> pd.Series:
    """Unique flight-line id. Line names repeat between surveys, so prefix the survey area."""
    return d.survey_area + "|" + d["LINE_NUMB-SEG"]
