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


def line_key(d: pd.DataFrame) -> pd.Series:
    """Unique flight-line id. Line names repeat between surveys, so prefix the survey area."""
    return d.survey_area + "|" + d["LINE_NUMB-SEG"]
