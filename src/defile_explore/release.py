"""Reading the defile-dataset release tables.

A copy, not a shared module: taken from defile-migration-forecast `src/data/counts.py` (commit
13c5107, 2026-10-09) when Explore moved to its own repo. Only what Explore reads is kept; the
forecast's survey-period processing stays there. The two repos are deliberately independent: a
change to one copy is carried to the other by hand, as a decision (DECISIONS.md -> Repository).

The release (`count.csv`, `survey.csv`, `taxonomy.csv`, `report_text.csv`, `paper_text.csv`,
`README.md`, `datapackage.json`, and the build's `metadata.json`) is copied into
`data/count/dataset/` by `scripts/build_explore.py --dataset <defile-dataset output folder>`. The
build reads the count tables only: the report and paper extracts are copied for the editors of
`content/accounts/` (`scripts/accounts/`), whose written accounts the page shows instead.
"""

import os

import numpy as np
import pandas as pd
from suncalc import get_position

# The site's local time: the dataset's dates are local calendar dates.
TIMEZONE = "Europe/Paris"
# Défilé de l'Écluse (lat, lon), as `LOCATIONS["Defile"]` in the forecast's weather layer.
SITE = (46.117215, 5.914877)

DATASET_DIR = os.path.join("count", "dataset")  # under the data dir
# Release tables, in the release's `dataset/` folder; `metadata.json` sits one level up. The text
# extracts are sources for `content/accounts/`, not read by the build.
DATASET_FILES = (
    "count.csv",
    "survey.csv",
    "taxonomy.csv",
    "report_text.csv",
    "paper_text.csv",
    "datapackage.json",
    "README.md",
)
METADATA_FILE = "metadata.json"

# The main migration direction (`count.csv` `count_category`); `reverse` and `local` are other
# quantities.
MAIN_CATEGORY = "normal"
# Presence without a number (`count_estimation`); its `count` is empty.
PRESENCE_ONLY = "x"
# The dataset keeps an entry timed outside its survey at day level, like an untimed one; only its
# `remark_processing` tells the two apart.
OUTSIDE_SURVEY_REMARK = "Entry time outside the native survey"
# Civil twilight: the sun's altitude (deg) at dawn and dusk.
NIGHT_SUN_ALTITUDE = -6.0


def local_date(t: pd.Series) -> pd.Series:
    """Local calendar date (a naive midnight timestamp) of UTC times."""
    return t.dt.tz_convert(TIMEZONE).dt.tz_localize(None).dt.normalize()


def parse_surveys(survey: pd.DataFrame) -> pd.DataFrame:
    """`survey.csv` plus `start`/`end` (UTC), local `date` and `source`.

    `source` is `trektellen` or `historical` (notebook, spreadsheet, Naturalist, and the curated
    non-counting periods).
    """
    s = survey.copy()
    bounds = s["datetime"].str.split("/", expand=True)
    s["start"] = pd.to_datetime(bounds[0], utc=True)
    s["end"] = pd.to_datetime(bounds[1], utc=True)
    s["date"] = local_date(s["start"])
    s["source"] = np.where(s["recording_era"] == "trektellen", "trektellen", "historical")
    return s


def outside_survey(counts: pd.DataFrame) -> pd.Series:
    """Entries the dataset marks as timed outside their survey."""
    return counts["remark_processing"].fillna("").str.startswith(OUTSIDE_SURVEY_REMARK)


def parse_counts(count: pd.DataFrame, surveys: pd.DataFrame, taxonomy: pd.DataFrame):
    """`count.csv` plus `species`, the `source` and `date` of its survey, and its own timestamp
    `datetime` (UTC; NaT when untimed).

    A count's `datetime` is either empty (it inherits its survey's interval), a local date (untimed,
    or timed outside its survey: the dataset keeps it at day level), or a UTC time. The dataset also
    allows a UTC interval, which no release has used; it is refused rather than guessed at.
    """
    c = count.copy()
    raw = c["datetime"].fillna("")
    interval = raw.str.contains("/")
    if interval.any():
        raise ValueError(
            f"{interval.sum()} count(s) with their own interval (e.g. "
            f"{c.loc[interval, 'count_id'].iloc[0]}): not handled."
        )
    date_only = raw.str.len() == 10
    timed = raw.str.contains("T")
    c["datetime"] = pd.to_datetime(raw.where(timed), utc=True, format="ISO8601")
    c = c.merge(surveys[["survey_id", "source", "date"]], on="survey_id", how="left")
    assert c["source"].notna().all(), "count without a released survey"
    c.loc[timed.values, "date"] = local_date(c.loc[timed.values, "datetime"])
    c.loc[date_only.values, "date"] = pd.to_datetime(raw[date_only]).values
    c["count"] = c["count"].astype(float)  # presence-only counts are empty
    names = taxonomy.set_index("taxon_id")["english_name"]
    c["species"] = c["taxon_id"].map(names)
    return c


def civil_twilight(dates: pd.Series) -> tuple[pd.Series, pd.Series]:
    """Civil dawn and dusk at Defile on each local date, UTC: the first minute of the day with the
    sun at or above -6 deg, and the minute after the last (as defile-dataset computed them)."""
    lat, lon = SITE
    days = pd.DatetimeIndex(pd.to_datetime(dates.unique()))
    midnight = days.tz_localize(TIMEZONE).tz_convert("UTC")
    minutes = np.arange(24 * 60).astype("timedelta64[m]")
    grid = pd.DatetimeIndex((midnight.values[:, None] + minutes[None, :]).ravel())
    altitude = np.degrees(np.asarray(get_position(grid, lon, lat)["altitude"]))
    day = altitude.reshape(len(days), -1) >= NIGHT_SUN_ALTITUDE
    first = day.argmax(axis=1)
    last = day.shape[1] - 1 - day[:, ::-1].argmax(axis=1)
    dawn = dict(zip(days, midnight + pd.to_timedelta(first, unit="min")))
    dusk = dict(zip(days, midnight + pd.to_timedelta(last + 1, unit="min")))
    return dates.map(dawn), dates.map(dusk)
