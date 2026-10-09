"""The season as counted: every year's passage day by day, its passage dates, and the chance of
seeing the taxon on a given date. Empirical: counts and coverage only, no model.

A day's rate is its birds per full counted day, `y / c` (`defile_explore.profile.coverage`), on
days counted for at least `COVERAGE_MIN` of the expected passage; other days are blank. A year's
share of a day divides that rate by the year's total, in which blank days are filled by linear
interpolation between counted days (so the total, not the shares shown, assumes the passage between
two counted days ran between their rates). The year's passage dates (10/50/90%) come from the same
filled series; `counted` is the share of the window's days actually counted, so a year with long
gaps can be shown as less certain.
"""

import numpy as np
import pandas as pd

from defile_explore.profile import COVERAGE_MIN

METHOD = "season@1"
PASSAGE_QUANTILES = (0.1, 0.5, 0.9)  # as the smooth season's (`trend.PASSAGE_QUANTILES`)
CHANCE_YEARS = 10  # chances are for the last this many complete seasons
CHANCE_THRESHOLDS = (1, 10, 100, 1000)  # at least this many birds in the day
CHANCE_MIN_SHARE = 0.02  # a threshold reached on fewer days than this anywhere is left out
PENTAD = 5  # days per period for the chances


def daily_rates(frame: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """`(rate, filled)`: birds per full counted day by year (rows) and day of year (columns), NaN
    where not counted (`c < COVERAGE_MIN`), and the same with gaps filled by linear interpolation
    (edges held at the nearest counted day).

    `frame`: `defile_explore.trend.model_frame` (one row per window day: `year`, `doy`, `y`, `c`).
    """
    f = frame.assign(rate=(frame["y"] / frame["c"]).where(frame["c"] >= COVERAGE_MIN))
    rate = f.pivot(index="year", columns="doy", values="rate")
    filled = rate.interpolate(axis=1, limit_direction="both").fillna(0)
    return rate, filled


def cumulative_quantiles(doy: np.ndarray, rate: np.ndarray, qs=PASSAGE_QUANTILES) -> list[float]:
    """Days of year by which shares `qs` of `rate`'s total have passed, interpolated between days;
    NaN for an empty series."""
    total = rate.sum()
    if total <= 0:
        return [np.nan] * len(qs)
    cum = np.cumsum(rate) / total
    return [float(np.interp(q, cum, doy)) for q in qs]


def year_passage(rate: pd.DataFrame, filled: pd.DataFrame) -> pd.DataFrame:
    """Per year, the 10/50/90% passage dates of the filled series and `counted`, the share of the
    window's days counted."""
    doy = filled.columns.to_numpy()
    rows = [
        [y, *cumulative_quantiles(doy, filled.loc[y].to_numpy()), rate.loc[y].notna().mean()]
        for y in filled.index
    ]
    cols = ["year", *(f"q{round(q * 100)}" for q in PASSAGE_QUANTILES), "counted"]
    return pd.DataFrame(rows, columns=cols).round(2)


def chances(frame: pd.DataFrame, last_year: int) -> dict:
    """Per period of `PENTAD` days, the share of well-counted days (`c >= COVERAGE_MIN`) of the
    last `CHANCE_YEARS` seasons with at least each of `CHANCE_THRESHOLDS` birds."""
    f = frame[(frame["year"] > last_year - CHANCE_YEARS) & (frame["c"] >= COVERAGE_MIN)]
    f = f.assign(period=(f["doy"] // PENTAD) * PENTAD)
    g = f.groupby("period")["y"]
    out = {"years": [last_year - CHANCE_YEARS + 1, last_year], "doy": g.size().index + PENTAD / 2}
    out["days"] = g.size().to_numpy()
    for k in CHANCE_THRESHOLDS:
        share = g.apply(lambda s, k=k: (s >= k).mean())
        if share.max() >= CHANCE_MIN_SHARE:
            out[f"at_least_{k}"] = share.round(3).to_numpy()
    return out


def season_block(frame: pd.DataFrame, last_year: int) -> dict:
    """The `season` block of a species file: `years` x `doy` shares of the year's birds (null =
    not counted), each year's passage dates, and the chances."""
    rate, filled = daily_rates(frame)
    share = rate.div(filled.sum(axis=1).where(lambda t: t > 0), axis=0)
    return {
        "method": METHOD,
        "years": rate.index.to_numpy(),
        "doy": rate.columns.to_numpy(),
        "share": share.to_numpy(),
        "passage": year_passage(rate, filled),
        "chances": chances(frame, last_year),
    }
