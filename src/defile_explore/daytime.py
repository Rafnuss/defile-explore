"""Passage through the day, as counted: birds per counted hour by solar hour, on the days timed to
the hour, with the smooth time-of-day profile (`defile_explore.profile`) beside it.

The rate of an hour is the birds timed to it over the hours counted in it, pooled over the days of
a period, so hours counted less often weigh less, not zero. Only days with birds and at least
`PROFILE_MIN_TIMED` of them timed are used (the profile's own sample), from whenever the counts are
timed (2014 on). Shares are of the period's summed hourly rates.

The block gives one histogram for the whole season (`hours`), and the date x hour matrix only where
the passage hour changes through the season (`change`): the late part's mean passage hour minus the
early part's (`PART_QUANTILES` of the passage), with a 95% interval from resampling days, so a
single flock cannot make a change look certain. `change.show` when the interval excludes zero and
the shift is at least `SHIFT_MIN_HOURS`. In solar time (`defile_explore.export.solar_shift`), like
every hour of the export: the end of summer time and the equation of time no longer make late taxa
look earlier, so a change shown is the birds' own.
"""

import numpy as np
import pandas as pd

from defile_explore.export import HOURS
from defile_explore.profile import PROFILE_DOY, PROFILE_MIN_TIMED

METHOD = "daytime@3"  # 3: solar hours
PENTAD = 5  # days per period of the date x hour grid
MIN_HOUR_EFFORT = 3.0  # counted hours below which an hour's rate in a period is not shown
MIN_PERIOD_BIRDS = 30  # timed birds below which a period's column is not shown
PART_QUANTILES = (0.1, 0.35, 0.65, 0.9)  # early / peak / late: between these passage shares
MIN_PART_DAYS = 5  # timed days in the early and the late part for the change to be tested
SHIFT_MIN_HOURS = 1.0  # a change in the mean passage hour smaller than this is not shown
BOOTSTRAP = 400  # resamples of days for the change's interval
SEED = 0


def timed_days(days: pd.DataFrame, effort: pd.DataFrame) -> pd.DataFrame:
    """The counted days with birds, mostly timed to the hour: `date` and `hourly` (hours counted in
    each solar hour)."""
    d = days.loc[(days["count"] > 0) & (days["timed"] >= PROFILE_MIN_TIMED), "date"]
    counted = effort[(effort["state"] == "counted") & effort["date"].isin(d)]
    return counted[["date", "hourly"]]


def hour_rates(hourly: pd.DataFrame, sample: pd.DataFrame, key: pd.Series) -> tuple:
    """`(birds, hours)`, both indexed by `key`'s values (one per sample day) x 24 hours."""
    effort = pd.DataFrame(np.stack(sample["hourly"].to_numpy()), index=sample.index)
    hours = effort.groupby(key.to_numpy()).sum()
    h = hourly[hourly["date"].isin(sample["date"])]
    k = h["date"].map(pd.Series(key.to_numpy(), index=sample["date"].to_numpy()))
    birds = (
        h.assign(k=k.to_numpy())
        .pivot_table(index="k", columns="hour", values="count", aggfunc="sum")
        .reindex(index=hours.index, columns=range(HOURS))
        .fillna(0)
    )
    return birds, hours


def pooled_share(birds: np.ndarray, hours: np.ndarray) -> np.ndarray | None:
    """Share per hour of the pooled rates `birds / hours` (hours under `MIN_HOUR_EFFORT`: 0)."""
    r = np.where(hours >= MIN_HOUR_EFFORT, birds / np.maximum(hours, 1e-9), 0.0)
    return r / r.sum() if r.sum() > 0 else None


def mean_hour(birds: np.ndarray, hours: np.ndarray) -> float:
    """Mean passage hour (the middle of the solar hour) of the pooled rates."""
    p = pooled_share(birds, hours)
    return float((p * (np.arange(HOURS) + 0.5)).sum()) if p is not None else np.nan


def seasonal_change(birds: np.ndarray, hours: np.ndarray, doy: np.ndarray, parts_doy) -> dict:
    """Late minus early mean passage hour with a 95% interval from resampling days within each
    part; `birds`, `hours`: one row per timed day."""
    early = np.flatnonzero((doy >= parts_doy[0]) & (doy < parts_doy[1]))
    late = np.flatnonzero((doy >= parts_doy[2]) & (doy < parts_doy[3]))
    out = {"days": [len(early), len(late)], "show": False}
    if min(len(early), len(late)) < MIN_PART_DAYS:
        return out
    m_early = mean_hour(birds[early].sum(0), hours[early].sum(0))
    m_late = mean_hour(birds[late].sum(0), hours[late].sum(0))
    rng = np.random.default_rng(SEED)
    shifts = []
    for _ in range(BOOTSTRAP):
        a, b = rng.choice(early, len(early)), rng.choice(late, len(late))
        shifts.append(
            mean_hour(birds[b].sum(0), hours[b].sum(0))
            - mean_hour(birds[a].sum(0), hours[a].sum(0))
        )
    lo, hi = np.nanpercentile(shifts, [2.5, 97.5])
    shift = m_late - m_early
    out |= {"early": m_early, "late": m_late, "shift": shift, "lo": float(lo), "hi": float(hi)}
    out["show"] = bool((lo > 0 or hi < 0) and abs(shift) >= SHIFT_MIN_HOURS)
    return out


def daytime_block(
    days: pd.DataFrame,
    hourly: pd.DataFrame,
    effort: pd.DataFrame,
    profile: np.ndarray,
    parts_doy: list[float] | None,
) -> dict | None:
    """The `daytime` block: the counted share per hour over the whole season (`hours`) beside the
    smooth profile averaged over the main passage (`parts_doy`: the days of year at
    `PART_QUANTILES` of the passage), the date x hour grid of shares, and whether it changes enough
    to be shown (`change`).

    None without timed days.
    """
    sample = timed_days(days, effort)
    if sample.empty:
        return None
    doy = sample["date"].dt.dayofyear
    per_day, day_hours = hour_rates(hourly, sample, pd.Series(np.arange(len(sample))))
    b, h = per_day.to_numpy(), day_hours.to_numpy()
    birds, hours = hour_rates(hourly, sample, (doy // PENTAD) * PENTAD)
    rate = birds / hours.where(hours >= MIN_HOUR_EFFORT)
    share = rate.div(rate.sum(axis=1), axis=0)[birds.sum(axis=1) >= MIN_PERIOD_BIRDS]
    total = pooled_share(b.sum(0), h.sum(0))
    out = {
        "method": METHOD,
        "years": [int(sample["date"].dt.year.min()), int(sample["date"].dt.year.max())],
        "days": len(sample),
        "birds": float(b.sum()),
        "hours": total,
        "profile": None,
        "doy": share.index + PENTAD / 2,
        "share": share.to_numpy(),
        "change": {"days": [0, 0], "show": False},
    }
    if parts_doy is None:
        return out
    rows = np.clip(np.arange(round(parts_doy[0]), round(parts_doy[-1]) + 1), *PROFILE_DOY)
    out["profile"] = profile[rows - PROFILE_DOY[0]].mean(axis=0)
    out["change"] = seasonal_change(b, h, doy.to_numpy(), parts_doy)
    return out
