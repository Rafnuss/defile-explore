"""Passage through the day, as counted: birds per counted hour by solar hour, on the days timed to
the hour, with the smooth time-of-day profile (`defile_explore.profile`) beside it.

The rate of an hour is the birds timed to it over the hours counted in it, pooled over the days of
a period, so hours counted less often weigh less, not zero. Only days with birds and at least
`PROFILE_MIN_TIMED` of them timed are used, and on each only the hours counted for at least
`PROFILE_MIN_HOUR_EFFORT` (the profile's own sample): a block's total entered at its last minute
cannot make a bar of an hour barely counted. Shares are of the period's summed hourly rates.

The block gives one histogram over the main passage (`hours`, the days between the first and last
`PART_QUANTILES`, else the whole season) beside what the profile predicts for the same days and
minutes counted (`expected`: each day's birds spread over its counted hours by the profile, pooled
the same way), so bars and line answer the same question; and the date x hour matrix only where the
passage hour changes through the season (`change`): the late part's mean passage hour minus the
early part's (`PART_QUANTILES` of the passage), with a 95% interval from resampling days, so a
single flock cannot make a change look certain. `change.show` when the interval excludes zero and
the shift is at least `SHIFT_MIN_HOURS`. In solar time (`defile_explore.export.solar_shift`), like
every hour of the export: the end of summer time and the equation of time no longer make late taxa
look earlier, so a change shown is the birds' own.
"""

import numpy as np
import pandas as pd

from defile_explore.export import HOURS, by_hour
from defile_explore.profile import PROFILE_DOY, PROFILE_MIN_HOUR_EFFORT, PROFILE_MIN_TIMED

METHOD = "daytime@4"  # 4: main passage, hours counted at least half, profile on the same days
PENTAD = 5  # days per period of the date x hour grid
MIN_HOUR_EFFORT = 3.0  # counted hours below which an hour's rate in a period is not shown
MIN_PERIOD_BIRDS = 30  # timed birds below which a period's column is not shown
PART_QUANTILES = (0.1, 0.35, 0.65, 0.9)  # early / peak / late: between these passage shares
MIN_PART_DAYS = 5  # timed days in the early and the late part for the change to be tested
SHIFT_MIN_HOURS = 1.0  # a change in the mean passage hour smaller than this is not shown
BOOTSTRAP = 400  # resamples of days for the change's interval
SEED = 0


def timed_days(days: pd.DataFrame, effort: pd.DataFrame) -> pd.DataFrame:
    """The counted days with birds, mostly timed to the hour: `date`, `hourly` (hours counted in
    each solar hour) and `slots`."""
    d = days.loc[(days["count"] > 0) & (days["timed"] >= PROFILE_MIN_TIMED), "date"]
    counted = effort[(effort["state"] == "counted") & effort["date"].isin(d)]
    return counted[["date", "hourly", "slots"]]


def day_hours(hourly: pd.DataFrame, sample: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    """`(birds, hours)`, one row per sample day x 24 hours, both zero in the hours counted for less
    than `PROFILE_MIN_HOUR_EFFORT`."""
    hours = np.stack(sample["hourly"].to_numpy())
    h = hourly[hourly["date"].isin(sample["date"])]
    row = h["date"].map(pd.Series(np.arange(len(sample)), index=sample["date"].to_numpy()))
    birds = np.zeros_like(hours, dtype=float)
    np.add.at(birds, (row.to_numpy(), h["hour"].to_numpy()), h["count"].to_numpy())
    used = hours >= PROFILE_MIN_HOUR_EFFORT
    return np.where(used, birds, 0.0), np.where(used, hours, 0.0)


def hour_rates(hourly: pd.DataFrame, sample: pd.DataFrame, key: pd.Series) -> tuple:
    """`(birds, hours)` (`day_hours`) summed by `key`'s values (one per sample day) x 24 hours."""
    birds, hours = day_hours(hourly, sample)
    k = key.to_numpy()
    return pd.DataFrame(birds).groupby(k).sum(), pd.DataFrame(hours).groupby(k).sum()


def expected_birds(sample: pd.DataFrame, birds: np.ndarray, profile: np.ndarray) -> np.ndarray:
    """What the profile (on `PROFILE_DOY`) predicts per hour on each sample day: the day's birds
    (`birds`, `day_hours`) spread over the hours it uses in proportion to the profile in the
    minutes counted."""
    doy = sample["date"].dt.dayofyear.clip(*PROFILE_DOY).to_numpy() - PROFILE_DOY[0]
    used = np.stack(sample["hourly"].to_numpy()) >= PROFILE_MIN_HOUR_EFFORT
    e = np.where(used, by_hour(profile[doy] * np.stack(sample["slots"].to_numpy())), 0.0)
    total = e.sum(axis=1, keepdims=True)
    share = np.divide(e, total, out=np.zeros_like(e), where=total > 0)
    return share * birds.sum(axis=1, keepdims=True)


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
    """The `daytime` block: the counted share per hour over the main passage (`hours`; `parts_doy`:
    the days of year at `PART_QUANTILES` of the passage, else the whole season) beside the
    profile's prediction for the same days and minutes (`expected`), the date x hour grid of
    shares, and whether it changes enough to be shown (`change`).

    None without timed days.
    """
    sample = timed_days(days, effort)
    if sample.empty:
        return None
    doy = sample["date"].dt.dayofyear
    b, h = day_hours(hourly, sample)
    birds, hours = hour_rates(hourly, sample, (doy // PENTAD) * PENTAD)
    rate = birds / hours.where(hours >= MIN_HOUR_EFFORT)
    share = rate.div(rate.sum(axis=1), axis=0)[birds.sum(axis=1) >= MIN_PERIOD_BIRDS]
    main = np.ones(len(sample), bool)
    if parts_doy is not None:
        main = ((doy >= parts_doy[0]) & (doy <= parts_doy[-1])).to_numpy()
    x = expected_birds(sample, b, profile)
    out = {
        "method": METHOD,
        "years": [int(sample["date"].dt.year.min()), int(sample["date"].dt.year.max())],
        "days": len(sample),
        "birds": float(b.sum()),
        "main_days": int(main.sum()),
        "main_doy": None if parts_doy is None else [float(parts_doy[0]), float(parts_doy[-1])],
        "hours": pooled_share(b[main].sum(0), h[main].sum(0)),
        "expected": pooled_share(x[main].sum(0), h[main].sum(0)),
        "doy": share.index + PENTAD / 2,
        "share": share.to_numpy(),
        "change": {"days": [0, 0], "show": False},
    }
    if parts_doy is not None:
        out["change"] = seasonal_change(b, h, doy.to_numpy(), parts_doy)
    return out
