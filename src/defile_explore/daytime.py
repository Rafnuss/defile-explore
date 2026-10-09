"""Passage through the day, as counted: birds per counted hour by local clock hour, on the days
timed to the hour, with the smooth time-of-day profile (`defile_explore.profile`) beside it.

The rate of an hour is the birds timed to it over the hours counted in it, pooled over the days of
a period, so hours counted less often weigh less, not zero. Only days with birds and at least
`PROFILE_MIN_TIMED` of them timed are used (the profile's own sample), from whenever the counts are
timed (2014 on). Shares are of the period's summed hourly rates.
"""

import numpy as np
import pandas as pd

from defile_explore.export import HOURS
from defile_explore.profile import PROFILE_DOY, PROFILE_MIN_TIMED

METHOD = "daytime@1"
PENTAD = 5  # days per period of the date x hour grid
MIN_HOUR_EFFORT = 3.0  # counted hours below which an hour's rate in a period is not shown
MIN_PERIOD_BIRDS = 30  # timed birds below which a period's column is not shown
PART_QUANTILES = (0.1, 0.35, 0.65, 0.9)  # early / peak / late: between these passage shares
PART_LABELS = ("early", "peak", "late")


def timed_days(days: pd.DataFrame, effort: pd.DataFrame) -> pd.DataFrame:
    """The counted days with birds, mostly timed to the hour: `date` and `hourly` (hours counted in
    each local clock hour)."""
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


def daytime_block(
    days: pd.DataFrame,
    hourly: pd.DataFrame,
    effort: pd.DataFrame,
    profile: np.ndarray,
    parts_doy: list[float] | None,
) -> dict | None:
    """The `daytime` block: the date x hour grid of shares and, for the early, peak and late season
    (`parts_doy`: the days of year at `PART_QUANTILES` of the passage), the counted shares per hour
    beside the profile averaged over those days.

    None without timed days.
    """
    sample = timed_days(days, effort)
    if sample.empty:
        return None
    doy = sample["date"].dt.dayofyear
    birds, hours = hour_rates(hourly, sample, (doy // PENTAD) * PENTAD)
    rate = birds / hours.where(hours >= MIN_HOUR_EFFORT)
    share = rate.div(rate.sum(axis=1), axis=0)[birds.sum(axis=1) >= MIN_PERIOD_BIRDS]
    out = {
        "method": METHOD,
        "years": [int(sample["date"].dt.year.min()), int(sample["date"].dt.year.max())],
        "days": len(sample),
        "birds": float(birds.to_numpy().sum()),
        "doy": share.index + PENTAD / 2,
        "share": share.to_numpy(),
        "parts": [],
    }
    if parts_doy is None:
        return out
    for label, lo, hi in zip(PART_LABELS, parts_doy[:-1], parts_doy[1:]):
        sel = (doy >= lo) & (doy < hi)
        if not sel.any():
            continue
        b, h = hour_rates(hourly, sample[sel.to_numpy()], pd.Series(0, index=sample.index[sel]))
        r = (b.iloc[0] / h.iloc[0].where(h.iloc[0] >= MIN_HOUR_EFFORT)).fillna(0)
        rows = np.clip(np.arange(round(lo), round(hi)), *PROFILE_DOY) - PROFILE_DOY[0]
        out["parts"].append(
            {
                "label": label,
                "doy": [lo, hi],
                "birds": float(b.to_numpy().sum()),
                "counted": (r / r.sum()).to_numpy() if r.sum() > 0 else None,
                "profile": profile[rows].mean(axis=0),
            }
        )
    return out
