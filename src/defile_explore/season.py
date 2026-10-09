"""The season day by day: every year's passage, its passage dates, and the chance of seeing the
taxon on a given date. One series per taxon, from its trend fit where it has one.

With a trend (`source: gam`), every day of the model window has the trend GAM's gap-filled full-day
birds (`trend.fill_draws`, the mean of the draws, as `trend.days`): the birds counted plus the
hours and days nobody counted. A year's share of a day divides it by the year's total over the
model window; the year's passage dates (10/50/90%) are taken in each draw, so each comes with the
median over draws and the median date with its 80% band (`passage_from_draws`); the chances are the
probability, over the draws, that a full day holds at least each threshold (`chance_from_draws`).
On the gap transplant this beats interpolation day by day, on the passage dates and on the totals,
the more so the gappier the year (`DECISIONS.md`). Where the trend's totals are hidden
(`reliability.fills_season`), or without a trend, the season is from the counts.

From the counts (`source: counts`): a day's rate is its birds per full counted day, `y / c`
(`defile_explore.profile.coverage`), on days counted for at least `COVERAGE_MIN` of the expected
passage; other days are blank. A year's total fills blank days by linear interpolation between
counted days, and its passage dates come from that filled series.

Either way, days are shown over the taxon's view window but totals, shares and passage dates are of
its model window (`defile_explore.window`); a view day beyond the model window has its rate from
the counts (`y / c`, blank below `COVERAGE_MIN`). `count` holds the birds counted on every day
counted at all, `c` its coverage (0 on a day not counted, or counted below `trend.MIN_COVERAGE`:
its birds are in `count`), and `counted` per year the share of the window's days counted.
"""

import numpy as np
import pandas as pd

from defile_explore.profile import COVERAGE_MIN

METHOD = "season@3"
PASSAGE_QUANTILES = (0.1, 0.5, 0.9)  # as the smooth season's (`trend.PASSAGE_QUANTILES`)
CHANCE_YEARS = 10  # chances are for the last this many complete seasons
CHANCE_THRESHOLDS = (1, 10, 100, 1000)  # at least this many birds in the day
CHANCE_MIN_SHARE = 0.02  # a threshold reached on fewer days than this anywhere is left out
PENTAD = 5  # days per period for the chances
PASSAGE_BAND = (0.1, 0.9)  # the median date's band over the draws (80%)


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


def quantiles_of_draws(doy: np.ndarray, draws: np.ndarray, qs=PASSAGE_QUANTILES) -> np.ndarray:
    """`cumulative_quantiles` of each row of `draws` (draws x days, in day order): `(n, len(qs))`,
    NaN for an empty draw."""
    total = draws.sum(axis=1, keepdims=True)
    cum = np.cumsum(draws, axis=1) / np.where(total > 0, total, 1)
    rows = np.arange(len(draws))
    out = np.empty((len(draws), len(qs)))
    for j, q in enumerate(qs):
        k = np.minimum(np.argmax(cum >= q - 1e-12, axis=1), len(doy) - 1)
        prev = np.where(k > 0, cum[rows, k - 1], 0.0)
        step = cum[rows, k] - prev
        frac = np.divide(q - prev, step, out=np.ones(len(rows)), where=step > 0)
        out[:, j] = np.where(k > 0, doy[k - 1] + frac * (doy[k] - doy[k - 1]), doy[0])
    out[total[:, 0] <= 0] = np.nan
    return out


def passage_from_draws(frame: pd.DataFrame, draws: np.ndarray) -> pd.DataFrame:
    """Per year of `frame` (`trend.model_frame` over the model window), the 10/50/90% passage dates
    of its gap-filled whole days (`draws`: `trend.fill_draws`, draws x rows of `frame`): the median
    over draws of each, and the 80% band of the median date (`q50_lo`, `q50_hi`, `PASSAGE_BAND`);
    `counted`, the share of the window's days counted.

    A year with no bird counted has no dates:
    they would be the model's alone.
    """
    rows = []
    for y, g in frame.groupby("year"):
        g = g.sort_values("doy")
        q = quantiles_of_draws(g["doy"].to_numpy(), draws[:, g.index.to_numpy()])
        if (g["y"] + g["extra"]).sum() <= 0:
            q[:] = np.nan
        mid = np.nanmedian(q, axis=0) if np.isfinite(q).any() else [np.nan] * len(q[0])
        band = np.nanquantile(q[:, 1], PASSAGE_BAND) if np.isfinite(q).any() else [np.nan] * 2
        counted = ((g["c"] > 0) | (g["extra"] > 0)).mean()
        rows.append([y, *mid, *band, counted])
    cols = ["year", *(f"q{round(q * 100)}" for q in PASSAGE_QUANTILES), "q50_lo", "q50_hi"]
    return pd.DataFrame(rows, columns=[*cols, "counted"]).round(2)


def chance_from_draws(frame: pd.DataFrame, draws: np.ndarray, last_year: int) -> pd.DataFrame:
    """Per day of the last `CHANCE_YEARS` seasons of `frame` (as `passage_from_draws`), the share
    of draws in which its whole day holds at least each of `CHANCE_THRESHOLDS` birds
    (`at_least_<k>`)."""
    on = (frame["year"] > last_year - CHANCE_YEARS).to_numpy()
    out = frame.loc[on, ["year", "doy"]].reset_index(drop=True)
    for k in CHANCE_THRESHOLDS:
        out[f"at_least_{k}"] = (draws[:, on] >= k).mean(axis=0)
    return out


def chances(frame: pd.DataFrame, last_year: int, filled: pd.DataFrame | None = None) -> dict:
    """Per period of `PENTAD` days, the share of days of the last `CHANCE_YEARS` seasons with at
    least each of `CHANCE_THRESHOLDS` birds: from `filled` (`chance_from_draws`) on its days, each
    day weighing its probability; elsewhere from the well-counted days (`c >= COVERAGE_MIN`), the
    birds counted.

    `days` is the days each period averages.
    """
    f = frame[frame["year"] > last_year - CHANCE_YEARS]
    cols = [f"at_least_{k}" for k in CHANCE_THRESHOLDS]
    ok = f["c"] >= COVERAGE_MIN
    day = pd.DataFrame(
        {c: (f["y"] >= k).astype(float).where(ok) for c, k in zip(cols, CHANCE_THRESHOLDS)}
    )
    day[["year", "doy"]] = f[["year", "doy"]]
    if filled is not None:
        day = day.set_index(["year", "doy"])
        p = filled.set_index(["year", "doy"])[cols]
        day.loc[day.index.isin(p.index), cols] = p.reindex(day.index[day.index.isin(p.index)])
        day = day.reset_index()
    day = day.dropna(subset=cols)
    g = day.assign(period=(day["doy"] // PENTAD) * PENTAD).groupby("period")
    out = {"years": [last_year - CHANCE_YEARS + 1, last_year], "doy": g.size().index + PENTAD / 2}
    out["days"] = g.size().to_numpy()
    for c in cols:
        share = g[c].mean()
        if share.max() >= CHANCE_MIN_SHARE:
            out[c] = share.round(3).to_numpy()
    return out


def season_block(
    frame: pd.DataFrame, last_year: int, model_window=None, fill: dict | None = None
) -> dict:
    """The `season` block of a species file: `years` x `doy` shares of the year's birds, the birds
    counted (`count`, null on a day not counted) and the coverage (`c`), each year's passage dates,
    and the chances.

    `frame` spans the view window; the year's total and its passage dates are of the model window
    (`model_window`, season days; all of `frame` if None), so a day shown beyond it has its share
    of that total. `fill`: the trend's gap-filled series (`trend.taxon_trend`'s `season_fill`:
    `days` with the mean whole day of every model-window day, `passage`, `chance`); from the counts
    if None, and for a year with no bird counted in the model window.
    """
    inside = (
        frame["doy"].between(*model_window) if model_window is not None else frame["doy"] > -np.inf
    )
    rate, filled = daily_rates(frame[inside])
    view_rate, _ = daily_rates(frame)
    if fill is None:
        passage = year_passage(view_rate.loc[:, rate.columns], filled)
        chance = chances(frame, last_year)
    else:
        model = fill["days"].pivot(index="year", columns="doy", values="total")
        seen = frame[inside].groupby("year")[["y", "extra"]].sum().sum(axis=1) > 0
        model = model[model.index.isin(seen.index[seen])]  # a year with no bird: from the counts
        filled = model.reindex(index=view_rate.index, columns=rate.columns).fillna(filled)
        view_rate.loc[:, rate.columns] = filled.where(model.reindex_like(filled).notna(), rate)
        passage = fill["passage"]
        chance = chances(frame, last_year, fill["chance"])
    share = view_rate.div(filled.sum(axis=1).where(lambda t: t > 0), axis=0)
    birds = frame.assign(n=frame["y"] + frame["extra"]).where(
        (frame["c"] > 0) | (frame["extra"] > 0)
    )
    birds[["year", "doy"]] = frame[["year", "doy"]]
    count = birds.pivot(index="year", columns="doy", values="n").reindex_like(view_rate)
    c = frame.pivot(index="year", columns="doy", values="c").reindex_like(view_rate)
    return {
        "method": METHOD,
        "source": "counts" if fill is None else "gam",
        "years": view_rate.index.to_numpy(),
        "doy": view_rate.columns.to_numpy(),
        "share": share.to_numpy(),
        "count": count.to_numpy(),
        "c": c.round(2).to_numpy(),
        "passage": passage,
        "chances": chance,
    }
