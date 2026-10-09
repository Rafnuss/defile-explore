"""Effort adjustment for the Explore page: a time-of-day profile per taxon, and coverage.

Part of `defile_explore`. Variant A of the baseline, provisional until compared with a Gaussian-
process variant (DECISIONS.md -> Explore): a taxon's profile p(h | doy) (its own, else its group's,
else uniform over daylight), fitted with the time-of-day GAM
(`defile_explore.timeofday.fit_ratio_surface`, a copy of the forecast's), gives each counted day a
coverage `c`, the share of that day's expected passage in the hours counted. The adjusted day is
count / c, and the annual index is Σ birds / Σ c over counted days in the window: birds per full-
day equivalent. Neither is given below `COVERAGE_MIN`, nor before the taxon's `start_year`.
"""

import numpy as np
import pandas as pd

from defile_explore.export import HOURS, hours_per_solar_hour, in_window, solar_shift, to_local
from defile_explore.release import civil_twilight
from defile_explore.timeofday import RATIO_WEIGHT_POWER, fit_ratio_surface

# A species' profile p(h | doy) is the share of its full-day passage in each solar hour,
# fitted on days timed to the hour (`profile_samples`) and used to say how much of that day's
# passage the hours counted covered (`coverage`). See DECISIONS.md -> Explore.
PROFILE_MIN_TIMED = 0.9  # a sample day has at least this share of its birds timed to the hour
PROFILE_MIN_HOUR_EFFORT = 0.5  # a sample hour was counted for at least half its length
PROFILE_MIN_HOURS = 2  # ...and the day has at least this many such hours
PROFILE_MIN_BIRDS = 500  # timed birds for an own profile; fewer -> the group's, then uniform
# Groups for the fallback profile, by checklist order; other orders are "other".
PROFILE_GROUPS = {
    "Accipitriformes": "raptors",
    "Falconiformes": "raptors",
    "Columbiformes": "pigeons",
    "Passeriformes": "passerines",
}
PROFILE_DOY = (182, 336)  # 1 Jul - 2 Dec (non-leap): the counted season, with a margin
# Below this coverage, a day or a year gets no adjusted value: its raw count and `c` only. Most
# pre-1993 days fall below it (see DECISIONS.md -> Explore).
COVERAGE_MIN = 0.5
PROFILE_EXPORT_STEP = 7  # days between the exported profile's rows


def daylight(doy: np.ndarray) -> np.ndarray:
    """Fraction of each solar hour between civil dawn and dusk, shape (len(doy), 24).

    Days of year are taken in a non-leap year (2025): a day's difference is negligible here.
    """
    dates = pd.Series(pd.Timestamp("2024-12-31") + pd.to_timedelta(doy, unit="D"))
    dawn, dusk = civil_twilight(dates)
    shift = solar_shift(dates)
    return np.array(
        [
            hours_per_solar_hour([(a, b)], s)
            for a, b, s in zip(to_local(dawn), to_local(dusk), shift)
        ]
    )


def profile_samples(
    days: pd.DataFrame, hourly: pd.DataFrame, effort: pd.DataFrame
) -> pd.DataFrame:
    """Hour-of-day ratio samples (`doy`, `hour`, `ratio`, `weight`) for one or more taxa.

    For each day with migrating birds, mostly timed to the hour (`PROFILE_MIN_TIMED`) and counted
    in at least `PROFILE_MIN_HOURS` solar hours: every hour counted for at least
    `PROFILE_MIN_HOUR_EFFORT` gives `ratio` = its rate / the day's rate over those hours,
    zero when no bird was timed to it. `weight` = the day's birds ** `RATIO_WEIGHT_POWER`, as in
    the model's phenology.
    """
    d = days[(days["count"] > 0) & (days["timed"] >= PROFILE_MIN_TIMED)][["taxon_id", "date"]]
    eff = effort.set_index("date")["hourly"]
    d = d[d["date"].isin(eff.index)]
    if d.empty:
        return pd.DataFrame(columns=["doy", "hour", "ratio", "weight"])
    cover = np.stack(d["date"].map(eff).to_numpy())  # (n, 24)
    birds = np.zeros_like(cover)
    pos = pd.Series(np.arange(len(d)), index=pd.MultiIndex.from_frame(d))
    h = hourly.set_index(["taxon_id", "date"])
    h = h[h.index.isin(pos.index)]
    np.add.at(birds, (pos[h.index].to_numpy(), h["hour"].to_numpy()), h["count"].to_numpy())
    used = cover >= PROFILE_MIN_HOUR_EFFORT
    birds_used = np.where(used, birds, 0).sum(axis=1)
    ok = (used.sum(axis=1) >= PROFILE_MIN_HOURS) & (birds_used > 0)
    day_rate = birds_used / np.where(used, cover, 0).sum(axis=1).clip(min=1e-9)
    rows, hours = np.nonzero(used & ok[:, None])
    return pd.DataFrame(
        {
            "doy": d["date"].dt.dayofyear.to_numpy()[rows],
            "hour": hours,
            "ratio": birds[rows, hours] / cover[rows, hours] / day_rate[rows],
            "weight": birds_used[rows] ** RATIO_WEIGHT_POWER,
        }
    )


def fit_profile(samples: pd.DataFrame, doy_grid: np.ndarray, light: np.ndarray, label=""):
    """P(h | doy) on `doy_grid` x 24 solar hours: the fitted ratio surface on the hours the samples
    span, times the `daylight` fraction, normalised to sum to 1 on each day."""
    if samples.empty:
        raise ValueError(f"{label}: no profile samples")
    hours = np.arange(samples["hour"].min(), samples["hour"].max() + 1)
    ratio = fit_ratio_surface(
        samples["doy"].to_numpy(),
        samples["hour"].to_numpy(),
        samples["ratio"].to_numpy(),
        samples["weight"].to_numpy(),
        doy_grid,
        hours,
        label=label,
    )
    p = np.zeros((len(doy_grid), HOURS))
    p[:, hours] = np.clip(ratio, 0, None)
    p *= light
    return p / p.sum(axis=1, keepdims=True)


def uniform_profile(light: np.ndarray) -> np.ndarray:
    """Every daylight minute equal: coverage is then the share of daylight counted."""
    return light / light.sum(axis=1, keepdims=True)


def build_profiles(
    taxa: pd.DataFrame,
    days: pd.DataFrame,
    hourly: pd.DataFrame,
    effort: pd.DataFrame,
    only=None,
) -> tuple[dict, pd.Series]:
    """`(profiles, source)`: one p(h | doy) per taxon id, group and "uniform" (on `PROFILE_DOY`),
    and each taxon's profile source (`own`, a group name, or `uniform`).

    With `only` (taxon ids), just the profiles those taxa use, and `source` for them alone: a
    group's profile is fitted on all its members, but only if one of them lacks its own (the
    passerine and "other" groups are most of the time here).
    """
    doy_grid = np.arange(PROFILE_DOY[0], PROFILE_DOY[1] + 1)
    light = daylight(doy_grid)
    profiles = {"uniform": uniform_profile(light)}
    timed = hourly.groupby("taxon_id")["count"].sum()
    group = taxa.set_index("taxon_id")["order"].map(PROFILE_GROUPS).fillna("other")
    wanted = group.index if only is None else pd.Index(only)
    source = pd.Series("uniform", index=wanted)
    for taxon_id in timed.index[(timed >= PROFILE_MIN_BIRDS) & timed.index.isin(wanted)]:
        samples = profile_samples(
            days[days["taxon_id"] == taxon_id], hourly[hourly["taxon_id"] == taxon_id], effort
        )
        if len(samples):
            profiles[taxon_id] = fit_profile(samples, doy_grid, light, label=taxon_id)
            source[taxon_id] = "own"
    needed = group if only is None else group[source.index[source != "own"]]
    for g in sorted(needed.unique()):
        members = group.index[group == g]
        samples = profile_samples(
            days[days["taxon_id"].isin(members)], hourly[hourly["taxon_id"].isin(members)], effort
        )
        if len(samples):
            profiles[g] = fit_profile(samples, doy_grid, light, label=g)
            source[source.index.isin(members) & (source != "own")] = g
    return profiles, source


def coverage(effort: pd.DataFrame, profile: np.ndarray) -> pd.Series:
    """C per counted day: the share of the day's expected passage (under `profile`, on
    `PROFILE_DOY`) in the hours counted.

    Days outside `PROFILE_DOY` use its nearest day.
    """
    doy = effort["date"].dt.dayofyear.clip(*PROFILE_DOY) - PROFILE_DOY[0]
    cover = np.stack(effort["hourly"].to_numpy())
    return pd.Series((profile[doy.to_numpy()] * cover).sum(axis=1), index=effort.index)


def annual_index(
    days: pd.DataFrame, effort: pd.DataFrame, c: pd.Series, start_year: int = 0
) -> pd.DataFrame:
    """Per year, over counted days in the window: birds, Σc and the index Σbirds / Σc (birds per
    full-day equivalent), none below `COVERAGE_MIN` or before `start_year`.

    Counted days without a record of the taxon count as zero birds.
    """
    e = effort.assign(c=c)
    e = e[(e["state"] == "counted") & in_window(e["date"])]
    birds = e["date"].map(days.set_index("date")["count"]).fillna(0)
    g = e.assign(birds=birds, year=e["date"].dt.year).groupby("year")
    a = g.agg(birds=("birds", "sum"), c=("c", "sum"), hours=("hours", "sum"), days=("c", "size"))
    a["c_mean"] = a["c"] / a["days"]
    ok = (a["c_mean"] >= COVERAGE_MIN) & (a.index >= start_year)
    a["index"] = (a["birds"] / a["c"]).where(ok)
    return a.reset_index()


def adjust_days(
    days: pd.DataFrame, effort: pd.DataFrame, c: pd.Series, start_year: int = 0
) -> pd.DataFrame:
    """`days` plus `c` (that day's coverage; NaN if not a counted day) and `adjusted` = count / c
    when c >= `COVERAGE_MIN`, from `start_year`."""
    counted = effort[effort["state"] == "counted"]
    d = days.assign(c=days["date"].map(c[counted.index].set_axis(counted["date"])))
    ok = (d["c"] >= COVERAGE_MIN) & (d["date"].dt.year >= start_year)
    d["adjusted"] = (d["count"] / d["c"]).where(ok)
    return d


def profile_table(profile: np.ndarray) -> dict:
    """The profile every `PROFILE_EXPORT_STEP` days, as `{"doy": [...], "p": [[24 values]...]}`."""
    rows = np.arange(0, profile.shape[0], PROFILE_EXPORT_STEP)
    return {"doy": (rows + PROFILE_DOY[0]).tolist(), "p": np.round(profile[rows], 4).tolist()}
