"""Each taxon's windows: the days of the season its model is fitted on, and the days its season
panels show.

The default window (`export.WINDOW`, 18 Jul - 18 Nov) is where the count is near-complete: most
days are counted in 80-100% of years, against about a third from 20 Nov to 1 Dec. Late taxa still
pass beyond it (Red Kite: still around 120 birds a day on 30 Nov in the years counted then), so a
taxon's windows are extended, by rule:

- `envelope`: the season days counted (`c >= COVERAGE_MIN`) in at least a share of the taxon's
  years (from its start year), smoothed over `ENVELOPE_SMOOTH` days so a one-off gap does not cut
  it: `MODEL_MIN_YEARS` for the model, which must fill the rest, `VIEW_MIN_YEARS` for the panels,
  which only show what was counted.
- `passage_tails`: from the taxon's mean birds per full counted day by date over its years
  (`climatology`), the dates by which `tails` of the passage within an envelope have passed.
- `windows`: the model window is the default one extended to `MODEL_TAILS` of the passage plus
  `MARGIN_DAYS`, within the model envelope; the view window is the model window extended the same
  way to `VIEW_TAILS`, within the view envelope. Only extended, never cut: a taxon whose passage
  sits inside the default window keeps the fit benchmarked on it.
- `edge_shares`: the share of the passage within the model envelope that falls in its first and
  last `EDGE_DAYS`. At `EDGE_SHARE` or more, much of the passage likely falls where counting thins
  out (`beyond_counting`): its totals and dates are those of the counted season, and the page says
  so. A test on the tails alone (`MODEL_TAILS` within `MARGIN_DAYS` of the edge) flagged 51 taxa,
  most for a trickle of late birds.

Windows are pairs of season days (`trend.season_day`: the day of year of a non-leap year, so one
calendar date is one day every year), inclusive. `settings` can override both.
"""

import numpy as np
import pandas as pd

from defile_explore.export import WINDOW
from defile_explore.profile import COVERAGE_MIN

METHOD = "window@2"
MODEL_MIN_YEARS = 1 / 3  # share of years a date must be counted in to be modelled
VIEW_MIN_YEARS = 0.1  # ... to be shown
ENVELOPE_SMOOTH = 7  # days of the rolling mean of that share (one-off gaps do not cut it)
CLIMATOLOGY_SMOOTH = 7  # days of the rolling mean of the mean daily rate
MODEL_TAILS = (0.005, 0.995)  # passage shares the model window covers
VIEW_TAILS = (0.001, 0.999)  # ... and the view window
MARGIN_DAYS = 5  # added beyond those dates
EDGE_DAYS = 7  # days at each end of the model envelope...
EDGE_SHARE = 0.1  # ... holding this share of the passage or more: `beyond_counting`
NON_LEAP = 2001  # a non-leap year, to turn (month, day) into a season day and back


def season_day_of(month: int, day: int) -> int:
    return pd.Timestamp(NON_LEAP, month, day).dayofyear


def default_window() -> tuple[int, int]:
    return season_day_of(*WINDOW[0]), season_day_of(*WINDOW[1])


def as_dates(window: tuple[int, int]) -> list[str]:
    """A window as `MM-DD` strings."""
    base = pd.Timestamp(NON_LEAP, 1, 1)
    return [(base + pd.Timedelta(days=d - 1)).strftime("%m-%d") for d in window]


def from_dates(dates: list[str]) -> tuple[int, int]:
    """A window from `MM-DD` strings (as in `overrides.yaml`)."""
    return tuple(season_day_of(*map(int, d.split("-"))) for d in dates)


def envelope(frame: pd.DataFrame, min_years: float) -> tuple[int, int]:
    """The run of season days around the default window counted (`c >= COVERAGE_MIN`) in at least
    `min_years` of the years of `frame` (a `trend.model_frame` over the whole counted season); the
    default window's own days always belong to it."""
    share = (frame["c"] >= COVERAGE_MIN).groupby(frame["doy"]).mean()
    share = share.rolling(ENVELOPE_SMOOTH, center=True, min_periods=1).mean()
    lo, hi = default_window()
    ok = (share >= min_years) | ((share.index >= lo) & (share.index <= hi))
    days = ok.index.to_numpy()
    a = days[(days <= lo) & ~ok.to_numpy()]
    b = days[(days >= hi) & ~ok.to_numpy()]
    return int(a.max() + 1 if len(a) else days.min()), int(b.min() - 1 if len(b) else days.max())


def climatology(frame: pd.DataFrame) -> pd.Series:
    """Mean birds per full counted day by season day, over the years it was counted, smoothed over
    `CLIMATOLOGY_SMOOTH` days; NaN where never counted."""
    rate = (frame["y"] / frame["c"]).where(frame["c"] >= COVERAGE_MIN)
    mean = rate.groupby(frame["doy"]).mean()
    return mean.rolling(CLIMATOLOGY_SMOOTH, center=True, min_periods=1).mean()


def passage_tails(clim: pd.Series, env: tuple[int, int], tails) -> tuple[float, float] | None:
    """Season days by which `tails` of the passage within `env` have passed; None without birds."""
    c = clim.loc[env[0] : env[1]].fillna(0)
    if c.sum() <= 0:
        return None
    cum = np.cumsum(c.to_numpy()) / c.sum()
    return tuple(float(np.interp(q, cum, c.index.to_numpy())) for q in tails)


def edge_shares(clim: pd.Series, env: tuple[int, int]) -> tuple[float, float] | None:
    """Shares of the passage within `env` in its first and last `EDGE_DAYS`; None without birds."""
    c = clim.loc[env[0] : env[1]].fillna(0)
    if c.sum() <= 0:
        return None
    return float(c.iloc[:EDGE_DAYS].sum() / c.sum()), float(c.iloc[-EDGE_DAYS:].sum() / c.sum())


def extend(window, tails, env) -> tuple[int, int]:
    """`window`, extended towards `tails` plus `MARGIN_DAYS`, but not beyond `env`."""
    if tails is None:
        return window
    lo = min(window[0], max(env[0], int(np.floor(tails[0])) - MARGIN_DAYS))
    hi = max(window[1], min(env[1], int(np.ceil(tails[1])) + MARGIN_DAYS))
    return lo, hi


def windows(frame: pd.DataFrame) -> dict:
    """The rule's windows from `frame` (a `trend.model_frame` over the whole counted season, from
    the taxon's start year): `model` and `view`, the envelopes and passage tails they come from,
    the shares of the passage at the model envelope's edges (`edge_shares`) and `beyond_counting`
    (`"early"`, `"late"`, both, or neither: `EDGE_SHARE` or more at that edge)."""
    clim = climatology(frame)
    env_m, env_v = envelope(frame, MODEL_MIN_YEARS), envelope(frame, VIEW_MIN_YEARS)
    tails_m = passage_tails(clim, env_m, MODEL_TAILS)
    model = extend(default_window(), tails_m, env_m)
    view = extend(model, passage_tails(clim, env_v, VIEW_TAILS), env_v)
    edges = edge_shares(clim, env_m)
    beyond = [k for k, e in zip(("early", "late"), edges or ()) if e >= EDGE_SHARE]
    return {
        "method": METHOD,
        "model": model,
        "view": view,
        "model_envelope": env_m,
        "view_envelope": env_v,
        "tails": tails_m,
        "edge_shares": edges,
        "beyond_counting": beyond,
    }
