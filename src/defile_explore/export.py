"""Explore export: the defile-dataset release tables -> the JSON files defileViz's Explore page
reads.

Unlike `src/data/counts.py`, nothing here is a model choice: every value is a plain aggregation
of the release tables (`count.csv`, `survey.csv`, `taxonomy.csv`), following
the dataset README's daily-total rules (a count's own timing, else its survey's; the local
Europe/Paris day). No count is moved, dropped, imputed or redistributed, so raw daily and annual
totals reconcile with `count.csv` (`tests/test_explore.py`). Entries the release keeps at day
level (untimed, or timed outside their survey) stay untimed here.

Files, written by `scripts/build_explore.py` (one `build_*` function each):

- `manifest.json`: when and from which dataset build it was made, the date range, partial years,
  the default window and the tier thresholds.
- `taxa.json`: one entry per taxon, with names (French from the eBird taxonomy), rank,
  occurrence tier, totals, and `start_year`, the first year its counts are comparable from. A group
  taxon ("harrier sp.") is read as everything below it (`add_rollups`), with its `members`.
- `effort.json`: one entry per local day with any survey: hours counted (union of complete
  surveys, `survey_complete`, weather stops apart), first/last time, hours covered per solar hour
  (weather stops included: no bird passed), and the hours of weather stops and incomplete surveys.
- `species/<taxon_id>.json`: per day, the main-direction count with its qualifiers, the reverse
  and local counts, the share of birds timed to the hour, the day's coverage `c` and the
  effort-adjusted count; per solar hour, the birds timed to it; per year, totals and the
  effort-adjusted index; the time-of-day profile used; and for
  full-tier taxa with a trend, `trend` (`defile_explore.trend.taxon_trend`): gap-filled
  annual totals with intervals, the smooth trend, the median passage date per year, the smooth
  season in the first and last year, and each year's season with and without its weather episodes.
  `defile_explore.pipeline.build_taxon` adds the blocks derived from the same data and fit
  (`settings`, `window`, `links`, `key_numbers`, `season`, `daytime`, `age`, `sex`, `records`,
  `accounts`, `benchmark`, `reliability`, `diagnostics`); each of those modules documents its block.

Tables are columnar (`{"date": [...], "count": [...]}`) to keep the files small.

Hours of the day are local apparent solar time (`solar_shift`): hour 12 starts at the sun's transit
over Defile, whatever the clock says. Days are local calendar dates. A day's passage keeps its place
in solar time through the season, where on the clock it moves by an hour at the end of summer time
and by 20 min with the equation of time (DECISIONS.md -> Explore).
"""

import datetime as dt
import json
import os

import numpy as np
import pandas as pd

from defile_explore.release import (
    DATASET_DIR,
    MAIN_CATEGORY,
    METADATA_FILE,
    PRESENCE_ONLY,
    SITE,
    TIMEZONE,
    parse_counts,
    parse_surveys,
)

# Bump when a field's meaning changes, so defileViz can refuse an export it does not understand.
DEFINITIONS_VERSION = 2  # 2: hours of the day in solar time

# Default comparison window (month, day), inclusive. A window for annual totals, not a filter:
# every day is exported.
WINDOW = ((7, 18), (11, 18))

# Occurrence tiers. "full": a dashboard; "short": records, annual totals and report links;
# "rare": the rare-records table only. On days with migrating birds (`tier_of`): a seasonal curve
# needs days, not birds (one large flock is still a rare record). Provisional, to tune on the page.
TIER_FULL_MIN_DAYS = 50
TIER_FULL_MIN_YEARS = 5
TIER_RARE_MAX_DAYS = 10

# Start year: the first year a taxon's counts are comparable from (`start_year`); earlier years
# are exported raw but get no adjusted value, and defileViz crops its annual series there.
# Systematic daily counting began in 1993. Its targets were always raptors, herons and egrets,
# storks, pigeons and corvids, and other large birds counted individually (cranes, geese, ducks); passerines were hardly
# recorded before 2007, when the number of passerine taxa doubles and their birds rise 15-fold
# (defile-dataset docs/sampling-history.md). Targets start in 1993 even when rare then (Peregrine
# recovered); any other taxon starts in 1993 only if recorded in at least `START_MIN_SHARE` as
# large a share of the 1993-2006 years as of the later ones, else in 2007. A year counts as
# recorded only with at least `START_MIN_YEAR_BIRDS` of the taxon's median year since 2007: a few
# birds noted in a year it was not counted (swallows 2000-2006: under 1% of a year since) are not
# a series. A floor relative to recent years also penalises a real increase, so it stays low.
SYSTEMATIC_FROM = 1993  # daily systematic counting starts (defile-dataset sampling history)
ALL_GROUPS_FROM = 2007
TARGET_ORDERS = {
    "Accipitriformes",
    "Falconiformes",
    "Pelecaniformes",
    "Ciconiiformes",
    "Columbiformes",
    "Gruiformes",
    "Anseriformes",
}
TARGET_FAMILIES = {"Corvidae"}
START_MIN_SHARE = 0.75
START_MIN_YEAR_BIRDS = 0.02

# Roll-up series: a group taxon ("harrier sp.", "Red/Black Kite", "swallow sp.") is read as all the
# birds counted under it: its own count plus those of every taxon below it in the release's
# `parent_taxon_id` tree (defile-dataset `taxonomy/parent_taxa.csv`). This makes a group's series
# independent of how far the birds were identified (swallows were "swallow sp." in 1993-1999 and
# are increasingly identified since 2021; "Columba sp." appears only from 2014). Each is exported
# like a taxon (`add_rollups`), with its `members`, and gets a start year by the same rule. A group
# adds nothing to read when it is alone (fewer than `ROLLUP_MIN_BELOW` taxa below it: "Short-eared
# Owl" is all "Asio sp." can be) or when one taxon is more than `ROLLUP_MAX_SHARE` of its birds
# (the group is that taxon), or when it is too large to say anything (more than `ROLLUP_MAX_BELOW`
# taxa below it: "bird sp.", "passerine sp."): it is `EXCLUDED_TIER`, with no page, and keeps its
# own series. A species is always its own birds plus those of its subspecies, whatever their share.
ROLLUP_MIN_BELOW = 2
ROLLUP_MAX_BELOW = 100
ROLLUP_MAX_SHARE = 0.9
EXCLUDED_TIER = "excluded"

# French names: the eBird taxonomy in French (France), joined on `ebird_code`; no API key needed.
# Downloaded once into the data dir by `scripts/build_explore.py` (`--refresh-names` to update).
EBIRD_TAXONOMY_URL = "https://api.ebird.org/v2/ref/taxonomy/ebird?fmt=csv&locale=fr_FR"
EBIRD_TAXONOMY_FILE = os.path.join("count", "ebird_taxonomy_fr_FR.csv")
# Taxa whose `ebird_code` is no longer in the eBird taxonomy.
FRENCH_NAME_OVERRIDES = {"avibase-81B32602": "Corneille mantelée"}  # hoocro4, Corvus c. cornix

# A count is "timed" when it has its own time, or inherits a survey at most this long.
TIMED_MAX_DURATION = pd.Timedelta(hours=1)


HOURS = 24
# A solar hour is cut into this many slots: survey edges, dawn and dusk fall within a slot, and the
# time-of-day profile is read at this resolution, so coverage integrates it over the minutes counted.
STEPS_PER_HOUR = 12  # 5-minute slots
SLOTS = HOURS * STEPS_PER_HOUR
DECIMALS = 3  # hours and fractions


# --- local and solar time --------------------------------------------------


def to_local(t: pd.Series) -> pd.Series:
    """UTC times -> naive local wall-clock times."""
    return t.dt.tz_convert(TIMEZONE).dt.tz_localize(None)


def solar_shift(dates: pd.Series) -> pd.Series:
    """Hours from local apparent solar time to local clock time on each local date: clock time less
    this is solar time, 12:00 at the sun's transit over Defile (`SITE`).

    The zone's offset from UTC at local noon (summer time included), less the longitude's offset
    and the equation of time (NOAA's Fourier series, within a minute). Constant through a day:
    summer time changes at night.
    """
    days = pd.DatetimeIndex(pd.to_datetime(dates.unique()))
    noon = days + pd.Timedelta(hours=12)
    utc_noon = noon.tz_localize(TIMEZONE).tz_convert("UTC").tz_localize(None)
    zone = (noon - utc_noon) / pd.Timedelta(hours=1)
    g = 2 * np.pi / 365 * (days.dayofyear - 1)
    eot = 229.18 * (  # minutes
        0.000075
        + 0.001868 * np.cos(g)
        - 0.032077 * np.sin(g)
        - 0.014615 * np.cos(2 * g)
        - 0.040849 * np.sin(2 * g)
    )
    shift = np.asarray(zone) - SITE[1] / 15 - np.asarray(eot) / 60
    return pd.Series(dates).map(dict(zip(days, shift)))


def split_at_midnight(start: pd.Series, end: pd.Series) -> pd.DataFrame:
    """Local intervals -> pieces within one local day (`index` points back to the input row)."""
    pieces = []
    for i, s, e in zip(start.index, start, end):
        while s < e:
            cut = min(e, s.normalize() + pd.Timedelta(days=1))
            pieces.append((i, s, cut))
            s = cut
    return pd.DataFrame(pieces, columns=["index", "start", "end"]).set_index("index")


def merge_intervals(start: np.ndarray, end: np.ndarray) -> list[tuple]:
    """Union of intervals, as sorted non-overlapping (start, end) pairs."""
    out = []
    for s, e in sorted(zip(start, end)):
        if out and s <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], e))
        else:
            out.append((s, e))
    return out


def hours_per_slot(intervals: list[tuple], shift: float) -> np.ndarray:
    """Hours covered in each of the `SLOTS` solar slots of a day by (start, end) intervals in local
    clock time within that day, `shift` the day's `solar_shift`; what falls outside the solar day
    (at night) is left out."""
    cover = np.zeros(SLOTS)
    for s, e in intervals:
        day = s.normalize()
        a = ((s - day) / pd.Timedelta(hours=1) - shift) * STEPS_PER_HOUR
        b = ((e - day) / pd.Timedelta(hours=1) - shift) * STEPS_PER_HOUR
        for k in range(max(int(np.floor(a)), 0), min(int(np.ceil(b)), SLOTS)):
            cover[k] += max(0.0, min(b, k + 1) - max(a, k))
    return cover / STEPS_PER_HOUR


def by_hour(x: np.ndarray) -> np.ndarray:
    """Slots (last axis, `SLOTS`) summed into the 24 solar hours."""
    return x.reshape(*x.shape[:-1], HOURS, STEPS_PER_HOUR).sum(axis=-1)


# --- effort -----------------------------------------------------------------


def _flag(values: pd.Series, default: bool) -> pd.Series:
    """A release boolean column (`true`/`false`, blank: `default`) as bool."""
    return values.map(
        lambda v: default if pd.isna(v) or v == "" else str(v).strip().lower() == "true"
    )


def _hours(day: pd.DataFrame) -> float:
    union = merge_intervals(day["start"].tolist(), day["end"].tolist())
    return sum((e - b) / pd.Timedelta(hours=1) for b, e in union)


def build_effort(surveys: pd.DataFrame) -> pd.DataFrame:
    """One row per local day with any survey.

    A survey with `survey_complete` holds every bird that passed, weather stops included (counting
    impossible, no bird: `weather_stop`). `slots` is the share of each solar slot covered by the
    union of complete surveys (`STEPS_PER_HOUR`), weather stops included, so their zeros count;
    `hourly` the same per solar hour; `first`/`last` its local bounds (clock time); `periods` the
    number of complete surveys. `hours` is the union of complete surveys without the weather stops
    (hours actually counted), `weather_hours` and `incomplete_hours` the others'. `state`:
    `counted` (some complete survey), else `uncertain` (only incomplete ones).
    """
    s = surveys.assign(lstart=to_local(surveys["start"]), lend=to_local(surveys["end"]))
    s["complete"] = _flag(s.get("survey_complete", pd.Series(index=s.index, dtype=object)), True)
    s["weather"] = _flag(s.get("weather_stop", pd.Series(index=s.index, dtype=object)), False)
    pieces = split_at_midnight(s["lstart"], s["lend"])
    pieces = pieces.join(s[["complete", "weather", "recording_era"]])
    pieces["date"] = pieces["start"].dt.normalize()
    shift = solar_shift(pd.Series(pieces["date"].unique()))
    shift = dict(zip(pieces["date"].unique(), shift))

    rows = []
    for date, day in pieces.groupby("date"):
        row = {"date": date}
        complete = day[day["complete"]]
        union = merge_intervals(complete["start"].tolist(), complete["end"].tolist())
        slots = hours_per_slot(union, shift[date])
        row["hours"] = _hours(complete[~complete["weather"]])
        row["weather_hours"] = _hours(complete[complete["weather"]])
        row["incomplete_hours"] = _hours(day[~day["complete"]])
        row["first"] = union[0][0] if union else None
        row["last"] = union[-1][1] if union else None
        row["periods"] = len(complete)
        row["slots"] = slots * STEPS_PER_HOUR
        row["hourly"] = by_hour(slots)
        eras = day.loc[day["recording_era"] != "curated", "recording_era"]
        row["era"] = eras.mode().iloc[0] if len(eras) else "curated"
        rows.append(row)
    e = pd.DataFrame(rows)
    e["state"] = np.where(e["periods"] > 0, "counted", "uncertain")
    return e


def annual_effort(effort: pd.DataFrame) -> pd.DataFrame:
    """Per year: days by `state`, counted hours, and the same in the default window.

    `weather_days`: days wholly stopped by weather (counted, no bird), within `days`.
    """
    e = effort.assign(year=effort["date"].dt.year, win=in_window(effort["date"]))
    counted = e["state"] == "counted"
    return (
        e.assign(
            days=counted,
            uncertain_days=e["state"] == "uncertain",
            weather_days=(e["weather_hours"] > 0) & (e["hours"] == 0),
            window_days=counted & e["win"],
            window_hours=e["hours"].where(e["win"], 0),
        )
        .groupby("year")[
            ["days", "hours", "window_days", "window_hours", "uncertain_days", "weather_days"]
        ]
        .sum()
        .reset_index()
    )


# --- counts -----------------------------------------------------------------


def timed_span(counts: pd.DataFrame, surveys: pd.DataFrame) -> pd.DataFrame:
    """`solar_start`, `solar_end`: the solar hours of its day each count is timed to, NaN if
    untimed.

    Timed: the count's own time (a point, start = end), or, without one, its survey's interval when
    at most `TIMED_MAX_DURATION` long (historical hour-by-hour records). Day-level entries are
    untimed.
    """
    sv = surveys.set_index("survey_id")
    start = to_local(counts["survey_id"].map(sv["start"]))
    end = to_local(counts["survey_id"].map(sv["end"]))
    own = counts["datetime"].notna()
    inherits = counts["raw_datetime"].isna() & ((end - start) <= TIMED_MAX_DURATION) & ~own
    shift = solar_shift(counts["date"]).set_axis(counts.index)

    def solar(t):
        return (t - counts["date"]) / pd.Timedelta(hours=1) - shift

    span = pd.DataFrame(np.nan, index=counts.index, columns=["solar_start", "solar_end"])
    point = solar(to_local(counts["datetime"].where(own)))
    span.loc[own, "solar_start"] = span.loc[own, "solar_end"] = point[own]
    span.loc[inherits, "solar_start"] = solar(start)[inherits]
    span.loc[inherits, "solar_end"] = solar(end)[inherits]
    return span


def daily_counts(counts: pd.DataFrame) -> pd.DataFrame:
    """One row per taxon and local day with any count row.

    `count`: main-direction birds (NaN if only presence-only records); `reverse`, `local`: NaN when
    no row (unknown, not zero); `qualifiers`: the day's distinct `count_estimation` codes on main-
    direction rows, sorted (e.g. `"~>"`); `timed`: the share of main-direction birds timed to the
    hour (NaN without birds).
    """
    normal = counts[counts["count_category"] == MAIN_CATEGORY]
    key = ["taxon_id", "date"]
    d = normal.groupby(key).agg(
        count=("count", lambda x: x.sum(min_count=1)),
        timed_birds=("timed_count", "sum"),
        qualifiers=("count_estimation", lambda x: "".join(sorted(set(x.dropna())))),
    )
    for cat in ("reverse", "local"):
        rows = counts[counts["count_category"] == cat]
        d = d.join(rows.groupby(key)["count"].sum().rename(cat), how="outer")
    d["qualifiers"] = d["qualifiers"].fillna("")
    d["timed"] = (d["timed_birds"] / d["count"]).where(d["count"] > 0)
    return d.drop(columns="timed_birds").reset_index()


def hourly_counts(counts: pd.DataFrame) -> pd.DataFrame:
    """Main-direction birds per taxon, local day and solar hour, from timed counts only.

    A count timed to an interval is spread over the solar hours it overlaps, in proportion (an
    hour-by-hour sheet's clock hour falls across two solar hours); one timed to a point falls in
    that point's hour. What falls outside the solar day is left out.
    """
    t = counts[(counts["count_category"] == MAIN_CATEGORY) & counts["solar_start"].notna()]
    a, b = t["solar_start"].to_numpy()[:, None], t["solar_end"].to_numpy()[:, None]
    edge = np.arange(HOURS)[None, :]
    overlap = np.clip(np.minimum(b, edge + 1) - np.maximum(a, edge), 0, None)
    point = (b <= a).ravel()
    share = np.where(
        point[:, None], np.floor(a) == edge, overlap / np.where(point, 1, (b - a).ravel())[:, None]
    )
    i, hour = np.nonzero(share)
    h = pd.DataFrame(
        {
            "taxon_id": t["taxon_id"].to_numpy()[i],
            "date": t["date"].to_numpy()[i],
            "hour": hour,
            "count": t["count"].to_numpy()[i] * share[i, hour],
        }
    )
    return h.groupby(["taxon_id", "date", "hour"])["count"].sum().reset_index()


def in_window(dates: pd.Series) -> pd.Series:
    """Whether local dates fall in the default comparison window `WINDOW`."""
    md = dates.dt.month * 100 + dates.dt.day
    (m0, d0), (m1, d1) = WINDOW
    return (md >= m0 * 100 + d0) & (md <= m1 * 100 + d1)


def annual(days: pd.DataFrame) -> pd.DataFrame:
    """Per year: total birds, birds in the window, days with a record, and the top day."""
    d = days.assign(year=days["date"].dt.year, win=in_window(days["date"]))
    d["win_count"] = d["count"].where(d["win"])
    g = d.groupby("year")
    a = g.agg(
        total=("count", "sum"),
        window=("win_count", "sum"),
        days=("date", "size"),
    )
    top = d.loc[d["count"].fillna(-1).groupby(d["year"]).idxmax()]
    a["max"] = top.set_index("year")["count"]
    a["max_date"] = top.set_index("year")["date"]
    return a.reset_index()


# --- taxa -------------------------------------------------------------------


def tier_of(years: int, days: int) -> str:
    """Occurrence tier from the years and days with migrating birds (thresholds above)."""
    if days <= TIER_RARE_MAX_DAYS:
        return "rare"
    if days >= TIER_FULL_MIN_DAYS and years >= TIER_FULL_MIN_YEARS:
        return "full"
    return "short"


def start_year(order: str, family: str, birds: dict[int, float], last_year: int) -> int:
    """First comparable year of a taxon (see `SYSTEMATIC_FROM`), from its order, family and birds
    per year (years with migrating birds only); `last_year` is the last year of the data."""
    if order in TARGET_ORDERS or family in TARGET_FAMILIES:
        return SYSTEMATIC_FROM
    after = [n for y, n in birds.items() if ALL_GROUPS_FROM <= y <= last_year]
    floor = START_MIN_YEAR_BIRDS * np.median(after) if after else 0
    years = [y for y, n in birds.items() if n >= floor]
    before = sum(SYSTEMATIC_FROM <= y < ALL_GROUPS_FROM for y in years)
    after = sum(ALL_GROUPS_FROM <= y <= last_year for y in years)
    share_before = before / (ALL_GROUPS_FROM - SYSTEMATIC_FROM)
    share_after = after / (last_year - ALL_GROUPS_FROM + 1)
    systematic = before > 0 and share_before >= START_MIN_SHARE * share_after
    return SYSTEMATIC_FROM if systematic else ALL_GROUPS_FROM


def french_names(taxonomy: pd.DataFrame, ebird: pd.DataFrame) -> pd.Series:
    """French name of each taxon (eBird `COMMON_NAME` by `ebird_code`), first letter capitalised
    ("labbe sp." -> "Labbe sp."), with `FRENCH_NAME_OVERRIDES`."""
    names = taxonomy["ebird_code"].map(ebird.set_index("SPECIES_CODE")["COMMON_NAME"])
    overrides = FRENCH_NAME_OVERRIDES
    names = names.fillna(taxonomy["taxon_id"].map(overrides))
    return names.str[:1].str.upper() + names.str[1:]


def build_taxa(
    taxonomy: pd.DataFrame, days: pd.DataFrame, ebird: pd.DataFrame | None = None
) -> pd.DataFrame:
    """Taxonomy plus French name and occurrence: years and days with migrating birds (or a
    presence-only record), first/last year, total birds, tier and `start_year`."""
    d = days[(days["count"] > 0) | days["qualifiers"].str.contains(PRESENCE_ONLY)]
    occ = d.groupby("taxon_id").agg(
        years=("date", lambda x: x.dt.year.nunique()),
        first_year=("date", lambda x: x.dt.year.min()),
        last_year=("date", lambda x: x.dt.year.max()),
        days=("date", "size"),
        birds=("count", "sum"),
    )
    cols = ["taxon_id", "english_name", "scientific_name", "taxon_rank", "order", "family"]
    extra = [c for c in ("members", "rollup_excluded") if c in taxonomy]
    t = taxonomy[cols + extra].join(occ, on="taxon_id")
    if ebird is not None:
        t.insert(1, "french_name", french_names(taxonomy, ebird))
    t["tier"] = [tier_of(y, n) if pd.notna(n) else "rare" for y, n in zip(t["years"], t["days"])]
    if "rollup_excluded" in t:
        t.loc[t["rollup_excluded"].fillna(False).astype(bool), "tier"] = EXCLUDED_TIER
        t = t.drop(columns="rollup_excluded")
    per_year = d.groupby(["taxon_id", d["date"].dt.year])["count"].sum()
    recorded = {i: g.droplevel(0).to_dict() for i, g in per_year.groupby(level=0)}
    last = int(days["date"].max().year)
    t["start_year"] = [
        start_year(o, f, recorded.get(i, {}), last)
        for i, o, f in zip(t["taxon_id"], t["order"], t["family"])
    ]
    return t.sort_values(["order", "family", "scientific_name"], na_position="last")


def descendants(parents: pd.Series) -> dict[str, list[str]]:
    """`{taxon_id: every taxon below it}` from `parents` (`parent_taxon_id` by `taxon_id`), for the
    taxa that have any."""
    children: dict[str, list[str]] = {}
    for child, parent in parents.dropna().items():
        children.setdefault(parent, []).append(child)
    out = {}
    for taxon_id in children:
        found, todo = [], list(children[taxon_id])
        while todo:
            node = todo.pop()
            found.append(node)
            todo.extend(children.get(node, []))
        out[taxon_id] = found
    return out


def add_rollups(
    taxonomy: pd.DataFrame, days: pd.DataFrame, hourly: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """`taxonomy`, `days` and `hourly` with each group taxon's series replaced by the sum over the
    group and everything below it (`parent_taxon_id`): summed counts (`count`, `reverse`, `local`
    NaN only when no member has one), merged qualifiers and the birds-weighted `timed` share.

    `members` (taxon ids, the group first) in `taxonomy` for a group, None otherwise (a species
    with subspecies below it is always one), and `rollup_excluded` for a group that is alone or
    dominated by one taxon (`ROLLUP_MIN_BELOW`, `ROLLUP_MAX_BELOW`, `ROLLUP_MAX_SHARE`: share of
    all birds counted under it), whose own series is left as it was. Sums of all taxa no longer
    reconcile with `count.csv`, so the reconciliation uses the release taxa.
    """
    below = descendants(taxonomy.set_index("taxon_id")["parent_taxon_id"])
    birds = days.groupby("taxon_id")["count"].sum()
    rank = taxonomy.set_index("taxon_id")["taxon_rank"]
    members, excluded = {}, set()
    for taxon_id, kids in below.items():
        own = birds.reindex([taxon_id, *kids]).fillna(0)
        if rank[taxon_id] == "species":  # its subspecies are always part of it
            members[taxon_id] = [taxon_id, *sorted(kids)]
        elif (
            not ROLLUP_MIN_BELOW <= len(kids) <= ROLLUP_MAX_BELOW
            or own.max() > ROLLUP_MAX_SHARE * own.sum()
        ):
            excluded.add(taxon_id)
        else:
            members[taxon_id] = [taxon_id, *sorted(kids)]
    d_parts = [days[~days["taxon_id"].isin(members)]]
    h_parts = [hourly[~hourly["taxon_id"].isin(members)]]
    for taxon_id, ids in members.items():
        dm = days[days["taxon_id"].isin(ids)].assign(timed_birds=lambda x: x["timed"] * x["count"])
        g = dm.groupby("date")
        dc = g[["count", "reverse", "local", "timed_birds"]].sum(min_count=1)
        dc["qualifiers"] = g["qualifiers"].agg(lambda x: "".join(sorted(set("".join(x)))))
        dc["timed"] = (dc["timed_birds"].fillna(0) / dc["count"]).where(dc["count"] > 0)
        d_parts.append(dc.drop(columns="timed_birds").reset_index().assign(taxon_id=taxon_id))
        hm = hourly[hourly["taxon_id"].isin(ids)]
        h_parts.append(
            hm.groupby(["date", "hour"])["count"].sum().reset_index().assign(taxon_id=taxon_id)
        )
    t = taxonomy.assign(
        members=taxonomy["taxon_id"].map(members).astype(object),
        rollup_excluded=taxonomy["taxon_id"].isin(excluded),
    )
    return (
        t,
        pd.concat(d_parts, ignore_index=True)[days.columns],
        pd.concat(h_parts, ignore_index=True)[hourly.columns],
    )


def has_members(members) -> bool:
    """Whether a taxon's `members` value is a roll-up's list (not None or NaN)."""
    return isinstance(members, list)


# --- assembly ---------------------------------------------------------------


def release_counts(count: pd.DataFrame, surveys: pd.DataFrame, taxonomy: pd.DataFrame):
    """`count.csv` parsed (`counts.parse_counts`, without entry times), plus its own
    `raw_datetime`, the solar hours it is timed to (`timed_span`) and its `timed_count`."""
    c = parse_counts(count, surveys, taxonomy)
    c["raw_datetime"] = count["datetime"].to_numpy()
    c[["solar_start", "solar_end"]] = timed_span(c, surveys)
    c["timed_count"] = c["count"].where(c["solar_start"].notna(), 0).fillna(0)
    return c


def read_release(data_dir: str):
    """`(surveys, counts, taxonomy, metadata)` from `<data_dir>/count/dataset/`.

    The release tables only: `entry_times.csv`, the model's extra, is not read. `counts` is
    `release_counts`.
    """
    folder = os.path.join(data_dir, DATASET_DIR)
    surveys = parse_surveys(pd.read_csv(os.path.join(folder, "survey.csv"), low_memory=False))
    taxonomy = pd.read_csv(os.path.join(folder, "taxonomy.csv"))
    count = pd.read_csv(os.path.join(folder, "count.csv"), low_memory=False)
    c = release_counts(count, surveys, taxonomy)
    path = os.path.join(folder, METADATA_FILE)
    metadata = json.load(open(path)) if os.path.exists(path) else {}
    return surveys, c, taxonomy, metadata


def partial_years(days: pd.DataFrame) -> list[int]:
    """The last year, if its last record is before the window's end (a season in progress)."""
    last = days["date"].max()
    (m1, d1) = WINDOW[1]
    return [int(last.year)] if (last.month, last.day) < (m1, d1) else []


def manifest(metadata: dict, days: pd.DataFrame, effort: pd.DataFrame, git_sha: str) -> dict:
    """When and from what the export was built, and the definitions it uses."""
    return {
        "definitions_version": DEFINITIONS_VERSION,
        "built_at": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        "git_sha": git_sha,
        "dataset_built_at": metadata.get("built_at"),
        "dataset_git_sha": metadata.get("git_sha"),
        "first_date": effort["date"].min(),
        "last_date": max(effort["date"].max(), days["date"].max()),
        "partial_years": partial_years(days),
        "window": [f"{m:02d}-{d:02d}" for m, d in WINDOW],
        "tiers": {
            "full_min_days": TIER_FULL_MIN_DAYS,
            "full_min_years": TIER_FULL_MIN_YEARS,
            "rare_max_days": TIER_RARE_MAX_DAYS,
        },
        "start_years": [SYSTEMATIC_FROM, ALL_GROUPS_FROM],
        "timezone": TIMEZONE,
    }


# --- JSON -------------------------------------------------------------------


def _value(v):
    if isinstance(v, pd.DataFrame):
        return records(v)
    if isinstance(v, (pd.Index, pd.Series)):
        return _value(v.to_numpy())
    if isinstance(v, np.bool_):
        return bool(v)
    if v is None or (isinstance(v, float) and np.isnan(v)) or v is pd.NaT:
        return None
    if isinstance(v, pd.Timestamp):
        return v.strftime("%Y-%m-%d") if v == v.normalize() else v.strftime("%H:%M")
    if isinstance(v, (np.integer,)):
        return int(v)
    if isinstance(v, (float, np.floating)):
        f = round(float(v), DECIMALS)
        return int(f) if f.is_integer() else f
    if isinstance(v, np.ndarray):
        return [_value(x) for x in v]
    return v


def columns(df: pd.DataFrame) -> dict:
    """A frame as `{column: [values]}`, with dates as `YYYY-MM-DD`, times as `HH:MM`, NaN as
    null."""
    return {c: [_value(v) for v in df[c]] for c in df.columns}


def records(df: pd.DataFrame) -> list[dict]:
    """A frame as a list of `{column: value}`, same conversions as `columns`."""
    return [{c: _value(v) for c, v in row.items()} for row in df.to_dict("records")]


def _plain(obj):
    """Containers walked, every float through `_value`: `json` writes a float (numpy's too) as is,
    NaN included, which is not JSON, and never calls `default` for it."""
    if isinstance(obj, dict):
        return {k: _plain(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_plain(v) for v in obj]
    if isinstance(obj, float):
        return None if not np.isfinite(obj) else _value(obj)
    return obj


def dumps(obj) -> str:
    """Strict JSON (no NaN or Infinity: null instead), dates and frames as `_value` writes them."""
    return json.dumps(
        _plain(obj), ensure_ascii=False, separators=(",", ":"), default=_value, allow_nan=False
    )
