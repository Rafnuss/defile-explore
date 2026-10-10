"""Hourly Trektellen observations for the weather pilot, preserving unsupported timing."""

import pandas as pd

from defile_explore.release import MAIN_CATEGORY, TIMEZONE, TREKTELLEN_ERA


def survey_hours(surveys: pd.DataFrame) -> pd.DataFrame:
    """Split complete observed surveys at UTC hour boundaries, including empty/partial hours."""
    observed = surveys[
        surveys["recording_era"].eq(TREKTELLEN_ERA)
        & surveys["survey_complete"].eq(True)
        & ~surveys["weather_stop"].eq(True)
    ]
    rows = []
    for s in observed.itertuples():
        edges = [s.start, *pd.date_range(s.start.ceil("h"), s.end, freq="h"), s.end]
        for start, end in zip(edges[:-1], edges[1:]):
            if end > start:
                rows.append((s.survey_id, start.floor("h"), start, end))
    out = pd.DataFrame(rows, columns=["survey_id", "slot_start", "start", "end"])
    out["date"] = out["start"].dt.tz_convert(TIMEZONE).dt.date
    out["exposure_hours"] = (out["end"] - out["start"]) / pd.Timedelta(hours=1)
    return out


def count_timing(counts: pd.DataFrame, surveys: pd.DataFrame) -> pd.DataFrame:
    """Assign only supported point/single-hour intervals; retain every main-direction record."""
    out = counts[
        counts["survey_id"].isin(surveys["survey_id"]) & counts["count_category"].eq(MAIN_CATEGORY)
    ].merge(surveys[["survey_id", "start", "end"]], on="survey_id", validate="many_to_one")
    raw = out["datetime"].fillna("")
    point = raw.str.contains("T") & ~raw.str.contains("/")
    timestamp = pd.to_datetime(raw.where(point), utc=True, format="ISO8601")
    out["slot_start"] = pd.Series(pd.NaT, index=out.index, dtype="datetime64[ns, UTC]")
    out["timing_status"] = "untimed"
    inside = point & timestamp.ge(out["start"]) & timestamp.lt(out["end"])
    out.loc[inside, "slot_start"] = timestamp[inside].dt.floor("h")
    out.loc[inside, "timing_status"] = "point"
    out.loc[point & ~inside, "timing_status"] = "outside_survey"

    interval = raw.str.contains("/")
    bounds = raw.where(interval).str.split("/", expand=True).reindex(columns=[0, 1])
    start = pd.to_datetime(bounds[0], utc=True, format="ISO8601").where(~raw.eq(""), out["start"])
    end = pd.to_datetime(bounds[1], utc=True, format="ISO8601").where(~raw.eq(""), out["end"])
    supported = (
        (interval | raw.eq(""))
        & start.ge(out["start"])
        & end.le(out["end"])
        & start.dt.floor("h").eq((end - pd.Timedelta(nanoseconds=1)).dt.floor("h"))
    )
    out.loc[supported, "slot_start"] = start[supported].dt.floor("h")
    out.loc[supported, "timing_status"] = "single_hour_interval"
    out.loc[out["count"].isna(), ["slot_start", "timing_status"]] = [pd.NaT, "presence_only"]
    return out


def species_hours(hours: pd.DataFrame, timing: pd.DataFrame, taxa: pd.DataFrame) -> pd.DataFrame:
    """Species by effort grid, with zeros only where the species' hourly count is complete."""
    frame = hours.merge(taxa[["taxon_id", "english_name"]], how="cross")
    timed = timing[timing["slot_start"].notna()]
    keys = ["survey_id", "slot_start", "taxon_id"]
    totals = timed.groupby(keys, as_index=False)["count"].sum()
    frame = frame.merge(totals, on=keys, how="left", validate="one_to_one")
    frame["count"] = frame["count"].fillna(0)
    frame["timed_count"] = frame["count"]
    unresolved = timing[
        timing["slot_start"].isna()
        & (timing["count"].gt(0) | timing["timing_status"].eq("presence_only"))
    ]
    unknown = unresolved.groupby(["survey_id", "taxon_id"]).size().rename("unresolved_entries")
    frame = frame.merge(unknown, on=["survey_id", "taxon_id"], how="left")
    frame["unresolved_entries"] = frame["unresolved_entries"].fillna(0).astype(int)
    frame["hourly_complete"] = frame["unresolved_entries"].eq(0)
    frame["count"] = frame["count"].where(frame["hourly_complete"])
    return frame
