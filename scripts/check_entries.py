#!/usr/bin/env python3
"""List the entries of the release that look like recording errors, for the editors to correct at
the source (defile-dataset or Trektellen); nothing here changes the build.

Checks, each a row per suspect count (or survey):

- `block_edge`: a survey of at least `BLOCK_HOURS` whose birds timed to a point are bunched in its
  first or last `EDGE_MINUTES` (at least `EDGE_SHARE` of them, `EDGE_BIRDS` or more): a block's
  total entered at one time, which then reads as a flock at that minute;
- `outside_survey`: a count timed more than `OUTSIDE_MINUTES` before its survey's start or after
  its end;
- `day_shifted`: a day with at least `SHIFTED_SHARE` of its timed birds in the dark: its clock
  looks an hour off (seen on days summer time ends), every count of it listed;
- `dark`: other counts timed more than `DARK_MINUTES` before civil dawn or after civil dusk (sun
  below -6 deg): a typo in the hour, a block's total entered late, or a real night record;
- `duplicate`: the same taxon, category, count and time entered twice in a survey;
- `survey_overlap`: two surveys of the same day overlapping by more than `OUTSIDE_MINUTES`, so the
  effort of those minutes is counted twice.

Writes `logs/qa/explore/entry_errors.csv` and prints a summary.

Usage:
    python scripts/check_entries.py [--data-dir data]
"""

import argparse
import os

import pandas as pd

from defile_explore import export as E
from defile_explore.release import civil_twilight

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BLOCK_HOURS = 2.0
EDGE_MINUTES = 30
EDGE_SHARE = 0.5
EDGE_BIRDS = 20
OUTSIDE_MINUTES = 5
DARK_MINUTES = 15  # a record this close to twilight is plausible
SHIFTED_SHARE = 0.2  # a day with this share of its timed birds in the dark: its clock is off


def local(t: pd.Series) -> pd.Series:
    return E.to_local(t).dt.strftime("%Y-%m-%d %H:%M")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--data-dir", default=os.path.join(ROOT, "data"))
    ap.add_argument("--out", default=os.path.join(ROOT, "logs", "qa", "explore"))
    args = ap.parse_args(argv)

    surveys, counts, taxonomy, _ = E.read_release(args.data_dir)
    names = taxonomy.set_index("taxon_id")["english_name"]
    sv = surveys.set_index("survey_id")
    c = counts.assign(
        taxon=counts["taxon_id"].map(names),
        s_start=counts["survey_id"].map(sv["start"]),
        s_end=counts["survey_id"].map(sv["end"]),
    )
    hours = (c["s_end"] - c["s_start"]) / pd.Timedelta(hours=1)
    pt = c[c["datetime"].notna() & (c["count"] > 0)]
    rows = []

    def add(check, d, detail):
        rows.append(d.assign(check=check, detail=detail))

    # block_edge
    p = pt[hours[pt.index] >= BLOCK_HOURS]
    edge = pd.Timedelta(minutes=EDGE_MINUTES)
    p = p.assign(
        at_start=p["datetime"] <= p["s_start"] + edge, at_end=p["datetime"] >= p["s_end"] - edge
    )
    g = p.groupby("survey_id")
    birds = g["count"].sum()
    for side in ("at_start", "at_end"):
        share = p["count"].where(p[side], 0).groupby(p["survey_id"]).sum() / birds
        bad = share.index[(share >= EDGE_SHARE) & (birds >= EDGE_BIRDS)]
        d = p[p["survey_id"].isin(bad) & p[side]]
        add(
            "block_edge",
            d,
            d["survey_id"].map(
                lambda s, side=side: f"{share[s]:.0%} of the survey's {birds[s]:.0f} timed birds "
                f"in its {'first' if side == 'at_start' else 'last'} {EDGE_MINUTES} min"
            ),
        )

    # outside_survey
    tol = pd.Timedelta(minutes=OUTSIDE_MINUTES)
    d = pt[(pt["datetime"] < pt["s_start"] - tol) | (pt["datetime"] > pt["s_end"] + tol)]
    add("outside_survey", d, "timed outside its survey")

    # dark
    dawn, dusk = civil_twilight(pt["date"])
    margin = pd.Timedelta(minutes=DARK_MINUTES)
    dark = (pt["datetime"] < dawn - margin) | (pt["datetime"] >= dusk + margin)
    share = (
        pt["count"].where(dark, 0).groupby(pt["date"]).sum() / pt.groupby("date")["count"].sum()
    )
    shifted = pt["date"].map(share) >= SHIFTED_SHARE
    d = pt[dark & shifted]
    add(
        "day_shifted",
        d,
        f"{SHIFTED_SHARE:.0%}+ of the day's timed birds in the dark; civil dawn "
        + local(dawn[d.index]).str[-5:]
        + ", dusk "
        + local(dusk[d.index]).str[-5:],
    )
    d = pt[dark & ~shifted]
    add(
        "dark",
        d,
        "civil dawn " + local(dawn[d.index]).str[-5:] + ", dusk " + local(dusk[d.index]).str[-5:],
    )

    # duplicate
    key = ["survey_id", "taxon_id", "count_category", "count", "datetime", "age", "sex", "plumage"]
    key += ["count_estimation", "remark"]
    d = pt[pt.duplicated(key, keep=False) & pt["count"].notna()]
    add("duplicate", d, "same taxon, category, count, time, age, sex and remark in the survey")

    # survey_overlap
    s = surveys.sort_values(["date", "start"])
    prev_end = s.groupby("date")["end"].shift()
    prev_id = s.groupby("date")["survey_id"].shift()
    over = s[(prev_end - s["start"]) > tol]
    survey_rows = over.assign(
        check="survey_overlap",
        detail="overlaps "
        + prev_id[over.index].map(sv["source_survey_id"]).astype(str)
        + " until "
        + local(prev_end[over.index]).str[-5:],
    )

    cols = ["check", "date", "time", "taxon", "count", "count_category", "source"]
    cols += ["source_count_id", "trektellen_data_id", "source_survey_id", "survey", "detail"]
    out = pd.concat(rows, ignore_index=True)
    out = out.assign(
        time=local(out["datetime"]).str[-5:],
        source_survey_id=out["survey_id"].map(sv["source_survey_id"]),
        survey=local(out["s_start"]) + "-" + local(out["s_end"]).str[-5:],
    )
    survey_rows = survey_rows.assign(
        time="", survey=local(survey_rows["start"]) + "-" + local(survey_rows["end"]).str[-5:]
    )
    out = pd.concat([out, survey_rows], ignore_index=True).reindex(columns=cols)
    out["date"] = pd.to_datetime(out["date"]).dt.strftime("%Y-%m-%d")
    out = out.sort_values(["check", "date", "time"])
    os.makedirs(args.out, exist_ok=True)
    path = os.path.join(args.out, "entry_errors.csv")
    out.to_csv(path, index=False)
    summary = out.groupby("check").agg(
        rows=("date", "size"), days=("date", "nunique"), birds=("count", "sum")
    )
    print(summary.to_string())
    print(f"-> {path}")


if __name__ == "__main__":
    main()
