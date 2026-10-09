"""Check source annual totals against the current release for website editing.

Run from the repository root: python3 scripts/accounts/check_website_account_totals.py Totals cover
whole local calendar years and main-direction numeric counts, without effort correction. They are
editorial checks, not population trend estimates.
"""

import csv
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from zoneinfo import ZoneInfo

root = Path(__file__).resolve().parents[2]
destination = root / "content/accounts/review"
destination.mkdir(parents=True, exist_ok=True)
surveys = {
    r["survey_id"]: r["datetime"]
    for r in csv.DictReader((root / "data/count/dataset/survey.csv").open(encoding="utf-8-sig"))
}
totals = defaultdict(int)
presence = set()
for row in csv.DictReader((root / "data/count/dataset/count.csv").open(encoding="utf-8-sig")):
    if row["count_category"] != "normal":
        continue
    timing = (row["datetime"] or surveys[row["survey_id"]]).split("/")[0]
    date = (
        timing
        if len(timing) == 10
        else datetime.fromisoformat(timing.replace("Z", "+00:00"))
        .astimezone(ZoneInfo("Europe/Paris"))
        .date()
        .isoformat()
    )
    key = (date[:4], row["taxon_id"])
    if row["count"]:
        presence.add(key)
        totals[key] += int(row["count"])

names = {
    r["avibase_id"]: r["taxon_name_original"]
    for r in csv.DictReader((destination / "reference/source_taxa.csv").open())
    if r["source"] == "historical"
}
published = {
    (r["year"], r["species"]): int(r["count"])
    for r in csv.DictReader((destination / "reference/annual-totals.csv").open())
    if r["count"]
}
fields = [
    "year",
    "key",
    "species",
    "dataset_total",
    "published_total",
    "difference",
    "rank_since_1993",
    "rank_years",
    "scope",
]
with (destination / "data-checks.csv").open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=fields)
    writer.writeheader()
    for row in csv.DictReader(
        (root / "data/count/dataset/report_text.csv").open(encoding="utf-8-sig")
    ):
        if row["category"] != "species":
            continue
        key = (row["year"], row["key"])
        name = names.get(row["key"], "")
        value = totals[key] if key in presence else None
        reference = published.get((row["year"], name))
        comparison = [
            totals[k]
            for k in presence
            if k[1] == row["key"] and 1993 <= int(k[0]) <= int(row["year"])
        ]
        writer.writerow(
            dict(
                year=row["year"],
                key=row["key"],
                species=name,
                dataset_total=value,
                published_total=reference,
                difference=value - reference
                if value is not None and reference is not None
                else "",
                rank_since_1993=1 + sum(v > value for v in comparison)
                if value is not None
                else "",
                rank_years=len(comparison),
                scope="calendar year; normal direction; numeric values; no effort correction; absent entries unknown",
            )
        )

print(f"Wrote {destination / 'data-checks.csv'}")
