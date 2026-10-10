"""Scientific accounting checks for the independent hourly weather pilot."""

import pandas as pd

from defile_explore import hourly as H


def test_split_preserves_empty_partial_hours_and_excludes_unobserved_statuses():
    surveys = pd.DataFrame(
        {
            "survey_id": ["observed", "rain", "incomplete"],
            "recording_era": ["trektellen"] * 3,
            "start": pd.to_datetime(["2023-09-01T09:20Z"] * 3),
            "end": pd.to_datetime(["2023-09-01T11:10Z"] * 3),
            "survey_complete": [True, True, False],
            "weather_stop": [None, True, None],
        }
    )
    hours = H.survey_hours(surveys)
    assert hours["survey_id"].tolist() == ["observed"] * 3
    assert hours["exposure_hours"].tolist() == [40 / 60, 1, 10 / 60]
    assert hours["exposure_hours"].sum() == 110 / 60


def test_untimed_birds_never_make_false_species_zeros_and_all_birds_reconcile():
    surveys = pd.DataFrame(
        {
            "survey_id": ["observed"],
            "recording_era": ["trektellen"],
            "start": pd.to_datetime(["2023-09-01T09:20Z"]),
            "end": pd.to_datetime(["2023-09-01T11:10Z"]),
            "survey_complete": [True],
            "weather_stop": [None],
        }
    )
    counts = pd.DataFrame(
        {
            "count_id": ["a", "b", "c", "d", "e"],
            "survey_id": ["observed"] * 5,
            "taxon_id": ["kite", "kite", "buzzard", "buzzard", "buzzard"],
            "datetime": [
                "2023-09-01T09:30Z",
                "2023-09-01",
                "2023-09-01T10:00Z/2023-09-01T11:00Z",
                "2023-09-01T09:40Z/2023-09-01T10:10Z",
                "2023-09-01T11:10Z",
            ],
            "count": [7, 13, 5, 11, 2],
            "count_category": ["normal"] * 5,
        }
    )
    timing = H.count_timing(counts, surveys)
    assert timing["timing_status"].tolist() == [
        "point",
        "untimed",
        "single_hour_interval",
        "untimed",
        "outside_survey",
    ]
    assert timing.loc[timing["slot_start"].notna(), "count"].sum() == 12
    assert timing.loc[timing["slot_start"].isna(), "count"].sum() == 26
    taxa = pd.DataFrame(
        {"taxon_id": ["kite", "buzzard", "osprey"], "english_name": ["Kite", "Buzzard", "Osprey"]}
    )
    frame = H.species_hours(H.survey_hours(surveys), timing, taxa)
    assert frame.loc[frame["taxon_id"].isin(["kite", "buzzard"]), "count"].isna().all()
    assert frame.loc[frame["taxon_id"].eq("osprey"), "count"].eq(0).all()
    assert frame["timed_count"].sum() == 12


def test_inherited_short_interval_is_supported_and_presence_is_not_a_zero():
    surveys = pd.DataFrame(
        {
            "survey_id": ["short"],
            "recording_era": ["trektellen"],
            "start": pd.to_datetime(["2023-09-01T09:20Z"]),
            "end": pd.to_datetime(["2023-09-01T09:50Z"]),
            "survey_complete": [True],
            "weather_stop": [None],
        }
    )
    counts = pd.DataFrame(
        {
            "count_id": ["a", "b"],
            "survey_id": ["short"] * 2,
            "taxon_id": ["kite", "osprey"],
            "datetime": [None, None],
            "count": [4, None],
            "count_category": ["normal"] * 2,
        }
    )
    timing = H.count_timing(counts, surveys)
    assert timing["timing_status"].tolist() == ["single_hour_interval", "presence_only"]
    taxa = pd.DataFrame({"taxon_id": ["kite", "osprey"], "english_name": ["Kite", "Osprey"]})
    frame = H.species_hours(H.survey_hours(surveys), timing, taxa)
    assert frame.loc[frame["taxon_id"].eq("kite"), "count"].iloc[0] == 4
    assert frame.loc[frame["taxon_id"].eq("osprey"), "count"].isna().all()
