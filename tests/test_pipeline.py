"""Tests for the per-taxon blocks (season, daytime, demography, settings, key numbers), on small
synthetic tables."""

import numpy as np
import pandas as pd
import pytest

from defile_explore import accounts as A
from defile_explore import daytime as Y
from defile_explore import demography as G
from defile_explore import pipeline as L
from defile_explore import remarks as M
from defile_explore import season as S
from defile_explore import settings as X
from defile_explore import window as W
from defile_explore.profile import PROFILE_DOY


def frame(rates: dict, c: float = 1.0) -> pd.DataFrame:
    """A `model_frame`-like table: {year: [birds per day over doy 200..]}; NaN = not counted."""
    rows = []
    for year, birds in rates.items():
        for i, b in enumerate(birds):
            counted = not np.isnan(b)
            rows.append(
                {
                    "year": year,
                    "doy": 200 + i,
                    "c": c if counted else 0.0,
                    "y": b if counted else 0,
                }
            )
    return pd.DataFrame(rows)


# --- season ------------------------------------------------------------------------------------


def test_passage_quantiles_of_a_flat_season_are_evenly_spaced():
    q = S.cumulative_quantiles(np.arange(100), np.ones(100))
    assert q == pytest.approx([9, 49, 89], abs=1)


def test_daily_rates_divide_by_coverage_and_fill_gaps_between_counted_days():
    f = frame({2020: [10, np.nan, 30]}, c=0.5)
    rate, filled = S.daily_rates(f)
    assert rate.loc[2020].tolist()[0] == 20 and np.isnan(rate.loc[2020].tolist()[1])
    assert filled.loc[2020].tolist() == [20, 40, 60]


def test_a_day_below_the_coverage_threshold_is_not_counted():
    f = frame({2020: [10, 10]}, c=0.3)
    rate, _ = S.daily_rates(f)
    assert rate.isna().all().all()


def test_season_shares_sum_to_one_where_every_day_is_counted():
    b = S.season_block(frame({2020: [1, 2, 3, 4], 2021: [4, 3, 2, 1]}), 2021)
    assert np.nansum(b["share"], axis=1) == pytest.approx([1, 1])
    assert b["passage"]["counted"].tolist() == [1, 1]


def test_days_beyond_the_model_window_are_shown_but_not_in_the_total():
    f = frame({2020: [10, 10, np.nan, np.nan, 50]})  # doy 200..204; 202-203 not counted
    b = S.season_block(f, 2020, model_window=(200, 201))
    assert b["doy"].tolist() == [200, 201, 202, 203, 204]
    assert b["share"][0].tolist()[:2] == [0.5, 0.5] and b["share"][0][4] == 2.5
    assert b["passage"]["counted"].tolist() == [1]


def test_chances_count_days_reaching_each_threshold():
    birds = [0, 5, 50, 500] * 5  # 20 days: a quarter at each level
    ch = S.chances(frame({2020: birds}), 2020)
    days = np.asarray(ch["days"])
    at1 = (np.asarray(ch["at_least_1"]) * days).sum() / days.sum()
    assert at1 == pytest.approx(0.75)
    assert "at_least_1000" not in ch


# --- window ------------------------------------------------------------------------------------


def season_frame(rate, counted_share, years: int = 10) -> pd.DataFrame:
    """Season days 182..336 over `years` years: `rate(doy)` birds on counted days, each day counted
    in the first `counted_share(doy)` share of the years."""
    rows = []
    for i in range(years):
        for doy in range(182, 337):
            counted = i < round(counted_share(doy) * years)
            rows.append(
                {"year": 2000 + i, "doy": doy, "c": float(counted), "y": rate(doy) * counted}
            )
    return pd.DataFrame(rows)


def test_a_passage_inside_the_default_window_keeps_it():
    lo, hi = W.default_window()
    f = season_frame(lambda d: 100.0 if 240 <= d <= 280 else 0.0, lambda d: 1.0)
    w = W.windows(f)
    assert w["model"] == (lo, hi) and w["view"] == (lo, hi) and w["beyond_counting"] == []


def test_a_late_passage_extends_to_where_enough_years_were_counted():
    lo, hi = W.default_window()  # 18 Nov = 322
    share = lambda d: 1.0 if d <= hi else (0.5 if d <= 328 else 0.2)  # noqa: E731
    w = W.windows(season_frame(lambda d: 100.0 if d >= 280 else 0.0, share))
    assert hi < w["model"][1] <= 328 + W.ENVELOPE_SMOOTH // 2  # modelled while 1/3 counted
    assert w["view"][1] == 336  # shown while 10% counted
    assert w["model"][0] == lo and w["beyond_counting"] == ["late"]


def test_the_default_window_is_never_cut():
    lo, hi = W.default_window()
    share = lambda d: 0.0 if d < lo + 10 else 1.0  # noqa: E731  (start of the window not counted)
    w = W.windows(season_frame(lambda d: 10.0, share))
    assert w["model"][0] <= lo and w["model"][1] >= hi


def test_a_window_override_is_read_as_dates(tmp_path):
    text = "t1: {model_window: ['07-18', '11-25'], reason: test}"
    o = X.load_overrides(write_overrides(tmp_path, text))
    rule = {"model": (199, 322), "view": (199, 330)}
    s = X.resolve_windows("t1", rule, o)
    assert W.as_dates(s["model_window"].value) == ["07-18", "11-25"]
    assert s["view_window"].value == (199, 330) and s["view_window"].source == "rule"


# --- daytime -----------------------------------------------------------------------------------


def test_hour_rates_divide_birds_by_the_hours_counted():
    dates = pd.to_datetime(["2020-09-01", "2020-09-02"])
    effort = np.zeros((2, 24))
    effort[0, 10:12] = 1  # day 1: 10-12 counted
    effort[1, 10] = 1  # day 2: 10 only
    sample = pd.DataFrame({"date": dates, "hourly": list(effort)})
    hourly = pd.DataFrame({"date": dates[[0, 0, 1]], "hour": [10, 11, 10], "count": [4, 6, 2]})
    birds, hours = Y.hour_rates(hourly, sample, pd.Series([0, 0]))
    rate = birds.iloc[0] / hours.iloc[0]
    assert rate[10] == 3 and rate[11] == 6


def timed(n_days: int, early_hour: int, late_hour: int):
    """`n_days` early (doy 200..) and `n_days` late (doy 300..) timed days, counted 6-18 h, the
    birds at one hour, with a little noise one hour later."""
    rng = np.random.default_rng(1)
    hours = np.zeros((2 * n_days, 24))
    hours[:, 6:18] = 1
    birds = np.zeros((2 * n_days, 24))
    for i in range(2 * n_days):
        h = early_hour if i < n_days else late_hour
        birds[i, h], birds[i, h + 1] = rng.integers(20, 40), rng.integers(0, 5)
    doy = np.r_[200 + np.arange(n_days), 300 + np.arange(n_days)]
    return birds, hours, doy


def test_a_clear_change_in_the_passage_hour_shows_the_matrix():
    c = Y.seasonal_change(*timed(20, 8, 11), [200, 250, 300, 350])
    assert c["show"] and c["shift"] == pytest.approx(3, abs=0.2) and c["lo"] > 0


def test_the_same_passage_hour_all_season_shows_one_histogram():
    c = Y.seasonal_change(*timed(20, 9, 9), [200, 250, 300, 350])
    assert not c["show"] and abs(c["shift"]) < 0.3


def test_too_few_days_are_not_tested():
    c = Y.seasonal_change(*timed(3, 8, 12), [200, 250, 300, 350])
    assert not c["show"] and "shift" not in c


# --- demography --------------------------------------------------------------------------------


def demo_rows(year: int, codes: dict, field: str = "age") -> pd.DataFrame:
    return pd.DataFrame(
        {
            "date": pd.Timestamp(f"{year}-09-15"),
            "count": list(codes.values()),
            "age": list(codes) if field == "age" else None,
            "sex": list(codes) if field == "sex" else None,
        }
    )


def test_wilson_interval_matches_the_textbook_value():
    lo, hi = G.wilson(5, 10)
    assert (round(float(lo), 3), round(float(hi), 3)) == (0.237, 0.763)


def test_non_adult_codes_are_one_class():
    rows = pd.concat([demo_rows(y, {"1": 10, "I": 20, "A": 30}) for y in (2020, 2021, 2022)])
    counted = pd.Series(100.0, index=[2020, 2021, 2022])
    b = G.age_block(rows, counted, [2020, 2021, 2022])
    assert b["share"].tolist() == [0.5, 0.5, 0.5]
    assert b["overall"]["n"] == 180
    assert b["overall"]["counted"] == 300


def test_sex_is_the_male_share_against_female_types():
    rows = pd.concat(
        [demo_rows(y, {"M": 30, "F": 10, "FC": 20}, "sex") for y in (2020, 2021, 2022)]
    )
    b = G.sex_block(rows, pd.Series(100.0, index=[2020, 2021, 2022]), [2020, 2021, 2022])
    assert b["classes"] == ["male", "female_type"]
    assert b["share"].tolist() == [0.5, 0.5, 0.5]
    assert b["codes"]["FC"].tolist() == [20, 20, 20]
    assert b["timing"]["male_birds"] == 90


def test_years_need_enough_birds_aged_and_enough_years():
    rows = pd.concat([demo_rows(2020, {"A": 30, "1": 20}), demo_rows(2021, {"A": 5})])
    counted = pd.Series([100.0, 100.0], index=[2020, 2021])
    assert G.usable_years(rows, "age", counted) == []  # one usable year, MIN_YEARS needed
    rows = pd.concat([demo_rows(y, {"A": 30, "1": 20}) for y in (2020, 2021, 2022)])
    assert G.usable_years(rows, "age", pd.Series(100.0, index=[2020, 2021, 2022])) == [
        2020,
        2021,
        2022,
    ]


def test_a_year_with_only_one_class_recorded_is_not_usable():
    # juveniles tagged, adults left blank: the share among aged birds would be ~100%
    years = (2020, 2021, 2022)
    rows = pd.concat([demo_rows(y, {"J": 400, "A": 2}) for y in years])
    assert G.usable_years(rows, "age", pd.Series(1000.0, index=list(years))) == []


# --- settings ----------------------------------------------------------------------------------


def write_overrides(tmp_path, text: str) -> str:
    path = tmp_path / "overrides.yaml"
    path.write_text(text)
    return str(path)


def test_an_override_without_a_reason_is_refused(tmp_path):
    with pytest.raises(ValueError, match="reason"):
        X.load_overrides(write_overrides(tmp_path, "t1: {start_year: 2007}"))


def test_an_unknown_setting_is_refused(tmp_path):
    with pytest.raises(ValueError, match="unknown"):
        X.load_overrides(write_overrides(tmp_path, "t1: {window: 3, reason: x}"))


def test_links_need_no_reason(tmp_path):
    o = X.load_overrides(write_overrides(tmp_path, "t1: {links: {vogelwarte: red-kite}}"))
    assert o["t1"]["links"]["vogelwarte"] == "red-kite"


def taxon(**kw) -> dict:
    return {
        "taxon_id": "t1",
        "tier": "full",
        "taxon_rank": "species",
        "start_year": 1993,
        "scientific_name": "Milvus milvus",
        "ebird_code": "redkit1",
        "trektellen_species_id": 101.0,
    } | kw


def test_overrides_win_over_rules_and_say_so():
    rows = pd.concat([demo_rows(y, {"M": 3}, "sex") for y in (2001, 2005, 2030)])
    o = {"t1": {"start_year": 2007, "sex_years": {"from": 2002}, "reason": "test"}}
    s, links = X.resolve(taxon(), rows, pd.Series(dtype=float), o, last_year=2025)
    assert s["start_year"].as_dict() == {"value": 2007, "source": "override", "reason": "test"}
    assert s["sex_years"].value == [2005]  # 2030 is after the last complete season
    assert s["trend"].source == "rule" and s["trend"].value is True
    assert links["ebba2"].endswith("/Milvus-milvus/ebba2/occurrence/")
    assert "/2422/101/" in links["trektellen"]


def test_no_trend_below_the_full_tier():
    s, _ = X.resolve(taxon(tier="short"), demo_rows(2020, {}), pd.Series(dtype=float), {}, 2025)
    assert s["trend"].value is False


# --- key numbers -------------------------------------------------------------------------------


def test_best_hours_are_the_fewest_holding_the_share():
    profile = np.zeros((PROFILE_DOY[1] - PROFILE_DOY[0] + 1, 24))
    profile[:, 10], profile[:, 11], profile[:, 12] = 0.5, 0.3, 0.2
    assert L.best_hours(profile, 250) == {"from": 10, "to": 12, "share": 0.8}


# --- accounts ----------------------------------------------------------------------------------


def write_accounts(tmp_path, sections: str, years: str) -> str:
    (tmp_path / A.SECTIONS_FILE).write_text("key\tsection\ttext_fr\ttext_en\n" + sections)
    (tmp_path / A.YEARS_FILE).write_text("year\tkey\ttext_fr\ttext_en\n" + years)
    return str(tmp_path)


def test_the_authored_accounts_load_and_check():
    a = A.load_accounts()
    assert len(a["sections"]) and len(a["years"])
    assert set(a["sections"]["section"]) <= set(A.SECTIONS)


def test_an_account_missing_one_language_is_refused(tmp_path):
    folder = write_accounts(tmp_path, "t1\tpassage\tTexte\t \n", "")
    with pytest.raises(ValueError, match="one language"):
        A.load_accounts(folder)


def test_an_unknown_section_or_taxon_is_refused(tmp_path):
    with pytest.raises(ValueError, match="unknown section"):
        A.load_accounts(write_accounts(tmp_path, "t1\thabitat\tT\tT\n", ""))
    folder = write_accounts(tmp_path, "", "2020\tt9\tT\tT\n")
    with pytest.raises(ValueError, match="unknown taxa"):
        A.load_accounts(folder, taxon_ids=["t1"])


def test_accounts_block_orders_sections_and_years(tmp_path):
    sections = "t1\tparticularites\tP\tP\nt1\tpassage\tA\tA\n"
    years = "2020\tt1\tx\tx\n2024\tt1\ty\ty\n2026\tt1\tz\tz\n2021\tt2\tw\tw\n"
    a = A.load_accounts(write_accounts(tmp_path, sections, years))
    b = A.accounts_block(A.taxon_rows(a, "t1"), last_year=2025)
    assert list(b["general"]) == ["passage", "particularites"]
    assert [r["year"] for r in b["years"]] == [2024, 2020]  # newest first, none after 2025
    assert A.accounts_block(A.taxon_rows(a, "t3"), 2025) is None


# --- remarks -----------------------------------------------------------------------------------


def test_a_remark_is_split_into_paragraphs_with_their_source():
    text = (
        "details: \n\n[Rapport annuel 2020; Contexte journalier 2020-10-14, Milan royal] Un pic."
    )
    assert M.split(text) == [("Rapport annuel 2020", "Un pic.")]
    assert M.split("detail: 2x mâles adultes") == [(None, "2x mâles adultes")]
    assert M.split("An explicit 'no species' entry was recorded in the source.") == []


def test_a_day_lists_survey_notes_first_and_each_text_once():
    d = pd.Timestamp("2020-10-14")
    remarks = pd.DataFrame(
        [
            ("t", d, "count", "R", "Un pic."),
            ("t", d, "count", "R", "Un pic."),
            (None, d, "survey", None, "Quelle journée !"),
        ],  # fmt: skip
        columns=M.COLUMNS,
    )
    assert [n["text"] for n in M.notes_of(remarks, d)] == ["Quelle journée !", "Un pic."]
