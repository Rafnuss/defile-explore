"""How far each claim a taxon's trend makes can be trusted: the `reliability` block.

Defined here, applied by defileViz: the page filters on the classes and can say why from the
reasons. Three claims, each classed `show`, `caveat` (shown, with a sentence on why) or `hide`:

- `totals`: the gap-filled annual totals. Judged mainly by the benchmark (`benchmark.py`): how far
  the fill misses birds actually counted (`FILL_ERROR`, median |log(estimate / counted)|), and how
  often its 80% intervals hold them (`MIN_COVER80`); then by how much of each year is filled and
  how wide the intervals are.
- `trend`: the smooth total over the years (the long-term change). Needs the totals, and a smooth
  whose 95% band is not too wide (`SMOOTH_BAND`, median over years of q97.5 / q2.5).
- `season`: the passage dates. Judged by the width of the median passage date's 80% band
  (`PASSAGE_BAND`, days) and by whether the passage runs past the counted days.

Every reason is a code with the level it sets; a claim takes its worst reason's level. `elements`
applies the classes to what a page draws (`ELEMENTS`: each field of the export that rests on a
claim, and the claim): defileViz shows, caveats or hides each one by it, with no rule of its own.
Fields not listed (the birds counted, the empirical season, the chances, the daytime, age and sex,
the records) are counts, shown whatever the trend's classes.
`estimated_years` lists the years mostly filled (`ESTIMATED_SHARE`): a page can mark them. The
numbers behind the reasons are in `inputs`. The recent-years test of the benchmark is not used: a
season's level is not predictable from the earlier ones, which is the forecast's problem, not
Explore's.
"""

import numpy as np
import pandas as pd

METHOD = "reliability@2"
CLASSES = ("show", "caveat", "hide")  # in order of severity
FILL_ERROR = (0.15, 0.4)  # benchmark gap error: caveat above the first, hide above the second
MIN_TRIALS = 6  # gap trials for the benchmark to count as a test
MIN_COVER80 = 0.6  # ... and its 80% intervals to hold the truth at least this often (9+ trials)
MIN_COVER_TRIALS = 9
MIN_OBSERVED = 0.6  # median share of a year's total actually counted
INTERVAL_RATIO = 2.0  # median q90 / q10 of the gap-filled totals (`wide_intervals`)
SMOOTH_BAND = (4.0, 20.0)  # median q97.5 / q2.5 of the smooth total: caveat, hide
PASSAGE_BAND = (10.0, 30.0)  # median width of the passage date's 80% band, days: caveat, hide
ESTIMATED_SHARE = 0.5  # a year with less than this share counted is mostly an estimate
ELEMENTS = {  # field of the species file: the claim it rests on
    "trend.annual.total": "totals",  # the gap-filled season totals, with q10/q90
    "key_numbers.typical_season": "totals",
    "trend.annual.smooth": "trend",  # the smooth total, with its 95% band
    "key_numbers.trend": "trend",
    "trend.passage": "season",  # each year's modelled passage dates
    "trend.passage_q": "season",  # the smooth passage dates (lines on the phenology)
    "key_numbers.passage": "season",  # when its `source` is `smooth`
}


def _ratio(hi: pd.Series, lo: pd.Series) -> float:
    return float(np.median(hi / lo.clip(lower=1)))


def inputs(trend: dict, bench: dict | None, window: dict, profile: str) -> dict:
    """The numbers the rules read."""
    a = pd.DataFrame(trend["annual"])
    p = pd.DataFrame(trend["passage"])
    gap = (bench or {}).get("gap") or {}
    return {
        "gap_error": gap.get("abs_log_err"),
        "gap_cover80": gap.get("cover80"),
        "gap_trials": int(gap.get("n", 0)),
        "observed_share": float(a["observed_share"].median()),
        "interval_ratio": _ratio(a["q90"], a["q10"]),
        "smooth_band": _ratio(a["smooth_q97.5"], a["smooth_q2.5"]),
        "passage_band": float((p["hi"] - p["lo"]).median()),
        "beyond_counting": list(window["beyond_counting"]),
        "profile": profile,
    }


def _level(lo: float, hi: float, x) -> str | None:
    if x is None or np.isnan(x) or x <= lo:
        return None
    return "hide" if x > hi else "caveat"


def _claim(reasons: dict) -> dict:
    """`{class, reasons}` from `{code: level or None}`."""
    kept = {k: v for k, v in reasons.items() if v}
    worst = max(kept.values(), key=CLASSES.index, default="show")
    return {"class": worst, "reasons": [{"code": k, "level": v} for k, v in kept.items()]}


def classes(x: dict) -> dict:
    """The three claims' classes from `inputs`."""
    tested = x["gap_error"] is not None and x["gap_trials"] >= MIN_TRIALS
    cover = x["gap_cover80"]
    totals = _claim(
        {
            "untested": None if tested else "caveat",
            "fill_error": _level(*FILL_ERROR, x["gap_error"]) if tested else None,
            "intervals_too_narrow": (
                "caveat"
                if x["gap_trials"] >= MIN_COVER_TRIALS
                and cover is not None
                and cover < MIN_COVER80
                else None
            ),
            "mostly_filled": "caveat" if x["observed_share"] < MIN_OBSERVED else None,
            "wide_intervals": "caveat" if x["interval_ratio"] > INTERVAL_RATIO else None,
            # the benchmark measures a borrowed profile's error; untested, it is a doubt
            "borrowed_profile": "caveat" if x["profile"] != "own" and not tested else None,
        }
    )
    trend = _claim(
        {
            "totals": None if totals["class"] == "show" else totals["class"],
            "smooth_band": _level(*SMOOTH_BAND, x["smooth_band"]),
        }
    )
    season = _claim(
        {
            "passage_band": _level(*PASSAGE_BAND, x["passage_band"]),
            "beyond_counting": "caveat" if x["beyond_counting"] else None,
        }
    )
    return {"totals": totals, "trend": trend, "season": season}


def reliability_block(trend: dict | None, bench, window: dict, profile: str) -> dict | None:
    """The `reliability` block of a taxon with a trend; None without one."""
    if trend is None:
        return None
    x = inputs(trend, bench, window, profile)
    a = pd.DataFrame(trend["annual"])
    return {
        "method": METHOD,
        **(c := classes(x)),
        "elements": {k: c[claim]["class"] for k, claim in ELEMENTS.items()},
        "estimated_years": [int(y) for y in a.loc[a["observed_share"] < ESTIMATED_SHARE, "year"]],
        "inputs": x,
    }
