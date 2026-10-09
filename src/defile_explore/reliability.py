"""How far each claim a taxon's trend makes can be trusted: the `reliability` block.

Defined here, applied by defileViz: the page filters on the classes and can say why from the
reasons. Three claims, each classed `show`, `caveat` (shown, with a sentence on why) or `hide`:

- `totals`: the gap-filled annual totals. Judged mainly by the benchmark (`benchmark.py`): how far
  the fill misses birds actually counted (`FILL_ERROR`, median |log(estimate / counted)|), and how
  often its 80% intervals hold them (`MIN_COVER80`); then by how much of each year is filled and
  how wide the intervals are.
- `trend`: the smooth total over the years (the long-term change). Needs the totals, and a smooth
  whose 95% band is not too wide (`SMOOTH_BAND`, median over years of q97.5 / q2.5).
- `season`: the passage dates. Judged by the width of the 80% bands of the median passage date
  and of the 10% and 90% dates (the wider of the two), against `PASSAGE_BAND` (days); by the birds
  behind them (`SEASON_BIRDS`: the model's smooth timing keeps those bands narrow even for a few
  birds a season); and by whether much of the passage falls where counting thins out
  (`window.EDGE_SHARE`).

Every reason is a code with the level it sets; a claim takes its worst reason's level. `elements`
applies the classes to what a page draws (`ELEMENTS`: each field of the export that rests on a
claim, and the claim): defileViz shows, caveats or hides each one by it, with no rule of its own.
A field in `NO_CAVEAT` is hidden rather than caveated: lines drawn over the counts, where a caveat
has no room and a doubtful line reads as a finding. A field in `COUNTS_FALLBACK` (the season block,
from the trend's gap-filled days) is built from the counts instead where its claim is hidden
(`fills_season`), and then shown.
Fields not listed (the birds counted, the daytime, age and sex, the records) are counts, shown
whatever the trend's classes.
`estimated_years` lists the years mostly filled (`ESTIMATED_SHARE`): a page can mark them. The
numbers behind the reasons are in `inputs`. The recent-years test of the benchmark is not used: a
season's level is not predictable from the earlier ones, which is the forecast's problem, not
Explore's.
"""

import numpy as np
import pandas as pd

METHOD = "reliability@3"
CLASSES = ("show", "caveat", "hide")  # in order of severity
FILL_ERROR = (0.15, 0.4)  # benchmark gap error: caveat above the first, hide above the second
MIN_TRIALS = 6  # gap trials for the benchmark to count as a test
MIN_COVER80 = 0.6  # ... and its 80% intervals to hold the truth at least this often (9+ trials)
MIN_COVER_TRIALS = 9
MIN_OBSERVED = 0.6  # median share of a year's total actually counted
INTERVAL_RATIO = 2.0  # median q90 / q10 of the gap-filled totals (`wide_intervals`)
SMOOTH_BAND = (4.0, 20.0)  # median q97.5 / q2.5 of the smooth total: caveat, hide
PASSAGE_BAND = (10.0, 30.0)  # median width of the passage date's 80% band, days: caveat, hide
SEASON_BIRDS = 50  # median birds counted per season, below which the dates are a caveat
ESTIMATED_SHARE = 0.5  # a year with less than this share counted is mostly an estimate
ELEMENTS = {  # field of the species file: the claim it rests on
    "trend.annual.total": "totals",  # the gap-filled season totals, with q10/q90
    "trend.days.total": "totals",  # each counted day's gap-filled total, with q10/q90
    "key_numbers.typical_season": "totals",
    "trend.annual.smooth": "trend",  # the smooth total, with its 95% band
    "key_numbers.trend": "trend",
    "trend.passage": "season",  # each year's modelled passage dates
    "trend.passage_q": "season",  # the smooth passage dates (lines on the phenology)
    "key_numbers.passage": "season",  # when its `source` is `smooth`
    "season.share": "totals",  # every day's share of the year, gap-filled (`season.source: gam`)
    "season.passage": "totals",  # every year's passage dates, over the gap-filled days
    "season.chances": "totals",  # the chances of a full day, over the gap-filled days
}
NO_CAVEAT = ("trend.passage_q",)  # elements shown only when their claim is `show`
COUNTS_FALLBACK = ("season.share", "season.passage", "season.chances")  # from the counts if hidden


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
        "season_birds": float(a["observed"].median()),
        "interval_ratio": _ratio(a["q90"], a["q10"]),
        "smooth_band": _ratio(a["smooth_q97.5"], a["smooth_q2.5"]),
        "passage_band": float((p["hi"] - p["lo"]).median()),
        "tails_band": (
            float(np.maximum(p["q10_hi"] - p["q10_lo"], p["q90_hi"] - p["q90_lo"]).median())
            if "q10_lo" in p
            else None
        ),
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
            "tails_band": _level(*PASSAGE_BAND, x["tails_band"]),
            "few_birds": "caveat" if x["season_birds"] < SEASON_BIRDS else None,
            "beyond_counting": "caveat" if x["beyond_counting"] else None,
        }
    )
    return {"totals": totals, "trend": trend, "season": season}


def element(field: str, cls: str) -> str:
    """The class of a drawn field from its claim's: a caveat hides a field in `NO_CAVEAT`; a field
    in `COUNTS_FALLBACK` whose claim is hidden is built from the counts, and shown."""
    if field in COUNTS_FALLBACK and cls == "hide":
        return "show"
    return "hide" if field in NO_CAVEAT and cls == "caveat" else cls


def fills_season(block: dict | None) -> bool:
    """Whether the season block is built from the trend's gap-filled days (`season.season_block`):

    with a trend whose totals are not hidden.
    """
    return block is not None and block["totals"]["class"] != "hide"


def reliability_block(trend: dict | None, bench, window: dict, profile: str) -> dict | None:
    """The `reliability` block of a taxon with a trend; None without one."""
    if trend is None:
        return None
    x = inputs(trend, bench, window, profile)
    a = pd.DataFrame(trend["annual"])
    return {
        "method": METHOD,
        **(c := classes(x)),
        "elements": {k: element(k, c[claim]["class"]) for k, claim in ELEMENTS.items()},
        "estimated_years": [int(y) for y in a.loc[a["observed_share"] < ESTIMATED_SHARE, "year"]],
        "inputs": x,
    }
