"""Age and sex of the birds counted, where the counters recorded them.

Only some birds are aged or sexed, and the share changes from year to year, so every value is a
share among the birds given an age (or sex) in that year, with a 95% Wilson interval. That assumes
the birds aged are representative of those passing. Counters often tag only one class and leave the
other blank (Red Kite since 2022: thousands of juveniles a year, almost no adults), which makes the
share among aged birds meaningless. So a year counts as usable only if enough birds were aged
(`MIN_BIRDS`, and `MIN_SHARE` of those counted) and both classes were recorded (the smaller at
least `MIN_CLASS_SHARE` of the aged birds), and a taxon gets the panel only with `MIN_YEARS` usable
years (`usable_years`; `settings` can override). A taxon where one class is genuinely rare loses
years to that rule: add them back as an override, with the reason.

Age codes: `A` adult; `1` (first calendar year), `J` (juvenile), `I` (immature) and `2` (second
calendar year) are one class, non-adult, because the codes used switch between years (Red Kite:
mostly `1` in 2019, mostly `I` in 2022), which says more about the recording than the birds. Sex:
`M`, `F`, and `FC` (female-coloured: a female or a juvenile, not separable in the field); `F` and
`FC` are one class, female-type, against males.

Age and sex make the same block (`share_block`): the share of the first class (non-adult, male) per
year, the birds classed and counted behind it, and when each class passes.
"""

import numpy as np
import pandas as pd

METHOD = "demography@3"
NON_ADULT = ("1", "J", "I", "2")
ADULT = ("A",)
SEXES = ("M", "F", "FC")
CLASSES = {"age": (NON_ADULT, ADULT), "sex": (("M",), ("F", "FC"))}  # the two classes compared
NAMES = {"age": ("non_adult", "adult"), "sex": ("male", "female_type")}  # ... in the block
CODES = {"age": NON_ADULT + ADULT, "sex": SEXES}  # codes counted per year in the block
MIN_BIRDS = 20  # aged (or sexed) birds in a year for the year to count
MIN_SHARE = 0.05  # ... and this share of the year's birds counted
MIN_CLASS_SHARE = 0.1  # ... and the smaller class this share of them (both classes recorded)
MIN_YEARS = 3  # usable years for a taxon to get the panel
Z95 = 1.96


def wilson(k, n, z: float = Z95):
    """Wilson score interval of a binomial share k / n, `(lo, hi)`; (nan, nan) where n = 0."""
    k, n = np.asarray(k, float), np.asarray(n, float)
    with np.errstate(invalid="ignore", divide="ignore"):
        p = k / n
        d = 1 + z**2 / n
        c = (p + z**2 / (2 * n)) / d
        h = z * np.sqrt(p * (1 - p) / n + z**2 / (4 * n**2)) / d
    return c - h, c + h


def per_year(rows: pd.DataFrame, field: str, counted: pd.Series) -> pd.DataFrame:
    """Birds per year and code of `field` (`age` or `sex`), with `counted` (all birds that
    year)."""
    r = rows[rows[field].notna()]
    g = r.pivot_table(
        index=r["date"].dt.year, columns=field, values="count", aggfunc="sum"
    ).fillna(0)
    g.index.name = "year"
    g.columns = [str(c) for c in g.columns]
    return g.join(counted.rename("counted"), how="left").fillna({"counted": 0})


def usable_years(rows: pd.DataFrame, field: str, counted: pd.Series) -> list[int]:
    """Years with at least `MIN_BIRDS` birds given a `field`, `MIN_SHARE` of those counted, and
    each of the field's two `CLASSES` at least `MIN_CLASS_SHARE` of them; empty if fewer than
    `MIN_YEARS` such years."""
    g = per_year(rows, field, counted)
    a, b = (g[[c for c in cls if c in g]].sum(axis=1) for cls in CLASSES[field])
    n = a + b
    ok = (n >= MIN_BIRDS) & (n >= MIN_SHARE * g["counted"])
    ok &= np.minimum(a, b) >= MIN_CLASS_SHARE * n
    years = [int(y) for y in g.index[ok]]
    return years if len(years) >= MIN_YEARS else []


def share_block(rows: pd.DataFrame, counted: pd.Series, years: list[int], field: str):
    """The share of the field's first class (`CLASSES`, named by `NAMES`) per usable year and over
    them all, with the birds classed (`n`, either class) and counted (`counted`) behind each, and
    when each class passes (cumulative share by day of year); None without usable years."""
    if not years:
        return None
    first, second = CLASSES[field]
    g = per_year(rows, field, counted).loc[years]
    k = g[[c for c in first if c in g]].sum(axis=1)
    n = k + g[[c for c in second if c in g]].sum(axis=1)
    lo, hi = wilson(k, n)
    K, N = float(k.sum()), float(n.sum())
    LO, HI = wilson(K, N)
    names = NAMES[field]
    r = rows[rows["date"].dt.year.isin(years) & rows[field].isin(first + second)]
    r = r.assign(cls=np.where(r[field].isin(first), names[0], names[1]))
    by_day = r.pivot_table(
        index=r["date"].dt.dayofyear, columns="cls", values="count", aggfunc="sum"
    ).fillna(0)
    timing = {"doy": by_day.index.to_numpy()}
    for cls in names:
        if cls in by_day:
            timing[cls] = (by_day[cls].cumsum() / by_day[cls].sum()).to_numpy()
            timing[f"{cls}_birds"] = float(by_day[cls].sum())
    return {
        "method": METHOD,
        "classes": list(names),
        "years": years,
        "n": n.to_numpy(),
        "counted": g["counted"].to_numpy(),
        "k": k.to_numpy(),
        "share": (k / n).to_numpy(),
        "lo": lo,
        "hi": hi,
        "codes": {c: g[c].to_numpy() if c in g else np.zeros(len(g)) for c in CODES[field]},
        "overall": {
            "share": K / N,
            "lo": float(LO),
            "hi": float(HI),
            "n": N,
            "counted": float(g["counted"].sum()),
        },
        "timing": timing,
    }


def age_block(rows: pd.DataFrame, counted: pd.Series, years: list[int]) -> dict | None:
    """The non-adult share (`share_block`)."""
    return share_block(rows, counted, years, "age")


def sex_block(rows: pd.DataFrame, counted: pd.Series, years: list[int]) -> dict | None:
    """The male share (`share_block`)."""
    return share_block(rows, counted, years, "sex")
