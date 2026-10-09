"""The Explore build: one shared stage, then one function per taxon that makes all its blocks.

`load_shared` reads the release and computes what every taxon needs (effort, daily and hourly
counts, taxa, the time-of-day profiles, age and sex rows), cached on disk under `data/cache/` keyed
by the release and the code that produced it. `taxon_jobs` cuts it into one job per taxon, and
`build_taxon` turns a job into the taxon's species
file:
the raw tables as before (`days`,
`hourly`, `annual`, `profile`), `trend` (the GAM, cached per taxon), and the blocks
derived from the same data and fit: `season`, `daytime`, `age`, `sex`, `records`,
`key_numbers`, the written `accounts`, with the taxon's `settings`, `window`, `links` and
`diagnostics`. The trend is fitted
on the taxon's model window, the season panels show its view window (`defile_explore.window`).

Each block names its `method` (`<name>@<version>`), matched by a section of defileViz's method
page; bump the version when what a block means changes.
"""

import hashlib
import inspect
import os
import pickle
from dataclasses import dataclass

import numpy as np
import pandas as pd

from defile_explore import accounts as A
from defile_explore import daytime as Y
from defile_explore import demography as G
from defile_explore import export as E
from defile_explore import profile as P
from defile_explore import release as R
from defile_explore import remarks as M
from defile_explore import season as S
from defile_explore import settings as X
from defile_explore import timeofday
from defile_explore import trend as T
from defile_explore import window as W

CACHE_DIR = os.path.join("cache")  # under the data dir
RECORD_DAYS = 15  # top days listed
KEY_YEARS = 10  # "typical season": the last this many complete seasons
BEST_HOURS_SHARE = 0.6  # best hours: the fewest hours holding this share of a peak day
DEMOGRAPHY_COLUMNS = ["taxon_id", "date", "count", "age", "sex"]


def _source_hash(*modules) -> str:
    h = hashlib.sha1()
    for m in modules:
        h.update(inspect.getsource(m).encode())
    return h.hexdigest()[:10]


def release_key(data_dir: str) -> str:
    """The release (file names, sizes, modification times) and the code of the shared stage."""
    folder = os.path.join(data_dir, R.DATASET_DIR)
    h = hashlib.sha1()
    for f in sorted(os.listdir(folder)):
        st = os.stat(os.path.join(folder, f))
        h.update(f"{f}:{st.st_size}:{st.st_mtime_ns}".encode())
    h.update(_source_hash(R, E, P, M, timeofday).encode())
    return h.hexdigest()[:12]


@dataclass
class Shared:
    key: str
    metadata: dict
    taxonomy: pd.DataFrame
    taxa: pd.DataFrame
    effort: pd.DataFrame
    days: pd.DataFrame
    hourly: pd.DataFrame
    profiles: dict
    source: pd.Series
    demography: pd.DataFrame
    remarks: pd.DataFrame
    last_year: int


def load_shared(data_dir: str, ebird: pd.DataFrame, use_cache: bool = True) -> Shared:
    """The shared stage, from `<data_dir>/cache/shared-<key>.pkl` when there."""
    key = release_key(data_dir)
    path = os.path.join(data_dir, CACHE_DIR, f"shared-{key}.pkl")
    if use_cache and os.path.exists(path):
        with open(path, "rb") as f:
            return pickle.load(f)
    surveys, counts, taxonomy, metadata = E.read_release(data_dir)
    effort = E.build_effort(surveys)
    days = E.daily_counts(counts)
    hourly = E.hourly_counts(counts)
    taxonomy, days, hourly = E.add_combined(taxonomy, days, hourly)
    taxa = E.build_taxa(taxonomy, days, ebird)
    profiles, source = P.build_profiles(taxa, days, hourly, effort)
    taxa["profile"] = taxa["taxon_id"].map(source)
    demo = counts.loc[
        (counts["count_category"] == R.MAIN_CATEGORY)
        & (counts["age"].notna() | counts["sex"].notna()),
        DEMOGRAPHY_COLUMNS,
    ]
    shared = Shared(
        key=key,
        metadata=metadata,
        taxonomy=taxonomy,
        taxa=taxa,
        effort=effort,
        days=days,
        hourly=hourly,
        profiles=profiles,
        source=source,
        demography=demo.astype({"age": "string", "sex": "string"}),
        remarks=M.day_remarks(surveys, counts),
        last_year=int(days["date"].max().year) - len(E.partial_years(days)),
    )
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "wb") as f:
        pickle.dump(shared, f)
    return shared


def taxon_jobs(
    shared: Shared, taxon_ids, overrides: dict, accounts: dict, data_dir: str, use_cache=True
):
    """One job per taxon: everything `build_taxon` needs, cut from the shared stage and the
    `accounts` (`accounts.load_accounts`)."""
    tax = shared.taxonomy.set_index("taxon_id")
    for taxon_id in taxon_ids:
        t = shared.taxa.set_index("taxon_id").loc[taxon_id]
        members = t["members"] if isinstance(t["members"], list) else [taxon_id]
        src = shared.source[taxon_id]
        yield {
            "taxon": {
                "taxon_id": taxon_id,
                **t.to_dict(),
                "ebird_code": tax["ebird_code"].get(taxon_id),
                "trektellen_species_id": tax["trektellen_species_id"].get(taxon_id),
            },
            "days": shared.days[shared.days["taxon_id"] == taxon_id].drop(columns="taxon_id"),
            "hourly": shared.hourly[shared.hourly["taxon_id"] == taxon_id].drop(
                columns="taxon_id"
            ),
            "effort": shared.effort,
            "profile": shared.profiles[taxon_id if src == "own" else src],
            "profile_source": src,
            "accounts": A.taxon_rows(accounts, taxon_id),
            "demography": shared.demography[shared.demography["taxon_id"].isin(members)],
            "remarks": shared.remarks[
                shared.remarks["taxon_id"].isna() | shared.remarks["taxon_id"].isin(members)
            ],
            "overrides": overrides,
            "last_year": shared.last_year,
            "cache": os.path.join(data_dir, CACHE_DIR, "trend") if use_cache else None,
            "shared_key": shared.key,
        }


# --- trend ------------------------------------------------------------------------------------


def trend_key(job: dict, start: int, window: tuple[int, int]) -> str:
    h = hashlib.sha1()
    for part in (job["days"], job["hourly"]):
        h.update(pd.util.hash_pandas_object(part, index=False).to_numpy().tobytes())
    h.update(np.ascontiguousarray(job["profile"]).tobytes())
    h.update(f"{start}:{window}:{job['last_year']}:{job['shared_key']}".encode())
    h.update(_source_hash(T).encode())
    return h.hexdigest()[:16]


def trend_block(job: dict, start: int, window: tuple[int, int]) -> dict:
    """`taxon_trend`, from the cache when the taxon's data, profile, years, window and the trend
    code are unchanged."""
    path = None
    if job["cache"]:
        path = os.path.join(
            job["cache"], f"{job['taxon']['taxon_id']}-{trend_key(job, start, window)}.pkl"
        )
        if os.path.exists(path):
            with open(path, "rb") as f:
                return pickle.load(f)
    t = T.taxon_trend(
        days=job["days"],
        hourly=job["hourly"],
        effort=job["effort"],
        profile=job["profile"],
        first_year=start,
        last_year=job["last_year"],
        window=window,
    )
    out = {
        "method": "trend-gam@1",
        "model": t["model"],
        "first_year": t["first_year"],
        "last_year": t["last_year"],
        "window": W.as_dates(window),
        "theta": t["theta"],
        "kappa": t["kappa"],
        "annual": E.records(t["annual"]),
        "passage": E.records(t["passage"]),
        "passage_q": E.records(t["passage_q"]),
        "season": t["season"],
        "episodes": t["episodes"],
    }
    if path:
        os.makedirs(job["cache"], exist_ok=True)
        with open(path, "wb") as f:
            pickle.dump(out, f)
    return out


# --- derived blocks ---------------------------------------------------------------------------


def records_block(days: pd.DataFrame, remarks: pd.DataFrame) -> dict:
    """The `RECORD_DAYS` days with the most birds, all years, with what was written about each
    (`remarks.notes_of`)."""
    top = days.dropna(subset=["count"]).nlargest(RECORD_DAYS, "count")[["date", "count"]]
    top["notes"] = [M.notes_of(remarks, d) for d in top["date"]]
    return {"method": "records@2", "top_days": top}


def best_hours(profile: np.ndarray, doy: float) -> dict:
    """The fewest local clock hours holding `BEST_HOURS_SHARE` of the profile on day `doy`, as a
    range from the first to the last of them."""
    p = profile[int(np.clip(round(doy), *P.PROFILE_DOY)) - P.PROFILE_DOY[0]]
    order = np.argsort(p)[::-1]
    top = order[: np.searchsorted(np.cumsum(p[order]), BEST_HOURS_SHARE) + 1]
    return {"from": int(top.min()), "to": int(top.max()) + 1, "share": float(p[top].sum())}


def main_passage(trend: dict | None, season: dict, last_year: int) -> dict:
    """The 10/50/90% passage dates: of the smooth season in the last year if there is a trend, else
    of the counted seasons of the last `KEY_YEARS` pooled."""
    if trend:
        q = trend["passage_q"][-1]
        return {"q10": q["q10"], "q50": q["q50"], "q90": q["q90"], "source": "smooth"}
    years = np.asarray(season["years"])
    share = np.asarray(season["share"], float)[years > last_year - KEY_YEARS]
    pooled = np.nan_to_num(share).sum(axis=0)
    q = S.cumulative_quantiles(np.asarray(season["doy"]), pooled)
    return {"q10": q[0], "q50": q[1], "q90": q[2], "source": "counted"}


def key_numbers(job, frame, trend, season, records) -> dict:
    last = job["last_year"]
    passage = main_passage(trend, season, last)
    out = {"passage": passage}
    if not np.isnan(passage["q50"]):
        out["best_hours"] = best_hours(job["profile"], passage["q50"])
        f = frame[
            (frame["year"] > last - KEY_YEARS)
            & frame["doy"].between(passage["q10"], passage["q90"])
            & (frame["c"] >= P.COVERAGE_MIN)
        ]
        if len(f):
            out["chance"] = {
                "at_least_1": float((f["y"] >= 1).mean()),
                "at_least_10": float((f["y"] >= 10).mean()),
                "days": len(f),
            }
    if trend:
        a = pd.DataFrame(trend["annual"])
        recent = a[a["year"] > last - KEY_YEARS]["total"]
        out["typical_season"] = {
            "median": float(recent.median()),
            "min": float(recent.min()),
            "max": float(recent.max()),
            "years": [last - KEY_YEARS + 1, last],
        }
        out["trend"] = {
            "change": float(a["smooth"].iloc[-1] / a["smooth"].iloc[0] - 1),
            "from": int(a["year"].iloc[0]),
            "to": int(a["year"].iloc[-1]),
        }
    top = records["top_days"]
    if len(top):
        out["record"] = {"date": top["date"].iloc[0], "count": float(top["count"].iloc[0])}
    return out


def daytime_parts(trend: dict | None, season: dict, last_year: int) -> list[float] | None:
    """Days of year splitting the passage into early / peak / late (`Y.PART_QUANTILES`): of the
    smooth season of the last year, else of the pooled counted seasons."""
    if trend:
        doy, mu = np.asarray(trend["season"]["doy"]), np.asarray(trend["season"][str(last_year)])
    else:
        doy = np.asarray(season["doy"])
        mu = np.nan_to_num(np.asarray(season["share"], float)).sum(axis=0)
    q = S.cumulative_quantiles(doy, mu, Y.PART_QUANTILES)
    return None if np.isnan(q).any() else q


def window_block(rule: dict, settings: dict) -> dict:
    """The taxon's windows as dates, the envelopes they come from, and where the passage reaches
    the edge of counting."""
    return {
        "method": rule["method"],
        "model": W.as_dates(settings["model_window"].value),
        "view": W.as_dates(settings["view_window"].value),
        "default": W.as_dates(W.default_window()),
        "model_envelope": W.as_dates(rule["model_envelope"]),
        "view_envelope": W.as_dates(rule["view_envelope"]),
        "beyond_counting": rule["beyond_counting"],
    }


def diagnostics(job, frame, trend, daytime, age, window) -> dict:
    """Numbers and flags for the QA review: what may make this taxon's page misleading."""
    last = job["last_year"]
    recent = frame[(frame["year"] > last - KEY_YEARS) & (frame["c"] > 0)]
    d = {
        "profile": job["profile_source"],
        "coverage_recent": float(recent["c"].mean()) if len(recent) else None,
        "timed_days": daytime["days"] if daytime else 0,
        "daytime_shift": daytime["change"].get("shift") if daytime else None,
        "age_years": len(age["years"]) if age else 0,
    }
    flags = []
    if job["profile_source"] != "own":
        flags.append("borrowed_profile")
    if window["beyond_counting"] and job["taxon"]["tier"] == "full":
        flags.append("passage_beyond_counting")
    if trend:
        a = pd.DataFrame(trend["annual"])
        d["interval_ratio"] = float((a["q90"] / a["q10"].clip(lower=1)).median())
        d["theta"], d["kappa"] = trend["theta"], trend["kappa"]
        if d["interval_ratio"] > 2:
            flags.append("wide_intervals")
    if daytime is None:
        flags.append("no_timed_days")
    d["flags"] = flags
    return d


# --- one taxon --------------------------------------------------------------------------------


def build_taxon(job: dict) -> tuple[str, dict]:
    """A taxon's species file (see the module docstring)."""
    taxon = job["taxon"]
    taxon_id = taxon["taxon_id"]
    effort, profile, last = job["effort"], job["profile"], job["last_year"]
    days = job["days"]
    counted = days.groupby(days["date"].dt.year)["count"].sum()
    settings, links = X.resolve(taxon, job["demography"], counted, job["overrides"], last)
    start = settings["start_year"].value

    c = P.coverage(effort, profile)
    d = P.adjust_days(days, effort, c, start)
    annual = E.annual(d).merge(
        P.annual_index(d, effort, c, start)[["year", "c_mean", "index"]], on="year", how="left"
    )
    season_frame = T.model_frame(days, effort, c, start, last, P.PROFILE_DOY)
    rule = W.windows(season_frame)
    settings |= X.resolve_windows(taxon_id, rule, job["overrides"])
    model_w, view_w = settings["model_window"].value, settings["view_window"].value
    window = window_block(rule, settings)
    trend = trend_block(job, start, model_w) if settings["trend"].value else None
    frame = season_frame[season_frame["doy"].between(*view_w)]
    season = S.season_block(frame, last, model_w)
    daytime = Y.daytime_block(
        days, job["hourly"], effort, profile, daytime_parts(trend, season, last)
    )
    age = G.age_block(job["demography"], counted, settings["age_years"].value)
    sex = G.sex_block(job["demography"], counted, settings["sex_years"].value)
    records = records_block(d, job["remarks"])
    return taxon_id, {
        "taxon_id": taxon_id,
        "settings": {k: v.as_dict() for k, v in settings.items()},
        "window": window,
        "links": links,
        "key_numbers": key_numbers(job, frame, trend, season, records),
        "days": E.columns(d),
        "hourly": E.columns(job["hourly"]),
        "annual": E.records(annual),
        "profile": {"source": job["profile_source"], **P.profile_table(profile)},
        "accounts": A.accounts_block(job["accounts"], last),
        "trend": trend,
        "season": season,
        "daytime": daytime,
        "age": age,
        "sex": sex,
        "records": records,
        "diagnostics": diagnostics(job, frame, trend, daytime, age, window),
    }
