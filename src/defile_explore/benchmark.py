"""How well a taxon's trend model is known, tested on its own data: the `benchmark` block.

Two tests, on the taxon's model frame (its years and model window, `trend.model_frame`):

1. Gap transplant (gap filling, what Explore shows): the well-counted recent years (mean coverage
   of counted days >= `TARGET_MIN_COVERAGE`, from `TARGET_FROM`, when counts are timed so the birds
   of a hidden hour are known) are given older years' gaps, hour by hour and day by day
   (`DONOR_YEARS`), all at once, with the donors rotating over `ROTATIONS` refits (the full fit's
   smoothing parameters). Each hidden total is filled as the export fills it (posterior predictive
   median and intervals) and compared with the birds the year actually counted. The ratio
   estimator (birds kept / coverage kept, applied to the hidden coverage) is the reference.
2. Recent years (extrapolation of the trend): fitted on the seasons before `RECENT_FROM`
   (smoothing parameters refitted), the counted birds of the later seasons predicted from the model
   alone.

`score` sums each up as the median absolute log error (|log(estimate / counted)|, one bird added to
both so an estimate of 0 stays finite; 0.05 is 5%), its
median sign (bias) and how often the 80% and 95% intervals hold the truth. A taxon without target
or donor years gets no gap test (`taxon_benchmark`'s `gap` is None).

`scripts/benchmark_trend.py` runs the same tests on chosen taxa with the season-only model beside
the GAM, and draws them.
"""

import numpy as np
import pandas as pd
from scipy.special import logsumexp

from defile_explore import export as E
from defile_explore import profile as P
from defile_explore import trend as T

METHOD = "benchmark@1"
TARGET_MIN_COVERAGE = 0.8
TARGET_FROM = 2014  # timed counts, so the birds of hidden hours are known
DONOR_YEARS = (1993, 2013)  # their gaps (hours and days not counted) are transplanted
RECENT_FROM = 2023
ROTATIONS = 3  # refits per taxon in test 1, each hiding every target year
KINDS = ("hours", "day")  # hidden hours of a day still partly counted; whole days hidden
QUANTILE_COLUMNS = ["q2.5", "q10", "q90", "q97.5"]


def hourly_effort(effort: pd.DataFrame) -> pd.Series:
    """Hours counted per solar hour, by date, on counted days."""
    e = effort[effort["state"] == "counted"]
    return pd.Series(list(e["hourly"]), index=e["date"])


# --- test 1: gap transplant -------------------------------------------------


def transplant(target: pd.DataFrame, donor_year: int, ctx: dict) -> pd.DataFrame:
    """The target year's rows with the donor year's gaps: `c` and `y` become the coverage and birds
    kept, `c_target` the coverage the target counted.

    Hidden birds are known from their hour; untimed birds are kept in proportion to the coverage
    kept.
    """
    eff = hourly_effort(ctx["effort"])
    hourly = {d: h.set_index("hour")["count"] for d, h in ctx["hourly"].groupby("date")}
    zero = np.zeros(E.HOURS)
    rows = []
    for _, r in target.iterrows():
        tgt = eff.get(r["date"], zero) if r["c"] > 0 else zero
        donor_date = r["date"].replace(year=donor_year)
        kept = np.minimum(tgt, eff.get(donor_date, zero))
        doy = r["date"].dayofyear  # as `profile.coverage` indexes the profile
        p = ctx["profile"][int(np.clip(doy, *P.PROFILE_DOY)) - P.PROFILE_DOY[0]]
        c_kept = float((p * kept).sum())
        birds_h = np.zeros(E.HOURS)
        if r["y"] > 0 and r["date"] in hourly:
            h = hourly[r["date"]]
            birds_h[h.index.to_numpy()] = h.to_numpy()
        share = np.divide(kept, tgt, out=np.zeros(E.HOURS), where=tgt > 0)
        untimed = max(r["y"] - birds_h.sum(), 0.0)
        y_kept = (birds_h * share).sum() + (untimed * c_kept / r["c"] if r["c"] > 0 else 0.0)
        extra = r["extra"]
        c_target = r["c"]
        hidden = r["y"] - y_kept
        if c_kept < T.MIN_COVERAGE:  # as `model_frame` would: kept as counted, not modelled
            extra, c_target = extra + y_kept, c_target - c_kept
            c_kept = y_kept = 0.0
        rows.append(
            {
                "c": c_kept,
                "y": y_kept,
                "extra": extra,
                "c_target": c_target,
                "hidden": hidden,
                "kind": "hours" if c_kept > 0 else "day",
            }
        )
    out = target.drop(columns=["c", "y", "extra"]).reset_index(drop=True)
    return pd.concat([out, pd.DataFrame(rows)], axis=1)


def targets_and_donors(frame: pd.DataFrame, start: int) -> tuple[list[int], list[int]]:
    counted = frame[frame["c"] > 0].groupby("year")["c"].mean()
    targets = [
        int(y)
        for y in counted.index
        if y >= max(TARGET_FROM, start) and counted[y] >= TARGET_MIN_COVERAGE
    ]
    donors = [y for y in range(DONOR_YEARS[0], DONOR_YEARS[1] + 1) if y >= start]
    return targets, donors


def gap_trials(ctx: dict, frame: pd.DataFrame, fits: dict, kappa: float) -> pd.DataFrame:
    """`ROTATIONS` refits per fitted variant (`fits`: `{variant: trend.Fit}` on `frame`), each with
    every target year given a different donor's gaps at once (the donors rotate), so a taxon has
    `ROTATIONS` x targets trials per method; plus the ratio estimator.

    Empty without targets or donors.
    """
    years = np.arange(ctx["start"], ctx["last"] + 1)
    doy_range = (int(frame["doy"].min()), int(frame["doy"].max()))
    targets, donors = targets_and_donors(frame, ctx["start"])
    if not targets or not donors:
        return pd.DataFrame()
    rows = []
    for k in range(ROTATIONS):
        donor = {
            ty: donors[(i + k * len(donors) // ROTATIONS) % len(donors)]
            for i, ty in enumerate(targets)
        }
        masked = {ty: transplant(frame[frame["year"] == ty], donor[ty], ctx) for ty in targets}
        train = pd.concat(
            [frame[~frame["year"].isin(targets)], *[m[frame.columns] for m in masked.values()]],
            ignore_index=True,
        )
        refits = {v: T.fit(v, train, years, doy_range, hyper=f.hyper) for v, f in fits.items()}
        for i, (ty, m) in enumerate(masked.items()):
            target = frame[frame["year"] == ty]
            base = {
                "rotation": k,
                "target": ty,
                "donor": donor[ty],
                "truth": target["y"].sum() + target["extra"].sum(),
                "kept": m["y"].sum() + m["extra"].sum(),
                "hidden_coverage": 1 - m["c"].sum() / m["c_target"].sum(),
            }
            for kind in KINDS:
                base[f"truth_{kind}"] = m.loc[m["kind"] == kind, "hidden"].sum()
            ratio = m["y"].sum() / max(m["c"].sum(), 1e-9)
            filled = ratio * (m["c_target"] - m["c"])
            rows.append(
                {
                    **base,
                    "method": "ratio",
                    "estimate": base["kept"] + filled.sum(),
                    **{f"est_{kind}": filled[m["kind"] == kind].sum() for kind in KINDS},
                }
            )
            for v, f in refits.items():
                draws = T.fill_draws(f, m, m["c_target"].to_numpy(), seed=i, kappa=kappa)
                q = np.quantile(draws.sum(axis=1), [0.025, 0.1, 0.5, 0.9, 0.975])
                missed = draws - (m["y"] + m["extra"]).to_numpy()
                rows.append(
                    {
                        **base,
                        "method": v,
                        "estimate": q[2],
                        **dict(zip(QUANTILE_COLUMNS, q[[0, 1, 3, 4]])),
                        **{
                            f"est_{kind}": np.median(
                                missed[:, (m["kind"] == kind).to_numpy()].sum(axis=1)
                            )
                            for kind in KINDS
                        },
                    }
                )
    return pd.DataFrame(rows)


# --- test 2: recent years ---------------------------------------------------


def recent_trials(ctx: dict, frame: pd.DataFrame, kappa: float, variants=("gam",)) -> pd.DataFrame:
    """Each season from `RECENT_FROM`, predicted by `variants` fitted on the seasons before it:

    counted birds against the predicted total of the same days, and the days' mean log predictive
    density. Empty if no season from `RECENT_FROM` was counted or none before it.
    """
    years = np.arange(ctx["start"], ctx["last"] + 1)
    doy_range = (int(frame["doy"].min()), int(frame["doy"].max()))
    future = frame["year"] >= RECENT_FROM
    test = frame[future & (frame["c"] > 0)]
    if test.empty or ctx["start"] >= RECENT_FROM:
        return pd.DataFrame()
    train = frame.assign(c=frame["c"].where(~future, 0.0), y=frame["y"].where(~future, 0.0))
    rows = []
    for v in variants:
        f = T.fit(v, train, years, doy_range)
        unseen = test.assign(c=0.0, y=0.0)
        draws = T.fill_draws(f, unseen, test["c"].to_numpy(), seed=1, kappa=kappa)
        beta = f.draws(T.DRAWS, np.random.default_rng(2))
        mu = np.exp(T.eta_draws(f, test["year"], test["doy"], beta)) * test["c"].to_numpy()
        ll = logsumexp(T.nb_loglik(test["y"].to_numpy(), mu, f.theta), axis=0) - np.log(len(mu))
        for yr in sorted(test["year"].unique()):
            sel = (test["year"] == yr).to_numpy()
            q = np.quantile(draws[:, sel].sum(axis=1), [0.025, 0.1, 0.5, 0.9, 0.975])
            rows.append(
                {
                    "method": v,
                    "year": int(yr),
                    "truth": test["y"].to_numpy()[sel].sum(),
                    "estimate": q[2],
                    **dict(zip(QUANTILE_COLUMNS, q[[0, 1, 3, 4]])),
                    "log_score": ll[sel].mean(),
                    "days": int(sel.sum()),
                }
            )
    return pd.DataFrame(rows)


# --- scores and the block ---------------------------------------------------


def score(df: pd.DataFrame, by=("method",)) -> pd.DataFrame:
    """Per `by`: median |log(estimate / truth)| (`abs_log_err`), its median sign (`bias`), how
    often the 80% and 95% intervals hold the truth (`cover80`, `cover95`; NaN without intervals)
    and the trials (`n`).

    Trials with no bird counted are left out.
    """
    d = df[df["truth"] > 0].assign(err=lambda x: np.log((x["estimate"] + 1) / (x["truth"] + 1)))
    if "q10" not in d:
        d = d.assign(**{c: np.nan for c in QUANTILE_COLUMNS})
    d["in80"] = ((d["truth"] >= d["q10"]) & (d["truth"] <= d["q90"])).where(d["q10"].notna())
    d["in95"] = ((d["truth"] >= d["q2.5"]) & (d["truth"] <= d["q97.5"])).where(d["q10"].notna())
    g = d.groupby(list(by))
    return pd.DataFrame(
        {
            "abs_log_err": g["err"].apply(lambda x: x.abs().median()),
            "bias": g["err"].median(),
            "cover80": g["in80"].mean(),
            "cover95": g["in95"].mean(),
            "n": g.size(),
        }
    )


def _summary(df: pd.DataFrame, method: str) -> dict | None:
    if df.empty:
        return None
    s = score(df)
    if method not in s.index:
        return None
    return {k: float(v) for k, v in s.loc[method].items()}


def taxon_benchmark(
    days: pd.DataFrame,
    hourly: pd.DataFrame,
    effort: pd.DataFrame,
    profile: np.ndarray,
    first_year: int,
    last_year: int,
    window: tuple[int, int] | None = None,
) -> dict:
    """Both tests on the GAM as exported (same years, window, profile): the `benchmark` block.

    `gap`: test 1's scores of the GAM and, as `gap_ratio`, of the ratio estimator, with the median
    share of the target years' coverage hidden; `gap_trials`: one row per trial (target, donor,
    counted, estimate and intervals). `recent`, `recent_trials`: test 2 the same way, with the mean
    log predictive density per day. Each is None (or empty) where the test cannot run.
    """
    ctx = {
        "start": first_year,
        "last": last_year,
        "effort": effort,
        "hourly": hourly,
        "profile": profile,
    }
    frame = T.model_frame(days, effort, P.coverage(effort, profile), first_year, last_year, window)
    years = np.arange(first_year, last_year + 1)
    doy_range = (int(frame["doy"].min()), int(frame["doy"].max()))
    kappa = T.hour_dispersion(days, hourly, effort, profile)
    fit = T.fit("gam", frame, years, doy_range)
    gaps = gap_trials(ctx, frame, {"gam": fit}, kappa)
    recent = recent_trials(ctx, frame, kappa)
    gap = _summary(gaps, "gam")
    if gap is not None:
        gap["hidden_coverage"] = float(gaps["hidden_coverage"].median())
    recent_score = _summary(recent, "gam")
    if recent_score is not None:
        recent_score["log_score"] = float(recent["log_score"].mean())
    trial_cols = ["target", "donor", "truth", "estimate", *QUANTILE_COLUMNS]
    return {
        "method": METHOD,
        "gap": gap,
        "gap_ratio": _summary(gaps, "ratio"),
        "gap_trials": gaps[gaps["method"] == "gam"][trial_cols] if gap else None,
        "recent": recent_score,
        "recent_trials": (
            recent[["year", "truth", "estimate", *QUANTILE_COLUMNS]] if recent_score else None
        ),
    }
