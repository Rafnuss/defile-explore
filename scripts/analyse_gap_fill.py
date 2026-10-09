#!/usr/bin/env python3
"""Day by day gap filling: the trend GAM against the season panel's linear interpolation.

The gap transplant of `defile_explore.benchmark` (well-counted recent years given older years'
gaps, hour by hour and day by day, `ROTATIONS` refits with the donors rotating), scored on each
day instead of the season: does a complete series from the GAM beat the interpolation the season
panel fills its gaps with (`season.daily_rates`)? Each target day is filled over the coverage the
target counted (`c_target`), so the truth is the birds it actually counted:

- `gam`: the posterior predictive (`trend.fill_draws`), the mean per day as `trend.days` exports
  it; per-draw passage dates and totals give their median and 80% interval;
- `interp`: the season panel as built: birds / c on days kept at `COVERAGE_MIN` or more, linear
  interpolation over the other days of the year (edges held), times `c_target`: a day kept below
  `COVERAGE_MIN` loses its count;
- `interp_kept`: the same rate, but a partly kept day keeps its birds and only its hidden coverage
  is filled at the interpolated rate (a fairer baseline).

Scores, over the target-year trials:

- `day_err`: median |log((estimate + 1) / (truth + 1))| over the days with any coverage hidden
  and birds counted, split into `hours` (partly kept) and `day` (wholly hidden); `day_werr` the
  same as a mean weighted by the day's birds; `cover80` (GAM): truth within the day's 80% interval;
  `q50_cover80` (GAM): the true median date within the draws' 80% interval (`_1d`: give or take a
  day, as dates are whole days);
- `q10`, `q50`, `q90`: |error| in days of the year's 10/50/90% passage dates, from the cumulative
  series of the target's counted time (estimate against truth); `q50_cover80` (GAM);
- `total_err`: |log| error of the year's total (the existing benchmark's measure).

Writes `logs/qa/explore/gap_fill_days.csv` and `gap_fill_years.csv` and prints a summary.

Usage:
    python scripts/analyse_gap_fill.py [--taxa "Red Kite" ...] [--workers 8]
"""

import argparse
import json
import os
import time

import numpy as np
import pandas as pd

from defile_explore import accounts as A
from defile_explore import benchmark as B
from defile_explore import export as E
from defile_explore import pipeline as L
from defile_explore import profile as P
from defile_explore import settings as X
from defile_explore import trend as T
from defile_explore import window as W

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
METHODS = ("gam", "interp", "interp_kept")
QS = (0.1, 0.5, 0.9)
INTERVAL = (0.1, 0.9)


def passage_dates(doy: np.ndarray, series: np.ndarray) -> np.ndarray:
    """First day of year by which each share in `QS` of the series' total has passed; `series` is
    (n_days,) or (n_draws, n_days); NaN for an empty series."""
    s = np.atleast_2d(series)
    cum = np.cumsum(s, axis=1)
    tot = cum[:, -1:]
    with np.errstate(invalid="ignore", divide="ignore"):
        share = cum / tot
    out = np.stack([doy[np.argmax(share >= q, axis=1)] for q in QS], axis=1).astype(float)
    out[tot[:, 0] <= 0] = np.nan
    return out


def interp_rate(m: pd.DataFrame) -> np.ndarray:
    """Birds per full day of each row as the season panel has it: y / c where c >= COVERAGE_MIN,
    linear interpolation over the year's days between (edges held), 0 if none."""
    rate = (m["y"] / m["c"]).where(m["c"] >= P.COVERAGE_MIN)
    s = pd.Series(rate.to_numpy(), index=m["doy"].to_numpy()).sort_index()
    return s.interpolate(limit_direction="both").fillna(0).reindex(m["doy"]).to_numpy()


def run_taxon(args) -> dict:
    name, job, start, window = args
    t0 = time.time()
    ctx = {
        "start": start,
        "last": job["last_year"],
        "effort": job["effort"],
        "hourly": job["hourly"],
        "profile": job["profile"],
    }
    c = P.coverage(job["effort"], job["profile"])
    frame = T.model_frame(job["days"], job["effort"], c, start, job["last_year"], window)
    years = np.arange(start, job["last_year"] + 1)
    doy_range = (int(frame["doy"].min()), int(frame["doy"].max()))
    targets, donors = B.targets_and_donors(frame, start)
    if not targets or not donors:
        print(f"  {name}: no target or donor years", flush=True)
        return {"days": pd.DataFrame(), "years": pd.DataFrame()}
    kappa = T.hour_dispersion(job["days"], job["hourly"], job["effort"], job["profile"])
    full = T.fit("gam", frame, years, doy_range)
    day_rows, year_rows = [], []
    for k in range(B.ROTATIONS):
        donor = {
            ty: donors[(i + k * len(donors) // B.ROTATIONS) % len(donors)]
            for i, ty in enumerate(targets)
        }
        masked = {ty: B.transplant(frame[frame["year"] == ty], donor[ty], ctx) for ty in targets}
        train = pd.concat(
            [frame[~frame["year"].isin(targets)], *[m[frame.columns] for m in masked.values()]],
            ignore_index=True,
        )
        f = T.fit("gam", train, years, doy_range, hyper=full.hyper)
        for i, (ty, m) in enumerate(masked.items()):
            orig = frame[frame["year"] == ty].reset_index(drop=True)
            truth = (orig["y"] + orig["extra"]).to_numpy(float)
            kept = (m["y"] + m["extra"]).to_numpy(float)
            hidden = (m["c_target"] - m["c"]).to_numpy() > 1e-9
            draws = T.fill_draws(f, m, m["c_target"].to_numpy(), seed=k * 100 + i, kappa=kappa)
            rate = interp_rate(m)
            gap = np.clip(m["c_target"] - m["c"], 0, None).to_numpy()
            est = {
                "gam": draws.mean(axis=0),
                "interp": np.where(
                    m["c"] >= P.COVERAGE_MIN,
                    kept + rate * gap,
                    m["extra"].to_numpy() + rate * m["c_target"].to_numpy(),
                ),
                "interp_kept": kept + rate * gap,
            }
            lo, hi = np.quantile(draws, INTERVAL, axis=0)
            doy = m["doy"].to_numpy()
            true_q = passage_dates(doy, truth)[0]
            gam_q = passage_dates(doy, draws)
            base = {"taxon": name, "rotation": k, "target": ty, "donor": donor[ty]}
            for meth, e in est.items():
                for j in np.flatnonzero(hidden):
                    day_rows.append(
                        {
                            **base,
                            "method": meth,
                            "doy": int(doy[j]),
                            "kind": m["kind"].iat[j],
                            "truth": truth[j],
                            "estimate": e[j],
                            "in80": (lo[j] <= truth[j] <= hi[j]) if meth == "gam" else np.nan,
                        }
                    )
                q = np.nanmedian(gam_q, axis=0) if meth == "gam" else passage_dates(doy, e)[0]
                row = {
                    **base,
                    "method": meth,
                    "hidden_coverage": 1 - m["c"].sum() / m["c_target"].sum(),
                    "truth_total": truth.sum(),
                    "estimate_total": (np.median(draws.sum(axis=1)) if meth == "gam" else e.sum()),
                }
                for n, qn in enumerate(("q10", "q50", "q90")):
                    row[qn] = abs(q[n] - true_q[n])
                if meth == "gam":
                    qlo, qhi = np.nanquantile(gam_q[:, 1], INTERVAL)
                    row["q50_in80"] = qlo <= true_q[1] <= qhi
                    row["q50_in80_1d"] = qlo - 1 <= true_q[1] <= qhi + 1
                year_rows.append(row)
    print(f"  {name}: {len(targets)} target years, {time.time() - t0:.0f} s", flush=True)
    return {"days": pd.DataFrame(day_rows), "years": pd.DataFrame(year_rows)}


def summary(days: pd.DataFrame, years: pd.DataFrame) -> None:
    d = days[days["truth"] > 0].assign(
        err=lambda x: np.abs(np.log((x["estimate"] + 1) / (x["truth"] + 1)))
    )

    def day_scores(g):
        return pd.Series(
            {
                "day_err": g["err"].median(),
                "day_werr": np.average(g["err"], weights=g["truth"]),
                "cover80": g["in80"].astype(float).mean(),
                "n": len(g),
            }
        )

    print("\nDAYS with coverage hidden and birds counted (|log| error; lower is better)")
    print(
        d.groupby(["kind", "method"]).apply(day_scores, include_groups=False).round(3).to_string()
    )
    print("\nall days")
    print(d.groupby("method").apply(day_scores, include_groups=False).round(3).to_string())
    # per taxon, so a few taxa with many days do not decide it
    per_taxon = d.groupby(["taxon", "method"])["err"].median().unstack()
    for meth in ("interp", "interp_kept"):
        print(
            f"  gam better than {meth} in {(per_taxon['gam'] < per_taxon[meth]).sum()} / "
            f"{per_taxon[[meth, 'gam']].dropna().shape[0]} taxa (median day error)"
        )

    y = years[years["truth_total"] > 0].assign(
        total_err=lambda x: np.abs(np.log((x["estimate_total"] + 1) / (x["truth_total"] + 1)))
    )
    print("\nYEARS: passage dates (|error| in days) and total (|log| error), medians")
    agg = y.groupby("method")[["q10", "q50", "q90", "total_err"]].median()
    agg["q50_mean"] = y.groupby("method")["q50"].mean()
    for col in ("q50_in80", "q50_in80_1d"):
        agg[col.replace("in80", "cover80")] = y[col].astype(float).groupby(y["method"]).mean()
    agg["n"] = y.groupby("method").size()
    print(agg.round(3).to_string())
    pt = y.groupby(["taxon", "method"])["q50"].mean().unstack()
    for meth in ("interp", "interp_kept"):
        print(
            f"  gam's median date closer than {meth} in {(pt['gam'] < pt[meth]).sum()}, "
            f"farther in {(pt['gam'] > pt[meth]).sum()} / {len(pt)} taxa"
        )
    print(
        "  median hidden coverage of the target years: "
        f"{y[y['method'] == 'gam']['hidden_coverage'].median():.2f}"
    )


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--data-dir", default=os.path.join(ROOT, "data"))
    ap.add_argument("--export", default=os.path.join(ROOT, "data", "explore"))
    ap.add_argument("--out", default=os.path.join(ROOT, "logs", "qa", "explore"))
    ap.add_argument("--taxa", nargs="+", help="English names (default: every taxon with a trend)")
    ap.add_argument("--workers", type=int, default=8)
    args = ap.parse_args(argv)

    ebird = pd.read_csv(os.path.join(args.data_dir, E.EBIRD_TAXONOMY_FILE))
    shared = L.load_shared(args.data_dir, ebird)
    ids, names = [], []
    for tid, name in zip(shared.taxa["taxon_id"], shared.taxa["english_name"]):
        path = os.path.join(args.export, "species", f"{tid}.json")
        if (args.taxa and name not in args.taxa) or not os.path.exists(path):
            continue
        ids.append(tid)
        names.append(name)
    jobs = L.taxon_jobs(
        shared, ids, X.load_overrides(), A.load_accounts(), args.data_dir, False, False
    )
    work = []
    for name, job in zip(names, jobs):
        sp = json.load(
            open(os.path.join(args.export, "species", f"{job['taxon']['taxon_id']}.json"))
        )
        if not sp.get("trend"):
            continue
        work.append((name, job, sp["trend"]["first_year"], W.from_dates(sp["window"]["model"])))
    print(f"{len(work)} taxa on {args.workers} workers", flush=True)
    if args.workers > 1:
        with T.worker_pool(args.workers) as ex:
            results = list(ex.map(run_taxon, work))
    else:
        results = [run_taxon(w) for w in work]

    days = pd.concat([r["days"] for r in results], ignore_index=True)
    years = pd.concat([r["years"] for r in results], ignore_index=True)
    os.makedirs(args.out, exist_ok=True)
    days.to_csv(os.path.join(args.out, "gap_fill_days.csv"), index=False)
    years.to_csv(os.path.join(args.out, "gap_fill_years.csv"), index=False)
    pd.set_option("display.width", 200)
    summary(days, years)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
