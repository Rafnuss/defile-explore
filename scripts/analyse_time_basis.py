#!/usr/bin/env python3
"""Held-out test of the time-of-day profile's time axis: solar hours against day-length units.

`solar` is the profile as built until 2026-10-09 (the ratio GAM on solar hours, cut by daylight);
`tau_day16` is the profile as built since (`profile.fit_profile`). `tau` fits the same GAM on tau = (t - civil dawn) / (civil dusk - civil dawn), the
share of the day's light elapsed, so dawn and dusk are the surface's edges on every date; it is
read back on the same 5-minute solar slots and cut by daylight the same way. `solar_day` and
`tau_day` (`tau_day16`: 16 splines on the time axis) fit only the samples whose hour overlaps civil
daylight, so a few night hours counted on 24-hour days no longer stretch the axis; `tau_day`'s
axis is exactly civil dawn (0) to civil dusk (1). `uniform` is the floor. Each taxon with its
own profile is fitted on even years and scored on odd years' timed days, then the other way round.

Scores, on each held-out day timed to the hour (`PROFILE_MIN_TIMED`) with birds in at least two
counted hours, from the profile's expected share of the day's counted passage in each counted
hour (`by_hour(p * slots)`, as coverage reads it):

- `l1`: half the L1 distance between the observed and expected shares, days weighted by their
  birds ** `RATIO_WEIGHT_POWER` (the forecast's held-out test of the same GAM);
- `edge`: the same, on the first and the last counted hour only, on days counted from within an
  hour of dawn or to within an hour of dusk: where the two axes differ most;
- `fill`: |log| error of the birds of the second half of the counted hours predicted from the
  first half (birds kept / coverage kept x coverage hidden), +1 on both, days weighted as `l1`:
  gap filling, what coverage is for.

Writes `logs/qa/explore/time_basis.csv` (one row per taxon, fold and axis) and prints a summary.

Usage:
    python scripts/analyse_time_basis.py [--taxa "Red Kite" ...] [--workers 8]
"""

import argparse
import multiprocessing
import os
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd

from defile_explore import export as E
from defile_explore import pipeline as L
from defile_explore import profile as P
from defile_explore import trend as T
from defile_explore.timeofday import RATIO_WEIGHT_POWER, fit_ratio_surface

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BASES = ("solar", "tau", "solar_day", "tau_day", "tau_day16")
FOLDS = (0, 1)  # fitted on years of this parity, scored on the others
TAU_GRID = 241  # points of the tau axis the GAM is read at, then interpolated to slots
EDGE_HOURS = 1.0  # a day counted from within this of dawn (or to within it of dusk)
MIN_HOURS = 2  # a scored day has birds and at least this many counted hours
SHARE_FLOOR = 1e-4  # an expected share below this is taken as this (log scores)

SHARED = None  # the shared stage, loaded once in each worker (`load`)


def load(data_dir: str) -> None:
    """Worker initializer: the shared stage, from its cache."""
    global SHARED
    ebird = pd.read_csv(os.path.join(data_dir, E.EBIRD_TAXONOMY_FILE))
    SHARED = L.load_shared(data_dir, ebird)


def fit_tau(samples, doy_grid, light, dawn, dusk, label=""):
    """P(t | doy) on solar slots from the ratio GAM fitted on tau, cut by daylight."""
    row = np.clip(samples["doy"].to_numpy(), *P.PROFILE_DOY) - P.PROFILE_DOY[0]
    length = dusk - dawn
    tau = (samples["hour"].to_numpy() + 0.5 - dawn[row]) / length[row]
    grid = np.linspace(tau.min(), tau.max(), TAU_GRID)
    ratio = fit_ratio_surface(
        samples["doy"].to_numpy(),
        tau,
        samples["ratio"].to_numpy(),
        samples["weight"].to_numpy(),
        doy_grid,
        grid,
        label=label,
    )
    t = (np.arange(E.SLOTS) + 0.5) / E.STEPS_PER_HOUR
    p = np.zeros((len(doy_grid), E.SLOTS))
    for i in range(len(doy_grid)):
        p[i] = np.interp((t - dawn[i]) / length[i], grid, np.clip(ratio[i], 0, None))
    p *= light
    return p / p.sum(axis=1, keepdims=True)


def in_daylight(samples, dawn, dusk):
    """The samples whose hour overlaps civil daylight on its day (`dawn`, `dusk` on
    `PROFILE_DOY`)."""
    row = np.clip(samples["doy"].to_numpy(), *P.PROFILE_DOY) - P.PROFILE_DOY[0]
    h = samples["hour"].to_numpy()
    return samples[(h + 1 > dawn[row]) & (h < dusk[row])]


def fit_solar(samples, doy_grid, light, label=""):
    """The profile as built until 2026-10-09: the ratio GAM on solar hours, each sample at its
    hour's midpoint, read at each slot's midpoint (flat beyond the first and last hour sampled),
    cut by daylight."""
    mid = samples["hour"].to_numpy() + 0.5
    t = (np.arange(E.SLOTS) + 0.5) / E.STEPS_PER_HOUR
    span = (t > samples["hour"].min()) & (t < samples["hour"].max() + 1)
    ratio = fit_ratio_surface(
        samples["doy"].to_numpy(),
        mid,
        samples["ratio"].to_numpy(),
        samples["weight"].to_numpy(),
        doy_grid,
        np.clip(t[span], mid.min(), mid.max()),
        label=label,
    )
    p = np.zeros((len(doy_grid), E.SLOTS))
    p[:, span] = np.clip(ratio, 0, None)
    p *= light
    return p / p.sum(axis=1, keepdims=True)


def fit_solar_day(samples, doy_grid, light, dawn, dusk, label=""):
    """`fit_solar` on the samples within daylight only: the hour axis spans the counted day."""
    return fit_solar(in_daylight(samples, dawn, dusk), doy_grid, light, label=label)


def fit_tau_day(samples, doy_grid, light, dawn, dusk, label="", k1=12):
    """The ratio GAM on tau within daylight, its axis exactly civil dawn (0) to dusk (1)."""
    s = in_daylight(samples, dawn, dusk)
    row = np.clip(s["doy"].to_numpy(), *P.PROFILE_DOY) - P.PROFILE_DOY[0]
    length = dusk - dawn
    # the hour's midpoint, kept within its daylight part
    lo = np.maximum(s["hour"].to_numpy(), dawn[row])
    hi = np.minimum(s["hour"].to_numpy() + 1, dusk[row])
    tau = ((lo + hi) / 2 - dawn[row]) / length[row]
    grid = np.linspace(0, 1, TAU_GRID)
    ratio = fit_ratio_surface(
        s["doy"].to_numpy(),
        tau,
        s["ratio"].to_numpy(),
        s["weight"].to_numpy(),
        doy_grid,
        grid,
        k1=k1,
        label=label,
    )
    t = (np.arange(E.SLOTS) + 0.5) / E.STEPS_PER_HOUR
    p = np.zeros((len(doy_grid), E.SLOTS))
    for i in range(len(doy_grid)):
        p[i] = np.interp((t - dawn[i]) / length[i], grid, np.clip(ratio[i], 0, None))
    p *= light
    return p / p.sum(axis=1, keepdims=True)


def score(profile, days, hourly, effort, dawn, dusk) -> dict:
    """`l1`, `edge` and `fill` (module docstring) of `profile` on the timed days in `days`."""
    d = days[(days["count"] > 0) & (days["timed"] >= P.PROFILE_MIN_TIMED)]
    eff = effort[effort["state"] == "counted"].set_index("date")
    d = d[d["date"].isin(eff.index)]
    birds = {k: g.groupby("hour")["count"].sum() for k, g in hourly.groupby("date")}
    out = {"l1": [], "l1_w": [], "edge": [], "edge_w": [], "fill": [], "fill_w": []}
    for date in d["date"]:
        if date not in birds:
            continue
        slots = eff.at[date, "slots"]
        hours = E.by_hour(slots) / E.STEPS_PER_HOUR
        b = np.zeros(E.HOURS)
        h = birds[date]
        b[h.index.to_numpy()] = h.to_numpy()
        on = hours > 0
        b = np.where(on, b, 0)
        total = b.sum()
        if on.sum() < MIN_HOURS or total <= 0:
            continue
        row = int(np.clip(date.dayofyear, *P.PROFILE_DOY)) - P.PROFILE_DOY[0]
        expected = E.by_hour(profile[row] * slots)
        if expected.sum() <= 0:
            continue
        share = expected / expected.sum()
        w = total**RATIO_WEIGHT_POWER
        out["l1"].append(0.5 * np.abs(b / total - share).sum())
        out["l1_w"].append(w)
        counted = np.flatnonzero(on)
        first, last = counted[0], counted[-1]
        if first <= dawn[row] + EDGE_HOURS or last + 1 >= dusk[row] - EDGE_HOURS:
            out["edge"].append(np.abs(b / total - share)[[first, last]].sum() / 2)
            out["edge_w"].append(w)
        half = counted[: len(counted) // 2]
        rest = counted[len(counted) // 2 :]
        kept, hidden = expected[half].sum(), expected[rest].sum()
        if kept > 0 and b[half].sum() > 0:
            pred = b[half].sum() / kept * hidden
            out["fill"].append(abs(np.log((pred + 1) / (b[rest].sum() + 1))))
            out["fill_w"].append(w)
    res = {}
    for k in ("l1", "edge", "fill"):
        v, w = np.array(out[k]), np.array(out[k + "_w"])
        res[k] = float((v * w).sum() / w.sum()) if len(v) else np.nan
        res[f"n_{k}"] = len(v)
    return res


def run_taxon(args) -> list[dict]:
    taxon_id, name = args
    s = SHARED
    t0 = time.time()
    days = s.days[s.days["taxon_id"] == taxon_id]
    hourly = s.hourly[s.hourly["taxon_id"] == taxon_id]
    doy_grid = np.arange(P.PROFILE_DOY[0], P.PROFILE_DOY[1] + 1)
    light = P.daylight(doy_grid)
    dawn, dusk = P.twilight(doy_grid)
    rows = []
    for fold in FOLDS:
        fit_years = days["date"].dt.year % 2 == fold
        test_years = hourly["date"].dt.year % 2 != fold
        samples = P.profile_samples(
            days[fit_years], hourly[hourly["date"].dt.year % 2 == fold], s.effort
        )
        if len(samples) == 0:
            continue
        fits = {
            "solar": fit_solar(samples, doy_grid, light, label=name),
            "tau": fit_tau(samples, doy_grid, light, dawn, dusk, label=name),
            "solar_day": fit_solar_day(samples, doy_grid, light, dawn, dusk, label=name),
            "tau_day": fit_tau_day(samples, doy_grid, light, dawn, dusk, label=name),
            "tau_day16": P.fit_profile(samples, doy_grid, light, name, (dawn, dusk)),
            "uniform": P.uniform_profile(light),
        }
        for basis, profile in fits.items():
            sc = score(profile, days[~fit_years], hourly[test_years], s.effort, dawn, dusk)
            rows.append({"taxon": name, "fold": fold, "basis": basis, **sc})
    print(f"  {name}: {len(samples)} samples, {time.time() - t0:.0f} s", flush=True)
    return rows


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--data-dir", default=os.path.join(ROOT, "data"))
    ap.add_argument("--out", default=os.path.join(ROOT, "logs", "qa", "explore"))
    ap.add_argument("--taxa", nargs="+", help="English names (default: every own profile)")
    ap.add_argument("--workers", type=int, default=8)
    args = ap.parse_args(argv)

    load(args.data_dir)  # builds the cache if missing, before the workers read it
    taxa = SHARED.taxa.set_index("taxon_id")
    own = SHARED.source.index[SHARED.source == "own"]
    names = taxa.loc[own, "english_name"]
    if args.taxa:
        names = names[names.isin(args.taxa)]
    work = list(names.items())
    print(f"{len(work)} taxa on {args.workers} workers", flush=True)
    for var in T.WORKER_THREAD_VARS:  # one BLAS thread per worker, as `trend.worker_pool`
        os.environ[var] = "1"
    ctx = multiprocessing.get_context("spawn")
    pool = ProcessPoolExecutor(args.workers, ctx, initializer=load, initargs=(args.data_dir,))
    with pool as ex:
        rows = [r for rs in ex.map(run_taxon, work) for r in rs]

    df = pd.DataFrame(rows)
    os.makedirs(args.out, exist_ok=True)
    df.to_csv(os.path.join(args.out, "time_basis.csv"), index=False)
    wide = df.pivot_table(index=["taxon", "fold"], columns="basis", values=["l1", "edge", "fill"])
    print(f"{df['taxon'].nunique()} taxa, {len(wide)} taxon-folds; lower is better")
    for k in ("l1", "edge", "fill"):
        w = wide[k].dropna()
        for b in w.columns.drop("solar"):
            diff = w[b] - w["solar"]
            print(
                f"  {k:5s} {b:10s} {w[b].mean():.4f} (solar {w['solar'].mean():.4f})  paired "
                f"median {diff.median():+.4f}  better {(diff < 0).sum()} / worse {(diff > 0).sum()}"
            )
    for b in ("solar_day", "tau_day"):
        per_taxon = (wide["edge"][b] - wide["edge"]["solar"]).groupby("taxon").mean().sort_values()
        print(f"  edge, {b} - solar, by taxon (mean of folds): best for {b}")
        print(per_taxon.head(8).round(4).to_string())
        print("  ... worst")
        print(per_taxon.tail(4).round(4).to_string())


if __name__ == "__main__":
    main()
