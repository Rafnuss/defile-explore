#!/usr/bin/env python3
"""Compares trend models on the benchmark (`defile_explore.benchmark`): the season-only model
beside the GAM, and the ratio estimator, on chosen taxa (defile-migration-forecast issue #55).

The build already benchmarks the GAM on every taxon with a trend (the `benchmark` block, drawn by
the viewer); this script is for trying a model change. It runs on the export's own inputs: the
cached shared stage, and each taxon's start year and model window read from its species file, so
build first (`scripts/build_explore.py`). Writes the trials as CSV to `logs/qa/explore/` and prints
the scores per taxon and method, and the error split by kind of gap.

Usage:
    python scripts/benchmark_trend.py
    python scripts/benchmark_trend.py --taxa "Red Kite" --workers 1
"""

import argparse
import json
import os

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
TAXA = (
    "Black Kite",
    "European Honey Buzzard",
    "Red Kite",
    "Common Buzzard",
    "Osprey",
    "All pigeons (Columba)",
    "Eurasian Chaffinch",
)


def run_taxon(args) -> dict:
    name, job, start, window = args
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
    kappa = T.hour_dispersion(job["days"], job["hourly"], job["effort"], job["profile"])
    fits = {v: T.fit(v, frame, years, doy_range) for v in T.VARIANTS}
    gaps = B.gap_trials(ctx, frame, fits, kappa).assign(taxon=name)
    recent = B.recent_trials(ctx, frame, kappa, T.VARIANTS).assign(taxon=name)
    print(f"  {name} done", flush=True)
    return {"gaps": gaps, "recent": recent}


def split_score(df: pd.DataFrame) -> pd.DataFrame:
    """The gap test's error by kind of gap (`benchmark.KINDS`), per method: the hidden birds' share
    of the year's total, and the median absolute and signed error of their estimate, as a share of
    the year's total."""
    d = df[(df["truth"] > 0) & df["est_hours"].notna()]
    out = {}
    for k in B.KINDS:
        err = (d[f"est_{k}"] - d[f"truth_{k}"]) / d["truth"]
        out[f"{k}_share"] = (d[f"truth_{k}"] / d["truth"]).groupby(d["method"]).median()
        out[f"{k}_abs_err"] = err.abs().groupby(d["method"]).median()
        out[f"{k}_bias"] = err.groupby(d["method"]).median()
    return pd.DataFrame(out)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--data-dir", default=os.path.join(ROOT, "data"))
    ap.add_argument("--export", default=os.path.join(ROOT, "data", "explore"))
    ap.add_argument("--out", default=os.path.join(ROOT, "logs", "qa", "explore"))
    ap.add_argument("--taxa", nargs="+", default=list(TAXA))
    ap.add_argument("--workers", type=int, default=4)
    args = ap.parse_args(argv)

    ebird = pd.read_csv(os.path.join(args.data_dir, E.EBIRD_TAXONOMY_FILE))
    shared = L.load_shared(args.data_dir, ebird)
    by_name = dict(zip(shared.taxa["english_name"], shared.taxa["taxon_id"]))
    ids = [by_name.get(n, n) for n in args.taxa]
    jobs = L.taxon_jobs(
        shared, ids, X.load_overrides(), A.load_accounts(), args.data_dir, False, False
    )
    work = []
    for name, job in zip(args.taxa, jobs):
        sp = json.load(
            open(os.path.join(args.export, "species", f"{job['taxon']['taxon_id']}.json"))
        )
        if not sp["trend"]:
            print(f"  {name}: no trend in the export, skipped")
            continue
        work.append((name, job, sp["trend"]["first_year"], W.from_dates(sp["window"]["model"])))
    if args.workers > 1:
        with T.worker_pool(args.workers) as ex:
            results = list(ex.map(run_taxon, work))
    else:
        results = [run_taxon(w) for w in work]

    os.makedirs(args.out, exist_ok=True)
    gaps = pd.concat([r["gaps"] for r in results], ignore_index=True)
    recent = pd.concat([r["recent"] for r in results], ignore_index=True)
    gaps.to_csv(os.path.join(args.out, "trend_benchmark_gaps.csv"), index=False)
    recent.to_csv(os.path.join(args.out, "trend_benchmark_recent.csv"), index=False)
    pd.set_option("display.width", 200)
    for title, df in (("GAP TRANSPLANT", gaps), ("RECENT SEASONS", recent)):
        s = B.score(df, by=("taxon", "method"))
        print(f"\n{title}\n", s.round(3).to_string())
        print(
            "\n",
            s.groupby("method")[["abs_log_err", "bias", "cover80", "cover95"]].mean().round(3),
        )
    print("\nby kind of gap (share of the year's total)\n", split_score(gaps).round(3).to_string())
    print(f"Wrote the trials to {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
