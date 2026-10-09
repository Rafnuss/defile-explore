#!/usr/bin/env python3
"""Builds the Explore export for defileViz from the defile-dataset release tables.

Reads the release in `data/count/dataset/` (copied there with `--dataset <dir>`) and writes
`data/explore/`: `manifest.json`, `taxa.json`, `effort.json` and `species/<taxon_id>.json`
(`defile_explore.pipeline.build_taxon`), with French names from the eBird taxonomy (downloaded
once, no key) and the written accounts from `content/accounts/`. The shared stage (~80 s of time-of-day profiles) and
each taxon's trend fit are cached in `data/cache/`, so a rebuild after a change to one block
refits nothing. Copy `data/explore/` to defileViz's `public/data/explore/` to publish it.

Usage:
    python scripts/build_explore.py --dataset ../defile-dataset/output   # copy a release, then build
    python scripts/build_explore.py                       # rebuild from data/count/dataset/
    python scripts/build_explore.py --taxa "Red Kite" avibase-ED5A7E8F   # only these species files
    python scripts/build_explore.py --skip-trend          # without the trend model
    python scripts/build_explore.py --no-cache            # recompute everything
"""

import argparse
import os
import shutil
import subprocess
import time
import urllib.request

import pandas as pd

from defile_explore import accounts as A
from defile_explore import export as E
from defile_explore import pipeline as L
from defile_explore import release as R
from defile_explore import settings as X
from defile_explore import trend as T

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def write(path: str, obj) -> int:
    with open(path, "w", encoding="utf-8") as f:
        f.write(E.dumps(obj))
    return os.path.getsize(path)


def copy_release(release: str, data_dir: str) -> None:
    """Copy a defile-dataset release (its `dataset/` tables and `metadata.json`) into
    `<data_dir>/count/dataset/`."""
    folder = os.path.join(data_dir, R.DATASET_DIR)
    os.makedirs(folder, exist_ok=True)
    for f in R.DATASET_FILES:
        shutil.copy2(os.path.join(release, "dataset", f), os.path.join(folder, f))
    shutil.copy2(os.path.join(release, R.METADATA_FILE), os.path.join(folder, R.METADATA_FILE))
    print(f"Copied the release from {release} to {folder}")


def select(taxa: pd.DataFrame, wanted: list[str] | None) -> list[str]:
    """Taxon ids from ids or English names (all taxa when `wanted` is None)."""
    if not wanted:
        return taxa["taxon_id"].tolist()
    by_name = dict(zip(taxa["english_name"], taxa["taxon_id"]))
    ids = [by_name.get(w, w) for w in wanted]
    unknown = [w for w, i in zip(wanted, ids) if i not in set(taxa["taxon_id"])]
    if unknown:
        raise SystemExit(f"Unknown taxa: {unknown}")
    return ids


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--data-dir", default=os.path.join(ROOT, "data"))
    ap.add_argument("--out", default=os.path.join(ROOT, "data", "explore"))
    ap.add_argument("--dataset", help="defile-dataset output folder to copy the release from")
    ap.add_argument("--taxa", nargs="+", help="only these taxa (ids or English names)")
    ap.add_argument(
        "--refresh-names", action="store_true", help="download the eBird taxonomy again"
    )
    ap.add_argument("--skip-trend", action="store_true", help="no trend model (fast build)")
    ap.add_argument("--no-cache", action="store_true", help="ignore and rewrite no cache")
    ap.add_argument("--workers", type=int, default=os.cpu_count() or 1)
    args = ap.parse_args(argv)
    t0 = time.time()

    if args.dataset:
        copy_release(args.dataset, args.data_dir)
    names_path = os.path.join(args.data_dir, E.EBIRD_TAXONOMY_FILE)
    if args.refresh_names or not os.path.exists(names_path):
        urllib.request.urlretrieve(E.EBIRD_TAXONOMY_URL, names_path)
        print(f"Downloaded the eBird taxonomy (French names) to {names_path}")
    ebird = pd.read_csv(names_path)

    shared = L.load_shared(args.data_dir, ebird, use_cache=not args.no_cache)
    print(f"Shared stage ready in {time.time() - t0:.0f} s (cache key {shared.key})")
    taxa, days, effort = shared.taxa, shared.days, shared.effort
    overrides = X.load_overrides()
    accounts = A.load_accounts(taxon_ids=taxa["taxon_id"])
    ids = select(taxa, args.taxa)
    jobs = list(L.taxon_jobs(shared, ids, overrides, accounts, args.data_dir, not args.no_cache))
    if args.skip_trend:
        for j in jobs:
            j["overrides"] = j["overrides"] | {
                j["taxon"]["taxon_id"]: {"trend": False, "reason": "--skip-trend"}
            }
    git_sha = subprocess.run(
        ["git", "rev-parse", "--short", "HEAD"], capture_output=True, text=True, cwd=ROOT
    ).stdout.strip()

    if not args.taxa and os.path.isdir(args.out):
        shutil.rmtree(args.out)
    os.makedirs(os.path.join(args.out, "species"), exist_ok=True)
    sizes = {}
    if not args.taxa:
        sizes = {
            "manifest.json": write(
                os.path.join(args.out, "manifest.json"),
                E.manifest(shared.metadata, days, effort, git_sha)
                | {"accounts_sha256": A.fingerprint()},
            ),
            "taxa.json": write(os.path.join(args.out, "taxa.json"), E.records(taxa)),
            "effort.json": write(
                os.path.join(args.out, "effort.json"),
                {"days": E.columns(effort), "annual": E.records(E.annual_effort(effort))},
            ),
        }
    species_bytes, n_trend, flags = 0, 0, {}
    with T.worker_pool(min(args.workers, len(jobs))) as ex:
        for taxon_id, out in ex.map(L.build_taxon, jobs, chunksize=4):
            species_bytes += write(os.path.join(args.out, "species", f"{taxon_id}.json"), out)
            n_trend += out["trend"] is not None
            for f in out["diagnostics"]["flags"]:
                flags[f] = flags.get(f, 0) + 1
    sizes[f"species/ ({len(jobs)} files)"] = species_bytes

    print(f"Explore export -> {args.out} in {time.time() - t0:.0f} s")
    for name, n in sizes.items():
        print(f"  {name:28s} {n / 1e3:9.1f} kB")
    print(f"  trends: {n_trend} taxa, to {shared.last_year}")
    print("  flags:", flags)
    if not args.taxa:
        print("  tiers:", taxa["tier"].value_counts().to_dict())
        print("  profiles:", taxa["profile"].value_counts().to_dict())
        missing = taxa.loc[taxa["french_name"].isna(), "english_name"].tolist()
        if missing:
            print(f"  no French name ({len(missing)}): {', '.join(missing)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
