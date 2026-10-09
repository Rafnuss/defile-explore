# defile-explore

The statistics behind the **Explore** page of [defileViz](https://github.com/Rafnuss/defileViz):
long-term trends, phenology, time of day and the other per-species views of the bird migration
counted at Défilé de l'Écluse (France) since 1966.

It reads the release of the **defile-dataset** (survey and count tables) and writes one JSON file
per taxon, which defileViz renders. The data are counts, and the page shows them as counts first:
model output (a GAM trend, gap-filled totals, smooth phenology) is drawn on top, never instead.

Sister repos:

- **defile-dataset**: the counts, their corrections and documentation (the input here).
- **defile-migration-forecast**: the daily migration forecast (deep learning on weather). Fully
  independent of this repo; see `DECISIONS.md` -> Repository.
- **defileViz**: the website (Vue), which shows both.

## Setup

[uv](https://docs.astral.sh/uv/) with Python 3.11; versions are pinned in `pyproject.toml` and
locked in `uv.lock`.

```bash
uv sync
uv run pytest
```

## Build the export

```bash
# copy a defile-dataset release into data/count/dataset/, then build data/explore/
uv run python scripts/build_explore.py --dataset ../defile-dataset/output

# rebuild from the copied release; --skip-trend for a build without the GAM (seconds)
uv run python scripts/build_explore.py
uv run python scripts/build_explore.py --out ../defileViz/public/data/explore
```

`data/` and `logs/` are generated and never committed.

## Layout

```
src/defile_explore/
  release.py     reading the release tables
  export.py      raw aggregation: daily totals, effort, taxa, reports -> JSON
  timeofday.py   the time-of-day GAM
  profile.py     effort adjustment: time-of-day profile, coverage, annual index
  trend.py       the trend GAM: smooth trend, season, timing shift, gap-filled totals
scripts/
  build_explore.py           the export
  benchmark_trend.py         gap-filling benchmark of the trend model (PDF + CSV)
  analyse_explore_effort.py  comparison of effort normalisations (PDF)
tests/
```

`DECISIONS.md` records what was decided and why, `DEVELOPMENT.md` what is still open.
