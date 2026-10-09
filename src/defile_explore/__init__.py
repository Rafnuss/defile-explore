"""Defile-explore: the data and statistics behind defileViz's Explore page.

- `release`: reading the defile-dataset release tables.
- `export`: a raw aggregation of those tables (daily totals, effort, taxa, reports), written as
  JSON by `scripts/build_explore.py`. No model processing.
- `timeofday`: the time-of-day GAM, a smooth (doy, hour) surface of an hour's share of the day.
- `profile`: the effort adjustment (time-of-day profile, coverage, annual index).
- `trend`: a taxon's smooth trend, season and phenology shift, and gap-filled annual totals, from
  a GAM benchmarked by `scripts/benchmark_trend.py`.

Independent of defile-migration-forecast: `release` and `timeofday` are copies of code there, and
nothing is shared or imported across the two repos. Carrying a change from one to the other is a
manual decision, recorded in DECISIONS.md -> Repository.
"""
