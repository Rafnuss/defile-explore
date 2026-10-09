# AGENTS.md

Guidance for AI coding agents working in this repo.

## What this project does

Computes the statistics behind defileViz's Explore page from the defile-dataset release: a JSON
file per taxon (daily counts, effort, time-of-day profile, GAM trend and gap-filled totals, ...)
plus `taxa.json`, `effort.json`, `reports.json` and `manifest.json`. `README.md` has the layout
and the commands.

Principles:

- **Empirical first.** Every panel shows the counts; a model (the trend GAM, the time-of-day
  GAM) is drawn on top and labelled as such.
- **One fit per taxon.** Panels derive from the taxon's fit and its data; nothing refits a model
  to draw a panel.
- **Effort is explicit.** A count is read with its coverage `c` (the share of the day's expected
  passage counted). No adjusted value below `COVERAGE_MIN`, nothing adjusted before a taxon's
  start year.

## Independence from the forecast

defile-migration-forecast and this repo share no code and import nothing from each other.
`release.py` and `timeofday.py` are copies of forecast code (origin in their docstrings). Never
add a dependency on the forecast repo or a shared package; carrying a change across is a decision
for the user, logged in the table in `DECISIONS.md` -> Repository.

## The output is a contract with defileViz

defileViz (`src/services/explore.js`, `src/components/explore/`) reads `taxa.json` and
`species/<taxon_id>.json`, including `trend.annual`, `trend.passage`, `trend.season`,
`trend.episodes`, `trend.theta`, `trend.kappa`, `profile` and `days`. Renaming or reshaping a
field breaks the page: change both repos together.

## Conventions

- Fixed values are module-level `UPPER_SNAKE_CASE` constants, defined once in the module that owns
  the concept and imported elsewhere. Grep before writing a literal.
- `DECISIONS.md` is the log of settled calls and of what was tried and rejected: read it before
  re-proposing something, and add to it when a question is settled. `DEVELOPMENT.md` holds only
  what is open.
- Before treating a change as done: `uv run pytest` and `pre-commit run --all-files` (black and
  isort at line length 99, docformatter, mdformat, codespell).
- `data/` and `logs/` are generated: never commit them. The release is copied into
  `data/count/dataset/` by `scripts/build_explore.py --dataset <dir>`.
- Builds are verbose (pygam warnings, per-taxon progress): redirect to `logs/` and grep.

## Environment gotcha

This repo sits in a OneDrive-synced folder. A cloud-only placeholder can make a read hang at 0%
CPU or fail with "Resource deadlock avoided"; if a build or `uv run` stalls, suspect that first
and ask the user to make the folder available offline.
