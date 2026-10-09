# Development roadmap

What is still open, in phases. Settled calls go to `DECISIONS.md`, not here.

## Phase 1: one pipeline per taxon

Every Explore panel comes from one per-taxon function, so each taxon is fitted once and every
panel is derived from that fit. A prototype of all candidate panels (6 species, Plotly, local)
fixed which panels to build: key numbers, trend, daily phenology heatmap with passage dates,
chances, time of day, weather, age/sex, records.

- `settings.py`: the settings of a taxon, each with a default computed by rule (start year, tier
  and profile source already are) and exceptions in `overrides.yaml`, each with a reason. The JSON
  records where each value came from.
- `pipeline.py`: `build_taxon(shared, settings)` returns the taxon's blocks: `trend`
  (`taxon_trend`, plus the smooth 10/50/90% passage dates and the key numbers), `season` (daily
  heatmap, each year's passage dates, chances), `daytime`, `demography` (age, sex), `records`.
- JSON contract for `species/<taxon_id>.json`: `taxon`, `settings`, `blocks` (each naming its
  method and version, which maps to a section of defileViz's method modal), `diagnostics`.
- Speed: the shared stage (release, effort, profiles) cached on disk by release version, the fit
  cached per taxon by a hash of its data and settings, `--taxa` and `--only <block>`. Targets: a
  cold build under 6 min on 11 cores, one taxon and one block under 15 s. The time-of-day
  profiles (~80 s, serial) are the first thing to cache.
- `scripts/explore_viewer.py`: the prototype turned into a viewer of the JSON, with no maths of
  its own, so the contract is tested before defileViz uses it.

## Phase 2: review and refine

- A QA index of every taxon sorted by diagnostics: low coverage, a borrowed profile, passage
  running into the window's edge, wide intervals, few aged birds.
- Review in batches (raptors; pigeons and other large birds; passerines). A fix is a better rule
  (preferred) or an override with its reason.
- Per-species windows: late species (Red Kite, Common Buzzard) run into the 18 Nov window end.
  A window per taxon changes the benchmarked fit, so it needs a re-benchmark and a decision.
- Age: the codes switch between `1` and `I` across years (Red Kite 2019 vs 2022), so non-adults
  are one class. The usable years are a manual guess for now; a rule (share of birds aged) is
  to be found.

## Phase 3: methods

defileViz's `ExploreMethod.vue` (a modal with KaTeX equations and widgets drawn from the selected
species) is where readers get the methods. To grow it:

- Split it into a shell plus one section component per method (coverage, model, fill, ...,
  then phenology, time of day, age/sex, weather).
- Each JSON block names its method id and version; the plot's info button opens that section.
- The benchmark's numbers (gap error, interval coverage) written by `scripts/benchmark_trend.py`
  into `manifest.json` and read by the modal, instead of copied into its text.
- `DECISIONS.md` stays the developers' record of why; the modal is the readers' account.

## Phase 4: weather

The prototype binned days of the main passage by one weather variable (12 variables: tailwind
and crosswind, wind speed, rain on the day and the day before, cloud, low cloud, temperature
anomaly and change, pressure change, sunshine, upstream tailwind) and showed birds relative to
the GAM's expectation without weather episodes, with a keep flag (best and worst bins' 95%
intervals apart). Hard to read, and the keep rule passes too much for rare taxa.

- Design first: the binned ratio vs. weather terms in the GAM; which 6-10 variables; a keep rule
  with a minimum effect and a minimum number of birds; a simple "good day / bad day" summary.
- The baseline is already exported: `trend.episodes.base` (expected birds without episodes).
- Weather data: Explore needs daily ERA5 at Défilé (and one upstream site). Getting it is a
  decision of its own, since nothing is shared with the forecast's weather cache.
