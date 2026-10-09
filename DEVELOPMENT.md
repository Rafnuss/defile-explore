# Development roadmap

What is still open, in phases. Settled calls go to `DECISIONS.md`, not here. Phase 1 (one
pipeline per taxon, settings, caches, QA viewer) is done: `DECISIONS.md` -> Pipeline.

## Phase 2: review and refine

Go through `logs/viewer/index.html` (`scripts/explore_viewer.py`) in batches: raptors; pigeons,
storks, cranes and other large birds; passerines. A fix is a better rule (preferred) or an
override with its reason. What the first full build already shows:

- `passage_at_window_edge` (15 taxa): the smooth 10% or 90% passage date within 5 days of the
  window's ends (7 Jul, 18 Nov). Late species (Red Kite, Common Buzzard, Hen Harrier) are cut off.
  A window per taxon changes the benchmarked fit, so it needs a re-benchmark and a decision.
- `wide_intervals` (11 taxa): median q90/q10 of the gap-filled totals above 2.
- `borrowed_profile` (205 taxa, most of them rare): the time of day of a group, not the taxon's.
- Age: whether aged birds represent those passing. Red Kite is 97% non-adult, carried by
  2024-2025 when many juveniles were aged; the usable-year rule (20 birds, 5% of those counted, 3
  years) may need to weigh years equally or ask for a steadier aged share.
- Records and the record day include years before the start year (all pigeons: 116 340 on
  20 Oct 1975). Keep the all-time record, or count from the start year?
- Time of day uses every day timed to the hour, so for pigeons it pools the hour-by-hour
  notebooks of the late 1960s with recent days. Check that the two agree before pooling them.
- `links`: Vogelwarte and Migration Atlas pages go in `overrides.yaml` by hand.
- Combined series have no links and take the report paragraphs of their first member.

## Phase 3: methods

defileViz's `ExploreMethod.vue` (a modal with KaTeX equations and widgets drawn from the selected
species) is where readers get the methods. To grow it:

- Split it into a shell plus one section component per method (coverage, model, fill, ...,
  then season, time of day, age/sex, weather).
- Each JSON block names its method and version; the plot's info button opens that section.
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

## Publishing

How the export reaches defileViz (copy into `public/data/explore/`, a release asset, or the GCE
host the forecasts use) is open, and so is the license of the exported data.
