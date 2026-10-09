# Development roadmap

What is still open, in phases. Settled calls go to `DECISIONS.md`, not here. Phase 1 (one
pipeline per taxon, settings, caches, QA viewer) is done: `DECISIONS.md` -> Pipeline.

## Phase 2: review and refine

Go through `logs/viewer/index.html` (`scripts/explore_viewer.py`) in batches: raptors; pigeons,
storks, cranes and other large birds; passerines. A fix is a better rule (preferred) or an
override with its reason. The viewer stays a model-inspection tool (everything that helps judge a
fit, in detail); the visitors' page is designed in defileViz. What the first full build shows:

- Windows (`DECISIONS.md` -> Pipeline): re-run `scripts/benchmark_trend.py` on the extended model
  windows (it uses the default one), and decide how the page says "still passing when counting
  thins out" for the 51 `passage_beyond_counting` taxa.
- `wide_intervals` (12 taxa): median q90/q10 of the gap-filled totals above 2.
- `borrowed_profile` (205 taxa, most of them rare): the time of day of a group, not the taxon's.
- Age: the both-classes rule drops Red Kite, Black Kite, Grey Heron and gull sp. Their juveniles
  are still informative as a minimum share of the birds counted (Red Kite 2024: at least 23%),
  if a panel for that is wanted.
- Time of day uses every day timed to the hour, so for pigeons it pools the hour-by-hour
  notebooks of the late 1960s with recent days. Check that the two agree before pooling them.
- `links`: Vogelwarte and Migration Atlas pages go in `overrides.yaml` by hand.
- Combined series have no links and no written account of their own.
- Written accounts: 35 full-tier taxa have none (most passerines, the unidentified groups, the
  combined series); the season context of the reports (weather, monitoring, results) has no
  account and is no longer exported (`reports.json` is gone).

### Trend reliability

The benchmark runs on every trend taxon and each claim is classed (`DECISIONS.md` -> Pipeline,
`reliability@2`); the viewer flags on every page and figure what defileViz will caveat or hide.
Open:

- Review the thresholds on the viewer (index columns totals / trend / season), now that the
  blown-up bands are fixed (`DECISIONS.md` -> Pipeline, information floor).
- Passerine trends from 2007 rise 20-140x in some taxa (Common Reed Bunting, Common Linnet):
  check that identification and recording effort did not grow with them before showing a trend.
- defileViz: implement the filter from `reliability.elements` and the caveat sentences from the
  reason codes.

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

defileViz's `explore.js` still mentions `reports` in a comment; the species files now carry
`accounts` (French and English) instead. How the export reaches defileViz (copy into `public/data/explore/`, a release asset, or the GCE
host the forecasts use) is open, and so is the license of the exported data.
