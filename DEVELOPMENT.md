# Development roadmap

What is still open, in phases. Settled calls go to `DECISIONS.md`, not here. Phase 1 (one
pipeline per taxon, settings, caches, QA viewer) is done: `DECISIONS.md` -> Pipeline.

## Phase 2: review and refine

Go through `logs/viewer/index.html` (`scripts/explore_viewer.py`) in batches: raptors; pigeons,
storks, cranes and other large birds; passerines. A fix is a better rule (preferred) or an
override with its reason. The viewer stays a model-inspection tool (everything that helps judge a
fit, in detail); the visitors' page is designed in defileViz. What the first full build shows:

- Windows (`DECISIONS.md` -> Pipeline): re-run `scripts/benchmark_trend.py` on the extended model
  windows (it uses the default one).
- `wide_intervals` (12 taxa): median q90/q10 of the gap-filled totals above 2.
- `borrowed_profile` (205 taxa, most of them rare): the time of day of a group, not the taxon's.
- Age and sex (`demography@4`): Red Kite's single pooled year (2018, 489 birds, 66% young) and
  those of Great Cormorant, Grey Heron and Eurasian Sparrowhawk rest on one season: check they are
  representative. Black Kite and Common Buzzard still get nothing (juveniles tagged, adults left
  blank every year).
- Time of day uses every day timed to the hour, so for pigeons it pools the hour-by-hour
  notebooks of the late 1960s with recent days. Check that the two agree before pooling them.
- Entry errors (`scripts/check_entries.py`, `logs/qa/explore/entry_errors.csv`): correct them at
  the source (defile-dataset, Trektellen), then rebuild. Mornings of the days summer time ends
  (2017, 2019, 2020) are timed before dawn; 6 Nov 2021 has a block's birds at its last minutes.
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
- Filled days from the hours within the day: `trend.days` uses only a day's count and coverage.
  Measure first whether passage is correlated from hour to hour beyond the profile (residuals of
  neighbouring hours on timed days; in the gap benchmark, does the last counted hour predict the
  hidden ones better than count / c?). Only if it is, a field over (hour, day) per year.

## Phase 3: page and methods

The species picker: `taxa.json` carries `taxon_order` (eBird's sequence), `group` (raptors,
waterbirds, passerines, other, by order), `season_birds` (median birds counted a season over the
last 10), `story` (`full`, `caveat`, `counts`: what the trend's totals class lets the page show)
and `highlight`, all from `catalogue.py`. The highlights are curated in `content/highlights.tsv`,
seeded with the 32 species of the 2019 and 2020 papers that have a page. Open: whether to add the
most counted passerines the papers leave out (Eurasian Chaffinch, Common Starling, Barn Swallow,
Western House Martin, Common Swift, the finches); the group names; whether "other" should split
off the pigeons.

The page's structure is drafted in both repos (`DECISIONS.md` -> Pipeline, page and methods), for
review. Open:

- Review the page and the five modals species by species (raptors, a flocking species with hidden
  totals such as Purple Heron, a taxon without a trend such as finch sp., one with age and sex
  such as Western Marsh Harrier), then the French and German labels.
- The modals are English only, and quote the benchmark's medians by hand (`MethodReliability.vue`
  `STATS`): written into `manifest.json` by the build instead, if they drift too often.
- `ExploreMethod.vue.orig` in defileViz is the single modal before the split, kept for the review:
  delete it once the split is accepted.
- The timing panel could select a season the same way as the totals (a click on its row).

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

The species files grew from 8 to 18 MB in all (benchmark trials, episodes): check what the page
needs before publishing. The license of the exported data is open (defile-dataset's is too).
