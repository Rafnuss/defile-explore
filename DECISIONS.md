# Défilé Explore — decisions log

A running log of settled calls for the Explore statistics: what was tried, what worked, what was
rejected and why. Not a design spec; read the code for that. See `DEVELOPMENT.md` for what is
still open.

## Repository

**Narrative accounts belong in Explore (2026-10-09).** The bilingual species and species-year accounts moved from defile-dataset to `content/accounts/`, with their authored TSV inputs, CSV exports, reading copies, editorial coverage and generation instructions. Their scripts live in `scripts/accounts/` and read the unchanged release tables in `data/count/dataset/`, including the original annual and paper text extracts. Published annual totals and the original-name crosswalk are retained as supporting snapshots under `content/accounts/reference/`. Authored content is versioned outside ignored `data/`; exporting it does not change the statistical JSON contract or generate new prose.

**Explore has its own repo, split from defile-migration-forecast (2026-10-09).** It started there
as `src/explore/` (that repo's issue #55) and grew into its own project: per-species statistics,
a QA viewer and method pages, with none of the forecast's dependencies (torch, Lightning, Hydra).
The history of its files came along (`git filter-repo`), so the commits below the split are the
forecast repo's, with paths renamed from `src/explore/` to `src/defile_explore/`.

**Nothing is shared with the forecast, by design.** Two pieces of forecast code are copied here,
each with its origin in its docstring: the release reader (`release`, from `src/data/counts.py`)
and the time-of-day GAM (`timeofday`, from `src/phenology.py`), both at forecast commit 13c5107.
In the forecast repo, Explore's results were meant to cross over only as named files (a profile,
a population trend); that still holds, now across repos. Carrying a change from one copy to the
other is a decision made by hand and logged here, never an import or a shared package.

| Date       | Direction           | What                                | Why       |
| ---------- | ------------------- | ----------------------------------- | --------- |
| 2026-10-09 | forecast -> explore | `release.py`, `timeofday.py` copied | the split |

## Explore

**The Explore export is a separate, raw aggregation of the release** (`src/defile_explore/`,
defile-migration-forecast #55), not the model's counts: no count is moved, dropped or imputed, and daily totals reconcile
with `count.csv` and the dataset's own daily totals (`tests/test_explore.py`). Effort is the union
of `complete` survey intervals per local day. `species_doy_statistics.json` stays the model's
contract, unchanged.

**Effort-adjusted values use a time-of-day profile, not birds per hour.** A taxon's profile
p(h | doy) is fitted with the model's own ratio GAM (`defile_explore.timeofday.fit_ratio_surface`, one
implementation) on solar hours (local clock hours until 2026-10-09, see below), from days timed to the hour; a counted day's coverage `c` is
the share of the profile in the hours counted, the adjusted day is count / c, and the annual index
is Σ birds / Σ c over counted days in the window (a ratio estimator, so low-coverage days do not
dominate). Taxa with fewer than 500 timed birds use a group profile (raptors, pigeons,
passerines, other), then a uniform one. On the 2026-10 release, 66 taxa have their own profile.
`scripts/analyse_explore_effort.py` compares the options for Black Kite, Honey Buzzard, Red Kite,
Common Buzzard, Wood Pigeon and Chaffinch:

- Birds per hour and the uniform (daylight) index are the same curve. Neither is offered.
- The profile matters where the species has a peaked day: before 1993 Common Buzzard is at 0.25
  of its 2010-2025 level per hour and 0.82 with its profile (midday passage, morning counts).
- Profiles from the hourly sheets (2014-2020) and from Trektellen entry times (2021-2026) agree
  (total variation 0.05-0.12; Wood Pigeon 0.22, shaped by a few huge mornings), so Trektellen
  times are usable as passage times. Nothing is timed to the hour before 2014, so stability across
  decades cannot be tested.

**No adjusted value below c = 0.5** (`COVERAGE_MIN`), for a day or for a year's mean: the raw
count and `c` are shown instead. Pre-1993 counts cover a median 0.05-0.11 of a raptor's day, and
the step at 1993 survives the adjustment (Kestrel x14, Honey Buzzard x5; Common Buzzard x2.4 ->
x1.45), so effort is not what separates those years. They stay visible, unadjusted, and out of
reference bands. Taxa not counted systematically in some years (passerines before ~2007) are a
protocol matter no effort metric fixes; the dataset README has that history.

**Explore lives in this repo, contained in `src/defile_explore/`.** It shares the forecast's release
reader and time-of-day GAM, and a population trend fitted for Explore may later help the
forecast (Red Kite's under-prediction of recent years is a missing trend). But the forecast never
imports it (`tests/test_explore.py` enforces the direction); a result crosses over only as a named
file the forecast reads, as `species_doy_statistics.json` does, decided here first. The GAM
profile is variant A of the Explore baseline, provisional until compared with a hierarchical
Gaussian process on a chronological holdout.

**Each taxon has a start year, 1993 or 2007**: its counts are compared from then on, and nothing
is adjusted before it (earlier years stay in the export, raw). Daily systematic counting began in
1993 and always targeted raptors, herons and egrets, storks, pigeons and corvids, and other large
birds counted individually (cranes, geese, ducks), which start in 1993 even if rare then (Peregrine
recovered). Passerines were hardly recorded before 2007, when their taxa double and their birds
rise 15-fold (defile-dataset `docs/sampling-history.md`): any other taxon starts in 1993 only if
recorded in at least 0.75 times as large a share of the 1993-2006 years as of the later ones. A
year counts as recorded only with at least 2% of the taxon's median year since 2007, so that a
trickle noted while the taxon was not counted is not a series ("swallow sp." 2000-2006: ~600 birds
a year against ~150 000 since). The floor stays low because it also penalises a real increase
(Common Crane, ~30 a year before 2007 and 380 since). Full tier: 45 from 1993, 41 from 2007.

**Combined series where names were split differently over the years** (`COMBINED`): all
Columba pigeons (Wood Pigeon, Stock Dove, "Columba sp."; not the local Feral and Rock Pigeons), and
all swallows and martins. "Columba sp." appears only in 2014, when the hourly sheets begin, and has
been a quarter of the pigeons since, so Wood Pigeon's own series drops there for a recording
reason; swallows were "swallow sp." in 1993-1999 and are increasingly identified since 2021. A
combined series is exported like a taxon (rank `combined`, its `members`, its own profile and start
year: pigeons 1993, swallows 2007), beside its members, which stay as they are. Its trend is the
one to read.

**Tiers count days with migrating birds**, not birds: `full` from 50 days over 5 years, `rare`
at 10 days or fewer (84 / 45 / 142). Provisional, to tune once the page exists. French names come
from the eBird taxonomy (fr_FR) by `ebird_code`.

**Trends and annual totals come from a GAM, not a Gaussian process** (`src/defile_explore/trend.py`,
`scripts/benchmark_trend.py`). A day's count is negative binomial around coverage x exp(intercept
\+ trend(year) + season(doy) + shift(year, doy) + year level + episode), where `episode` is each
year's own short-range curve (runs of good or bad migration days, the weather's share); hours not
counted are filled with an hourly over-dispersion `kappa` fitted from the timed days (flocks). The
annual total is the birds counted plus the posterior predictive of the coverage missed. Both
priors were fitted on one engine and benchmarked on Black Kite, Honey Buzzard, Red Kite, Common
Buzzard, Osprey, all pigeons and Chaffinch:

- Gap filling (well-counted years 2014-2025 given an old year's gaps): GAM and GP both 5.0%
  error, the GP better in 49% of paired trials; the ratio index 14.6% (40% too high for pigeons),
  the season alone 5.7%.
- Calibration came from the episode term and `kappa`, not from the prior: without them both 80%
  intervals held the truth 54% of the time, with them 80% (95%: 94%).
- Predicting 2023-2025 from the years before (extrapolation, which Explore never does): GP 35%
  error, GAM 49%, almost all of it Chaffinch, where the P-spline extends its last slope; daily log
  scores equal (-3.170 / -3.171). A trend fed to the forecast would need a flat extrapolation.
- The Hilbert-space GP misbehaved at long lengthscales (prior variance 0.17 instead of 1 at the
  series' ends for a 30-year lengthscale, with the domain at 1.5 x the data), and full Bayes
  (NumPyro NUTS: 3-5 min a fit; PyMC's C backend does not compile here) bought nothing once
  intervals were calibrated. The GAM refits in 0.2 s.

Every smooth sums to zero over its grid, so the intercept, trend, season, `shift` (an interaction
only, as mgcv's `ti`) and the year level each own their part: unconstrained, the smooth trend moved
with the optimiser's stopping point (Black Kite 2025/1993: 0.80 to 1.21). Smoothing parameters
come from the Fellner-Schall update (mgcv's `efs`): the same optimum from any start. `theta` then
alternates with them as the maximum of the Laplace marginal likelihood: its likelihood given the
fitted means was 15-35% too large (days looked less variable than they are). About 3 s a taxon
(median 2.9 s of CPU over the 74 exported trends).

The constraints separate the terms on the log scale only, so `trend` and `shift` are not
abundance and timing on their own, and are never reported separately. `shift` sums to zero over
the window's days with each day weighted equally, near-empty tails included: when the tails fill
up, `trend` rises and `shift` lowers the peak, and in birds the two cancel. Change from the first
season to the last (1993-2025; Chaffinch from 2007):

| Taxon          | Totals (trend + season + shift) | `trend` alone | `trend` with `shift` recentred on passage | Largest gap |
| -------------- | ------------------------------- | ------------- | ----------------------------------------- | ----------- |
| Black Kite     | x0.82                           | x2.95         | x0.83                                     | 6%          |
| Honey Buzzard  | x0.30                           | x1.67         | x0.30                                     | 2%          |
| Red Kite       | x9.6                            | x4.9          | x9.8                                      | 9%          |
| Common Buzzard | x0.79                           | x5.66         | x0.78                                     | 1%          |
| Osprey         | x0.83                           | x1.02         | x0.83                                     | 0%          |
| All pigeons    | x0.94                           | x2.98         | x0.92                                     | 2%          |
| Chaffinch      | x16                             | x9.5          | x16                                       | 2%          |

Abundance is therefore read from the smooth total in birds and timing from the median passage
date, which is what the export and defileViz show. Weighting `shift`'s constraint by the season's
passage, inside the fit, was considered and not done:

- It adds no output: the two quantities it would approximate are already exported exactly.
- It is not exact: a pure shift of δ days on a season of spread σ still moves the weighted `trend`
  by about δ²/2σ² (Red Kite, 12 days later since 1993, has the table's largest gap, 9%).
- It changes what the trend and shift penalties smooth, so the totals would need re-benchmarking.
- Its weights come from the fitted season (circular) or the raw counts (noisy tails).

If a split of abundance and timing is ever needed, for example for the forecast, recentre `shift`
on the passage after the fit (the table's fourth column) rather than refitting.

The benchmark (three refits per taxon, each hiding every target year; 40 s for seven taxa) splits
the error by kind of gap: hours hidden on days still partly counted are 15% of a year's birds and
3% of error, whole days hidden 6% and 2%. So the day's own count and its time-of-day profile matter
most, and the smooth terms least. What was tried on that split (`fe09ac4` has the code):

- `kappa` on blocks of 4 clock hours (`KAPPA_BLOCK`), kept: a day's passage shifts as a whole, so
  neighbouring hours are correlated and `kappa` halves from single hours to 3-4 hour blocks, then
  levels off. It fixed most of the under-coverage (Honey Buzzard 55% -> 70% of 80% intervals).
- Counted hours as flocks too (y negative binomial of size kappa x c in the day's rate update),
  dropped: 3% too low on gap-filled totals and no more accurate.
- Daily ERA5 weather at Défilé (wind, rain, cloud, temperature smooths), dropped for Explore: the
  error on whole days fell from 2.1% to 1.7% of a year but the total's did not (5.4% vs 5.3%), at
  2.5 x the time. It improves daily log scores (Black Kite -2.94 vs -3.02): a forecast matter.
- Finer episodes (`GAM_EPISODE_K` 35 and 50 instead of 25), kept at 25. At 25 knots the
  first-difference penalty sits at its lower bound for six of the seven taxa, so an episode's
  length is set by the knot spacing (fitted correlation halves after ~6 days for every taxon), not
  by the data. The data's own departures from the smooth model are shorter: on well-counted core
  days, with each year's mean removed, their autocorrelation halves from lag 1 to lag 2 and is gone
  by lag 3-4 (raptors, 1993-2025). Fitted on 2014-2025, more knots raise the marginal likelihood a
  little (+2 to +11 at 50 knots) and let the penalty act again, but even 50 knots (one every 2.6
  days) only bring the halving to ~4 days. The gap benchmark at 35 knots matched 25 (5.3% error,
  whole days 2.0% vs 2.1%, intervals 90% / 97%), with recent years slightly worse (0.62 vs 0.54).
  Fewer is worse: at 15 knots (episodes correlated for ~10 days) the error rose to 5.9% (Honey
  Buzzard 8.1%, Chaffinch 9.4%), intervals narrowed to 87% / 96%, and the marginal likelihood fell
  by 12-14, for no gain in time.
  On 1993-2025, now that the episode block is solved per year (below), 50 knots costs 2x the time
  of 25 and still does not help: the marginal likelihood rises by 0-25 per taxon, but the gap
  error is 5.4% vs 5.3% and recent years 0.60 vs 0.55. Variation faster than an episode (a day or
  two of weather) is left to the negative binomial's independent days, which the benchmark says is
  enough for totals.

Result: 5.3% gap-filling error, bias +0.5%, intervals slightly wide (91% / 97% for 80% / 95%),
Honey Buzzard still narrow (70% / 91%). Trends are exported for full-tier species and combined
series, not for unidentified birds ("falcon sp."), whose numbers follow identification effort.

**Fitting speed and convergence (2026-10).** A fit took 15-25 s and the full build ~4 min, and
the Fellner-Schall loop did not always converge:

- H is solved by blocks (`Factor`): the `episode` columns of different years never share a row,
  so the episode corner of H is block diagonal (33 blocks of 24) and the 891-column system reduces
  to the Schur complement on the other 99 columns. The results are unchanged (log marginal within
  0.06, totals within Monte Carlo noise; `tests/test_explore.py` checks the algebra against
  the dense matrix).
- The smoothing-parameter loop hit its 100-iteration cap for Black Kite: a fully penalised trend
  (`lambda` near its upper bound, a straight line) has a vanishing wiggle and room, and its update
  oscillated by +-0.1 forever. A penalty leaving less than 0.01 degrees of freedom now counts as
  converged, and repeated steps in the same direction are lengthened (Red Kite's episode penalty
  crept to its lower bound by 0.11 a step). Neither moves the fixed point.
- Worker processes ran one OpenBLAS thread per core each (load 110 on 12 cores): `worker_pool`
  gives each worker one thread.
- `doy` was `dayofyear`, so every window day of a leap year sat one day later in the season,
  `shift` and `episode` than the same date in other years: now the day of a non-leap year
  (`season_day`).

Together about 10x: a fit ~2-3 s, the benchmark 40 s with loading, the full build ~2.5 min, of
which ~80 s is the time-of-day profiles (`build_profiles`, serial).

Settings tried on the benchmark (GAM alone, seven taxa) beside the episode knots: none moved the
gap error off 5.3-5.4%, the benchmark's noise level, nor fixed Honey Buzzard's coverage. Season
12 or 30 knots, `shift` 4x6 or 8x12, trend 20 knots, a lower `lambda` bound (-10), observed instead
of expected information in the Laplace approximation (+20 log marginal over seven taxa, same
totals), and `MIN_COVERAGE` 0.05 or 0.2 (0.05 widens nothing and loses Osprey's coverage, 0.78)
all landed within +-0.002 of the base error. Dropping the year level or `shift` cost 111 and 56
of log marginal. Honey Buzzard's misses are whole hidden days at the late-August peak, estimated
~30% low; its day-level dispersion is not peakier than the season's (theta 1.2 in the tails, 0.9
at the peak), and with 11 target years an 80% coverage of 0.70 is within one standard error.
Across the 74 exported trends, `trend` is at its upper bound (a straight line on the log scale)
for 27 and `shift_year` (timing constant over the years) for 26: the year levels take the
multi-year waves, which the marginal likelihood prefers.

## Pipeline

**One function per taxon builds every block, from one fit (2026-10-09).** `pipeline.build_taxon`
takes a taxon's slice of a shared stage (release, effort, daily and hourly counts, profiles, age
and sex rows) and returns its species file: the raw tables, the trend GAM, and the blocks derived
from them (`season`, `daytime`, `age`, `sex`, `records`, `key_numbers`, `diagnostics`). Nothing
refits a model to draw a panel; a prototype that did (six species, every candidate panel) fixed
which blocks exist. The viewer (`scripts/explore_viewer.py`) only draws the export, so what it
shows is what defileViz would get.

**Settings by rule, exceptions by hand with a reason.** About 90 taxa get a full page and 74 a
trend, too many to configure one by one. Every setting has a rule (`settings.py`); a wrong value
is fixed by a better rule where possible, else by an entry in `overrides.yaml`, which refuses an
override without a `reason`. The species file records each value's source. Links follow the
same scheme: eBird, Birds of the World, EBBA2 and Trektellen by rule from the taxonomy's codes;
pages without a code to build them from (Vogelwarte, Migration Atlas) only by hand.

**Caches keyed by content, not by flags.** The shared stage is cached per release (file names,
sizes, times) and per source code of the modules that compute it; each trend fit per hash of the
taxon's data, profile, years and `trend.py`'s source. A code change invalidates exactly what it
touches, without remembering to clear anything. A rebuild after a change to a derived block takes
seconds (four taxa: 5 s instead of ~2 min), so a `--only <block>` option was not needed. A cold
build of all 273 taxa (74 trends) takes ~3 min on 12 cores: 109 s of shared stage, mostly the
time-of-day profiles, then 63 s for the taxa.

**Blocks name their method and version** (`season@1`, `trend-gam@1`, ...): defileViz's method
page will have one section per method, and a version bump says that what a block means changed.

**The season as counted: shares of each year's birds per day, gaps filled only for the total.** A
day's rate is birds per full counted day (`y / c`, not below `COVERAGE_MIN`); its share divides by
the year's total, where uncounted days are linearly interpolated between counted ones. The cells
shown are counted days only, each with its birds counted (`count`) for the hover. Each year's
10/50/90% passage dates come from the same filled series, with the share of the window counted, so
the viewer can fade years with long gaps. The
smooth dates (the GAM's season without the year's level and episodes, `trend.passage_quantiles`)
are drawn over them.

**Time of day: one histogram for the season, the date x hour matrix only where the hour changes**
(`daytime@2`, decided by the user: simpler, and a change shown only when there is one). Bars are
birds per counted hour as shares, pooled over every timed day (days with birds and at least
`PROFILE_MIN_TIMED` of them timed, so hours counted less often weigh less, not zero; an hour
counted under 3 h is left out), beside the smooth profile averaged over the main passage. The
change is the late part's mean passage hour minus the early part's (10-35% and 65-90% of the
passage), with a 95% interval from resampling days within each part, so one large flock cannot
make it look certain. The matrix is shown when the interval excludes zero and the shift is at
least 1 h: 19 of the 86 full-tier taxa, nearly all late taxa passing earlier on the clock
(Chaffinch -1.9 h, Red Kite -1.6 h), plus Great Cormorant, Barn Swallow and Black-headed Gull
passing later. Superseded by solar time (`daytime@3`, next entry). Measured in clock time on purpose: the end of summer time in late October and later
sunrises explain much of it, but a visitor reads a clock. A shift that is clear but under an hour
(Common Buzzard -0.7 h, Sparrowhawk -0.8 h) stays one histogram.

**Every hour of the day is solar time, not clock time (2026-10-09, decided by the user;
`DEFINITIONS_VERSION` 2, `daytime@3`).** Hour 12 starts at the sun's transit over Defile
(`export.solar_shift`: the zone's offset, less the longitude and the equation of time); days stay
local calendar dates. On the clock the transit moves from 13:43 (1 Aug) to 12:20 (3 Nov), with a
1 h step when summer time ends, in late October, at the peak of Wood Pigeon and Chaffinch passage.
The profile GAM, with 4 splines over the season, cannot draw that step. Comparing the mean passage
time of each timed day in four time bases (share of its spread explained by the season, lower is
steadier): Red Kite 0.18 by clock, 0.03 solar, 0.21 hours after sunrise, 0.03 share of daylight;
Common Buzzard 0.21 / 0.03 / 0.23 / 0.04; Sparrowhawk 0.12 / 0.03 / 0.14 / 0.02; Wood Pigeon
0.06 / 0.01 / 0.04 / 0.02; Chaffinch 0.13 / 0.09 / 0.10 / 0.08; only Great Cormorant is steadier
after sunrise. Solar time and the share of daylight tie; solar time keeps hours as units.

- Effort (`hourly`) is split into solar hours from the survey intervals; a count timed to an
  interval (an hour-by-hour sheet) is spread over the solar hours it overlaps, a point time falls
  in its hour, so counts and effort are binned alike. `hours` is the counted duration, night
  included. A survey of up to an hour now times its counts even across a clock hour.
- On the 2026-10 release, against a clock-time build of the same code: the gap benchmark's paired
  median change is -0.003 (mean 0.184 -> 0.170; 27 taxa better, 18 worse; Eurasian Curlew 1.17 ->
  0.75, Greylag Goose 0.57 -> 0.23, Alpine Swift 0.59 -> 0.32; Whimbrel 1.01 -> 1.18). Annual
  indices move 1.5% for the median taxon, up to 15-19% for Common Crane, Skylark, Greylag Goose,
  Blue Tit, Chaffinch. Totals hidden 7 -> 5, trends shown 38 -> 41; Yellow-legged Gull's season
  becomes hidden.
- `daytime` still shows 15 taxa by date, a different 15: the summer-time artefacts go (Wood Pigeon
  -1.5 -> -0.4 h, Merlin, Siskin, Hawfinch, Blue Tit, Greenfinch, Honey Buzzard), and morning
  migrants passing later in solar time as sunrise moves later appear (finch sp. +2.6 h, Redwing
  +1.7, Grey Wagtail +1.6, Hen Harrier +1.1). Best hours move about an hour earlier in their
  numbers (Red Kite 11-16 solar, was 12-17 clock).
- A page wanting the clock adds the day's `solar_shift`: about 1.6 h in summer time, 0.5 h in
  winter time.

**Age and sex: shares among the birds given one, only in usable years.** A year is usable with at
least 20 birds aged (or sexed), 5% of those counted, and both classes recorded (`demography@2`),
and a taxon needs 3 such years. The age codes `1`, `J`, `I`, `2` are one class, non-adult: the
codes used switch between years (Red Kite: mostly `1` in 2019, `I` in 2022). Every share has a 95%
Wilson interval; a pooled share over the usable years is given too.

**Both classes recorded: the smaller at least 10% of the birds aged or sexed.** Counters often tag
one class and leave the other blank. Red Kite came out 97% non-adult: in 2024, 4 030 non-adults
and 8 adults among 17 372 birds counted, in 2025 3 452 and none. Black Kite (2 334 non-adult, 0
adult in 2025), Grey Heron and gull sp. do the same; Honey Buzzard in 2024 the reverse (644 adults,
52 non-adults). The share among aged birds then measures what was tagged, not what passed. The
rule drops those years (Red Kite, Black Kite, Grey Heron and gull sp. lose the panel; 7 taxa keep
it). It also drops a year where one class is genuinely under 10%: such a taxon gets its years back
as an override, with the reason. Sex uses males against females and female-coloured birds.

**Sex is drawn like age** (`demography@3`, one `share_block` for both): the male share per usable
year, with its interval, sized by the birds sexed, the pooled share as a band, and when males and
female-types pass. The title gives the birds sexed (or aged) against all birds counted in those
years, since only part of them are.

**Records are all-time**, historical counts included (all pigeons: 116 340 on 20 Oct 1975),
decided by the user: a record is a count, not an estimate, so the start year does not apply.
Each top day carries what was written about it in the release (`records@2`, `remarks.py`): the
counters' survey remark of the day and the taxon's count remarks, including the report paragraphs
the dataset attached to that day (Red Kite, 14 Oct 2020: the site record of 2 112). Field labels
(`details:`) and the dataset's placeholders are dropped; the texts stay in French. Coverage is no
longer listed: a record is what was counted. The best day of each year is gone too.

**Every trend is benchmarked in the build, and each of its claims classed (`benchmark@1`,
`reliability@1`).** The gap-transplant test that judged the model on seven taxa now runs on all 74
trend taxa, on the same years, window and profile as the export (`benchmark.py`, cached with the
trend; ~1 min more on a full build). Its median error is 0.10 (0.07-0.20 between quartiles); seven
taxa miss by more than 0.4, all flocking or rare species whose hidden days hold a few big flocks
(Purple Heron 1.8, Eurasian Curlew 1.4, Whimbrel 1.2, Mediterranean Gull 1.1, Alpine Swift, Greylag
Goose, Northern Lapwing ~0.6), always too low. Grey Heron has no well-counted target year and is
untested. Three claims are classed, decided by the user on the proposal:

- `totals`: hidden above 0.4 of gap error; a caveat above 0.15, with 80% intervals holding the truth
  under 60% of 9+ trials, under 60% of a typical year counted, q90/q10 above 2, or untested (with a
  borrowed profile as a second reason: the benchmark measures a borrowed profile's error).
- `trend`: the totals' class, and the smooth's 95% band (median q97.5/q2.5 over the years): a
  caveat above 4, hidden above 20.
- `season`: the median passage date's 80% band: a caveat above 10 days, hidden above 30; a caveat
  when the passage runs past the counted days.

With the posterior's information floor (next entry): totals 48 show / 19 caveat / 7 hide; trend
37 / 30 / 7 (the seven whose totals are hidden); season 27 / 47 / 0. Before it, the blown-up bands
hid 22 trends and 4 seasons.
The recent-seasons test is exported and drawn but not used: a season's level is not predictable
from earlier seasons (median error 0.8), which is a forecast question. The year-by-year marking
(`estimated_years`, under half counted) is for the page to show, not a class.

**The classes are applied to what a page draws in the export, not in defileViz
(`reliability@2`).** `reliability.ELEMENTS` names each field that rests on a claim (the season
totals and the typical-season key number on `totals`; the smooth and the trend key number on
`trend`; the modelled passage dates, their smooth and the main-passage key number on `season`) and
the block's `elements` gives each its class, so defileViz looks a field up rather than knowing
which claim it belongs to. Counts (birds counted, the empirical phenology, chances, daytime, age,
sex, records) rest on no claim and are always shown. The viewer draws everything and flags what
defileViz will do: a summary at the top of each page, a strip over each figure, hidden traces
greyed and tagged `[HIDE]` with a watermark, the mostly estimated years ringed. On the 2025
release: totals 49 / 18 / 7, trend 38 / 29 / 7, season 27 / 47 / 0.

**A count with its own interval of up to an hour is timed at its midpoint
(`release.MAX_COUNT_INTERVAL`); a longer one is kept at day level.** Trektellen's hourly entries
reach the release from 2025 (30 counts on 23 Nov 2025, each 1 h); refusing them stopped the build,
and the hour they fall in is all the pipeline uses. Decided by the user: no need for more.

**The posterior's information is evaluated at a rate of at least 0.1 birds per full day
(`trend.MIN_RATE`, `information_weights`).** Where a taxon is counted but never seen (Common Wood
Pigeon on 18-25 July: 33 years of zeros), the season spline runs towards minus infinity and the
negative binomial's information, proportional to the expected count, with it. The log-likelihood
is flat below the mode but steep above it (zeros were counted), and the Laplace posterior, a
Gaussian with the curvature at the mode, was flat both ways: log-scale sd 12-16 on those days
against 0.4-0.5 in the passage, so a few draws held millions of birds. The smooth's q97.5 was 13
million for Wood Pigeon in 2009 against a median of 29 000 (now 47 000), its median passage
date's 80% band up to 100 days wide (Song Thrush, Redwing; now at most 29). Taking each day's
information as at least that of 0.1 birds per full day bounds such days by a rate the zeros still
allow, and leaves days with more birds, and the posterior mode, unchanged. It applies to the
posterior and the Laplace marginal likelihood, not to the IRLS steps, which it slowed. Across the
74 trend taxa: gap-filled totals within 1% (Water Pipit -8%), benchmark error 0.101 -> 0.097,
80% / 95% interval cover 0.81 / 0.94 -> 0.79 / 0.93, smooth band median 3.6 -> 2.7 (max 5e8 ->
14), same build time. `tests/test_explore.py` has a synthetic case (24x without the floor).
Rejected: summing only days with expected birds (hides the problem in one output), a prior on the
spline's edge (would need tuning per taxon).

**Windows per taxon: fitted on a model window, shown over a wider view window, by rule**
(`window@1`). The default window (18 Jul - 18 Nov) is where 80-100% of years are counted; from
20 Nov to 1 Dec only about a third are. Late taxa still pass then: Red Kite averages around 120
birds a day on 30 Nov in the years counted, and 2025 had 2 600 counted on 19-23 Nov. So each taxon
gets:

- a model window: the default, extended to 0.5% / 99.5% of its mean passage plus 5 days, but only
  over days counted in at least a third of its years, since the trend fills the rest. Never cut
  below the default, so a taxon passing inside it keeps the benchmarked fit. 41 of the 74 trend taxa
  are extended, by at most 6 days (to 20-24 Nov; Black Kite and Sand Martin to 14-15 Jul). Red
  Kite's 2025 total goes from 20 700 to 24 700.
- a view window for the season panels, extended the same way to 0.1% / 99.9% over days counted in
  at least 10% of years (to 2 Dec for most late taxa). Totals, shares and the yearly passage dates
  stay those of the model window: beyond it, the gaps to fill are too long, and filling them held
  the last counted day flat to 2 Dec.
- `passage_beyond_counting` (51 full-tier taxa): the passage still reaches the model window's
  outer limit, so the totals are totals up to that date. A property of the count, for the page to
  state, not a fault to fix.

Both windows are settings, overridable as `[MM-DD, MM-DD]` with a reason.

**Key numbers** are read from the blocks, not computed apart: the main passage is the smooth
season's 10-90% in the last year (else the pooled counted seasons of the last 10), best hours the
fewest solar hours holding 60% of the profile on the median passage day, the typical season the
median of the last 10 gap-filled totals, the trend the smooth's change from the start year, the
chance the share of well-counted days in the main passage with at least 1 or 10 birds.

**Written accounts replace the report extracts** (`accounts@1`). The page showed the paragraphs
of the annual reports as extracted (French, unedited, and for a combined series those of its first
member). It now shows authored accounts from `content/accounts/`: a general account in up to three
sections and a short account per season, in French and English, edited from the reports and the
2019 paper with the numbers left to the page's figures. The build reads only the two authored TSVs
and checks them; `report_text.csv` and `paper_text.csv` are still copied with the release, for the
editors' coverage and totals checks. `reports.json` (the non-species report texts) is no longer
written: defileViz never read it. Location: `content/accounts/` in this repo, versioned and next to
the code that reads it, rather than in defile-dataset (the accounts are written for this page, not
data) or `data/` (generated, never committed); the authored files at its top, `reading/` and
`review/` for what the editorial scripts regenerate.
