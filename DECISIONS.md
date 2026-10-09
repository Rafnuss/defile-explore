# Défilé Explore — decisions log

A running log of settled calls for the Explore statistics: what was tried, what worked, what was
rejected and why. Not a design spec; read the code for that. See `DEVELOPMENT.md` for what is
still open.

## Repository

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
implementation) on local clock hours, from days timed to the hour; a counted day's coverage `c` is
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
