# Hourly weather pilot

This experiment splits complete observed Trektellen surveys independently of bird entries, fits raw hourly counts with real-duration exposure, and tests whether weather improves prediction of contiguous withheld intervals. It reads the current defile-dataset release and makes a bounded snapshot of the local ERA5 cache. It imports no forecast code and does not change Explore's production models or exports.

Run from the defile-explore root:

```bash
uv run --with pyarrow==13.0.0 python scripts/research/build_hourly_weather.py --dataset ../defile-dataset/output/dataset --weather-cache ../defile-migration-forecast/data/weather/era5_hourly
Rscript scripts/research/analyse_hourly_weather.R
```

The builder defaults to completed seasons 2022–2025 and four pilot taxa: Red Kite, Common Buzzard, European Honey Buzzard and Common Wood Pigeon. Override `--years`, `--taxa` and `--out` to extend the pilot. R requires mgcv, dplyr, tidyr, readr, ggplot2, here, cli and MASS; run the script section by section to inspect its intermediate objects. Optional R arguments override the input and output folders, in that order.

## Observation accounting

- `hourly_effort.csv` contains each complete observed survey cut at UTC clock-hour boundaries, including empty and partial hours. Calendar-day weather closures are not observed effort.
- `hourly_counts.csv` contains main-direction counts with supported hourly placement, for all taxa in the selected seasons. Exact timestamps are assigned to their hour. Own or inherited intervals are assigned only when wholly contained within one clock hour and their survey.
- `unresolved_counts.csv` retains daily/untimed counts, intervals crossing hour boundaries and timestamps outside the released survey. Entries are never moved, redistributed or silently discarded.
- `model_frame.csv` joins effort to the selected species, including species-specific zeros. An unresolved positive or presence-only record makes all hours of that species' survey unsuitable for hourly fitting. Their `count` is missing, while `timed_count` retains the supported lower bound. An untimed explicit zero does not invalidate the survey.
- `weather_stops.csv` keeps documented weather closures separately. Their zeros are assumptions and are not training or validation observations.
- `timing_reconciliation.csv`, `species_summary.csv` and `metadata.json` record conservation of birds and effort, exclusions, units, source paths and SHA-256 hashes.

Survey completeness is used under the dataset's assumption of systematic species recording. Exact survey edges are retained. The first fit uses only intervals wholly within civil dawn/dusk (sun altitude −6°); boundary and night intervals remain in the prepared data and their exclusions are reported in `eligibility.csv`. It does not interpret missing dates or unobserved hours as zero.

## Weather and time

The predictor `tau` is the fraction of that local day's civil daylight elapsed at the observation interval's midpoint. Duration remains in real hours. The GAM is an hourly approximation: it evaluates the rate at the interval midpoint rather than integrating its varying diurnal curve within the hour. The time-coordinate transformation does not normalize predicted counts or force equal daily totals.

The weather snapshot uses the forecast cache's existing converted units. Precipitation in metres is converted to millimetres; accumulated solar radiation in J/m² is converted to W/m². [Open-Meteo defines rain and radiation at timestamp t over the preceding hour](https://open-meteo.com/en/docs/historical-weather-api), so the builder aligns the value at t+1 hour with count interval [t, t+1 hour). The previous-24-hour rain summary excludes the current interval. Wind components and low cloud are approximated by the mean of the instantaneous values at the two endpoints. Lags are computed on a complete hourly weather index before selecting observation times, retaining missingness across any cache gaps.

## First model and comparison

Both models are negative-binomial GAMs with `offset(log(exposure_hours))`, separate year levels, a seasonal curve, a noncyclic daylight-time curve, a season × daylight-time interaction and an independent random effect for each day. The weather model adds current precipitation, preceding-24-hour precipitation, wind u/v, solar radiation and low cloud. Precipitation predictors use log1p. Smoothing and dispersion are estimated separately within each training fold using mgcv::bam. Neither model includes the current episode smooth. With four recent seasons, this pilot estimates year levels rather than a long-term trend or gradual phenology shift.

Three folds withhold seven-day calendar blocks, balanced across years, and test each eligible interval once. A fourth fold hides up to three adjacent intervals around the wettest observed hour, or the central hour on dry days, retaining other hours on the same day. Masks depend on dates, weather and effort, never bird counts. A partly observed day is conditioned on its remaining hours through its random effect. Wholly withheld days receive draws from the learned distribution of day effects, rather than their full-data fitted effects. Daily totals sum the same posterior predictive draws, preserving shared-day and coefficient uncertainty. Intervals condition on estimated smoothing and dispersion parameters; they do not account for weather measurement error.

`logs/hourly-weather/` contains fitted models, hourly and daily held-out predictions, wet/dry scores, dispersion and residual diagnostics, partial precipitation-effect curves, figures and the R session information. Scores concern the observed intervals that were hidden, not an independently known dawn-to-dusk count. The rain effect is conditional on the other predictors and is not a causal estimate. Diagnostics include basis checks and concurvity because radiation, rain and cloud may overlap strongly.

This is a pilot for deciding whether to develop the interval-integrated model. It does not yet fit historical daily totals, fill production gaps, estimate residual correlation across days, establish zero passage in unobserved severe weather, or provide a count distribution invariant to arbitrary temporal aggregation. Negative-binomial overdispersion is per observation; posterior predictive daily sums share day effects but retain conditional hourly independence. Review blocked-validation bias and interval coverage before interpreting any model as an improvement.

## First run (10 October 2026)

The 2022–2025 release produces 6,421 observed hourly/partial-hour intervals, totalling 5,953.82 hours. All 2,022,996 main-direction birds reconcile between supported timed counts and retained unresolved records; no weather data are missing for the selected years. Each pilot species has approximately 1,500 observed wet intervals. All 85 Python tests pass; saved-model checks verify real-hour exposure scaling, conservation of withheld counts and independently recomputed summary metrics.

Weather improves the mean held-out wet-hour log predictive density for all four pilot species. On wholly withheld wet days, Common Buzzard mean absolute error drops from 88.6 to 42.7 birds and Red Kite from 111.3 to 58.0 birds. Their aggregate wet-period bias drops from +137% to +14% and from +114% to +13%, respectively. These are exploratory differences from the same completed validation masks, not significance tests or a comparison against the production episode model.

The Gaussian coefficient draws are unstable in zero-count seasonal tails for Common Wood Pigeon and European Honey Buzzard: a few extreme draws dominate the predictive mean even when upper count quantiles remain near zero. Their mean-based gap totals are unusable. The existing trend engine already addresses this low-information problem (`trend.MIN_RATE`, `information_weights`); applying an appropriate treatment to this mgcv pilot is the next uncertainty task. Until that is resolved and the pilot is compared against the current production model, retain this experiment for inspection rather than using it to replace gap-filled totals.

The first execution saved all fits and outputs but encountered a trailing parser error after a source comment was edited while R was running. The source then passed a complete syntax check, and summaries and plots were regenerated from the completed validation predictions; `logs/hourly-weather/verification.log` records successful model, exposure, conservation and score checks.
