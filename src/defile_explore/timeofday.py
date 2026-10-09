"""The time-of-day GAM: a smooth surface of (an hour's rate / its day's rate) over (doy, hour).

A copy, not a shared module: taken from defile-migration-forecast `src/phenology.py` (commit
13c5107, 2026-10-09), where the same fit gives the forecast's phenology prior. The two repos are
deliberately independent: a change to one copy is carried to the other by hand, as a decision
(DECISIONS.md -> Repository).
"""

import numpy as np

# GAM spline counts for the doy/hour ratio surface.
DOY_SPLINES = 4
HOUR_SPLINES = 12

# Regularization ladder for the PIRLS fit, weakest first. `ratio` is heavily zero-inflated and
# long-tailed, and PoissonGAM's PIRLS diverges for some taxa at (100, 10): progressively stronger
# regularization is tried and the first that converges kept.
GAM_LAM_LADDER = [(100, 10), (1000, 100), (10_000, 1_000), (100_000, 10_000)]

# Each day's hourly-ratio samples are weighted by (that day's bird count) ** this. 0.5 was best on
# held-out years in the forecast's test (fit on even years, L1 to odd years' count-weighted hourly
# profile, 7 species): 0.253 / 0.213 / 0.229 for 0 / 0.5 / 1.
RATIO_WEIGHT_POWER = 0.5


def fit_ratio_surface(
    doy: np.ndarray,
    hour: np.ndarray,
    ratio: np.ndarray,
    weights: np.ndarray,
    doy_grid: np.ndarray,
    hours: np.ndarray,
    k0: int = DOY_SPLINES,
    k1: int = HOUR_SPLINES,
    interaction: bool = True,
    lam_ladder=GAM_LAM_LADDER,
    label: str = "",
) -> np.ndarray:
    """Poisson GAM of `ratio` (a period's rate / its day's rate) on (doy, hour), predicted on
    `doy_grid` x `hours`: shape `(len(doy_grid), len(hours))`.

    `interaction` fits `te(doy, hour)`, else `s(doy) + s(hour)`. If PIRLS diverges, the fit is
    retried with progressively stronger regularization from `lam_ladder`, and so is one whose
    prediction overflows.
    """
    from pygam import PoissonGAM, s, te
    from pygam.utils import OptimizationError

    term = te(0, 1, n_splines=[k0, k1]) if interaction else s(0, n_splines=k0) + s(1, n_splines=k1)
    X = np.column_stack([doy, hour])
    grid = np.column_stack([np.repeat(doy_grid, len(hours)), np.tile(hours, len(doy_grid))])
    for lam in lam_ladder:
        try:
            gam = PoissonGAM(term, lam=list(lam)).fit(X, ratio, weights=weights)
        except OptimizationError:
            continue
        surface = gam.predict(grid)
        if np.isfinite(surface).all():  # PIRLS can also overflow without raising
            break
    else:
        raise OptimizationError(f"{label}: PIRLS did not converge even at lam={lam_ladder[-1]}")
    if lam != lam_ladder[0]:
        print(f"  {label}: PIRLS needed lam={lam} to converge")
    return surface.reshape(len(doy_grid), len(hours))
