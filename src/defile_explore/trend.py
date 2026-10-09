"""Trend model for the Explore page: one taxon's smooth trend, season and phenology shift, and
gap-filled annual totals. Part of `defile_explore`.

A GAM, benchmarked by `scripts/benchmark_trend.py` (`DECISIONS.md` -> Explore has why a GAM):

    count_day ~ NegBin(mean = c_day * exp(eta), theta)
    eta = a + trend(year) + year_effect(year) + season(doy) + shift(year, doy)
          + episode_year(doy)

`c_day` is the day's coverage (`defile_explore.profile.coverage`: the share of the day's expected
passage in the hours counted), `year_effect` an independent level per year (good and bad years),
`shift` a smooth change of the season's timing or width over the years, and `episode` a short-range
curve of its own each year: runs of good or bad migration days, the weather's share, which a
weather model would one day explain. Without it the shift absorbed them (lengthscales of 3 years by
6 days), and filled days were treated as independent of their neighbours. The negative binomial is a
day's gamma-distributed rate seen through a Poisson count; within the day, passage comes in flocks,
so the hours not counted are drawn with an hourly over-dispersion `kappa` (`hour_dispersion`,
`fill_draws`).

Every component is a basis times coefficients with a Gaussian prior (a penalty): P-splines, cubic
B-spline bases with second-difference penalties (a tensor product for `shift`, first differences
plus a ridge for `episode`), a smoothing parameter per penalty, as mgcv would. The posterior mode is
found by penalised IRLS, and the smoothing parameters and `theta` by the Laplace approximation of
the marginal likelihood (`fit`). Variants:

- `gam`: the model above.
- `season`: the season alone, a fixed phenology without trend or year effect, as the forecast's
  baseline has: the reference the GAM must beat.

`trend` alone is not a change in abundance, nor `shift` alone a change in timing. The constraints
hold on the log scale with every day of the window weighted equally, near-empty tails included:
tails filling up raise `trend`, and `shift` lowers the peak to keep its mean at zero, so the two
cancel in birds (Common Buzzard 1993-2025: `trend` x5.7, totals x0.8). Read abundance from the
smooth total in birds (`annual_totals`, which keeps `shift`) and timing from the median passage
date (`peak_doy`). Recentring `shift` on the passage-weighted mean after the fit gives a `trend` within
0-9% of the totals' change, but a pure shift of the season still leaks into it, and moving that
weighting into the fit would change what the penalties smooth (`DECISIONS.md` -> Explore).

Beyond the data (a year not yet counted), the trend extrapolates its last slope: Explore never asks
for it, and a forecast that did would need a flat extrapolation instead.

Uncertainty is the posterior given the fitted hyperparameters (empirical Bayes); intervals for annual
totals are calibrated by the benchmark's gap-transplant test.
"""

from __future__ import annotations

import multiprocessing
import os
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from scipy import optimize, sparse
from scipy.interpolate import BSpline
from scipy.linalg import cho_factor, cho_solve, solve_triangular
from scipy.special import gammaln

from defile_explore.export import HOURS, by_hour, in_window
from defile_explore.profile import PROFILE_DOY, PROFILE_MIN_TIMED, coverage
from defile_explore.season import chance_from_draws, passage_from_draws

VARIANTS = ("season", "gam")
# Ranks whose trend is exported: a species, or a group read as everything below it (`members`).
# A group with nothing below it ("Acrocephalus sp.") trends with how hard observers try to
# identify it, not with the birds.
TREND_RANKS = ("species",)

# Prior precision on coefficients no penalty reaches (intercept, spline null spaces): effectively
# flat on the log scale (sd 100), and what keeps the overlapping constants identifiable.
RIDGE = 1e-4
ETA_MAX = 20.0  # log birds per full-day equivalent; guards exp() in early IRLS steps

# P-spline bases (cubic, second-order difference penalties).
GAM_TREND_K = 12
GAM_SEASON_K = 20
GAM_SHIFT_K = (6, 8)  # year x doy
# Per year, over the window's ~124 days: a knot every ~5.6 days, which (not the data) sets how long
# an episode lasts; 15 is worse, 35 and 50 no better on totals (DECISIONS.md -> Explore).
GAM_EPISODE_K = 25
GAM_LOG_LAMBDA_BOUNDS = (-6.0, 16.0)
LOG_LAMBDA_START = 2.0
EFS_MAX_ITER = 100
EFS_MAX_GAIN = 8.0  # longest step, in multiples of the plain update (`_efs`)
THETA_ROUNDS = 5  # alternations of theta (marginal likelihood) and the smoothing parameters
EFS_TOL = 1e-2  # largest change of a log smoothing parameter or log theta: 1%

THETA_BOUNDS = (0.02, 200.0)  # negative binomial shape; small = over-dispersed
KAPPA_BOUNDS = (0.01, 1e4)  # hourly over-dispersion, per full day of coverage; large = Poisson
# Hours summed into blocks before `kappa` is fitted: a day's passage shifts as a whole (later on a
# slow morning), so neighbouring hours are not independent and `kappa` falls from single hours to
# blocks of 3-4 hours, then levels off (Honey Buzzard 4.7 -> 2.6). Gaps are blocks of hours too.
KAPPA_BLOCK = 4
SMOOTH = ("year", "episode")  # left out of the smooth trend, season curves and passage dates

# Below this coverage a day says nothing about its rate (75 Black Kites in 11 minutes on
# 2025-07-27, c = 0.0004, would be a 190 000-bird day): it is neither fitted nor conditioned on,
# and its birds are added to the filled total as they are (`extra`).
MIN_COVERAGE = 0.1
# birds per full day at which a day's information is evaluated, at least (`information_weights`)
MIN_RATE = 0.1

IRLS_MAX_ITER = 100
IRLS_TOL = 1e-8
DRAWS = 1000
QUANTILES = (0.025, 0.1, 0.5, 0.9, 0.975)
FILLED_DAY_INTERVAL = (0.1, 0.9)  # a filled day: its mean and this 80% interval
PASSAGE_QUANTILES = (0.1, 0.5, 0.9)  # passage dates: 10%, half and 90% of the season passed

# One BLAS thread per worker (`worker_pool`): numpy's OpenBLAS otherwise starts one per core in every
# worker, and a pool of one worker per core ran ~12x oversubscribed (load 110 on 12 cores).
WORKER_THREAD_VARS = ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS")


def worker_pool(workers: int) -> ProcessPoolExecutor:
    """A process pool for fitting taxa in parallel, each worker on one BLAS thread: the workers are
    spawned (not forked, which would inherit the parent's threads) with `WORKER_THREAD_VARS`
    set."""
    for var in WORKER_THREAD_VARS:
        os.environ[var] = "1"
    return ProcessPoolExecutor(workers, mp_context=multiprocessing.get_context("spawn"))


# --- data -------------------------------------------------------------------


def season_day(dates: pd.DatetimeIndex) -> np.ndarray:
    """Day of year in a non-leap year, so that a calendar date has the same `doy` every year: with
    `dayofyear`, every window day of a leap year was one day later in the season, `shift` and
    `episode` than the same date in other years."""
    return dates.dayofyear.to_numpy() - (dates.is_leap_year & (dates.month > 2))


def model_frame(
    days: pd.DataFrame,
    effort: pd.DataFrame,
    c: pd.Series,
    first_year: int,
    last_year: int,
    window: tuple[int, int] | None = None,
) -> pd.DataFrame:
    """One row per window day of `first_year`..`last_year`: `date`, `year`, `doy` (`season_day`),
    `c`, `y` and `extra`. `window`: first and last season day (`defile_explore.window`); the
    default window `WINDOW` if None.

    Counted days (`state == "counted"`) carry their coverage `c` and count `y` (0 without a record);
    every other day, and a counted day with only a presence record, has `c = 0` and `y = 0`: nothing
    observed, all of it to fill. A counted day below `MIN_COVERAGE` is treated as not counted, its
    birds kept in `extra`.
    """
    dates = pd.date_range(f"{first_year}-01-01", f"{last_year}-12-31", freq="D")
    if window is None:
        dates = dates[in_window(pd.Series(dates)).to_numpy()]
    else:
        sd = season_day(dates)
        dates = dates[(sd >= window[0]) & (sd <= window[1])]
    f = pd.DataFrame({"date": dates, "year": dates.year, "doy": season_day(dates)})
    counted = effort["state"] == "counted"
    f["c"] = f["date"].map(c[counted].set_axis(effort.loc[counted, "date"])).fillna(0.0)
    f["y"] = f["date"].map(days.set_index("date")["count"])
    no_count = (f["c"] > 0) & f["y"].isna() & f["date"].isin(days["date"])  # presence only
    f.loc[no_count, "c"] = 0.0
    f["y"] = f["y"].where(f["c"] > 0, 0.0).fillna(0.0)
    low = f["c"] < MIN_COVERAGE
    f["extra"] = f["y"].where(low, 0.0)
    f.loc[low, ["c", "y"]] = 0.0
    return f


# --- bases and priors -------------------------------------------------------


def bspline_basis(x: np.ndarray, lo: float, hi: float, k: int) -> np.ndarray:
    """Cubic B-spline basis with `k` functions on equally spaced knots over [lo, hi]."""
    inner = np.linspace(lo, hi, k - 2)
    step = inner[1] - inner[0]
    knots = np.r_[lo - 3 * step + step * np.arange(3), inner, hi + step * np.arange(1, 4)]
    return BSpline.design_matrix(np.clip(x, lo, hi), knots, 3).toarray()


def difference_penalty(k: int, order: int = 2) -> np.ndarray:
    d = np.diff(np.eye(k), order, axis=0)
    return d.T @ d


@dataclass
class Block:
    """A model component: its columns in the design and its penalties, each a smoothing parameter's
    name and the matrix it multiplies (the prior precision is their sum).

    With `groups` > 1 the columns are that many equal groups (`episode`: one per year), each with
    the same penalties.
    """

    name: str
    columns: slice
    penalties: list[tuple[str, np.ndarray]]
    groups: int = 1

    @property
    def size(self) -> int:
        """Columns of one group."""
        return (self.columns.stop - self.columns.start) // self.groups


@dataclass
class Design:
    variant: str
    years: np.ndarray  # every year of the model's range, observed or not
    doy_range: tuple[int, int]
    blocks: list[Block]
    _basis: callable = field(repr=False)

    def X(self, year: np.ndarray, doy: np.ndarray) -> sparse.csr_matrix:
        return self._basis(np.asarray(year), np.asarray(doy))

    @property
    def n_columns(self) -> int:
        return self.blocks[-1].columns.stop

    @property
    def grouped(self) -> Block | None:
        """The last block when it is grouped: its columns are kept per group (`Matrix`)."""
        return self.blocks[-1] if self.blocks[-1].groups > 1 else None

    @property
    def hyper_names(self) -> list[str]:
        return [name for b in self.blocks for name, _ in b.penalties]

    def block_precision(self, b: Block, h: dict) -> np.ndarray:
        """Prior precision of one group of the block's columns."""
        return sum((np.exp(h[n]) * S for n, S in b.penalties), RIDGE * np.eye(b.size))

    def prior(self, h: dict) -> Prior:
        g = self.grouped
        split = g.columns.start if g else self.n_columns
        Q = np.zeros((split, split))
        for b in self.blocks:
            if b is not g:
                Q[b.columns, b.columns] = self.block_precision(b, h)
        return Prior(
            Q, self.block_precision(g, h) if g else np.zeros((0, 0)), g.groups if g else 0
        )

    def columns_of(self, *names: str) -> np.ndarray:
        """Indices of the named blocks' columns."""
        return np.concatenate(
            [np.arange(b.columns.start, b.columns.stop) for b in self.blocks if b.name in names]
            + [np.zeros(0, int)]
        )


def design(variant: str, years: np.ndarray, doy_range: tuple[int, int]) -> Design:
    """The bases and penalties of a variant over `years` (all of them, including years without
    data, whose coefficients then stay at their prior) and the window's days of year.

    Every smooth sums to zero over its grid (`sum_to_zero`), so each level has one owner: the
    intercept, then `year` for a year's level (`episode` has zero mean over each season), and the
    trend and the season for the main effects of `shift` (an interaction only, as mgcv's `ti`).
    """
    y0, y1 = int(years.min()), int(years.max())
    d0, d1 = doy_range
    n_years = y1 - y0 + 1
    year_grid, doy_grid = np.arange(y0, y1 + 1), np.arange(d0, d1 + 1)
    parts, blocks = [], []
    col = 0

    def add(name, basis, penalties, groups=1):
        nonlocal col
        k = basis(np.array([y0]), np.array([d0])).shape[1]
        parts.append(basis)
        blocks.append(Block(name, slice(col, col + k), penalties, groups))
        col += k

    def centred(lo, hi, k, grid, order=2):
        """A B-spline basis that sums to zero over `grid`, and its difference penalty."""
        Z = sum_to_zero(bspline_basis(grid, lo, hi, k))
        return (lambda v: bspline_basis(v, lo, hi, k) @ Z), Z.T @ difference_penalty(k, order) @ Z

    add("intercept", lambda y, d: np.ones((len(y), 1)), [])
    season, S = centred(d0, d1, GAM_SEASON_K, doy_grid)
    add("season", lambda y, d: season(d), [("season", S)])
    if variant == "gam":
        trend, St = centred(y0, y1, GAM_TREND_K, year_grid)
        add("trend", lambda y, d: trend(y), [("trend", St)])
        ky, kd = GAM_SHIFT_K
        shift_y, Sy = centred(y0, y1, ky, year_grid)
        shift_d, Sd = centred(d0, d1, kd, doy_grid)
        add(
            "shift",
            lambda y, d: _row_kron(shift_y(y), shift_d(d)),
            [
                ("shift_year", np.kron(Sy, np.eye(kd - 1))),
                ("shift_doy", np.kron(np.eye(ky - 1), Sd)),
            ],
        )
        add(
            "year",
            lambda y, d: np.eye(n_years)[np.clip(y - y0, 0, n_years - 1)],
            [("year", np.eye(n_years))],
        )
        episode, Se = centred(d0, d1, GAM_EPISODE_K, doy_grid, order=1)
        add(
            "episode",
            lambda y, d: _per_year(y - y0, n_years, episode(d)),
            [("episode", Se), ("episode_size", np.eye(GAM_EPISODE_K - 1))],
            groups=n_years,
        )

    def basis(y, d):
        return sparse.hstack([sparse.csr_matrix(p(y, d)) for p in parts], format="csr")

    return Design(variant, year_grid, doy_range, blocks, basis)


def sum_to_zero(B: np.ndarray) -> np.ndarray:
    """`Z` (k, k - 1) such that `B @ Z @ b` sums to zero over the rows of `B` for any `b`: the
    centring constraint that separates a smooth from the intercept and from other smooths."""
    q, _ = np.linalg.qr(B.sum(axis=0)[:, None], mode="complete")
    return q[:, 1:]


def _per_year(index: np.ndarray, n_years: int, b: np.ndarray) -> sparse.csr_matrix:
    """`b`'s row placed in the columns of its year: a separate curve per year."""
    n, k = b.shape
    index = np.clip(index, 0, n_years - 1)
    cols = (index[:, None] * k + np.arange(k)).ravel()
    rows = np.repeat(np.arange(n), k)
    return sparse.csr_matrix((b.ravel(), (rows, cols)), shape=(n, n_years * k))


def _row_kron(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Row-wise Kronecker product: the tensor-product basis."""
    return (a[:, :, None] * b[:, None, :]).reshape(len(a), -1)


# --- fitting ----------------------------------------------------------------


def nb_loglik(y: np.ndarray, mu: np.ndarray, theta: float) -> np.ndarray:
    return (
        gammaln(y + theta)
        - gammaln(theta)
        - gammaln(y + 1)
        + theta * np.log(theta / (theta + mu))
        + y * np.log(mu / (theta + mu) + 1e-300)
    )


class Matrix:
    """The design at the observed rows, kept as its dense columns and, for a grouped last block
    (`episode`, one group per year), each row's own group and that group's `k` values: a row of
    `episode` is zero outside its year.

    X'WX is then a dense part, a cross term and a block diagonal (`information`), which `Factor`
    solves by blocks.
    """

    def __init__(self, X: sparse.csr_matrix, split: int, groups: int):
        n = X.shape[0]
        self.split, self.groups = split, groups
        self.dense = X[:, :split].toarray()
        if groups:
            Xg = X[:, split:].toarray().reshape(n, groups, -1)
            self.group = np.abs(Xg).sum(axis=2).argmax(axis=1)
            self.values = Xg[np.arange(n), self.group]  # (n, k)
            self.rows = [np.flatnonzero(self.group == g) for g in range(groups)]

    def __matmul__(self, b: np.ndarray) -> np.ndarray:
        out = self.dense @ b[: self.split]
        if self.groups:
            bg = b[self.split :].reshape(self.groups, -1)
            out += (self.values * bg[self.group]).sum(axis=1)
        return out

    def rmatvec(self, r: np.ndarray) -> np.ndarray:
        top = self.dense.T @ r
        if not self.groups:
            return top
        grouped = [self.values[i].T @ r[i] for i in self.rows]
        return np.concatenate([top, np.concatenate(grouped)])

    def information(self, w: np.ndarray) -> tuple:
        """X'WX as `(A, C, D)`: the dense block, the cross terms per group `(G, k, p)` and the
        diagonal blocks per group `(G, k, k)`."""
        dw = self.dense * w[:, None]
        A = self.dense.T @ dw
        if not self.groups:
            return A, np.zeros((0, 0, len(A))), np.zeros((0, 0, 0))
        C = np.stack([self.values[i].T @ dw[i] for i in self.rows])
        D = np.stack([self.values[i].T @ (self.values[i] * w[i, None]) for i in self.rows])
        return A, C, D


@dataclass
class Prior:
    """Prior precision: `dense` on the first columns, then `group` repeated `groups` times."""

    dense: np.ndarray
    group: np.ndarray
    groups: int

    def __matmul__(self, b: np.ndarray) -> np.ndarray:
        p = len(self.dense)
        top = self.dense @ b[:p]
        if not self.groups:
            return top
        return np.concatenate([top, (b[p:].reshape(self.groups, -1) @ self.group).ravel()])

    def logdet(self) -> float:
        return _logdet(self.dense) + (self.groups * _logdet(self.group) if self.groups else 0.0)


class Factor:
    """H = [[A, C'], [C, blockdiag(D)]] (A the dense columns, D one block per group) by its blocks:
    each D's Cholesky factor and the Schur complement S = A - C' D^-1 C. Solves, the
    log-determinant, the diagonal blocks of H^-1 and draws then cost O(G k^3 + p^3) instead of
    O((p + G k)^3): for 1993-2025 with 25 episode knots, 33 blocks of 24 and 99 columns against 891.
    """

    def __init__(self, A: np.ndarray, C: np.ndarray, D: np.ndarray):
        self.p, self.G = len(A), len(D)
        self.C = C
        self.LD = np.linalg.cholesky(D) if self.G else D
        self.DiC = np.linalg.solve(D, C) if self.G else C  # D^-1 C, (G, k, p)
        S = A - C.reshape(-1, self.p).T @ self.DiC.reshape(-1, self.p)
        self.LS = np.linalg.cholesky(S)

    def _group_solve(self, e: np.ndarray) -> np.ndarray:
        """D^-1 e per group, e `(G, k)`."""
        return np.stack([cho_solve((L, True), v) for L, v in zip(self.LD, e)])

    def solve(self, b: np.ndarray) -> np.ndarray:
        a = b[: self.p]
        if not self.G:
            return cho_solve((self.LS, True), a)
        Die = self._group_solve(b[self.p :].reshape(self.G, -1))
        xa = cho_solve((self.LS, True), a - self.C.reshape(-1, self.p).T @ Die.ravel())
        return np.concatenate([xa, (Die - self.DiC @ xa).ravel()])

    def logdet(self) -> float:
        d = 2 * np.log(np.diagonal(self.LS)).sum()
        return d + (2 * np.log(np.diagonal(self.LD, axis1=1, axis2=2)).sum() if self.G else 0.0)

    def inverse_blocks(self) -> tuple[np.ndarray, np.ndarray]:
        """The dense block of H^-1 `(p, p)` and its diagonal blocks per group `(G, k, k)`."""
        Si = cho_solve((self.LS, True), np.eye(self.p))
        if not self.G:
            return Si, np.zeros((0, 0, 0))
        k = self.LD.shape[1]
        Di = np.stack([cho_solve((L, True), np.eye(k)) for L in self.LD])
        return Si, Di + (self.DiC @ Si) @ self.DiC.transpose(0, 2, 1)

    def draws(self, z: np.ndarray) -> np.ndarray:
        """`z` `(n_columns, n)` standard normal to draws of N(0, H^-1): with the groups ordered
        first, H = L L' with L = [[LD, 0], [C' LD^-T, LS]], and x = L^-T z."""
        xa = solve_triangular(self.LS, z[: self.p], lower=True, trans="T")
        if not self.G:
            return xa
        ze = z[self.p :].reshape(self.G, -1, z.shape[1])
        xe = [
            solve_triangular(
                L,
                zg - solve_triangular(L, Cg @ xa, lower=True),
                lower=True,
                trans="T",
            )
            for L, Cg, zg in zip(self.LD, self.C, ze)
        ]
        return np.concatenate([xa, np.concatenate(xe)])


def penalised_information(X: Matrix, w: np.ndarray, Q: Prior) -> Factor:
    A, C, D = X.information(w)
    return Factor(A + Q.dense, C, D + Q.group)


def information_weights(mu: np.ndarray, offset: np.ndarray, theta: float) -> np.ndarray:
    """The negative binomial's expected information per row on the log scale, with each day's mean
    taken as at least `MIN_RATE` birds per full day (times its coverage).

    Where a taxon is counted but never seen (Common Wood Pigeon on 18-25 July), its rate runs to
    zero and the information `mu` with it: the log-likelihood is flat below the mode but rises
    steeply above it, since zeros were counted. The Gaussian (Laplace) posterior is then nearly
    flat there (log-scale sd 12-16 against 0.4-0.5 in the passage) and a few draws held millions of
    birds. Evaluated at `MIN_RATE`, those days are as uncertain as a rate the zeros still allow;
    days carrying more birds are unchanged, and the mode, where the gradient is zero, too. Used for
    the posterior (and the Laplace marginal likelihood), not for the IRLS steps.
    """
    return np.maximum(mu, MIN_RATE * np.exp(offset)) * theta / (theta + mu)


def posterior_mode(X: Matrix, offset, y, Q: Prior, theta, beta):
    """Penalised IRLS (Fisher scoring with step halving) for the negative binomial log link:

    `(beta, H, penalised log-likelihood)`, H the penalised expected information (a `Factor`).
    """

    def objective(b):
        eta = np.minimum(offset + X @ b, ETA_MAX)
        return nb_loglik(y, np.exp(eta), theta).sum() - 0.5 * b @ (Q @ b)

    current = objective(beta)
    for _ in range(IRLS_MAX_ITER):
        mu = np.exp(np.minimum(offset + X @ beta, ETA_MAX))
        w = mu * theta / (theta + mu)
        grad = X.rmatvec((y - mu) * theta / (theta + mu)) - Q @ beta
        step = penalised_information(X, w, Q).solve(grad)
        for _ in range(30):
            new = objective(beta + step)
            if new >= current - 1e-10:
                break
            step /= 2
        beta = beta + step
        done = abs(new - current) < IRLS_TOL * (abs(current) + 1)
        current = new
        if done:
            break
    mu = np.exp(np.minimum(offset + X @ beta, ETA_MAX))
    return beta, penalised_information(X, information_weights(mu, offset, theta), Q), current


@dataclass
class Fit:
    design: Design
    hyper: dict
    beta: np.ndarray
    H: Factor
    log_marginal: float

    @property
    def theta(self) -> float:
        return float(np.exp(self.hyper["theta"]))

    def draws(self, n: int, rng: np.random.Generator) -> np.ndarray:
        """Posterior draws of the coefficients, `(n, k)`."""
        z = rng.standard_normal((len(self.beta), n))
        return self.beta + self.H.draws(z).T


def _logdet(A: np.ndarray) -> float:
    return 2 * np.log(np.diag(cho_factor(A, lower=True)[0])).sum()


def _laplace(design, X, offset, y, h, beta):
    Q = design.prior(h)
    beta, H, pll = posterior_mode(X, offset, y, Q, np.exp(h["theta"]), beta)
    log_marginal = pll + 0.5 * Q.logdet() - 0.5 * H.logdet()
    return log_marginal, beta, H


def _efs(d: Design, X, offset, y, hyper: dict, beta, theta_ml: bool):
    """Fellner-Schall iterations of the smoothing parameters to convergence, `theta` fixed or (with
    `theta_ml`) its maximum likelihood given the fitted means after each update.

    Two departures from the plain update, neither of which moves its fixed point: a step in the
    same direction as the last is lengthened (x2 each time, up to `EFS_MAX_GAIN`), so that a
    parameter heading for a bound gets there in a few iterations rather than tens; and a penalty
    that leaves its smooth less than `EFS_TOL` degrees of freedom (`lambda * room`: a fully
    penalised trend, a straight line) has converged whatever its `lambda` does, which near the
    upper bound is a ratio of two vanishing numbers and oscillates.
    """
    gain, last = dict.fromkeys(hyper, 1.0), dict.fromkeys(hyper, 0.0)
    for _ in range(EFS_MAX_ITER):
        beta, H, _ = posterior_mode(X, offset, y, d.prior(hyper), np.exp(hyper["theta"]), beta)
        dense_inv, group_inv = H.inverse_blocks()
        new, change = dict(hyper), 0.0
        for b in d.blocks:
            if not b.penalties:
                continue
            Q_inv = np.linalg.inv(d.block_precision(b, hyper))
            bb = beta[b.columns].reshape(b.groups, b.size)
            Hb = group_inv if b.groups > 1 else dense_inv[None, b.columns, b.columns]
            for name, S in b.penalties:
                # tr(Q^-1 S) - tr(H^-1 S), over every group
                room = b.groups * np.sum(Q_inv * S) - np.sum(Hb * S)
                wiggle = np.sum((bb @ S) * bb)
                step = np.log(max(room, 1e-12) / max(wiggle, 1e-12))
                gain[name] = min(2 * gain[name], EFS_MAX_GAIN) if step * last[name] > 0 else 1.0
                last[name] = step
                new[name] = float(np.clip(hyper[name] + gain[name] * step, *GAM_LOG_LAMBDA_BOUNDS))
                if np.exp(hyper[name]) * room >= EFS_TOL:
                    change = max(change, abs(new[name] - hyper[name]))
        if theta_ml:
            mu = np.exp(np.minimum(offset + X @ beta, ETA_MAX))
            new["theta"] = float(
                optimize.minimize_scalar(
                    lambda t: -nb_loglik(y, mu, np.exp(t)).sum(),
                    bounds=np.log(THETA_BOUNDS),
                    method="bounded",
                ).x
            )
            change = max(change, abs(new["theta"] - hyper["theta"]))
        hyper = new
        if change < EFS_TOL:
            break
    return hyper, beta


def fit(variant: str, frame: pd.DataFrame, years: np.ndarray, doy_range, hyper=None) -> Fit:
    """Fit a variant to the observed rows of `frame` (`c > 0`), over the year range `years`.

    The smoothing parameters by the extended Fellner-Schall update (Wood and Fasiolo 2017; mgcv's
    `optimizer = "efs"`), a fixed point of the Laplace marginal likelihood: each penalty's
    `lambda` is scaled by `(tr(Q^-1 S) - tr(H^-1 S)) / (beta' S beta)`, the degrees of freedom its
    prior allows minus those the data use, over the wiggliness it is asked to explain. `theta`
    starts at its maximum likelihood given the fitted means, which over-fit means make too large
    (days look less variable than they are), then alternates with the smoothing parameters as the
    maximum of the Laplace marginal likelihood. With `hyper` given, only the posterior mode is
    found: how the benchmark refits a year with part of it hidden.
    """
    d = design(variant, years, doy_range)
    obs = frame[frame["c"] > 0]
    g = d.grouped
    X = Matrix(
        d.X(obs["year"].to_numpy(), obs["doy"].to_numpy()),
        g.columns.start if g else d.n_columns,
        g.groups if g else 0,
    )
    offset = np.log(obs["c"].to_numpy())
    y = obs["y"].to_numpy(float)
    beta = np.zeros(d.n_columns)
    beta[0] = np.log(max(y.sum() / obs["c"].sum(), 1e-3))
    if hyper is None:
        hyper = {name: LOG_LAMBDA_START for name in d.hyper_names} | {"theta": 0.0}
        hyper, beta = _efs(d, X, offset, y, hyper, beta, theta_ml=True)
        for _ in range(THETA_ROUNDS):
            state = {"beta": beta}

            def negative(t):
                lm, b, _ = _laplace(d, X, offset, y, hyper | {"theta": t}, state["beta"])
                state["beta"] = b
                return -lm

            t = hyper["theta"]
            lo, hi = np.log(THETA_BOUNDS)
            theta = optimize.minimize_scalar(
                negative,
                bounds=(max(lo, t - 2), min(hi, t + 1)),
                method="bounded",
                options={"xatol": EFS_TOL},
            ).x
            hyper, beta = _efs(d, X, offset, y, hyper | {"theta": float(theta)}, beta, False)
            if abs(theta - t) < EFS_TOL:
                break
    lm, beta, H = _laplace(d, X, offset, y, hyper, beta)
    return Fit(d, dict(hyper), beta, H, lm)


# --- predictions ------------------------------------------------------------


def eta_draws(f: Fit, year, doy, beta_draws, without=()) -> np.ndarray:
    """Log rate per full-day equivalent, `(n_draws, n_rows)`, optionally leaving blocks out."""
    X = f.design.X(np.asarray(year), np.asarray(doy)).tocsc()
    if without:
        keep = np.ones(X.shape[1])
        keep[f.design.columns_of(*without)] = 0
        X = X @ sparse.diags(keep)
    return np.minimum(np.asarray((X @ beta_draws.T).T), ETA_MAX)


def hour_dispersion(
    days: pd.DataFrame,
    hourly: pd.DataFrame,
    effort: pd.DataFrame,
    profile: np.ndarray,
    block: int = KAPPA_BLOCK,
) -> float:
    """Hourly over-dispersion `kappa` of one taxon: within a day, the passage of a coverage `m` is.

    gamma(kappa * m, kappa / rate) a priori, so given the day's birds, the birds of its counted
    hours are Dirichlet-multinomial with weights kappa * the profile's share in the minutes counted
    of each hour. Maximum likelihood over the days timed to the hour (`PROFILE_MIN_TIMED`) with at
    least two counted hours; birds timed to an hour with no coverage are left out. Without such
    days: `KAPPA_BOUNDS[1]`, Poisson. Fitted on blocks of `block` solar hours (`KAPPA_BLOCK`), at
    least two of them counted.
    """
    d = days[(days["count"] > 0) & (days["timed"] >= PROFILE_MIN_TIMED)]
    counted = effort[effort["state"] == "counted"].set_index("date")["slots"]
    d = d[d["date"].isin(counted.index)].reset_index(drop=True)
    if d.empty:
        return KAPPA_BOUNDS[1]
    pos = pd.Series(d.index, index=d["date"])
    h = hourly[hourly["date"].isin(pos.index)]
    birds = np.zeros((len(d), HOURS))
    np.add.at(birds, (pos[h["date"]].to_numpy(), h["hour"].to_numpy()), h["count"].to_numpy())
    doy = d["date"].dt.dayofyear.clip(*PROFILE_DOY) - PROFILE_DOY[0]
    a = by_hour(profile[doy.to_numpy()] * np.stack(d["date"].map(counted).to_numpy()))
    on = a > 0
    birds = np.where(on, birds, 0)
    if block > 1:  # hours summed into blocks: a Dirichlet's kappa is unchanged by this
        a = a.reshape(len(a), -1, block).sum(axis=2)
        birds = birds.reshape(len(birds), -1, block).sum(axis=2)
        on = a > 0
    ok = (on.sum(axis=1) >= 2) & (birds.sum(axis=1) > 0)
    a, birds, on = a[ok], birds[ok], on[ok]
    if not len(a):
        return KAPPA_BOUNDS[1]
    total = birds.sum(axis=1)

    def nll(log_k):
        alpha = np.exp(log_k) * a
        A = alpha.sum(axis=1)
        per_hour = np.where(on, gammaln(birds + alpha) - gammaln(np.where(on, alpha, 1)), 0)
        return -(gammaln(A) - gammaln(total + A) + per_hour.sum(axis=1)).sum()

    res = optimize.minimize_scalar(nll, bounds=np.log(KAPPA_BOUNDS), method="bounded")
    return float(np.exp(res.x))


def fill_draws(
    f: Fit, frame: pd.DataFrame, c_total, n: int = DRAWS, seed: int = 0, kappa: float | None = None
) -> np.ndarray:
    """Draws of each row's birds over coverage `c_total` (1: the whole day), `(n, n_rows)`.

    The observed `y` (over `c`) is kept and the rest, `c_total - c`, drawn from the posterior
    predictive: the day's rate is gamma(theta, theta / mu) a priori, so given y birds in c it is
    gamma(theta + y, theta / mu + c), and the birds missed are Poisson of that rate times the
    coverage missed, or with `kappa` (`hour_dispersion`) of a gamma(kappa * missed, kappa / rate)
    passage: flocks. A day with nothing counted is a plain negative binomial draw (with `kappa`, a
    little wider). `extra` birds (counted below `MIN_COVERAGE`) are added as they are.

    The day's rate given y treats the counted hours as Poisson. Treating them as flocks too (y
    negative binomial of size kappa * c, which shrinks y and c by kappa c / (kappa c + y)) was
    tried: it biased gap-filled totals 3% low (`DECISIONS.md` -> Explore).
    """
    rng = np.random.default_rng(seed)
    beta = f.draws(n, rng)
    mu = np.exp(eta_draws(f, frame["year"], frame["doy"], beta))
    y, c = frame["y"].to_numpy(float), frame["c"].to_numpy(float)
    missing = np.clip(np.asarray(c_total, float) - c, 0, None)
    rate = rng.gamma(f.theta + y, 1 / (f.theta / mu + c))
    passage = rate * missing
    if kappa is not None:
        passage = rng.gamma(kappa * missing + 1e-300, passage / (kappa * missing + 1e-300))
    extra = frame["extra"].to_numpy(float) if "extra" in frame else 0.0
    return y + extra + rng.poisson(passage)


def annual_totals(
    f: Fit,
    frame: pd.DataFrame,
    n: int = DRAWS,
    seed: int = 0,
    kappa: float | None = None,
    totals: np.ndarray | None = None,
) -> pd.DataFrame:
    """Gap-filled window total per year (observed birds + the posterior predictive of the hours and
    days not counted), its quantiles, the observed share, and the smooth expected total: the trend
    in a typical year, every hour counted.

    Typical means without the year's own level but with its
    weather episodes on average: leaving them out altogether would lower the total, a sum of
    exponentials (Red Kite 2025: 13 500 against totals of 15 000-21 000), so the smooth is scaled
    by their mean effect, the geometric mean over years of total with episodes / total without.
    `totals`: the draws of `fill_draws` (whole days) if already made.
    """
    if totals is None:
        totals = fill_draws(f, frame, 1.0, n, seed, kappa)
    years = frame["year"].to_numpy()
    uy = np.unique(years)
    per_year = np.stack([totals[:, years == y].sum(axis=1) for y in uy], axis=1)
    beta = f.draws(n, np.random.default_rng(seed + 1))

    def per_year_sum(without):
        mu = np.exp(eta_draws(f, frame["year"], frame["doy"], beta, without=without))
        return np.stack([mu[:, years == y].sum(axis=1) for y in uy], axis=1)

    smooth = per_year_sum(SMOOTH)
    with_episodes = per_year_sum(("year",))
    smooth *= np.exp(np.log(with_episodes / smooth).mean(axis=1, keepdims=True))
    counted = frame["y"] + frame.get("extra", 0.0)
    out = pd.DataFrame({"year": uy, "observed": counted.groupby(frame["year"]).sum().to_numpy()})
    out["total"] = per_year.mean(axis=0)
    for q in QUANTILES:
        out[f"q{q * 100:g}"] = np.quantile(per_year, q, axis=0)
    out["observed_share"] = out["observed"] / out["total"]
    out["smooth"] = np.median(smooth, axis=0)
    out["smooth_q2.5"], out["smooth_q97.5"] = np.quantile(smooth, [0.025, 0.975], axis=0)
    return out


def season_curve(f: Fit, year: int, n: int = 400, seed: int = 0) -> pd.DataFrame:
    """Expected birds per full day by day of year in `year` (median and 95% band)."""
    doy = np.arange(f.design.doy_range[0], f.design.doy_range[1] + 1)
    beta = f.draws(n, np.random.default_rng(seed))
    mu = np.exp(eta_draws(f, np.full(len(doy), year), doy, beta, without=SMOOTH))
    q = np.quantile(mu, [0.025, 0.5, 0.975], axis=0)
    return pd.DataFrame({"doy": doy, "lo": q[0], "mid": q[1], "hi": q[2]})


def episode_curves(f: Fit, years: np.ndarray) -> dict:
    """Each year's season with and without its weather episodes, at the posterior mode: `base`,

    expected birds per full day without the episodes (trend, year level, season and shift), and
    `episode`, the episodes' log multiplier, so that `base * exp(episode)` is the full fit.
    """
    doy = np.arange(f.design.doy_range[0], f.design.doy_range[1] + 1)
    base, episode = {}, {}
    for yr in years:
        yd = (np.full(len(doy), yr), doy)
        full = eta_draws(f, *yd, f.beta[None])[0]
        without = eta_draws(f, *yd, f.beta[None], without=("episode",))[0]
        base[str(yr)] = np.round(np.exp(without), 1)
        episode[str(yr)] = np.round(full - without, 3)
    return {"doy": doy, "base": base, "episode": episode}


def peak_doy(f: Fit, years: np.ndarray, n: int = 400, seed: int = 0) -> pd.DataFrame:
    """Median passage date per year (the day by which half the season's birds have passed), of the
    smooth season (without the year's weather episodes), with its 80% band (`lo`, `hi`), and the
    80% bands of the dates by which the first and last `PASSAGE_QUANTILES` have passed (`q10_lo`,
    `q10_hi`, `q90_lo`, `q90_hi`)."""
    doy = np.arange(f.design.doy_range[0], f.design.doy_range[1] + 1)
    beta = f.draws(n, np.random.default_rng(seed))
    lo, hi = PASSAGE_QUANTILES[0], PASSAGE_QUANTILES[-1]
    rows = []
    for yr in years:
        mu = np.exp(eta_draws(f, np.full(len(doy), yr), doy, beta, without=SMOOTH))
        cum = mu.cumsum(axis=1) / mu.sum(axis=1, keepdims=True)
        date = {q: doy[np.minimum((cum < q).sum(axis=1), len(doy) - 1)] for q in (lo, 0.5, hi)}
        rows.append(
            {
                "year": yr,
                "lo": np.quantile(date[0.5], 0.1),
                "mid": np.median(date[0.5]),
                "hi": np.quantile(date[0.5], 0.9),
                **{
                    f"q{round(q * 100)}_{k}": np.quantile(date[q], b)
                    for q in (lo, hi)
                    for k, b in (("lo", 0.1), ("hi", 0.9))
                },
            }
        )
    return pd.DataFrame(rows)


def passage_quantiles(f: Fit, years: np.ndarray, qs=PASSAGE_QUANTILES) -> pd.DataFrame:
    """Day of year by which each share `qs` of the smooth season (without the year's level and
    weather episodes) has passed, per year, at the posterior mode: interpolated between days."""
    doy = np.arange(f.design.doy_range[0], f.design.doy_range[1] + 1)
    rows = []
    for yr in years:
        mu = np.exp(eta_draws(f, np.full(len(doy), yr), doy, f.beta[None], without=SMOOTH))[0]
        cum = np.cumsum(mu) / mu.sum()
        rows.append([yr, *(float(np.interp(q, cum, doy)) for q in qs)])
    return pd.DataFrame(rows, columns=["year", *(f"q{round(q * 100)}" for q in qs)]).round(1)


# --- export -----------------------------------------------------------------


def filled_days(frame: pd.DataFrame, totals: np.ndarray, c: pd.Series) -> pd.DataFrame:
    """Every day of `frame`: its coverage `c` (by date; 0 if not counted), the birds counted, and
    its gap-filled full-day total from the `fill_draws` draws `totals`: the mean and an 80%
    interval.

    The model's value for the day, counted or not: the birds counted plus the posterior predictive
    of the hours not counted, which leans on the day's own count when much of it was counted and
    on the neighbouring days and the season when little or none was. Days counted below
    `MIN_COVERAGE` are filled as days not counted, their birds added as they are. The mean, not
    the median: the birds missed come in flocks, so their median is low (Common Wood Pigeon, c
    0.5-0.9: medians sum to 0.87 of Σ count / c, means to 1.00), and means add up over days.
    """
    q = np.quantile(totals, FILLED_DAY_INTERVAL, axis=0)
    return pd.DataFrame(
        {
            "date": frame["date"].to_numpy(),
            "c": frame["date"].map(c).fillna(0.0).round(3).to_numpy(),
            "count": (frame["y"] + frame["extra"]).to_numpy(),
            "total": totals.mean(axis=0).round(1),
            "q10": q[0],
            "q90": q[1],
        }
    )


def taxon_trend(
    days: pd.DataFrame,
    hourly: pd.DataFrame,
    effort: pd.DataFrame,
    profile: np.ndarray,
    first_year: int,
    last_year: int,
    variant: str = "gam",
    seed: int = 0,
    window: tuple[int, int] | None = None,
) -> dict:
    """One taxon's trend for `species/<taxon_id>.json`, `first_year`..`last_year` (complete seasons
    from its start year), over `window` (season days; the default window if None).

    `annual`: per year, birds counted, the gap-filled total (median and 80%/95% intervals), the
    share counted, and the smooth expected total (trend without the year's level and weather
    episodes, every hour counted) with its 95% band. `passage`: the smooth season's median passage
    date (day of year, 80% band), and the 80% bands of the 10% and 90% dates. `passage_q`: the
    smooth season's 10/50/90% passage dates (`passage_quantiles`). `season`: expected birds per
    full day on each day of year, smooth, in the first and the last year. `episodes`: every year's
    season with and without its weather episodes (`episode_curves`). `days`: every window day's
    gap-filled full-day total (`filled_days`), from the same draws as `annual`. `season_fill`: the
    same series for the season block (`season.season_block`): each day's mean, every year's passage
    dates over the draws (`season.passage_from_draws`) and the chances
    (`season.chance_from_draws`); not exported as such.
    """
    c = coverage(effort, profile)
    frame = model_frame(days, effort, c, first_year, last_year, window)
    years = np.arange(first_year, last_year + 1)
    doy_range = (int(frame["doy"].min()), int(frame["doy"].max()))
    kappa = hour_dispersion(days, hourly, effort, profile)
    f = fit(variant, frame, years, doy_range)
    totals = fill_draws(f, frame, 1.0, seed=seed, kappa=kappa)
    a = annual_totals(f, frame, seed=seed, kappa=kappa, totals=totals)
    a = a.drop(columns=["total"]).rename(columns={"q50": "total"})
    a["observed_share"] = a["observed"] / a["total"]
    birds = [c for c in a.columns if c not in ("year", "observed_share")]
    a[birds] = a[birds].round(0)
    curves = {
        str(y): season_curve(f, y, seed=seed)["mid"].to_numpy() for y in (first_year, last_year)
    }
    return {
        "model": variant,
        "first_year": first_year,
        "last_year": last_year,
        "theta": f.theta,
        "kappa": kappa,
        "annual": a,
        "passage": peak_doy(f, years, seed=seed),
        "passage_q": passage_quantiles(f, years),
        "season": {"doy": np.arange(doy_range[0], doy_range[1] + 1), **curves},
        "episodes": episode_curves(f, years) if variant == "gam" else None,
        "days": filled_days(frame, totals, c.set_axis(effort["date"])),
        "season_fill": {
            "days": frame[["year", "doy"]].assign(total=totals.mean(axis=0)),
            "passage": passage_from_draws(frame, totals),
            "chance": chance_from_draws(frame, totals, last_year),
        },
    }
