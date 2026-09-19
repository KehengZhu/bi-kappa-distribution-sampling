"""Frozen statistical machinery for Experiment 6.

Every threshold, transform and interval rule in this file was fixed before any production
number was looked at, and is not to be retuned afterwards.  The repository has no shared
statistics module, so these are written here rather than borrowed.

Three house rules from the rest of the experiment suite apply throughout and are the
reason several functions look indirect:

*  never form ``R**2``.  It overflows for ``R > 1.3e154`` while ``R`` itself is perfectly
   representable, which silently discards exactly the heavy-tail draws a radial diagnostic
   exists to test.  The radial variable is ``W = 1/(1+R^2)``, computed as
   ``expit(-2 log R)`` and carried in the log domain as ``log W = -logaddexp(0, 2 log R)``.
*  never use ``Y = R^2/(1+R^2) ~ Beta(3/2, a)``.  At small ``a`` a large fraction of ``Y``
   rounds to exactly 1 and the test statistic becomes meaningless.  ``W ~ Beta(a, 3/2)``
   does not have that problem: its mass sits away from both endpoints.
*  never compute a rejected fraction as ``1 - P(accept)``.

The upper-tail probability of the radius is
``q = Pr(R > r) = Pr(W < w) = I_w(a, 3/2)``, and ``Z = -log q`` is exactly ``Exp(1)`` under
the target.  ``Z`` is the diagnostic of record at the smallest shapes, because ``W`` itself
rounds to zero there while ``log W`` does not.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy import special, stats

# ---------------------------------------------------------------------------
# Predeclared constants.  Fixed before any result was seen.  Not to be adjusted.
# ---------------------------------------------------------------------------
FAMILYWISE_ALPHA = 0.01          # Holm-corrected familywise type-I rate per configuration
INTERVAL_CONF = 0.95             # two-sided interval level for failure probabilities
TAIL_INTERVAL_CONF = 0.99        # bootstrap level for tail-ratio and method differences
BOOTSTRAP_RESAMPLES = 10000
QUANTILE_LEVELS = (0.5, 0.9, 0.99, 0.999, 0.9999)
TAIL_LEVELS = (0.99, 0.999, 0.9999)
LOG_W_SWITCH = -30.0             # below this, the small-W log form is used for log q
LOG_Q_ORACLE_TOL = 1e-10         # required agreement with the arbitrary-precision oracle
BETA_B = 1.5                     # the second Beta shape is always 3/2


# ---------------------------------------------------------------------------
# Binomial intervals
# ---------------------------------------------------------------------------
def clopper_pearson(k: int, n: int, conf: float = INTERVAL_CONF) -> tuple[float, float]:
    """Exact two-sided Clopper-Pearson interval for a binomial proportion."""
    if n == 0:
        return (float("nan"), float("nan"))
    alpha = 1.0 - conf
    lo = 0.0 if k == 0 else float(special.betaincinv(k, n - k + 1, alpha / 2))
    hi = 1.0 if k == n else float(special.betaincinv(k + 1, n - k, 1 - alpha / 2))
    return (lo, hi)


def zero_count_upper(n: int, conf: float = INTERVAL_CONF) -> float:
    """One-sided upper limit for a rate after observing zero events in ``n`` trials.

    ``1 - (1-conf)^(1/n)``, the exact form of the rule of three.  A zero count is reported
    as this bound, never as "zero probability".
    """
    if n <= 0:
        return float("nan")
    return 1.0 - (1.0 - conf) ** (1.0 / n)


def rate_with_interval(k: int, n: int, conf: float = INTERVAL_CONF) -> dict:
    """Point estimate plus interval, with zero counts reported as an upper bound."""
    if n == 0:
        return {"rate": float("nan"), "lo": float("nan"), "hi": float("nan"),
                "is_upper_bound": False, "n": 0, "k": 0}
    if k == 0:
        hi = zero_count_upper(n, conf)
        return {"rate": 0.0, "lo": 0.0, "hi": hi, "is_upper_bound": True, "n": n, "k": 0}
    lo, hi = clopper_pearson(k, n, conf)
    return {"rate": k / n, "lo": lo, "hi": hi, "is_upper_bound": False, "n": n, "k": k}


# ---------------------------------------------------------------------------
# The radial transforms
# ---------------------------------------------------------------------------
def log_w_from_log_r(log_r: np.ndarray) -> np.ndarray:
    """``log W`` from ``log R`` without ever forming ``R`` or ``R**2``.

    ``W = 1/(1+R^2)`` so ``log W = -log(1 + exp(2 log R))``, evaluated with logaddexp so
    that it is exact to rounding for every representable ``log R``.
    """
    log_r = np.asarray(log_r, dtype=float)
    return -np.logaddexp(0.0, 2.0 * log_r)


def log_r_from_log_w(log_w: np.ndarray) -> np.ndarray:
    """Inverse of :func:`log_w_from_log_r`:  ``log R = (log(1-W) - log W)/2``."""
    log_w = np.asarray(log_w, dtype=float)
    return 0.5 * (np.log(-np.expm1(log_w)) - log_w)


def log_beta(a: float, b: float = BETA_B) -> float:
    return float(special.betaln(a, b))


def log_q_from_log_w(log_w: np.ndarray, a: float, b: float = BETA_B) -> np.ndarray:
    """``log I_w(a, b)``: the log upper-tail probability of the radius.

    The evaluator is frozen in two pieces, switching at ``LOG_W_SWITCH``:

    * ``log W > -30``: form ``W`` and take ``log(betainc(a, b, W))``.  Any library result
      that is not strictly positive and finite is rejected rather than patched.
    * ``log W <= -30``: use the small-``W`` log form
      ``log q = a log W - log a - log B(a, b)``.  The omitted correction is ``O(W)``, i.e.
      below ``1e-13`` at the switch.

    Both pieces, and the switch itself, are validated against an arbitrary-precision
    incomplete beta before any sampling result is used; see
    :func:`validate_log_q_evaluator`.
    """
    log_w = np.asarray(log_w, dtype=float)
    out = np.empty_like(log_w)
    hi = log_w > LOG_W_SWITCH
    if np.any(hi):
        w = np.exp(log_w[hi])
        val = special.betainc(a, b, w)
        with np.errstate(divide="ignore"):
            out[hi] = np.where(val > 0, np.log(val), -np.inf)
    lo = ~hi
    if np.any(lo):
        out[lo] = a * log_w[lo] - math.log(a) - log_beta(a, b)
    return out


def z_from_log_w(log_w: np.ndarray, a: float, b: float = BETA_B) -> np.ndarray:
    """``Z = -log I_w(a, b)``, which is exactly ``Exp(1)`` under the target law."""
    return -log_q_from_log_w(log_w, a, b)


def validate_log_q_evaluator(shapes, b: float = BETA_B, span: float = 10.0,
                             points: int = 41, tol: float = LOG_Q_ORACLE_TOL) -> list[dict]:
    """Check the frozen evaluator against ``mpmath`` on both sides of the switch.

    The grid spans at least ``span`` natural-log units on each side of ``LOG_W_SWITCH``
    for every shape in ``shapes``.  Run before sampling results are inspected; if it fails
    the evaluator is replaced, never the switch retuned.
    """
    import mpmath as mp

    mp.mp.dps = 60
    rows = []
    grid = np.linspace(LOG_W_SWITCH - span, LOG_W_SWITCH + span, points)
    for a in shapes:
        ours = log_q_from_log_w(grid, a, b)
        worst = 0.0
        worst_at = float("nan")
        for lw, mine in zip(grid, ours):
            w = mp.e ** mp.mpf(float(lw))
            exact = mp.log(mp.betainc(mp.mpf(float(a)), mp.mpf(float(b)), 0, w,
                                      regularized=True))
            d = abs(float(exact) - float(mine))
            if d > worst:
                worst, worst_at = d, float(lw)
        rows.append({
            "shape_a": a,
            "log_w_min": float(grid[0]),
            "log_w_max": float(grid[-1]),
            "points": points,
            "max_abs_log_q_error": worst,
            "worst_at_log_w": worst_at,
            "tolerance": tol,
            "pass": bool(worst < tol),
        })
    return rows


# ---------------------------------------------------------------------------
# Goodness of fit
# ---------------------------------------------------------------------------
@dataclass
class GofResult:
    name: str
    statistic: float
    pvalue: float
    n: int

    def as_row(self) -> dict:
        return {"test": self.name, "statistic": self.statistic, "pvalue": self.pvalue,
                "n": self.n}


def gof_exponential(z: np.ndarray, label: str = "Z") -> list[GofResult]:
    """KS and Cramer-von Mises against the unit exponential."""
    z = np.asarray(z, dtype=float)
    z = z[np.isfinite(z)]
    if z.size < 2:
        return []
    ks = stats.kstest(z, "expon")
    cvm = stats.cramervonmises(z, "expon")
    return [GofResult(f"{label}:KS", float(ks.statistic), float(ks.pvalue), z.size),
            GofResult(f"{label}:CvM", float(cvm.statistic), float(cvm.pvalue), z.size)]


def gof_beta(w: np.ndarray, a: float, b: float = BETA_B, label: str = "W") -> list[GofResult]:
    """KS and Cramer-von Mises against ``Beta(a, b)``, on resolved values only."""
    w = np.asarray(w, dtype=float)
    w = w[np.isfinite(w) & (w > 0) & (w < 1)]
    if w.size < 2:
        return []
    ks = stats.kstest(w, "beta", args=(a, b))
    cvm = stats.cramervonmises(w, "beta", args=(a, b))
    return [GofResult(f"{label}:KS", float(ks.statistic), float(ks.pvalue), w.size),
            GofResult(f"{label}:CvM", float(cvm.statistic), float(cvm.pvalue), w.size)]


def gof_uniform(x: np.ndarray, lo: float, hi: float, label: str) -> list[GofResult]:
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if x.size < 2:
        return []
    ks = stats.kstest(x, "uniform", args=(lo, hi - lo))
    cvm = stats.cramervonmises(x, "uniform", args=(lo, hi - lo))
    return [GofResult(f"{label}:KS", float(ks.statistic), float(ks.pvalue), x.size),
            GofResult(f"{label}:CvM", float(cvm.statistic), float(cvm.pvalue), x.size)]


def holm(pvalues, alpha: float = FAMILYWISE_ALPHA) -> list[bool]:
    """Holm step-down.  Returns True where the null is rejected at familywise ``alpha``."""
    p = np.asarray(list(pvalues), dtype=float)
    m = p.size
    if m == 0:
        return []
    order = np.argsort(p)
    reject = np.zeros(m, dtype=bool)
    for rank, idx in enumerate(order):
        if p[idx] <= alpha / (m - rank):
            reject[idx] = True
        else:
            break
    return [bool(x) for x in reject]


# ---------------------------------------------------------------------------
# Quantiles
# ---------------------------------------------------------------------------
def order_statistic_interval(n: int, p: float, conf: float = INTERVAL_CONF,
                             bonferroni: int = 1) -> tuple[int, int] | None:
    """Exact binomial order-statistic index bracket for the ``p`` quantile.

    Returns 0-based indices ``(lo, hi)`` into the sorted sample such that the interval
    ``[x_(lo), x_(hi)]`` covers the true quantile with at least ``conf`` probability,
    Bonferroni-corrected over ``bonferroni`` simultaneous levels.  ``None`` when the
    sample cannot resolve the level at all, which is reported rather than extrapolated.
    """
    if n <= 0 or not (0.0 < p < 1.0):
        return None
    alpha = (1.0 - conf) / bonferroni
    lo = int(stats.binom.ppf(alpha / 2, n, p)) - 1
    hi = int(stats.binom.isf(alpha / 2, n, p))
    lo = max(lo, 0)
    hi = min(hi, n - 1)
    if hi <= lo or hi >= n:
        return None
    # An upper index that coincides with the largest order statistic cannot bound the
    # quantile from above; that is the unresolved case.
    if hi == n - 1 and p > 1 - 3.0 / n:
        return None
    return (lo, hi)


def quantile_with_interval(sample: np.ndarray, p: float, conf: float = INTERVAL_CONF,
                           bonferroni: int = 1) -> dict:
    x = np.sort(np.asarray(sample, dtype=float))
    x = x[np.isfinite(x)]
    n = x.size
    if n == 0:
        return {"p": p, "estimate": float("nan"), "lo": float("nan"), "hi": float("nan"),
                "n": 0, "resolved": False}
    est = float(np.quantile(x, p))
    br = order_statistic_interval(n, p, conf, bonferroni)
    if br is None:
        return {"p": p, "estimate": est, "lo": float("nan"), "hi": float("nan"),
                "n": n, "resolved": False}
    return {"p": p, "estimate": est, "lo": float(x[br[0]]), "hi": float(x[br[1]]),
            "n": n, "resolved": True}


def exact_log_r_quantile(p: float, a: float, b: float = BETA_B) -> float:
    """The target ``log R`` quantile, obtained in the log domain and never from ``R**2``.

    ``Pr(R <= r) = p``  <=>  ``Pr(W < w) = 1 - p``  with  ``w = I^{-1}_{1-p}(a, b)``.

    Inverting the incomplete beta directly is only usable while ``w`` is representable.  At
    the shapes this experiment lives at it is not: for ``a = 1e-4`` the median of ``W`` is
    ``exp(-6931)``, so ``betaincinv`` returns a denormal or zero and every quantile from
    ``p = 0.5`` upwards collapses onto the same saturated value.  Below the same
    ``LOG_W_SWITCH`` used by :func:`log_q_from_log_w`, the small-``W`` form is inverted in
    closed form instead:

        ``log q = a log w - log a - log B(a, b)``  =>
        ``log w = (log q + log a + log B(a, b)) / a``,

    whose omitted ``O(w)`` correction shifts ``log w`` by about ``w/2``, i.e. below ``1e-13``
    at the switch.
    """
    if not (0.0 < p < 1.0) or a <= 0.0:
        return float("nan")
    log_q = math.log1p(-p)
    log_w_asym = (log_q + math.log(a) + log_beta(a, b)) / a
    if log_w_asym <= LOG_W_SWITCH:
        log_w = log_w_asym
    else:
        w = float(special.betaincinv(a, b, 1.0 - p))
        if not (0.0 < w < 1.0):
            return float("nan")
        log_w = math.log(w)
    return 0.5 * (math.log1p(-math.exp(log_w)) - log_w)


def mad_ratio_interval(x: np.ndarray, y: np.ndarray, conf: float = INTERVAL_CONF,
                       resamples: int = 1000, max_n: int = 50000,
                       rng_seed: int = 60003) -> dict:
    """Bootstrap percentile interval for ``MAD(x)/MAD(y)``.

    A fixed relative tolerance cannot serve here: the sampling error of a MAD ratio scales
    as ``n^-1/2``, so any constant would be far too tight at one sample size and far too
    loose at another.  The bootstrap is taken on a deterministic subsample of at most
    ``max_n`` pairs, because a MAD is a full pass over the data and the samples run to
    5e5 values per configuration; the subsample size is recorded alongside the interval.
    """
    rng = np.random.default_rng(rng_seed)
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    n = x.size
    if n < 100:
        return {"ratio": float("nan"), "lo": float("nan"), "hi": float("nan"),
                "n": int(n), "subsample": 0, "resamples": 0}
    if n > max_n:
        idx = rng.choice(n, size=max_n, replace=False)
        x, y = x[idx], y[idx]
    draws = np.empty(resamples, dtype=float)
    for i in range(resamples):
        j = rng.integers(0, x.size, x.size)
        my = mad(y[j])
        draws[i] = mad(x[j]) / my if my else np.nan
    alpha = 1.0 - conf
    lo, hi = np.nanpercentile(draws, [100 * alpha / 2, 100 * (1 - alpha / 2)])
    my = mad(y)
    return {"ratio": mad(x) / my if my else float("nan"), "lo": float(lo), "hi": float(hi),
            "n": int(n), "subsample": int(x.size), "resamples": resamples}


# ---------------------------------------------------------------------------
# Bootstrap
# ---------------------------------------------------------------------------
def seed_stratified_bootstrap(values_by_seed, statistic, resamples: int = BOOTSTRAP_RESAMPLES,
                              conf: float = TAIL_INTERVAL_CONF, rng_seed: int = 60001) -> dict:
    """Percentile interval from a seed-stratified bootstrap.

    ``values_by_seed`` is a sequence of per-seed arrays; each resample draws, within every
    seed, a sample of that seed's own size.  Preserving the seed strata is what keeps the
    interval honest when the seeds differ systematically.
    """
    rng = np.random.default_rng(rng_seed)
    arrays = [np.asarray(v, dtype=float) for v in values_by_seed]
    arrays = [v[np.isfinite(v)] for v in arrays]
    arrays = [v for v in arrays if v.size]
    if not arrays:
        return {"estimate": float("nan"), "lo": float("nan"), "hi": float("nan"),
                "resamples": 0}
    point = float(statistic([np.concatenate(arrays)]))
    draws = np.empty(resamples, dtype=float)
    for i in range(resamples):
        boot = [v[rng.integers(0, v.size, v.size)] for v in arrays]
        draws[i] = statistic(boot)
    alpha = 1.0 - conf
    lo, hi = np.nanpercentile(draws, [100 * alpha / 2, 100 * (1 - alpha / 2)])
    return {"estimate": point, "lo": float(lo), "hi": float(hi), "resamples": resamples}


def seed_stratified_bootstrap_indicator(counts, sizes, scale: float = 1.0,
                                        resamples: int = BOOTSTRAP_RESAMPLES,
                                        conf: float = TAIL_INTERVAL_CONF,
                                        rng_seed: int = 60001) -> dict:
    """Seed-stratified bootstrap of a pooled proportion, done exactly rather than by
    resampling the data.

    The statistics this is used for are means of 0/1 indicators over samples of 10^6 draws
    per seed.  Resampling those arrays 10000 times would be 5x10^10 element copies per
    interval and would never finish.  It is also unnecessary: the bootstrap distribution of
    a mean of indicators within a stratum of size ``n`` and observed proportion ``p`` is
    exactly ``Binomial(n, p)/n``, so the resampling can be done on the counts.  The result
    is the same estimator, the same strata and the same percentile interval.

    ``scale`` divides the pooled proportion, which is how the tail-retention ratios turn a
    probability into a ratio against ``1-p``.
    """
    rng = np.random.default_rng(rng_seed)
    counts = np.asarray(list(counts), dtype=np.int64)
    sizes = np.asarray(list(sizes), dtype=np.int64)
    keep = sizes > 0
    counts, sizes = counts[keep], sizes[keep]
    if sizes.size == 0 or sizes.sum() == 0:
        return {"estimate": float("nan"), "lo": float("nan"), "hi": float("nan"),
                "resamples": 0}
    total = float(sizes.sum())
    point = float(counts.sum()) / total / scale
    draws = np.zeros(resamples, dtype=float)
    for n_i, k_i in zip(sizes, counts):
        draws += rng.binomial(int(n_i), float(k_i) / float(n_i), size=resamples)
    draws = draws / total / scale
    alpha = 1.0 - conf
    lo, hi = np.percentile(draws, [100 * alpha / 2, 100 * (1 - alpha / 2)])
    return {"estimate": point, "lo": float(lo), "hi": float(hi), "resamples": resamples}


# ---------------------------------------------------------------------------
# Independence
# ---------------------------------------------------------------------------
def binned_independence_test(x: np.ndarray, y: np.ndarray, nbins: int = 8) -> GofResult:
    """Chi-square of independence on equal-count bins of ``x`` and equal-width bins of ``y``.

    Equal-count bins in ``x`` (the radius, whose scale spans hundreds of log units) and
    equal-width bins in ``y`` (a direction cosine on a bounded range) keep every cell
    populated without tuning.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    if x.size < 100:
        return GofResult("independence:chi2", float("nan"), float("nan"), int(x.size))
    xedges = np.quantile(x, np.linspace(0, 1, nbins + 1))
    xedges[0] -= 1e-9
    xedges[-1] += 1e-9
    yedges = np.linspace(y.min() - 1e-12, y.max() + 1e-12, nbins + 1)
    table, _, _ = np.histogram2d(x, y, bins=(xedges, yedges))
    keep_r = table.sum(axis=1) > 0
    keep_c = table.sum(axis=0) > 0
    table = table[np.ix_(keep_r, keep_c)]
    if table.shape[0] < 2 or table.shape[1] < 2:
        return GofResult("independence:chi2", float("nan"), float("nan"), int(x.size))
    chi2, p, _, _ = stats.chi2_contingency(table)
    return GofResult("independence:chi2", float(chi2), float(p), int(x.size))


def mad(x: np.ndarray) -> float:
    """Median absolute deviation about the median.

    Moment-free on purpose: for ``kappa <= 3/2`` the loaded population has no finite
    variance, so a scale comparison built on second moments would be estimating something
    that does not exist.
    """
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return float("nan")
    return float(np.median(np.abs(x - np.median(x))))
