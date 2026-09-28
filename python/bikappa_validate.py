#!/usr/bin/env python3
"""Distributional validation of a bi-Kappa velocity sample.

This module tests a *sample*, not a sampler.  Nothing in it assumes the
Gamma-ratio construction, so it applies to velocities produced by any bi-Kappa
loader -- inverse transform, acceptance-rejection, or scale mixture.

Input contract
--------------
``v``            (N, 3) array of velocities, in the simulation's global frame,
                 with N >= MIN_SAMPLES (1000).
``kappa``        shape index, > 1/2.
``theta_perp``   perpendicular thermal speed, > 0, in the units of ``v``.
``theta_par``    parallel thermal speed, > 0, in the units of ``v``.
``bhat``         magnetic-field direction, or ``None`` if the sample is already
                 field-aligned (i.e. the third column is the parallel one).
``n_attempts``   total draws attempted, if the loader enforced a velocity bound
                 by discarding and redrawing.  Supplying it lets the report
                 quote the rejected fraction, which for a hard bound equals the
                 total-variation distance from the uncapped distribution.
``alpha``        family-wise significance level, 0 < alpha < 1.

The thermal speeds follow the convention of Eq. (2) of the accompanying paper,

    f(v) ~ [1 + v_perp^2/(kappa theta_perp^2) + v_par^2/(kappa theta_par^2)]^-(kappa+1),

normalizable for every kappa > 1/2.  A thermal speed is a speed, not a
temperature.  A temperature-linked convention such as
theta^2 = 2[(kappa-3/2)/kappa](k_B T/m) is defined only for kappa > 3/2 and is
*not* interchangeable with this one; convert before calling.

The target is the uncapped distribution above.  A sample from a loader that
rejects draws outside a velocity bound follows a different, capped distribution;
the tests then measure the effect of the bound.

What the tests establish, and what they do not
----------------------------------------------
The cell-count test of the accompanying paper (Sec. VI B) and the radial,
directional and independence tests below together pin down the
three-dimensional distribution, because a spherically symmetric density is
determined by its radial distribution plus a uniform, radius-independent
direction.  They are distribution-free and use no moments, so they remain valid
on 1/2 < kappa <= 3/2 where the second moment does not exist.  A pass is
evidence that the sample is consistent with the target at the given size; it is
not a proof of correctness, and it says nothing about a bound the loader may
have applied -- for that, supply ``n_attempts`` and read the tail quantiles.

Verdicts
--------
The nine distributional tests (seven when the sample is too small for the
cell-count test) form one family.  Each has a p-value, and the family is
corrected with the Holm-Bonferroni step-down procedure: a test fails when its
Holm-adjusted p-value is at most ``alpha``, which keeps the probability that a
correct sample fails any of them at or below ``alpha`` (to the accuracy of the
asymptotic p-values).  Two further tests have
their own criteria: the anisotropy ratio fails when it is more than four
standard errors from its expected value, and the finiteness test fails on any
non-finite draw.  The tail quantiles and the largest normalized speeds are
diagnostics and carry no verdict.

Usage
-----
    from bikappa_validate import validate_sample, format_report
    report = validate_sample(v, kappa=2.0, theta_perp=1.0, theta_par=2.0)
    print(format_report(report))

or from the command line (``--help`` lists the options and input formats):

    python bikappa_validate.py sample.npy --kappa 2 --theta-perp 1 --theta-par 2

Exit status is 0 if every test passes, 1 if at least one test fails, and 2 if
the input or the options are invalid.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
import warnings

import numpy as np
from scipy import special, stats

__all__ = ["validate_sample", "format_report", "load_sample", "field_basis",
           "log_radius_quantile", "MIN_SAMPLES"]

# Smallest sample the tests accept.
MIN_SAMPLES = 1000

# Quantile probes for the radial tail.  The two extreme ones are where a radial
# error shows up first and where Cartesian marginals are least sensitive.
QUANTILE_PROBES = (0.5, 0.9, 0.99, 0.999)

# Direction cells for the joint-direction and independence tests.  The solid-angle
# element is d(cos Theta) d(Phi), so equal intervals of cos(Theta) crossed with
# equal intervals of Phi give cells of equal solid angle.  Binning both angles is
# what makes the independence test cover the whole direction: with cos(Theta)
# alone, a dependence between radius and azimuth would pass unseen.
N_RADIAL_BINS = 4
N_COS_BINS = 5
N_PHI_BINS = 8

# The accompanying paper's test (Sec. VI B): Pearson's chi^2 over cells of equal
# probability under the target -- N_SHELLS radius shells bounded by the exact
# deciles of R, crossed with the equal-solid-angle direction cells above.  Every
# one of the N_SHELLS * 40 cells then has probability 1/400.
N_SHELLS = 10
MIN_EXPECTED_PER_CELL = 5.0

# The test of the radius alone divides the outermost decile shell, R > r_9, further
# at the radii that R exceeds with these probabilities (its 99th ... 99.999th
# percentiles), so that the fastest 10% of draws are resolved in radius rather than
# counted in one shell.  An edge is used only if the shell beyond it expects at least
# MIN_EXPECTED_PER_CELL draws; at the paper's 5 x 10^7 draws all four are used.
TAIL_EXCEEDANCE = (1e-2, 1e-3, 1e-4, 1e-5)


def log_radius_quantile(kappa: float, exceed) -> np.ndarray:
    """log of the radius r with P(R > r) = exceed under the target, elementwise.

    P(R > r) = P(W < w) with W = 1/(1+R^2) ~ Beta(kappa-1/2, 3/2), and R^2 = (1-w)/w.
    Near kappa = 1/2 the edge w lies below the smallest double -- w ~ 10^-500 for
    the 10^-5 edge at kappa = 0.51, and even the deciles underflow once
    kappa - 1/2 < 3e-3 -- so Beta.ppf returns 0 there.  For such w the leading term
    of I_w(a, b) = w^a / (a B(a, b)) (1 + O(w)) is exact to double precision and
    gives log w directly.  Comparing log R with these edges then needs no bounded
    transform of either, so nothing underflows however large R is.
    """
    a = kappa - 0.5
    exceed = np.atleast_1d(np.asarray(exceed, dtype=float))
    with np.errstate(divide="ignore"):
        log_w = np.log(stats.beta.ppf(exceed, a, 1.5))
    tiny = ~(log_w > math.log(1e-280))
    log_w[tiny] = (np.log(exceed[tiny]) + math.log(a) + special.betaln(a, 1.5)) / a
    return 0.5 * (np.log1p(-np.exp(log_w)) - log_w)


def log_shell_edges(kappa: float) -> np.ndarray:
    """log of the radii splitting R into N_SHELLS shells of equal probability."""
    return log_radius_quantile(kappa, 1.0 - np.arange(1, N_SHELLS) / N_SHELLS)


def radius_shells(kappa: float, n: int) -> tuple[np.ndarray, np.ndarray]:
    """log edges and probabilities of the shells for the test of the radius alone.

    The N_SHELLS decile shells, with the outermost divided at those TAIL_EXCEEDANCE
    edges beyond which n draws expect at least MIN_EXPECTED_PER_CELL.
    """
    tail = [t for t in TAIL_EXCEEDANCE if n * t >= MIN_EXPECTED_PER_CELL]
    edges = np.concatenate([log_shell_edges(kappa), log_radius_quantile(kappa, tail)])
    exceed = np.concatenate([1.0 - np.arange(N_SHELLS) / N_SHELLS, tail, [0.0]])
    return edges, exceed[:-1] - exceed[1:]


# ----------------------------------------------------------------------------
# frame handling
# ----------------------------------------------------------------------------

def field_basis(bhat: np.ndarray) -> np.ndarray:
    """Orthonormal basis (e1, e2, e3) with e3 = bhat, as columns.

    Built the same way as the released C++ loader: start from (1,1,1), replace
    the component of largest |b_k| so the result is orthogonal to bhat.  Picking
    the largest component keeps the division safe -- for a unit vector
    |b_k| >= 1/sqrt(3) -- including the axis-aligned cases.

    Any orthonormal basis with e3 = bhat serves for the tests below, since they
    depend on the perpendicular directions only through rotationally invariant
    statistics.
    """
    b = np.asarray(bhat, dtype=float)
    norm = np.linalg.norm(b)
    if not np.isfinite(norm) or norm == 0.0:
        raise ValueError("bhat must be a nonzero finite vector")
    e3 = b / norm

    e2 = np.ones(3)
    k = int(np.argmax(np.abs(e3)))
    e2[k] = 1.0 - e3.sum() / e3[k]
    e2 /= np.linalg.norm(e2)

    e1 = np.cross(e2, e3)
    return np.column_stack([e1, e2, e3])


# ----------------------------------------------------------------------------
# the tests
# ----------------------------------------------------------------------------

def mad_relative_sigma(kappa: float, n: int) -> float:
    """Asymptotic relative standard error of a sample MAD, per component.

    Each Cartesian marginal of Eq. (2) is a scaled Student-t distribution with
    nu = 2 kappa - 1 degrees of freedom, so in units of its own scale the MAD is
    t0 = ppf(3/4, nu).  For a symmetric distribution the MAD is the median of |X|, whose
    density at t0 is 2 f(t0), giving

        sqrt(n) (MAD_n - t0) -> N(0, 1 / (16 f(t0)^2)).

    The relative error therefore scales as 1/(4 sqrt(n) t0 f(t0)), which blows up
    as kappa -> 1/2 because the marginal flattens: a fixed percentage tolerance
    would be far too tight there and far too loose at large kappa.
    """
    nu = 2.0 * kappa - 1.0
    t0 = float(stats.t.ppf(0.75, nu))
    f0 = float(stats.t.pdf(t0, nu))
    return 1.0 / (4.0 * math.sqrt(n) * t0 * f0)


def _holm_adjusted(pvalues: dict) -> dict:
    """Holm-Bonferroni adjusted p-values over a family of p-values.

    With the m p-values sorted as p_(1) <= ... <= p_(m), the adjusted value of
    p_(i) is max_{j <= i} min(1, (m - j + 1) p_(j)).  Rejecting every test whose
    adjusted p-value is at most alpha is the Holm step-down procedure at
    family-wise level alpha.
    """
    items = sorted(pvalues.items(), key=lambda kv: kv[1])
    m = len(items)
    adjusted, running = {}, 0.0
    for i, (name, p) in enumerate(items):
        running = max(running, min(1.0, (m - i) * p))
        adjusted[name] = running
    return adjusted


def _finite_number(name: str, x) -> float:
    try:
        x = float(x)
    except (TypeError, ValueError):
        raise ValueError(f"{name} must be a number; got {x!r}") from None
    if not math.isfinite(x):
        raise ValueError(f"{name} must be finite; got {x}")
    return x


def validate_sample(v,
                    kappa: float,
                    theta_perp: float,
                    theta_par: float,
                    bhat=None,
                    n_attempts: int | None = None,
                    alpha: float = 0.01) -> dict:
    """Run the tests on an (N, 3) velocity sample.  See the module docstring.

    Returns a nested dict of per-test statistics and verdicts, plus a top-level
    ``passed`` flag.  Raises ``ValueError`` on a malformed input contract.
    """
    v = np.asarray(v, dtype=float)
    if v.ndim != 2 or v.shape[1] != 3:
        raise ValueError(f"the sample must have shape (N, 3); got {v.shape}")
    kappa = _finite_number("kappa", kappa)
    if not kappa > 0.5:
        raise ValueError(f"kappa must exceed 1/2; got {kappa}")
    theta_perp = _finite_number("theta_perp", theta_perp)
    theta_par = _finite_number("theta_par", theta_par)
    if not (theta_perp > 0.0 and theta_par > 0.0):
        raise ValueError(f"the thermal speeds must be positive; got theta_perp = "
                         f"{theta_perp}, theta_par = {theta_par}")
    alpha = _finite_number("alpha", alpha)
    if not 0.0 < alpha < 1.0:
        raise ValueError(f"alpha must lie strictly between 0 and 1; got {alpha}")
    if bhat is not None:
        bhat = np.asarray(bhat, dtype=float)
        if bhat.shape != (3,):
            raise ValueError(f"bhat must have three components; got shape {bhat.shape}")
        field_basis(bhat)              # raises on a zero or non-finite vector

    n_total = v.shape[0]
    if n_total < MIN_SAMPLES:
        raise ValueError(f"the tests need at least {MIN_SAMPLES} draws; got {n_total}")
    if n_attempts is not None:
        if isinstance(n_attempts, bool) or int(n_attempts) != n_attempts:
            raise ValueError(f"n_attempts must be an integer; got {n_attempts!r}")
        n_attempts = int(n_attempts)
        if n_attempts < n_total:
            raise ValueError(f"the number of attempts ({n_attempts}) cannot be smaller "
                             f"than the number of draws in the sample ({n_total})")

    report: dict = {
        "inputs": {
            "n_samples": n_total,
            "kappa": float(kappa),
            "theta_perp": float(theta_perp),
            "theta_par": float(theta_par),
            "bhat": None if bhat is None else [float(x) for x in np.asarray(bhat, float)],
            "n_attempts": None if n_attempts is None else int(n_attempts),
            "alpha": float(alpha),
        },
        "tests": {},
        "notes": [],
    }

    # --- rotate into the field-aligned frame -------------------------------
    if bhat is None:
        local = v
        report["inputs"]["frame"] = "assumed field-aligned (third column parallel)"
    else:
        M = field_basis(bhat)
        local = v @ M              # local = M^T v
        report["inputs"]["frame"] = "rotated into the field-aligned frame using bhat"
        report["tests"]["frame_basis"] = {
            "orthonormality_error": float(np.abs(M.T @ M - np.eye(3)).max()),
            "description": "max |M^T M - I| for the reconstructed basis",
        }

    # --- normalized coordinates: the isotropic core ------------------------
    scales = math.sqrt(kappa) * np.array([theta_perp, theta_perp, theta_par])
    u = local / scales

    # Never form T = R^2: it overflows for R > 1.3e154 while R itself is
    # representable, which would silently drop exactly the heavy-tail draws the
    # radial test exists to examine.  np.linalg.norm squares internally; hypot
    # does not.
    R = np.hypot(np.hypot(u[:, 0], u[:, 1]), u[:, 2])
    finite = np.isfinite(R) & np.isfinite(local).all(axis=1)
    n_nonfinite = int((~finite).sum())
    u, R, local_f = u[finite], R[finite], local[finite]
    n = int(R.size)
    if n == 0:
        raise ValueError("every draw is non-finite; nothing to test")

    report["tests"]["finiteness"] = {
        "n_nonfinite": n_nonfinite,
        "fraction": n_nonfinite / n_total,
        "passed": n_nonfinite == 0,
        "description": "non-finite draws, which indicate the loader's precision limit rather than a distributional error",
    }

    b = kappa - 0.5
    pvalues: dict = {}

    # --- 1. radial law -----------------------------------------------------
    # W = 1/(1+T) ~ Beta(kappa-1/2, 3/2), formed without materializing T.  The
    # complementary variable Y = T/(1+T) carries identical information in exact
    # arithmetic but piles its mass against 1.0, where doubles have no relative
    # resolution left; at kappa = 0.55 that alone makes a KS test fail on a
    # perfectly good sample.  The orientation is not cosmetic.
    #
    # For R > 1 the form (1/R)^2 / (1 + (1/R)^2) keeps W resolvable down to the
    # subnormal range, R up to about 6e161, whereas 1/(1 + R^2) or
    # 1/(1 + exp(2 log R)) returns W = 0 already for R > 1.3e154.
    with np.errstate(over="ignore", divide="ignore", under="ignore"):
        small = R <= 1.0
        W = np.empty_like(R)
        W[small] = 1.0 / (1.0 + R[small] ** 2)
        inv2 = (1.0 / R[~small]) ** 2
        W[~small] = inv2 / (1.0 + inv2)

    # The radial tests run on the probability integral transform F_W(W), which
    # is U(0,1) under the target.  Where W underflows to 0 its CDF is still
    # resolvable: for W -> 0, F_W(W) = W^b / (b B(b, 3/2)) (1 + O(W)), and
    # log W = -2 log R to within a relative 1/R^2.  Near kappa = 1/2 these draws
    # carry a measurable share of the probability (about 6e-4 at kappa = 0.51),
    # so mapping them all to F = 0 displaces them and fails an exact sample of a
    # few million draws.
    with np.errstate(divide="ignore"):
        U = stats.beta.cdf(W, b, 1.5)
        under = W <= 0.0
        U[under] = np.exp(-2.0 * b * np.log(R[under])
                          - np.log(b) - special.betaln(b, 1.5))
    n_w_underflow = int(under.sum())

    ks_radial = stats.kstest(U, "uniform")
    cvm_radial = stats.cramervonmises(U, "uniform")
    pvalues["radial_ks"] = float(ks_radial.pvalue)
    pvalues["radial_cvm"] = float(cvm_radial.pvalue)
    report["tests"]["radial_law"] = {
        "statistic_sqrtn_D": float(ks_radial.statistic * math.sqrt(n)),
        "ks_pvalue": float(ks_radial.pvalue),
        "cvm_statistic": float(cvm_radial.statistic),
        "cvm_pvalue": float(cvm_radial.pvalue),
        "n_w_underflow": n_w_underflow,
        "description": ("scaled radius R = |v| in normalized coordinates against "
                        "T = R^2 ~ BetaPrime(3/2, kappa-1/2), tested through the "
                        "CDF of W = 1/(1+T) ~ Beta(kappa-1/2, 3/2); draws whose W "
                        "underflows use the small-W expansion of that CDF"),
    }

    # --- 2. directional uniformity ----------------------------------------
    cos_theta = u[:, 2] / R
    phi = np.arctan2(u[:, 1], u[:, 0])
    ks_cos = stats.kstest(cos_theta, "uniform", args=(-1.0, 2.0))
    ks_phi = stats.kstest(phi, "uniform", args=(-np.pi, 2.0 * np.pi))

    # Equal-solid-angle direction cells, shared with the independence test.
    r_edges = np.quantile(R, np.linspace(0.0, 1.0, N_RADIAL_BINS + 1))
    r_bin = np.clip(np.searchsorted(r_edges[1:-1], R, side="right"), 0, N_RADIAL_BINS - 1)
    c_bin = np.clip(np.floor((cos_theta + 1.0) / 2.0 * N_COS_BINS).astype(int),
                    0, N_COS_BINS - 1)
    p_bin = np.clip(np.floor((phi + np.pi) / (2.0 * np.pi) * N_PHI_BINS).astype(int),
                    0, N_PHI_BINS - 1)
    table = np.zeros((N_RADIAL_BINS, N_COS_BINS * N_PHI_BINS))
    np.add.at(table, (r_bin, c_bin * N_PHI_BINS + p_bin), 1.0)
    # The two KS tests check each angle on its own; uniform marginals make the
    # direction uniform on the sphere only if the angles are also independent.
    # Equal occupancy of the cells, summed over radius, tests the joint distribution.
    cells = stats.chisquare(table.sum(axis=0))

    # The paper's cell-count test: radius shells x direction cells, and the radius
    # alone on the decile shells with the outermost divided into the tail.  Skipped,
    # with a note, when a cell would expect fewer than MIN_EXPECTED_PER_CELL draws,
    # where the chi^2 approximation is poor.  Non-finite draws were removed above
    # and are reported by the finiteness test: a validator cannot tell a component
    # that overflowed at the end of the calculation from one that went infinite
    # through an intermediate underflow, so it does not assign them a radius.
    n_cells = N_SHELLS * N_COS_BINS * N_PHI_BINS
    if n / n_cells >= MIN_EXPECTED_PER_CELL:
        with np.errstate(divide="ignore"):
            log_R = np.log(R)
        shell = np.searchsorted(log_shell_edges(kappa), log_R, side="right")
        cell_counts = np.zeros((N_SHELLS, N_COS_BINS * N_PHI_BINS))
        np.add.at(cell_counts, (shell, c_bin * N_PHI_BINS + p_bin), 1.0)
        cells_all = stats.chisquare(cell_counts.ravel())
        r_log_edges, r_probs = radius_shells(kappa, n)
        r_counts = np.bincount(np.searchsorted(r_log_edges, log_R, side="right"),
                               minlength=r_probs.size)
        cells_radius = stats.chisquare(r_counts, n * r_probs)
        # Pearson's chi^2 against fully specified cell probabilities: k - 1 dof.
        pvalues["cells_radius"] = float(cells_radius.pvalue)
        pvalues["cells_all"] = float(cells_all.pvalue)
        report["tests"]["cells"] = {
            "radius_pvalue": float(cells_radius.pvalue),
            "radius_statistic": float(cells_radius.statistic),
            "radius_dof": int(r_probs.size - 1),
            "radius_shells": int(r_probs.size),
            "radius_outermost_percentile": float(100.0 * (1.0 - r_probs[-1])),
            "radius_expected_outermost": float(n * r_probs[-1]),
            "all_pvalue": float(cells_all.pvalue),
            "all_statistic": float(cells_all.statistic),
            "all_dof": int(n_cells - 1),
            "expected_per_cell": n / n_cells,
            "min_count": int(cell_counts.min()),
            "description": (f"Pearson chi^2 against equal occupancy of {N_SHELLS} "
                            "equal-probability radius shells x "
                            f"{N_COS_BINS * N_PHI_BINS} equal-solid-angle direction "
                            "cells; 'radius' uses the radius shells alone, the "
                            "outermost divided at the 99th to 99.999th percentiles "
                            "of R as far as the sample size allows"),
        }
    else:
        report["notes"].append(
            f"cell-count test skipped: {n} draws give fewer than "
            f"{MIN_EXPECTED_PER_CELL:g} expected per cell over {n_cells} cells")

    pvalues["cos_theta"] = float(ks_cos.pvalue)
    pvalues["phi"] = float(ks_phi.pvalue)
    pvalues["direction_cells"] = float(cells.pvalue)
    report["tests"]["direction"] = {
        "cos_theta_sqrtn_D": float(ks_cos.statistic * math.sqrt(n)),
        "cos_theta_pvalue": float(ks_cos.pvalue),
        "phi_sqrtn_D": float(ks_phi.statistic * math.sqrt(n)),
        "phi_pvalue": float(ks_phi.pvalue),
        "cells_chi2_statistic": float(cells.statistic),
        "cells_chi2_dof": int(N_COS_BINS * N_PHI_BINS - 1),
        "cells_chi2_pvalue": float(cells.pvalue),
        "description": ("cos(Theta) ~ U(-1,1) and Phi ~ U(-pi,pi) about bhat, and "
                        f"equal occupancy of {N_COS_BINS * N_PHI_BINS} equal-solid-angle "
                        f"cells ({N_COS_BINS} cos(Theta) x {N_PHI_BINS} Phi intervals)"),
    }

    # --- 3. radius-direction independence ----------------------------------
    # A direction cell that no draw reaches has no expected count and would stop
    # chi2_contingency; drop it here.  An empty cell already fails the
    # equal-occupancy test above, so nothing is hidden.
    occupied = table.sum(axis=0) > 0
    chi2 = stats.chi2_contingency(table[:, occupied])
    inner = cos_theta[R <= r_edges[1]]
    outer = cos_theta[R > r_edges[-2]]
    ks_io = stats.ks_2samp(inner, outer)
    pvalues["independence_chi2"] = float(chi2.pvalue)
    pvalues["independence_inner_outer"] = float(ks_io.pvalue)
    report["tests"]["independence"] = {
        "chi2_statistic": float(chi2.statistic),
        "chi2_pvalue": float(chi2.pvalue),
        "chi2_dof": int(chi2.dof),
        "min_cell_count": int(table.min()),
        "inner_outer_ks_D": float(ks_io.statistic),
        "inner_outer_ks_pvalue": float(ks_io.pvalue),
        "description": ("chi^2 contingency over radial quartiles x the equal-solid-"
                        "angle direction cells, plus a two-sample KS of cos(Theta) "
                        "between the innermost and outermost radial quartiles"),
    }

    # --- 4. moment-free anisotropy ratio -----------------------------------
    # MAD needs no moments, so this is the anisotropy check that survives
    # 1/2 < kappa <= 3/2 where the variance does not exist.
    mad = stats.median_abs_deviation(local_f, axis=0)
    ratio = float(mad[2] / mad[0])
    expected = float(theta_par / theta_perp)
    perp_ratio = float(mad[1] / mad[0])
    # A ratio of two MADs, each carrying the relative error above; the tolerance
    # is four of those combined standard errors, chosen so that this test does
    # not dominate the false-alarm rate of the family it sits in.  The sqrt(2)
    # treats the two MADs as independent, which they are not: both are driven by
    # the same radial variate, and the correlation grows as the tail heavies.
    # The tolerance is therefore conservative, and increasingly so as
    # kappa -> 1/2 -- at kappa = 0.55 it is about twice the spread actually
    # observed over replicates, and at kappa = 0.51 about ten times it.  This
    # test will not raise a false alarm at small kappa; it will also not catch
    # much there.
    sigma = math.sqrt(2.0) * mad_relative_sigma(kappa, n)
    tol = 4.0 * sigma
    report["tests"]["anisotropy"] = {
        "mad_ratio_par_perp": ratio,
        "expected": expected,
        "relative_error": ratio / expected - 1.0,
        "mad_ratio_perp_perp": perp_ratio,
        "relative_sigma": sigma,
        "tolerance": tol,
        "passed": bool(abs(ratio / expected - 1.0) <= tol and abs(perp_ratio - 1.0) <= tol),
        "description": ("MAD(v_par)/MAD(v_x) against theta_par/theta_perp, and "
                        "MAD(v_y)/MAD(v_x) against 1; the tolerance is four "
                        "asymptotic standard errors and widens as kappa -> 1/2"),
    }

    # --- 5. tail quantiles --------------------------------------------------
    # Compared in log R: with kappa - 1/2 < 1 the upper quantiles span hundreds
    # of decades, so a relative error on R is dominated by tail sampling noise.
    with np.errstate(divide="ignore"):
        theory = 0.5 * np.log(stats.betaprime.ppf(QUANTILE_PROBES, 1.5, b))
        empirical = np.log(np.quantile(R, QUANTILE_PROBES))
    log_err = np.where(np.isfinite(theory), empirical - theory, np.nan)
    report["tests"]["tail_quantiles"] = {
        "probes": list(QUANTILE_PROBES),
        "log_error": [float(x) for x in log_err],
        "ratio": [float(math.exp(x)) if np.isfinite(x) else None for x in log_err],
        "description": ("diagnostic without a verdict: empirical minus exact "
                        "quantiles of log R; ratio is the empirical quantile of R "
                        "divided by the exact one, so 1 means an undistorted tail"),
    }

    # --- 6. an applied velocity bound ---------------------------------------
    max_component = float(np.abs(u).max() * math.sqrt(kappa))
    max_speed_norm = float((R * math.sqrt(kappa)).max())
    bound: dict = {
        "max_normalized_component": max_component,
        "max_normalized_radius": max_speed_norm,
        "description": ("diagnostic without a verdict: largest |v_i|/theta_i and "
                        "largest sqrt(sum (v_i/theta_i)^2) in the sample; a bound "
                        "applied by the loader shows up as a hard edge in whichever "
                        "of these the bounding region uses"),
    }
    if n_attempts is not None:
        rejected = 1.0 - n_total / n_attempts
        bound["rejected_fraction"] = rejected
        bound["total_variation_from_unbounded"] = rejected
        bound["note"] = ("for a bound enforced by discarding and redrawing, the rejected "
                         "fraction equals the total-variation distance from the uncapped "
                         "bi-Kappa distribution exactly; it bounds how far any probability "
                         "can move and does not bound a tail quantile")
        report["notes"].append(
            "A bound was declared. The tests compare the sample with the uncapped "
            "distribution, so the radial tests and the tail quantiles show the effect "
            "of the bound; a failure there measures the bound, not a defect of the "
            "loader.")
    report["tests"]["bound"] = bound

    # --- overall verdict ----------------------------------------------------
    adjusted = _holm_adjusted(pvalues)
    family_tests = {
        name: {"pvalue": p,
               "holm_adjusted_pvalue": adjusted[name],
               "rejected_at_alpha": bool(adjusted[name] <= alpha)}
        for name, p in pvalues.items()
    }
    report["family"] = {
        "alpha": alpha,
        "correction": ("Holm-Bonferroni over the distributional tests; a test is "
                       "rejected when its Holm-adjusted p-value is at most alpha"),
        "n_tests": len(family_tests),
        "tests": family_tests,
        "n_rejected": int(sum(d["rejected_at_alpha"] for d in family_tests.values())),
    }
    report["passed"] = bool(
        report["family"]["n_rejected"] == 0
        and report["tests"]["anisotropy"]["passed"]
        and report["tests"]["finiteness"]["passed"]
    )
    return report


# ----------------------------------------------------------------------------
# reporting and I/O
# ----------------------------------------------------------------------------

def format_report(report: dict) -> str:
    """Render a report as a short human-readable table.

    The distributional tests are listed with their statistic, raw p-value and
    Holm-adjusted p-value; each verdict is FAIL exactly when the adjusted
    p-value is at most alpha.  The tests with their own criteria and the
    diagnostics without a verdict follow.
    """
    i = report["inputs"]
    t = report["tests"]
    fam = report["family"]
    alpha = fam["alpha"]
    out = [
        "bi-Kappa sample validation",
        f"  N = {i['n_samples']}, kappa = {i['kappa']:g}, "
        f"theta_perp = {i['theta_perp']:g}, theta_par = {i['theta_par']:g}",
        f"  frame: {i['frame']}",
        "",
        f"  Distributional tests: {fam['n_tests']} tests, Holm-Bonferroni correction at "
        f"family-wise alpha = {alpha:g}.",
        "  A test fails when its Holm-adjusted p-value is at most alpha.",
        "",
        f"  {'test':36s} {'statistic':>22s} {'p-value':>9s} {'Holm p':>9s}  verdict",
        "  " + "-" * 88,
    ]

    def pfmt(p):
        return f"{p:.3f}" if p >= 1e-3 else f"{p:.1e}"

    def row(name, key, stat):
        d = fam["tests"][key]
        out.append(f"  {name:36s} {stat:>22s} {pfmt(d['pvalue']):>9s} "
                   f"{pfmt(d['holm_adjusted_pvalue']):>9s}  "
                   f"{'FAIL' if d['rejected_at_alpha'] else 'PASS'}")

    if "cells" in t:
        c = t["cells"]
        row(f"cell counts: radius, {c['radius_shells']} shells", "cells_radius",
            f"chi2 = {c['radius_statistic']:.1f} ({c['radius_dof']} dof)")
        row("cell counts: radius x direction", "cells_all",
            f"chi2 = {c['all_statistic']:.1f} ({c['all_dof']} dof)")
    rl = t["radial_law"]
    row("radius, Kolmogorov-Smirnov", "radial_ks",
        f"sqrt(n) D = {rl['statistic_sqrtn_D']:.3f}")
    row("radius, Cramer-von Mises", "radial_cvm", f"T = {rl['cvm_statistic']:.4f}")
    d = t["direction"]
    row("direction: cos(Theta), KS", "cos_theta", f"sqrt(n) D = {d['cos_theta_sqrtn_D']:.3f}")
    row("direction: Phi, KS", "phi", f"sqrt(n) D = {d['phi_sqrtn_D']:.3f}")
    row("direction: equal-solid-angle cells", "direction_cells",
        f"chi2 = {d['cells_chi2_statistic']:.1f} ({d['cells_chi2_dof']} dof)")
    ind = t["independence"]
    row("radius-direction independence", "independence_chi2",
        f"chi2 = {ind['chi2_statistic']:.1f} ({ind['chi2_dof']} dof)")
    row("cos(Theta), inner vs outer quartile", "independence_inner_outer",
        f"D = {ind['inner_outer_ks_D']:.4f}")

    an, fin = t["anisotropy"], t["finiteness"]
    out += [
        "",
        "  Tests with their own criteria:",
        f"  {'anisotropy':36s} {'':>42s}  {'PASS' if an['passed'] else 'FAIL'}",
        f"    MAD(v_par)/MAD(v_perp1)   = {an['mad_ratio_par_perp']:.4f}, expected "
        f"{an['expected']:.4f} (error {100 * an['relative_error']:+.2f}%)",
        f"    MAD(v_perp2)/MAD(v_perp1) = {an['mad_ratio_perp_perp']:.4f}, expected 1",
        f"    each must lie within {100 * an['tolerance']:.2f}% of its expected value "
        "(four standard errors)",
        f"  {'finite draws':36s} {'':>42s}  {'PASS' if fin['passed'] else 'FAIL'}",
        f"    {fin['n_nonfinite']} of {i['n_samples']} draws are non-finite; any is a failure",
        "",
        "  Diagnostics without a verdict:",
        "  tail quantiles of the radius, empirical / exact:",
    ]
    for p, r in zip(t["tail_quantiles"]["probes"], t["tail_quantiles"]["ratio"]):
        out.append(f"    p = {p:<6} {r:.4f}" if r is not None else f"    p = {p:<6} n/a")
    bd = t["bound"]
    out.append(f"  largest |v_i|/theta_i: {bd['max_normalized_component']:.4g}; "
               f"largest normalized speed: {bd['max_normalized_radius']:.4g}")
    if "frame_basis" in t:
        out.append(f"  basis orthonormality error: {t['frame_basis']['orthonormality_error']:.2e}")
    if "rejected_fraction" in bd:
        out.append(f"  declared bound: rejected fraction {bd['rejected_fraction']:.3e}, "
                   "which is the total-variation distance from the uncapped distribution")
    out.append("")
    for note in report["notes"]:
        out.append(f"  note: {note}")
    out.append(f"  overall: {'PASS' if report['passed'] else 'FAIL'}")
    return "\n".join(out)


RAW_EXTENSIONS = (".bin", ".f64")


def _read_text(path: str, delimiter) -> tuple[np.ndarray, str | None]:
    """Read a text table, skipping one non-numeric header line if there is one.

    Blank lines and lines starting with '#' are ignored.  If the first other line
    does not parse as numbers it is taken as a column header and skipped; any
    later non-numeric line is an error.
    """
    header, skip = None, 0
    with open(path, "r") as fh:
        for lineno, raw in enumerate(fh, start=1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            try:
                [float(tok) for tok in line.split(delimiter)]
            except ValueError:
                header, skip = line, lineno
            break
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)        # "input contained no data"
        try:
            v = np.loadtxt(path, delimiter=delimiter, skiprows=skip, ndmin=2)
        except ValueError as exc:
            # NumPy's message may end with advice on its own keyword arguments.
            raise ValueError(f"{path}: {str(exc).split('; use')[0]}") from None
    return v, header


def _load(path: str, binary: bool) -> tuple[np.ndarray, str | None]:
    ext = os.path.splitext(path)[1].lower()
    header = None
    if binary or ext in RAW_EXTENSIONS:
        v = np.fromfile(path, dtype=np.float64)
        if v.size % 3:
            raise ValueError(f"{path}: raw float64 input holds {v.size} values, "
                             "which is not a multiple of 3")
        v = v.reshape(-1, 3)
    elif ext == ".npy":
        v = np.load(path, allow_pickle=False)
    else:
        v, header = _read_text(path, "," if ext == ".csv" else None)
    v = np.asarray(v, dtype=float)
    if v.size == 0:
        raise ValueError(f"{path}: the file holds no data")
    if v.ndim != 2 or v.shape[1] != 3:
        raise ValueError(f"{path}: expected three columns (vx vy vz), got an array "
                         f"of shape {v.shape}")
    return v, header


def load_sample(path: str, binary: bool = False) -> np.ndarray:
    """Read an (N, 3) velocity sample from a file.

    The format follows the extension: ``.bin`` and ``.f64`` are raw float64
    values in the machine's byte order, three per draw (``binary=True`` forces this for any
    extension); ``.npy`` is a NumPy array of shape (N, 3); anything else is text
    with three columns, separated by whitespace, or by commas for ``.csv``.  In
    text, blank lines and lines starting with '#' are ignored, and a first line
    that is not numeric (a header such as ``vx vy vz``) is skipped.

    Raises ``ValueError`` if the file does not hold an (N, 3) array of numbers,
    and ``OSError`` if it cannot be read.
    """
    return _load(path, binary)[0]


EPILOG = """\
input formats:
  .bin, .f64   raw float64 values, three per draw (vx vy vz); --binary forces
               this for any extension
  .npy         NumPy array of shape (N, 3)
  other        text, three columns per line, separated by whitespace (by
               commas for .csv); blank lines and lines starting with '#' are
               ignored, and a non-numeric first line such as 'vx vy vz' is
               skipped as a header

The target is the uncapped bi-Kappa distribution of Eq. (2) of the paper,
  f(v) ~ [1 + v_perp^2/(kappa theta_perp^2) + v_par^2/(kappa theta_par^2)]^-(kappa+1).
A sample from a sampler that rejects draws outside a velocity bound follows a
different, capped distribution; pass --attempts to have the rejected fraction
reported.

The distributional tests are corrected for multiple testing with the
Holm-Bonferroni procedure at family-wise level --alpha: a test fails when its
Holm-adjusted p-value is at most alpha.  The anisotropy test fails beyond four
standard errors, and any non-finite draw is a failure.  Tail quantiles are
reported without a verdict.

exit status:
  0  every test passes
  1  at least one test fails
  2  invalid input or options (message on stderr)

example:
  python bikappa_validate.py samples_bikappa.txt --kappa 2 --theta-perp 1 --theta-par 2
"""


def _parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        prog="bikappa_validate.py",
        description=("Test an (N, 3) velocity sample from any bi-Kappa loader against\n"
                     f"the bi-Kappa distribution.  N must be at least {MIN_SAMPLES}."),
        epilog=EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("sample", help="file holding the sample; see 'input formats' below")
    ap.add_argument("--kappa", type=float, required=True,
                    help="spectral index kappa of the target distribution, > 1/2")
    ap.add_argument("--theta-perp", type=float, required=True, metavar="THETA",
                    help=("perpendicular thermal speed, > 0: a speed in the units of the "
                          "sample, not a temperature (the convention of Eq. (2) of the "
                          "paper and of the C++ class)"))
    ap.add_argument("--theta-par", type=float, required=True, metavar="THETA",
                    help="parallel thermal speed, > 0, in the same convention and units")
    ap.add_argument("--bhat", type=float, nargs=3, default=None, metavar=("BX", "BY", "BZ"),
                    help=("magnetic-field direction in the frame of the sample, need not "
                          "be normalized; omit if the third column is already the "
                          "parallel component"))
    ap.add_argument("--attempts", type=int, default=None, metavar="N",
                    help=("total draws attempted, for a sampler that rejects draws "
                          "outside a velocity bound and redraws them; must be at least "
                          "the number of draws in the sample"))
    ap.add_argument("--alpha", type=float, default=0.01,
                    help="family-wise significance level, 0 < alpha < 1 (default 0.01)")
    ap.add_argument("--binary", action="store_true",
                    help="read the file as raw float64 whatever its extension")
    ap.add_argument("--json", action="store_true", help="print the full report as JSON")
    return ap


def main(argv=None) -> int:
    args = _parser().parse_args(argv)
    prog = "bikappa_validate.py"
    try:
        v, header = _load(args.sample, args.binary)
        if header is not None:
            print(f"{prog}: note: skipped the header line {header!r}", file=sys.stderr)
        report = validate_sample(v, args.kappa, args.theta_perp, args.theta_par,
                                 bhat=args.bhat, n_attempts=args.attempts, alpha=args.alpha)
    except (ValueError, OSError) as exc:
        msg = " ".join(str(exc).split()) or type(exc).__name__
        print(f"{prog}: error: {msg}", file=sys.stderr)
        return 2
    print(json.dumps(report, indent=2) if args.json else format_report(report))
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
