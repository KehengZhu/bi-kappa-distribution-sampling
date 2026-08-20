#!/usr/bin/env python3
"""Distributional validation of a bi-Kappa velocity sample.

This module tests a *sample*, not a sampler.  Nothing in it assumes the
Gamma-ratio construction, so it applies to velocities produced by any bi-Kappa
loader -- inverse transform, acceptance-rejection, or scale mixture.

Input contract
--------------
``v``            (N, 3) array of velocities, in the simulation's global frame.
``kappa``        shape index, > 1/2.
``theta_perp``   perpendicular thermal speed, > 0, in the units of ``v``.
``theta_par``    parallel thermal speed, > 0, in the units of ``v``.
``bhat``         magnetic-field direction, or ``None`` if the sample is already
                 field-aligned (i.e. the third column is the parallel one).
``n_attempts``   total draws attempted, if the loader enforced a velocity bound
                 by discarding and redrawing.  Supplying it lets the report
                 quote the rejected fraction, which for a hard bound equals the
                 total-variation distance from the unbounded law.

The thermal speeds follow the convention of Eq. (2) of the accompanying paper,

    f(v) ~ [1 + v_perp^2/(kappa theta_perp^2) + v_par^2/(kappa theta_par^2)]^-(kappa+1),

normalizable for every kappa > 1/2.  A temperature-linked convention such as
theta^2 = 2[(kappa-3/2)/kappa](k_B T/m) is defined only for kappa > 3/2 and is
*not* interchangeable with this one; convert before calling.

What the tests establish, and what they do not
----------------------------------------------
The radial, directional and independence tests together pin down the
three-dimensional law, because a spherically symmetric density is determined by
its radial law plus a uniform, radius-independent direction.  They are
distribution-free and use no moments, so they remain valid on
1/2 < kappa <= 3/2 where the second moment does not exist.  A pass is evidence
that the sample is consistent with the target at the given size; it is not a
proof of correctness, and it says nothing about a bound the loader may have
applied -- for that, supply ``n_attempts`` and read the tail-quantile block.

Usage
-----
    from bikappa_validate import validate_sample, format_report
    report = validate_sample(v, kappa=2.0, theta_perp=1.0, theta_par=2.0)
    print(format_report(report))

or from the command line:

    python bikappa_validate.py sample.npy --kappa 2 --theta-perp 1 --theta-par 2

Exit status is 0 if every test passes and 1 otherwise.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys

import numpy as np
from scipy import stats

__all__ = ["validate_sample", "format_report", "load_sample"]

# Quantile probes for the radial tail.  The two extreme ones are where a radial
# error shows up first and where Cartesian marginals are least sensitive.
QUANTILE_PROBES = (0.5, 0.9, 0.99, 0.999)

N_RADIAL_BINS = 4
N_COS_BINS = 10


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
# the battery
# ----------------------------------------------------------------------------

def _ks_critical(alpha: float) -> float:
    """Asymptotic critical value of sqrt(n) D for a two-sided KS test."""
    return math.sqrt(-0.5 * math.log(alpha / 2.0))


def mad_relative_sigma(kappa: float, n: int) -> float:
    """Asymptotic relative standard error of a sample MAD, per component.

    Each Cartesian marginal of Eq. (2) is a scaled Student-t law with
    nu = 2 kappa - 1 degrees of freedom, so in units of its own scale the MAD is
    t0 = ppf(3/4, nu).  For a symmetric law the MAD is the median of |X|, whose
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


def _holm(pvalues: dict, alpha: float) -> dict:
    """Holm-Bonferroni step-down over a family of p-values."""
    items = sorted(pvalues.items(), key=lambda kv: kv[1])
    m = len(items)
    rejected, still = {}, True
    for i, (name, p) in enumerate(items):
        thresh = alpha / (m - i)
        if still and p <= thresh:
            rejected[name] = True
        else:
            still = False
            rejected[name] = False
    return rejected


def validate_sample(v,
                    kappa: float,
                    theta_perp: float,
                    theta_par: float,
                    bhat=None,
                    n_attempts: int | None = None,
                    alpha: float = 0.01) -> dict:
    """Run the battery on an (N, 3) velocity sample.  See the module docstring.

    Returns a nested dict of per-test statistics and verdicts, plus a top-level
    ``passed`` flag.  Raises ``ValueError`` on a malformed input contract.
    """
    v = np.asarray(v, dtype=float)
    if v.ndim != 2 or v.shape[1] != 3:
        raise ValueError(f"v must have shape (N, 3); got {v.shape}")
    if not kappa > 0.5:
        raise ValueError(f"kappa must exceed 1/2; got {kappa}")
    if not (theta_perp > 0 and theta_par > 0):
        raise ValueError("thermal speeds must be positive")

    n_total = v.shape[0]
    if n_total < 1000:
        raise ValueError(f"need at least 1000 draws for the tests to mean anything; got {n_total}")

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
    ks_crit = _ks_critical(alpha)
    pvalues: dict = {}

    # --- 1. radial law -----------------------------------------------------
    # W = 1/(1+T) ~ Beta(kappa-1/2, 3/2), computed from log R so that T is never
    # materialized.  The complementary variable Y = T/(1+T) carries identical
    # information in exact arithmetic but piles its mass against 1.0, where
    # doubles have no relative resolution left; at kappa = 0.55 that alone makes
    # a KS test fail on a perfectly good sample.  The orientation is not
    # cosmetic.
    with np.errstate(over="ignore", divide="ignore", under="ignore"):
        logR = np.log(R)
        W = 1.0 / (1.0 + np.exp(2.0 * logR))       # = expit(-2 log R)
        W = np.where(R == 0.0, 1.0, W)

    ks_radial = stats.kstest(W, "beta", args=(b, 1.5))
    cvm_radial = stats.cramervonmises(W, "beta", args=(b, 1.5))
    pvalues["radial_ks"] = float(ks_radial.pvalue)
    pvalues["radial_cvm"] = float(cvm_radial.pvalue)
    report["tests"]["radial_law"] = {
        "statistic_sqrtn_D": float(ks_radial.statistic * math.sqrt(n)),
        "critical_sqrtn_D": ks_crit,
        "ks_pvalue": float(ks_radial.pvalue),
        "cvm_statistic": float(cvm_radial.statistic),
        "cvm_pvalue": float(cvm_radial.pvalue),
        "description": ("scaled radius R = |v| in normalized coordinates against "
                        "T = R^2 ~ BetaPrime(3/2, kappa-1/2), tested through "
                        "W = 1/(1+T) ~ Beta(kappa-1/2, 3/2)"),
    }

    # --- 2. directional uniformity ----------------------------------------
    cos_theta = u[:, 2] / R
    phi = np.arctan2(u[:, 1], u[:, 0])
    ks_cos = stats.kstest(cos_theta, "uniform", args=(-1.0, 2.0))
    ks_phi = stats.kstest(phi, "uniform", args=(-np.pi, 2.0 * np.pi))
    pvalues["cos_theta"] = float(ks_cos.pvalue)
    pvalues["phi"] = float(ks_phi.pvalue)
    report["tests"]["direction"] = {
        "cos_theta_sqrtn_D": float(ks_cos.statistic * math.sqrt(n)),
        "cos_theta_pvalue": float(ks_cos.pvalue),
        "phi_sqrtn_D": float(ks_phi.statistic * math.sqrt(n)),
        "phi_pvalue": float(ks_phi.pvalue),
        "critical_sqrtn_D": ks_crit,
        "description": "cos(Theta) ~ U(-1,1) and Phi ~ U(-pi,pi) about bhat",
    }

    # --- 3. radius-direction independence ----------------------------------
    r_edges = np.quantile(R, np.linspace(0.0, 1.0, N_RADIAL_BINS + 1))
    r_edges[0], r_edges[-1] = -np.inf, np.inf
    c_edges = np.linspace(-1.0, 1.0, N_COS_BINS + 1)
    c_edges[0], c_edges[-1] = -np.inf, np.inf
    table = np.histogram2d(R, cos_theta, bins=[r_edges, c_edges])[0]
    chi2 = stats.chi2_contingency(table)
    inner = cos_theta[R <= r_edges[1]]
    outer = cos_theta[R > r_edges[-2]]
    ks_io = stats.ks_2samp(inner, outer)
    pvalues["independence_chi2"] = float(chi2.pvalue)
    pvalues["independence_inner_outer"] = float(ks_io.pvalue)
    report["tests"]["independence"] = {
        "chi2_pvalue": float(chi2.pvalue),
        "inner_outer_ks_pvalue": float(ks_io.pvalue),
        "description": ("chi^2 contingency over radial quartiles x direction-cosine "
                        "deciles, plus a two-sample KS of cos(Theta) between the "
                        "innermost and outermost radial quartiles"),
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
        "description": ("empirical minus exact quantiles of log R; ratio is the "
                        "empirical quantile of R divided by the exact one, so 1 "
                        "means an undistorted tail"),
    }

    # --- 6. an applied velocity bound ---------------------------------------
    max_component = float(np.abs(u).max() * math.sqrt(kappa))
    max_speed_norm = float((R * math.sqrt(kappa)).max())
    bound: dict = {
        "max_normalized_component": max_component,
        "max_normalized_radius": max_speed_norm,
        "description": ("largest |v_i|/theta_i and largest sqrt(sum (v_i/theta_i)^2) "
                        "in the sample; a bound applied by the loader shows up as a "
                        "hard edge in whichever of these the bounding region uses"),
    }
    if n_attempts is not None:
        if n_attempts < n_total:
            raise ValueError("n_attempts cannot be smaller than the number of draws returned")
        rejected = 1.0 - n_total / n_attempts
        bound["rejected_fraction"] = rejected
        bound["total_variation_from_unbounded"] = rejected
        bound["note"] = ("for a bound enforced by discarding and redrawing, the rejected "
                         "fraction equals the total-variation distance from the unbounded "
                         "bi-Kappa law exactly; it bounds how far any probability can move "
                         "and does not bound a tail quantile")
        report["notes"].append(
            "A bound was declared. The radial and tail tests above compare against the "
            "*unbounded* law and are expected to fail in the tail by construction; read "
            "them as a measurement of the bound's effect, not as a defect.")
    report["tests"]["bound"] = bound

    # --- overall verdict ----------------------------------------------------
    rejected_flags = _holm(pvalues, alpha)
    for name, rej in rejected_flags.items():
        pvalues[name] = {"pvalue": pvalues[name], "rejected_at_alpha": rej}
    report["family"] = {
        "alpha": alpha,
        "correction": "Holm-Bonferroni over the distributional tests",
        "tests": pvalues,
        "n_rejected": int(sum(rejected_flags.values())),
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
    """Render a report as a short human-readable table."""
    i = report["inputs"]
    out = [
        "bi-Kappa sample validation",
        f"  N = {i['n_samples']}, kappa = {i['kappa']}, "
        f"theta_perp = {i['theta_perp']}, theta_par = {i['theta_par']}",
        f"  frame: {i['frame']}",
        "",
    ]
    t = report["tests"]
    crit = t["radial_law"]["critical_sqrtn_D"]
    out.append(f"  {'test':32s} {'statistic':>14s}  verdict")
    out.append("  " + "-" * 62)

    def line(name, stat, ok):
        out.append(f"  {name:32s} {stat:>14s}  {'PASS' if ok else 'FAIL'}")

    fam = report["family"]["tests"]
    line("radial law (sqrt(n) D)",
         f"{t['radial_law']['statistic_sqrtn_D']:.3f}", not fam["radial_ks"]["rejected_at_alpha"])
    line("radial law (Cramer-von Mises)",
         f"{t['radial_law']['cvm_statistic']:.4f}", not fam["radial_cvm"]["rejected_at_alpha"])
    line("direction cos(Theta)",
         f"{t['direction']['cos_theta_sqrtn_D']:.3f}", not fam["cos_theta"]["rejected_at_alpha"])
    line("direction Phi",
         f"{t['direction']['phi_sqrtn_D']:.3f}", not fam["phi"]["rejected_at_alpha"])
    line("radius-direction independence",
         f"p={t['independence']['chi2_pvalue']:.3f}",
         not fam["independence_chi2"]["rejected_at_alpha"])
    line("  inner vs outer quartile",
         f"p={t['independence']['inner_outer_ks_pvalue']:.3f}",
         not fam["independence_inner_outer"]["rejected_at_alpha"])
    line("anisotropy MAD ratio",
         f"{t['anisotropy']['mad_ratio_par_perp']:.4f}", t["anisotropy"]["passed"])
    line("finite draws",
         f"{t['finiteness']['n_nonfinite']} bad", t["finiteness"]["passed"])
    if "frame_basis" in t:
        out.append(f"  {'basis orthonormality':32s} "
                   f"{t['frame_basis']['orthonormality_error']:>14.2e}")
    out.append("")
    out.append(f"  critical sqrt(n) D at alpha = {report['inputs']['alpha']}: {crit:.3f}"
               "   (about 0.87 on average when the law is correct)")
    out.append("")
    out.append("  tail quantiles of the radius, empirical / exact:")
    for p, r in zip(t["tail_quantiles"]["probes"], t["tail_quantiles"]["ratio"]):
        out.append(f"    p = {p:<6} {r:.4f}" if r is not None else f"    p = {p:<6} n/a")
    bd = t["bound"]
    if "rejected_fraction" in bd:
        out.append("")
        out.append(f"  declared bound: rejected fraction {bd['rejected_fraction']:.3e}, "
                   f"which is the total-variation distance from the unbounded law")
    out.append("")
    for note in report["notes"]:
        out.append(f"  note: {note}")
    out.append(f"  overall: {'PASS' if report['passed'] else 'FAIL'}")
    return "\n".join(out)


def load_sample(path: str, binary: bool = False) -> np.ndarray:
    """Read an (N, 3) sample from .npy, whitespace/comma text, or raw float64."""
    ext = os.path.splitext(path)[1].lower()
    if ext == ".npy":
        v = np.load(path)
    elif binary or ext in (".bin", ".dat", ".f64"):
        v = np.fromfile(path, dtype=np.float64)
    else:
        v = np.loadtxt(path, delimiter="," if ext == ".csv" else None)
    v = np.asarray(v, dtype=float)
    if v.ndim == 1:
        if v.size % 3:
            raise ValueError("flat input length is not a multiple of 3")
        v = v.reshape(-1, 3)
    return v


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description="Validate an (N, 3) bi-Kappa velocity sample from any loader.")
    ap.add_argument("sample", help="path to .npy, text (.txt/.csv), or raw float64 (.bin)")
    ap.add_argument("--kappa", type=float, required=True)
    ap.add_argument("--theta-perp", type=float, required=True)
    ap.add_argument("--theta-par", type=float, required=True)
    ap.add_argument("--bhat", type=float, nargs=3, default=None,
                    help="magnetic-field direction; omit if the sample is field-aligned")
    ap.add_argument("--attempts", type=int, default=None,
                    help="total draws attempted, if a velocity bound was enforced by redrawing")
    ap.add_argument("--alpha", type=float, default=0.01)
    ap.add_argument("--binary", action="store_true", help="force raw float64 input")
    ap.add_argument("--json", action="store_true", help="emit the full report as JSON")
    args = ap.parse_args(argv)

    v = load_sample(args.sample, binary=args.binary)
    report = validate_sample(v, args.kappa, args.theta_perp, args.theta_par,
                             bhat=args.bhat, n_attempts=args.attempts, alpha=args.alpha)
    print(json.dumps(report, indent=2) if args.json else format_report(report))
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
