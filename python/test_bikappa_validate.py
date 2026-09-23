#!/usr/bin/env python3
"""Regression tests for bikappa_validate.

A validation battery that never fails is not a validation battery, so the
negative controls below matter more than the positive one: each injects a
specific, plausible loader bug and asserts that the battery names it.

Run with:  uv run --project python python python/test_bikappa_validate.py
Exit status is 0 on success.
"""

from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from bikappa_validate import validate_sample, format_report, field_basis  # noqa: E402

N = 100_000


def draw(n, kappa, theta_perp, theta_par, rng, uniform_theta=False, bhat=None):
    """A reference bi-Kappa sample, independent of the C++ implementation."""
    # sqrt(X1)/sqrt(X2), never sqrt(X1/X2): the quotient overflows at small
    # kappa in a regime where the radius itself is perfectly representable.
    R = np.sqrt(rng.gamma(1.5, 1.0, n)) / np.sqrt(rng.gamma(kappa - 0.5, 1.0, n))
    if uniform_theta:                        # the classic direction bug
        c = np.cos(rng.uniform(0.0, np.pi, n))
    else:
        c = rng.uniform(-1.0, 1.0, n)
    phi = rng.uniform(0.0, 2.0 * np.pi, n)
    s = np.sqrt(1.0 - c * c)
    u = np.column_stack([R * s * np.cos(phi), R * s * np.sin(phi), R * c])
    v = u * (np.sqrt(kappa) * np.array([theta_perp, theta_perp, theta_par]))
    if bhat is not None:
        v = v @ field_basis(bhat).T
    return v


def draw_log_radius(n, kappa, theta_perp, theta_par, rng):
    """An exact reference sample at kappa near 1/2, with the radius formed in logs.

    X2 ~ Gamma(kappa - 1/2) is drawn as Gamma(kappa + 1/2) * U^(1/(kappa - 1/2)),
    so log X2 is available even where X2 itself underflows.  The radius then stays
    finite up to exp(709), far past the point where W = 1/(1+R^2) underflows.
    """
    b = kappa - 0.5
    log_x1 = np.log(rng.gamma(1.5, 1.0, n))
    log_x2 = np.log(rng.gamma(b + 1.0, 1.0, n)) + np.log(rng.uniform(size=n)) / b
    with np.errstate(over="ignore"):
        R = np.exp(0.5 * (log_x1 - log_x2))
    c = rng.uniform(-1.0, 1.0, n)
    phi = rng.uniform(0.0, 2.0 * np.pi, n)
    s = np.sqrt(1.0 - c * c)
    u = np.column_stack([R * s * np.cos(phi), R * s * np.sin(phi), R * c])
    return u * (np.sqrt(kappa) * np.array([theta_perp, theta_perp, theta_par]))


def failing_tests(report):
    bad = [t for t, d in report["family"]["tests"].items() if d["rejected_at_alpha"]]
    if not report["tests"]["anisotropy"]["passed"]:
        bad.append("anisotropy")
    if not report["tests"]["finiteness"]["passed"]:
        bad.append("finiteness")
    return bad


def main() -> int:
    rng = np.random.default_rng(20260820)
    failures = []

    def check(name, ok, detail=""):
        print(f"  [{'ok' if ok else 'FAILED'}] {name}{'  ' + detail if detail else ''}")
        if not ok:
            failures.append(name)

    print("positive control")
    for kappa in (0.55, 0.75, 1.0, 1.5, 2.0, 5.0, 10.0):
        r = validate_sample(draw(N, kappa, 1.0, 2.0, rng), kappa, 1.0, 2.0)
        check(f"correct sample passes at kappa = {kappa}", r["passed"],
              f"failing: {failing_tests(r)}")

    # At kappa = 0.51 the reference generator above underflows in double
    # precision -- it draws a Gamma variate of shape 0.01 -- so the battery
    # should report non-finite draws rather than a distributional failure.  The
    # fraction it sees is the same order as the released C++ loader's.
    r = validate_sample(draw(N, 0.51, 1.0, 2.0, rng), 0.51, 1.0, 2.0)
    check("non-finite draws are reported, not misread as a bad law",
          failing_tests(r) == ["finiteness"],
          f"fraction = {r['tests']['finiteness']['fraction']:.1e}")

    # With the radius formed in logs the sample is exact, and about 6e-4 of its
    # probability lies where W = 1/(1+R^2) underflows to 0.  Release 2.2.0 mapped
    # those draws to F_W = 0 and failed an exact sample of this size
    # (sqrt(n) D = 2.33 at n = 8e6).
    n_big = 8_000_000
    v = draw_log_radius(n_big, 0.51, 1.0, 2.0, rng)
    finite = np.isfinite(v).all(axis=1)
    r = validate_sample(v[finite], 0.51, 1.0, 2.0)
    check("exact kappa = 0.51 sample passes where W underflows", r["passed"],
          f"sqrt(n) D = {r['tests']['radial_law']['statistic_sqrtn_D']:.3f}, "
          f"W underflows: {r['tests']['radial_law']['n_w_underflow']}, "
          f"failing: {failing_tests(r)}")
    check("the underflow branch is exercised",
          r["tests"]["radial_law"]["n_w_underflow"] > 1000)

    print("positive control, rotated frame")
    bhat = np.array([0.3, -0.5, 0.8])
    bhat = bhat / np.linalg.norm(bhat)
    r = validate_sample(draw(N, 2.0, 1.0, 2.0, rng, bhat=bhat), 2.0, 1.0, 2.0, bhat=bhat)
    check("rotated sample passes when bhat is supplied", r["passed"], f"failing: {failing_tests(r)}")
    check("basis orthonormality is at round-off",
          r["tests"]["frame_basis"]["orthonormality_error"] < 1e-14)

    r = validate_sample(draw(N, 2.0, 1.0, 2.0, rng, bhat=bhat), 2.0, 1.0, 2.0, bhat=None)
    check("rotated sample fails when bhat is withheld", not r["passed"],
          f"caught by: {failing_tests(r)}")

    print("negative controls")
    r = validate_sample(draw(N, 2.0, 1.0, 2.0, rng, uniform_theta=True), 2.0, 1.0, 2.0)
    check("uniform-Theta direction bug is caught",
          "cos_theta" in failing_tests(r) and "cells_all" in failing_tests(r),
          f"caught by: {failing_tests(r)}")

    r = validate_sample(draw(N, 3.0, 1.0, 2.0, rng), 2.0, 1.0, 2.0)
    check("wrong kappa is caught by the radial test",
          "radial_ks" in failing_tests(r) and "cells_radius" in failing_tests(r),
          f"caught by: {failing_tests(r)}")

    r = validate_sample(draw(N, 2.0, 1.0, 2.0, rng), 2.0, 2.0, 1.0)
    check("swapped thermal speeds are caught", "anisotropy" in failing_tests(r),
          f"caught by: {failing_tests(r)}")

    # A radius-direction coupling that leaves both marginals correct: rotate the
    # direction of the outer half of the sample onto a preferred axis by pairing
    # sorted radii with sorted |cos Theta|.  Marginals are untouched by
    # construction; only the joint law is wrong.
    v = draw(N, 2.0, 1.0, 1.0, rng)
    R = np.hypot(np.hypot(v[:, 0], v[:, 1]), v[:, 2])
    order_r = np.argsort(R)
    c = v[:, 2] / R
    order_c = np.argsort(np.abs(c))
    v2 = np.empty_like(v)
    v2[order_r] = (v[order_c] / np.linalg.norm(v[order_c], axis=1)[:, None]) * R[order_r, None]
    r = validate_sample(v2, 2.0, 1.0, 1.0)
    check("radius-direction coupling is caught with marginals intact",
          "independence_chi2" in failing_tests(r) or "independence_inner_outer" in failing_tests(r),
          f"caught by: {failing_tests(r)}")

    # The same kind of coupling in the azimuth: pair sorted radii with sorted Phi.
    # cos(Theta) is untouched and independent of the radius, so a test binned in
    # cos(Theta) alone cannot see this.
    v = draw(N, 2.0, 1.0, 1.0, rng)
    R = np.hypot(np.hypot(v[:, 0], v[:, 1]), v[:, 2])
    c = v[:, 2] / R
    s = np.sqrt(1.0 - c * c)
    phi = np.empty(N)
    phi[np.argsort(R)] = np.sort(rng.uniform(-np.pi, np.pi, N))
    v3 = R[:, None] * np.column_stack([s * np.cos(phi), s * np.sin(phi), c])
    r = validate_sample(v3, 2.0, 1.0, 1.0)
    check("radius-azimuth coupling is caught with cos(Theta) independent of R",
          "independence_chi2" in failing_tests(r) and "cells_all" in failing_tests(r),
          f"caught by: {failing_tests(r)}")

    # Both angles uniform on their own but tied to each other, Phi = pi cos(Theta):
    # the two KS tests pass, and only the joint direction test can fail.
    v = draw(N, 2.0, 1.0, 1.0, rng)
    R = np.hypot(np.hypot(v[:, 0], v[:, 1]), v[:, 2])
    c = rng.uniform(-1.0, 1.0, N)
    phi = np.pi * c
    s = np.sqrt(1.0 - c * c)
    v4 = R[:, None] * np.column_stack([s * np.cos(phi), s * np.sin(phi), c])
    r = validate_sample(v4, 2.0, 1.0, 1.0)
    bad = failing_tests(r)
    check("angle-angle coupling is caught with both angle marginals uniform",
          "direction_cells" in bad and "cells_all" in bad
          and "cos_theta" not in bad and "phi" not in bad,
          f"caught by: {bad}")

    print("bound accounting")

    lam = 5.0
    v = draw(4 * N, 2.0, 1.0, 2.0, rng)
    keep = (np.abs(v[:, 0]) <= lam) & (np.abs(v[:, 1]) <= lam) & (np.abs(v[:, 2]) <= 2.0 * lam)
    vc, attempts = v[keep], v.shape[0]
    r = validate_sample(vc, 2.0, 1.0, 2.0, n_attempts=attempts)
    rf = r["tests"]["bound"]["rejected_fraction"]
    check("rejected fraction is reported and equals the TV distance",
          abs(rf - (1 - vc.shape[0] / attempts)) < 1e-12
          and r["tests"]["bound"]["total_variation_from_unbounded"] == rf,
          f"rejected fraction = {rf:.3e}")
    check("the capped tail is reported as shortened",
          r["tests"]["tail_quantiles"]["ratio"][-1] < 1.0,
          f"p99.9 ratio = {r['tests']['tail_quantiles']['ratio'][-1]:.4f}")

    print("input contract")
    for bad_input, why in (
        (dict(v=np.zeros((10, 3)), kappa=2.0, theta_perp=1.0, theta_par=1.0), "too few draws"),
        (dict(v=np.zeros((N, 2)), kappa=2.0, theta_perp=1.0, theta_par=1.0), "wrong shape"),
        (dict(v=np.zeros((N, 3)), kappa=0.4, theta_perp=1.0, theta_par=1.0), "kappa <= 1/2"),
        (dict(v=np.zeros((N, 3)), kappa=2.0, theta_perp=0.0, theta_par=1.0), "zero thermal speed"),
    ):
        try:
            validate_sample(**bad_input)
            check(f"rejects {why}", False)
        except ValueError:
            check(f"rejects {why}", True)

    print("report rendering")
    r = validate_sample(draw(N, 2.0, 1.0, 2.0, rng), 2.0, 1.0, 2.0)
    text = format_report(r)
    check("report renders", "overall: PASS" in text and "radial law" in text)

    print()
    if failures:
        print(f"FAILED: {len(failures)} check(s): {failures}")
        return 1
    print("all checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
