#!/usr/bin/env python3
"""Regression tests for bikappa_validate.

A set of tests that never fails tests nothing, so the
negative controls below matter more than the positive one: each injects a
specific, plausible loader bug and asserts that the tests catch it.

Run with:  uv run --project python python python/test_bikappa_validate.py
Exit status is 0 on success.  The same checks run under pytest, one test per
section:  uv run --project python --with pytest pytest python/test_bikappa_validate.py
"""

from __future__ import annotations

import contextlib
import io
import os
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import bikappa_validate  # noqa: E402
from bikappa_validate import (validate_sample, format_report, field_basis,  # noqa: E402
                              log_radius_quantile, MIN_SAMPLES)

N = 100_000
SEED = 20260820


def draw(n, kappa, theta_perp, theta_par, rng, uniform_theta=False, bhat=None):
    """A reference bi-Kappa sample, independent of the C++ implementation."""
    # sqrt(X1)/sqrt(X2), never sqrt(X1/X2): the quotient overflows at small
    # kappa in a regime where the radius itself is perfectly representable.
    # Near kappa = 1/2, X2 can underflow to 0; the division then gives R = inf,
    # which the validator is expected to report as a non-finite draw.
    with np.errstate(divide="ignore"):
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


class Checks:
    """Collects named checks, printing each as it runs."""

    def __init__(self):
        self.failures = []

    def __call__(self, name, ok, detail=""):
        print(f"  [{'ok' if ok else 'FAILED'}] {name}{'  ' + detail if detail else ''}")
        if not ok:
            self.failures.append(name)


# ----------------------------------------------------------------------------
# sections: each takes a random generator and a Checks object
# ----------------------------------------------------------------------------

def positive_control(rng, check):
    print("positive control")
    for kappa in (0.55, 0.75, 1.0, 1.5, 2.0, 5.0, 10.0):
        r = validate_sample(draw(N, kappa, 1.0, 2.0, rng), kappa, 1.0, 2.0)
        check(f"correct sample passes at kappa = {kappa}", r["passed"],
              f"failing: {failing_tests(r)}")

    # At kappa = 0.51 the reference generator above underflows in double
    # precision -- it draws a Gamma variate of shape 0.01 -- and the draws it
    # loses are those with R above about 10^161, beyond the 99.9th percentile of R.
    # The validator should report the non-finite draws, and the radius shells beyond
    # the 99.9th percentile should see that the finite draws lack that tail; the
    # remaining tests, which do not resolve the fastest 0.1%, should pass.
    r = validate_sample(draw(N, 0.51, 1.0, 2.0, rng), 0.51, 1.0, 2.0)
    check("underflow losses are reported and the depleted tail is caught",
          sorted(failing_tests(r)) == ["cells_radius", "finiteness"],
          f"fraction = {r['tests']['finiteness']['fraction']:.1e}, "
          f"caught by: {failing_tests(r)}")

    # With the radius formed in logs the sample is exact, and about 6e-4 of its
    # probability lies where W = 1/(1+R^2) underflows to 0.  Mapping those draws
    # to F_W = 0 would fail an exact sample of this size (sqrt(n) D = 2.33 at
    # n = 8e6); the validator evaluates F_W from log R there instead.
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


def rotated_frame(rng, check):
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


def negative_controls(rng, check):
    print("negative controls")
    r = validate_sample(draw(N, 2.0, 1.0, 2.0, rng, uniform_theta=True), 2.0, 1.0, 2.0)
    check("uniform-Theta direction bug is caught",
          "cos_theta" in failing_tests(r) and "cells_all" in failing_tests(r),
          f"caught by: {failing_tests(r)}")

    r = validate_sample(draw(N, 3.0, 1.0, 2.0, rng), 2.0, 1.0, 2.0)
    check("wrong kappa is caught by the radial test",
          "radial_ks" in failing_tests(r) and "cells_radius" in failing_tests(r),
          f"caught by: {failing_tests(r)}")

    # A tail error confined to the fastest 1%: beyond the 99th percentile r_99 the
    # radius is compressed, R -> r_99 (R/r_99)^(1/2).  The draws stay beyond r_99,
    # so every decile-shell count, and with it the 400-cell test, is unchanged; the
    # CDF moves by at most 10^-2, below what a KS test resolves at this size.  Only
    # the radius shells beyond r_99 see it.
    kappa = 2.0
    v = draw(N, kappa, 1.0, 2.0, rng)
    u = v / (np.sqrt(kappa) * np.array([1.0, 1.0, 2.0]))
    R = np.hypot(np.hypot(u[:, 0], u[:, 1]), u[:, 2])
    r99 = np.exp(log_radius_quantile(kappa, 1e-2)[0])
    far = R > r99
    v[far] *= ((r99 * np.sqrt(R[far] / r99)) / R[far])[:, None]
    r = validate_sample(v, kappa, 1.0, 2.0)
    check("tail compression beyond the 99th percentile is caught by the radius shells",
          "cells_radius" in failing_tests(r) and "cells_all" not in failing_tests(r),
          f"caught by: {failing_tests(r)}")

    r = validate_sample(draw(N, 2.0, 1.0, 2.0, rng), 2.0, 2.0, 1.0)
    check("swapped thermal speeds are caught", "anisotropy" in failing_tests(r),
          f"caught by: {failing_tests(r)}")

    # A radius-direction coupling that leaves both marginals correct: rotate the
    # direction of the outer half of the sample onto a preferred axis by pairing
    # sorted radii with sorted |cos Theta|.  Marginals are untouched by
    # construction; only the joint distribution is wrong.
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


def bound_accounting(rng, check):
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


def input_contract(rng, check):
    print("input contract")
    ok = dict(kappa=2.0, theta_perp=1.0, theta_par=1.0)
    v = np.zeros((N, 3))
    for bad_input, why in (
        (dict(ok, v=np.zeros((10, 3))), "too few draws"),
        (dict(ok, v=np.zeros((N, 2))), "wrong shape"),
        (dict(ok, v=v, kappa=0.4), "kappa <= 1/2"),
        (dict(ok, v=v, kappa=float("nan")), "kappa = nan"),
        (dict(ok, v=v, kappa=float("inf")), "kappa = inf"),
        (dict(ok, v=v, theta_perp=0.0), "zero thermal speed"),
        (dict(ok, v=v, theta_par=float("inf")), "infinite thermal speed"),
        (dict(ok, v=v, alpha=5.0), "alpha >= 1"),
        (dict(ok, v=v, alpha=0.0), "alpha = 0"),
        (dict(ok, v=v, bhat=[0.0, 0.0, 0.0]), "zero bhat"),
        (dict(ok, v=v, n_attempts=N - 1), "fewer attempts than draws"),
    ):
        try:
            validate_sample(**bad_input)
            check(f"rejects {why}", False)
        except ValueError:
            check(f"rejects {why}", True)


def report_rendering(rng, check):
    print("report rendering")
    r = validate_sample(draw(N, 2.0, 1.0, 2.0, rng), 2.0, 1.0, 2.0)
    text = format_report(r)
    check("report renders", "overall: PASS" in text and "Holm" in text
          and "Kolmogorov-Smirnov" in text)

    # Every family verdict must follow from the printed Holm-adjusted p-value, and
    # the adjusted p-values must reproduce the Holm step-down procedure.
    fam = r["family"]
    check("each verdict is 'Holm-adjusted p <= alpha'",
          all(d["rejected_at_alpha"] == (d["holm_adjusted_pvalue"] <= fam["alpha"])
              for d in fam["tests"].values()))
    p = {f"t{k}": x for k, x in enumerate(rng.uniform(0, 0.02, 12) ** 2)}
    adj = bikappa_validate._holm_adjusted(p)
    for alpha in (1e-4, 1e-3, 1e-2):
        stepdown, still = {}, True
        for i, (name, x) in enumerate(sorted(p.items(), key=lambda kv: kv[1])):
            still = still and x <= alpha / (len(p) - i)
            stepdown[name] = still
        check(f"adjusted p-values reproduce the step-down at alpha = {alpha:g}",
              all(stepdown[k] == (adj[k] <= alpha) for k in p))


def command_line(rng, check):
    print("command line")
    n = 20_000
    v = draw(n, 2.0, 1.0, 2.0, rng)
    params = ["--kappa", "2", "--theta-perp", "1", "--theta-par", "2"]

    def run(argv):
        out, err = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
            try:
                code = bikappa_validate.main(argv)
            except SystemExit as exc:          # argparse usage errors
                code = exc.code
        return code, out.getvalue(), err.getvalue()

    with tempfile.TemporaryDirectory() as tmp:
        def path(name):
            return os.path.join(tmp, name)

        np.savetxt(path("s.txt"), v, fmt="%.17g")
        np.savetxt(path("s.dat"), v, fmt="%.17g", header="vx vy vz", comments="")
        np.savetxt(path("s.csv"), v, fmt="%.17g", delimiter=",", header="vx,vy,vz",
                   comments="")
        np.savetxt(path("commented.txt"), v, fmt="%.17g", header="written by a test")
        np.save(path("s.npy"), v)
        v.tofile(path("s.bin"))
        v.tofile(path("s.raw"))
        v.ravel()[:-1].tofile(path("short.bin"))
        np.savetxt(path("two.txt"), v[:, :2])
        np.savetxt(path("small.txt"), v[:MIN_SAMPLES - 1])
        with open(path("broken.txt"), "w") as fh:
            fh.write("1 2 3\n4 5 6\n7 x 9\n")
        open(path("empty.txt"), "w").close()

        for name, extra in (("s.txt", []), ("s.dat", []), ("s.csv", []),
                            ("commented.txt", []), ("s.npy", []), ("s.bin", []),
                            ("s.raw", ["--binary"])):
            code, out, err = run([path(name)] + params + extra)
            check(" ".join(["correct sample in", name] + extra + ["exits 0"]),
                  code == 0 and "overall: PASS" in out, f"exit {code}, stderr {err!r}")
        code, _, err = run([path("s.dat")] + params)
        check("a header line in a text .dat file is skipped with a note",
              "skipped the header line 'vx vy vz'" in err)

        code, out, _ = run([path("s.txt"), "--kappa", "1", "--theta-perp", "1",
                            "--theta-par", "2"])
        check("wrong kappa exits 1", code == 1 and "overall: FAIL" in out, f"exit {code}")

        code, out, _ = run([path("s.txt"), "--json"] + params)
        check("--json prints the report as JSON",
              code == 0 and '"holm_adjusted_pvalue"' in out)

        for argv, why in (
            ([path("two.txt")] + params, "two-column file"),
            ([path("small.txt")] + params, f"fewer than {MIN_SAMPLES} draws"),
            ([path("broken.txt")] + params, "non-numeric line after the first"),
            ([path("empty.txt")] + params, "empty file"),
            ([path("short.bin")] + params, "raw file of the wrong length"),
            ([path("missing.txt")] + params, "missing file"),
            ([path("s.txt")] + params + ["--alpha", "5"], "alpha = 5"),
            ([path("s.txt"), "--kappa", "0.5", "--theta-perp", "1", "--theta-par", "2"],
             "kappa = 1/2"),
            ([path("s.txt"), "--kappa", "nan", "--theta-perp", "1", "--theta-par", "2"],
             "kappa = nan"),
            ([path("s.txt"), "--kappa", "2", "--theta-perp", "-1", "--theta-par", "2"],
             "negative theta_perp"),
            ([path("s.txt"), "--kappa", "2", "--theta-perp", "1", "--theta-par", "inf"],
             "infinite theta_par"),
            ([path("s.txt")] + params + ["--bhat", "0", "0", "0"], "zero bhat"),
            ([path("s.txt")] + params + ["--attempts", str(n - 1)], "attempts < N"),
        ):
            code, out, err = run(argv)
            lines = err.strip().splitlines()
            check(f"{why} exits 2 with a one-line message",
                  code == 2 and out == "" and len(lines) == 1
                  and lines[0].startswith("bikappa_validate.py: error:"),
                  f"exit {code}, stderr {err.strip()!r}")

        code, _, _ = run([path("s.txt"), "--theta-perp", "1", "--theta-par", "2"])
        check("a missing required option exits 2", code == 2, f"exit {code}")


SECTIONS = (positive_control, rotated_frame, negative_controls, bound_accounting,
            input_contract, report_rendering, command_line)


def run_section(section) -> list:
    check = Checks()
    section(np.random.default_rng([SEED, SECTIONS.index(section)]), check)
    return check.failures


# pytest entry points, one per section
def test_positive_control():
    assert not run_section(positive_control)


def test_rotated_frame():
    assert not run_section(rotated_frame)


def test_negative_controls():
    assert not run_section(negative_controls)


def test_bound_accounting():
    assert not run_section(bound_accounting)


def test_input_contract():
    assert not run_section(input_contract)


def test_report_rendering():
    assert not run_section(report_rendering)


def test_command_line():
    assert not run_section(command_line)


def main() -> int:
    failures = []
    for section in SECTIONS:
        failures += run_section(section)
    print()
    if failures:
        print(f"FAILED: {len(failures)} check(s): {failures}")
        return 1
    print("all checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
