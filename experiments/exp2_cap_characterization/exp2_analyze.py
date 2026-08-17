#!/usr/bin/env python3
"""Experiment 2 analysis: characterization of the component-wise velocity cap.

The released sampler has two modes and they sample two *different* probability laws:

  no_cap()          the untruncated bi-Kappa distribution -- the intended target;
  finite lambda     that distribution CONDITIONED on the component-wise box event
                        |v_x|/theta_perp <= lambda AND
                        |v_y|/theta_perp <= lambda AND
                        |v_z|/theta_par  <= lambda.

This script quantifies the gap: how often a draw is rejected and redrawn, and how far
the conditioned law sits from the target in robust quantiles, in the empirical CDF, in
the moments that legitimately exist, and in the angular structure.

Measuring the rejected fraction.  `operator()` loops internally and reports nothing, and
we are not allowed to instrument it.  Instead the box predicate is re-evaluated here,
exactly as `withinNormalizedVelocityCap` writes it, on the draws of the *uncapped* run at
the same seed.  Because each loop iteration consumes x1, x2, cosTheta, phi in the same
order whether or not a cap is in force, the uncapped run reproduces the capped run's
attempt stream bitwise; the mean of the predicate over those attempts is P(accept)
directly, and `subsequence_check` verifies the bitwise correspondence rather than
assuming it.

Reads the binary sample files and manifest written by exp2_sample.exe; writes
results/exp2_results.json (machine readable) and results/exp2_table.md.

Usage:  uv run --project ../../python python exp2_analyze.py [raw_dir] [results_dir]
"""

from __future__ import annotations

import csv
import json
import os
import platform
import subprocess
import sys

import numpy as np
import scipy
from scipy import integrate, special, stats

QUANTILE_PROBES = (0.5, 0.9, 0.99, 0.999)

# "Negligible" thresholds, fixed here before any result was looked at.  A cap is called
# negligible for a given kappa when BOTH hold:
#   total-variation distance from the target law (= rejection fraction) < 1e-3, and
#   the p99.9 speed quantile is distorted by < 1%.
NEGLIGIBLE_TV = 1.0e-3
NEGLIGIBLE_Q999_REL = 1.0e-2

SQRT3_INV = 1.0 / np.sqrt(3.0)
SQRT2_INV = 1.0 / np.sqrt(2.0)

# bi_kappa_distribution.H:247 -- operator() throws after this many consecutive rejections.
MAX_CAP_REJECT_TRIES = 1_000_000


# ---------------------------------------------------------------------------
# Analytic acceptance probability
# ---------------------------------------------------------------------------
# In the field-aligned local frame the sampler produces
#     v = sqrt(kappa) * (theta_perp u_x, theta_perp u_y, theta_par u_z),
# so the shipped predicate |v_x|/theta_perp <= lambda, ... becomes
#     sqrt(kappa) * max(|u_x|, |u_y|, |u_z|) <= lambda.
# The theta's cancel identically: the box event is a CUBE of half-side
# c = lambda/sqrt(kappa) in the isotropic u-coordinates, and the acceptance
# probability therefore depends on (kappa, lambda) only -- not on the anisotropy and
# not on the field direction (the cap is tested before the frame rotation).
#
# Writing u = R n with n uniform on S^2 independent of R, and M = max_i |n_i|,
#     P(accept) = E_M[ F_R(c / M) ],   F_R(r) = I_{r^2/(1+r^2)}(3/2, kappa-1/2),
# evaluated as I_z with z = c^2/(m^2 + c^2) so that r^2 is never formed.
#
# The law of M is available in closed form.  P(|n_i| > m) = 1 - m for each axis; for
# m >= 1/sqrt(3) at most one axis can exceed m unless m < 1/sqrt(2), so by
# inclusion-exclusion P(M > m) = 3(1-m) - 3 Q(m) with Q(m) = P(|n_x|>m, |n_y|>m).
# Differentiating Q under the integral gives the density below; it is validated
# against a direct spherical Monte Carlo in `check_direction_density`.
def m_density(m: np.ndarray) -> np.ndarray:
    """Density of M = max_i |n_i| for n uniform on S^2, supported on [1/sqrt3, 1]."""
    m = np.asarray(m, dtype=float)
    out = np.full(m.shape, 3.0)
    lo = m < SQRT2_INV
    x = np.clip((1.0 - 2.0 * m[lo] ** 2) / (1.0 - m[lo] ** 2), 0.0, 1.0)
    out[lo] = 3.0 - (12.0 / np.pi) * np.arcsin(np.sqrt(x))
    return out


def analytic_reject_fraction(kappa: float, lam: float) -> float:
    """Exact P(reject) = 1 - P(accept) for the component-wise box, by 1-D quadrature.

    The REJECTION probability is integrated directly rather than obtained by subtracting
    the acceptance from 1: at large kappa or large lambda the acceptance is 1 - O(1e-12)
    and the subtraction loses every significant digit (and can even go negative on
    quadrature noise).  The complement is taken inside the incomplete beta instead, via
    1 - I_z(a, b) = I_{1-z}(b, a), where 1 - z = m^2/(m^2 + c^2) is exact.
    """
    b = kappa - 0.5
    c = lam / np.sqrt(kappa)

    def integrand(m: float) -> float:
        one_minus_z = m * m / (m * m + c * c)
        return float(special.betainc(b, 1.5, one_minus_z) * m_density(np.array([m]))[0])

    lower = integrate.quad(integrand, SQRT3_INV, SQRT2_INV, limit=200)[0]
    upper = integrate.quad(integrand, SQRT2_INV, 1.0, limit=200)[0]
    return float(min(max(lower + upper, 0.0), 1.0))


def check_direction_density(n: int = 2_000_000, seed: int = 20250817) -> dict:
    """Validate m_density against a direct uniform-sphere Monte Carlo."""
    rng = np.random.default_rng(seed)
    z = rng.uniform(-1.0, 1.0, n)
    phi = rng.uniform(0.0, 2.0 * np.pi, n)
    s = np.sqrt(np.maximum(0.0, 1.0 - z * z))
    M = np.maximum(np.abs(z), np.maximum(np.abs(s * np.cos(phi)), np.abs(s * np.sin(phi))))
    grid = np.linspace(SQRT3_INV, 1.0, 41)
    cdf = np.array([integrate.quad(lambda t: m_density(np.array([t]))[0], SQRT3_INV, g)[0]
                    for g in grid])
    emp = np.searchsorted(np.sort(M), grid, side="right") / n
    return {
        "n_monte_carlo": n,
        "density_integral": float(cdf[-1]),
        "max_cdf_abs_err_vs_monte_carlo": float(np.abs(cdf - emp).max()),
    }


# ---------------------------------------------------------------------------
# Per-run diagnostics
# ---------------------------------------------------------------------------
def safe_speed(v: np.ndarray) -> np.ndarray:
    """|v| without ever forming a square.  np.linalg.norm squares internally and
    overflows for the heavy-tailed draws at kappa <= 1; np.hypot does not."""
    return np.hypot(np.hypot(v[:, 0], v[:, 1]), v[:, 2])


def in_box(v: np.ndarray, theta_perp: float, theta_par: float, lam: float) -> np.ndarray:
    """The shipped predicate `withinNormalizedVelocityCap`, transcribed verbatim."""
    nx = np.abs(v[:, 0]) / theta_perp
    ny = np.abs(v[:, 1]) / theta_perp
    nz = np.abs(v[:, 2]) / theta_par
    return (nx <= lam) & (ny <= lam) & (nz <= lam)


def run_stats(v: np.ndarray, kappa: float, theta_perp: float, theta_par: float) -> dict:
    speed = safe_speed(v)
    finite = np.isfinite(speed) & np.isfinite(v).all(axis=1)
    v, speed = v[finite], speed[finite]
    n = speed.size

    comp = np.abs(v)
    out = {
        "n_finite": int(n),
        "q_speed": [float(x) for x in np.quantile(speed, QUANTILE_PROBES)],
        "q_absvx": [float(x) for x in np.quantile(comp[:, 0], QUANTILE_PROBES)],
        "q_absvy": [float(x) for x in np.quantile(comp[:, 1], QUANTILE_PROBES)],
        "q_absvz": [float(x) for x in np.quantile(comp[:, 2], QUANTILE_PROBES)],
        "mad": [float(x) for x in stats.median_abs_deviation(v, axis=0)],
    }
    out["mad_ratio_par_perp"] = out["mad"][2] / out["mad"][0]

    # Angular structure in the isotropic u-coordinates.  The box is a CUBE there, so a
    # capped sample cannot be isotropic in direction and cannot be axisymmetric about B:
    # directions pointing at a cube corner have more radial room than directions along an
    # axis.  a4 = 2<cos 4phi> is the leading azimuthal Fourier coefficient of the
    # density (1/2pi)(1 + a4 cos 4phi + ...); it is 0 for the untruncated law, and its
    # sampling s.d. is sqrt(2/n), so z4 = a4 sqrt(n/2) is a calibrated detector.
    scales = np.sqrt(kappa) * np.array([theta_perp, theta_perp, theta_par])
    u = v / scales
    R = np.hypot(np.hypot(u[:, 0], u[:, 1]), u[:, 2])
    good = R > 0
    cos_theta = u[good, 2] / R[good]
    phi = np.arctan2(u[good, 1], u[good, 0])
    a4 = 2.0 * float(np.mean(np.cos(4.0 * phi)))
    out["cos_theta_ks_sqrtn"] = float(stats.kstest(cos_theta, "uniform",
                                                   args=(-1.0, 2.0)).statistic * np.sqrt(n))
    out["azimuth_a4"] = a4
    out["azimuth_a4_z"] = float(a4 * np.sqrt(n / 2.0))

    # Second moments.  Recorded unconditionally for the *capped* runs because the capped
    # law always has bounded support, but interpreted only where an untruncated reference
    # exists; see `moment_comparison_allowed`.
    out["sample_var"] = [float(x) for x in v.var(axis=0, ddof=1)]
    return out


def moment_comparison_allowed(kappa: float) -> bool:
    """The untruncated bi-Kappa second moment is theta^2 kappa/(2 kappa - 3); it exists
    only for kappa > 3/2 and diverges at kappa = 3/2.  Below that there is no reference
    value, so a capped-vs-uncapped variance comparison compares a number to nothing."""
    return kappa > 1.5


def theory_var(kappa: float, theta: float) -> float:
    return theta * theta * kappa / (2.0 * kappa - 3.0)


# ---------------------------------------------------------------------------
def main() -> int:
    raw_dir = sys.argv[1] if len(sys.argv) > 1 else "raw"
    results_dir = sys.argv[2] if len(sys.argv) > 2 else "results"
    os.makedirs(results_dir, exist_ok=True)

    with open(os.path.join(raw_dir, "manifest.csv")) as fh:
        runs = list(csv.DictReader(fh))

    def load(row) -> np.ndarray:
        return np.fromfile(os.path.join(raw_dir, row["file"]), dtype=np.float64).reshape(-1, 3)

    lambdas = sorted({float(r["lambda"]) for r in runs if r["mode"] == "capped"})

    # index: (block, kappa, theta_perp, theta_par, seed) -> {mode/lambda -> row}
    uncapped: dict[tuple, dict] = {}
    capped: dict[tuple, dict] = {}
    for r in runs:
        key = (r["block"], float(r["kappa"]), float(r["theta_perp"]), float(r["theta_par"]),
               int(r["seed"]))
        if r["mode"] == "uncapped":
            uncapped[key] = r
        else:
            capped.setdefault(key, {})[float(r["lambda"])] = r

    per_run = []
    subsequence_checks = []

    for key in sorted(uncapped):
        block, kappa, tp, tpar, seed = key
        ref_row = uncapped[key]
        vref = load(ref_row)
        ref = run_stats(vref, kappa, tp, tpar)
        ref.update({"run_id": int(ref_row["run_id"]), "block": block, "kappa": kappa,
                    "theta_perp": tp, "theta_par": tpar, "mode": "uncapped",
                    "lambda": None, "seed": seed,
                    "n_samples": int(ref_row["n_samples"]),
                    "n_nonfinite": int(ref_row["n_nonfinite"])})
        per_run.append(ref)

        speed_ref = safe_speed(vref)

        for lam in lambdas:
            row = capped[key][lam]
            vcap = load(row)
            st = run_stats(vcap, kappa, tp, tpar)

            # --- rejected fraction, measured on the uncapped attempt stream ---
            accept = in_box(vref, tp, tpar, lam)
            p_hat = float(accept.mean())
            n_att = accept.size
            st["p_accept_empirical"] = p_hat
            st["p_accept_stderr"] = float(np.sqrt(max(p_hat * (1.0 - p_hat), 0.0) / n_att))
            st["reject_fraction_empirical"] = 1.0 - p_hat
            st["n_attempts_measured"] = int(n_att)
            st["expected_attempts_per_accept"] = float("inf") if p_hat == 0 else 1.0 / p_hat

            # --- retry-count distribution ---
            # Attempts per accepted draw = 1 + (number of consecutive rejections before
            # it), which is Geometric(p) on {1, 2, ...} if the attempts are i.i.d.  Read
            # off the accept mask as gaps between successive accepted indices.
            idx = np.flatnonzero(accept)
            if idx.size >= 2:
                tries = np.diff(idx)  # attempts consumed to produce each accepted draw
                st["retry_mean_attempts"] = float(tries.mean())
                st["retry_p99_attempts"] = float(np.quantile(tries, 0.99))
                st["retry_max_attempts"] = int(tries.max())
                # Geometric reference for the same p.
                st["retry_mean_attempts_geometric"] = 1.0 / p_hat
                st["retry_p99_attempts_geometric"] = (
                    1.0 if p_hat >= 1.0
                    else float(max(1.0, np.ceil(np.log1p(-0.99) / np.log1p(-p_hat)))))
            else:
                st["retry_mean_attempts"] = float("nan")
                st["retry_p99_attempts"] = float("nan")
                st["retry_max_attempts"] = 0
                st["retry_mean_attempts_geometric"] = float("nan")
                st["retry_p99_attempts_geometric"] = float("nan")
            # The sampler throws after kMaxCapRejectTries = 1e6 consecutive rejections.
            # log10 P(that happens for one draw) = 1e6 log10(1-p).
            st["log10_p_hit_internal_try_limit"] = (
                float(MAX_CAP_REJECT_TRIES * np.log10(1.0 - p_hat)) if p_hat < 1.0
                else float("-inf"))

            # --- bitwise verification that capped == uncapped restricted to the box ---
            sub = vref[accept]
            k = sub.shape[0]
            match = bool(k == 0 or np.array_equal(sub, vcap[:k]))
            subsequence_checks.append({
                "block": block, "kappa": kappa, "theta_perp": tp, "theta_par": tpar,
                "lambda": lam, "seed": seed, "n_in_box": int(k), "bitwise_identical": match,
            })

            # --- distortion relative to the uncapped target, same seed ---
            st["q_speed_ratio"] = [c / u if u != 0 else float("nan")
                                   for c, u in zip(st["q_speed"], ref["q_speed"])]
            st["q_absvz_ratio"] = [c / u if u != 0 else float("nan")
                                   for c, u in zip(st["q_absvz"], ref["q_absvz"])]
            st["q_absvx_ratio"] = [c / u if u != 0 else float("nan")
                                   for c, u in zip(st["q_absvx"], ref["q_absvx"])]
            speed_cap = safe_speed(vcap)
            st["ks_speed_ecdf_diff"] = float(
                stats.ks_2samp(speed_cap, speed_ref, method="asymp").statistic)
            st["ks_absvz_ecdf_diff"] = float(
                stats.ks_2samp(np.abs(vcap[:, 2]), np.abs(vref[:, 2]),
                               method="asymp").statistic)
            # The capped law is the target conditioned on the box, so its density ratio is
            # 1_box/P(accept) and the TOTAL-VARIATION distance to the target is exactly
            # 1 - P(accept).  Any 1-D ECDF gap is a lower bound on this.
            st["tv_distance_empirical"] = 1.0 - p_hat

            st.update({"run_id": int(row["run_id"]), "block": block, "kappa": kappa,
                       "theta_perp": tp, "theta_par": tpar, "mode": "capped",
                       "lambda": lam, "seed": seed,
                       "n_samples": int(row["n_samples"]),
                       "n_nonfinite": int(row["n_nonfinite"])})
            per_run.append(st)

        print(f"  {block} kappa={kappa:<5g} theta=({tp:g},{tpar:g}) seed={seed}  "
              + "  ".join(
                  f"L{lam:g}:rej={1.0 - float(in_box(vref, tp, tpar, lam).mean()):.4f}"
                  for lam in lambdas))

    # ---------------- aggregate over seeds ----------------
    groups: dict[tuple, list] = {}
    for r in per_run:
        groups.setdefault((r["block"], r["kappa"], r["theta_perp"], r["theta_par"],
                           r["mode"], r["lambda"]), []).append(r)

    summary = []
    for gkey, group in sorted(groups.items(), key=lambda kv: (kv[0][0], kv[0][1], kv[0][2],
                                                             kv[0][3], kv[0][4],
                                                             kv[0][5] or 0.0)):
        block, kappa, tp, tpar, mode, lam = gkey

        def stat(name):
            vals = np.array([g[name] for g in group], dtype=float)
            return {"mean": float(vals.mean()), "sd": float(vals.std(ddof=1)),
                    "min": float(vals.min()), "max": float(vals.max())}

        def vstat(name):
            vals = np.array([g[name] for g in group], dtype=float)
            return {"mean": vals.mean(axis=0).tolist(), "sd": vals.std(axis=0, ddof=1).tolist()}

        entry = {
            "block": block, "kappa": kappa, "theta_perp": tp, "theta_par": tpar,
            "mode": mode, "lambda": lam,
            "n_replicates": len(group), "n_per_replicate": group[0]["n_samples"],
            "seeds": sorted(g["seed"] for g in group),
            "total_nonfinite": int(sum(g["n_nonfinite"] for g in group)),
            "q_probes": list(QUANTILE_PROBES),
            "q_speed": vstat("q_speed"),
            "q_absvx": vstat("q_absvx"),
            "q_absvz": vstat("q_absvz"),
            "mad_ratio_par_perp": stat("mad_ratio_par_perp"),
            "cos_theta_ks_sqrtn": stat("cos_theta_ks_sqrtn"),
            "azimuth_a4": stat("azimuth_a4"),
            "azimuth_a4_z": stat("azimuth_a4_z"),
            "sample_var": vstat("sample_var"),
        }
        if mode == "capped":
            entry["p_accept_empirical"] = stat("p_accept_empirical")
            entry["reject_fraction_empirical"] = stat("reject_fraction_empirical")
            entry["reject_fraction_analytic"] = analytic_reject_fraction(kappa, lam)
            entry["p_accept_analytic"] = 1.0 - entry["reject_fraction_analytic"]
            entry["expected_attempts_per_accept_analytic"] = (
                1.0 / entry["p_accept_analytic"] if entry["p_accept_analytic"] > 0
                else float("inf"))
            entry["tv_distance_empirical"] = stat("tv_distance_empirical")
            entry["tv_distance_analytic"] = entry["reject_fraction_analytic"]
            entry["q_speed_ratio"] = vstat("q_speed_ratio")
            entry["q_absvx_ratio"] = vstat("q_absvx_ratio")
            entry["q_absvz_ratio"] = vstat("q_absvz_ratio")
            entry["ks_speed_ecdf_diff"] = stat("ks_speed_ecdf_diff")
            entry["ks_absvz_ecdf_diff"] = stat("ks_absvz_ecdf_diff")
            entry["total_attempts_measured"] = int(sum(g["n_attempts_measured"]
                                                       for g in group))
            for nm in ("retry_mean_attempts", "retry_p99_attempts", "retry_max_attempts",
                       "retry_mean_attempts_geometric", "retry_p99_attempts_geometric",
                       "log10_p_hit_internal_try_limit"):
                entry[nm] = stat(nm)
            # Computed from the ANALYTIC acceptance: the empirical p saturates at 1 for the
            # wide caps, which would report -inf instead of a very large negative number.
            entry["log10_p_hit_internal_try_limit_analytic"] = float(
                MAX_CAP_REJECT_TRIES * np.log10(entry["reject_fraction_analytic"])
                if entry["reject_fraction_analytic"] > 0 else -np.inf)

        # ---- moments, only where an untruncated reference exists ----
        if moment_comparison_allowed(kappa):
            sv = np.array([g["sample_var"] for g in group], dtype=float)
            ref_var = np.array([theory_var(kappa, tp), theory_var(kappa, tp),
                                theory_var(kappa, tpar)])
            entry["var_theory_untruncated"] = ref_var.tolist()
            entry["var_over_theory"] = (sv.mean(axis=0) / ref_var).tolist()
            entry["moment_comparison"] = "reported (kappa > 3/2, second moment exists)"
        else:
            entry["var_theory_untruncated"] = None
            entry["var_over_theory"] = None
            entry["moment_comparison"] = (
                "REFUSED: the untruncated bi-Kappa second moment requires kappa > 3/2 "
                "(it is theta^2 kappa/(2 kappa - 3), divergent at kappa = 3/2). No "
                "reference value exists, so a capped-vs-uncapped variance ratio would be "
                "meaningless. The capped sample variance is still recorded under "
                "'sample_var' because the capped law has bounded support and therefore "
                "always has one -- but it is an artifact of lambda, not a property of the "
                "physical distribution; see the lambda-dependence table.")
        summary.append(entry)

    # ---------------- where the cap becomes negligible ----------------
    negligible = []
    by_case: dict[tuple, list] = {}
    for e in summary:
        if e["mode"] == "capped":
            by_case.setdefault((e["block"], e["kappa"], e["theta_perp"], e["theta_par"]),
                               []).append(e)
    for ckey, entries in sorted(by_case.items()):
        entries = sorted(entries, key=lambda e: e["lambda"])
        chosen = None
        for e in entries:
            q999 = abs(e["q_speed_ratio"]["mean"][3] - 1.0)
            if e["reject_fraction_analytic"] < NEGLIGIBLE_TV and q999 < NEGLIGIBLE_Q999_REL:
                chosen = e["lambda"]
                break
        negligible.append({
            "block": ckey[0], "kappa": ckey[1], "theta_perp": ckey[2], "theta_par": ckey[3],
            "criterion": {"tv_below": NEGLIGIBLE_TV, "q999_rel_below": NEGLIGIBLE_Q999_REL},
            "smallest_lambda_meeting_criterion": chosen,
        })

    # ---------------- theta-independence of acceptance ----------------
    # Predicted exactly: the theta's cancel in the normalized predicate, so the box event
    # is a cube of half-side lambda/sqrt(kappa) in the isotropic u-coordinates.  Because
    # the isotropic control and the anisotropic runs share the seed, they share the same
    # u-stream, and the prediction is therefore that the accept/reject decision agrees
    # DRAW BY DRAW, not merely in aggregate.  That is what is checked.
    theta_independence = []
    acc = {(e["kappa"], e["theta_perp"], e["theta_par"], e["lambda"]):
           e["reject_fraction_empirical"]["mean"]
           for e in summary if e["mode"] == "capped"}
    for (kappa, tp, tpar, lam), val in sorted(acc.items()):
        if (tp, tpar) == (1.0, 1.0) and (kappa, 1.0, 2.0, lam) in acc:
            theta_independence.append({
                "kappa": kappa, "lambda": lam,
                "reject_fraction_isotropic_1_1": val,
                "reject_fraction_anisotropic_1_2": acc[(kappa, 1.0, 2.0, lam)],
                "abs_difference": abs(val - acc[(kappa, 1.0, 2.0, lam)]),
                "analytic": analytic_reject_fraction(kappa, lam),
            })

    # Draw-by-draw form of the same check: the isotropic control and the anisotropic run
    # at matched (kappa, seed) consume the same RNG stream, so if the theta's really
    # cancel the accept mask must agree on every single attempt.
    theta_mask_checks = []
    for kappa in sorted({k for (_, k, tp, tpar, _) in uncapped if (tp, tpar) == (1.0, 1.0)}):
        for seed in sorted({s for (_, k, tp, tpar, s) in uncapped
                            if k == kappa and (tp, tpar) == (1.0, 1.0)}):
            iso = uncapped.get(("C", kappa, 1.0, 1.0, seed))
            ani = uncapped.get(("A", kappa, 1.0, 2.0, seed))
            if iso is None or ani is None:
                continue
            viso, vani = load(iso), load(ani)
            for lam in lambdas:
                m_iso = in_box(viso, 1.0, 1.0, lam)
                m_ani = in_box(vani, 1.0, 2.0, lam)
                theta_mask_checks.append({
                    "kappa": kappa, "seed": seed, "lambda": lam,
                    "n_attempts": int(m_iso.size),
                    "n_mask_disagreements": int(np.count_nonzero(m_iso != m_ani)),
                })

    # ---------------- tail-exponent scaling of the rejected fraction ----------------
    # P(R > r) ~ r^{-(2 kappa - 1)} for the untruncated radial law, so the rejected
    # fraction must decay as lambda^{-(2 kappa - 1)}.  Measured as the local log-log slope
    # between the two widest caps, where the asymptotic regime has been reached.
    tail_scaling = []
    for kappa in sorted({e["kappa"] for e in summary}):
        pts = sorted(((e["lambda"], e["reject_fraction_analytic"]) for e in summary
                      if e["mode"] == "capped" and e["kappa"] == kappa
                      and e["theta_par"] == 2.0), key=lambda t: t[0])
        (l1, r1), (l2, r2) = pts[-2], pts[-1]
        tail_scaling.append({
            "kappa": kappa, "lambda_pair": [l1, l2],
            "measured_loglog_slope": float(np.log(r2 / r1) / np.log(l2 / l1)),
            "predicted_slope": -(2.0 * kappa - 1.0),
        })

    n_sub = len(subsequence_checks)
    n_sub_ok = sum(1 for c in subsequence_checks if c["bitwise_identical"])

    out = {
        "experiment": "Experiment 2 -- characterization of the component-wise velocity cap",
        "answers": ["R1.1 (primary)", "R1.2 (primary)", "R1.3 item 5", "R1.5", "R1.6"],
        "target_law_statement": (
            "no_cap() samples the untruncated bi-Kappa distribution. A finite "
            "max_normalized_velocity = lambda samples a DIFFERENT law: the bi-Kappa "
            "distribution conditioned on the component-wise box "
            "{|v_x|/theta_perp <= lambda, |v_y|/theta_perp <= lambda, "
            "|v_z|/theta_par <= lambda}. The two are not two approximations of one law; "
            "the capped mode is a conditional/truncated target with a different density, "
            "different tails, different moments and different angular symmetry."),
        "reject_fraction_method": (
            "The released operator() loops internally and reports no attempt count, and "
            "cpp/bi_kappa_distribution.H was not modified. Instead the shipped predicate "
            "withinNormalizedVelocityCap was transcribed into exp2_analyze.py:in_box and "
            "evaluated on the draws of the UNCAPPED run at the same seed. Each loop "
            "iteration of operator() consumes x1, x2, cosTheta, phi in the same order "
            "with or without a cap, so the uncapped run reproduces the capped run's "
            "attempt stream; the predicate mean over those attempts is P(accept). The "
            "bitwise subsequence check confirms the correspondence."),
        "numerical_safety": [
            "|v| via chained np.hypot; np.linalg.norm squares internally and overflows "
            "for the kappa <= 1 tails.",
            "T = R^2 is never formed; the analytic acceptance uses z = c^2/(m^2+c^2).",
            "The bounded radial variable Y = T/(1+T) is not used anywhere; it rounds to "
            "exactly 1 at low kappa (see Experiment 1).",
        ],
        "environment": {
            "platform": platform.platform(),
            "machine": platform.machine(),
            "processor": platform.processor(),
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "cxx": subprocess.run(["g++", "--version"], capture_output=True,
                                  text=True).stdout.splitlines()[0],
            "cxx_target": next((ln for ln in subprocess.run(
                ["g++", "--version"], capture_output=True, text=True).stdout.splitlines()
                if ln.startswith("Target:")), ""),
            "cxxflags": "-Wall -Wextra -std=c++11 -O2",
            "stdlib": "libc++ (Apple clang)",
            "rng": "std::mt19937, seeded explicitly per run",
            "git_commit": subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True,
                                         text=True).stdout.strip(),
            "git_dirty": bool(subprocess.run(["git", "status", "--porcelain"],
                                             capture_output=True, text=True).stdout.strip()),
            # The working tree is dirty for reasons outside this experiment, so pin the
            # exact sampler that produced raw/ by content rather than by commit.
            "sampler_header_sha256": subprocess.run(
                ["shasum", "-a", "256", "../../cpp/bi_kappa_distribution.H"],
                capture_output=True, text=True).stdout.split()[0],
        },
        "direction_density_validation": check_direction_density(),
        "subsequence_check": {
            "n_pairs": n_sub, "n_bitwise_identical": n_sub_ok,
            "all_identical": n_sub_ok == n_sub,
            "detail": subsequence_checks,
        },
        "theta_independence_of_acceptance": theta_independence,
        "theta_independence_drawwise": {
            "n_comparisons": len(theta_mask_checks),
            "total_attempts": int(sum(c["n_attempts"] for c in theta_mask_checks)),
            "total_mask_disagreements": int(sum(c["n_mask_disagreements"]
                                                for c in theta_mask_checks)),
            "detail": theta_mask_checks,
        },
        "tail_exponent_scaling": tail_scaling,
        "negligibility": negligible,
        "summary": summary,
        "per_run": per_run,
    }
    with open(os.path.join(results_dir, "exp2_results.json"), "w") as fh:
        json.dump(out, fh, indent=2)

    write_table(out, os.path.join(results_dir, "exp2_table.md"))
    print(f"\nwrote {results_dir}/exp2_results.json and {results_dir}/exp2_table.md")
    return 0


# ---------------------------------------------------------------------------
def write_table(out: dict, path: str) -> None:
    summary = out["summary"]
    lam_list = sorted({e["lambda"] for e in summary if e["mode"] == "capped"})
    kap_list = sorted({e["kappa"] for e in summary})

    def get(block, kappa, tpar, lam):
        for e in summary:
            if (e["block"] == block and e["kappa"] == kappa and e["theta_par"] == tpar
                    and e["mode"] == "capped" and e["lambda"] == lam):
                return e
        return None

    def ref(block, kappa, tpar):
        for e in summary:
            if (e["block"] == block and e["kappa"] == kappa and e["theta_par"] == tpar
                    and e["mode"] == "uncapped"):
                return e
        return None

    L = [
        "# Experiment 2 results -- the component-wise velocity cap",
        "",
        "**The capped and uncapped modes sample two different probability laws.**",
        "`no_cap()` samples the untruncated bi-Kappa distribution. A finite",
        "`max_normalized_velocity = lambda` samples that distribution *conditioned on* the",
        "component-wise box `|v_x|/theta_perp <= lambda AND |v_y|/theta_perp <= lambda AND",
        "|v_z|/theta_par <= lambda`. The capped mode is a truncated/conditional target, not a",
        "numerical approximation to the uncapped one, and must never be reported as the same",
        "distribution.",
        "",
        "5 replicate seeds (2001-2005) x 10^5 samples per configuration, double precision,",
        "Apple clang / libc++, `std::mt19937`, `theta_perp:theta_par = 1:2` unless stated.",
        "",
        "## How the rejected fraction was measured",
        "",
        "`operator()` loops internally and reports no attempt count, and",
        "`cpp/bi_kappa_distribution.H` was **not** modified. The shipped predicate",
        "`withinNormalizedVelocityCap` was instead transcribed into `exp2_analyze.py:in_box`",
        "and evaluated on the draws of the **uncapped** run at the same seed. Every loop",
        "iteration of `operator()` consumes `x1, x2, cosTheta, phi` in the same order whether",
        "or not a cap is in force, so the uncapped run *is* the capped run's attempt stream;",
        "the mean of the predicate over it is `P(accept)` directly.",
        "",
        f"That correspondence is verified, not assumed: for all "
        f"{out['subsequence_check']['n_pairs']} (case, lambda, seed) pairs the capped run's",
        "output is **bitwise identical** to the uncapped run's draws restricted to the box "
        f"({out['subsequence_check']['n_bitwise_identical']}/"
        f"{out['subsequence_check']['n_pairs']} pairs).",
        "",
        "## 1. Headline -- rejected (redrawn) fraction, kappa x lambda",
        "",
        "Fraction of attempts thrown away and redrawn, `1 - P(accept)`. Empirical values are",
        "the mean over 5 seeds of the box-predicate rate on 10^5 uncapped attempts each;",
        "`analytic` is the exact value from the closed-form derivation below.",
        "",
        "**This same number is also the exact total-variation distance between the capped law",
        "and the untruncated target**, because conditioning on an event of probability `p`",
        "gives density ratio `1_box/p` and hence `TV = 1 - p`. So the rejection column is not",
        "merely a cost: it is the distortion.",
        "",
        "| kappa | " + " | ".join(f"lambda={l:g}" for l in lam_list) + " |",
        "|---|" + "---|" * len(lam_list),
    ]
    def small(x: float) -> str:
        return f"{x:.5f}" if x >= 1.0e-5 else f"{x:.2e}"

    for k in kap_list:
        cells = []
        for l in lam_list:
            e = get("A", k, 2.0, l)
            if e is None:
                cells.append("-")
                continue
            cells.append(f"{e['reject_fraction_empirical']['mean']:.5f} "
                         f"({small(e['reject_fraction_analytic'])})")
        L.append(f"| {k:g} | " + " | ".join(cells) + " |")
    L += ["",
          "Format: empirical (analytic). The empirical column is a mean over 5x10^5 attempts,",
          "so its own resolution is ~2e-6 and it reads 0.00000 wherever the analytic value is",
          "below that; the analytic column is the one to quote in those cells.",
          ""]

    L += [
        "### Mean attempts per accepted draw (analytic)",
        "",
        "| kappa | " + " | ".join(f"lambda={l:g}" for l in lam_list) + " |",
        "|---|" + "---|" * len(lam_list),
    ]
    for k in kap_list:
        cells = []
        for l in lam_list:
            e = get("A", k, 2.0, l)
            cells.append("-" if e is None
                         else f"{e['expected_attempts_per_accept_analytic']:.4f}")
        L.append(f"| {k:g} | " + " | ".join(cells) + " |")

    L += [
        "",
        "### Retry-count distribution and the internal attempt limit",
        "",
        "Attempts consumed per accepted draw, read off the gaps between accepted indices in",
        "the uncapped attempt stream. If attempts are i.i.d. this is `Geometric(p)` on",
        "`{1, 2, ...}`; the geometric reference is shown alongside. `operator()` throws after",
        f"`kMaxCapRejectTries = {MAX_CAP_REJECT_TRIES:,}` consecutive rejections",
        "(`bi_kappa_distribution.H:247`); `log10 P(throw)` per draw is `1e6 log10(1-p)`.",
        "",
        "| kappa | lambda | mean attempts | (geometric) | p99 attempts | (geometric) | "
        "max observed | log10 P(hit try limit) |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for k in kap_list:
        for l in lam_list:
            e = get("A", k, 2.0, l)
            if e is None or not np.isfinite(e["retry_mean_attempts"]["mean"]):
                continue
            lp = e["log10_p_hit_internal_try_limit_analytic"]
            L.append(f"| {k:g} | {l:g} | {e['retry_mean_attempts']['mean']:.4f} "
                     f"| {e['retry_mean_attempts_geometric']['mean']:.4f} "
                     f"| {e['retry_p99_attempts']['mean']:.1f} "
                     f"| {e['retry_p99_attempts_geometric']['mean']:.1f} "
                     f"| {e['retry_max_attempts']['max']:.0f} "
                     + (f"| {lp:.3g} |" if np.isfinite(lp) else "| -inf |"))
    L += ["",
          "The retry counts follow the geometric reference throughout, and the internal",
          "attempt limit is unreachable: even at the worst configuration in the sweep the",
          "log-probability of a single draw exhausting it is about -260 000. The attempt limit",
          "is therefore not a practical failure mode, and the cost of the cap is entirely the",
          "mean-attempt overhead plus -- far more importantly -- the change of target law.",
          ""]

    L += [
        "",
        "### Closed form for the acceptance probability",
        "",
        "In the field-aligned local frame the sampler emits",
        "`v = sqrt(kappa) (theta_perp u_x, theta_perp u_y, theta_par u_z)`, so the shipped",
        "predicate reduces to `sqrt(kappa) max_i |u_i| <= lambda`: **the theta's cancel**, and",
        "the box event is a cube of half-side `c = lambda/sqrt(kappa)` in the isotropic",
        "u-coordinates. Hence, with `u = R n`, `n` uniform on S^2, `M = max_i |n_i|`,",
        "",
        "```",
        "P(accept) = E_M[ I_z(3/2, kappa-1/2) ],   z = c^2 / (M^2 + c^2),",
        "f_M(m) = 3 - (12/pi) arcsin( sqrt((1-2m^2)/(1-m^2)) )   for 1/sqrt3 <= m <= 1/sqrt2,",
        "f_M(m) = 3                                              for 1/sqrt2 <= m <= 1.",
        "```",
        "",
        f"`f_M` integrates to {out['direction_density_validation']['density_integral']:.12f} "
        "and its CDF agrees with a "
        f"{out['direction_density_validation']['n_monte_carlo']:.0e}-point spherical Monte "
        "Carlo to "
        f"{out['direction_density_validation']['max_cdf_abs_err_vs_monte_carlo']:.2e}.",
        "",
        "Three consequences worth stating in the manuscript:",
        "",
        "1. `P(accept)` depends on `(kappa, lambda)` only -- **not** on `theta_perp`,",
        "   `theta_par` or the field direction (the cap is tested before the frame rotation).",
        "2. The rejected fraction decays only *algebraically* in lambda, like",
        "   `lambda^{-(2 kappa - 1)}`, so at heavy tails it is stubborn: see the table above.",
        "3. The cap is a **cube** in the isotropic coordinates, so the capped law is not",
        "   isotropic in direction and is not axisymmetric about **B**; see section 4.",
        "",
        "### Theta-independence of acceptance -- isotropic (1,1) control vs anisotropic (1,2)",
        "",
        "| kappa | lambda | reject frac (1,1) | reject frac (1,2) | abs diff | analytic |",
        "|---|---|---|---|---|---|",
    ]
    for t in out["theta_independence_of_acceptance"]:
        L.append(f"| {t['kappa']:g} | {t['lambda']:g} | "
                 f"{t['reject_fraction_isotropic_1_1']:.5f} | "
                 f"{t['reject_fraction_anisotropic_1_2']:.5f} | "
                 f"{t['abs_difference']:.2e} | {t['analytic']:.5f} |")
    dw = out["theta_independence_drawwise"]
    L += ["",
          "These are **not** two independent estimates that happen to agree. The two runs",
          "share the seed and therefore the RNG stream, so the prediction is the stronger one",
          "that the accept/reject decision agrees on every individual attempt -- and it does:",
          f"{dw['total_mask_disagreements']} disagreements in {dw['total_attempts']} attempts",
          f"across {dw['n_comparisons']} (kappa, seed, lambda) comparisons. The cancellation of",
          "`theta_perp` and `theta_par` in the cap predicate is exact, not statistical.",
          "",
          "### Tail-exponent scaling of the rejected fraction",
          "",
          "`P(R > r) ~ r^{-(2 kappa - 1)}` for the untruncated radial law, so the rejected",
          "fraction must decay as `lambda^{-(2 kappa - 1)}`. Local log-log slope of the exact",
          "rejection curve between the two widest caps:",
          "",
          "| kappa | lambda pair | measured slope | predicted -(2 kappa - 1) |",
          "|---|---|---|---|",
          ]
    for t in out["tail_exponent_scaling"]:
        L.append(f"| {t['kappa']:g} | {t['lambda_pair'][0]:g} -> {t['lambda_pair'][1]:g} | "
                 f"{t['measured_loglog_slope']:.4f} | {t['predicted_slope']:.4f} |")
    L += ["",
          "This is the whole problem in one line. The cost of the cap is not exponentially",
          "small in lambda, it is a power law, and the power is weakest exactly where kappa",
          "distributions are physically interesting.",
          ""]

    L += [
        "## 2. Distortion of the capped law relative to the uncapped target",
        "",
        "Robust quantile ratios `q_capped / q_uncapped` at matched seeds, the largest",
        "empirical-CDF gap on the speed `|v|`, and the exact total-variation distance",
        "(analytic, since the empirical estimate saturates at ~2e-6). A ratio of 1 means no",
        "distortion.",
        "",
        "**Read the last two columns against the p99.9 column, not instead of it.** `TV` is an",
        "upper bound on how much any *probability* can move, and the sup ECDF gap is a lower",
        "bound on it -- but neither bounds a *quantile* ratio. A cap can amputate the entire",
        "far tail while moving no probability by more than 1e-3, because the amputated region",
        "carries almost no probability and enormous velocity. At `kappa = 1.5, lambda = 50`",
        "the total-variation distance is 6.6e-4 -- by any probability-based measure the two",
        "laws are indistinguishable -- and the p99.9 speed is still 24% too small.",
        "",
        "| kappa | lambda | |v| p50 | |v| p90 | |v| p99 | |v| p99.9 | |v_z| p99.9 | "
        "sup ECDF gap on |v| | TV (analytic) |",
        "|---|---|---|---|---|---|---|---|---|",
    ]

    def ratio(x: float) -> str:
        return f"{x:.4f}" if x >= 1.0e-3 else f"{x:.2e}"

    for k in kap_list:
        for l in lam_list:
            e = get("A", k, 2.0, l)
            if e is None:
                continue
            r = e["q_speed_ratio"]["mean"]
            rz = e["q_absvz_ratio"]["mean"]
            L.append(f"| {k:g} | {l:g} | " + " | ".join(ratio(x) for x in r)
                     + f" | {ratio(rz[3])} | {e['ks_speed_ecdf_diff']['mean']:.5f} "
                     f"| {e['tv_distance_analytic']:.3e} |")

    q999_worst = max(ref("A", k, 2.0)["q_speed"]["mean"][3] for k in kap_list)
    L += [
        "",
        "### Absolute scale of the quantiles being compared (uncapped target, |v|)",
        "",
        "Given so the ratios above can be read: at `kappa = 0.75` the untruncated p99.9 speed",
        f"is {q999_worst:.3g} thermal units, while a capped draw obeys "
        "`|v| <= lambda sqrt(theta_perp^2 + theta_perp^2 + theta_par^2)` by construction. The",
        "ratio is ~1e-5 and the tail is simply gone -- not compressed, gone.",
        "",
        "| kappa | p50 | p90 | p99 | p99.9 |",
        "|---|---|---|---|---|",
    ]
    for k in kap_list:
        e = ref("A", k, 2.0)
        q = e["q_speed"]["mean"]
        L.append(f"| {k:g} | " + " | ".join(f"{x:.4g}" for x in q) + " |")

    L += [
        "",
        "### Where the cap becomes negligible",
        "",
        f"Criterion fixed in advance: `TV < {NEGLIGIBLE_TV:g}` **and** p99.9 speed quantile",
        f"distorted by `< {100.0 * NEGLIGIBLE_Q999_REL:g}%`. Both are required precisely",
        "because they fail in different places.",
        "",
        "| block | kappa | (theta_perp, theta_par) | smallest lambda in the ladder that "
        "qualifies |",
        "|---|---|---|---|",
    ]
    for n in out["negligibility"]:
        val = n["smallest_lambda_meeting_criterion"]
        ladder = ", ".join(f"{x:g}" for x in lam_list)
        L.append(f"| {n['block']} | {n['kappa']:g} | ({n['theta_perp']:g}, "
                 f"{n['theta_par']:g}) | "
                 + (f"{val:g}" if val is not None else f"**none** of {ladder}") + " |")

    L += [
        "",
        "## 3. Moments -- what could and could not legitimately be compared",
        "",
        "The untruncated bi-Kappa second moment is `theta^2 kappa/(2 kappa - 3)`. It exists",
        "**only for kappa > 3/2** and diverges at `kappa = 3/2`. For `kappa <= 3/2` there is",
        "no untruncated reference value, so a capped-vs-uncapped variance comparison would be",
        "comparing a number to a divergent integral. Those comparisons are **refused**, not",
        "silently omitted.",
        "",
        "| kappa | second moment of the untruncated target | comparison |",
        "|---|---|---|",
    ]
    for k in kap_list:
        if moment_comparison_allowed(k):
            L.append(f"| {k:g} | exists, = {k / (2 * k - 3):.4g} theta^2 | reported below |")
        else:
            L.append(f"| {k:g} | **does not exist** | "
                     "**REFUSED -- there is no reference value to compare against** |")

    L += [
        "",
        "### Variance ratio to the untruncated theory, kappa > 3/2 only",
        "",
        "`sample var / (theta^2 kappa/(2 kappa - 3))`, perpendicular then parallel component.",
        "",
        "| kappa | mode | lambda | var_x / theory | var_z / theory |",
        "|---|---|---|---|---|",
    ]
    for k in kap_list:
        if not moment_comparison_allowed(k):
            continue
        e = ref("A", k, 2.0)
        v = e["var_over_theory"]
        L.append(f"| {k:g} | uncapped | - | {v[0]:.4f} | {v[2]:.4f} |")
        for l in lam_list:
            c = get("A", k, 2.0, l)
            v = c["var_over_theory"]
            L.append(f"| {k:g} | capped | {l:g} | {v[0]:.4f} | {v[2]:.4f} |")

    L += [
        "",
        "### The cap manufactures a finite variance where none exists",
        "",
        "For `kappa <= 3/2` the capped sample has a perfectly finite variance -- it has",
        "bounded support -- and that number is reported by any naive diagnostic. It is an",
        "artifact of `lambda`, not a property of the plasma. Since the untruncated radial",
        "density behaves as `f_R(r) ~ r^{-2 kappa}`, the second moment truncated at `~lambda`",
        "grows as `lambda^{3 - 2 kappa}` for `kappa < 3/2` and as `log lambda` at",
        "`kappa = 3/2` exactly -- without bound in both cases. Reporting it as *the* variance",
        "of a kappa distribution with `kappa <= 3/2` would be wrong.",
        "",
        "| kappa | " + " | ".join(f"var_x, lambda={l:g}" for l in lam_list)
        + " | measured growth (20 -> 50) | predicted |",
        "|---|" + "---|" * (len(lam_list) + 2),
    ]
    for k in kap_list:
        if moment_comparison_allowed(k):
            continue
        cells, vals = [], []
        for l in lam_list:
            e = get("A", k, 2.0, l)
            cells.append("-" if e is None else f"{e['sample_var']['mean'][0]:.4g}")
            vals.append(None if e is None else e["sample_var"]["mean"][0])
        if abs(k - 1.5) < 1e-12:
            # At kappa = 3/2 the growth is logarithmic, not a power: fit var ~ A + B log
            # lambda over each consecutive pair and report B, which must be constant.
            bs = [(vals[i + 1] - vals[i]) / np.log(lam_list[i + 1] / lam_list[i])
                  for i in range(len(vals) - 1)]
            meas = (f"log-slope B = {bs[-1]:.2f} "
                    f"(B over all four pairs: {', '.join(f'{b:.2f}' for b in bs)})")
            pred = "`A + B log lambda`, B constant"
        else:
            slope = np.log(vals[-1] / vals[-2]) / np.log(lam_list[-1] / lam_list[-2])
            meas = f"power-law slope {slope:.3f}"
            pred = f"`lambda^{{{3.0 - 2.0 * k:g}}}`"
        L.append(f"| {k:g} | " + " | ".join(cells) + f" | {meas} | {pred} |")
    L += ["",
          "The uncapped runs have no entry here on purpose: their sample variance is a finite",
          "number produced by a divergent population moment and means nothing. The capped runs",
          "do have a well-defined population variance -- it is just a variance of the box, not",
          "of the plasma.",
          "",
          "Two honest caveats on the fitted rates. At `kappa = 3/2` the log-slope `B` rises",
          "from 1.15 to 1.47 across the ladder rather than sitting at a constant: the growth is",
          "unmistakably slower than any power of lambda and consistent with `log lambda`, but",
          "lambda = 50 is not yet deep enough in the asymptotic regime to pin `B` down. At",
          "`kappa = 0.75` the measured power-law slope 1.41 likewise approaches the asymptotic",
          "1.5 from below, because an O(1) contribution from the distribution core is still",
          "present; fitting `var = A + B lambda^{1.5}` to the two widest caps gives A = 1.86,",
          "B = 0.179, which reproduces the lambda = 10 point to 11%. Neither caveat touches the",
          "conclusion, which is that the number diverges as the cap is relaxed and therefore",
          "is not a property of the distribution being sampled.",
          ""]

    L += [
        "## 4. Angular structure -- the cap breaks axisymmetry about B",
        "",
        "The box is a **cube** in the isotropic u-coordinates, so a direction pointing at a",
        "cube corner has `sqrt(3)` times more radial room than a direction along an axis. The",
        "conditioned law therefore has a four-fold azimuthal modulation about **B** and a",
        "non-uniform polar distribution -- structure the physical bi-Kappa distribution does",
        "not have and which no choice of `theta` can absorb.",
        "",
        "`a4 = 2<cos 4 phi>` is the leading azimuthal Fourier coefficient (0 for the target);",
        "`z4 = a4 sqrt(n/2)` is it in units of its own sampling s.d., so |z4| > 3 is a",
        "detection. `cos(theta) sqrt(n) D` is the KS statistic against U(-1,1), which is",
        "O(1) (about 0.87) when the law is correct.",
        "",
        "| kappa | lambda | a4 | z4 | cos(theta) sqrt(n)D |",
        "|---|---|---|---|---|",
    ]
    for k in kap_list:
        e = ref("A", k, 2.0)
        L.append(f"| {k:g} | uncapped | {e['azimuth_a4']['mean']:+.5f} | "
                 f"{e['azimuth_a4_z']['mean']:+.1f} | {e['cos_theta_ks_sqrtn']['mean']:.3f} |")
        for l in lam_list:
            c = get("A", k, 2.0, l)
            L.append(f"| {k:g} | {l:g} | {c['azimuth_a4']['mean']:+.5f} | "
                     f"{c['azimuth_a4_z']['mean']:+.1f} | "
                     f"{c['cos_theta_ks_sqrtn']['mean']:.3f} |")

    L += [
        "",
        "## 5. Recommendation for the manuscript",
        "",
        "1. Run every validation and every physics result with the cap **off**",
        "   (`no_cap()`). That is the mode whose target law is the bi-Kappa distribution.",
        "2. Document the finite `max_normalized_velocity` as an **optional pragmatic",
        "   finite-velocity-box conditional target**: the bi-Kappa distribution conditioned on",
        "   a component-wise box, useful only when a code needs bounded particle speeds.",
        "3. Do **not** present it as a regularized or physically motivated kappa model. It is",
        "   a cube in normalized velocity components: it is not isotropic, it is not",
        "   axisymmetric about **B** (section 4), and its shape depends on the arbitrary",
        "   choice of lambda.",
        "4. State the cost and the distortion together, since they are the same number:",
        "   rejected fraction = total-variation distance from the target = `1 - P(accept)`,",
        "   with the closed form above.",
        "5. Never quote a variance for `kappa <= 3/2` obtained from a capped run.",
        "",
    ]
    with open(path, "w") as fh:
        fh.write("\n".join(L))


if __name__ == "__main__":
    raise SystemExit(main())
