#!/usr/bin/env python3
"""Experiment 1 analysis: radial law, directional uniformity, independence, frame invariance.

Answers reviewer comment R1.3 by testing the construction that the sampler is actually
built on -- T = R^2 ~ BetaPrime(3/2, kappa-1/2) -- rather than Cartesian marginals alone,
and by covering the low-kappa range 1/2 < kappa <= 3/2 where the second moment does not
exist and variance-based checks are unavailable.

Reads the binary sample files and manifest written by exp1_sample.exe; writes
results/exp1_results.json (machine readable, full per-run detail) and
results/exp1_table.md (the summary table for the manuscript).

Usage:  python/.venv/bin/python exp1_analyze.py [raw_dir] [results_dir]
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
from scipy import special, stats

# Q-Q / quantile probes.  The two extreme ones probe the heavy tail, which is where a
# radial-law error would show up first and where Cartesian marginals are least sensitive.
QUANTILE_PROBES = (0.5, 0.9, 0.99, 0.999)

# Independence test binning: radial quartiles against 40 direction cells.  The solid-angle
# element is d(cos theta) d(phi), so equal intervals of cos(theta) crossed with equal
# intervals of phi give cells of equal solid angle -- no sphere-pixelization scheme needed.
# Binning both angles is what makes this a test of the whole direction: with cos(theta)
# alone, a dependence between radius and azimuth would pass unseen.
N_RADIAL_BINS = 4
N_COS_BINS = 5
N_PHI_BINS = 8

# The manuscript's test: Pearson's chi^2 over cells of equal probability under the target.
# N_SHELLS radius shells of equal probability, bounded by the exact deciles of R, crossed
# with the N_COS_BINS x N_PHI_BINS equal-solid-angle direction cells above.  Because the
# radius and the direction are independent under the target and the direction is uniform,
# every one of the N_SHELLS * 40 cells has probability 1/400.  The shells are compared in
# log R against log edges, so no bounded transform of R is formed and nothing underflows.
N_SHELLS = 10


def log_shell_edges(kappa: float) -> np.ndarray:
    """log of the radii that split R into N_SHELLS shells of equal probability.

    P(T <= t) = P(W >= 1/(1+t)) with W = 1/(1+T) ~ Beta(kappa-1/2, 3/2); the edge at
    cumulative probability q therefore has W_q = Beta.ppf(1-q) and R_q^2 = (1-W_q)/W_q.
    W_q stays well inside the double range for every kappa tested (R_q <= 1.4e50 at 0.51).
    """
    q = np.arange(1, N_SHELLS) / N_SHELLS
    w = stats.beta.ppf(1.0 - q, kappa - 0.5, 1.5)
    return 0.5 * (np.log1p(-w) - np.log(w))


def cell_pvalues(counts: np.ndarray) -> tuple[float, float]:
    """Pearson chi^2 p-values of (radius shells alone, all cells) against equal occupancy."""
    return (float(stats.chisquare(counts.sum(axis=1)).pvalue),
            float(stats.chisquare(counts.ravel()).pvalue))


def field_basis(ub: np.ndarray) -> np.ndarray:
    """Rebuild the orthonormal field-aligned basis exactly as bi_kappa_distribution does.

    Returns M with columns (e1, e2, e3) so that v_global = M @ v_local.  Reproducing the
    C++ construction here rather than assuming a canonical basis is the point: it is what
    makes the recovered local coordinates a genuine test of the shipped rotation.
    """
    ub = np.asarray(ub, dtype=float)
    e3 = ub / np.linalg.norm(ub)

    e2 = np.ones(3)
    maxcomp = 0
    if abs(e3[1]) > abs(e3[maxcomp]):
        maxcomp = 1
    if abs(e3[2]) > abs(e3[maxcomp]):
        maxcomp = 2
    e2[maxcomp] = 1.0 - e3.sum() / e3[maxcomp]
    e2 /= np.linalg.norm(e2)

    e1 = np.cross(e2, e3)
    return np.column_stack([e1, e2, e3])


def is_axis_aligned_z(ub: np.ndarray) -> bool:
    """Mirror the C++ short-circuit: for ub == +z the sampler skips the rotation entirely."""
    eps = np.finfo(float).eps

    def approx(a: float, b: float) -> bool:
        return abs(a - b) <= eps * max(1.0, abs(a), abs(b))

    return approx(ub[0], 0.0) and approx(ub[1], 0.0) and approx(ub[2], 1.0)


def analyze_run(path: str, kappa: float, theta_perp: float, theta_par: float,
                ub: np.ndarray) -> dict:
    v = np.fromfile(path, dtype=np.float64).reshape(-1, 3)
    n_total = v.shape[0]

    # --- recover the field-aligned frame ---
    if is_axis_aligned_z(ub):
        local = v
        basis_orthonormality_err = 0.0
    else:
        M = field_basis(ub)
        local = v @ M  # local = M^T v
        basis_orthonormality_err = float(np.abs(M.T @ M - np.eye(3)).max())

    # --- normalized (dimensionless) coordinates: the isotropic core ---
    scales = np.sqrt(kappa) * np.array([theta_perp, theta_perp, theta_par])
    u = local / scales

    # Work in R, never in T = R^2.  Forming T first overflows for R > ~1.3e154 even
    # though R itself is perfectly representable -- the same trap that the sampler used
    # to fall into at cpp/bi_kappa_distribution.H:273, and it would silently discard
    # valid heavy-tail draws from the diagnostic exactly where the diagnostic matters.
    # np.linalg.norm squares internally and so overflows here; hypot does not.
    R = np.hypot(np.hypot(u[:, 0], u[:, 1]), u[:, 2])
    finite = np.isfinite(R) & np.isfinite(local).all(axis=1)
    n_nonfinite = int((~finite).sum())

    u, R, local_f = u[finite], R[finite], local[finite]
    n = R.size

    b = kappa - 0.5  # BetaPrime(3/2, kappa-1/2); Beta(3/2, kappa-1/2) after Y = T/(1+T)

    # --- radial law ---
    # PRIMARY diagnostic: W = 1/(1+T) ~ Beta(kappa-1/2, 3/2), the *complement* of the
    # usual bounded transform.  Both W and Y = T/(1+T) = 1-W are exact bijections of T,
    # so they carry identical information in exact arithmetic -- but not in doubles.
    # For small b the mass of Y piles up against 1, where doubles have no relative
    # resolution left (see n_y_saturated below), and a KS test then measures rounding
    # rather than the sampler.  W puts that same mass near 0, where relative resolution
    # is full down to ~1e-308.  Orientation of the bounded transform is not cosmetic.
    with np.errstate(over="ignore", divide="ignore", under="ignore"):
        small = R <= 1.0
        W = np.empty_like(R)
        W[small] = 1.0 / (1.0 + R[small] ** 2)
        inv2 = (1.0 / R[~small]) ** 2
        W[~small] = inv2 / (1.0 + inv2)

        Y = 1.0 - W  # kept only to quantify the failure of the naive orientation
        Y_direct = np.empty_like(R)
        Y_direct[small] = R[small] ** 2 / (1.0 + R[small] ** 2)
        Y_direct[~small] = 1.0 / (1.0 + inv2)

    n_y_saturated = int((Y_direct >= 1.0).sum())  # draws beyond double resolution of Y
    n_w_underflow = int((W <= 0.0).sum())         # draws beyond double resolution of W

    # The radial tests run on the probability integral transform F_W(W), which is
    # U(0,1) under the target.  Where W underflows to 0 (R above ~1e161, reached at
    # kappa near 1/2) its CDF is still resolvable: for W -> 0,
    # F_W(W) = W^b / (b B(b, 3/2)) (1 + O(W)), and log W = -2 log R to within a
    # relative 1/R^2.  At b = 0.01 these draws carry about 6e-4 of the probability,
    # so mapping them all to F = 0 would displace them rather than drop them, and
    # would add up to ~0.2 to sqrt(n) D on its own.
    with np.errstate(divide="ignore"):
        U = stats.beta.cdf(W, b, 1.5)
        under = W <= 0.0
        U[under] = np.exp(-2.0 * b * np.log(R[under])
                          - np.log(b) - special.betaln(b, 1.5))
    ks_radial = stats.kstest(U, "uniform")
    cvm_radial = stats.cramervonmises(U, "uniform")
    ks_radial_y = stats.kstest(Y_direct, "beta", args=(1.5, b))

    # Quantiles are compared in log R.  With b < 1 the upper quantiles of T span
    # hundreds of decades, so a relative error on R is dominated by tail sampling noise
    # and says nothing; a difference in log R is the meaningful scale.
    with np.errstate(divide="ignore"):
        tq = 0.5 * np.log(stats.betaprime.ppf(QUANTILE_PROBES, 1.5, b))
        eq = np.log(np.quantile(R, QUANTILE_PROBES))
    quantile_log_err = np.where(np.isfinite(tq), eq - tq, np.nan)

    # --- directional uniformity ---
    cos_theta = u[:, 2] / R
    phi = np.arctan2(u[:, 1], u[:, 0])
    ks_cos = stats.kstest(cos_theta, "uniform", args=(-1.0, 2.0))
    ks_phi = stats.kstest(phi, "uniform", args=(-np.pi, 2.0 * np.pi))

    # --- radial-direction independence ---
    r_edges = np.quantile(R, np.linspace(0.0, 1.0, N_RADIAL_BINS + 1))
    r_bin = np.clip(np.searchsorted(r_edges[1:-1], R, side="right"), 0, N_RADIAL_BINS - 1)
    c_bin = np.clip(np.floor((cos_theta + 1.0) / 2.0 * N_COS_BINS).astype(int),
                    0, N_COS_BINS - 1)
    p_bin = np.clip(np.floor((phi + np.pi) / (2.0 * np.pi) * N_PHI_BINS).astype(int),
                    0, N_PHI_BINS - 1)
    cell = c_bin * N_PHI_BINS + p_bin
    table = np.zeros((N_RADIAL_BINS, N_COS_BINS * N_PHI_BINS))
    np.add.at(table, (r_bin, cell), 1.0)
    chi2 = stats.chi2_contingency(table)
    # --- the manuscript's cell-count test ------------------------------------
    shell = np.searchsorted(log_shell_edges(kappa), np.log(R), side="right")
    cell_counts = np.zeros((N_SHELLS, N_COS_BINS * N_PHI_BINS), dtype=np.int64)
    np.add.at(cell_counts, (shell, cell), 1)
    cells_radius_p, cells_all_p = cell_pvalues(cell_counts)

    # The same counts summed over radius test the joint direction law: the two KS tests
    # above check cos(theta) and phi separately, and two uniform marginals do not make
    # the direction uniform on the sphere unless the angles are also independent.
    direction_gof = stats.chisquare(table.sum(axis=0))
    # Direct check too: is the direction law the same in the innermost and outermost
    # radial quartiles?  A contingency test can wash out a localized tail effect.
    inner = cos_theta[R <= r_edges[1]]
    outer = cos_theta[R > r_edges[-2]]
    ks_inner_outer = stats.ks_2samp(inner, outer)

    result = {
        "n_samples": n_total,
        "n_nonfinite": n_nonfinite,
        "nonfinite_fraction": n_nonfinite / n_total,
        "n_y_saturated": n_y_saturated,
        "n_w_underflow": n_w_underflow,
        "basis_orthonormality_err": basis_orthonormality_err,
        # sqrt(n)*D is reported because it is O(1) when the law is correct, independent
        # of n -- unlike the p-value, which collapses to 0 for any n large enough.
        "radial_ks_stat": float(ks_radial.statistic),
        "radial_ks_sqrtn": float(ks_radial.statistic * np.sqrt(n)),
        "radial_ks_pvalue": float(ks_radial.pvalue),
        "radial_cvm_stat": float(cvm_radial.statistic),
        "radial_cvm_pvalue": float(cvm_radial.pvalue),
        # The same test run on the naive orientation Y = T/(1+T), for comparison only.
        "radial_ks_sqrtn_naive_Y": float(ks_radial_y.statistic * np.sqrt(n)),
        "quantile_probes": list(QUANTILE_PROBES),
        "quantile_log_err": [float(x) for x in quantile_log_err],
        "cos_theta_ks_sqrtn": float(ks_cos.statistic * np.sqrt(n)),
        "cos_theta_ks_pvalue": float(ks_cos.pvalue),
        "phi_ks_sqrtn": float(ks_phi.statistic * np.sqrt(n)),
        "phi_ks_pvalue": float(ks_phi.pvalue),
        "independence_chi2_pvalue": float(chi2.pvalue),
        "independence_chi2_dof": int(chi2.dof),
        "independence_min_cell_count": int(table.min()),
        "direction_cells_chi2_pvalue": float(direction_gof.pvalue),
        "independence_inner_outer_ks_pvalue": float(ks_inner_outer.pvalue),
        "cells_radius_pvalue": cells_radius_p,
        "cells_all_pvalue": cells_all_p,
        "cell_counts": cell_counts.tolist(),
    }

    # --- anisotropy alignment, valid at every kappa ---
    # MAD needs no moments, so this works in the 1/2 < kappa <= 3/2 range where the
    # variance does not exist.  A misaligned frame rotation shows up here immediately.
    mad = stats.median_abs_deviation(local_f, axis=0)
    result["mad_ratio_par_perp"] = float(mad[2] / mad[0])
    result["mad_ratio_expected"] = float(theta_par / theta_perp)

    # --- moments, only where they exist ---
    if kappa > 1.5:
        var_theory_perp = theta_perp**2 * kappa / (2.0 * kappa - 3.0)
        var_theory_par = theta_par**2 * kappa / (2.0 * kappa - 3.0)
        sample_var = local_f.var(axis=0, ddof=1)
        result["var_rel_err"] = [
            float((sample_var[0] - var_theory_perp) / var_theory_perp),
            float((sample_var[1] - var_theory_perp) / var_theory_perp),
            float((sample_var[2] - var_theory_par) / var_theory_par),
        ]
        cov = np.cov(local_f, rowvar=False)
        d = np.sqrt(np.diag(cov))
        corr = cov / np.outer(d, d)
        result["max_offdiag_corr"] = float(np.abs(corr - np.diag(np.diag(corr))).max())
    else:
        result["var_rel_err"] = None
        result["max_offdiag_corr"] = None

    return result


def load_normalized(raw_dir: str, row: dict) -> np.ndarray:
    """Recover the normalized velocity u = A^-1 Q^T v of every draw in a run."""
    v = np.fromfile(os.path.join(raw_dir, row["file"]), dtype=np.float64).reshape(-1, 3)
    ub = np.array([float(row["ub_x"]), float(row["ub_y"]), float(row["ub_z"])])
    local = v if is_axis_aligned_z(ub) else v @ field_basis(ub)
    kappa = float(row["kappa"])
    return local / (np.sqrt(kappa) * np.array(
        [float(row["theta_perp"]), float(row["theta_perp"]), float(row["theta_par"])]))


def load_radii(raw_dir: str, row: dict, label: str) -> np.ndarray:
    """Recover the normalized radius R from a run, undoing the field-frame rotation."""
    v = np.fromfile(os.path.join(raw_dir, row["file"]), dtype=np.float64).reshape(-1, 3)
    ub = np.array([float(row["ub_x"]), float(row["ub_y"]), float(row["ub_z"])])
    local = v if is_axis_aligned_z(ub) else v @ field_basis(ub)
    kappa = float(row["kappa"])
    scales = np.sqrt(kappa) * np.array(
        [float(row["theta_perp"]), float(row["theta_perp"]), float(row["theta_par"])]
    )
    u = local / scales
    return np.hypot(np.hypot(u[:, 0], u[:, 1]), u[:, 2])


def main() -> int:
    raw_dir = sys.argv[1] if len(sys.argv) > 1 else "raw"
    results_dir = sys.argv[2] if len(sys.argv) > 2 else "results"
    os.makedirs(results_dir, exist_ok=True)

    with open(os.path.join(raw_dir, "manifest.csv")) as fh:
        runs = list(csv.DictReader(fh))

    per_run = []
    for row in runs:
        ub = np.array([float(row["ub_x"]), float(row["ub_y"]), float(row["ub_z"])])
        res = analyze_run(
            os.path.join(raw_dir, row["file"]),
            float(row["kappa"]), float(row["theta_perp"]), float(row["theta_par"]), ub,
        )
        res.update({
            "run_id": int(row["run_id"]), "block": row["block"],
            "kappa": float(row["kappa"]), "theta_perp": float(row["theta_perp"]),
            "theta_par": float(row["theta_par"]), "ub_label": row["ub_label"],
            "ub": ub.tolist(), "mode": row["mode"], "seed": int(row["seed"]),
        })
        # The sampler counts a draw bad when a velocity component is non-finite; the
        # analysis counts it bad when the recovered radius is.  They must agree, and
        # they only do because both work in R rather than in R^2.
        res["n_nonfinite_sampler"] = int(row["n_nonfinite"])
        assert res["n_nonfinite"] == res["n_nonfinite_sampler"], (
            f"run {row['run_id']}: sampler saw {row['n_nonfinite']} non-finite draws, "
            f"analysis saw {res['n_nonfinite']}"
        )
        per_run.append(res)
        print(f"  run {res['run_id']:3d}  block {res['block']}  kappa={res['kappa']:<6g} "
              f"{res['ub_label']:<8s} seed={res['seed']}  "
              f"radial sqrt(n)D={res['radial_ks_sqrtn']:.3f}")

    # --- aggregate across replicate seeds ---
    configs: dict[tuple, list] = {}
    for r in per_run:
        key = (r["block"], r["kappa"], r["theta_perp"], r["theta_par"], r["ub_label"])
        configs.setdefault(key, []).append(r)

    summary = []
    for key, group in configs.items():
        block, kappa, tp, tpar, ub_label = key

        def stat(name):
            vals = np.array([g[name] for g in group], dtype=float)
            return {"mean": float(vals.mean()), "sd": float(vals.std(ddof=1)),
                    "min": float(vals.min()), "max": float(vals.max())}

        entry = {
            "block": block, "kappa": kappa, "theta_perp": tp, "theta_par": tpar,
            "ub_label": ub_label, "n_replicates": len(group),
            "n_per_replicate": group[0]["n_samples"],
            "seeds": sorted(g["seed"] for g in group),
            "total_nonfinite": int(sum(g["n_nonfinite"] for g in group)),
            "total_y_saturated": int(sum(g["n_y_saturated"] for g in group)),
            "total_w_underflow": int(sum(g["n_w_underflow"] for g in group)),
            "radial_ks_sqrtn": stat("radial_ks_sqrtn"),
            "radial_ks_sqrtn_naive_Y": stat("radial_ks_sqrtn_naive_Y"),
            "radial_cvm_stat": stat("radial_cvm_stat"),
            "cos_theta_ks_sqrtn": stat("cos_theta_ks_sqrtn"),
            "phi_ks_sqrtn": stat("phi_ks_sqrtn"),
            "independence_chi2_pvalue": stat("independence_chi2_pvalue"),
            "independence_min_cell_count": stat("independence_min_cell_count"),
            "direction_cells_chi2_pvalue": stat("direction_cells_chi2_pvalue"),
            "independence_inner_outer_ks_pvalue": stat("independence_inner_outer_ks_pvalue"),
            "mad_ratio_par_perp": stat("mad_ratio_par_perp"),
            "basis_orthonormality_err": stat("basis_orthonormality_err"),
            "max_tail_quantile_log_err": float(
                np.nanmax([abs(np.array(g["quantile_log_err"])).max() for g in group])
            ),
        }
        # Cell test on the replicates pooled: 5 x 10^5 draws, 1250 expected per cell.
        pooled = np.sum([np.array(g["cell_counts"]) for g in group], axis=0)
        entry["cells_pooled_radius_pvalue"], entry["cells_pooled_all_pvalue"] = cell_pvalues(pooled)
        entry["cells_pooled_expected"] = float(pooled.sum() / pooled.size)
        entry["cells_pooled_min_count"] = int(pooled.min())
        entry["cells_pooled_shell_counts"] = pooled.sum(axis=1).tolist()
        if kappa > 1.5:
            varr = np.array([g["var_rel_err"] for g in group], dtype=float)
            entry["var_rel_err_mean"] = varr.mean(axis=0).tolist()
            entry["var_rel_err_sd"] = varr.std(axis=0, ddof=1).tolist()
            entry["max_offdiag_corr"] = stat("max_offdiag_corr")
        summary.append(entry)

    summary.sort(key=lambda e: (e["block"], e["kappa"], e["theta_par"], e["ub_label"]))

    # --- frame invariance, stated as a direct comparison rather than inferred ---
    # Same (kappa, theta, seed) driven through different field directions must produce
    # the *same* normalized radial sample, because the field transform is an orthogonal
    # rotation applied after the radial draw and the RNG stream does not depend on ub.
    # Comparing the recovered radii bitwise is a far sharper test than comparing summary
    # statistics, and it is what makes the C2 claim a validated one rather than an API one.
    frame_checks = []
    by_phys: dict[tuple, dict] = {}
    for r in per_run:
        if r["block"] != "B":
            continue
        by_phys.setdefault(
            (r["kappa"], r["theta_perp"], r["theta_par"], r["seed"]), {}
        )[r["ub_label"]] = r["run_id"]
    for key, runs_by_dir in sorted(by_phys.items()):
        if "z" not in runs_by_dir:
            continue
        ref = load_radii(raw_dir, runs[runs_by_dir["z"]], "z")
        for label, rid in sorted(runs_by_dir.items()):
            if label == "z":
                continue
            got = load_radii(raw_dir, runs[rid], label)
            m = np.isfinite(ref) & np.isfinite(got)
            rel = np.abs(got[m] - ref[m]) / np.maximum(np.abs(ref[m]), 1e-300)
            frame_checks.append({
                "kappa": key[0], "theta_perp": key[1], "theta_par": key[2], "seed": key[3],
                "direction": label,
                "max_rel_radius_diff": float(rel.max()),
                "n_bitwise_identical": int((got[m] == ref[m]).sum()),
                "n_compared": int(m.sum()),
            })

    # --- the same comparison on the whole normalized velocity -------------------
    # Every block-B run shares its random numbers with the (theta_par = 2, B || z) run at
    # the same kappa and seed, so after undoing its own scaling and rotation it must
    # return the same normalized velocity u, draw by draw.  Comparing the vector rather
    # than only |u| also covers the direction, and covering theta_par = 1 checks that the
    # thermal speeds enter exactly as specified.
    ref_rows = {(float(r["kappa"]), int(r["seed"])): r for r in runs
                if r["block"] == "B" and float(r["theta_par"]) == 2.0 and r["ub_label"] == "z"}
    vector_checks = []
    for r in runs:
        if r["block"] != "B":
            continue
        ref_row = ref_rows[(float(r["kappa"]), int(r["seed"]))]
        if ref_row is r:
            continue
        ref, got = load_normalized(raw_dir, ref_row), load_normalized(raw_dir, r)
        m = np.isfinite(ref).all(axis=1) & np.isfinite(got).all(axis=1)
        ref_norm = np.hypot(np.hypot(ref[m, 0], ref[m, 1]), ref[m, 2])
        d = got[m] - ref[m]
        rel = np.hypot(np.hypot(d[:, 0], d[:, 1]), d[:, 2]) / np.maximum(ref_norm, 1e-300)
        vector_checks.append({
            "kappa": float(r["kappa"]), "seed": int(r["seed"]),
            "theta_par": float(r["theta_par"]), "direction": r["ub_label"],
            "max_rel_vector_diff": float(rel.max()), "n_compared": int(m.sum()),
        })

    out = {
        "experiment": "Experiment 1 -- radial, directional, anisotropic and frame validation",
        "answers": ["R1.3 (primary)", "R1.5", "R2.A2", "R2.A3", "R2.C1", "R2.C3"],
        "mode": "uncapped (no_cap()); the sampled law is the full bi-Kappa distribution",
        "environment": {
            "platform": platform.platform(),
            "machine": platform.machine(),
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "cxx": os.popen("g++ --version").read().splitlines()[0],
            "cxx_target": next((ln for ln in os.popen("g++ --version").read().splitlines()
                                if ln.startswith("Target:")), ""),
            "cxxflags": "-Wall -Wextra -std=c++11 -O2",
            "stdlib": "libc++ (Apple clang)",
            "rng": "std::mt19937, seeded explicitly per run",
            "git_commit": subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True,
                                         text=True).stdout.strip(),
            "git_dirty": bool(subprocess.run(["git", "status", "--porcelain"],
                                             capture_output=True, text=True).stdout.strip()),
            # Pin the exact sampler that produced raw/ by content rather than by commit,
            # so the result stays traceable even when the working tree is dirty for
            # reasons outside this experiment.
            "sampler_header_sha256": subprocess.run(
                ["shasum", "-a", "256", "../../cpp/bi_kappa_distribution.H"],
                capture_output=True, text=True).stdout.split()[0],
        },
        "summary": summary,
        "frame_invariance": frame_checks,
        "normalized_vector_invariance": vector_checks,
        "per_run": per_run,
    }
    with open(os.path.join(results_dir, "exp1_results.json"), "w") as fh:
        json.dump(out, fh, indent=2)

    write_table(summary, frame_checks, os.path.join(results_dir, "exp1_table.md"), vector_checks)
    print(f"\nwrote {results_dir}/exp1_results.json and {results_dir}/exp1_table.md")
    return 0


def write_table(summary, frame_checks, path, vector_checks=()):
    lines = [
        "# Experiment 1 results",
        "",
        "Uncapped sampler (`no_cap()`), 5 replicate seeds x 10^5 samples per configuration,",
        "double precision, Apple clang / libc++, `std::mt19937`.",
        "",
        "`sqrt(n)D` is the KS statistic scaled by sqrt(n). It is O(1) -- around 0.87 on average --",
        "when the law is correct, and does not shrink with n, unlike a p-value, which collapses to",
        "0 for any large enough sample regardless of how small the discrepancy is.",
        "Entries are mean +/- sd over the 5 replicates.",
        "",
        "## Block A -- radial law and directional uniformity (theta_perp:theta_par = 1:2, B || z)",
        "",
        "| kappa | non-finite / 5x10^5 | radial sqrt(n)D | CvM | cos(theta) sqrt(n)D | phi sqrt(n)D "
        "| indep. chi2 p | direction cells chi2 p | max |dlog R| |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    for e in summary:
        if e["block"] != "A":
            continue
        lines.append(
            f"| {e['kappa']:g} | {e['total_nonfinite']} "
            f"| {e['radial_ks_sqrtn']['mean']:.3f} +/- {e['radial_ks_sqrtn']['sd']:.3f} "
            f"| {e['radial_cvm_stat']['mean']:.3f} "
            f"| {e['cos_theta_ks_sqrtn']['mean']:.3f} +/- {e['cos_theta_ks_sqrtn']['sd']:.3f} "
            f"| {e['phi_ks_sqrtn']['mean']:.3f} +/- {e['phi_ks_sqrtn']['sd']:.3f} "
            f"| {e['independence_chi2_pvalue']['mean']:.3f} "
            f"| {e['direction_cells_chi2_pvalue']['mean']:.3f} "
            f"| {e['max_tail_quantile_log_err']:.3f} |"
        )

    lines += [
        "",
        f"## Cell-count test (the manuscript's test), replicates pooled",
        "",
        f"Pearson chi^2 against equal occupancy of {N_SHELLS} equal-probability radius shells x "
        f"{N_COS_BINS * N_PHI_BINS} equal-solid-angle",
        f"direction cells ({N_COS_BINS} cos(theta) x {N_PHI_BINS} phi intervals). `radius` sums "
        "over direction and tests the radial law alone;",
        "`all cells` tests radius, direction and their independence together.",
        "",
        "| block | kappa | theta_par | B direction | draws | expected per cell | min cell | radius p | all cells p |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    for e in summary:
        lines.append(
            f"| {e['block']} | {e['kappa']:g} | {e['theta_par']:g} | {e['ub_label']} "
            f"| {e['n_replicates'] * e['n_per_replicate']} | {e['cells_pooled_expected']:.0f} "
            f"| {e['cells_pooled_min_count']} | {e['cells_pooled_radius_pvalue']:.3f} "
            f"| {e['cells_pooled_all_pvalue']:.3f} |"
        )

    lines += [
        "",
        "## Orientation of the bounded radial diagnostic",
        "",
        "Both W = 1/(1+T) ~ Beta(kappa-1/2, 3/2) and Y = T/(1+T) ~ Beta(3/2, kappa-1/2) are exact",
        "bijections of T and carry identical information in exact arithmetic. In doubles they do",
        "not. The table above uses W; the same test on Y is shown here for comparison.",
        "",
        "| kappa | radial sqrt(n)D using W | radial sqrt(n)D using Y | Y values rounded to exactly 1 "
        "| W values underflowed to 0 |",
        "|---|---|---|---|---|",
    ]
    for e in summary:
        if e["block"] != "A":
            continue
        n_tot = e["n_replicates"] * e["n_per_replicate"]
        lines.append(
            f"| {e['kappa']:g} | {e['radial_ks_sqrtn']['mean']:.3f} "
            f"| {e['radial_ks_sqrtn_naive_Y']['mean']:.3f} "
            f"| {e['total_y_saturated']} / {n_tot} "
            f"({100.0 * e['total_y_saturated'] / n_tot:.2f}%) "
            f"| {e['total_w_underflow']} / {n_tot} |"
        )

    lines += [
        "",
        "## Block B -- anisotropy and arbitrary field frame",
        "",
        "| kappa | theta_par | B direction | radial sqrt(n)D | MAD ratio (expect theta_par/theta_perp) "
        "| basis orthonormality err | max off-diag corr |",
        "|---|---|---|---|---|---|---|",
    ]
    for e in summary:
        if e["block"] != "B":
            continue
        offdiag = (f"{e['max_offdiag_corr']['mean']:.2e}" if "max_offdiag_corr" in e else "n/a")
        lines.append(
            f"| {e['kappa']:g} | {e['theta_par']:g} | {e['ub_label']} "
            f"| {e['radial_ks_sqrtn']['mean']:.3f} +/- {e['radial_ks_sqrtn']['sd']:.3f} "
            f"| {e['mad_ratio_par_perp']['mean']:.4f} (expect {e['theta_par'] / e['theta_perp']:g}) "
            f"| {e['basis_orthonormality_err']['max']:.1e} | {offdiag} |"
        )

    worst = max(frame_checks, key=lambda c: c["max_rel_radius_diff"]) if frame_checks else None
    n_ident = sum(c["n_bitwise_identical"] for c in frame_checks)
    n_cmp = sum(c["n_compared"] for c in frame_checks)
    lines += [
        "",
        "## Frame invariance (direct comparison)",
        "",
        "For matched (kappa, theta, seed), the normalized radius recovered after rotating into an",
        "arbitrary field direction is compared against the axis-aligned run draw by draw.",
        "",
        f"- comparisons: {len(frame_checks)} run pairs, {n_cmp} draws",
        f"- bitwise identical radii: {n_ident} / {n_cmp} ({100.0 * n_ident / n_cmp:.4f}%)",
        f"- largest relative difference: {worst['max_rel_radius_diff']:.2e} "
        f"(kappa={worst['kappa']:g}, {worst['direction']}, seed={worst['seed']})",
        "",
        "Whole normalized velocity: every block-B run against the (theta_par = 2, B || z) run",
        "with the same kappa and seed, after undoing each run's own scaling and rotation.",
        "",
        f"- comparisons: {len(vector_checks)} run pairs, "
        f"{sum(c['n_compared'] for c in vector_checks)} draws",
        f"- largest relative difference |u - u_ref| / |u_ref|: "
        f"{max((c['max_rel_vector_diff'] for c in vector_checks), default=float('nan')):.2e}",
        "",
        "Notes.",
        "",
        f"- Independence: chi^2 contingency of {N_RADIAL_BINS} radial quartiles x "
        f"{N_COS_BINS * N_PHI_BINS} equal-solid-angle direction",
        f"  cells ({N_COS_BINS} cos(theta) x {N_PHI_BINS} phi intervals). The direction-cell "
        "column tests the same cell",
        "  counts, summed over radius, against equal occupancy: the joint direction law.",
        "- Variance-based checks are reported only for kappa > 3/2, where the second moment",
        "  exists. Below that the MAD ratio and the radial/quantile diagnostics carry the test.",
        "- `max off-diag corr` is the largest off-diagonal correlation of the sample covariance",
        "  in the recovered field-aligned frame; a misaligned rotation would inflate it. It is",
        "  itself noisy at kappa = 2, where the fourth moment does not exist and the correlation",
        "  estimator therefore has no finite variance.",
        "",
    ]
    with open(path, "w") as fh:
        fh.write("\n".join(lines))


if __name__ == "__main__":
    raise SystemExit(main())
