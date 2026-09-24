#!/usr/bin/env python3
"""Experiment 1 large-sample cell test: the cell-count test at all nine kappa on 5 x 10^7 draws each.

This is the paper's cell-count test (Table III, label "table:validation").  At each of the
nine kappa in KAPPAS it draws 100 seeds x 5 x 10^5 = 5 x 10^7 velocities,
(theta_perp, theta_par) = (1, 2), field along z, uncapped, double precision, and applies the
Pearson chi^2 test over the 10 radius shells alone and over all 400 cells
(exp1_common.log_shell_edges, cell_pvalues).  The kappa at index i of KAPPAS uses seeds
100001 + 1000 i ... 100100 + 1000 i (seeds_for).

Draws with a component beyond the largest double (expected only near kappa = 1/2, at a
rate of about 7 x 10^-7 at kappa = 0.51) have R beyond every shell edge.  They are
counted in the outermost shell of the radius test; their direction cannot be recovered,
so the 400-cell test uses the finite draws only.  Draws with an undefined (NaN) component
would indicate a defect and are counted separately.

It also checks the scatter between seeds: if successive draws were correlated, or the
counts otherwise more variable than a multinomial allows, the per-seed chi^2 would
average above its degrees of freedom and the per-seed p-values would not be uniform.

Writes results/exp1_cells.json and results/exp1_cells.md.  The raw draws of each seed
are written to a temporary file, counted, and deleted; only the counts are kept.

Usage:  uv run --project ../../python python exp1_cells.py
"""

from __future__ import annotations

import json
import os
import platform
import subprocess
import sys
import tempfile

import numpy as np
import scipy
from scipy import stats

from exp1_common import (N_COS_BINS, N_PHI_BINS, N_SHELLS, cell_pvalues,
                          log_shell_edges)

KAPPAS = [0.51, 0.55, 0.75, 1.0, 1.25, 1.5, 2.0, 5.0, 10.0]
N_SEEDS = 100
N_PER_SEED = 500_000
EXE = "./exp1_draw.exe"


def seeds_for(i: int) -> range:
    return range(100001 + 1000 * i, 100001 + 1000 * i + N_SEEDS)


def cell_counts(v: np.ndarray, kappa: float):
    """Cell counts for (1, 2) with the field along z, plus the overflow and NaN tallies.

    Returns (counts[N_SHELLS, 40] of the finite draws, n_overflow, n_nan).
    """
    u = v / (np.sqrt(kappa) * np.array([1.0, 1.0, 2.0]))
    nan = np.isnan(u).any(axis=1)
    finite = np.isfinite(u).all(axis=1)
    overflow = ~finite & ~nan
    u = u[finite]
    R = np.hypot(np.hypot(u[:, 0], u[:, 1]), u[:, 2])
    cos_theta = u[:, 2] / R
    phi = np.arctan2(u[:, 1], u[:, 0])
    c_bin = np.clip(np.floor((cos_theta + 1.0) / 2.0 * N_COS_BINS).astype(int), 0, N_COS_BINS - 1)
    p_bin = np.clip(np.floor((phi + np.pi) / (2.0 * np.pi) * N_PHI_BINS).astype(int), 0, N_PHI_BINS - 1)
    shell = np.searchsorted(log_shell_edges(kappa), np.log(R), side="right")
    counts = np.zeros((N_SHELLS, N_COS_BINS * N_PHI_BINS), dtype=np.int64)
    np.add.at(counts, (shell, c_bin * N_PHI_BINS + p_bin), 1)
    return counts, int(overflow.sum()), int(nan.sum())


def run_kappa(i: int, kappa: float, tmp: str) -> dict:
    path = os.path.join(tmp, "draws.bin")
    per_seed, n_overflow, n_nan = [], 0, 0
    per_seed_overflow = []
    for seed in seeds_for(i):
        subprocess.run([EXE, repr(kappa), str(seed), str(N_PER_SEED), path], check=True)
        v = np.fromfile(path, dtype=np.float64).reshape(-1, 3)
        os.remove(path)
        counts, over, bad = cell_counts(v, kappa)
        per_seed.append(counts)
        per_seed_overflow.append(over)
        n_overflow += over
        n_nan += bad
    per_seed = np.array(per_seed)
    pooled = per_seed.sum(axis=0)

    # Radius test: overflowing draws lie beyond every shell edge, i.e. in the outermost shell.
    shell = pooled.sum(axis=1)
    shell[-1] += n_overflow
    radius_p = float(stats.chisquare(shell).pvalue)
    _, all_p = cell_pvalues(pooled)

    # Scatter between seeds, from the radius shells of each seed (9 degrees of freedom).
    seed_shells = per_seed.sum(axis=2)
    seed_shells[:, -1] += np.array(per_seed_overflow)
    chi2 = np.array([stats.chisquare(s).statistic for s in seed_shells])
    seed_p = stats.chi2.sf(chi2, N_SHELLS - 1)
    n = int(shell.sum())
    s = seeds_for(i)
    rec = {
        "kappa": kappa, "seeds": [s.start, s.stop - 1], "n_per_seed": N_PER_SEED,
        "n_draws": n, "n_overflow": n_overflow, "n_nan": n_nan,
        "expected_per_cell": n / pooled.size,
        "radius_pvalue": radius_p, "all_cells_pvalue": all_p,
        "shell_counts": shell.tolist(),
        "shell_z": ((shell - n / N_SHELLS) / np.sqrt(n * (1 / N_SHELLS) * (1 - 1 / N_SHELLS))).tolist(),
        "per_seed_radius_chi2_over_dof_mean": float(chi2.mean() / (N_SHELLS - 1)),
        "per_seed_radius_chi2_over_dof_sem": float(chi2.std(ddof=1) / (N_SHELLS - 1) / np.sqrt(len(chi2))),
        "per_seed_pvalue_uniformity_ks_pvalue": float(stats.kstest(seed_p, "uniform").pvalue),
        "per_seed_fraction_p_below_0.05": float((seed_p < 0.05).mean()),
        "pooled_cell_counts": pooled.tolist(),
        "per_seed_shell_counts": seed_shells.tolist(),
    }
    print(f"kappa={kappa:g}: radius p={radius_p:.3f}, all cells p={all_p:.3f}, "
          f"overflow={n_overflow}, nan={n_nan}", flush=True)
    return rec


def main() -> int:
    if not os.path.exists(EXE):
        print(f"{EXE} not found; run `make cells`", file=sys.stderr)
        return 1
    with tempfile.TemporaryDirectory() as tmp:
        records = [run_kappa(i, k, tmp) for i, k in enumerate(KAPPAS)]
    out = {
        "experiment": "Experiment 1 large-sample cell test -- nine kappa, 5 x 10^7 draws each",
        "theta_perp": 1.0, "theta_par": 2.0, "field": "z", "mode": "uncapped",
        "precision": "double", "n_shells": N_SHELLS, "n_direction_cells": N_COS_BINS * N_PHI_BINS,
        "summary": records,
        "environment": {
            "platform": platform.platform(), "python": platform.python_version(),
            "numpy": np.__version__, "scipy": scipy.__version__,
            "sampler_header_sha256": subprocess.run(
                ["shasum", "-a", "256", "../../cpp/bi_kappa_distribution.H"],
                capture_output=True, text=True).stdout.split()[0],
            "git_commit": subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True,
                                         text=True).stdout.strip(),
        },
    }
    os.makedirs("results", exist_ok=True)
    with open("results/exp1_cells.json", "w") as fh:
        json.dump(out, fh, indent=1)

    lines = [
        "# Experiment 1 large-sample cell test -- nine kappa, 5 x 10^7 draws each",
        "",
        f"{N_SEEDS} seeds x {N_PER_SEED} draws per kappa (kappa index i uses seeds 100001 + 1000 i "
        "onward); (theta_perp, theta_par) = (1, 2), B || z, "
        "uncapped, double precision. Expected per cell: "
        f"{records[0]['expected_per_cell']:.0f}. Overflow = draws with a component beyond "
        "the largest double (counted in the outermost shell of the radius test, excluded "
        "from the 400-cell test).",
        "",
        "| kappa | radius p (10 shells) | all cells p (400) | overflow | NaN | "
        "seed chi^2/dof | seed KS p |",
        "|---|---|---|---|---|---|---|",
    ]
    for r in records:
        lines.append(
            f"| {r['kappa']:g} | {r['radius_pvalue']:.3f} | {r['all_cells_pvalue']:.3f} | "
            f"{r['n_overflow']} | {r['n_nan']} | "
            f"{r['per_seed_radius_chi2_over_dof_mean']:.3f} +/- "
            f"{r['per_seed_radius_chi2_over_dof_sem']:.3f} | "
            f"{r['per_seed_pvalue_uniformity_ks_pvalue']:.3f} |")
    lines += ["", "Shell deviations in binomial standard deviations:", ""]
    for r in records:
        lines.append(f"- kappa = {r['kappa']:g}: " + " ".join(f"{z:+.1f}" for z in r["shell_z"]))
    lines.append("")
    with open("results/exp1_cells.md", "w") as fh:
        fh.write("\n".join(lines))
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
