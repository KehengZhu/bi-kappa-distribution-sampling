#!/usr/bin/env python3
"""Experiment 1 large-sample variances: the field-aligned second moments at kappa = 2, 5, 10.

This is the paper's variance table (Table IV, label "table:moments").  It redraws, seed for seed, the
large-sample cell-test runs of exp1_cells.py at the three kappa > 3/2 where the variance
exists (100 seeds x 5 x 10^5 = 5 x 10^7 draws per kappa, (theta_perp, theta_par) = (1, 2),
field along z, uncapped, double precision), and records each run's sample variance of the
three components.  The table reports, per component, the mean of the 100 run variances and
its standard error (standard deviation of the run variances / sqrt(100)).

The deviation of the mean from theory is compared with that standard error, not with the
spread of single runs.  At kappa = 2 the fourth moment does not exist, so the sample
variance converges only as about N^(-1/3) and a large sample is needed.  Using the cell-test
seeds makes the variance sample the same draws the cell test counted; the script checks that
by recounting the cells and comparing with results/exp1_cells.json.

Writes results/exp1_moments.json and results/exp1_moments.md.  The raw draws are deleted
after use.

Usage:  uv run --project ../../python python exp1_moments.py   (after `make exp1_draw.exe`)
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

from exp1_cells import EXE, KAPPAS, N_PER_SEED, cell_counts, seeds_for

MOMENT_KAPPAS = [2.0, 5.0, 10.0]
THETA = np.array([1.0, 1.0, 2.0])
NAMES = ["v_perp1", "v_perp2", "v_par"]


def main() -> int:
    if not os.path.exists(EXE):
        print(f"{EXE} not found; run `make exp1_draw.exe`", file=sys.stderr)
        return 1
    with open("results/exp1_cells.json") as fh:
        cells = {r["kappa"]: r for r in json.load(fh)["summary"]}
    records = []
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "draws.bin")
        for kappa in MOMENT_KAPPAS:
            i = KAPPAS.index(kappa)
            run_var, pooled = [], np.zeros((10, 40), dtype=np.int64)
            for seed in seeds_for(i):
                subprocess.run([EXE, repr(kappa), str(seed), str(N_PER_SEED), path], check=True)
                v = np.fromfile(path, dtype=np.float64).reshape(-1, 3)
                os.remove(path)
                if not np.isfinite(v).all():
                    raise RuntimeError(f"non-finite draw at kappa={kappa}, seed={seed}")
                run_var.append(np.var(v, axis=0, ddof=1))
                counts, _, _ = cell_counts(v, kappa)
                pooled += counts
            same_draws = pooled.tolist() == cells[kappa]["pooled_cell_counts"]
            run_var = np.array(run_var)
            expected = THETA**2 * kappa / (2.0 * kappa - 3.0)
            mean = run_var.mean(axis=0)
            sd = run_var.std(axis=0, ddof=1)
            se = sd / np.sqrt(len(run_var))
            rec = {
                "kappa": kappa, "seeds": [seeds_for(i).start, seeds_for(i).stop - 1],
                "n_runs": len(run_var), "n_per_run": N_PER_SEED,
                "same_draws_as_cell_test": same_draws,
                "components": {
                    n: {"expected": float(expected[j]), "mean_run_variance": float(mean[j]),
                        "sd_run_variance": float(sd[j]), "standard_error": float(se[j]),
                        "diff_percent": float(100 * (mean[j] - expected[j]) / expected[j]),
                        "diff_over_se": float((mean[j] - expected[j]) / se[j])}
                    for j, n in enumerate(NAMES)},
                "run_variances": run_var.tolist(),
            }
            records.append(rec)
            print(f"kappa={kappa:g}: same draws as cell test: {same_draws}", flush=True)
    out = {
        "experiment": "Experiment 1 large-sample variances -- kappa = 2, 5, 10, 5 x 10^7 draws each",
        "theta_perp": 1.0, "theta_par": 2.0, "field": "z", "mode": "uncapped",
        "precision": "double", "summary": records,
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
    with open("results/exp1_moments.json", "w") as fh:
        json.dump(out, fh, indent=1)
    lines = [
        "# Experiment 1 large-sample variances -- kappa = 2, 5, 10",
        "",
        f"100 runs x {N_PER_SEED} draws per kappa, the exp1_cells.py seeds; (theta_perp, "
        "theta_par) = (1, 2), B || z, uncapped, double precision. Mean of the run variances "
        "+/- standard error (sd / sqrt(100)); sd = spread of single runs.",
        "",
        "| kappa | component | expected | mean +/- SE | sd of runs | diff % | diff / SE |",
        "|---|---|---|---|---|---|---|",
    ]
    for r in records:
        for n, c in r["components"].items():
            lines.append(
                f"| {r['kappa']:g} | {n} | {c['expected']:.4f} | {c['mean_run_variance']:.4f} "
                f"+/- {c['standard_error']:.4f} | {c['sd_run_variance']:.4f} | "
                f"{c['diff_percent']:+.2f} | {c['diff_over_se']:+.2f} |")
    lines += ["", "Same draws as the cell test: "
              + ", ".join(f"kappa = {r['kappa']:g}: {r['same_draws_as_cell_test']}" for r in records),
              ""]
    with open("results/exp1_moments.md", "w") as fh:
        fh.write("\n".join(lines))
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
