#!/usr/bin/env python3
"""Experiment 1 tail extension of the cell-count test: the radius test resolved into the tail.

The ten equal-probability shells of exp1_cells.py leave the fastest 10% of draws in one
open shell, R > r_9, so the test counts those draws but does not resolve their radii.
This run divides that shell further at the radii that R exceeds with probability
10^-2, 10^-3, 10^-4 and 10^-5, giving 14 radius shells (13 degrees of freedom), and
applies Pearson's chi^2 test with the exact shell probabilities.

The draws are the same as the cell test's: the same executable, kappa values and seed
blocks (exp1_cells.seeds_for), so every test in the paper's validation section uses one
set of 5 x 10^7 draws per kappa.  As a check that the draws are the same, the 400-cell
counts and the overflow tallies are recomputed and compared, cell by cell, with
results/exp1_cells.json; any difference stops the run.

Passing kappa values on the command line runs only those kappa and writes no results.

Tail edges.  With W = 1/(1+R^2) ~ Beta(kappa-1/2, 3/2), P(R > r) = P(W < w) = I_w(a, 3/2),
a = kappa - 1/2.  Near kappa = 1/2 the edge w is far below the smallest double
(w ~ 10^-500 for the 10^-5 edge at kappa = 0.51), so scipy's Beta quantile cannot return
it.  There the leading term of the series I_w(a, b) = w^a / (a B(a, b)) (1 + O(w)) is
exact to double precision, and log w = (log t + log a + log B(a, b)) / a.  Every edge
is checked against the regularized incomplete Beta function where w is representable.

Draws with a component beyond the largest double have R above every edge (their
components exceed 10^308, the 10^-5 edge is at most ~10^250), so the radius test
counts them in the outermost shell, as the cell test does.

Writes results/exp1_tail.json and results/exp1_tail.md.

Usage:  uv run --project ../../python python exp1_tail.py
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
from scipy import special, stats

from exp1_common import log_shell_edges
from exp1_cells import EXE, KAPPAS, N_PER_SEED, N_SEEDS, cell_counts, seeds_for

TAIL_EXCEEDANCE = (1e-2, 1e-3, 1e-4, 1e-5)
B = 1.5


def log_tail_edges(kappa: float) -> np.ndarray:
    """log of the radii that R exceeds with the probabilities in TAIL_EXCEEDANCE."""
    a = kappa - 0.5
    out = []
    for t in TAIL_EXCEEDANCE:
        w = float(stats.beta.ppf(t, a, B))
        if w > 1e-280:
            # Representable: check the quantile against the CDF directly.
            if abs(special.betainc(a, B, w) / t - 1.0) > 1e-9:
                raise RuntimeError(f"Beta quantile inaccurate at kappa={kappa}, t={t}")
            log_w = np.log(w)
        else:
            log_w = (np.log(t) + np.log(a) + special.betaln(a, B)) / a
        out.append(0.5 * (np.log1p(-np.exp(log_w)) - log_w))
    return np.array(out)


def check_asymptote() -> None:
    """The series form must agree with the Beta quantile where both apply."""
    for kappa in (0.51, 0.55):
        a = kappa - 0.5
        for t in (1e-3, 1e-2):
            w = float(stats.beta.ppf(t, a, B))
            if not w > 1e-300:
                continue
            log_w = (np.log(t) + np.log(a) + special.betaln(a, B)) / a
            if abs(log_w - np.log(w)) > 1e-8 * abs(np.log(w)):
                raise RuntimeError(f"series form disagrees at kappa={kappa}, t={t}")


def radius_edges(kappa: float) -> tuple[np.ndarray, np.ndarray]:
    """log edges and probabilities of the 14 radius shells."""
    edges = np.concatenate([log_shell_edges(kappa), log_tail_edges(kappa)])
    if not np.all(np.diff(edges) > 0):
        raise RuntimeError(f"shell edges not increasing at kappa={kappa}")
    exceed = np.concatenate([1.0 - np.arange(10) / 10.0, TAIL_EXCEEDANCE, [0.0]])
    # exceed = [1, 0.9, ..., 0.1, 1e-2, ..., 1e-5, 0]; shell j lies between exceed[j], exceed[j+1].
    probs = exceed[:-1] - exceed[1:]
    return edges, probs


def run_kappa(i: int, kappa: float, tmp: str, ref: dict) -> dict:
    path = os.path.join(tmp, "draws.bin")
    edges, probs = radius_edges(kappa)
    n_shells = len(probs)
    per_seed, pooled_cells, n_overflow, n_nan = [], None, 0, 0
    for seed in seeds_for(i):
        subprocess.run([EXE, repr(kappa), str(seed), str(N_PER_SEED), path], check=True)
        v = np.fromfile(path, dtype=np.float64).reshape(-1, 3)
        os.remove(path)
        cells, over, bad = cell_counts(v, kappa)
        pooled_cells = cells if pooled_cells is None else pooled_cells + cells
        n_overflow += over
        n_nan += bad

        u = v / (np.sqrt(kappa) * np.array([1.0, 1.0, 2.0]))
        u = u[np.isfinite(u).all(axis=1)]
        R = np.hypot(np.hypot(u[:, 0], u[:, 1]), u[:, 2])
        shell = np.searchsorted(edges, np.log(R), side="right")
        counts = np.bincount(shell, minlength=n_shells)
        counts[-1] += over
        per_seed.append(counts)

    # The draws must be exactly those of the cell test.
    if (pooled_cells.tolist() != ref["pooled_cell_counts"] or n_overflow != ref["n_overflow"]
            or n_nan != ref["n_nan"]):
        raise RuntimeError(f"kappa={kappa}: draws differ from results/exp1_cells.json")

    per_seed = np.array(per_seed)
    shell = per_seed.sum(axis=0)
    n = int(shell.sum())
    expected = n * probs
    chi2, p = stats.chisquare(shell, expected)
    z = (shell - expected) / np.sqrt(n * probs * (1.0 - probs))
    s = seeds_for(i)
    print(f"kappa={kappa:g}: 14-shell radius p={p:.3f}, max|z|={np.abs(z).max():.2f}, "
          f"tail counts={shell[9:].tolist()}", flush=True)
    return {
        "kappa": kappa, "seeds": [s.start, s.stop - 1], "n_per_seed": N_PER_SEED,
        "n_draws": n, "n_overflow": n_overflow, "n_nan": n_nan,
        "log_edges": edges.tolist(), "log10_edges": (edges / np.log(10)).tolist(),
        "shell_probabilities": probs.tolist(),
        "shell_counts": shell.tolist(), "shell_expected": expected.tolist(),
        "shell_z": z.tolist(),
        "radius_chi2": float(chi2), "radius_dof": n_shells - 1, "radius_pvalue": float(p),
        "cells_match_exp1_cells": True,
        "per_seed_shell_counts": per_seed.tolist(),
    }


def main() -> int:
    if not os.path.exists(EXE):
        print(f"{EXE} not found; run `make cells`", file=sys.stderr)
        return 1
    check_asymptote()
    with open("results/exp1_cells.json") as fh:
        ref_all = json.load(fh)
    ref = {r["kappa"]: r for r in ref_all["summary"]}
    only = [float(x) for x in sys.argv[1:]]
    with tempfile.TemporaryDirectory() as tmp:
        records = [run_kappa(i, k, tmp, ref[k]) for i, k in enumerate(KAPPAS)
                   if not only or k in only]
    if only:
        return 0
    header = subprocess.run(["shasum", "-a", "256", "../../cpp/bi_kappa_distribution.H"],
                            capture_output=True, text=True).stdout.split()[0]
    if header != ref_all["environment"]["sampler_header_sha256"]:
        print("warning: sampler header differs from the cell test's", file=sys.stderr)
    out = {
        "experiment": "Experiment 1 tail extension -- 14 radius shells, nine kappa, "
                      "the cell test's 5 x 10^7 draws each",
        "theta_perp": 1.0, "theta_par": 2.0, "field": "z", "mode": "uncapped",
        "precision": "double", "tail_exceedance": list(TAIL_EXCEEDANCE),
        "summary": records,
        "environment": {
            "platform": platform.platform(), "python": platform.python_version(),
            "numpy": np.__version__, "scipy": scipy.__version__,
            "sampler_header_sha256": header,
            "git_commit": subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True,
                                         text=True).stdout.strip(),
        },
    }
    with open("results/exp1_tail.json", "w") as fh:
        json.dump(out, fh, indent=1)

    lines = [
        "# Experiment 1 tail extension -- 14 radius shells",
        "",
        "The cell test's draws (same seeds; 400-cell counts reproduced exactly). The "
        "outermost of the ten equal-probability shells is divided at P(R > r) = "
        + ", ".join(f"{t:g}" for t in TAIL_EXCEEDANCE)
        + "; Pearson chi^2 with the exact shell probabilities, 13 degrees of freedom.",
        "",
        "| kappa | radius p (14 shells) | log10 r at 1e-5 | count beyond 1e-5 (E = "
        f"{records[0]['shell_expected'][-1]:.0f}) | overflow |",
        "|---|---|---|---|---|",
    ]
    for r in records:
        lines.append(f"| {r['kappa']:g} | {r['radius_pvalue']:.3f} | "
                     f"{r['log10_edges'][-1]:.1f} | {r['shell_counts'][-1]} | "
                     f"{r['n_overflow']} |")
    lines += ["", "Shell deviations in binomial standard deviations (14 shells):", ""]
    for r in records:
        lines.append(f"- kappa = {r['kappa']:g}: " + " ".join(f"{z:+.1f}" for z in r["shell_z"]))
    lines.append("")
    with open("results/exp1_tail.md", "w") as fh:
        fh.write("\n".join(lines))
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
