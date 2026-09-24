#!/usr/bin/env python3
"""Experiment 1 large-sample marginals: component histograms at kappa = 1, 2, 10.

This supplies the paper's component-distribution figure (Fig. 5).  It redraws the exp1_cells.py
runs at kappa = 1, 2 and 10 (100 seeds x 5 x 10^5 = 5 x 10^7 draws per kappa,
(theta_perp, theta_par) = (1, 2), field along z) and histograms each component divided by
its thermal speed, v_i / theta_i, on the figure's bins: 48 bins of width 0.25 on [-6, 6].
Counts outside [-6, 6] are recorded but not drawn.  The draws are the cell-test draws; the
script checks that by recounting the cells and comparing with results/exp1_cells.json.

Writes results/exp1_marginals.json.  Raw draws are deleted after use.

Usage:  uv run --project ../../python python exp1_marginals.py   (after `make exp1_draw.exe`)
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile
from concurrent.futures import ProcessPoolExecutor

import numpy as np

from exp1_cells import EXE, KAPPAS, N_PER_SEED, cell_counts, seeds_for

MARGINAL_KAPPAS = [1.0, 2.0, 10.0]
THETA = np.array([1.0, 1.0, 2.0])
BINS = np.linspace(-6.0, 6.0, 49)


def one_seed(args):
    kappa, seed = args
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "draws.bin")
        subprocess.run([EXE, repr(kappa), str(seed), str(N_PER_SEED), path], check=True)
        v = np.fromfile(path, dtype=np.float64).reshape(-1, 3)
    if not np.isfinite(v).all():
        raise RuntimeError(f"non-finite draw at kappa={kappa}, seed={seed}")
    s = v / THETA
    hist = np.array([np.histogram(s[:, j], bins=BINS)[0] for j in range(3)])
    counts, _, _ = cell_counts(v, kappa)
    return kappa, hist, counts, len(v)


def main() -> int:
    if not os.path.exists(EXE):
        print(f"{EXE} not found; run `make exp1_draw.exe`", file=sys.stderr)
        return 1
    with open("results/exp1_cells.json") as fh:
        cells = {r["kappa"]: r for r in json.load(fh)["summary"]}
    jobs = [(k, s) for k in MARGINAL_KAPPAS for s in seeds_for(KAPPAS.index(k))]
    acc = {k: {"hist": np.zeros((3, len(BINS) - 1), dtype=np.int64),
               "cells": np.zeros((10, 40), dtype=np.int64), "n": 0} for k in MARGINAL_KAPPAS}
    with ProcessPoolExecutor(max_workers=max(1, (os.cpu_count() or 2) // 2)) as pool:
        for kappa, hist, counts, n in pool.map(one_seed, jobs):
            acc[kappa]["hist"] += hist
            acc[kappa]["cells"] += counts
            acc[kappa]["n"] += n
    records = []
    for k in MARGINAL_KAPPAS:
        same = acc[k]["cells"].tolist() == cells[k]["pooled_cell_counts"]
        print(f"kappa={k:g}: n={acc[k]['n']}, same draws as cell test: {same}")
        records.append({"kappa": k, "n_total": acc[k]["n"], "same_draws_as_cell_test": same,
                        "seeds": [seeds_for(KAPPAS.index(k)).start,
                                  seeds_for(KAPPAS.index(k)).stop - 1],
                        "counts": {"v_perp1": acc[k]["hist"][0].tolist(),
                                   "v_perp2": acc[k]["hist"][1].tolist(),
                                   "v_par": acc[k]["hist"][2].tolist()}})
    out = {"experiment": "Experiment 1 large-sample marginals -- kappa = 1, 2, 10",
           "bin_edges": BINS.tolist(), "variable": "v_i / theta_i",
           "theta_perp": 1.0, "theta_par": 2.0, "field": "z", "summary": records,
           "sampler_header_sha256": subprocess.run(
               ["shasum", "-a", "256", "../../cpp/bi_kappa_distribution.H"],
               capture_output=True, text=True).stdout.split()[0]}
    with open("results/exp1_marginals.json", "w") as fh:
        json.dump(out, fh, indent=1)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
