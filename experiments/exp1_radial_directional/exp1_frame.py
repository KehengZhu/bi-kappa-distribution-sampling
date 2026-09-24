#!/usr/bin/env python3
"""Experiment 1 large-sample frame invariance: kappa = 0.55, 2, 10 on the cell-test seeds.

This is the paper's scaling-and-rotation invariance test (Sec. VI C).  For each of the 100
exp1_cells.py seeds at kappa = 0.55, 2 and 10 it draws 5 x 10^5 velocities for each of the
six combinations of (theta_perp, theta_par) in {(1, 1), (1, 2)} and field directions z,
(1,1,1)/sqrt(3) and (0.3,-0.5,0.8)/|.|.  The (1, 2), z run is the cell-test run itself.
Each of the other five runs is mapped back to the normalized velocity u' = A^-1 Q^T v
(exp1_common.field_basis, the C++ basis construction) and compared, draw by draw, with u
of the (1, 2), z run at the same seed: max |u' - u| / |u| over all draws.  The test thus
uses the same 5 x 10^7 draws per kappa as the other validation tests.

Writes results/exp1_frame.json and results/exp1_frame.md.  Raw draws are deleted after use.

Usage:  uv run --project ../../python python exp1_frame.py   (after `make exp1_draw.exe`)
"""

from __future__ import annotations

import json
import os
import platform
import subprocess
import sys
import tempfile
from concurrent.futures import ProcessPoolExecutor

import numpy as np

from exp1_common import field_basis, is_axis_aligned_z
from exp1_cells import EXE, KAPPAS, N_PER_SEED, seeds_for

FRAME_KAPPAS = [0.55, 2.0, 10.0]


def unit(x: float, y: float, z: float) -> np.ndarray:
    n = np.sqrt(x * x + y * y + z * z)
    return np.array([x / n, y / n, z / n])


DIRECTIONS = {"z": np.array([0.0, 0.0, 1.0]), "diag111": unit(1.0, 1.0, 1.0),
              "oblique": unit(0.3, -0.5, 0.8)}
THETAS = [(1.0, 1.0), (1.0, 2.0)]
REFERENCE = ((1.0, 2.0), "z")


def normalized(v: np.ndarray, kappa: float, theta, ub: np.ndarray) -> np.ndarray:
    local = v if is_axis_aligned_z(ub) else v @ field_basis(ub)
    return local / (np.sqrt(kappa) * np.array([theta[0], theta[0], theta[1]]))


def one_seed(args):
    kappa, seed = args
    out = {}
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "draws.bin")
        u = {}
        for theta in THETAS:
            for label, ub in DIRECTIONS.items():
                cmd = [EXE, repr(kappa), str(seed), str(N_PER_SEED), path,
                       repr(theta[0]), repr(theta[1])] + [repr(float(c)) for c in ub]
                subprocess.run(cmd, check=True)
                v = np.fromfile(path, dtype=np.float64).reshape(-1, 3)
                os.remove(path)
                u[(theta, label)] = normalized(v, kappa, theta, ub)
    ref = u[REFERENCE]
    for key, got in u.items():
        if key == REFERENCE:
            continue
        m = np.isfinite(ref).all(axis=1) & np.isfinite(got).all(axis=1)
        ref_norm = np.hypot(np.hypot(ref[m, 0], ref[m, 1]), ref[m, 2])
        d = got[m] - ref[m]
        rel = np.hypot(np.hypot(d[:, 0], d[:, 1]), d[:, 2]) / ref_norm
        out[f"theta_par={key[0][1]:g},{key[1]}"] = (float(rel.max()), int(m.sum()),
                                                    int((~m).sum()))
    return kappa, seed, out


def main() -> int:
    if not os.path.exists(EXE):
        print(f"{EXE} not found; run `make exp1_draw.exe`", file=sys.stderr)
        return 1
    jobs = [(k, s) for k in FRAME_KAPPAS for s in seeds_for(KAPPAS.index(k))]
    results = {}
    with ProcessPoolExecutor(max_workers=max(1, (os.cpu_count() or 2) // 2)) as pool:
        for kappa, seed, out in pool.map(one_seed, jobs):
            for combo, (mx, n, bad) in out.items():
                r = results.setdefault((kappa, combo), {"max": 0.0, "n": 0, "nonfinite": 0})
                r["max"] = max(r["max"], mx)
                r["n"] += n
                r["nonfinite"] += bad
    summary = [{"kappa": k, "combination": c, "max_rel_vector_diff": r["max"],
                "n_compared": r["n"], "n_nonfinite": r["nonfinite"]}
               for (k, c), r in sorted(results.items())]
    worst = max(s["max_rel_vector_diff"] for s in summary)
    n_all = sum(s["n_compared"] for s in summary)
    out = {
        "experiment": "Experiment 1 large-sample frame invariance -- kappa = 0.55, 2, 10",
        "n_per_seed": N_PER_SEED, "reference": "theta = (1, 2), field along z",
        "summary": summary, "max_rel_vector_diff": worst, "n_compared": n_all,
        "environment": {
            "platform": platform.platform(), "python": platform.python_version(),
            "numpy": np.__version__,
            "sampler_header_sha256": subprocess.run(
                ["shasum", "-a", "256", "../../cpp/bi_kappa_distribution.H"],
                capture_output=True, text=True).stdout.split()[0],
            "git_commit": subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True,
                                         text=True).stdout.strip(),
        },
    }
    with open("results/exp1_frame.json", "w") as fh:
        json.dump(out, fh, indent=1)
    lines = ["# Experiment 1 large-sample frame invariance -- kappa = 0.55, 2, 10", "",
             f"100 cell-test seeds x {N_PER_SEED} draws per combination; each run compared "
             "draw by draw with the (1, 2), B || z run at the same seed.", "",
             "| kappa | combination | max rel. diff | draws compared | non-finite |",
             "|---|---|---|---|---|"]
    for s in summary:
        lines.append(f"| {s['kappa']:g} | {s['combination']} | {s['max_rel_vector_diff']:.2e} "
                     f"| {s['n_compared']} | {s['n_nonfinite']} |")
    lines += ["", f"Overall: {n_all} draws compared, largest relative difference {worst:.2e}.", ""]
    with open("results/exp1_frame.md", "w") as fh:
        fh.write("\n".join(lines))
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
