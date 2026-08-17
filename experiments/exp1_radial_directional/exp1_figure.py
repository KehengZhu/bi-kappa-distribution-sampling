#!/usr/bin/env python3
"""Experiment 1 validation figure.

Four panels, each answering a distinct part of R1.3:
  (a) the radial law itself, as a P-P plot of W = 1/(1+T) against Beta(kappa-1/2, 3/2);
  (b) directional uniformity of cos(Theta);
  (c) radial-direction independence, cos(Theta) split by radial quartile;
  (d) the KS statistic across the whole kappa sweep, including the low-kappa range.

Captions are written to be self-contained (R1.6): every panel states kappa, theta_perp,
theta_par, N, the sampling mode and the curve identities.

Usage:  uv run --project ../../python python exp1_figure.py [raw_dir] [results_dir]

Note: this is the working validation figure. The manuscript version should be regenerated
through the project's figure workflow once the final section/figure numbering is fixed.
"""

from __future__ import annotations

import csv
import os
import sys

import matplotlib as mpl
import numpy as np
from matplotlib import pyplot as plt
from scipy import stats

from exp1_analyze import field_basis, is_axis_aligned_z

mpl.rcParams.update({
    "font.size": 8,
    "axes.labelsize": 8,
    "axes.titlesize": 8,
    "legend.fontsize": 7,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "axes.linewidth": 0.6,
    "lines.linewidth": 1.0,
    "figure.dpi": 200,
})

PANEL_KAPPAS = [0.55, 1.0, 2.0, 10.0]
COLORS = ["#1b6ca8", "#2a9d5c", "#d1701c", "#8b3a62"]


def load_run(raw_dir, row):
    v = np.fromfile(os.path.join(raw_dir, row["file"]), dtype=np.float64).reshape(-1, 3)
    ub = np.array([float(row["ub_x"]), float(row["ub_y"]), float(row["ub_z"])])
    local = v if is_axis_aligned_z(ub) else v @ field_basis(ub)
    kappa = float(row["kappa"])
    scales = np.sqrt(kappa) * np.array(
        [float(row["theta_perp"]), float(row["theta_perp"]), float(row["theta_par"])]
    )
    u = local / scales
    R = np.hypot(np.hypot(u[:, 0], u[:, 1]), u[:, 2])
    ok = np.isfinite(R) & (R > 0)
    u, R = u[ok], R[ok]
    with np.errstate(over="ignore", divide="ignore", under="ignore"):
        small = R <= 1.0
        W = np.empty_like(R)
        W[small] = 1.0 / (1.0 + R[small] ** 2)
        inv2 = (1.0 / R[~small]) ** 2
        W[~small] = inv2 / (1.0 + inv2)
    return R, W, u[:, 2] / R


def main() -> int:
    raw_dir = sys.argv[1] if len(sys.argv) > 1 else "raw"
    results_dir = sys.argv[2] if len(sys.argv) > 2 else "results"

    with open(os.path.join(raw_dir, "manifest.csv")) as fh:
        runs = [r for r in csv.DictReader(fh) if r["block"] == "A"]
    by_kappa_seed = {(float(r["kappa"]), int(r["seed"])): r for r in runs}

    theta_perp, theta_par = 1.0, 2.0
    n_per_run = int(runs[0]["n_samples"])

    fig, axes = plt.subplots(2, 2, figsize=(7.0, 5.6))
    ax_a, ax_b, ax_c, ax_d = axes.ravel()

    # ---- (a) radial law: P-P residual of W against its exact Beta law ----
    # The residual, not the raw P-P curve: on a raw P-P plot a correct sampler and a
    # sampler wrong by 1% both look like the diagonal.
    probs = np.linspace(0.001, 0.999, 400)
    for color, kappa in zip(COLORS, PANEL_KAPPAS):
        _, W, _ = load_run(raw_dir, by_kappa_seed[(kappa, 1001)])
        theo = stats.beta.cdf(np.quantile(W, probs), kappa - 0.5, 1.5)
        ax_a.plot(probs, 1e3 * (theo - probs), color=color, label=rf"$\kappa={kappa:g}$")
    band = 1e3 * 1.36 / np.sqrt(n_per_run)  # 95% KS acceptance band
    ax_a.axhspan(-band, band, color="0.90", zorder=0)
    ax_a.axhline(0.0, color="0.35", ls="--", lw=0.8)
    ax_a.text(0.02, band * 1.12, "95% KS band", fontsize=6.5, color="0.35")
    ax_a.set_xlabel("empirical quantile level (dimensionless)")
    ax_a.set_ylabel(r"$10^{3}\,[\,F_{\mathrm{theory}}-F_{\mathrm{empirical}}\,]$"
                    "\n(dimensionless)")
    ax_a.set_title(r"(a) radial law residual: $W=1/(1+R^{2})$", loc="left")
    ax_a.legend(frameon=False, loc="lower right", ncol=2)
    ax_a.set_xlim(0, 1)
    ax_a.set_ylim(-9, 9)

    # ---- (b) directional uniformity ----
    edges = np.linspace(-1, 1, 41)
    for color, kappa in zip(COLORS, PANEL_KAPPAS):
        _, _, cos_theta = load_run(raw_dir, by_kappa_seed[(kappa, 1001)])
        ax_b.hist(cos_theta, bins=edges, density=True, histtype="step",
                  color=color, label=rf"$\kappa={kappa:g}$")
    ax_b.axhline(0.5, color="0.35", ls="--", lw=0.8, label=r"uniform on $S^2$")
    ax_b.set_xlabel(r"$\cos\Theta$ (dimensionless)")
    ax_b.set_ylabel("probability density (dimensionless)")
    ax_b.set_title(r"(b) direction: $\cos\Theta$", loc="left")
    ax_b.set_ylim(0.40, 0.60)
    ax_b.legend(frameon=False, ncol=2, loc="upper center")

    # ---- (c) radial-direction independence ----
    R, _, cos_theta = load_run(raw_dir, by_kappa_seed[(2.0, 1001)])
    q = np.quantile(R, [0.0, 0.25, 0.5, 0.75, 1.0])
    for i in range(4):
        lo, hi = q[i], q[i + 1]
        sel = (R >= lo) & (R <= hi) if i == 3 else (R >= lo) & (R < hi)
        ax_c.hist(cos_theta[sel], bins=edges, density=True, histtype="step",
                  color=COLORS[i], label=f"radial quartile {i + 1}")
    ax_c.axhline(0.5, color="0.35", ls="--", lw=0.8)
    ax_c.set_xlabel(r"$\cos\Theta$ (dimensionless)")
    ax_c.set_ylabel("probability density (dimensionless)")
    ax_c.set_title(r"(c) independence of $R$ and direction, $\kappa=2$", loc="left")
    ax_c.set_ylim(0.35, 0.65)
    ax_c.legend(frameon=False, ncol=2, loc="upper center")

    # ---- (d) KS statistic across the sweep ----
    kappas = sorted({float(r["kappa"]) for r in runs})
    means, mins, maxs = [], [], []
    for kappa in kappas:
        vals = []
        for seed in (1001, 1002, 1003, 1004, 1005):
            _, W, _ = load_run(raw_dir, by_kappa_seed[(kappa, seed)])
            d = stats.kstest(W, "beta", args=(kappa - 0.5, 1.5)).statistic
            vals.append(d * np.sqrt(W.size))
        means.append(np.mean(vals))
        mins.append(np.min(vals))
        maxs.append(np.max(vals))
    means, mins, maxs = map(np.array, (means, mins, maxs))
    ax_d.errorbar(kappas, means, yerr=[means - mins, maxs - means],
                  fmt="o-", color=COLORS[0], ms=3.5, capsize=2,
                  label="observed (mean, range over 5 seeds)")
    ax_d.axhline(0.87, color="0.35", ls="--", lw=0.8,
                 label=r"expected under the exact law ($\approx 0.87$)")
    ax_d.axvspan(0.5, 1.5, color="0.90", zorder=0)
    ax_d.text(0.62, 2.05, r"$1/2<\kappa\leq 3/2$" "\n" "no second moment",
              fontsize=6.5, color="0.35")
    ax_d.set_xscale("log")
    ax_d.set_xlabel(r"$\kappa$ (dimensionless)")
    ax_d.set_ylabel(r"$\sqrt{n}\,D_{\mathrm{KS}}$ (dimensionless)")
    ax_d.set_title(r"(d) radial KS statistic vs $\kappa$", loc="left")
    ax_d.set_ylim(0, 2.4)
    ax_d.legend(frameon=False, loc="lower right")

    for ax in axes.ravel():
        ax.tick_params(direction="in", top=True, right=True)

    fig.suptitle(
        "Uncapped bi-Kappa sampler, "
        rf"$\theta_\perp={theta_perp:g}$, $\theta_\parallel={theta_par:g}$, "
        rf"$N={n_per_run // 1000}\,000$ samples per seed, "
        r"$\mathbf{B}\parallel\hat{z}$",
        fontsize=8, y=0.985,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.965))

    os.makedirs(results_dir, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(results_dir, f"exp1_validation.{ext}"), bbox_inches="tight")
    print(f"wrote {results_dir}/exp1_validation.pdf and .png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
