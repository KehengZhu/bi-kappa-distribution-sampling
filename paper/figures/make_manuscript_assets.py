#!/usr/bin/env python3
"""Regenerate every figure and numerical table in the manuscript.

Everything the manuscript prints is produced here, from the committed experiment
outputs in ``experiments/``.  Nothing is transcribed by hand: the LaTeX tables are
emitted as ``\\input``-able fragments so that no number can drift between the
evidence and the paper.

Sources:

  exp4  figures/fp1_failure_envelope.pdf        -> Fig. 2 (copied)
  exp2  results/exp2_results.json               -> Fig. 3, Table II
  exp1  results/exp1_tail.json                  -> Fig. 4, Table III (radius column)
  exp1  results/exp1_cells.json                 -> Table III (400-cell column)
  exp1  results/exp1_marginals.json             -> Fig. 5
  exp1  results/exp1_moments.json               -> Table IV
  verify_cap_geometry.py (closed form)          -> macros in tables/capgeom.tex

Usage
-----
    uv run --project ../../python python make_manuscript_assets.py

Writes into ``paper/overleaf/figures/`` and ``paper/overleaf/tables/``.
"""

from __future__ import annotations

import json
import os
import shutil
import sys

import matplotlib as mpl

mpl.use("Agg")

import numpy as np
from matplotlib import pyplot as plt
from scipy import stats
from scipy.special import gammaln

import verify_cap_geometry
from verify_cap_geometry import cap_for_tv_target

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
EXP = os.path.join(ROOT, "experiments")
OUT_FIG = os.path.join(ROOT, "paper", "overleaf", "figures")
OUT_TAB = os.path.join(ROOT, "paper", "overleaf", "tables")

# Manuscript figures are two-column REVTeX widetext panels.
mpl.rcParams.update({
    "font.size": 8,
    "axes.labelsize": 8,
    "axes.titlesize": 8,
    "legend.fontsize": 7,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "axes.linewidth": 0.6,
    "lines.linewidth": 1.1,
    "figure.dpi": 200,
    "savefig.bbox": "tight",
    "pdf.fonttype": 42,
})

C_HIST = "0.72"
C_KAPPA = "#1f4e9c"
C_MAXW = "#c1272d"

# The (kappa, TV target) pair Sec. V B uses to state how wide a component-wise
# cap would have to be before the capped conditional law is within a negligible
# total-variation distance of the intended one.  kappa = 0.75 is the same
# low-kappa case; 10^-3 is the negligibility threshold Experiment 2 fixed before
# any result was looked at.  ``verify_cap_geometry.py`` solves for the width and
# checks it two independent ways, so the number below cannot drift from the
# evidence.
TV_TARGET_EXAMPLE = verify_cap_geometry.TV_TARGET_EXAMPLE


# --------------------------------------------------------------------------
# Analytic reference laws
# --------------------------------------------------------------------------

def bikappa_marginal_pdf(s, kappa):
    """Exact marginal of the bi-Kappa law in the normalised variable s = v_i/theta_i.

    Integrating the 3-D density over the two remaining components gives
    ``p(v) = Gamma(k)/(sqrt(pi) Gamma(k-1/2) sqrt(k) theta) (1+v^2/(k theta^2))^-k``,
    so in ``s`` the thermal speed drops out entirely.  Every Cartesian component
    therefore has the *same* normalised marginal -- which is what makes the
    anisotropic scaling itself testable in these panels.
    """
    logn = gammaln(kappa) - gammaln(kappa - 0.5) - 0.5 * np.log(np.pi)
    return np.exp(logn) / np.sqrt(kappa) * (1.0 + s**2 / kappa) ** (-kappa)


def maxwellian_marginal_pdf(s):
    """Normalised bi-Maxwellian marginal, the kappa -> infinity limit."""
    return np.exp(-(s**2)) / np.sqrt(np.pi)


# --------------------------------------------------------------------------
# Shared loading
# --------------------------------------------------------------------------

def load_json(exp_dir, name):
    with open(os.path.join(exp_dir, "results", name)) as fh:
        return json.load(fh)


# --------------------------------------------------------------------------
# Figure: Cartesian marginals and the bi-Maxwellian limit
# --------------------------------------------------------------------------

def figure_marginals(exp1_dir):
    """One single-column figure: the three normalised components against the
    exact marginal and the bi-Maxwellian limit, on a logarithmic density axis.

    In s = v_i/theta_i all three components share one marginal, so overlaying
    them tests the anisotropic rescaling; the logarithmic axis shows the
    suprathermal tail, where the bi-Kappa and bi-Maxwellian laws differ, rather
    than only the core.  kappa = 1 has no finite variance; kappa = 2 and 10
    carry the approach to the bi-Maxwellian limit.
    """
    res = load_json(exp1_dir, "exp1_marginals.json")
    by_kappa = {r["kappa"]: r for r in res["summary"]}
    keys = ["v_perp1", "v_perp2", "v_par"]
    kappas = [1.0, 2.0, 10.0]
    thetas = [1.0, 1.0, 2.0]          # theta_perp, theta_perp, theta_par
    comps = [
        (r"$v_{\perp1}/\theta_\perp$", "o", "0.15"),
        (r"$v_{\perp2}/\theta_\perp$", "s", "0.45"),
        (r"$v_\parallel/\theta_\parallel$", "^", "0.70"),
    ]

    half = 6.0
    bins = np.array(res["bin_edges"])
    width = bins[1] - bins[0]
    centers = 0.5 * (bins[1:] + bins[:-1])
    grid = np.linspace(-half, half, 800)

    fig, axes = plt.subplots(3, 1, figsize=(3.4, 4.4), sharex=True, sharey=True)
    meta = {}

    for i, kappa in enumerate(kappas):
        ax = axes[i]
        rec = by_kappa[kappa]
        n_total = rec["n_total"]
        meta[kappa] = {"n_total": n_total, "seeds": rec["seeds"]}

        ax.plot(grid, bikappa_marginal_pdf(grid, kappa), color=C_KAPPA,
                label="bi-Kappa marginal", zorder=2)
        ax.plot(grid, maxwellian_marginal_pdf(grid), color=C_MAXW,
                ls="--", label="bi-Maxwellian limit", zorder=1)
        for j, (lab, mk, col) in enumerate(comps):
            # Counts of v_j / theta_j per bin, normalised by the FULL sample size,
            # so the markers are the probability density itself, not a density
            # renormalised over the displayed window.  Empty bins are not drawn.
            counts = np.array(rec["counts"][keys[j]])
            dens = counts / (n_total * width)
            keep = counts > 0
            ax.plot(centers[keep], dens[keep], mk, ms=2.4, mfc="none", mew=0.6,
                    color=col, label=lab, zorder=3)

        ax.set_yscale("log")
        ax.set_xlim(-half, half)
        ax.set_ylim(1e-8, 1.0)
        ax.set_ylabel("Probability density")
        ax.text(0.03, 0.93, rf"$\kappa={kappa:g}$", transform=ax.transAxes,
                va="top", ha="left")
    axes[-1].set_xlabel(r"$v_i/\theta_i$ (dimensionless)")

    handles, lab = axes[0].get_legend_handles_labels()
    fig.legend(handles, lab, frameon=False, ncol=2, fontsize=6.5,
               loc="upper center", bbox_to_anchor=(0.55, 1.06))
    fig.tight_layout(pad=0.4, rect=(0, 0, 1, 0.95))
    path = os.path.join(OUT_FIG, "validation-marginals.pdf")
    fig.savefig(path)
    plt.close(fig)
    print(f"  wrote {path}")
    return meta


# --------------------------------------------------------------------------
# Figure: what the component-wise cap does
# --------------------------------------------------------------------------

def capped_summary(exp2_dir):
    """Block-A capped runs, keyed by (kappa, lambda).

    Everything the cap figure and table need lives in ``summary``: the analytic
    rejected fraction (which *is* the TV distance) and the tail-quantile ratios.
    """
    res = load_json(exp2_dir, "exp2_results.json")
    table = {}
    for rec in res["summary"]:
        if rec.get("block") != "A" or rec.get("mode") != "capped":
            continue
        table[(float(rec["kappa"]), float(rec["lambda"]))] = rec
    kappas = sorted({k for k, _ in table})
    lams = sorted({l for _, l in table})
    probes = res["summary"][0]["q_probes"]
    return res, table, kappas, lams, probes.index(0.999)


def figure_cap(exp2_dir):
    _, table, kappas, lams, i999 = capped_summary(exp2_dir)
    qtab = table

    cmap = plt.get_cmap("viridis")
    colors = {k: cmap(t) for k, t in zip(kappas, np.linspace(0.05, 0.85, len(kappas)))}

    fig, axes = plt.subplots(1, 2, figsize=(7.1, 2.5))

    ax = axes[0]
    for k in kappas:
        xs = [l for l in lams if (k, l) in table]
        ys = [table[(k, l)]["reject_fraction_analytic"] for l in xs]
        ax.loglog(xs, ys, "o-", ms=3, color=colors[k], label=rf"$\kappa={k:g}$")
    ax.axhline(1e-3, color="0.4", ls=":", lw=0.8)
    ax.text(3.2, 1.3e-3, r"$\mathrm{TV}=10^{-3}$", fontsize=6, color="0.3")
    ax.set_xlabel(r"cap $\lambda$")
    ax.set_ylabel(r"rejected fraction $=\mathrm{TV}(f_\lambda,f_\kappa)$")
    ax.set_ylim(1e-8, 1.5)
    ax.legend(frameon=False, ncol=2, fontsize=6, loc="lower left")

    ax = axes[1]
    for k in kappas:
        xs = [l for l in lams if (k, l) in qtab]
        ys = [qtab[(k, l)]["q_speed_ratio"]["mean"][i999] for l in xs]
        ax.semilogx(xs, ys, "o-", ms=3, color=colors[k])
    ax.axhline(1.0, color="0.35", lw=0.7)
    ax.axhline(0.99, color="0.4", ls=":", lw=0.8)
    ax.set_xlabel(r"cap $\lambda$")
    ax.set_ylabel(r"$q_{99.9}(|v|)$ ratio, capped / uncapped")
    ax.set_ylim(-0.03, 1.08)

    fig.tight_layout(pad=0.4)
    path = os.path.join(OUT_FIG, "cap-characterization.pdf")
    fig.savefig(path)
    plt.close(fig)
    print(f"  wrote {path}")


# --------------------------------------------------------------------------
# Figure: radius-shell counts into the tail
# --------------------------------------------------------------------------

def figure_radius_shells(exp1_dir):
    """One single-column figure: the standardized deviation of each radius-shell count.

    The horizontal axis is the shell, placed by category: each of the ten shells of
    probability 1/10 gets one slot, and each of the five pieces of the outermost one a
    wider slot, since the tail is what the subdivision adds.  The tick labels give the
    shell edges as percentiles of R under the target (r_9 is the 90th), so the axis names
    positions in R, in the same terms as the 99.9th-percentile speed of Sec. V B, and
    shows how far into the tail the test resolves.  Each kappa is offset within the slot.  Under the target every deviation is
    approximately standard normal.
    """
    res = load_json(exp1_dir, "exp1_tail.json")
    rows = sorted(res["summary"], key=lambda r: float(r["kappa"]))
    n_shells = len(rows[0]["shell_z"])
    n_bulk = 9                       # shells inside r_9; the rest subdivide R > r_9
    widths = np.array([0.7] * n_bulk + [2.2] * (n_shells - n_bulk))
    edges = np.concatenate([[0.0], np.cumsum(widths)])
    centers = 0.5 * (edges[:-1] + edges[1:])

    fig, ax = plt.subplots(figsize=(3.4, 2.3))
    ax.axhspan(-2.0, 2.0, color="0.92", lw=0, zorder=0)
    ax.axhline(0.0, color="0.55", lw=0.6, zorder=1)
    ax.axvline(edges[n_bulk], color="0.35", lw=0.6, ls=":", zorder=1)

    cmap = mpl.colormaps["viridis"]
    frac = np.linspace(-0.36, 0.36, len(rows))
    for j, rec in enumerate(rows):
        k = float(rec["kappa"])
        col = cmap(0.05 + 0.85 * j / (len(rows) - 1))
        # Filled markers for kappa <= 3/2, where the second moments do not exist.
        filled = k <= 1.5
        ax.plot(centers + frac[j] * widths, rec["shell_z"], "o", ms=2.6, mew=0.7,
                color=col, mfc=col if filled else "white", label=rf"${k:g}$", zorder=3)

    # Ticks at the shell edges, labelled by the percentile of R there.
    major = {0: "$0$", 5: "$50$", 9: "$90$", 10: "$99$", 11: "$99.9$",
             12: "$99.99$", 13: "$99.999$", 14: "$100$"}
    ax.set_xticks([edges[i] for i in major])
    ax.set_xticklabels(list(major.values()))
    ax.set_xticks(edges, minor=True)
    ax.tick_params(axis="x", which="minor", length=1.5)
    ax.tick_params(axis="x", which="major", labelsize=6.5)
    ax.set_xlim(edges[0], edges[-1])
    ax.set_ylim(-3.6, 3.6)
    ax.set_yticks([-3, -2, -1, 0, 1, 2, 3])
    ax.set_xlabel(r"Percentile of $R$ at the shell edge")
    ax.set_ylabel(r"$(O-E)/\sigma$")
    ax.text(edges[n_bulk] + 0.2, 3.35, r"$R>r_9$", ha="left", va="top",
            fontsize=7, color="0.3")

    # Two legend rows: the six kappa <= 3/2 (filled) above, the three kappa > 3/2
    # (open) below.  Matplotlib fills columns first, so the handles are interleaved.
    handles, labels = ax.get_legend_handles_labels()
    order = [0, 6, 1, 7, 2, 8, 3, 4, 5]
    ax.legend([handles[i] for i in order], [labels[i] for i in order],
              title=r"$\kappa$", frameon=False, ncol=6, fontsize=6.5,
              title_fontsize=7, loc="lower center", bbox_to_anchor=(0.5, 1.0),
              handletextpad=0.1, columnspacing=0.7, borderaxespad=0.2)
    fig.tight_layout(pad=0.3)
    path = os.path.join(OUT_FIG, "validation-shells.pdf")
    fig.savefig(path)
    plt.close(fig)
    print(f"  wrote {path}")


# --------------------------------------------------------------------------
# Tables
# --------------------------------------------------------------------------

def fmt(x, nd=3):
    return f"{x:.{nd}f}"


def _math_sci(x, nd=1):
    """Scientific notation as a bare math-mode fragment, e.g. ``9.4\\times 10^{5}``.

    Unlike ``sci`` this emits no ``$``, so the manuscript can drop it inside an
    existing math environment; and a unit mantissa is suppressed, so a threshold
    of 10^-3 reads as a power of ten rather than as ``1.0 x 10^-3``.
    """
    mant, exp = f"{x:.{nd}e}".split("e")
    if float(mant) == 1.0:
        return rf"10^{{{int(exp)}}}"
    return rf"{mant}\times 10^{{{int(exp)}}}"


def sci(x, nd=1):
    """LaTeX scientific notation, e.g. 6.3\\times 10^{-4}."""
    if x == 0:
        return "$0$"
    if 1e-3 <= abs(x) < 1e3:
        return f"${x:.{nd + 2}g}$"
    mant, exp = f"{x:.{nd}e}".split("e")
    return rf"${mant}\times 10^{{{int(exp)}}}$"


def table_validation(exp1_dir):
    """The Experiment 1 cell-count test: chi^2 p-values.

    One row per kappa, on the large-sample draws (100 runs of 5 x 10^5 pooled,
    5 x 10^7 per kappa).  "Radius" is the test of the radius alone over 14 shells,
    the ten deciles with the outermost divided at the 99th to 99.999th percentiles
    (exp1_tail.py, on the same draws); "all cells" tests the radius, the direction
    and their independence together over 400 cells (exp1_cells.py).
    """
    res = load_json(exp1_dir, "exp1_cells.json")
    rows = res["summary"]
    radius = {float(r["kappa"]): r["radius_pvalue"]
              for r in load_json(exp1_dir, "exp1_tail.json")["summary"]}
    lines = []
    prev_k = None
    for rec in sorted(rows, key=lambda r: float(r["kappa"])):
        k = float(rec["kappa"])
        # A rule separates kappa <= 3/2, where the second moments do not exist.
        if prev_k is not None and prev_k <= 1.5 < k:
            lines.append("\\hline")
        prev_k = k
        lines.append(
            f"${k:g}$ & {fmt(radius[k])} & {fmt(rec['all_cells_pvalue'])} \\\\"
        )
    body = "\n".join(lines)
    path = os.path.join(OUT_TAB, "validation-summary.tex")
    with open(path, "w") as fh:
        fh.write(body + "\n")
    print(f"  wrote {path}")


def table_moments(exp1_dir):
    """Table IV: second moments where they exist.

    From the large-sample run (exp1_moments.py: the 100 cell-test runs of 5 x 10^5 at
    kappa = 2, 5, 10).  "Sample" is the mean of the 100 run variances +/- its standard
    error (sd of the run variances / 10).  The caveat that matters: the variance of the
    sample variance needs a finite fourth moment, i.e. kappa > 5/2, so at kappa = 2 the
    standard error has no finite population value and is only a rough scale.
    """
    res = load_json(exp1_dir, "exp1_moments.json")
    names = {"v_perp1": "$V_{\\perp1}$", "v_perp2": "$V_{\\perp2}$",
             "v_par": "$V_\\parallel$"}
    out_lines = []
    check = {}
    for rec in sorted(res["summary"], key=lambda r: r["kappa"]):
        kappa = rec["kappa"]
        check[kappa] = {}
        for key in ("v_perp1", "v_perp2", "v_par"):
            c = rec["components"][key]
            check[kappa][key] = {"theory": c["expected"], "mean": c["mean_run_variance"],
                                 "se": c["standard_error"], "rel_pct": c["diff_percent"],
                                 "diff_over_se": c["diff_over_se"]}
            diff = f"{c['diff_percent']:+.2f}"
            if float(diff) == 0.0:
                diff = "0.00"   # no sign on a difference that rounds to zero
            out_lines.append(
                f"${kappa:g}$ & {names[key]} & {fmt(c['expected'], 4)} & "
                f"{fmt(c['mean_run_variance'], 4)} $\\pm$ {fmt(c['standard_error'], 4)} & "
                f"{diff} \\\\"
            )
        out_lines.append("\\hline")
    body = "\n".join(out_lines[:-1])
    path = os.path.join(OUT_TAB, "moments.tex")
    with open(path, "w") as fh:
        fh.write(body + "\n")
    print(f"  wrote {path}")
    return check


def table_cap(exp2_dir):
    """Table II: the removed probability mass and the 99.9th-percentile ratio."""
    _, table, kappas, _, i999 = capped_summary(exp2_dir)
    lams = [3.0, 5.0, 10.0, 20.0, 50.0, 100.0]
    lines = []
    for k in kappas:
        cells = [f"${k:g}$"]
        for lam in lams:
            rec = table.get((k, lam))
            if rec is None:
                cells.append("---")
                continue
            tv = rec["reject_fraction_analytic"]
            q = rec["q_speed_ratio"]["mean"][i999]
            qs = "$<\\!0.001$" if q < 5e-4 else f"${q:.3f}$"
            cells.append(f"{sci(tv)} / {qs}")
        lines.append(" & ".join(cells) + " \\\\")
    body = "\n".join(lines)
    path = os.path.join(OUT_TAB, "cap.tex")
    with open(path, "w") as fh:
        fh.write(body + "\n")
    print(f"  wrote {path}")


def macro_tv_threshold():
    """Sec. V B velocity-bound numbers, as LaTeX macros.

    The half-width at which the capped conditional law comes within a negligible total-variation
    distance of the intended one.  It is solved for, not tabulated, because the
    answer falls between the entries of any practical lambda ladder, and it is the
    number that makes the cost of the heavy-tailed cases concrete -- at
    kappa = 3/4 the required width is five orders of magnitude beyond the widest
    cap anyone would set.  ``verify_cap_geometry.py``'s check 8 pins it from both
    directions, so it has the same provenance guarantee as every table entry.
    """
    path = os.path.join(OUT_TAB, "capgeom.tex")
    with open(path, "w") as fh:
        fh.write("% generated by paper/figures/make_manuscript_assets.py "
                 "-- do not edit\n")
        tv_kappa, tv_target = TV_TARGET_EXAMPLE
        lam = cap_for_tv_target(tv_kappa, tv_target)
        fh.write(f"\\newcommand{{\\TVThreshLambda}}{{{_math_sci(lam)}}}\n")
        fh.write(f"\\newcommand{{\\TVThreshKappa}}{{{tv_kappa:g}}}\n")
        fh.write(f"\\newcommand{{\\TVThreshTarget}}{{{_math_sci(tv_target)}}}\n")

    print(f"  wrote {path}")
    return lam


def figure_failure_envelope(exp4_dir):
    """Fig. 2: the finite-precision failure envelope, drawn by the experiment itself."""
    src = os.path.join(exp4_dir, "figures", "fp1_failure_envelope.pdf")
    path = os.path.join(OUT_FIG, "fp1_failure_envelope.pdf")
    shutil.copyfile(src, path)
    print(f"  wrote {path}")


def main() -> int:
    os.makedirs(OUT_FIG, exist_ok=True)
    os.makedirs(OUT_TAB, exist_ok=True)

    # The Sec. V B velocity-bound claims are analytic, so they are checked
    # before anything is written.  A drift between the closed forms in the
    # manuscript and independent quadrature must stop asset generation rather
    # than quietly emit a wrong number.
    print("velocity-bound closed forms (Sec. V B):")
    if verify_cap_geometry.main() != 0:
        print("cap-geometry verification failed; no assets written",
              file=sys.stderr)
        return 1
    print()

    exp1 = os.path.join(EXP, "exp1_radial_directional")
    exp2 = os.path.join(EXP, "exp2_cap_characterization")
    exp4 = os.path.join(EXP, "exp4_finite_precision")

    print("figures:")
    figure_failure_envelope(exp4)
    m1 = figure_marginals(exp1)
    figure_radius_shells(exp1)
    figure_cap(exp2)

    print("tables:")
    table_validation(exp1)
    check = table_moments(exp1)
    table_cap(exp2)
    tv_lam = macro_tv_threshold()

    print("\nprovenance for the captions:")
    for k, v in m1.items():
        print(f"  marginals kappa={k:g}: N={v['n_total']}, seeds={v['seeds']}")
    print(f"cap width reaching TV = {TV_TARGET_EXAMPLE[1]:g} at "
          f"kappa={TV_TARGET_EXAMPLE[0]:g}: lambda = {tv_lam:.6e}")
    print("\nmoment check (theory, mean, se, rel%, diff/se):")
    for k in sorted(check):
        for key, d in check[k].items():
            print(f"  kappa={k:g} {key}: {d['theory']:.4f} {d['mean']:.4f} "
                  f"{d['se']:.4f} {d['rel_pct']:+.2f}% {d['diff_over_se']:+.2f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
