#!/usr/bin/env python3
"""Regenerate every figure and numerical table in the JTJ1001 manuscript.

Everything the manuscript prints is produced here, from the frozen experiment
outputs in ``experiments/exp{1,2,3,4}_*``.  Nothing is transcribed by hand: the
LaTeX tables are emitted as ``\\input``-able fragments so that no number can
drift between the evidence and the paper.

Sources, and the reviewer comment each answers:

  exp1  raw/*.bin + results/exp1_results.json  -> marginal figure, validation and moment tables   R1.3, R1.5, R1.6
  exp2  results/exp2_results.json              -> Fig. 5, Table V            R1.1, R1.2
  exp3  results/exp3_results.json              -> Table VI                   R1.4
  exp4  results/exp4_results.json              -> Table VII                  R1.3

Usage
-----
    uv run --project ../../python python make_manuscript_assets.py

Writes into ``paper/overleaf/figures/`` and ``paper/overleaf/tables/``.  Those
live inside the git-ignored Overleaf project, so this script is the tracked
record of how they were made.
"""

from __future__ import annotations

import csv
import json
import os
import sys

import matplotlib as mpl

mpl.use("Agg")

import numpy as np
from matplotlib import pyplot as plt
from scipy import stats
from scipy.special import gammaln

import verify_cap_geometry
from verify_cap_geometry import cap_for_tv_target, speed_cap_anisotropy_limit

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

# The single (kappa, theta_par/theta_perp) pair Appendix A uses to illustrate the
# physical-speed cap's wide-cap anisotropy bias.  kappa = 0.75 is the low-kappa
# case already carried through Secs. IV and VI, and 2 is the anisotropy of every
# capped experiment, so the example costs the reader no new parameters.
SPEED_CAP_EXAMPLE = (0.75, 2.0)

# The (kappa, TV target) pair Sec. IV B uses to state how wide a component-wise
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

def load_manifest(exp_dir):
    with open(os.path.join(exp_dir, "raw", "manifest.csv"), newline="") as fh:
        return list(csv.DictReader(fh))


def load_runs(exp_dir, rows):
    """Concatenate the velocity vectors of the given manifest rows."""
    chunks = []
    for row in rows:
        path = os.path.join(exp_dir, "raw", row["file"])
        chunks.append(np.fromfile(path, dtype=np.float64).reshape(-1, 3))
    return np.vstack(chunks)


def load_json(exp_dir, name):
    with open(os.path.join(exp_dir, "results", name)) as fh:
        return json.load(fh)


def select(rows, **kw):
    def match(row):
        for key, val in kw.items():
            if abs(float(row[key]) - float(val)) > 1e-12 if key not in ("block", "mode", "ub_label") \
                    else row[key] != val:
                return False
        return True
    return [r for r in rows if match(r)]


def block_a_rows(rows, kappa):
    """Block A: theta = (1,2), B || z, uncapped, five seeds."""
    out = []
    for r in rows:
        if r["block"] != "A" or r["mode"] != "uncapped":
            continue
        if abs(float(r["kappa"]) - kappa) > 1e-12:
            continue
        out.append(r)
    return sorted(out, key=lambda r: int(r["seed"]))


# --------------------------------------------------------------------------
# Figure: Cartesian marginals and the bi-Maxwellian limit  (R1.6)
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
    rows = load_manifest(exp1_dir)
    kappas = [1.0, 2.0, 10.0]
    thetas = [1.0, 1.0, 2.0]          # theta_perp, theta_perp, theta_par
    comps = [
        (r"$v_{\perp1}/\theta_\perp$", "o", "0.15"),
        (r"$v_{\perp2}/\theta_\perp$", "s", "0.45"),
        (r"$v_\parallel/\theta_\parallel$", "^", "0.70"),
    ]

    half = 6.0
    bins = np.linspace(-half, half, 49)
    width = bins[1] - bins[0]
    centers = 0.5 * (bins[1:] + bins[:-1])
    grid = np.linspace(-half, half, 800)

    fig, axes = plt.subplots(3, 1, figsize=(3.4, 4.4), sharex=True, sharey=True)
    meta = {}

    for i, kappa in enumerate(kappas):
        ax = axes[i]
        sel = block_a_rows(rows, kappa)
        v = load_runs(exp1_dir, sel)
        n_total = v.shape[0]
        meta[kappa] = {"n_total": n_total, "seeds": [int(r["seed"]) for r in sel]}

        ax.plot(grid, bikappa_marginal_pdf(grid, kappa), color=C_KAPPA,
                label="bi-Kappa marginal", zorder=2)
        ax.plot(grid, maxwellian_marginal_pdf(grid), color=C_MAXW,
                ls="--", label="bi-Maxwellian limit", zorder=1)
        for j, (lab, mk, col) in enumerate(comps):
            s = v[:, j] / thetas[j]
            # Normalised by the FULL sample size, so the markers are the
            # probability density itself, not a density renormalised over
            # the displayed window.  Empty bins are not drawn.
            counts, _ = np.histogram(s, bins=bins)
            dens = counts / (n_total * width)
            keep = counts > 0
            ax.plot(centers[keep], dens[keep], mk, ms=2.4, mfc="none", mew=0.6,
                    color=col, label=lab, zorder=3)

        ax.set_yscale("log")
        ax.set_xlim(-half, half)
        ax.set_ylim(1e-6, 1.0)
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
# Figure: what the component-wise cap does  (R1.1, R1.2)
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
    """The Experiment 1 cell-count test (R1.3): chi^2 p-values over equal-probability cells.

    One row per kappa, the five replicate runs pooled.  "Radius" sums the cells over
    direction and tests the radial law alone; "all cells" tests the radius, the direction
    and their independence together.
    """
    res = load_json(exp1_dir, "exp1_results.json")
    rows = [r for r in res["summary"] if r["block"] == "A"]
    lines = []
    for rec in sorted(rows, key=lambda r: float(r["kappa"])):
        k = float(rec["kappa"])
        lines.append(
            f"${k:g}$ & {fmt(rec['cells_pooled_radius_pvalue'])} & "
            f"{fmt(rec['cells_pooled_all_pvalue'])} \\\\"
        )
    body = "\n".join(lines)
    path = os.path.join(OUT_TAB, "validation-summary.tex")
    with open(path, "w") as fh:
        fh.write(body + "\n")
    print(f"  wrote {path}")


def table_moments(exp1_dir):
    """Table IV -- second moments where they exist, with replicate spread (R1.5).

    Reports each replicate's sample variance, so statistical scatter is visible
    rather than inferred.  The caveat that matters: the variance of the sample
    variance needs a finite fourth moment, i.e. kappa > 5/2, so at kappa = 2 the
    spread below has no finite population value and must not be read as a
    standard error.
    """
    rows = load_manifest(exp1_dir)
    out_lines = []
    audit = {}
    for kappa in (2.0, 5.0, 10.0):
        sel = block_a_rows(rows, kappa)
        theta_perp, theta_par = 1.0, 2.0
        per_seed = {"vx": [], "vy": [], "vz": []}
        for r in sel:
            v = load_runs(exp1_dir, [r])
            per_seed["vx"].append(np.var(v[:, 0], ddof=1))
            per_seed["vy"].append(np.var(v[:, 1], ddof=1))
            per_seed["vz"].append(np.var(v[:, 2], ddof=1))

        second = kappa / (2.0 * kappa - 3.0)
        theory = {
            "vx": theta_perp**2 * second,
            "vy": theta_perp**2 * second,
            "vz": theta_par**2 * second,
        }
        names = {"vx": "$V_{\\perp1}$", "vy": "$V_{\\perp2}$", "vz": "$V_\\parallel$"}
        audit[kappa] = {}
        for key in ("vx", "vy", "vz"):
            vals = np.array(per_seed[key])
            mean, sd = vals.mean(), vals.std(ddof=1)
            rel = 100.0 * (mean - theory[key]) / theory[key]
            audit[kappa][key] = {"theory": theory[key], "mean": mean,
                                 "sd": sd, "rel_pct": rel}
            out_lines.append(
                f"${kappa:g}$ & {names[key]} & {fmt(theory[key], 4)} & "
                f"{fmt(mean, 4)} $\\pm$ {fmt(sd, 4)} & {rel:+.2f} \\\\"
            )
        out_lines.append("\\hline")
    body = "\n".join(out_lines[:-1])
    path = os.path.join(OUT_TAB, "moments.tex")
    with open(path, "w") as fh:
        fh.write(body + "\n")
    print(f"  wrote {path}")
    return audit


def table_cap(exp2_dir):
    """Table V -- the cap's cost and distortion, quoted as one number (R1.1, R1.2)."""
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


def macro_speed_cap_limit():
    """Appendix A and Sec. IV B velocity-bound numbers, as LaTeX macros.

    Appendix A states the wide-cap limit as an equation and illustrates it with a
    single value; it deliberately does not tabulate the (kappa, anisotropy) grid,
    because a sweep would present the bounding geometry as a study in its own
    right rather than as the analytic caution it is.  The number is emitted here
    anyway so that it keeps the same provenance rule as every other figure in the
    paper: nothing numeric is typed into the manuscript by hand.

    ``verify_cap_geometry.py`` checks this same function against an independent
    finite-cap quadrature carried out in velocity space; that check runs from
    ``main()`` below, so a drift in either one fails asset generation.

    The same fragment carries the component-wise cap's counterpart: the half-width
    at which the capped conditional law comes within a negligible total-variation
    distance of the intended one.  It is solved for, not tabulated, because the
    answer falls between the entries of any practical lambda ladder, and it is the
    number that makes the cost of the heavy-tailed cases concrete -- at
    kappa = 3/4 the required width is five orders of magnitude beyond the widest
    cap anyone would set.  ``verify_cap_geometry.py``'s check 8 pins it from both
    directions, so it has the same provenance guarantee as every table entry.
    """
    kappa, ratio = SPEED_CAP_EXAMPLE
    value = speed_cap_anisotropy_limit(kappa, ratio)
    path = os.path.join(OUT_TAB, "capgeom.tex")
    with open(path, "w") as fh:
        fh.write("% generated by paper/figures/make_manuscript_assets.py "
                 "-- do not edit\n")
        fh.write(f"\\newcommand{{\\SpeedCapLimit}}{{{value:.3f}}}\n")
        fh.write(f"\\newcommand{{\\SpeedCapLimitKappa}}{{{kappa:g}}}\n")
        fh.write(f"\\newcommand{{\\SpeedCapLimitRatio}}{{{ratio:g}}}\n")

        tv_kappa, tv_target = TV_TARGET_EXAMPLE
        lam = cap_for_tv_target(tv_kappa, tv_target)
        fh.write(f"\\newcommand{{\\TVThreshLambda}}{{{_math_sci(lam)}}}\n")
        fh.write(f"\\newcommand{{\\TVThreshKappa}}{{{tv_kappa:g}}}\n")
        fh.write(f"\\newcommand{{\\TVThreshTarget}}{{{_math_sci(tv_target)}}}\n")

    print(f"  wrote {path}")
    return value, lam


def table_performance(exp3_dir):
    """Table VI -- measured per-sample cost (R1.4)."""
    res = load_json(exp3_dir, "exp3_results.json")
    by = {}
    for rec in res["timing"]:
        if rec.get("variant") not in (None, "iso"):
            continue
        by.setdefault(float(rec["kappa"]), {})[rec["method"]] = rec
    order = ["gamma_ratio_spherical", "scale_mixture_normals", "pareto_rejection"]
    lines = []
    for k in sorted(by):
        cells = [f"${k:g}$"]
        for m in order:
            rec = by[k].get(m)
            if rec is None or not rec.get("applicable", True) \
                    or rec.get("ns_per_sample_median") is None:
                cells.append("n/a")
            else:
                cells.append(f"${rec['ns_per_sample_median']:.1f}$")
        base = by[k].get("gamma_ratio_spherical")
        par = by[k].get("pareto_rejection")
        if base and par and par.get("applicable", True) \
                and par.get("ns_per_sample_median"):
            ratio = par["ns_per_sample_median"] / base["ns_per_sample_median"]
            cells.append(f"${ratio:.2f}$")
        else:
            cells.append("n/a")
        lines.append(" & ".join(cells) + " \\\\")
    body = "\n".join(lines)
    path = os.path.join(OUT_TAB, "performance.tex")
    with open(path, "w") as fh:
        fh.write(body + "\n")
    print(f"  wrote {path}")


def table_precision(exp4_dir):
    """Table VII -- the measured operating envelope (R1.3, scoped per rule E)."""
    res = load_json(exp4_dir, "exp4_results.json")
    by = {}
    for r in res["released"]:
        by.setdefault(float(r["kappa"]), {}) \
          .setdefault(r["precision"], {})[r["stdlib"]] = r
    lines = []
    for k in sorted(by):
        cells = [f"${k:g}$"]
        for prec in ("double", "float"):
            per_lib = by[k].get(prec, {})
            if not per_lib:
                cells.append("---")
                continue
            # libc++ and libstdc++ agree to within seed noise, so the worst of
            # the two is quoted rather than an average that could hide one.
            frac = max(r["nonfinite_frac"] for r in per_lib.values())
            cells.append(sci(frac))
        lines.append(" & ".join(cells) + " \\\\")
    body = "\n".join(lines)
    path = os.path.join(OUT_TAB, "precision.tex")
    with open(path, "w") as fh:
        fh.write(body + "\n")
    print(f"  wrote {path}")


def main() -> int:
    os.makedirs(OUT_FIG, exist_ok=True)
    os.makedirs(OUT_TAB, exist_ok=True)

    # The Appendix A and Sec. IV B velocity-bound claims are analytic, so they are checked
    # before anything is written.  A drift between the closed forms in the
    # manuscript and independent quadrature must stop asset generation rather
    # than quietly emit a wrong number.
    print("velocity-bound closed forms (Appendix A and Sec. IV B):")
    if verify_cap_geometry.main() != 0:
        print("cap-geometry verification failed; no assets written",
              file=sys.stderr)
        return 1
    print()

    exp1 = os.path.join(EXP, "exp1_radial_directional")
    exp2 = os.path.join(EXP, "exp2_cap_characterization")
    exp3 = os.path.join(EXP, "exp3_benchmark")
    exp4 = os.path.join(EXP, "exp4_precision")

    print("figures:")
    m1 = figure_marginals(exp1)
    figure_cap(exp2)

    print("tables:")
    table_validation(exp1)
    audit = table_moments(exp1)
    table_cap(exp2)
    geom, tv_lam = macro_speed_cap_limit()
    table_performance(exp3)
    table_precision(exp4)

    print("\nprovenance for the captions:")
    for k, v in m1.items():
        print(f"  marginals kappa={k:g}: N={v['n_total']}, seeds={v['seeds']}")
    print(f"\nspeed-cap anisotropy limit at kappa={SPEED_CAP_EXAMPLE[0]:g}, "
          f"theta_par/theta_perp={SPEED_CAP_EXAMPLE[1]:g}: {geom:.6f}")
    print(f"cap width reaching TV = {TV_TARGET_EXAMPLE[1]:g} at "
          f"kappa={TV_TARGET_EXAMPLE[0]:g}: lambda = {tv_lam:.6e}")
    print("\nmoment audit (theory, mean, sd, rel%):")
    for k in sorted(audit):
        for key, d in audit[k].items():
            print(f"  kappa={k:g} {key}: {d['theory']:.4f} {d['mean']:.4f} "
                  f"{d['sd']:.4f} {d['rel_pct']:+.2f}%")
    return 0


if __name__ == "__main__":
    sys.exit(main())
