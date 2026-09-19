"""Experiment 6 figures.

Run with:  uv run --project ../../python python make_figures.py

This script plots and nothing else.  It reads only the committed summary CSVs under
``results/``, never the bulk binaries under ``raw/``, so every plotted point is recoverable
from a human-readable row.  No statistic is recomputed here.

Two main figures and three supplementary figures, each exported as PDF (canonical), SVG
(editable text) and a 600-dpi PNG preview, with ``figures/figure_manifest.json`` recording
the source CSVs and their hashes.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import sys
from datetime import datetime, timezone

import matplotlib as mpl

mpl.use("Agg")

import numpy as np                      # noqa: E402
from matplotlib import pyplot as plt    # noqa: E402
from matplotlib.lines import Line2D     # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "results")
FIGURES = os.path.join(HERE, "figures")

# Repository manuscript style, plus the two settings this experiment adds: editable text in
# SVG as well as PDF, because these figures are meant to be adjusted at proof stage.
mpl.rcParams.update({
    "font.size": 8,
    "axes.labelsize": 8,
    "axes.titlesize": 8,
    "legend.fontsize": 7,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "axes.linewidth": 0.7,
    "lines.linewidth": 1.1,
    "figure.dpi": 200,
    "figure.facecolor": "white",
    "savefig.facecolor": "white",
    "pdf.fonttype": 42,
    "svg.fonttype": "none",
    "axes.spines.top": False,
    "axes.spines.right": False,
})

MM = 1.0 / 25.4

# One restrained neutral/blue/vermilion palette.  Line style and marker shape carry the
# same distinction as colour, so every panel survives grayscale.
STYLE = {
    "QF":     dict(color="0.62", ls=":",  marker="o", mfc="none", ms=3.0, lw=0.9),
    "SPLIT":  dict(color="#1f4e9c", ls="-", marker="o", mfc="#1f4e9c", ms=3.4, lw=1.1),
    "LOG":    dict(color="#d55e00", ls="-", marker="s", mfc="#d55e00", ms=3.2, lw=1.1),
    "floor":  dict(color="black", ls="--", marker="D", mfc="none", ms=3.2, lw=0.9),
}
PASS_C, FAIL_C, GREY_C = "#1f4e9c", "#d55e00", "0.6"
PRIMARY_TAG = "libcxx"      # the declared primary environment

MANIFEST: list[dict] = []


# ---------------------------------------------------------------------------
def load(name: str) -> list[dict]:
    path = os.path.join(RESULTS, name)
    if not os.path.exists(path):
        return []
    with open(path, encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def num(row: dict, key: str, default=float("nan")) -> float:
    v = row.get(key)
    if v in (None, "", "NA"):
        return default
    if v == "inf":
        return float("inf")
    if v == "-inf":
        return float("-inf")
    try:
        return float(v)
    except ValueError:
        return default


def truthy(row: dict, key: str) -> bool:
    return str(row.get(key, "")).lower() == "true"


def sha256_of(name: str) -> str | None:
    path = os.path.join(RESULTS, name)
    if not os.path.exists(path):
        return None
    with open(path, "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()


def export(fig, basename: str, panels: int, conclusion: str, sources: list[str],
           width_mm: float, height_mm: float, interval_definition: str) -> None:
    os.makedirs(FIGURES, exist_ok=True)
    for ext, kw in (("pdf", {}), ("svg", {}), ("png", {"dpi": 600})):
        fig.savefig(os.path.join(FIGURES, f"{basename}.{ext}"), bbox_inches="tight", **kw)
    plt.close(fig)
    MANIFEST.append({
        "figure": basename, "panels": panels, "core_conclusion": conclusion,
        "source_csvs": sources,
        "source_csv_sha256s": {s: sha256_of(s) for s in sources},
        "script_sha256": hashlib.sha256(
            open(os.path.abspath(__file__), "rb").read()).hexdigest(),
        "width_mm": width_mm, "height_mm": height_mm,
        "interval_definition": interval_definition,
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
    })
    print(f"  wrote figures/{basename}.{{pdf,svg,png}}")


def panel_label(ax, letter: str, text: str = "") -> None:
    ax.set_title(f"{letter} {text}".rstrip(), loc="left", fontweight="bold", fontsize=8)


def plot_rate_series(ax, xs, rates, los, his, bounds, style, label):
    """A failure-probability curve on log axes, with zero counts drawn as upper bounds.

    A zero count has no place on a logarithmic axis: plotted at zero it disappears, and
    plotted at its estimate it asserts a precision the sample does not have.  It is drawn
    instead as a downward triangle sitting at the one-sided 95% upper limit.
    """
    xs = np.asarray(xs, dtype=float)
    rates = np.asarray(rates, dtype=float)
    obs = ~np.asarray(bounds, dtype=bool)
    if obs.any():
        ax.errorbar(xs[obs], rates[obs],
                    yerr=[rates[obs] - np.asarray(los)[obs],
                          np.asarray(his)[obs] - rates[obs]],
                    capsize=1.2, elinewidth=0.6, label=label, **style)
    if (~obs).any():
        ax.plot(xs[~obs], np.asarray(his)[~obs], marker="v", ls="none",
                color=style["color"], ms=3.4, mfc="none",
                label=None if obs.any() else label)
    return obs


# ---------------------------------------------------------------------------
# FP1 - failure envelope and mechanism decomposition
# ---------------------------------------------------------------------------
def figure_fp1() -> None:
    rows = [r for r in load("failure_envelope.csv")
            if r["scope"] == "pooled" and r["layer"] == "paired"
            and r.get("tag", PRIMARY_TAG) in (PRIMARY_TAG, "")]
    if not rows:
        print("  FP1: no data")
        return
    # failure_envelope.csv does not carry a tag column when pooled across seeds of one
    # toolchain, so the primary environment is selected by its standard library.
    rows = [r for r in rows if r["stdlib"] in ("libc++",)] or rows

    fig, axes = plt.subplots(1, 2, figsize=(183 * MM, 78 * MM))
    for ax, precision, letter in zip(axes, ("double", "float"), ("a", "b")):
        sub = [r for r in rows if r["precision"] == precision]
        if not sub:
            continue
        for method in ("QF", "SPLIT", "LOG"):
            ms = sorted([r for r in sub if r["method"] == method],
                        key=lambda r: num(r, "kappa"))
            xs = [num(r, "kappa") - 0.5 for r in ms]
            plot_rate_series(ax, xs, [num(r, "failure_rate") for r in ms],
                             [num(r, "failure_ci_lo") for r in ms],
                             [num(r, "failure_ci_hi") for r in ms],
                             [truthy(r, "failure_is_upper_bound") for r in ms],
                             STYLE[method], method)
        floor = sorted([r for r in sub if r["method"] == "SPLIT"],
                       key=lambda r: num(r, "kappa"))
        fx = [num(r, "kappa") - 0.5 for r in floor]
        plot_rate_series(ax, fx, [num(r, "honest_floor_rate") for r in floor],
                         [num(r, "honest_floor_ci_lo") for r in floor],
                         [num(r, "honest_floor_ci_hi") for r in floor],
                         [truthy(r, "honest_floor_is_upper_bound") for r in floor],
                         STYLE["floor"], "honest floor")

        # Shade the avoidable gap only where both curves are resolved observations.
        sp = {num(r, "kappa"): r for r in sub if r["method"] == "SPLIT"}
        ks = sorted(k for k in sp
                    if not truthy(sp[k], "failure_is_upper_bound")
                    and not truthy(sp[k], "honest_floor_is_upper_bound")
                    and num(sp[k], "failure_rate") > num(sp[k], "honest_floor_rate"))
        if len(ks) >= 2:
            ax.fill_between([k - 0.5 for k in ks],
                            [num(sp[k], "honest_floor_rate") for k in ks],
                            [num(sp[k], "failure_rate") for k in ks],
                            color="#1f4e9c", alpha=0.10, lw=0, zorder=0)

        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel(r"$\kappa-1/2$ (dimensionless)")
        panel_label(ax, letter, "double precision" if precision == "double"
                    else "single precision")
        ax.set_ylim(1.5e-7, 2.5)
    axes[0].set_ylabel("Failure probability per attempt\n(dimensionless)")
    handles = [Line2D([], [], label=k, **{kk: vv for kk, vv in STYLE[k].items()})
               for k in ("QF", "SPLIT", "LOG", "floor")]
    for h, lab in zip(handles, ("QF  $\\sqrt{X_1/X_2}$", "SPLIT  $\\sqrt{X_1}/\\sqrt{X_2}$",
                                "LOG  log domain", "honest overflow floor")):
        h.set_label(lab)
    from matplotlib.patches import Patch
    handles.append(Patch(facecolor="#1f4e9c", alpha=0.10,
                         label="avoidable loss (shaded)"))
    axes[1].legend(handles=handles, frameon=False, loc="lower left", fontsize=6.5)
    fig.tight_layout(pad=0.4)
    export(fig, "fp1_failure_envelope", 2,
           "The log-domain path removes the avoidable finite-precision loss of the released "
           "Gamma-ratio loader and reaches the honest final-output representability floor.",
           ["failure_envelope.csv"], 183, 78,
           "two-sided 95% Clopper-Pearson; zero counts shown as one-sided 95% upper limits "
           "(downward triangles)")


# ---------------------------------------------------------------------------
# FP2 - state-dependent failure and tail consequence
# ---------------------------------------------------------------------------
FP2_CASES = (("double", 0.51), ("float", 0.55))


def figure_fp2() -> None:
    bins = load("conditioning_bins.csv")
    tails = load("tail_metrics.csv")
    if not bins or not tails:
        print("  FP2: no data")
        return

    fig, axes = plt.subplots(2, 2, figsize=(183 * MM, 115 * MM))
    ratio_axes = []
    for col, (precision, kappa) in enumerate(FP2_CASES):
        ax = axes[0, col]
        for method in ("SPLIT", "LOG"):
            sub = [r for r in bins
                   if r["bin_kind"] == "target_upper_tail_q" and r["method"] == method
                   and r["precision"] == precision
                   and abs(num(r, "kappa") - kappa) < 1e-9
                   and r["tag"] == PRIMARY_TAG and num(r, "n_in_bin") >= 20]
            if not sub:
                continue
            sub.sort(key=lambda r: -num(r, "bin_hi"))
            x = [math.sqrt(max(num(r, "bin_hi"), 1e-300) * max(num(r, "bin_lo"), 1e-6))
                 for r in sub]
            y = [num(r, "success_rate") for r in sub]
            lo = [max(num(r, "success_rate") - num(r, "ci_lo"), 0) for r in sub]
            hi = [max(num(r, "ci_hi") - num(r, "success_rate"), 0) for r in sub]
            ax.errorbar(x, y, yerr=[lo, hi], capsize=1.2, elinewidth=0.6,
                        label=method, **STYLE[method])
        ax.set_xscale("log")
        ax.invert_xaxis()
        ax.set_ylim(-0.05, 1.05)
        ax.set_xlabel(r"target upper-tail probability $q=\Pr(R>r)$")
        panel_label(ax, "ab"[col], rf"{precision}, $\kappa={kappa:g}$")
        if col == 0:
            ax.set_ylabel("Pr(finite output | tail bin)\n(dimensionless)")
            ax.legend(frameon=False, loc="lower left", fontsize=6.5)

        ax = axes[1, col]
        ratio_axes.append(ax)
        levels = [0.99, 0.999, 0.9999]
        xpos = np.arange(len(levels))
        for k, method in enumerate(("SPLIT", "LOG")):
            sub = {num(r, "p"): r for r in tails
                   if r["method"] == method and r["precision"] == precision
                   and abs(num(r, "kappa") - kappa) < 1e-9 and r["tag"] == PRIMARY_TAG}
            off = (k - 0.5) * 0.22
            st = STYLE[method]
            for kind, mfc, dx in (("rho_attempt", st["color"], -0.055),
                                  ("rho_cond", "none", 0.055)):
                ys, los, his, xs = [], [], [], []
                for i, p in enumerate(levels):
                    r = sub.get(p)
                    if not r:
                        continue
                    xs.append(xpos[i] + off + dx)
                    ys.append(num(r, kind))
                    los.append(max(num(r, kind) - num(r, f"{kind}_lo"), 0))
                    his.append(max(num(r, f"{kind}_hi") - num(r, kind), 0))
                if xs:
                    ax.errorbar(xs, ys, yerr=[los, his], ls="none", capsize=1.2,
                                elinewidth=0.6, color=st["color"], marker=st["marker"],
                                mfc=mfc, ms=st["ms"])
        ax.axhline(1.0, color="black", lw=0.7, ls="-", zorder=0)
        ax.set_xticks(xpos)
        ax.set_xticklabels([f"{p:g}" for p in levels])
        ax.set_xlabel("target percentile $p$")
        panel_label(ax, "cd"[col], rf"{precision}, $\kappa={kappa:g}$")
        if col == 0:
            ax.set_ylabel("tail-retention ratio\n(dimensionless)")
            handles = [
                Line2D([], [], color=STYLE["SPLIT"]["color"], marker="o",
                       mfc=STYLE["SPLIT"]["color"], ls="none", ms=3.4, label="SPLIT"),
                Line2D([], [], color=STYLE["LOG"]["color"], marker="s",
                       mfc=STYLE["LOG"]["color"], ls="none", ms=3.2, label="LOG"),
                Line2D([], [], color="0.3", marker="o", mfc="0.3", ls="none", ms=3.4,
                       label=r"$\rho^{\rm attempt}$"),
                Line2D([], [], color="0.3", marker="o", mfc="none", ls="none", ms=3.4,
                       label=r"$\rho^{\rm cond}$"),
            ]
            ax.legend(handles=handles, frameon=False, loc="best", fontsize=6.5, ncol=2)
    # Panels c and d are the same metric at two precisions and are meant to be read
    # against each other, so they get one set of limits.
    if len(ratio_axes) == 2:
        lo = min(a.get_ylim()[0] for a in ratio_axes)
        hi = max(a.get_ylim()[1] for a in ratio_axes)
        for a in ratio_axes:
            a.set_ylim(lo, hi)
    fig.tight_layout(pad=0.4)
    export(fig, "fp2_conditioning_tail", 4,
           "Released-path failure is concentrated in the intended tail, so conditioning on "
           "finite generation changes tail observables; the log path restores the "
           "representable tail without hiding honest overflow.",
           ["conditioning_bins.csv", "tail_metrics.csv"], 183, 115,
           "binomial 95% for success rates; 99% seed-stratified bootstrap (10000 resamples) "
           "for tail ratios")


# ---------------------------------------------------------------------------
# SFP1 - scalar LOG validation
# ---------------------------------------------------------------------------
def figure_sfp1(selected: str) -> None:
    ecdf = [r for r in load("scalar_ecdf.csv")
            if r["method"] == selected and r["tag"] == PRIMARY_TAG]
    scal = [r for r in load("scalar_validation.csv")
            if r["method"] == selected and r["tag"] == PRIMARY_TAG]
    if not ecdf:
        print("  SFP1: no data")
        return

    panel_k = [0.501, 0.505, 0.51, 2.0]
    cmap = plt.get_cmap("viridis")
    kcol = {k: cmap(t) for k, t in zip(panel_k, np.linspace(0.05, 0.85, len(panel_k)))}

    fig, axes = plt.subplots(2, 2, figsize=(183 * MM, 115 * MM))
    for ax, precision, letter in zip(axes[0], ("double", "float"), ("a", "b")):
        for k in panel_k:
            sub = sorted([r for r in ecdf if r["transform"] == "Z"
                          and r["precision"] == precision
                          and abs(num(r, "kappa") - k) < 1e-9],
                         key=lambda r: num(r, "null_cdf"))
            if not sub:
                continue
            x = [num(r, "null_cdf") for r in sub]
            y = [num(r, "ecdf_residual") for r in sub]
            ax.plot(x, y, color=kcol[k], lw=0.9, label=rf"$\kappa={k:g}$")
            ax.axhline(num(sub[0], "band_hi"), color=kcol[k], lw=0.4, ls=":")
            ax.axhline(num(sub[0], "band_lo"), color=kcol[k], lw=0.4, ls=":")
        ax.axhline(0, color="black", lw=0.6)
        ax.set_xlabel(r"$F_0(Z)$, unit-exponential CDF")
        panel_label(ax, letter, f"{precision}: $Z=-\\log I_W(a,3/2)$")
        if letter == "a":
            ax.set_ylabel("ECDF residual\n(dimensionless)")
            ax.legend(frameon=False, fontsize=6.5, ncol=2)

    ax = axes[1, 0]
    # The two routes agree to the pixel, so drawing them with the same weight would hide
    # the agreement rather than show it: the W curve would vanish under the Z curve.  Z is
    # drawn as a wide pale band and W as a thin dashed line riding on it, so the reader can
    # see that the dashes track the band.
    for k in (0.55, 2.0):
        for transform, ls, lw, alpha in (("Z", "-", 2.2, 0.35), ("W", (0, (4, 3)), 0.9, 1.0)):
            sub = sorted([r for r in ecdf if r["transform"] == transform
                          and r["precision"] == "double"
                          and abs(num(r, "kappa") - k) < 1e-9],
                         key=lambda r: num(r, "null_cdf"))
            if not sub:
                continue
            xs = [num(r, "null_cdf_z_orientation") for r in sub]
            ys = [num(r, "ecdf_residual_z_orientation") for r in sub]
            order = np.argsort(xs)
            ax.plot(np.asarray(xs)[order], np.asarray(ys)[order],
                    color=kcol.get(k, "0.4"), ls=ls, lw=lw, alpha=alpha,
                    label=rf"$\kappa={k:g}$, {transform}")
    ax.axhline(0, color="black", lw=0.6)
    ax.set_xlabel(r"$F_0(Z)$, both routes on the $Z$ orientation")
    ax.set_ylabel("ECDF residual\n(dimensionless)")
    panel_label(ax, "c", "overlap: $W$ and $Z$ routes, double")
    ax.legend(frameon=False, fontsize=6.5, ncol=2)

    ax = axes[1, 1]
    levels = [0.5, 0.9, 0.99, 0.999, 0.9999]
    for j, precision in enumerate(("double", "float")):
        for i, k in enumerate(panel_k):
            rs = [r for r in scal if r["precision"] == precision
                  and abs(num(r, "kappa") - k) < 1e-9]
            if not rs:
                continue
            xs, ys, los, his = [], [], [], []
            for li, p in enumerate(levels):
                errs = [num(r, f"log_r_q{p}_error") for r in rs
                        if truthy(r, f"log_r_q{p}_resolved")]
                if not errs:
                    continue
                xs.append(li + (j - 0.5) * 0.3 + (i - 1.5) * 0.06)
                ys.append(float(np.mean(errs)))
                lo = [num(r, f"log_r_q{p}_lo") - num(r, f"log_r_q{p}_target") for r in rs
                      if truthy(r, f"log_r_q{p}_resolved")]
                hi = [num(r, f"log_r_q{p}_hi") - num(r, f"log_r_q{p}_target") for r in rs
                      if truthy(r, f"log_r_q{p}_resolved")]
                los.append(max(ys[-1] - float(np.mean(lo)), 0))
                his.append(max(float(np.mean(hi)) - ys[-1], 0))
            if xs:
                ax.errorbar(xs, ys, yerr=[los, his], ls="none", capsize=1.0,
                            elinewidth=0.5, ms=3.0, color=kcol[k],
                            marker="o" if precision == "double" else "^",
                            mfc=kcol[k] if precision == "double" else "none")
    ax.axhline(0, color="black", lw=0.7)
    # Symmetric log: the order-statistic interval at p = 0.9999 spans hundreds of log units
    # at the smallest shapes, where the radial scale itself is that large, and a linear axis
    # would compress every resolved error onto the zero line.
    ax.set_yscale("symlog", linthresh=1.0)
    ax.set_xticks(range(len(levels)))
    ax.set_xticklabels([f"{p:g}" for p in levels])
    ax.set_xlabel("target percentile $p$")
    ax.set_ylabel(r"$\widehat{\log r_p}-\log r_p$" "\n(natural log units, symlog)")
    panel_label(ax, "d", "log-radius quantile error")
    ax.legend(handles=[Line2D([], [], color="0.3", marker="o", mfc="0.3", ls="none",
                              ms=3, label="double"),
                       Line2D([], [], color="0.3", marker="^", mfc="none", ls="none",
                              ms=3, label="float")],
              frameon=False, fontsize=6.5, loc="best")
    fig.tight_layout(pad=0.4)
    export(fig, "sfp1_scalar_validation", 4,
           f"The selected primitive ({selected}) samples the intended radial law on a stable "
           "scale at the lower-limit stress cases and agrees with the ordinary W diagnostic "
           "where both are resolved.",
           ["scalar_ecdf.csv", "scalar_validation.csv"], 183, 115,
           "simultaneous Kolmogorov bands at familywise alpha 0.01; exact order-statistic "
           "quantile intervals, Bonferroni-corrected over five levels")


# ---------------------------------------------------------------------------
# SFP2 - complete-loader validation matrix
# ---------------------------------------------------------------------------
SFP2_COLS = ["radial", "direction", "independence", "anisotropy", "frame", "finiteness",
             "cap_law"]
SYMBOL = {"pass": ("o", PASS_C, PASS_C), "fail": ("X", FAIL_C, FAIL_C),
          "unresolved": ("^", GREY_C, "none"), "na": (".", GREY_C, GREY_C)}


def figure_sfp2() -> None:
    rows = [r for r in load("validation_matrix.csv") if r["tag"] == PRIMARY_TAG]
    if not rows:
        print("  SFP2: no data")
        return
    order = {c: i for i, c in enumerate(["C0", "C1", "C2", "C3", "C4", "C5", "C6",
                                         "NC1", "NC2"])}
    rows.sort(key=lambda r: (order.get(r["case"], 99), r["method"]))
    labels = [f"{r['case']} {r['method']}  (n={int(num(r, 'n')):,})" for r in rows]

    fig, ax = plt.subplots(figsize=(183 * MM, 92 * MM))
    for y, r in enumerate(rows):
        for x, col in enumerate(SFP2_COLS):
            marker, edge, face = SYMBOL.get(r.get(col, "na"), SYMBOL["na"])
            ms = 5.5 if r.get(col) != "na" else 2.0
            ax.plot(x, y, marker=marker, color=edge, mfc=face, ms=ms, mew=1.0, ls="none")
    ax.set_xticks(range(len(SFP2_COLS)))
    ax.set_xticklabels([c.replace("_", " ") for c in SFP2_COLS], rotation=30, ha="right")
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels(labels, fontsize=6.5)
    ax.set_xlim(-0.6, len(SFP2_COLS) - 0.4)
    ax.set_ylim(-0.7, len(rows) - 0.3)
    ax.invert_yaxis()
    ax.grid(True, color="0.92", lw=0.5)
    ax.set_axisbelow(True)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_visible(False)
    ax.tick_params(length=0)
    handles = [Line2D([], [], marker=m, color=e, mfc=f, ls="none", ms=5.5, label=lab)
               for lab, (m, e, f) in
               (("pass / detected as intended", SYMBOL["pass"]),
                ("fail", SYMBOL["fail"]),
                ("unresolved", SYMBOL["unresolved"]),
                ("not applicable", SYMBOL["na"]))]
    ax.legend(handles=handles, frameon=False, fontsize=6.5, ncol=4,
              loc="upper center", bbox_to_anchor=(0.5, 1.13))
    fig.tight_layout(pad=0.4)
    export(fig, "sfp2_loader_validation", 1,
           "The selected log primitive remains correct after direction sampling, anisotropic "
           "scaling, rotation, cap conditioning and final-component representability checks.",
           ["validation_matrix.csv"], 183, 92,
           "Holm-corrected decisions at familywise alpha 0.01; NC1 and NC2 pass only when "
           "the injected defect is detected")


# ---------------------------------------------------------------------------
# SFP3 - portability and performance
# ---------------------------------------------------------------------------
def figure_sfp3() -> None:
    port = [r for r in load("portability.csv")
            if r["scope"] == "pooled" and r["layer"] == "end_to_end"
            and str(r.get("completed", "")).lower() == "true"]
    perf = [r for r in load("performance.csv") if r.get("scope") == "ratio"]
    if not port and not perf:
        print("  SFP3: no data")
        return

    fig, axes = plt.subplots(1, 2, figsize=(183 * MM, 82 * MM))

    ax = axes[0]
    cases = [("C1", "o"), ("C3", "^")]
    platforms = sorted({(r["arch"], r["stdlib"]) for r in port})
    for pi, (arch, stdlib) in enumerate(platforms):
        for ci, (case, mk) in enumerate(cases):
            for mi, method in enumerate(("SPLIT", "LOG")):
                rs = [r for r in port if r["arch"] == arch and r["stdlib"] == stdlib
                      and r["case"] == case and r["method"] == method]
                if not rs:
                    continue
                r = rs[0]
                x = pi + (ci - 0.5) * 0.36 + (mi - 0.5) * 0.14
                st = STYLE[method]
                if truthy(r, "failure_is_upper_bound"):
                    ax.plot(x, num(r, "failure_ci_hi"), marker="v", ls="none",
                            color=st["color"], mfc="none", ms=4.0)
                else:
                    ax.errorbar([x], [num(r, "failure_rate")],
                                yerr=[[max(num(r, "failure_rate") - num(r, "failure_ci_lo"), 0)],
                                      [max(num(r, "failure_ci_hi") - num(r, "failure_rate"), 0)]],
                                ls="none", capsize=1.2, elinewidth=0.6, marker=mk,
                                color=st["color"], mfc=st["color"], ms=3.6)
    ax.set_yscale("log")
    ax.set_xticks(range(len(platforms)))
    ax.set_xticklabels([f"{a}\n{s}" for a, s in platforms], fontsize=6.5)
    ax.set_ylabel("Failure probability per attempt\n(dimensionless)")
    panel_label(ax, "a", "portability")
    ax.legend(handles=[
        Line2D([], [], color=STYLE["SPLIT"]["color"], marker="s", ls="none", ms=3.6,
               label="SPLIT"),
        Line2D([], [], color=STYLE["LOG"]["color"], marker="s", ls="none", ms=3.4,
               label="LOG"),
        Line2D([], [], color="0.3", marker="o", mfc="none", ls="none", ms=3.6,
               label=r"C1 (double, $\kappa=0.51$, rotated)"),
        Line2D([], [], color="0.3", marker="^", mfc="none", ls="none", ms=3.6,
               label=r"C3 (float, $\kappa=0.55$, rotated)"),
        Line2D([], [], color="0.3", marker="v", mfc="none", ls="none", ms=4.0,
               label="upper bound (no failure observed)")],
        frameon=False, fontsize=6.0, loc="best")

    ax = axes[1]
    tags = sorted({str(r.get("tag", "")) for r in perf})
    cases = sorted({str(r["case"]) for r in perf})
    for r in perf:
        i = cases.index(str(r["case"]))
        ti = tags.index(str(r.get("tag", "")))
        x = i + (ti - (len(tags) - 1) / 2) * 0.22
        med = num(r, "median_time_per_return_ratio")
        lo, hi = num(r, "ratio_ci_lo"), num(r, "ratio_ci_hi")
        ax.errorbar([x], [med], yerr=[[max(med - lo, 0)], [max(hi - med, 0)]],
                    ls="none", capsize=1.5, elinewidth=0.7, marker="s", ms=3.6,
                    color=STYLE["LOG"]["color"],
                    mfc=STYLE["LOG"]["color"] if ti == 0 else "none")
    labels = []
    for c in cases:
        r = next(r for r in perf if str(r["case"]) == c)
        cap = r.get("cap")
        labels.append(f"{c}\n{r['precision']}, $\\kappa$={num(r, 'kappa'):g}\n"
                      f"{'uncapped' if cap in (None, '', 'NA') else 'cap ' + str(cap)}")
    ax.axhline(1.0, color="black", lw=0.7)
    ax.set_xticks(range(len(cases)))
    ax.set_xticklabels(labels, fontsize=5.5)
    ax.legend(handles=[Line2D([], [], color=STYLE["LOG"]["color"], marker="s",
                              mfc=STYLE["LOG"]["color"] if i == 0 else "none",
                              ls="none", ms=3.6, label=t)
                       for i, t in enumerate(tags)],
              frameon=False, fontsize=6.0, loc="best")
    ax.set_ylabel("LOG / SPLIT time per returned sample\n(dimensionless)")
    panel_label(ax, "b", "cost")
    fig.tight_layout(pad=0.4)
    export(fig, "sfp3_portability_performance", 2,
           "The measured mitigation and its cost are reproducible across the tested "
           "environments, and any platform-specific limitation is visible rather than "
           "averaged away.",
           ["portability.csv", "performance.csv"], 183, 82,
           "two-sided 95% Clopper-Pearson for rates, one-sided 95% upper limits for zero "
           "counts; 95% bootstrap percentile interval for the paired blockwise time ratio")


# ---------------------------------------------------------------------------
def write_captions(selected: str) -> None:
    caps = {
        "fp1_failure_envelope":
            "Failure envelope and mechanism decomposition of the Gamma-ratio bi-Kappa "
            "radius. Probability that one attempt fails to return a finite three-vector, "
            "against $\\kappa-1/2$, in (a) double and (b) single precision, for the "
            "quotient-first formation $\\sqrt{X_1/X_2}$, the released split formation "
            "$\\sqrt{X_1}/\\sqrt{X_2}$, and the log-domain path. An attempt is one draw of "
            "the four high-level variates in the canonical isotropic, unrotated, uncapped "
            "configuration, $\\theta_\\perp=\\theta_\\parallel$, $\\hat b=\\hat z$; the "
            "curves are pooled over five frozen seeds of $10^6$ attempts each, with "
            "seed-level source rows retained. The honest overflow floor is the probability "
            "that at least one component of the intended returned vector is not "
            "representable in the run's floating-point type, judged against the exact "
            "threshold at which a real value rounds to infinity; no formation can go below "
            "it. Error bars are two-sided 95 per cent Clopper-Pearson intervals; a "
            "configuration with no observed failure is drawn as a downward triangle at its "
            "one-sided 95 per cent upper limit rather than at zero. Shading marks the "
            "avoidable loss, the gap between the released formation and the floor, where "
            "both are resolved observations.",
        "fp2_conditioning_tail":
            "State-dependent failure and its tail consequence. (a, b) Probability that an "
            "attempt returns a finite vector, as a function of the upper-tail probability "
            "$q=\\Pr(R>r)$ of the radius the attempt intended to produce, so that the "
            "horizontal axis runs from the body of the distribution on the left to the "
            "extreme tail on the right; bars are 95 per cent binomial intervals and bins "
            "holding fewer than twenty attempts are not drawn. (c, d) Tail retention at "
            "three target percentiles. Filled markers are the per-attempt ratio "
            "$\\rho^{\\rm attempt}(p)=\\Pr(R>r_p, S)/(1-p)$, in which a failed attempt "
            "stays outside the returned mass; open markers are the conditional ratio "
            "$\\rho^{\\rm cond}(p)=\\Pr(R>r_p \\mid S)/(1-p)$. For independent attempts the "
            "conditional law of the finite survivors is also the law produced by silently "
            "redrawing until success, so the open markers describe both. The failure atom "
            "is a mass, not a percentile error, and the two are measured separately. Bars "
            "are 99 per cent seed-stratified bootstrap intervals from 10\\,000 resamples; "
            "the black line marks the target itself.",
        "sfp1_scalar_validation":
            f"Scalar validation of the selected log primitive ({selected}). (a, b) Residuals "
            "of the empirical CDF of $Z=-\\log I_W(a,3/2)$ from the unit-exponential CDF, "
            "where $W=X_2/(X_1+X_2)\\sim\\mathrm{Beta}(a,3/2)$ and $a=\\kappa-1/2$; $Z$ is "
            "used because $W$ itself rounds to zero at the smallest shapes while $\\log W$ "
            "does not. Dotted lines are simultaneous Kolmogorov bands at the frozen "
            "familywise level. (c) Overlap check at the two cases where both routes are "
            "numerically resolved: the ordinary Beta diagnostic on $W$ and the $Z$ "
            "diagnostic must agree, and their test decisions do. (d) Error in the empirical "
            "$\\log R$ quantile relative to the exact target quantile, with exact "
            "order-statistic intervals Bonferroni-corrected over the five levels; the zero "
            "line is the target.",
        "sfp2_loader_validation":
            "Complete-loader validation matrix. Rows are the curated configurations C0-C6 "
            "and the two negative controls; columns are the predeclared test families. "
            "Decisions are Holm-corrected within each configuration at a familywise "
            "$\\alpha$ of 0.01. A filled circle is a pass, or for NC1 and NC2 a successful "
            "detection of the injected defect: NC1 couples radius to direction while leaving "
            "both one-dimensional marginals bit-identical, and NC2 tests the capped sample "
            "against the uncapped target. Those two rows pass only when the defect is "
            "detected. A cross is a failure, an open triangle an unresolved statistic, and a "
            "small grey dot a test that does not apply - direction uniformity and "
            "radius-direction independence are not properties of the capped conditional "
            "law, which is tested by its own conditional transform instead. Sample sizes are "
            "in the row labels. A cross in the radial column of a released-path row is not "
            "a defect in that row: it is the conditioning result. Where the released "
            "formation discards a non-negligible share of the intended draws, the draws it "
            "keeps are no longer distributed as the target, and the test detects that. The "
            "log path passes the same test on the same configuration.",
        "sfp3_portability_performance":
            "Portability and cost. (a) Per-attempt failure probability of each formation on "
            "every environment that was available, with SPLIT and LOG offset horizontally "
            "inside each platform and the two boundary cases distinguished by marker shape; "
            "downward triangles are one-sided 95 per cent upper limits where no failure was "
            "observed. Environments that were not available on this host are absent from the "
            "panel and recorded in portability.csv with the exact command to run them "
            "elsewhere; they are not drawn as zeros. (b) Paired blockwise ratio of LOG to "
            "SPLIT time per returned sample, with the 95 per cent bootstrap interval over "
            "ten timed blocks, by benchmark case. The black line is parity. Architectures "
            "are not pooled.",
    }
    os.makedirs(FIGURES, exist_ok=True)
    made = {m["figure"] for m in MANIFEST}
    lines = ["# Experiment 6 figure captions", "",
             f"Generated {datetime.now(timezone.utc).isoformat(timespec='seconds')}.", ""]
    for name, text in caps.items():
        if name not in made:
            continue
        lines += [f"## {name}", "", text, ""]
    with open(os.path.join(FIGURES, "captions.md"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    print("  wrote figures/captions.md")


def main() -> int:
    os.makedirs(FIGURES, exist_ok=True)
    selected = "LOG-ID"
    cfg = os.path.join(HERE, "config", "frozen.json")
    if os.path.exists(cfg):
        selected = json.load(open(cfg, encoding="utf-8")).get("log_primitive", selected)

    figure_fp1()
    figure_fp2()
    figure_sfp1(selected)
    figure_sfp2()
    figure_sfp3()
    write_captions(selected)

    with open(os.path.join(FIGURES, "figure_manifest.json"), "w", encoding="utf-8") as fh:
        json.dump(MANIFEST, fh, indent=2)
        fh.write("\n")
    print(f"  wrote figures/figure_manifest.json  ({len(MANIFEST)} figures)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
