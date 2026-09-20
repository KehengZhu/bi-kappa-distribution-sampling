#!/usr/bin/env python3
"""Experiment 7 figures.

Run with:  uv run --project ../../python python make_figures.py [--smoke]

This script plots and nothing else.  It reads only the committed CSVs under ``results/`` --
never a raw binary, never a JSONL counter file -- so every plotted point is recoverable from
a human-readable row, and no statistic is recomputed here.  The imports below are the whole
guarantee: ``csv``, ``json``, ``hashlib``, ``math``, ``os``, ``sys``, numpy and matplotlib.

Output is byte-identical between runs.  ``SOURCE_DATE_EPOCH`` is fixed before matplotlib is
imported, the PDF writer is told to omit its creation date and the SVG writer its date, and
``svg.hashsalt`` is pinned so that generated element ids do not depend on the process.  Six
of Experiment 6's tracked outputs carry a "Generated <utc>" line and its bundle can therefore
never re-verify; nothing written here does.

Two main figures and three supplementary figures, each exported as PDF (canonical), SVG
(editable text) and a 600-dpi PNG preview, with ``figures/captions.md`` and
``figures/figure_manifest.json`` recording the source CSVs and their hashes.  ``--smoke``
writes into ``results/smoke/figures/`` instead, which the repository excludes wholesale, so
a rehearsal can never reach the tracked figure set.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import sys

# Fixed before matplotlib is imported: the PDF and SVG writers read it when they stamp a
# document date, and an unset value would make every run a different file.
os.environ.setdefault("SOURCE_DATE_EPOCH", "1735689600")   # 2025-01-01T00:00:00Z

import matplotlib as mpl                                    # noqa: E402

mpl.use("Agg")

import numpy as np                                          # noqa: E402
from matplotlib import pyplot as plt                        # noqa: E402
from matplotlib.lines import Line2D                         # noqa: E402
from matplotlib.patches import Patch                        # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))

# Repository manuscript style.  Sans-serif at 7-8 pt, hairline axes, no top or right spine,
# text kept as text in both vector formats so that a proof-stage edit does not need the
# script.
mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "Helvetica", "Arial"],
    "font.size": 8,
    "axes.labelsize": 8,
    "axes.titlesize": 8,
    "legend.fontsize": 6.5,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "axes.linewidth": 0.7,
    "xtick.major.width": 0.7,
    "ytick.major.width": 0.7,
    "xtick.minor.width": 0.6,
    "ytick.minor.width": 0.6,
    "lines.linewidth": 1.1,
    "figure.dpi": 200,
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "savefig.facecolor": "white",
    "savefig.transparent": False,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
    "svg.hashsalt": "exp7-confirmatory",
    "axes.spines.top": False,
    "axes.spines.right": False,
})

MM = 1.0 / 25.4
WIDTH_MM = 183.0

# One restrained neutral/blue/vermilion palette.  Line style and marker shape carry the same
# distinction as colour, so every panel survives a grayscale print.
NEUTRAL, BLUE, VERMILION = "0.45", "#1f4e9c", "#d55e00"
STYLE = {
    "QF":        dict(color=NEUTRAL, ls=":", marker="o", mfc="none", ms=3.0, lw=0.9),
    "LEGACY":    dict(color=BLUE, ls="-", marker="o", mfc=BLUE, ms=3.4, lw=1.2),
    "CANDIDATE": dict(color=VERMILION, ls="-", marker="s", mfc=VERMILION, ms=3.2, lw=1.2),
    "floor":     dict(color="black", ls="--", marker="D", mfc="none", ms=3.0, lw=1.0),
}
PASS_C, FAIL_C, GREY_C = BLUE, VERMILION, "0.6"
PRIMARY_TAG = "libcxx"

MANIFEST: list[dict] = []


# ---------------------------------------------------------------------------
def load(ctx: dict, name: str) -> list[dict]:
    path = os.path.join(ctx["results"], name)
    if not os.path.exists(path):
        return []
    with open(path, encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh))


def num(row: dict, key: str, default=float("nan")) -> float:
    v = row.get(key)
    if v in (None, "", "NA"):
        return default
    if v == "inf":
        return math.inf
    if v == "-inf":
        return -math.inf
    if v == "nan":
        return math.nan
    try:
        return float(v)
    except ValueError:
        return default


def truthy(row: dict, key: str) -> bool:
    return str(row.get(key, "")).strip().lower() == "true"


def present(row: dict, key: str) -> bool:
    return row.get(key) not in (None, "", "NA")


def sha256_of(ctx: dict, name: str):
    path = os.path.join(ctx["results"], name)
    if not os.path.exists(path):
        return None
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def script_sha256() -> str:
    with open(os.path.abspath(__file__), "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()


def export(ctx: dict, fig, basename: str, panels: int, conclusion: str, sources: list[str],
           width_mm: float, height_mm: float, interval_definition: str) -> None:
    """PDF, SVG and a 600-dpi PNG, with every date field suppressed.

    ``metadata={"CreationDate": None}`` for the PDF and ``metadata={"Date": None}`` for the
    SVG are what make two runs produce the same bytes; without them each file carries the
    moment it was written and no byte-for-byte regeneration check is possible.
    """
    os.makedirs(ctx["figures"], exist_ok=True)
    for ext, kw in (("pdf", {"metadata": {"CreationDate": None}}),
                    ("svg", {"metadata": {"Date": None}}),
                    ("png", {"dpi": 600, "metadata": {"Software": None}})):
        fig.savefig(os.path.join(ctx["figures"], f"{basename}.{ext}"),
                    bbox_inches="tight", **kw)
    plt.close(fig)
    MANIFEST.append({
        "figure": basename,
        "panels": panels,
        "core_conclusion": conclusion,
        "source_csvs": sources,
        "source_csv_sha256s": {s: sha256_of(ctx, s) for s in sources},
        "script_sha256": script_sha256(),
        "protocol_sha256": ctx["protocol_sha256"],
        "width_mm": width_mm,
        "height_mm": height_mm,
        "interval_definition": interval_definition,
        "formats": ["pdf", "svg", "png@600dpi"],
    })
    print(f"  wrote {os.path.relpath(ctx['figures'], HERE)}/{basename}.{{pdf,svg,png}}")


def panel_label(ax, letter: str, text: str = "") -> None:
    ax.set_title(f"{letter} {text}".rstrip(), loc="left", fontweight="bold", fontsize=8)


def plot_rate_series(ax, xs, rates, los, his, bounds, style, label):
    """A failure-probability curve on log axes, with zero counts drawn as upper bounds.

    A zero count has no place on a logarithmic axis: at zero it disappears, and at its point
    estimate it asserts a precision the sample does not have.  It is drawn instead as a
    downward triangle sitting at the one-sided 95 per cent upper limit, which is what the
    data actually support.
    """
    xs = np.asarray(xs, dtype=float)
    rates = np.asarray(rates, dtype=float)
    los = np.asarray(los, dtype=float)
    his = np.asarray(his, dtype=float)
    obs = ~np.asarray(bounds, dtype=bool) & np.isfinite(rates) & (rates > 0)
    if obs.any():
        ax.errorbar(xs[obs], rates[obs],
                    yerr=[np.maximum(rates[obs] - los[obs], 0),
                          np.maximum(his[obs] - rates[obs], 0)],
                    capsize=1.2, elinewidth=0.6, label=label, **style)
    bound = np.asarray(bounds, dtype=bool) & np.isfinite(his)
    if bound.any():
        ax.plot(xs[bound], his[bound], marker="v", ls="none", color=style["color"],
                ms=3.6, mfc="none", mew=0.9,
                label=None if obs.any() else label)
    return obs


# ---------------------------------------------------------------------------
# FP1 - failure envelope and mechanism decomposition
# ---------------------------------------------------------------------------
def figure_fp1(ctx: dict) -> None:
    rows = [r for r in load(ctx, "failure_envelope.csv")
            if r["scope"] == "pooled" and r["layer"] == "paired"
            and r["tag"] == PRIMARY_TAG]
    floor_rows = load(ctx, "honest_floor.csv")
    if not rows or not floor_rows:
        print("  FP1: no data")
        return

    fig, axes = plt.subplots(1, 2, figsize=(WIDTH_MM * MM, 78 * MM))
    # The axis is set by what was observed.  The analytic floor runs to 1e-309 at the top of
    # the ladder, and letting it set the range would compress every measurement onto one
    # line; it is drawn where it is comparable to an observable rate and runs off the bottom
    # where it is not, which the caption states.
    observed = [num(r, "failure_rate") for r in rows if num(r, "failure_rate") > 0]
    observed += [num(r, "failure_ci_hi") for r in rows
                 if truthy(r, "failure_is_upper_bound") and num(r, "failure_ci_hi") > 0]
    ymin = (min(observed) / 10.0) if observed else 1e-8
    for ax, precision, letter in zip(axes, ("double", "float"), ("a", "b")):
        sub = [r for r in rows if r["precision"] == precision]
        if not sub:
            continue
        for method in ("QF", "LEGACY", "CANDIDATE"):
            ms = sorted([r for r in sub if r["method"] == method],
                        key=lambda r: num(r, "kappa"))
            if not ms:
                continue
            xs = [num(r, "kappa") - 0.5 for r in ms]
            obs = plot_rate_series(
                ax, xs, [num(r, "failure_rate") for r in ms],
                [num(r, "failure_ci_lo") for r in ms],
                [num(r, "failure_ci_hi") for r in ms],
                [truthy(r, "failure_is_upper_bound") for r in ms],
                STYLE[method], method)

        # The analytic floor, as a reference curve rather than as a measurement.  It is a
        # closed form of the type's own range, so it has no interval and is drawn as a line.
        fl = sorted([r for r in floor_rows if r["precision"] == precision],
                    key=lambda r: num(r, "kappa"))
        fx = [num(r, "kappa") - 0.5 for r in fl]
        fy = [num(r, "honest_floor_rate") for r in fl]
        ax.plot(fx, fy, label="analytic honest floor", **STYLE["floor"])
        if any(0 < v < ymin for v in fy):
            ax.annotate("floor continues below the axis",
                        xy=(0.98, 0.02), xycoords="axes fraction", fontsize=6.0,
                        ha="right", va="bottom", color="0.35")

        # Shade the avoidable gap: the distance between the released 1.0.0 form and the
        # floor, wherever both are resolved observations.
        leg = {round(num(r, "kappa"), 12): r for r in sub if r["method"] == "LEGACY"}
        flo = {round(num(r, "kappa"), 12): num(r, "honest_floor_rate") for r in fl}
        ks = sorted(k for k in leg
                    if k in flo and not truthy(leg[k], "failure_is_upper_bound")
                    and num(leg[k], "failure_rate") > flo[k] > 0)
        if len(ks) >= 2:
            ax.fill_between([k - 0.5 for k in ks], [flo[k] for k in ks],
                            [num(leg[k], "failure_rate") for k in ks],
                            color=BLUE, alpha=0.10, lw=0, zorder=0)

        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel(r"$\kappa-1/2$ (dimensionless)")
        panel_label(ax, letter, "double precision" if precision == "double"
                    else "single precision")
    for ax in axes:
        ax.set_ylim(max(ymin * 0.2, 1e-320), 3.0)
    axes[0].set_ylabel("Failure probability per attempt\n(dimensionless)")
    handles = [Line2D([], [], label="QF  $\\sqrt{X_1/X_2}$ (diagnostic)", **STYLE["QF"]),
               Line2D([], [], label="LEGACY  1.0.0, $\\sqrt{X_1}/\\sqrt{X_2}$",
                      **STYLE["LEGACY"]),
               Line2D([], [], label="CANDIDATE  2.0.0, log domain", **STYLE["CANDIDATE"]),
               Line2D([], [], label="analytic honest floor", **STYLE["floor"]),
               Line2D([], [], color="0.3", marker="v", mfc="none", ls="none", ms=3.6,
                      label="one-sided 95% upper limit\n(no failure observed)"),
               Patch(facecolor=BLUE, alpha=0.10, label="avoidable loss (shaded)")]
    axes[1].legend(handles=handles, frameon=False, loc="lower left")
    fig.tight_layout(pad=0.4)
    export(ctx, fig, "fp1_failure_envelope", 2,
           "The log-domain candidate removes the avoidable finite-precision loss of the "
           "released 1.0.0 Gamma-ratio loader and sits on the analytic representability "
           "floor, which no formation can go below.",
           ["failure_envelope.csv", "honest_floor.csv"], WIDTH_MM, 78,
           "two-sided 95% Clopper-Pearson; a zero count is the one-sided 95% upper limit "
           "1 - 0.05^(1/N), drawn as a downward triangle and never as zero; the floor is a "
           "closed form and carries no interval")


# ---------------------------------------------------------------------------
# FP2 - state-dependent failure and its tail consequence
# ---------------------------------------------------------------------------
def fp2_cases(bins: list[dict]) -> list[tuple]:
    """The two conditioning cells to draw, chosen from what P3 actually ran."""
    have = sorted({(r["precision"], round(num(r, "kappa"), 12)) for r in bins})
    want = [("double", 0.51), ("float", 0.55)]
    out = [w for w in want if (w[0], round(w[1], 12)) in
           {(p, k) for p, k in have}]
    if len(out) == 2:
        return out
    dbl = [h for h in have if h[0] == "double"]
    flt = [h for h in have if h[0] == "float"]
    return ([dbl[len(dbl) // 2]] if dbl else []) + ([flt[0]] if flt else [])


def figure_fp2(ctx: dict) -> None:
    bins = [r for r in load(ctx, "conditioning_bins.csv")
            if r["scope"] == "pooled" and r["tag"] == PRIMARY_TAG]
    tails = [r for r in load(ctx, "tail_metrics.csv") if r["tag"] == PRIMARY_TAG]
    if not bins or not tails:
        print("  FP2: no data")
        return
    cases = fp2_cases(bins)
    if not cases:
        print("  FP2: no conditioning cell to draw")
        return

    fig, axes = plt.subplots(2, 2, figsize=(WIDTH_MM * MM, 118 * MM))
    ratio_axes = []
    for col in range(2):
        precision, kappa = cases[min(col, len(cases) - 1)]
        ax = axes[0, col]
        for method in ("LEGACY", "CANDIDATE"):
            sub = [r for r in bins
                   if r["method"] == method and r["precision"] == precision
                   and abs(num(r, "kappa") - kappa) < 1e-9 and num(r, "n_in_bin") >= 20]
            if not sub:
                continue
            sub.sort(key=lambda r: -num(r, "bin_hi"))
            x = [math.sqrt(max(num(r, "bin_hi"), 1e-300)
                           * max(num(r, "bin_lo"), 1e-300 if num(r, "bin_lo") > 0
                                 else num(r, "bin_hi") * 1e-2)) for r in sub]
            y = [num(r, "success_rate") for r in sub]
            lo = [max(num(r, "success_rate") - num(r, "ci_lo"), 0) for r in sub]
            hi = [max(num(r, "ci_hi") - num(r, "success_rate"), 0) for r in sub]
            ax.errorbar(x, y, yerr=[lo, hi], capsize=1.2, elinewidth=0.6, label=method,
                        **STYLE[method])
        ax.set_xscale("log")
        ax.invert_xaxis()
        ax.set_ylim(-0.05, 1.05)
        ax.set_xlabel(r"intended upper-tail probability $q=\Pr(R>r)$")
        panel_label(ax, "ab"[col], rf"{precision}, $\kappa={kappa:g}$")
        if col == 0:
            ax.set_ylabel("Pr(a vector reaches the caller | bin)\n(dimensionless)")
            ax.legend(frameon=False, loc="lower left")

        ax = axes[1, col]
        ratio_axes.append(ax)
        levels = sorted({num(r, "p") for r in tails})
        xpos = np.arange(len(levels))
        for mi, method in enumerate(("LEGACY", "CANDIDATE")):
            sub = {num(r, "p"): r for r in tails
                   if r["method"] == method and r["precision"] == precision
                   and abs(num(r, "kappa") - kappa) < 1e-9}
            st = STYLE[method]
            off = (mi - 0.5) * 0.24
            for kind, mfc, dx in (("returned_tail_ratio_per_attempt", st["color"], -0.055),
                                  ("returned_tail_ratio_conditional_on_success", "none",
                                   0.055)):
                xs, ys, los, his = [], [], [], []
                for i, p in enumerate(levels):
                    r = sub.get(p)
                    if not r or not np.isfinite(num(r, kind)):
                        continue
                    xs.append(xpos[i] + off + dx)
                    ys.append(num(r, kind))
                    los.append(max(num(r, kind) - num(r, f"{kind}_lo"), 0))
                    his.append(max(num(r, f"{kind}_hi") - num(r, kind), 0))
                if xs:
                    ax.errorbar(xs, ys, yerr=[los, his], ls="none", capsize=1.2,
                                elinewidth=0.6, color=st["color"], marker=st["marker"],
                                mfc=mfc, ms=st["ms"], mew=0.9)
        ax.axhline(1.0, color="black", lw=0.7, zorder=0)
        ax.set_xticks(xpos)
        ax.set_xticklabels([f"{p:g}" for p in levels])
        ax.set_xlabel("target percentile $p$ of the intended radius")
        panel_label(ax, "cd"[col], rf"{precision}, $\kappa={kappa:g}$")
        if col == 0:
            ax.set_ylabel("returned tail mass / target tail mass\n(dimensionless)")
            ax.legend(handles=[
                Line2D([], [], color=STYLE["LEGACY"]["color"], marker="o",
                       mfc=STYLE["LEGACY"]["color"], ls="none", ms=3.4, label="LEGACY"),
                Line2D([], [], color=STYLE["CANDIDATE"]["color"], marker="s",
                       mfc=STYLE["CANDIDATE"]["color"], ls="none", ms=3.2,
                       label="CANDIDATE"),
                Line2D([], [], color="0.3", marker="o", mfc="0.3", ls="none", ms=3.4,
                       label="per attempt"),
                Line2D([], [], color="0.3", marker="o", mfc="none", ls="none", ms=3.4,
                       label="conditional on success")],
                frameon=False, loc="best", ncol=2)
    if len(ratio_axes) == 2:
        lo = min(a.get_ylim()[0] for a in ratio_axes)
        hi = max(a.get_ylim()[1] for a in ratio_axes)
        for a in ratio_axes:
            a.set_ylim(lo, hi)
    fig.tight_layout(pad=0.4)
    export(ctx, fig, "fp2_conditioning_tail", 4,
           "Released-path failure is concentrated in the intended tail, so conditioning on a "
           "returned vector changes tail observables; the candidate returns the representable "
           "tail and leaves the honest overflow visible instead of hiding it.",
           ["conditioning_bins.csv", "tail_metrics.csv"], WIDTH_MM, 118,
           "95% Clopper-Pearson for per-bin success rates (one-sided 95% upper limits where "
           "a bin saw no success); 99% cluster-bootstrap percentile intervals over seeds, "
           "10000 resamples, for the tail ratios")


# ---------------------------------------------------------------------------
# SFP1 - scalar validation
# ---------------------------------------------------------------------------
QUANTILE_KEYS = (("q0p5", 0.5), ("q0p9", 0.9), ("q0p99", 0.99), ("q0p999", 0.999),
                 ("q0p9999", 0.9999))
TAIL_KEYS = (("tail0p01", 1e-2), ("tail0p001", 1e-3), ("tail0p0001", 1e-4))


def figure_sfp1(ctx: dict) -> None:
    ecdf = [r for r in load(ctx, "scalar_ecdf.csv")
            if r["method"] == "CANDIDATE" and r["tag"] == PRIMARY_TAG]
    scal = [r for r in load(ctx, "scalar_validation.csv")
            if r["method"] == "CANDIDATE" and r["tag"] == PRIMARY_TAG]
    if not ecdf or not scal:
        print("  SFP1: no data")
        return

    available = sorted({round(num(r, "kappa"), 12) for r in ecdf})
    wanted = [0.5001, 0.51, 0.75, 2.0]
    panel_k = [k for k in wanted if k in available] or available[:4]
    cmap = plt.get_cmap("viridis")
    kcol = {k: cmap(t) for k, t in zip(panel_k, np.linspace(0.05, 0.82, max(len(panel_k),
                                                                            1)))}
    dashes = [(None, None), (4, 2), (1, 1.5), (5, 1.5, 1, 1.5)]
    kdash = {k: dashes[i % len(dashes)] for i, k in enumerate(panel_k)}

    fig, axes = plt.subplots(2, 2, figsize=(WIDTH_MM * MM, 118 * MM))
    for ax, precision, letter in zip(axes[0], ("double", "float"), ("a", "b")):
        for k in panel_k:
            sub = sorted([r for r in ecdf if r["precision"] == precision
                          and abs(num(r, "kappa") - k) < 1e-9],
                         key=lambda r: num(r, "grid_index"))
            if not sub:
                continue
            x = [num(r, "null_cdf") for r in sub]
            y = [num(r, "ecdf_residual") for r in sub]
            ln, = ax.plot(x, y, color=kcol[k], lw=1.0, label=rf"$\kappa={k:g}$")
            if kdash[k][0] is not None:
                ln.set_dashes(list(kdash[k]))
            ax.axhline(num(sub[0], "band_hi"), color=kcol[k], lw=0.5, ls=":")
            ax.axhline(num(sub[0], "band_lo"), color=kcol[k], lw=0.5, ls=":")
        ax.axhline(0.0, color="black", lw=0.7)
        ax.set_xlabel(r"$F_0(Z)$, unit-exponential CDF")
        panel_label(ax, letter, f"{precision}: $Z=-\\log I_W(a,3/2)$")
        if letter == "a":
            ax.set_ylabel("ECDF residual\n(dimensionless)")
            ax.legend(frameon=False, ncol=2)

    # (c) the signed quantile error, with the informativeness classification shown rather
    # than applied silently: a bracket wider than one natural-log unit is declared
    # non-informative in advance and is drawn open.
    ax = axes[1, 0]
    for j, precision in enumerate(("double", "float")):
        for i, k in enumerate(panel_k):
            rs = [r for r in scal if r["precision"] == precision
                  and abs(num(r, "kappa") - k) < 1e-9]
            if not rs:
                continue
            for li, (key, p) in enumerate(QUANTILE_KEYS):
                errs = [num(r, f"{key}_error") for r in rs]
                errs = [e for e in errs if np.isfinite(e)]
                if not errs:
                    continue
                infm = any(truthy(r, f"{key}_informative") for r in rs)
                x = li + (j - 0.5) * 0.30 + (i - (len(panel_k) - 1) / 2) * 0.055
                ax.plot([x], [float(np.mean(errs))], ls="none", ms=3.2,
                        color=kcol[k], marker="o" if precision == "double" else "^",
                        mfc=kcol[k] if infm else "none", mew=0.9)
    ax.axhline(0.0, color="black", lw=0.7)
    ax.set_yscale("symlog", linthresh=1.0)
    ax.set_xticks(range(len(QUANTILE_KEYS)))
    ax.set_xticklabels([f"{p:g}" for _, p in QUANTILE_KEYS])
    ax.set_xlabel("target percentile $p$")
    ax.set_ylabel(r"$\widehat{\log r_p}-\log r_p$" "\n(natural log units, symlog)")
    panel_label(ax, "c", "log-radius quantile error")
    ax.legend(handles=[
        Line2D([], [], color="0.3", marker="o", mfc="0.3", ls="none", ms=3.2,
               label="double"),
        Line2D([], [], color="0.3", marker="^", mfc="0.3", ls="none", ms=3.2,
               label="float"),
        Line2D([], [], color="0.3", marker="o", mfc="none", ls="none", ms=3.2,
               label="non-informative bracket")], frameon=False, loc="best")

    # (d) upper-tail mass: observed exceedances against the count the null predicts.  This is
    # the statistic with power against survivor conditioning where the bulk-weighted ones
    # have none, so it is shown rather than summarized.  The count is the one the family
    # tests: resolved exceedances PLUS the draws whose Z never resolved, over attempts, since
    # an unresolved draw lies above every threshold by construction.
    ax = axes[1, 1]
    for mi, (key, q0) in enumerate(TAIL_KEYS):
        xs, ys = [], []
        for r in sorted(scal, key=lambda r: num(r, "kappa")):
            if r["precision"] != "double":
                continue
            exp = num(r, f"{key}_expected_per_attempt")
            obs = num(r, f"{key}_observed_total")
            if not (np.isfinite(exp) and exp > 0 and np.isfinite(obs)):
                continue
            xs.append(num(r, "kappa") - 0.5)
            ys.append(obs / exp)
        if xs:
            ax.plot(xs, ys, ls="none", ms=3.2, mew=0.9,
                    marker=["o", "s", "^"][mi],
                    color=[NEUTRAL, BLUE, VERMILION][mi],
                    mfc=["none", BLUE, "none"][mi],
                    label=rf"$q_0={q0:g}$")
    ax.axhline(1.0, color="black", lw=0.7, zorder=0)
    ax.set_xscale("log")
    ax.set_xlabel(r"$\kappa-1/2$ (dimensionless), double")
    ax.set_ylabel("observed / expected exceedances\n(dimensionless)")
    panel_label(ax, "d", "upper-tail mass above $z_0=-\\log q_0$")
    ax.legend(frameon=False, loc="best")

    fig.tight_layout(pad=0.4)
    export(ctx, fig, "sfp1_scalar_validation", 4,
           "The candidate samples the intended radial law on a stable scale across the whole "
           "frozen ladder, including the shapes where the order-statistic brackets are too "
           "wide to be evidence and the upper-tail statistic carries the certification.",
           ["scalar_ecdf.csv", "scalar_validation.csv"], WIDTH_MM, 118,
           "dotted lines are simultaneous Kolmogorov bands at the frozen F1 level; quantile "
           "brackets are exact binomial order-statistic intervals, Bonferroni-corrected over "
           "the five levels")


# ---------------------------------------------------------------------------
# SFP2 - the complete-loader validation matrix
# ---------------------------------------------------------------------------
SYMBOL = {"pass": ("o", PASS_C, PASS_C), "fail": ("X", FAIL_C, FAIL_C),
          "ungated": ("o", GREY_C, "none"), "na": (".", GREY_C, GREY_C)}


def sfp2_tests(rows: list[dict]) -> list[str]:
    """The F5 members, in the order `loader_validation.csv` writes them.

    Read off the header rather than listed here: amendment 1.3.0 added one pair of upper-tail
    members per threshold in `F5_tail.q0`, and a figure that hard-codes the family's
    membership is a second place for the protocol to be recorded.
    """
    names = [k[len("applies_"):] for k in rows[0] if k.startswith("applies_")]
    return names


def sfp2_label(t: str) -> str:
    if t.startswith("tail_"):
        _, slug, kind = t.split("_", 2)
        return f"tail {slug.replace('m', '-')} {kind}"
    return t.replace("_", " ")


def figure_sfp2(ctx: dict) -> None:
    rows = [r for r in load(ctx, "loader_validation.csv")
            if r["scope"] == "pooled" and r["tag"] == PRIMARY_TAG]
    if not rows:
        print("  SFP2: no data")
        return
    order = {c: i for i, c in enumerate(["C0", "C1", "C2", "C3", "C4", "C5", "C6"])}
    rows.sort(key=lambda r: (order.get(r["case"], 99), r["method"]))
    SFP2_TESTS = sfp2_tests(rows)
    labels = []
    fig, ax = plt.subplots(figsize=(WIDTH_MM * MM, 104 * MM))
    for y, r in enumerate(rows):
        gated = truthy(r, "in_family")
        for x, t in enumerate(SFP2_TESTS):
            if not truthy(r, f"applies_{t}"):
                marker, edge, face = SYMBOL["na"]
                ms = 2.2
            else:
                # The Holm decision is the family's, taken in analyze.py over every cell and
                # test jointly and written into the row; re-deriving a threshold from an
                # alpha and a count of tests here would be a second copy of the rule.
                if not gated:
                    marker, edge, face = SYMBOL["ungated"]
                elif truthy(r, f"f5_holm_rejected_{t}"):
                    marker, edge, face = SYMBOL["fail"]
                else:
                    marker, edge, face = SYMBOL["pass"]
                ms = 5.5
            ax.plot(x, y, marker=marker, color=edge, mfc=face, ms=ms, mew=1.0, ls="none")
        cond = "conditional" if truthy(r, "conditional") else "unconditional"
        lf = num(r, "loss_fraction")
        lft = f", loss {lf:.2g}" if truthy(r, "conditional") and np.isfinite(lf) else ""
        labels.append(f"{r['case']} {r['method']}  (n={int(num(r, 'n_analyzed', 0)):,}; "
                      f"{cond}{lft})")
    ax.set_xticks(range(len(SFP2_TESTS)))
    ax.set_xticklabels([sfp2_label(t) for t in SFP2_TESTS], rotation=35, ha="right")
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels(labels, fontsize=6.5)
    ax.set_xlim(-0.6, len(SFP2_TESTS) - 0.4)
    ax.set_ylim(-0.7, len(rows) - 0.3)
    ax.invert_yaxis()
    ax.grid(True, color="0.92", lw=0.5)
    ax.set_axisbelow(True)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_visible(False)
    ax.tick_params(length=0)
    handles = [Line2D([], [], marker=m, color=e, mfc=f, ls="none", ms=5.5, label=lab)
               for lab, (m, e, f) in
               (("pass (candidate, in family F5)", SYMBOL["pass"]),
                ("fail", SYMBOL["fail"]),
                ("comparator, reported but not gated", SYMBOL["ungated"]),
                ("not a property of this cell's law", SYMBOL["na"]))]
    ax.legend(handles=handles, frameon=False, ncol=4, loc="upper center",
              bbox_to_anchor=(0.5, 1.11))
    fig.tight_layout(pad=0.4)
    export(ctx, fig, "sfp2_loader_validation", 1,
           "The candidate stays correct after direction sampling, anisotropic scaling, "
           "rotation into an arbitrary field frame and cap conditioning; every cell whose "
           "statistics ran on a radius-filtered subsample is marked conditional with its "
           "loss fraction.",
           ["loader_validation.csv"], WIDTH_MM, 104,
           "Holm-corrected jointly over all cells and tests at the frozen F5 familywise "
           "alpha; a conditional cell cannot support the fidelity claim on its own")


# ---------------------------------------------------------------------------
# SFP3 - portability and cost
# ---------------------------------------------------------------------------
def figure_sfp3(ctx: dict) -> None:
    port = load(ctx, "portability.csv")
    perf = [r for r in load(ctx, "performance.csv") if r["scope"] == "case_ratio"]
    if not port and not perf:
        print("  SFP3: no data")
        return

    fig, axes = plt.subplots(1, 2, figsize=(WIDTH_MM * MM, 84 * MM))

    # (a) Environment coverage and the bitwise comparison.  An environment that did not run
    # is drawn as an open triangle on a separate "not available" level, never as a zero and
    # never omitted: a missing platform must not be readable as a passing one.
    ax = axes[0]
    envs = [r for r in port if r["kind"] == "environment"]
    pairs = [r for r in port if r["kind"] == "cross_stdlib_bitwise"]
    envs.sort(key=lambda r: (not truthy(r, "completed"), str(r.get("arch")),
                             str(r.get("stdlib"))))
    for i, r in enumerate(envs):
        done = truthy(r, "completed") and truthy(r, "available")
        ax.plot([i], [1.0 if done else 0.0], ls="none", ms=5.0, mew=1.0,
                marker="o" if done else "v", color=PASS_C if done else GREY_C,
                mfc=PASS_C if done else "none")
    ax.set_yticks([0.0, 1.0])
    ax.set_yticklabels(["not available\n(gate stays open)", "native, completed"])
    ax.set_ylim(-0.6, 1.6)
    ax.set_xticks(range(len(envs)))
    ax.set_xticklabels([f"{r.get('arch')}\n{r.get('stdlib')}" for r in envs], fontsize=6.0)
    ax.set_xlim(-0.7, max(len(envs) - 0.3, 0.5))
    panel_label(ax, "a", "environment coverage")
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.legend(handles=[
        Line2D([], [], marker="o", color=PASS_C, mfc=PASS_C, ls="none", ms=5.0,
               label="native, completed"),
        Line2D([], [], marker="v", color=GREY_C, mfc="none", ls="none", ms=5.0,
               label="not available on this host")],
        frameon=False, loc="upper center", bbox_to_anchor=(0.5, 1.02), ncol=2)
    txt = []
    for r in pairs:
        n = int(num(r, "n_common", 0))
        pred = ("bitwise equality" if "bitwise" in str(r.get("prediction", ""))
                else "a difference")
        got = (f"all {n} configurations agree" if truthy(r, "identical")
               else f"{int(num(r, 'n_differences', 0))} of {n} configurations differ")
        txt.append(f"{r.get('label')}\n    predicted {pred}; observed {got}")
    if txt:
        ax.text(0.0, -0.40, "\n".join(txt), transform=ax.transAxes, fontsize=6.0,
                va="top", ha="left")

    # (b) Cost.
    ax = axes[1]
    perf = [r for r in perf if r["tag"] == PRIMARY_TAG] or perf
    perf.sort(key=lambda r: str(r["case"]))
    bound = num(perf[0], "pre_registered_bound") if perf else float("nan")
    for i, r in enumerate(perf):
        med = num(r, "ratio_candidate_over_legacy")
        lo, hi = num(r, "ratio_lo"), num(r, "ratio_hi")
        ax.errorbar([i], [med], yerr=[[max(med - lo, 0)], [max(hi - med, 0)]],
                    ls="none", capsize=1.5, elinewidth=0.7, marker="s", ms=3.6,
                    color=VERMILION, mfc=VERMILION)
    ax.axhline(1.0, color="black", lw=0.7, zorder=0)
    if np.isfinite(bound):
        ax.axhline(bound, color=NEUTRAL, lw=0.9, ls="--", zorder=0)
        ax.text(0.02, bound, f" pre-registered bound {bound:g}x", fontsize=6.0,
                va="bottom", ha="left", transform=ax.get_yaxis_transform())
    tick_labels = []
    for r in perf:
        cap_text = "uncapped" if not present(r, "cap") else f"cap {num(r, 'cap'):g}"
        tick_labels.append(f"{r['case']}\n{r.get('precision')}, "
                           f"$\\kappa$={num(r, 'kappa'):g}\n{cap_text}")
    ax.set_xticks(range(len(perf)))
    ax.set_xticklabels(tick_labels, fontsize=5.5)
    ax.set_ylabel("CANDIDATE / LEGACY time per returned sample\n(dimensionless)")
    panel_label(ax, "b", "cost")
    fig.tight_layout(pad=0.4)
    export(ctx, fig, "sfp3_portability_cost", 2,
           "Within one architecture the candidate's output is bit-for-bit independent of the "
           "standard library while the comparator's is not, and the mitigation costs less "
           "than the pre-registered bound; the architectures that were not available are "
           "shown as absent rather than as passing.",
           ["portability.csv", "performance.csv"], WIDTH_MM, 84,
           "99% cluster-bootstrap percentile interval over seeds, 10000 resamples, for the "
           "paired blockwise time ratio; the bitwise comparison is exact and has no interval")


# ---------------------------------------------------------------------------
def write_captions(ctx: dict) -> None:
    caps = {
        "fp1_failure_envelope":
            "Failure envelope and mechanism decomposition of the Gamma-ratio bi-Kappa "
            "radius. Probability that one attempt fails to return a finite three-vector, "
            "against $\\kappa-1/2$, in (a) double and (b) single precision, for the released "
            "1.0.0 split formation $\\sqrt{X_1}/\\sqrt{X_2}$ (LEGACY), the 2.0.0 log-domain "
            "candidate, and the quotient-first formation $\\sqrt{X_1/X_2}$, which is a "
            "diagnostic and carries no claim. An attempt is one draw in the canonical "
            "isotropic, unrotated, uncapped configuration with "
            "$\\theta_\\perp=\\theta_\\parallel$ and $\\hat b=\\hat z$; curves are pooled "
            "over the frozen seed block, with seed-level rows retained in "
            "`failure_envelope.csv`. The dashed black curve is the analytic honest "
            "representability floor, the closed-form probability that at least one component "
            "of the intended vector has no representation in the run's floating-point type; "
            "it is derived in `config/honest_floor.md`, is a property of the type rather "
            "than of a method, and no formation can go below it. Error bars are two-sided 95 "
            "per cent Clopper-Pearson intervals; a configuration with no observed failure is "
            "drawn as a downward triangle at its one-sided 95 per cent upper limit and never "
            "at zero. Shading marks the avoidable loss, the gap between the released "
            "formation and the floor, wherever both are resolved observations.",
        "fp2_conditioning_tail":
            "State-dependent failure and its tail consequence. (a, b) Probability that an "
            "attempt delivers a three-vector to the caller, as a function of the upper-tail "
            "probability $q=\\Pr(R>r)$ of the radius that attempt intended to produce, so "
            "that the horizontal axis runs from the body of the distribution on the left to "
            "the extreme tail on the right. Bars are 95 per cent Clopper-Pearson intervals "
            "and bins holding fewer than twenty attempts are not drawn. (c, d) Tail "
            "retention at the pre-registered percentiles of the intended radius. Filled "
            "markers are the ratio measured per attempt, in which an attempt that failed "
            "stays outside the returned mass; open markers are the ratio conditional on "
            "success, the law a caller actually receives, and for independent attempts it is "
            "also the law produced by silently redrawing until success. The two are "
            "different quantities and are labelled as such throughout "
            "`tail_metrics.csv`. Where a method's per-attempt law carries a failure atom its "
            "upper quantiles are undefined, so the figure and the source data report tail "
            "MASS above the target quantile rather than a quantile error. Bars are 99 per "
            "cent cluster-bootstrap percentile intervals over the seed block from 10\\,000 "
            "resamples; the black line marks the target itself.",
        "sfp1_scalar_validation":
            "Scalar validation of the 2.0.0 candidate. (a, b) Residuals of the empirical CDF "
            "of $Z=-\\log I_W(a,3/2)$ from the unit-exponential CDF, where "
            "$W=X_2/(X_1+X_2)\\sim\\mathrm{Beta}(a,3/2)$ and $a=\\kappa-1/2$; $Z$ is the "
            "diagnostic of record because $W$ itself rounds to zero at the smallest shapes "
            "while $\\log W$ does not. Dotted lines are simultaneous Kolmogorov bands at the "
            "frozen familywise level of family F1. (c) Error in the empirical $\\log R$ "
            "quantile relative to the exact target quantile, averaged over the seed block; "
            "open markers are cells whose order-statistic bracket is wider than one natural-"
            "log unit and which the protocol declares non-informative in advance, since an "
            "interval admitting a multiplicative error of $e^{8951}$ cannot distinguish a "
            "correct sampler from a wrong one. Those shapes are certified by the families "
            "that retain power at any shape, not by coverage. (d) Observed exceedances above "
            "$z_0=-\\log q_0$ divided by the count the null predicts, in double precision. "
            "The count is the one family F4 tests: the resolved draws above $z_0$ plus the "
            "draws whose $Z$ never resolved, over the attempts made, because a draw that "
            "could not be resolved lies above every threshold and counting only the resolved "
            "ones would condition on resolvability, which is itself monotone in the tail. "
            "This is the statistic that retains power against survivor conditioning at a "
            "loss fraction of $10^{-4}$, where the bulk-weighted statistics have none; the "
            "measured power of each is in `config/power_study.json` and was computed before "
            "the run.",
        "sfp2_loader_validation":
            "Complete-loader validation matrix. Rows are the curated configurations C0-C6 "
            "pooled over the frozen seed block; columns are the members of family F5 -- the "
            "five frozen with the protocol, then an exceedance count and an excess "
            "goodness-of-fit at each upper-tail threshold $q_0$. Decisions are "
            "Holm-corrected jointly over all cells and tests at the "
            "frozen familywise level, not within each cell. A small grey dot marks a test "
            "that is not a property of that cell's law rather than one that passed: under a "
            "cap the accepted set couples the radius to the direction and is not "
            "rotationally symmetric in the azimuth, so direction uniformity, independence, "
            "frame invariance and the upper-tail members do not apply there -- the accepted "
            "radius is truncated at a direction-dependent bound, so the count above $z_0$ is "
            "not binomial -- and the bounded law is tested by its "
            "own conditional transform instead. An exceedance count includes the draws the "
            "analysis could not resolve, each of which lies above every threshold, so the "
            "count also reads on a cell whose own loss fraction exceeds $q_0$; the resolved "
            "and unresolved parts of every count are published separately in "
            "`loader_validation.csv`. The frame test recovers the direction "
            "through a field-aligned basis re-derived independently of the loader's, so it "
            "is a genuine check of the rotation rather than a tautology. Open grey circles "
            "are the 1.0.0 comparator: its p-values are reported but it is not gated, "
            "because where the released formation discards a non-negligible share of the "
            "intended draws the draws it keeps are no longer distributed as the target, and "
            "a test detecting that is the conditioning result rather than a defect. Sample "
            "sizes, and the loss fraction of every cell whose statistics ran on a "
            "radius-filtered subsample, are in the row labels and in "
            "`loader_validation.csv`.",
        "sfp3_portability_cost":
            "Portability and cost. (a) Every environment the portability matrix names. An "
            "environment that ran natively and completed sits on the upper level; one that "
            "was not available on this host sits on the lower level as an open triangle, "
            "with the reason and the exact command to run it elsewhere recorded in "
            "`portability.csv` and `results/portability_remote.md`. An absent environment is "
            "never drawn as a zero and never omitted, so it cannot be read as a passing one; "
            "the gate stays open until it runs, and emulated execution is recorded as "
            "corroborating evidence and excluded from the decision. The annotation states "
            "the two predictions being tested within one architecture: the candidate's "
            "output is predicted to be bit-for-bit equal across standard libraries, because "
            "its Gamma, normal and uniform primitives live in the header and its stream is a "
            "function of the engine alone, while the comparator's is predicted to differ, "
            "because it draws from the standard library. (b) Paired blockwise ratio of "
            "candidate to comparator time per returned sample, by benchmark case, with the "
            "99 per cent cluster-bootstrap interval over the performance seed block. The "
            "solid line is parity and the dashed line the pre-registered bound. "
            "Architectures are never pooled.",
    }
    os.makedirs(ctx["figures"], exist_ok=True)
    L = ["# Experiment 7 figure captions", "",
         "Every figure is generated from the CSVs under `results/` by `make_figures.py`, "
         "which reads no raw binary and recomputes no statistic. This file carries no "
         "generation timestamp: `make reverify` regenerates every derived artifact and diffs "
         "it, which a wall-clock stamp would make impossible.", "",
         f"Protocol `config/protocol.json` SHA-256 `{ctx['protocol_sha256']}`.", "",
         "A caption is a property of the figure's definition rather than of the run that "
         "drew it, so all five are written here; `figure_manifest.json` records which "
         "figures this run actually produced, with the SHA-256 of every source CSV behind "
         "them.", ""]
    for name, text in caps.items():
        L += [f"## {name}", "", text, ""]
    with open(os.path.join(ctx["figures"], "captions.md"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(L) + "\n")
    print(f"  wrote {os.path.relpath(ctx['figures'], HERE)}/captions.md")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--smoke", action="store_true",
                    help="plot results/smoke/ into figures/smoke/")
    args = ap.parse_args()
    results = os.path.join(HERE, "results", "smoke") if args.smoke else os.path.join(
        HERE, "results")
    # A smoke rehearsal writes its figures inside results/smoke/ rather than into
    # figures/: the repository excludes experiments/*/results/smoke/ wholesale, so a
    # rehearsal cannot reach the tracked figure set, and `make SMOKE=1 reverify` -- which
    # diffs results/ -- then checks the figures for byte-identical regeneration as well.
    figures = os.path.join(HERE, "results", "smoke", "figures") if args.smoke \
        else os.path.join(HERE, "figures")
    if not os.path.isdir(results):
        print(f"make_figures.py: {os.path.relpath(results, HERE)} does not exist; run "
              f"`make {'SMOKE=1 ' if args.smoke else ''}analyze` first", file=sys.stderr)
        return 3
    sha = ""
    summary = os.path.join(results, "exp7_results.json")
    if os.path.exists(summary):
        with open(summary, encoding="utf-8") as fh:
            sha = json.load(fh).get("protocol_sha256", "")
    ctx = {"results": results, "figures": figures, "protocol_sha256": sha}
    expected = ("failure_envelope.csv", "conditioning_bins.csv", "tail_metrics.csv",
                "scalar_validation.csv", "scalar_ecdf.csv", "loader_validation.csv",
                "portability.csv", "performance.csv", "honest_floor.csv")
    absent = [n for n in expected if not os.path.exists(os.path.join(results, n))]
    if len(absent) == len(expected):
        print(f"make_figures.py: {os.path.relpath(results, HERE)} holds none of the source "
              f"CSVs the figures are drawn from; run "
              f"`make {'SMOKE=1 ' if args.smoke else ''}analyze` first", file=sys.stderr)
        return 3
    if absent:
        print(f"make_figures.py: missing source data {absent}; the figures that need them "
              "are not drawn", file=sys.stderr)
    os.makedirs(figures, exist_ok=True)

    figure_fp1(ctx)
    figure_fp2(ctx)
    figure_sfp1(ctx)
    figure_sfp2(ctx)
    figure_sfp3(ctx)
    write_captions(ctx)

    with open(os.path.join(figures, "figure_manifest.json"), "w", encoding="utf-8") as fh:
        json.dump(MANIFEST, fh, indent=2, sort_keys=True)
        fh.write("\n")
    print(f"  wrote {os.path.relpath(figures, HERE)}/figure_manifest.json "
          f"({len(MANIFEST)} figures)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
