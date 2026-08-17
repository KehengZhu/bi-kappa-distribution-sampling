"""Experiment 4 analysis — finite-precision / low-kappa audit.

Aggregates the JSONL emitted by ``exp4_probe.exe`` and validates the log-domain
radial path against the exact law.

Numerical conventions inherited from Experiment 1, and they are not optional:

* never form ``T = R**2`` at low kappa -- it overflows for ``R > 1.3e154`` even where
  ``R`` is perfectly representable, silently discarding the heavy-tail draws that the
  diagnostic exists to test;
* never use ``Y = T/(1+T)``; for small ``kappa - 1/2`` its mass piles up against 1.0
  where no relative resolution remains and a KS test measures rounding, not the sampler;
* use ``W = 1/(1+T) ~ Beta(kappa - 1/2, 3/2)``, computed here as ``expit(-2 log R)``
  straight from the log-domain radius, so ``T`` is never materialized at all.

Run:  uv run --project ../../python python exp4_analyze.py
"""

from __future__ import annotations

import json
import math
import pathlib
import platform
import subprocess
import sys
from collections import defaultdict

import numpy as np
import scipy
from scipy import stats
from scipy.special import expit

HERE = pathlib.Path(__file__).resolve().parent
RAW = HERE / "raw"
RESULTS = HERE / "results"

# Quantile probes, compared in log R so that no heavy-tail draw is ever squared.
QUANTILES = (0.5, 0.9, 0.99, 0.999)

# Fixed before any result was seen.  Not to be adjusted afterwards.
KS_ALPHA = 0.01


# --------------------------------------------------------------------------- #
# provenance
# --------------------------------------------------------------------------- #
def git_state() -> dict:
    def run(*args):
        try:
            return subprocess.check_output(args, cwd=HERE, text=True).strip()
        except Exception:
            return "unknown"

    return {
        "commit": run("git", "rev-parse", "HEAD"),
        "dirty": bool(run("git", "status", "--porcelain")),
    }


def environment() -> dict:
    return {
        "python": sys.version.split()[0],
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "platform": platform.platform(),
        "machine": platform.machine(),
        "git": git_state(),
    }


# --------------------------------------------------------------------------- #
# loading
# --------------------------------------------------------------------------- #
def load_jsonl(pattern: str) -> list[dict]:
    rows: list[dict] = []
    for path in sorted(RAW.glob(pattern)):
        with path.open() as fh:
            for line in fh:
                line = line.strip()
                if line:
                    rows.append(json.loads(line))
    return rows


def group(rows, keys):
    out = defaultdict(list)
    for r in rows:
        out[tuple(r[k] for k in keys)].append(r)
    return out


# --------------------------------------------------------------------------- #
# Q1 / Q2 / Q3 -- failure envelope, split-sqrt effect, log-domain headroom
# --------------------------------------------------------------------------- #
def summarize_variates(rows) -> list[dict]:
    """Per (stdlib, precision, kappa), pooled over seeds.

    The three loss categories are mutually exclusive by construction and answer
    three different questions:

    ``honest``          the mathematical radius exceeds the largest representable
                        number at this precision.  No reformulation can fix this.
    ``spurious_ratio``  the radius was representable, but ``sqrt(x1/x2)`` lost it
                        because the *ratio* overflowed before the square root.
                        Fixed by ``sqrt(x1)/sqrt(x2)``.
    ``spurious_split``  the radius was representable, but ``sqrt(x1)/sqrt(x2)``
                        still lost it because the denominator Gamma variate itself
                        underflowed to exactly zero.  Only a log-domain
                        construction can recover these.
    """
    out = []
    for (stdlib, precision, kappa), rs in sorted(
        group(rows, ["stdlib", "precision", "kappa"]).items(),
        key=lambda kv: (kv[0][0], kv[0][1], kv[0][2]),
    ):
        n = sum(r["n"] for r in rs)
        agg = {
            k: sum(r[k] for r in rs)
            for k in (
                "x2_zero",
                "x2_subnormal",
                "ratio_nonfinite",
                "split_nonfinite",
                "recovered_by_split",
                "honest_overflow",
                "spurious_ratio",
                "spurious_split",
            )
        }
        out.append(
            {
                "stdlib": stdlib,
                "precision": precision,
                "kappa": kappa,
                "shape_x2": kappa - 0.5,
                "n_draws": n,
                "n_seeds": len(rs),
                "counts": agg,
                "frac": {k: v / n for k, v in agg.items()},
                # finite fraction under each formation
                "finite_frac_ratio": 1.0 - agg["ratio_nonfinite"] / n,
                "finite_frac_split": 1.0 - agg["split_nonfinite"] / n,
                # what a log-domain construction could achieve: only the honest
                # overflows would remain
                "finite_frac_logdomain_bound": 1.0 - agg["honest_overflow"] / n,
                "max_finite_log_r": max(r["max_finite_log_r"] for r in rs),
                "max_log_r_true": max(r["max_log_r"] for r in rs),
            }
        )
    return out


def summarize_released(rows) -> list[dict]:
    out = []
    for (stdlib, precision, kappa), rs in sorted(
        group(rows, ["stdlib", "precision", "kappa"]).items(),
        key=lambda kv: (kv[0][0], kv[0][1], kv[0][2]),
    ):
        n = sum(r["n"] for r in rs)
        nf = sum(r["nonfinite"] for r in rs)
        out.append(
            {
                "stdlib": stdlib,
                "precision": precision,
                "kappa": kappa,
                "n_draws": n,
                "n_seeds": len(rs),
                "nonfinite": nf,
                "nonfinite_frac": nf / n,
                "thrown": sum(r["thrown"] for r in rs),
                "max_finite_log_r": max(r["max_finite_log_r"] for r in rs),
            }
        )
    return out


def summarize_capped(rows) -> list[dict]:
    out = []
    for (stdlib, precision, kappa, lam), rs in sorted(
        group(rows, ["stdlib", "precision", "kappa", "lambda"]).items(),
        key=lambda kv: (kv[0][0], kv[0][1], kv[0][2], kv[0][3]),
    ):
        n = sum(r["n"] for r in rs)
        att = sum(r["attempts"] for r in rs)
        out.append(
            {
                "stdlib": stdlib,
                "precision": precision,
                "kappa": kappa,
                "lambda": lam,
                "n_accepted": n,
                "attempts": att,
                "acceptance": n / att if att else float("nan"),
                # non-finite draws generated inside the rejection loop
                "nonfinite_redrawn": sum(r["nonfinite_redrawn"] for r in rs),
                "nonfinite_redrawn_per_attempt": sum(r["nonfinite_redrawn"] for r in rs) / att
                if att
                else float("nan"),
                # non-finite draws that actually reached the caller
                "nonfinite_returned": sum(r["nonfinite_returned"] for r in rs),
                "exhausted": sum(r["exhausted"] for r in rs),
            }
        )
    return out


# --------------------------------------------------------------------------- #
# Q3 validation -- does the log-domain path sample the intended law?
# --------------------------------------------------------------------------- #
def validate_logdomain() -> list[dict]:
    """More than one diagnostic, as required.

    ``W = 1/(1+T) ~ Beta(kappa-1/2, 3/2)`` is obtained directly from log R via the
    logistic function, so ``T`` is never formed.  Reported alongside it are robust
    quantiles of ``log R`` against their exact values, which are meaningful even
    where no moment exists.
    """
    out = []
    for path in sorted(RAW.glob("logr_k*_double.bin")):
        stem = path.stem  # logr_k0.55_s4001_double
        parts = stem.split("_")
        kappa = float(parts[1][1:])
        seed = int(parts[2][1:])
        log_r = np.fromfile(path, dtype=np.float64)
        if log_r.size == 0:
            continue

        a = kappa - 0.5
        # W = 1/(1+T) = 1/(1+exp(2 log R)) = expit(-2 log R).  Stable at both ends,
        # and it never materializes T.
        w = expit(-2.0 * log_r)

        ks = stats.kstest(w, "beta", args=(a, 1.5))
        cvm = stats.cramervonmises(w, "beta", args=(a, 1.5))

        # Exact log-R quantiles from the same Beta law, again without forming T:
        #   T = (1-W)/W,   log R = (log(1-W) - log(W)) / 2
        emp, exact = [], []
        for q in QUANTILES:
            wq = stats.beta.ppf(1.0 - q, a, 1.5)  # W decreases as R increases
            exact.append(0.5 * (math.log1p(-wq) - math.log(wq)))
            emp.append(float(np.quantile(log_r, q)))

        out.append(
            {
                "kappa": kappa,
                "seed": seed,
                "n": int(log_r.size),
                "ks_stat": float(ks.statistic),
                "ks_p": float(ks.pvalue),
                "ks_sqrt_n_d": float(ks.statistic * math.sqrt(log_r.size)),
                "cvm_stat": float(cvm.statistic),
                "cvm_p": float(cvm.pvalue),
                "quantiles": {
                    "levels": list(QUANTILES),
                    "empirical_log_r": emp,
                    "exact_log_r": exact,
                    "abs_diff": [abs(e - x) for e, x in zip(emp, exact)],
                },
                "pass": bool(ks.pvalue > KS_ALPHA),
            }
        )
    return out


# --------------------------------------------------------------------------- #
# Q4 -- supported-range table
# --------------------------------------------------------------------------- #
def supported_range(variates_summary) -> list[dict]:
    """(kappa, precision, implementation) -> status.

    The status is a statement about the CURRENT released code, i.e. about the
    ``sqrt(x1)/sqrt(x2)`` formation.  Thresholds are set here, before the numbers
    are read, and are not adjusted afterwards.
    """
    rows = []
    for s in variates_summary:
        loss = s["frac"]["split_nonfinite"]
        if loss == 0.0:
            status = "supported"
        elif loss < 1e-5:
            status = "supported (rare loss)"
        elif loss < 1e-2:
            status = "degraded"
        else:
            status = "unsupported"
        rows.append(
            {
                "kappa": s["kappa"],
                "precision": s["precision"],
                "stdlib": s["stdlib"],
                "finite_fraction": s["finite_frac_split"],
                "spurious_failure_fraction": s["frac"]["spurious_split"],
                "honest_overflow_fraction": s["frac"]["honest_overflow"],
                "recoverable_by_log_domain": s["frac"]["spurious_split"],
                "status": status,
            }
        )
    return rows


# --------------------------------------------------------------------------- #
# reporting
# --------------------------------------------------------------------------- #
def fmt(x: float) -> str:
    if x == 0:
        return "0"
    if x < 1e-4:
        return f"{x:.2e}"
    return f"{x:.5f}".rstrip("0").rstrip(".")


def write_table(res: dict) -> None:
    L: list[str] = []
    L.append("# Experiment 4 — supported numerical range of the released bi-Kappa sampler")
    L.append("")
    L.append(
        "All fractions are pooled over "
        f"{res['variates'][0]['n_seeds']} seeds "
        f"× {res['variates'][0]['n_draws'] // res['variates'][0]['n_seeds']:,} draws "
        "per configuration, uncapped."
    )
    L.append("")

    L.append("## Q1 + Q2 — failure envelope and the effect of `sqrt(x1)/sqrt(x2)`")
    L.append("")
    L.append(
        "`ratio` is the pre-fix formation `sqrt(x1/x2)`; `split` is the current "
        "formation `sqrt(x1)/sqrt(x2)`. Both columns are evaluated on **identical "
        "variates**, so the difference is the fix and nothing else."
    )
    L.append("")
    for stdlib in sorted({s["stdlib"] for s in res["variates"]}):
        for precision in ("float", "double"):
            sel = [
                s
                for s in res["variates"]
                if s["stdlib"] == stdlib and s["precision"] == precision
            ]
            if not sel:
                continue
            L.append(f"### {stdlib}, {precision}")
            L.append("")
            L.append(
                "| κ | shape(x2) | lost, `ratio` | lost, `split` | recovered by fix | "
                "honest overflow | still spurious |"
            )
            L.append("|---|---|---|---|---|---|---|")
            for s in sorted(sel, key=lambda r: r["kappa"]):
                f = s["frac"]
                L.append(
                    f"| {s['kappa']:g} | {s['shape_x2']:g} | "
                    f"{fmt(f['ratio_nonfinite'])} | {fmt(f['split_nonfinite'])} | "
                    f"{fmt(f['recovered_by_split'])} | {fmt(f['honest_overflow'])} | "
                    f"{fmt(f['spurious_split'])} |"
                )
            L.append("")

    L.append("## Released class, end to end (`bi_kappa_distribution<T>`, `no_cap()`)")
    L.append("")
    L.append("| stdlib | precision | κ | non-finite fraction | exceptions | max finite log R |")
    L.append("|---|---|---|---|---|---|")
    for s in res["released"]:
        L.append(
            f"| {s['stdlib']} | {s['precision']} | {s['kappa']:g} | "
            f"{fmt(s['nonfinite_frac'])} | {s['thrown']} | {s['max_finite_log_r']:.4g} |"
        )
    L.append("")

    L.append("## Q4 — supported range")
    L.append("")
    L.append("| stdlib | precision | κ | finite fraction | spurious | honest | status |")
    L.append("|---|---|---|---|---|---|---|")
    for s in res["supported_range"]:
        L.append(
            f"| {s['stdlib']} | {s['precision']} | {s['kappa']:g} | "
            f"{fmt(s['finite_fraction'])} | {fmt(s['spurious_failure_fraction'])} | "
            f"{fmt(s['honest_overflow_fraction'])} | {s['status']} |"
        )
    L.append("")

    L.append("## Q5 — interaction with the component-wise cap")
    L.append("")
    L.append(
        "`non-finite per attempt` counts draws generated inside the rejection loop that "
        "were not finite. `non-finite returned` counts those that reached the caller. "
        "A non-finite draw can never satisfy the box predicate, so the cap **redraws** it."
    )
    L.append("")
    L.append(
        "| stdlib | precision | κ | λ | acceptance | non-finite per attempt | "
        "non-finite returned | exhausted |"
    )
    L.append("|---|---|---|---|---|---|---|---|")
    for s in res["capped"]:
        L.append(
            f"| {s['stdlib']} | {s['precision']} | {s['kappa']:g} | {s['lambda']:g} | "
            f"{fmt(s['acceptance'])} | {fmt(s['nonfinite_redrawn_per_attempt'])} | "
            f"{s['nonfinite_returned']} | {s['exhausted']} |"
        )
    L.append("")

    L.append("## Q3 — log-domain construction, distributional validation")
    L.append("")
    L.append(
        "`W = 1/(1+T) ~ Beta(κ−1/2, 3/2)`, obtained as `expit(−2 log R)` so that `T` is "
        "never formed. Quantiles are compared in `log R`. Threshold α = "
        f"{KS_ALPHA} fixed in advance."
    )
    L.append("")
    L.append("| κ | seed | n | √n·D | KS p | CvM p | max |Δ log R| over p=0.5…0.999 | verdict |")
    L.append("|---|---|---|---|---|---|---|---|")
    for s in res["logdomain_validation"]:
        L.append(
            f"| {s['kappa']:g} | {s['seed']} | {s['n']:,} | {s['ks_sqrt_n_d']:.3f} | "
            f"{s['ks_p']:.3f} | {s['cvm_p']:.3f} | "
            f"{max(s['quantiles']['abs_diff']):.4f} | "
            f"{'PASS' if s['pass'] else 'FAIL'} |"
        )
    L.append("")

    RESULTS.mkdir(exist_ok=True)
    (RESULTS / "exp4_table.md").write_text("\n".join(L) + "\n")


def main() -> None:
    variates = load_jsonl("variates_*.jsonl")
    released = load_jsonl("released_*.jsonl")
    capped = load_jsonl("capped_*.jsonl")
    if not variates:
        sys.exit("no raw/ data — run `make run` first")

    vs = summarize_variates(variates)
    res = {
        "experiment": "exp4_precision",
        "question": (
            "What is the supported numerical range of the released implementation as a "
            "function of kappa, precision and standard library, and how much of the "
            "observed failure is avoidable?"
        ),
        "mode": "uncapped except the Q5 block",
        "environment": environment(),
        "toolchains": sorted(
            {(r["compiler"], r["stdlib"], r["arch"]) for r in variates}
        ),
        "variates": vs,
        "released": summarize_released(released),
        "capped": summarize_capped(capped),
        "supported_range": supported_range(vs),
        "logdomain_validation": validate_logdomain(),
    }
    RESULTS.mkdir(exist_ok=True)
    (RESULTS / "exp4_results.json").write_text(json.dumps(res, indent=2, default=str))
    write_table(res)
    print(f"wrote {RESULTS/'exp4_results.json'} and {RESULTS/'exp4_table.md'}")


if __name__ == "__main__":
    main()
