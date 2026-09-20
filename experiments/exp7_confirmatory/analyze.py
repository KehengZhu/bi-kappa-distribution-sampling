#!/usr/bin/env python3
"""Experiment 7 analysis driver: raw counters and per-draw records in, source data,
family decisions, gate table and verdict out.

Run with:  uv run --project ../../python python analyze.py [--smoke]

Four properties of this file are structural rather than stylistic, and each one exists
because its absence is a named defect in Experiment 6.

* **It refuses to run unless the protocol it judges by is the protocol the run was made
  against.**  ``config/protocol.json`` is hashed here and compared with the hash every probe
  process recorded in ``raw/manifest.csv`` and in every counter row.  A result that cannot be
  tied to the rules it was judged by is not a confirmatory result, so a mismatch is an error
  and not a warning.  The hash is stamped into every file written below.

* **A missing phase is an error that names the phase, not an empty result.**  Experiment 6's
  reader returned ``[]`` for a phase directory that did not exist and the analysis still
  printed a GO/PARTIAL/NO-GO line underneath.  Here every phase loader raises
  :class:`AnalysisError`, the driver prints it and exits non-zero, and no verdict is written.

* **Nothing tracked carries a wall clock, a host name or an absolute path.**  ``make
  reverify`` regenerates every derived artifact and diffs it; a "Generated <utc>" line makes
  that check impossible, and six of Experiment 6's tracked outputs carry one.  The single
  provenance file written here, ``provenance.md``, does record the run's wall clock and build
  identity -- but it copies them out of the recorded raw data rather than reading the clock,
  so it too regenerates byte for byte.

* **Every acceptance threshold comes from the protocol.**  The decision rules live in
  ``exp7_families`` and ``exp7_gates``; this file measures, and hands the measurements over.
  ``exp7_families.Protocol`` raises on a key that was not pre-registered, so a criterion that
  was not frozen cannot be applied from here.

Smoke output is isolated: ``--smoke`` reads ``raw/smoke/`` only and writes ``results/smoke/``
only, and the two modes refuse to read each other's data.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import sys
from collections import defaultdict

import numpy as np
from scipy import special, stats

import exp7_families as F
import exp7_gates as GATES
import exp7_censoring as CEN
import exp7_io as IO
import exp7_portability as PORT
import exp7_stats as S

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))

# The string a cell carries when a quantity does not apply to that row or was not measured.
# A blank cell is never written, so a blank can never be read as a zero.  A quantity that WAS
# computed and came out non-finite is written as `inf`, `-inf` or `nan`, which is a different
# statement from `NA` and is documented as such in results/source_data_README.md.
NA = "NA"

# The environment whose rows carry the single-environment figures and the F1-F5 decisions.
# Chosen by name, deterministically, rather than by whichever file happened to sort first.
PRIMARY_TAG = "libcxx"

# Per-seed resolution of the within-seed half of the cluster bootstrap.  The between-seed
# half -- the component Experiment 6's bootstrap lacked entirely -- is always exact, because
# it resamples the five seeds themselves.  The within-seed half is evaluated on at most this
# many draws per seed, which widens the interval relative to the full sample rather than
# narrowing it; `boot_n_per_seed` records what was used.  The decisive interval for a failure
# probability is the Clopper-Pearson one in the `*_ci_lo/_ci_hi` columns, per PROTOCOL.md 5.2.
BOOT_MAX_PER_SEED = 2000

# PROTOCOL.md 5 requires at least 2000 Monte-Carlo replicates for the F3 null.
F3_MC_REPLICATES = 2000

# Terminal categories, in the enum order src/exp7_common.H assigns them.  CondRecord.cat_* and
# LoaderRecord.status are indices into this list.  Taken from exp7_io rather than restated:
# a second copy of an enum is a second thing to drift out of step with the C++.
CATEGORIES = tuple(IO.CATEGORY_NAMES)
CAT_INDEX = {name: i for i, name in enumerate(CATEGORIES)}
# Categories on which a three-vector did reach the caller.  `finite_but_wrong` and
# `overflow_returned_finite` both return a vector; they are losses of accuracy and of
# representability, not of the return itself, so a tail-retention statement that asked "did a
# value enter the returned sample?" must count them.
RETURNED_FINITE_CATS = (CAT_INDEX["finite"], CAT_INDEX["finite_but_wrong"],
                        CAT_INDEX["overflow_returned_finite"])

QUANTILE_LEVELS = S.QUANTILE_LEVELS
TAIL_LEVELS = S.TAIL_LEVELS

# Fixed decade edges for the P3 conditioning bins, on the intended upper-tail probability of
# the radius.  Frozen here so that the binning is a property of the script and not of the
# data it happens to see.
Q_BIN_EDGES_LOG10 = (0.0, -1.0, -2.0, -3.0, -4.0, -5.0, -6.0, -8.0, -10.0, -15.0,
                     -20.0, -40.0, -100.0, -float("inf"))


class AnalysisError(RuntimeError):
    """A condition under which no verdict may be written."""


# ---------------------------------------------------------------------------
# Raw record readers
#
# exp7_io now reads Experiment 7's files itself: the `EXP7REC` magic, record kind 5, and the
# 32-byte CondRecord that drops Experiment 6's stored `log_w_ref` (which was exactly
# -logaddexp(0, 2 log R), a function of another field).  The layouts, the kind table and the
# category enum are therefore taken from it rather than restated here.  What this wrapper
# adds is the failure mode: a missing or malformed record file must stop the analysis with a
# message that names the file, because every call site below is on the path to a verdict.
# ---------------------------------------------------------------------------
EXP7_KIND_NAME = dict(IO.RECORD_KINDS)


def read_exp7_records(path: str, expect_kind: int, expect_schema: int = 1):
    """``(header, structured array)`` for one Experiment 7 bulk file.

    A truncated file, a wrong record size or a wrong kind is an error here rather than a
    plausible-looking array downstream.
    """
    if not os.path.exists(path):
        raise AnalysisError(f"record file {path} does not exist, but a counter row names it")
    try:
        return IO.read_records(path, expect_kind, expect_schema)
    except IO.SchemaError as exc:
        want = EXP7_KIND_NAME.get(expect_kind, "?")
        raise AnalysisError(f"{exc} (reading it as record kind {expect_kind}, {want})"
                            ) from exc


# ---------------------------------------------------------------------------
# Value helpers
# ---------------------------------------------------------------------------
def jnum(row: dict, key: str, default: float = float("nan")) -> float:
    """A numeric JSON field, honouring results/schema.md 0: non-finite values are the
    strings ``"inf"``, ``"-inf"`` and ``"nan"``, never bare tokens."""
    v = row.get(key)
    if v is None:
        return default
    if isinstance(v, bool):
        return float(v)
    if isinstance(v, str):
        s = v.strip().lower()
        if s == "inf":
            return math.inf
        if s == "-inf":
            return -math.inf
        if s in ("nan", "-nan"):
            return math.nan
        try:
            return float(v)
        except ValueError:
            return default
    return float(v)


def jint(row: dict, key: str, default=None):
    v = row.get(key)
    if v is None:
        return default
    if isinstance(v, bool):
        return int(v)
    try:
        return int(v)
    except (TypeError, ValueError):
        return default


def fmt(v) -> str:
    """One cell.  ``None`` is the missing-value code; a computed non-finite is written as
    such, because "not measured" and "measured, and infinite" are different statements."""
    if v is None:
        return NA
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, (int, np.integer)):
        return str(int(v))
    if isinstance(v, (float, np.floating)):
        x = float(v)
        if math.isnan(x):
            return "nan"
        if math.isinf(x):
            return "inf" if x > 0 else "-inf"
        return repr(x)
    return str(v)


def write_csv(path: str, columns, rows) -> None:
    """Deterministic CSV: fixed column order, ``\\n`` line endings, no locale, no clock."""
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(columns)
        for r in rows:
            unknown = set(r) - set(columns)
            if unknown:
                raise AnalysisError(f"{os.path.basename(path)}: row carries columns that the "
                                    f"header does not declare: {sorted(unknown)}")
            w.writerow([fmt(r.get(c)) for c in columns])


def sha256_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def nanmean(x) -> float:
    """``np.nanmean`` without the all-NaN warning: an empty conditional slice is NaN, which
    the bootstrap already treats as an unresolved resample."""
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    return float(np.mean(x)) if x.size else float("nan")


def stable_seed(*parts) -> int:
    """A reproducible RNG seed from a label.  ``hash()`` is salted per process and would make
    every interval a different one on the next run."""
    label = "|".join(str(p) for p in parts)
    return int.from_bytes(hashlib.sha256(label.encode("utf-8")).digest()[:4], "big")


# ---------------------------------------------------------------------------
# Intervals
# ---------------------------------------------------------------------------
def rate_columns(k: int, n: int, prefix: str, conf: float) -> dict:
    """A rate with its interval and, always, the KIND of interval it is.

    A zero count gets the one-sided upper limit ``1 - (1-conf)^(1/n)``; it never gets a
    zero-width interval.  Experiment 6 emits 28 rows whose interval is ``[0, 0]`` for a zero
    count, and its own README forbids exactly that.
    """
    out = {f"{prefix}_count": int(k), f"{prefix}_n": int(n)}
    if n <= 0:
        out.update({f"{prefix}_rate": None, f"{prefix}_ci_lo": None, f"{prefix}_ci_hi": None,
                    f"{prefix}_interval_kind": "none_no_trials",
                    f"{prefix}_is_upper_bound": None})
        return out
    if k == 0:
        hi = S.zero_count_upper(n, conf)
        out.update({f"{prefix}_rate": 0.0, f"{prefix}_ci_lo": 0.0, f"{prefix}_ci_hi": hi,
                    f"{prefix}_interval_kind": f"one_sided_upper_{conf:g}",
                    f"{prefix}_is_upper_bound": True})
        return out
    lo, hi = S.clopper_pearson(k, n, conf)
    out.update({f"{prefix}_rate": k / n, f"{prefix}_ci_lo": lo, f"{prefix}_ci_hi": hi,
                f"{prefix}_interval_kind": f"two_sided_clopper_pearson_{conf:g}",
                f"{prefix}_is_upper_bound": False})
    return out


def seed_indicator(k: int, n: int, cap: int = BOOT_MAX_PER_SEED) -> np.ndarray:
    """One seed's per-attempt 0/1 outcomes, at a bounded resolution.

    A seed contributes ``n`` Bernoulli outcomes of which ``k`` are ones.  Resampling ``n``
    of them with replacement is exactly ``Binomial(n, k/n)/n``, so the only thing ``cap``
    changes is the resolution of that inner draw: the interval comes out wider, never
    narrower.  ``boot_n_per_seed`` records the value used.
    """
    if n <= 0:
        return np.zeros(0, dtype=float)
    m = int(min(n, cap))
    ones = int(round(k * m / n))
    ones = max(0, min(m, ones))
    return np.concatenate([np.ones(ones, dtype=float), np.zeros(m - ones, dtype=float)])


def cluster_rate_interval(counts_by_seed, sizes_by_seed, proto: F.Protocol,
                          *label) -> dict:
    """Cluster bootstrap over seeds for a rate, with its own recorded RNG seed.

    ``exp7_stats.seed_stratified_bootstrap`` is deliberately not used: its own docstring
    records that its variance has no between-seed component, so no interval built on it can
    see the seed-level heterogeneity that an RNG-stream or portability defect produces.
    """
    resamples = int(proto.get("statistics", "bootstrap", "resamples"))
    conf = float(proto.get("statistics", "bootstrap", "conf"))
    rng_seed = stable_seed(*label)
    vals = {s: seed_indicator(int(counts_by_seed[s]), int(sizes_by_seed[s]))
            for s in sorted(counts_by_seed)}
    vals = {s: v for s, v in vals.items() if v.size}
    r = F.cluster_bootstrap(vals, np.mean, resamples, conf, rng_seed)
    return {"boot_estimate": r["estimate"], "boot_lo": r["lo"], "boot_hi": r["hi"],
            "boot_conf": conf, "boot_resamples": resamples, "boot_rng_seed": rng_seed,
            "boot_kind": "cluster_over_seeds", "boot_n_seeds": r["n_seeds"],
            "boot_n_per_seed": int(min(BOOT_MAX_PER_SEED,
                                       max([sizes_by_seed[s] for s in sizes_by_seed] or [0]))),
            "boot_resolved": r["resolved"]}


def cluster_value_interval(values_by_seed, statistic, proto: F.Protocol, *label) -> dict:
    """Cluster bootstrap over seeds for a statistic of per-seed value arrays."""
    resamples = int(proto.get("statistics", "bootstrap", "resamples"))
    conf = float(proto.get("statistics", "bootstrap", "conf"))
    rng_seed = stable_seed(*label)
    vals = {s: np.asarray(v, dtype=float) for s, v in sorted(values_by_seed.items())}
    vals = {s: v for s, v in vals.items() if v.size}
    r = F.cluster_bootstrap(vals, statistic, resamples, conf, rng_seed)
    return {"boot_estimate": r["estimate"], "boot_lo": r["lo"], "boot_hi": r["hi"],
            "boot_conf": conf, "boot_resamples": resamples, "boot_rng_seed": rng_seed,
            "boot_kind": "cluster_over_seeds", "boot_n_seeds": r["n_seeds"],
            "boot_resolved": r["resolved"]}


# ---------------------------------------------------------------------------
# Statistic -> p-value
#
# The probe reports A^2, D and W^2 rather than the sample they were computed on -- the sample
# is 10^6 values per configuration -- and results/schema.md fixes the estimators so that the
# probe's statistic and this module's agree to the last bit.  What is needed here is
# therefore the null tail of each statistic, at a fully specified null with no estimated
# parameter.
# ---------------------------------------------------------------------------
def ad_pvalue(a2: float) -> float:
    """Upper tail of ``A^2`` under the fully specified null.

    The limiting CDF of Marsaglia & Marsaglia (2004), J. Stat. Softw. 9(2) -- the same
    expression ``exp7_families.anderson_darling_exp1`` evaluates, restated here because that
    function takes a sample and what is available is the statistic.  Verified to agree with
    it to every printed digit on simulated samples.
    """
    a2 = float(a2)
    if not np.isfinite(a2):
        return float("nan")
    if a2 <= 0.0:
        return 1.0
    if a2 < 2.0:
        cdf = (a2 ** -0.5) * np.exp(-1.2337141 / a2) * (
            2.00012 + (0.247105 - (0.0649821 - (0.0347962 - (0.011672 - 0.00168691 * a2)
                                                * a2) * a2) * a2) * a2)
    else:
        cdf = np.exp(-np.exp(1.0776 - (2.30695 - (0.43424 - (0.082433 - (0.008056
                     - 0.0003146 * a2) * a2) * a2) * a2) * a2))
    return float(min(max(1.0 - float(cdf), 0.0), 1.0))


def ks_pvalue(d: float, n: int) -> float:
    """Upper tail of the two-sided Kolmogorov-Smirnov ``D`` at sample size ``n``."""
    if not np.isfinite(d) or n < 2:
        return float("nan")
    return float(stats.kstwo.sf(float(d), int(n)))


def cvm_pvalue(w2: float) -> float:
    """Upper tail of the Cramer-von Mises ``W^2``, asymptotic in ``n``.

    Csorgo & Faraway (1996), the series scipy evaluates for the limiting law.  The finite-n
    correction is negligible at the sample sizes here (n >= 10^5); at n = 2*10^5 the two
    agree to four significant figures on simulated samples.
    """
    w2 = float(w2)
    if not np.isfinite(w2) or w2 <= 0.0:
        return float("nan")
    tot = 0.0
    for k in range(10):
        u = np.exp(special.gammaln(k + 0.5) - special.gammaln(k + 1)) / (
            np.pi ** 1.5 * np.sqrt(w2))
        y = 4.0 * k + 1.0
        q = y * y / (16.0 * w2)
        tot += u * np.sqrt(y) * np.exp(-q) * special.kv(0.25, q)
    return float(min(max(1.0 - float(tot), 0.0), 1.0))


# ---------------------------------------------------------------------------
# Upper-tail exceedance, with the unresolved draws counted
#
# `config/protocol.json -> F5_tail.unresolved_Z_counts_toward_every_threshold` is true, and
# amendment 1.3.0 states the reason: a draw whose `Z` never resolved lies above every `z0`
# by construction, so the count a threshold sees is the resolved exceedances plus the
# unresolved draws, over the cell's total ATTEMPTS.  Testing the resolved subset alone
# against `Binomial(n_resolved, q0)` conditions on resolvability, which is monotone in the
# tail -- the exact bias these statistics exist to detect.  The rule is written for F5's new
# tail members and is applied to F4 as well, because F4 is the same statistic on the scalar
# phase and there is no reading under which one of them should condition and the other not.
#
# `F5_tail.unresolved_Z_note` fixes what "unresolved" means here, and the distinction is not
# a nicety: it is the difference between a sound count and a defective one.
#
#   * The `Z` TRANSFORM underflowed.  `Z = -log I_W(a, 3/2)` is evaluated on a radius the
#     loader did return, and no finite `Z` comes back when the intended tail probability
#     underflows the type.  That draw is further out than any threshold representable at
#     all, so it belongs above every `q0`.
#   * The VELOCITY overflowed.  The loader returned no number, but the draw has a perfectly
#     ordinary `Z`: its radius crossed the type's representable bound, which in the `Z`
#     scale is `z_f = -log f`, with `f` the cell's honest-overflow rate.  It lies above `z0`
#     when `z0 < z_f`, that is when `q0 > f`, and NOT otherwise.
#
# Both arrive at this function the same way -- as a count of draws with no measured `Z` --
# so it is the caller's applicability rule, not this function, that keeps the second from
#
# The excess Kolmogorov-Smirnov statistic has no such choice: an unresolved draw contributes
# no measured excess.  It is computed on the resolved excesses, by the frozen
# `exp7_families.exceedance_test`, and the count it was computed from is published beside it.
#
# `p_count_resolved_only` is carried through as published bookkeeping, never as a decision:
# it is what the count would have said had the unresolved draws been dropped, and the two
# together are what lets a reader see which of the two conventions a cell turns on.
# ---------------------------------------------------------------------------
def exceedance_test_counting_unresolved(z_resolved, n_unresolved: int, n_attempted: int,
                                        q0: float) -> dict:
    """The frozen exceedance test with the unresolved draws restored to its numerator."""
    base = F.exceedance_test(z_resolved, q0)
    if not base.get("resolved"):
        return {"resolved": False}
    k_res, n_res = int(base["observed"]), int(base["n"])
    k = k_res + int(n_unresolved)
    n = int(n_attempted)
    if n < k:
        raise AnalysisError(
            f"exceedance at q0={q0:g}: {k} draws are above z0 ({k_res} resolved, "
            f"{int(n_unresolved)} unresolved) but the cell records only {n} attempts. The "
            "count and the denominator do not come from the same cell.")
    z0 = -math.log(q0)
    zr = np.asarray(z_resolved, dtype=float)
    n_excess = int(np.sum(zr[np.isfinite(zr)] > z0))
    return {"resolved": True, "q0": q0, "n": n, "observed": k, "expected": n * q0,
            "p_count": float(stats.binomtest(k, n, q0).pvalue),
            "p_excess": base["p_excess"],
            "n_resolved": n_res, "observed_resolved": k_res,
            "n_unresolved": int(n_unresolved), "n_excess": n_excess,
            "p_count_resolved_only": base["p_count"], "z0": z0}


def tail_slug(q0: float) -> str:
    """``1em4`` for ``q0 = 1e-4``: a threshold names itself by its decade, in a form that is
    safe as a column name and reversible for a label."""
    e = round(math.log10(q0))
    if q0 > 0 and abs(q0 - 10.0 ** e) <= 1e-12 * q0:
        return f"1e{e}".replace("-", "m")
    return f"{q0:g}".replace("-", "m").replace(".", "p").replace("+", "")


def tail_label(q0: float) -> str:
    return f"1e-{-round(math.log10(q0))}" if q0 > 0 else f"{q0:g}"


# ---------------------------------------------------------------------------
# The analytic honest representability floor (config/honest_floor.md)
# ---------------------------------------------------------------------------
_DIRECTION_GRID = None


def _direction_moment(two_a: float) -> float:
    """``E[(max_j |n_j|)^(2a)]`` for ``n`` uniform on the sphere, by a fixed grid.

    Deterministic by construction: a midpoint rule on ``cos(theta)`` and ``phi``, not a Monte
    Carlo, so the floor is the same number on every run and on every host.
    """
    global _DIRECTION_GRID
    if _DIRECTION_GRID is None:
        nu, nphi = 512, 512
        u = (np.arange(nu) + 0.5) / nu                      # |cos theta| on (0, 1)
        phi = 2.0 * np.pi * (np.arange(nphi) + 0.5) / nphi
        st = np.sqrt(np.maximum(0.0, 1.0 - u[:, None] ** 2))
        m = np.maximum(np.maximum(np.abs(st * np.cos(phi[None, :])),
                                  np.abs(st * np.sin(phi[None, :]))),
                       u[:, None])
        _DIRECTION_GRID = m.ravel()
    m = _DIRECTION_GRID
    if two_a == 0.0:
        return 1.0
    return float(np.mean(np.exp(two_a * np.log(m))))


def honest_floor(kappa: float, precision: str) -> dict:
    """The probability that the intended vector has no representation in the run's type.

    Closed form of config/honest_floor.md, in the log domain throughout so that a floor of
    1e-309 is reported as a number rather than as an underflowed zero:

        log P = lgamma(3/2+a) - lgamma(3/2) - lgamma(a+1) - 2a log MAX
                + a log kappa + log E[(max_j |n_j|)^(2a)]

    in the isotropic, unrotated, theta = 1 configuration P1 and P2 run.  "No failures
    observed" means nothing without this number beside it: at double kappa = 0.55 the floor
    is 1.5e-31, and at float kappa = 0.55 it is 1.4e-4.
    """
    a = float(kappa) - 0.5
    max_val = float(np.finfo(np.float32 if precision == "float" else np.float64).max)
    if a <= 0.0:
        return {"log_floor": float("nan"), "floor": float("nan"), "log10_floor": float("nan")}
    log_floor = (special.gammaln(1.5 + a) - special.gammaln(1.5) - special.gammaln(a + 1.0)
                 - 2.0 * a * math.log(max_val) + a * math.log(float(kappa))
                 + math.log(_direction_moment(2.0 * a)))
    log_floor = float(log_floor)
    return {"log_floor": log_floor, "floor": float(math.exp(log_floor))
            if log_floor > -740.0 else 0.0, "log10_floor": log_floor / math.log(10.0)}


# ---------------------------------------------------------------------------
# Context: which tree is being read, which tree is being written, and the hash gate
# ---------------------------------------------------------------------------
class Context:
    def __init__(self, smoke: bool, proto: F.Protocol):
        self.smoke = smoke
        self.proto = proto
        self.raw_root = os.path.join(HERE, "raw", "smoke") if smoke else os.path.join(HERE,
                                                                                     "raw")
        self.out_dir = os.path.join(HERE, "results", "smoke") if smoke else os.path.join(
            HERE, "results")
        self.manifest_path = os.path.join(self.raw_root, "manifest.csv")
        self.manifest: list[dict] = []
        self.interval_conf = float(proto.get("statistics", "interval_conf"))
        self.notes: list[str] = []
        self._phase_cache: dict = {}

    # -- paths -------------------------------------------------------------
    def rel(self, path: str) -> str:
        return os.path.relpath(path, HERE)

    def out(self, name: str) -> str:
        os.makedirs(self.out_dir, exist_ok=True)
        return os.path.join(self.out_dir, name)

    def resolve_raw(self, recorded: str) -> str:
        """A path a counter row recorded, which is relative to the experiment directory."""
        return os.path.join(HERE, recorded)

    # -- requirement A -----------------------------------------------------
    def check_protocol_hash(self) -> None:
        if not os.path.exists(self.manifest_path):
            raise AnalysisError(
                f"the run manifest {self.rel(self.manifest_path)} does not exist. The "
                f"analysis will not judge a run whose rules it cannot tie to the run: "
                f"execute the phases first (`make {'SMOKE=1 ' if self.smoke else ''}p1 p2 p3 "
                f"p4 p5 p6`).")
        with open(self.manifest_path, encoding="utf-8", newline="") as fh:
            self.manifest = list(csv.DictReader(fh))
        if not self.manifest:
            raise AnalysisError(f"{self.rel(self.manifest_path)} has a header and no rows; no "
                                "output file was recorded, so there is nothing to analyse")
        recorded = sorted({r.get("protocol_sha256", "") for r in self.manifest})
        if recorded != [self.proto.sha256]:
            raise AnalysisError(
                "REFUSING TO RUN: config/protocol.json does not hash to what the run "
                f"recorded.\n  config/protocol.json now:  {self.proto.sha256}\n"
                f"  recorded in {self.rel(self.manifest_path)}: {recorded}\n"
                "A result that cannot be tied to the rules it was judged by is not a "
                "confirmatory result. Either restore the protocol the run was made against, "
                "or rerun the phases against this one.")
        # Smoke and production data are never mixed: they differ by three orders of magnitude
        # in sample size and would sit side by side in a pooled rate.
        smoke_flags = sorted({str(r.get("smoke", "")).lower() for r in self.manifest})
        want = ["true"] if self.smoke else ["false"]
        if smoke_flags != want:
            raise AnalysisError(
                f"{self.rel(self.manifest_path)} records smoke={smoke_flags} but this is a "
                f"{'smoke' if self.smoke else 'production'} analysis, which requires "
                f"smoke={want}. Smoke output never reaches a production result file.")

    # -- loud degradation --------------------------------------------------
    def load_phase(self, phase: str) -> list[dict]:
        """Every counter row a phase produced, and the protocol check they must survive.

        The directory walk is ``exp7_io.load_phase``, which raises
        :class:`exp7_io.MissingPhaseError` for an absent or empty phase rather than
        returning ``[]``; it is re-raised here as an :class:`AnalysisError` carrying the
        command that would produce the phase, because that is the form :func:`main` turns
        into an exit code with no verdict written.
        """
        if phase in self._phase_cache:
            return self._phase_cache[phase]
        try:
            rows = IO.load_phase(os.path.join(HERE, "raw"), phase, smoke=self.smoke)
        except IO.MissingPhaseError as exc:
            raise AnalysisError(
                f"phase {phase}: {exc} No verdict is computed. Run "
                f"`make {'SMOKE=1 ' if self.smoke else ''}{phase}`.") from exc
        bad = sorted({r.get("protocol_sha256") for r in rows} - {self.proto.sha256})
        if bad:
            raise AnalysisError(
                f"phase {phase}: counter rows record protocol hashes {bad}, but "
                f"config/protocol.json hashes to {self.proto.sha256}. The phase was run "
                "against a different protocol than the one this analysis applies.")
        self._phase_cache[phase] = rows
        return rows

    def oracle_rows(self) -> tuple[list[dict], list[dict]]:
        audit = os.path.join(self.raw_root, "oracle_audit.jsonl")
        disag = os.path.join(self.raw_root, "oracle_disagreements.jsonl")
        if not os.path.exists(audit):
            raise AnalysisError(
                f"{self.rel(audit)} does not exist: the arbitrary-precision oracle never "
                f"adjudicated the audit streams. G2 rests on that adjudication and cannot be "
                f"evaluated without it. Run `make {'SMOKE=1 ' if self.smoke else ''}oracle`.")
        if not os.path.exists(disag):
            raise AnalysisError(
                f"{self.rel(disag)} does not exist. PROTOCOL.md 2.3 requires a header record "
                "per audited file whether or not anything disagreed, so that an empty "
                "adjudication is distinguishable from an oracle that never ran.")
        return IO.read_jsonl(audit), IO.read_jsonl(disag)


# ---------------------------------------------------------------------------
# Coverage checks: the matrix that was frozen against the matrix that ran
# ---------------------------------------------------------------------------
def check_ladder_coverage(ctx: Context, phase: str, rows: list[dict], methods) -> None:
    ladder = [float(k) for k in ctx.proto.get("matrix", "kappa_ladder")]
    precisions = list(ctx.proto.get("matrix", "precisions"))
    seeds = sorted({jint(r, "seed") for r in rows})
    if not ctx.smoke:
        want = sorted(int(s) for s in ctx.proto.get("seeds", "production"))
        if seeds != want:
            raise AnalysisError(
                f"phase {phase}: seeds {seeds} ran, but PROTOCOL.md 3 freezes {want}. No "
                "result may be recomputed on a different seed block.")
    have = {(r.get("tag"), r.get("method"), r.get("precision"),
             round(jnum(r, "kappa"), 12), jint(r, "seed")) for r in rows}
    tags = sorted({r.get("tag") for r in rows})
    missing = []
    for tag in tags:
        for m in methods:
            for p in precisions:
                for k in ladder:
                    for s in seeds:
                        if (tag, m, p, round(k, 12), s) not in have:
                            missing.append(f"{tag}/{m}/{p}/kappa={k:g}/seed={s}")
    if missing:
        raise AnalysisError(
            f"phase {phase}: {len(missing)} configuration(s) of the frozen matrix produced "
            f"no row, for example {missing[:5]}. A partial matrix is not analysed into a "
            "verdict.")


# ---------------------------------------------------------------------------
# P1 -- scalar validation
# ---------------------------------------------------------------------------
P1_QUANTILE_FIELDS = ("lo_index", "hi_index", "achieved_coverage", "lo_value", "hi_value",
                      "point_value", "target", "error", "width", "resolved", "informative",
                      "covered")
P1_TAIL_FIELDS = ("z0", "observed_resolved", "observed_nonresolved", "observed_total",
                  "expected_resolved", "expected_per_attempt", "n_test", "n_excess",
                  "p_count", "p_count_resolved_only", "p_excess")


def qkey(p: float) -> str:
    return f"q{p:g}".replace(".", "p")


def q0key(q0: float) -> str:
    return f"tail{q0:g}".replace("-", "m").replace(".", "p").replace("+", "")


def analyse_p1(ctx: Context) -> dict:
    proto = ctx.proto
    rows = ctx.load_phase("p1")
    native = [r for r in rows if r.get("layer") == "native"]
    klass = [r for r in rows if r.get("layer") == "class"]
    if not native:
        raise AnalysisError("phase p1: no `native` counter rows; the verified replica layer "
                            "is where F1-F4 read their sample from")
    check_ladder_coverage(ctx, "p1", native, ("LEGACY", "CANDIDATE"))

    q0s = [float(q) for q in proto.get("matrix", "tail_q0")]
    alpha_F1 = proto.alpha("F1_radial_law")

    columns = ["phase", "tag", "stdlib", "arch", "execution", "layer", "method", "precision",
               "kappa", "shape_a", "seed", "n_attempted"]
    columns += [f"cat_{c}" for c in CATEGORIES]
    columns += ["n_finite", "n_avoidable", "n_honest", "nonfinite_output",
                "x2_zero", "x2_subnormal",
                "gamma_variates", "uniform_variates", "engine_calls", "seconds",
                "seconds_per_attempt", "max_finite_log_r",
                "n_resolved", "n_nonresolved", "loss_fraction", "conditional",
                "honest_floor_rate", "honest_floor_log10", "observed_over_floor",
                "ad_statistic", "ad_pvalue", "ks_statistic", "ks_pvalue",
                "cvm_statistic", "cvm_pvalue", "f1_cell_simes_p",
                "tail_records", "tail_file", "class_n_finite", "class_replica_agrees",
                "protocol_sha256"]
    for pref in ("failure", "avoidable", "honest"):
        columns += [f"{pref}_count", f"{pref}_n", f"{pref}_rate", f"{pref}_ci_lo",
                    f"{pref}_ci_hi", f"{pref}_interval_kind", f"{pref}_is_upper_bound"]
    for p in QUANTILE_LEVELS:
        columns += [f"{qkey(p)}_{f}" for f in P1_QUANTILE_FIELDS]
    for q0 in q0s:
        columns += [f"{q0key(q0)}_{f}" for f in P1_TAIL_FIELDS]

    class_by_key = {(r.get("tag"), r.get("method"), r.get("precision"),
                     round(jnum(r, "kappa"), 12), jint(r, "seed")): r for r in klass}

    out_rows: list[dict] = []
    ecdf_acc: dict = defaultdict(lambda: {"counts": None, "n": 0, "seeds": 0})
    f1_cells: dict = defaultdict(list)
    f2_cells: list[dict] = []
    f3_obs: dict = defaultdict(dict)
    f4_cells: dict = defaultdict(list)
    failure_by_config: dict = defaultdict(dict)
    size_by_config: dict = defaultdict(dict)
    benign = {}
    benign_digests: dict = {}
    candidate_avoidable_p1 = 0
    tail_mismatches: list[str] = []

    for r in sorted(native, key=lambda r: (r.get("tag"), r.get("method"), r.get("precision"),
                                           jnum(r, "kappa"), jint(r, "seed"))):
        tag, method = r.get("tag"), r.get("method")
        precision = r.get("precision")
        kappa, a = jnum(r, "kappa"), jnum(r, "shape_a")
        seed = jint(r, "seed")
        n = jint(r, "n_attempted", 0)
        n_resolved = jint(r, "n_resolved", 0)
        n_nonres = jint(r, "n_nonresolved", 0)
        floor = honest_floor(kappa, precision)
        nonfinite = jint(r, "nonfinite_output", 0)

        row = {"phase": "p1", "tag": tag, "stdlib": r.get("stdlib"), "arch": r.get("arch"),
               "execution": r.get("execution"), "layer": "native", "method": method,
               "precision": precision, "kappa": kappa, "shape_a": a, "seed": seed,
               "n_attempted": n, "protocol_sha256": proto.sha256}
        for c in CATEGORIES:
            row[f"cat_{c}"] = jint(r, f"cat_{c}", 0)
        for f in ("n_finite", "n_avoidable", "n_honest", "nonfinite_output", "x2_zero",
                  "x2_subnormal", "gamma_variates", "uniform_variates", "engine_calls",
                  "tail_records"):
            row[f] = jint(r, f, 0)
        row["seconds"] = jnum(r, "seconds")
        row["seconds_per_attempt"] = jnum(r, "seconds") / n if n else None
        row["max_finite_log_r"] = jnum(r, "max_finite_log_r")
        row["n_resolved"] = n_resolved
        row["n_nonresolved"] = n_nonres
        row["loss_fraction"] = jnum(r, "loss_fraction")
        row["conditional"] = bool(r.get("conditional"))
        row["honest_floor_rate"] = floor["floor"]
        row["honest_floor_log10"] = floor["log10_floor"]
        row["observed_over_floor"] = ((nonfinite / n) / floor["floor"]
                                      if n and floor["floor"] > 0 else None)
        row["tail_file"] = r.get("tail_file")

        row.update(rate_columns(nonfinite, n, "failure", ctx.interval_conf))
        row.update(rate_columns(jint(r, "n_avoidable", 0), n, "avoidable", ctx.interval_conf))
        row.update(rate_columns(jint(r, "n_honest", 0), n, "honest", ctx.interval_conf))

        if method == "CANDIDATE":
            candidate_avoidable_p1 += jint(r, "n_avoidable", 0)

        # --- F1: goodness of fit of Z against Exp(1) -------------------------------
        a2, d, w2 = (jnum(r, "ad_statistic"), jnum(r, "ks_statistic"),
                     jnum(r, "cvm_statistic"))
        p_ad, p_ks, p_cvm = ad_pvalue(a2), ks_pvalue(d, n_resolved), cvm_pvalue(w2)
        row.update({"ad_statistic": a2, "ad_pvalue": p_ad, "ks_statistic": d,
                    "ks_pvalue": p_ks, "cvm_statistic": w2, "cvm_pvalue": p_cvm})
        # corrections.F1_contents: F1 is Anderson-Darling, KS and Cramer-von Mises on Z and
        # nothing else.  The Beta test on W that PROTOCOL.md 5's table lists alongside them is
        # withdrawn -- the probe emits no W sample, Z is a monotone transform of log W, and G0
        # already checks the log-q evaluator against an arbitrary-precision incomplete beta to
        # 1e-10, which is a stronger check on the transform than comparing the two routes.
        cell_p = F.simes_global([p_ad, p_ks, p_cvm])
        row["f1_cell_simes_p"] = cell_p
        if tag == PRIMARY_TAG and method == "CANDIDATE":
            f1_cells[(precision, kappa)].extend([p_ad, p_ks, p_cvm])

        # --- F2 / F3: the five order-statistic brackets ----------------------------
        for qd in r.get("quantiles", []):
            p = float(qd["p"])
            lo_v, hi_v = jnum(qd, "lo_value"), jnum(qd, "hi_value")
            point = jnum(qd, "point_value")
            target = S.exact_log_r_quantile(p, a)
            resolved = bool(np.isfinite(lo_v) and np.isfinite(hi_v))
            inf_cell = F.informative(lo_v, hi_v, proto) if resolved else False
            covered = bool(resolved and lo_v <= target <= hi_v)
            k = qkey(p)
            row.update({f"{k}_lo_index": jint(qd, "lo_index"),
                        f"{k}_hi_index": jint(qd, "hi_index"),
                        f"{k}_achieved_coverage": jnum(qd, "achieved_coverage"),
                        f"{k}_lo_value": lo_v, f"{k}_hi_value": hi_v,
                        f"{k}_point_value": point, f"{k}_target": target,
                        f"{k}_error": point - target if np.isfinite(point) else point,
                        f"{k}_width": hi_v - lo_v, f"{k}_resolved": resolved,
                        f"{k}_informative": inf_cell,
                        f"{k}_covered": covered if resolved else None})
            if tag == PRIMARY_TAG and method == "CANDIDATE":
                f2_cells.append({"precision": precision, "kappa": kappa, "seed": seed,
                                 "p": p, "resolved": resolved, "informative": inf_cell,
                                 "covered": covered})
                f3_obs[(precision, kappa, p)][seed] = (point, target, a, n)

        # --- F4: exceedance counts, and the excesses from the tail file ------------
        tail_objs = {float(t["q0"]): t for t in r.get("tail", [])}
        z_tail = None
        tail_path = r.get("tail_file")
        if tail_path:
            _h, arr = read_exp7_records(ctx.resolve_raw(tail_path), 5)
            z_tail = np.asarray(arr["z"], dtype=float)
            lowest = tail_objs.get(max(q0s))
            if lowest is not None and z_tail.size != int(lowest["observed_resolved"]):
                tail_mismatches.append(
                    f"{tag}/{method}/{precision}/kappa={kappa:g}/seed={seed}: tail file holds "
                    f"{z_tail.size} records, the summary row says "
                    f"{int(lowest['observed_resolved'])} resolved exceedances of z0="
                    f"{float(lowest['z0']):.4f}")
        z_full = None
        if z_tail is not None and n_resolved >= z_tail.size:
            # The tail file holds every resolved Z above the lowest protocol threshold.  The
            # sub-threshold draws affect neither an exceedance count nor an excess, so the
            # exceedance test is evaluated on the tail file padded back to the resolved
            # sample size with a sub-threshold sentinel; k, n and both p-values are exact.
            z_full = np.concatenate([z_tail, np.zeros(n_resolved - z_tail.size)])
        tests = []
        for q0 in q0s:
            t = tail_objs.get(q0, {})
            kk = q0key(q0)
            obs_res = jint(t, "observed_resolved", 0)
            obs_non = jint(t, "observed_nonresolved", 0)
            if obs_non != n_nonres:
                raise AnalysisError(
                    f"{tag}/{method}/{precision}/kappa={kappa:g}/seed={seed}: the q0={q0:g} "
                    f"tail object records {obs_non} unresolved draws while the row records "
                    f"{n_nonres}. An unresolved draw is above every threshold, so the two "
                    "must be the same number; they are not, so the exceedance count has no "
                    "unambiguous numerator.")
            got = {}
            if z_full is not None:
                # protocol.json F5_tail.unresolved_Z_counts_toward_every_threshold.
                got = exceedance_test_counting_unresolved(z_full, obs_non, n, q0)
                tests.append(got)
                if got["observed_resolved"] != obs_res:
                    tail_mismatches.append(
                        f"{tag}/{method}/{precision}/kappa={kappa:g}/seed={seed}: the tail "
                        f"file holds {got['observed_resolved']} resolved draws above "
                        f"z0={got['z0']:.4f} (q0={q0:g}) but the summary row says "
                        f"{obs_res}")
            row.update({f"{kk}_z0": jnum(t, "z0"),
                        f"{kk}_observed_resolved": obs_res,
                        f"{kk}_observed_nonresolved": obs_non,
                        f"{kk}_observed_total": obs_res + obs_non,
                        f"{kk}_expected_resolved": n_resolved * q0,
                        f"{kk}_expected_per_attempt": n * q0,
                        f"{kk}_n_test": got.get("n"),
                        f"{kk}_n_excess": got.get("n_excess"),
                        f"{kk}_p_count": got.get("p_count"),
                        f"{kk}_p_count_resolved_only": got.get("p_count_resolved_only"),
                        f"{kk}_p_excess": got.get("p_excess")})
        if tag == PRIMARY_TAG and method == "CANDIDATE" and tests:
            f4_cells[(precision, kappa)].append(tests)

        # --- the released class beside its instrumented replica --------------------
        ck = (tag, method, precision, round(kappa, 12), seed)
        cr = class_by_key.get(ck)
        row["class_n_finite"] = jint(cr, "n_finite") if cr else None
        row["class_replica_agrees"] = (bool(jint(cr, "n_finite") == jint(r, "n_finite")
                                            and jint(cr, "nonfinite_output")
                                            == jint(r, "nonfinite_output"))
                                       if cr else None)

        # --- the benign control (PROTOCOL.md 2.1) ----------------------------------
        # The rejected primitive's failure mode was that at its own benign control every
        # draw fell through to std::gamma_distribution, so the control exercised the
        # standard library and not the candidate.  Three computed facts rule that out here:
        # exactly two Gamma variates per attempt, which is X1 and the *boosted* Gamma of the
        # shape-boosting identity and not a direct Gamma(a) draw; at least one further
        # uniform per Gamma, which is the boosting uniform; and every attempt returning.
        # The fourth fact, cross-library digest equality at this configuration, is collected
        # below: a standard-library sampler in the path could not produce it.
        if method == "CANDIDATE" and precision == "double" and abs(kappa - 2.0) < 1e-12:
            g, u = jint(r, "gamma_variates", 0), jint(r, "uniform_variates", 0)
            ok = bool(n > 0 and g == 2 * n and u > g and jint(r, "n_finite", 0) == n)
            benign.setdefault(tag, []).append(ok)
            if cr is not None:
                benign_digests.setdefault((precision, seed), {})[tag] = cr.get(
                    "digest_sha256")

        # --- per-configuration seed bookkeeping ------------------------------------
        cfg = (tag, method, precision, round(kappa, 12))
        failure_by_config[cfg][seed] = nonfinite
        size_by_config[cfg][seed] = n

        # --- the Z ECDF, pooled over seeds, for the figure -------------------------
        counts = np.asarray(r.get("z_ecdf_counts", []), dtype=float)
        acc = ecdf_acc[(tag, method, precision, round(kappa, 12))]
        if counts.size:
            acc["counts"] = counts if acc["counts"] is None else acc["counts"] + counts
            acc["n"] += n_resolved
            acc["seeds"] += 1
            acc["step"] = jnum(r, "z_grid_step", 0.05)

        out_rows.append(row)

    if tail_mismatches:
        raise AnalysisError(
            "the P1 tail files and the P1 summary rows disagree on the exceedance count, so "
            "F4 cannot be evaluated on them:\n  " + "\n  ".join(tail_mismatches[:5]))

    write_csv(ctx.out("scalar_validation.csv"), columns, out_rows)

    # --- scalar_ecdf.csv ---------------------------------------------------------
    ecdf_cols = ["tag", "method", "precision", "kappa", "n_seeds", "n_resolved_pooled",
                 "grid_index", "z", "null_cdf", "ecdf", "ecdf_residual", "band_lo",
                 "band_hi", "band_alpha", "protocol_sha256"]
    ecdf_rows = []
    for (tag, method, precision, kappa), acc in sorted(ecdf_acc.items(), key=lambda kv: str(
            kv[0])):
        if acc["counts"] is None or acc["n"] == 0:
            continue
        step = float(acc.get("step", 0.05))
        band = math.sqrt(-0.5 * math.log(alpha_F1 / 2.0) / acc["n"])
        for i, c in enumerate(acc["counts"]):
            z = step * i
            null = -math.expm1(-z)
            e = float(c) / acc["n"]
            ecdf_rows.append({"tag": tag, "method": method, "precision": precision,
                              "kappa": kappa, "n_seeds": acc["seeds"],
                              "n_resolved_pooled": acc["n"], "grid_index": i, "z": z,
                              "null_cdf": null, "ecdf": e, "ecdf_residual": e - null,
                              "band_lo": -band, "band_hi": band, "band_alpha": alpha_F1,
                              "protocol_sha256": proto.sha256})
    write_csv(ctx.out("scalar_ecdf.csv"), ecdf_cols, ecdf_rows)

    # --- honest_floor.csv --------------------------------------------------------
    floor_cols = ["precision", "kappa", "shape_a", "honest_floor_rate", "honest_floor_log10",
                  "type_max", "log_type_max", "configuration", "protocol_sha256"]
    floor_rows = []
    for precision in sorted(proto.get("matrix", "precisions")):
        for kappa in proto.get("matrix", "kappa_ladder"):
            fl = honest_floor(float(kappa), precision)
            mx = float(np.finfo(np.float32 if precision == "float" else np.float64).max)
            floor_rows.append({
                "precision": precision, "kappa": float(kappa),
                "shape_a": float(kappa) - 0.5, "honest_floor_rate": fl["floor"],
                "honest_floor_log10": fl["log10_floor"], "type_max": mx,
                "log_type_max": math.log(mx),
                "configuration": "isotropic, unrotated, theta_perp = theta_par = 1",
                "protocol_sha256": proto.sha256})
    write_csv(ctx.out("honest_floor.csv"), floor_cols, floor_rows)

    # --- families ----------------------------------------------------------------
    fam = {}
    fam["F1_radial_law"] = F.family_F1(proto, [
        {"label": f"{prec} kappa={k:g}", "pvalues": pv}
        for (prec, k), pv in sorted(f1_cells.items(), key=lambda kv: (kv[0][0], kv[0][1]))])

    # corrections.F2_informativeness: the miss count runs over every RESOLVED interval,
    # which is the set the frozen Poisson-binomial null was built over.  The informativeness
    # flag of PROTOCOL.md 5.1 is published as a map below and changes no count -- dropping
    # the non-informative cells from the count would leave F2 with no null to be tested
    # against.
    resolved_cells = [c for c in f2_cells if c["resolved"]]
    misses = sum(1 for c in resolved_cells if not c["covered"])
    fam["F2_quantile_coverage"] = F.family_F2(proto, misses, len(resolved_cells))

    f3_cells = []
    for (precision, kappa, p), by_seed in sorted(f3_obs.items(),
                                                 key=lambda kv: (kv[0][0], kv[0][1],
                                                                 kv[0][2])):
        zs = []
        for seed in sorted(by_seed):
            point, target, a, n = by_seed[seed]
            zs.append(f3_standardized(point, target, a, n, p,
                                      stable_seed("F3", precision, f"{kappa:.12g}",
                                                  f"{p:g}", seed)))
        f3_cells.append({"label": f"{precision} kappa={kappa:g} p={p:g}", "z": zs})
    fam["F3_quantile_direction"] = F.family_F3(proto, f3_cells)

    fam["F4_upper_tail_mass"] = F.family_F4(proto, [
        {"label": f"{prec} kappa={k:g}", "tests": [t for per_seed in lst for t in per_seed]}
        for (prec, k), lst in sorted(f4_cells.items(), key=lambda kv: (kv[0][0], kv[0][1]))])

    informativeness = [{"precision": c["precision"], "kappa": c["kappa"], "p": c["p"],
                        "seed": c["seed"], "resolved": c["resolved"],
                        "informative": c["informative"]} for c in f2_cells]

    benign_digest_agrees = bool(benign_digests) and all(
        len(set(d.values())) == 1 and all(v for v in d.values())
        for d in benign_digests.values())

    return {"families": fam, "scalar_rows": out_rows, "f2_cells": f2_cells,
            "informativeness": informativeness, "benign": benign,
            "benign_digest_agrees": benign_digest_agrees,
            "candidate_avoidable_p1": candidate_avoidable_p1,
            "failure_by_config": failure_by_config, "size_by_config": size_by_config,
            "class_rows": klass, "native_rows": native,
            "f2_misses": misses, "f2_resolved": len(resolved_cells)}


def f3_standardized(point: float, target: float, a: float, n: int, p: float,
                    rng_seed: int) -> float:
    """The signed quantile error, standardized against a parametric Monte-Carlo null.

    The empirical p-quantile of a sample of size ``n`` is the order statistic of rank
    ``k = ceil(pn)``, and under the exact law ``F(X_(k)) ~ Beta(k, n-k+1)``.  Drawing that
    Beta and pushing it through the exact quantile function is therefore a draw from the null
    distribution of the estimator at the same ``n`` and the same ``a`` -- the parametric Monte
    Carlo PROTOCOL.md 5 requires, evaluated exactly rather than by simulating 10^6 variates
    2000 times over.  Standardizing rather than using the raw error is what makes the null
    usable at small shape, where it is strongly asymmetric.
    """
    if not (np.isfinite(point) and np.isfinite(target) and n > 0 and a > 0):
        return float("nan")
    k = int(math.ceil(p * n))
    k = max(1, min(n, k))
    rng = np.random.default_rng(rng_seed)
    u = rng.beta(k, n - k + 1, size=F3_MC_REPLICATES)
    u = np.clip(u, 1e-15, 1.0 - 1e-15)
    null = np.array([S.exact_log_r_quantile(float(x), a) for x in u], dtype=float)
    null = null[np.isfinite(null)]
    if null.size < F3_MC_REPLICATES // 2:
        return float("nan")
    sd = float(np.std(null, ddof=1))
    if not np.isfinite(sd) or sd <= 0.0:
        return float("nan")
    return float((point - float(np.mean(null))) / sd)


# ---------------------------------------------------------------------------
# P2 -- the mechanism decomposition and the failure envelope
# ---------------------------------------------------------------------------
AUDIT_STRATA = ("within_margin_of_a_type_limit", "method_disagreement",
                "finite_but_wrong_candidate", "avoidable_loss_candidate",
                "subnormal_or_zero_denominator_representable_target",
                "unambiguous_failure_beyond_margin", "uniform_sample")


def analyse_p2(ctx: Context, p1: dict) -> dict:
    proto = ctx.proto
    rows = ctx.load_phase("p2")
    paired = [r for r in rows if r.get("layer") == "paired"]
    native = [r for r in rows if r.get("layer") == "native"]
    if not paired:
        raise AnalysisError("phase p2: no `paired` counter rows; the paired layer is the "
                            "authoritative mechanism decomposition and G2 rests on it")
    check_ladder_coverage(ctx, "p2", paired, ("LEGACY", "CANDIDATE", "QF"))

    columns = ["phase", "tag", "stdlib", "arch", "execution", "layer", "scope", "method",
               "diagnostic_only", "precision", "kappa", "shape_a", "seed", "n_attempted"]
    columns += [f"cat_{c}" for c in CATEGORIES]
    columns += ["n_finite", "n_avoidable", "n_honest", "nonfinite_output",
                "candidate_avoidable_count", "accounting_ok", "accounting_residual",
                "ref_nonrepresentable", "honest_minus_ref", "rotation_recoverable",
                "near_limit", "x2_zero", "x2_subnormal",
                "n_rel_err", "mean_rel_err_radius", "max_rel_err_radius",
                "max_rel_error_threshold", "max_log_r_ref", "max_finite_log_component",
                "honest_floor_rate", "honest_floor_log10", "observed_over_floor",
                "audit_margin_log_units", "audit_margin_assertions",
                "audit_margin_assertion_failures", "audit_min_unambiguous_margin",
                "audited", "seed_homogeneity_chi2", "seed_homogeneity_p",
                "cross_phase_native_agrees", "protocol_sha256"]
    for pref in ("failure", "avoidable", "honest"):
        columns += [f"{pref}_count", f"{pref}_n", f"{pref}_rate", f"{pref}_ci_lo",
                    f"{pref}_ci_hi", f"{pref}_interval_kind", f"{pref}_is_upper_bound"]
    columns += ["boot_estimate", "boot_lo", "boot_hi", "boot_conf", "boot_resamples",
                "boot_rng_seed", "boot_kind", "boot_n_seeds", "boot_n_per_seed",
                "boot_resolved"]
    for st in AUDIT_STRATA:
        columns += [f"audit_{st}_total", f"audit_{st}_audited", f"audit_{st}_rate"]

    p1_native = {(r.get("tag"), r.get("method"), r.get("precision"),
                  round(jnum(r, "kappa"), 12), jint(r, "seed")): r
                 for r in p1["native_rows"]}

    out_rows: list[dict] = []
    accounting_failures = 0
    candidate_avoidable_p2 = 0
    honest_minus_ref_total = 0
    stratum_totals = {st: [0, 0] for st in AUDIT_STRATA}
    per_cfg: dict = defaultdict(dict)
    cross_phase_disagreements = 0

    for r in sorted(rows, key=lambda r: (r.get("layer"), r.get("tag"), r.get("method"),
                                         r.get("precision"), jnum(r, "kappa"),
                                         jint(r, "seed"))):
        layer = r.get("layer")
        tag, method, precision = r.get("tag"), r.get("method"), r.get("precision")
        kappa, a, seed = jnum(r, "kappa"), jnum(r, "shape_a"), jint(r, "seed")
        n = jint(r, "n_attempted", 0)
        n_fin, n_av, n_hon = (jint(r, "n_finite", 0), jint(r, "n_avoidable", 0),
                              jint(r, "n_honest", 0))
        residual = n - (n_fin + n_av + n_hon)
        ok = bool(r.get("accounting_ok")) and residual == 0
        if not ok:
            accounting_failures += 1
        floor = honest_floor(kappa, precision)
        nonfinite = jint(r, "nonfinite_output", 0)
        ref_nonrep = jint(r, "ref_nonrepresentable")

        row = {"phase": "p2", "tag": tag, "stdlib": r.get("stdlib"), "arch": r.get("arch"),
               "execution": r.get("execution"), "layer": layer, "scope": "seed",
               "method": method, "diagnostic_only": bool(r.get("diagnostic_only")),
               "precision": precision, "kappa": kappa, "shape_a": a, "seed": seed,
               "n_attempted": n, "accounting_ok": ok, "accounting_residual": residual,
               "protocol_sha256": proto.sha256}
        for c in CATEGORIES:
            row[f"cat_{c}"] = jint(r, f"cat_{c}", 0)
        row.update({"n_finite": n_fin, "n_avoidable": n_av, "n_honest": n_hon,
                    "nonfinite_output": nonfinite,
                    "candidate_avoidable_count": n_av if method == "CANDIDATE" else None,
                    "ref_nonrepresentable": ref_nonrep,
                    "honest_minus_ref": (n_hon - ref_nonrep) if ref_nonrep is not None
                    else None,
                    "rotation_recoverable": jint(r, "rotation_recoverable"),
                    "near_limit": jint(r, "near_limit"),
                    "x2_zero": jint(r, "x2_zero"), "x2_subnormal": jint(r, "x2_subnormal"),
                    "n_rel_err": jint(r, "n_rel_err"),
                    "mean_rel_err_radius": jnum(r, "mean_rel_err_radius"),
                    "max_rel_err_radius": jnum(r, "max_rel_err_radius"),
                    "max_rel_error_threshold": jnum(r, "max_rel_error_threshold"),
                    "max_log_r_ref": jnum(r, "max_log_r_ref"),
                    "max_finite_log_component": jnum(r, "max_finite_log_component"),
                    "honest_floor_rate": floor["floor"],
                    "honest_floor_log10": floor["log10_floor"],
                    "observed_over_floor": ((nonfinite / n) / floor["floor"]
                                            if n and floor["floor"] > 0 else None),
                    "audit_margin_log_units": jnum(r, "audit_margin_log_units"),
                    "audit_margin_assertions": jint(r, "audit_margin_assertions"),
                    "audit_margin_assertion_failures":
                        jint(r, "audit_margin_assertion_failures"),
                    "audit_min_unambiguous_margin": jnum(r, "audit_min_unambiguous_margin"),
                    "audited": jint(r, "audited")})
        row.update(rate_columns(nonfinite, n, "failure", ctx.interval_conf))
        row.update(rate_columns(n_av, n, "avoidable", ctx.interval_conf))
        row.update(rate_columns(n_hon, n, "honest", ctx.interval_conf))
        for st in AUDIT_STRATA:
            t, au = jint(r, f"audit_{st}_total"), jint(r, f"audit_{st}_audited")
            row[f"audit_{st}_total"] = t
            row[f"audit_{st}_audited"] = au
            row[f"audit_{st}_rate"] = jnum(r, f"audit_{st}_rate")
            if layer == "paired" and method == "CANDIDATE" and t is not None:
                stratum_totals[st][0] += t
                stratum_totals[st][1] += au or 0

        if layer == "native":
            k1 = p1_native.get((tag, method, precision, round(kappa, 12), seed))
            agree = None
            if k1 is not None:
                agree = bool(jint(k1, "n_finite") == n_fin
                             and jint(k1, "nonfinite_output") == nonfinite
                             and jint(k1, "engine_calls") == jint(r, "engine_calls"))
                if not agree:
                    cross_phase_disagreements += 1
            row["cross_phase_native_agrees"] = agree

        if layer == "paired" and method == "CANDIDATE":
            candidate_avoidable_p2 += n_av
            if ref_nonrep is not None:
                honest_minus_ref_total += abs(n_hon - ref_nonrep)

        per_cfg[(layer, tag, method, precision, round(kappa, 12))][seed] = (nonfinite, n)
        out_rows.append(row)

    # --- pooled rows, with the seed-homogeneity diagnostic and a cluster bootstrap ---
    for key in sorted(per_cfg, key=str):
        layer, tag, method, precision, kappa = key
        by_seed = per_cfg[key]
        counts = {s: by_seed[s][0] for s in by_seed}
        sizes = {s: by_seed[s][1] for s in by_seed}
        n = sum(sizes.values())
        k = sum(counts.values())
        floor = honest_floor(kappa, precision)
        chi2, p_hom = F.seed_homogeneity([counts[s] for s in sorted(counts)],
                                         [sizes[s] for s in sorted(sizes)])
        row = {"phase": "p2", "tag": tag, "stdlib": None, "arch": None, "execution": None,
               "layer": layer, "scope": "pooled", "method": method,
               "diagnostic_only": method == "QF", "precision": precision, "kappa": kappa,
               "shape_a": kappa - 0.5, "seed": None, "n_attempted": n,
               "n_finite": None, "n_avoidable": None, "n_honest": None,
               "nonfinite_output": k, "accounting_ok": None, "accounting_residual": None,
               "honest_floor_rate": floor["floor"], "honest_floor_log10": floor["log10_floor"],
               "observed_over_floor": ((k / n) / floor["floor"]
                                       if n and floor["floor"] > 0 else None),
               "seed_homogeneity_chi2": chi2, "seed_homogeneity_p": p_hom,
               "protocol_sha256": proto.sha256}
        row.update(rate_columns(k, n, "failure", ctx.interval_conf))
        row.update(cluster_rate_interval(counts, sizes, proto, "p2", layer, tag, method,
                                         precision, f"{kappa:.12g}"))
        out_rows.append(row)

    write_csv(ctx.out("failure_envelope.csv"), columns, out_rows)

    coverage = {}
    for st in AUDIT_STRATA:
        total, audited = stratum_totals[st]
        coverage[st] = (audited / total) if total > 0 else float("nan")
    return {"rows": out_rows, "accounting_failures": accounting_failures,
            "candidate_avoidable_p2": candidate_avoidable_p2,
            "honest_minus_ref_total": honest_minus_ref_total,
            "stratum_totals": stratum_totals, "audit_coverage": coverage,
            "cross_phase_disagreements": cross_phase_disagreements,
            "paired_rows": paired, "native_rows": native}


# ---------------------------------------------------------------------------
# P3 -- conditioning and the tail consequence
# ---------------------------------------------------------------------------
def analyse_p3(ctx: Context, p2: dict) -> dict:
    proto = ctx.proto
    rows = ctx.load_phase("p3")
    paired = [r for r in rows if r.get("layer") == "paired"]
    if not paired:
        raise AnalysisError("phase p3: no `paired` counter rows")

    bin_cols = ["tag", "method", "precision", "kappa", "shape_a", "scope", "seed",
                "bin_kind", "bin_index", "bin_lo", "bin_hi", "n_in_bin", "n_success",
                "success_definition", "success_rate", "ci_lo", "ci_hi", "interval_kind",
                "is_upper_bound", "protocol_sha256"]
    tail_cols = ["tag", "method", "precision", "kappa", "shape_a", "scope", "p",
                 "log_r_target", "n_attempts", "n_success", "success_definition",
                 "n_above_target_all_attempts", "n_above_target_and_returned",
                 "failure_atom_mass_per_attempt",
                 "per_attempt_law_has_failure_atom", "per_attempt_quantile_reported",
                 "per_attempt_quantile_reason",
                 "tail_mass_per_attempt", "tail_mass_per_attempt_lo",
                 "tail_mass_per_attempt_hi",
                 "returned_tail_ratio_per_attempt", "returned_tail_ratio_per_attempt_lo",
                 "returned_tail_ratio_per_attempt_hi",
                 "tail_mass_conditional_on_success",
                 "tail_mass_conditional_on_success_lo", "tail_mass_conditional_on_success_hi",
                 "returned_tail_ratio_conditional_on_success",
                 "returned_tail_ratio_conditional_on_success_lo",
                 "returned_tail_ratio_conditional_on_success_hi",
                 "boot_conf", "boot_resamples", "boot_kind", "boot_n_seeds",
                 "boot_n_per_seed", "protocol_sha256"]

    bin_rows: list[dict] = []
    tail_rows: list[dict] = []
    stratum_totals = {st: [0, 0] for st in AUDIT_STRATA}
    candidate_avoidable_p3 = 0
    accounting_failures = 0
    honest_minus_ref_total = 0

    for r in paired:
        n = jint(r, "n_attempted", 0)
        residual = n - (jint(r, "n_finite", 0) + jint(r, "n_avoidable", 0)
                        + jint(r, "n_honest", 0))
        if not (bool(r.get("accounting_ok")) and residual == 0):
            accounting_failures += 1
        if r.get("method") == "CANDIDATE":
            candidate_avoidable_p3 += jint(r, "n_avoidable", 0)
            ref = jint(r, "ref_nonrepresentable")
            if ref is not None:
                honest_minus_ref_total += abs(jint(r, "n_honest", 0) - ref)
            for st in AUDIT_STRATA:
                t, au = jint(r, f"audit_{st}_total"), jint(r, f"audit_{st}_audited")
                if t is not None:
                    stratum_totals[st][0] += t
                    stratum_totals[st][1] += au or 0

    # The CondRecord streams are the only per-attempt record of a 10^6-attempt run, and
    # conditioning is a statement about the joint distribution of the intended state and
    # success, which no summary preserves.  One file is read at a time and reduced to
    # counts, so the memory cost is one configuration rather than the phase.
    files: dict = {}
    for r in paired:
        rf = r.get("raw_file")
        if not rf:
            continue
        key = (r.get("tag"), r.get("precision"), round(jnum(r, "kappa"), 12),
               jint(r, "seed"))
        files[key] = (rf, jnum(r, "shape_a"))
    if not files:
        raise AnalysisError("phase p3: no `raw_file` was named by any counter row, so the "
                            "per-attempt conditioning records cannot be located")

    acc_bins: dict = defaultdict(lambda: defaultdict(lambda: [0, 0]))
    acc_tail: dict = defaultdict(dict)
    boot_bits: dict = defaultdict(dict)
    edges = list(Q_BIN_EDGES_LOG10)

    for key in sorted(files, key=str):
        tag, precision, kappa, seed = key
        rf, a = files[key]
        _h, arr = read_exp7_records(ctx.resolve_raw(rf), 2)
        log_r = np.asarray(arr["log_r_ref"], dtype=float)
        log_w = S.log_w_from_log_r(log_r)
        with np.errstate(divide="ignore", invalid="ignore"):
            log_q = S.log_q_from_log_w(log_w, a)
        log10_q = log_q / math.log(10.0)
        cats = {"LEGACY": np.asarray(arr["cat_legacy"]),
                "CANDIDATE": np.asarray(arr["cat_candidate"]),
                "QF": np.asarray(arr["cat_qf"])}
        n_att = log_r.size
        stride = max(1, int(math.ceil(n_att / BOOT_MAX_PER_SEED)))
        for method, cat in cats.items():
            success = np.isin(cat, RETURNED_FINITE_CATS)
            for i in range(len(edges) - 1):
                sel = (log10_q <= edges[i]) & (log10_q > edges[i + 1])
                nb = int(sel.sum())
                if nb:
                    acc_bins[(tag, method, precision, kappa)][i][0] += nb
                    acc_bins[(tag, method, precision, kappa)][i][1] += int(
                        success[sel].sum())
                    bin_rows.append(_bin_row(ctx, tag, method, precision, kappa, a, seed, i,
                                             edges, nb, int(success[sel].sum())))
            for p in TAIL_LEVELS:
                target = S.exact_log_r_quantile(p, a)
                above = log_r > target
                cellkey = (tag, method, precision, kappa, p)
                cur = acc_tail[cellkey].setdefault(seed, [0, 0, 0, 0])
                cur[0] += n_att
                cur[1] += int(success.sum())
                cur[2] += int(above.sum())
                cur[3] += int((above & success).sum())
                bb = boot_bits[cellkey].setdefault(seed, {})
                bb.setdefault("mass", []).append(
                    (above & success)[::stride].astype(float))
                bb.setdefault("cond", []).append(
                    np.where(success[::stride], (above & success)[::stride], np.nan))

    for key in sorted(acc_bins, key=str):
        tag, method, precision, kappa = key
        for i, (nb, ns) in sorted(acc_bins[key].items()):
            bin_rows.append(_bin_row(ctx, tag, method, precision, kappa, kappa - 0.5, None,
                                     i, edges, nb, ns, scope="pooled"))
    write_csv(ctx.out("conditioning_bins.csv"), bin_cols, bin_rows)

    for cellkey in sorted(acc_tail, key=str):
        tag, method, precision, kappa, p = cellkey
        by_seed = acc_tail[cellkey]
        n_att = sum(v[0] for v in by_seed.values())
        n_suc = sum(v[1] for v in by_seed.values())
        n_above = sum(v[2] for v in by_seed.values())
        n_ret = sum(v[3] for v in by_seed.values())
        atom = 1.0 - (n_suc / n_att) if n_att else float("nan")
        mass_att = n_ret / n_att if n_att else float("nan")
        mass_cond = n_ret / n_suc if n_suc else float("nan")
        bmass = cluster_value_interval(
            {s: np.concatenate(boot_bits[cellkey][s]["mass"]) for s in boot_bits[cellkey]},
            np.mean, proto, "p3-mass", tag, method, precision, f"{kappa:.12g}", f"{p:g}")
        bcond = cluster_value_interval(
            {s: np.concatenate(boot_bits[cellkey][s]["cond"]) for s in boot_bits[cellkey]},
            nanmean, proto, "p3-cond", tag, method, precision, f"{kappa:.12g}", f"{p:g}")
        tail_rows.append({
            "tag": tag, "method": method, "precision": precision, "kappa": kappa,
            "shape_a": kappa - 0.5, "scope": "pooled", "p": p,
            "log_r_target": S.exact_log_r_quantile(p, kappa - 0.5),
            "n_attempts": n_att, "n_success": n_suc,
            "success_definition": "a three-vector reached the caller: terminal category "
                                  "finite, finite_but_wrong or overflow_returned_finite",
            "n_above_target_all_attempts": n_above,
            "n_above_target_and_returned": n_ret,
            "failure_atom_mass_per_attempt": atom,
            "per_attempt_law_has_failure_atom": bool(atom > 0.0),
            "per_attempt_quantile_reported": False,
            "per_attempt_quantile_reason":
                "the per-attempt law of this method puts an atom of mass "
                f"{atom:.6g} on failure, so its upper quantiles are undefined; the tail MASS "
                "above the target quantile is reported instead",
            "tail_mass_per_attempt": mass_att,
            "tail_mass_per_attempt_lo": bmass["boot_lo"],
            "tail_mass_per_attempt_hi": bmass["boot_hi"],
            "returned_tail_ratio_per_attempt": mass_att / (1.0 - p),
            "returned_tail_ratio_per_attempt_lo": bmass["boot_lo"] / (1.0 - p),
            "returned_tail_ratio_per_attempt_hi": bmass["boot_hi"] / (1.0 - p),
            "tail_mass_conditional_on_success": mass_cond,
            "tail_mass_conditional_on_success_lo": bcond["boot_lo"],
            "tail_mass_conditional_on_success_hi": bcond["boot_hi"],
            "returned_tail_ratio_conditional_on_success": mass_cond / (1.0 - p),
            "returned_tail_ratio_conditional_on_success_lo": bcond["boot_lo"] / (1.0 - p),
            "returned_tail_ratio_conditional_on_success_hi": bcond["boot_hi"] / (1.0 - p),
            "boot_conf": bmass["boot_conf"], "boot_resamples": bmass["boot_resamples"],
            "boot_kind": bmass["boot_kind"], "boot_n_seeds": bmass["boot_n_seeds"],
            "boot_n_per_seed": BOOT_MAX_PER_SEED,
            "protocol_sha256": proto.sha256})
    write_csv(ctx.out("tail_metrics.csv"), tail_cols, tail_rows)

    merged = {st: [p2["stratum_totals"][st][0] + stratum_totals[st][0],
                   p2["stratum_totals"][st][1] + stratum_totals[st][1]]
              for st in AUDIT_STRATA}
    declared = proto.get("audit", "strata_rates")
    cov_cols = ["stratum", "declared_rate", "attempts_in_stratum", "attempts_audited",
                "achieved_rate", "meets_declared_rate", "phases", "protocol_sha256"]
    cov_rows = []
    coverage = {}
    for st in AUDIT_STRATA:
        total, audited = merged[st]
        rate = (audited / total) if total > 0 else float("nan")
        coverage[st] = rate
        cov_rows.append({"stratum": st, "declared_rate": float(declared[st]),
                         "attempts_in_stratum": total, "attempts_audited": audited,
                         "achieved_rate": rate,
                         "meets_declared_rate": (None if total == 0
                                                 else bool(rate + 1e-12 >= declared[st])),
                         "phases": "p2+p3", "protocol_sha256": proto.sha256})
    write_csv(ctx.out("audit_coverage.csv"), cov_cols, cov_rows)

    return {"bin_rows": bin_rows, "tail_rows": tail_rows, "audit_coverage": coverage,
            "candidate_avoidable_p3": candidate_avoidable_p3,
            "accounting_failures": accounting_failures,
            "honest_minus_ref_total": honest_minus_ref_total}


def _bin_row(ctx, tag, method, precision, kappa, a, seed, i, edges, nb, ns,
             scope="seed") -> dict:
    r = {"tag": tag, "method": method, "precision": precision, "kappa": kappa,
         "shape_a": a, "scope": scope, "seed": seed, "bin_kind": "target_upper_tail_q",
         "bin_index": i, "bin_lo": 10.0 ** edges[i + 1] if np.isfinite(edges[i + 1]) else 0.0,
         "bin_hi": 10.0 ** edges[i], "n_in_bin": nb, "n_success": ns,
         "success_definition": "a three-vector reached the caller",
         "protocol_sha256": ctx.proto.sha256}
    ri = rate_columns(ns, nb, "s", ctx.interval_conf)
    r.update({"success_rate": ri["s_rate"], "ci_lo": ri["s_ci_lo"], "ci_hi": ri["s_ci_hi"],
              "interval_kind": ri["s_interval_kind"],
              "is_upper_bound": ri["s_is_upper_bound"]})
    return r


# ---------------------------------------------------------------------------
# P4 -- the complete three-dimensional loader
# ---------------------------------------------------------------------------
# The five structural members frozen with the protocol, then the upper-tail members
# amendment 1.3.0 added.  The tail names are built from `config/protocol.json -> F5_tail.q0`
# rather than written out, so the family is whatever the protocol says it is.
F5_STRUCTURAL_TESTS = ("direction_uniformity", "independence", "frame_invariance",
                       "anisotropy", "cap_law")


def f5_tail_names(q0s) -> tuple:
    return tuple(f"tail_{tail_slug(q0)}_{kind}"
                 for q0 in q0s for kind in ("count", "excess"))


def f5_test_names(q0s) -> tuple:
    return F5_STRUCTURAL_TESTS + f5_tail_names(q0s)


def f5_display(name: str) -> str:
    """The column heading a reader sees; `tail_1em4_count` reads as `tail 1e-4 count`."""
    if name.startswith("tail_"):
        _, slug, kind = name.split("_", 2)
        return f"tail {slug.replace('m', '-')} {kind}"
    return name.replace("_", " ")


def loader_tail_tests(z_intended, z_returned, c_returned, n_attempted: int, q0s,
                      label: str = "cell") -> dict:
    """The upper-tail members of F5, corrected for the representability boundary.

    The derivation, and why the frozen nulls were invalid, are in `exp7_censoring`.

    ``z_intended``   the intended ``Z`` of every recorded attempt, overflowed ones included.
    ``z_returned``   the ``Z`` of the draws the loader actually returned.
    ``c_returned``   each returned draw's own representability threshold ``C(n_i)``.
    ``n_attempted``  the cell's attempt count, which is the count member's denominator.

    COUNT: exact binomial on ``#{intended Z > z0}`` against ``Binomial(n_attempted, q0)``.
    Applies at every threshold of every uncapped cell -- the side condition amendment 1.3.0
    needed is gone, because the quantity it could not observe is observable after all.

    EXCESS: Kolmogorov-Smirnov of the per-draw conditional probability integral transform
    against ``Uniform(0,1)``.  Where nothing can overflow this is the frozen Exp(1) statistic,
    to the last bit.
    """
    zi = np.asarray(z_intended, dtype=float)
    zr = np.asarray(z_returned, dtype=float)
    cr = np.asarray(c_returned, dtype=float)
    out: dict = {}
    n_att = int(n_attempted)
    # A draw whose Z transform underflowed sits above every threshold by construction; it is
    # counted, never dropped, because dropping it would condition on resolvability, which is
    # monotone in the tail -- the exact bias these statistics exist to detect.
    zi_finite = np.isfinite(zi)
    n_zi_unresolved = int(np.sum(~zi_finite))
    zi_ok = zi[zi_finite]
    # The record count must equal the attempt count for an uncapped cell: the probe pushes
    # one record per attempt.  A shortfall means attempts went unrecorded, which is silent
    # conditioning, and it is reported rather than averaged into a rate.
    out["tail_records"] = int(zi.size)
    out["tail_records_missing"] = max(0, n_att - int(zi.size))

    for q0 in q0s:
        slug = tail_slug(q0)
        z0 = -math.log(q0)
        # ---- count, on the intended Z of every attempt --------------------------------
        k = int(np.sum(zi_ok > z0)) + n_zi_unresolved
        resolved = zi.size > 0 and n_att > 0
        out[f"p_tail_{slug}_count"] = (
            float(stats.binomtest(k, n_att, q0).pvalue) if resolved else None)
        out[f"stat_tail_{slug}_count"] = float(k) if resolved else None
        out[f"tail_{slug}_observed_intended"] = int(k) if resolved else None
        out[f"tail_{slug}_expected"] = n_att * q0 if resolved else None
        out[f"tail_{slug}_n"] = n_att if resolved else None
        out[f"tail_{slug}_z0"] = z0
        out[f"tail_{slug}_z_unresolved"] = n_zi_unresolved

        # ---- excess, on the returned draws, against their own censored null ------------
        sel = np.isfinite(zr) & (zr > z0)
        zz, cc = zr[sel], cr[sel]
        out[f"tail_{slug}_observed_returned"] = int(zz.size)
        # A NaN boundary is not the same as an absent one.  `+inf` means "this direction
        # cannot overflow", which is a real and common case; NaN means the geometry could not
        # be evaluated, and treating it as +inf would quietly test a censored draw against an
        # uncensored null -- the very substitution amendment 2.0.0 exists to undo.
        n_bad_c = int(np.sum(np.isnan(cc)))
        if n_bad_c:
            raise AnalysisError(
                f"{label}: {n_bad_c} of {cc.size} returned draws above z0={z0:.4f} have a "
                "NaN representability threshold, so their censored null is undefined. The "
                "cell's geometry (kappa, theta, ub, precision) could not be evaluated; this "
                "is a defect in the run record, not a result.")
        if zz.size >= 8:
            span = np.where(np.isinf(cc), np.inf, cc - z0)
            with np.errstate(divide="ignore", invalid="ignore"):
                u = -np.expm1(-(zz - z0)) / -np.expm1(-span)
            # A returned draw satisfies Z <= C by construction, so U cannot exceed 1 except
            # by the rounding of a recovered direction.  The count that does is published
            # rather than absorbed: if it is ever more than a handful, the boundary model and
            # the loader have stopped agreeing and the member should not be read.
            n_clip = int(np.sum(~((u >= 0.0) & (u <= 1.0))))
            u = np.clip(u, 0.0, 1.0)
            ks = stats.kstest(u, "uniform")
            out[f"p_tail_{slug}_excess"] = float(ks.pvalue)
            out[f"stat_tail_{slug}_excess"] = float(ks.statistic)
            out[f"tail_{slug}_excess_clipped"] = n_clip
            out[f"tail_{slug}_censoring_span_min"] = float(np.min(span))
            out[f"tail_{slug}_censoring_span_max"] = float(np.max(span))
        else:
            out[f"p_tail_{slug}_excess"] = None
            out[f"stat_tail_{slug}_excess"] = None
            out[f"tail_{slug}_excess_clipped"] = 0
            out[f"tail_{slug}_censoring_span_min"] = None
            out[f"tail_{slug}_censoring_span_max"] = None
    return out


TAIL_EXTRA_FIELDS = ("z0", "observed_intended", "observed_returned", "expected", "n",
                     "z_unresolved", "excess_clipped", "censoring_span_min",
                     "censoring_span_max")


def loader_sample(ctx: Context, row: dict):
    """``(log R, unit direction, status, intended log R)`` for one P4 file.

    The fourth value is the probe's ``log_r_ref``: the radius each attempt carried, recorded
    for every attempt including the ones whose velocity overflowed.  It is what makes the
    corrected tail count member possible -- see exp7_censoring -- and it was on disk all
    along; this analysis simply was not reading it.
    """
    rf = row.get("raw_file")
    if not rf:
        raise AnalysisError(f"phase p4: the {row.get('layer')} row for case "
                            f"{row.get('case')} / {row.get('method')} names no raw_file")
    _h, arr = read_exp7_records(ctx.resolve_raw(rf), 3)
    v = np.asarray(arr["v"], dtype=float)
    status = np.asarray(arr["status"])
    log_r_ref = np.asarray(arr["log_r_ref"], dtype=float)
    kappa = jnum(row, "kappa")
    ratio = jnum(row, "theta_ratio", 1.0)
    ub = row.get("ub") or [0.0, 0.0, 1.0]
    log_r, n_hat = IO.recover_radius_direction(v, kappa, 1.0, ratio, ub)
    return log_r, n_hat, status, log_r_ref


def loader_tests(log_r, n_hat, kappa, cap, n_attempted=None, q0s=(), z=None,
                 precision="double", theta_ratio=1.0, label="cell", log_r_ref=None,
                 ub=None, z_intended=None) -> dict:
    """The pre-registered F5 statistics, on one recovered sample.

    Which of them applies depends on the cell, and the ones that do not are reported as not
    applicable rather than as a pass.  Under a cap the accepted set is
    ``R sqrt(kappa) max_j|n_j| <= lambda``: it couples the radius to the direction and is not
    rotationally symmetric in the azimuth, so direction uniformity, independence and frame
    invariance are not properties of the capped conditional law and the cap law is tested by
    its own conditional transform instead.  The set is invariant under permuting the three
    components, so the anisotropy test -- which asks only that the descaled components stay
    exchangeable -- survives the cap and is applied to every cell.

    The upper-tail members are on the same footing: under a cap the accepted radius is
    truncated at a direction-dependent bound and the count above ``z0 = -log q0`` is not
    ``Binomial(n, q0)`` there, so they do not apply to a capped cell.  For an uncapped cell
    amendment 1.4.0 rebuilds both of them against the representability boundary; see
    `exp7_censoring` for the derivation.  ``n_attempted`` is the cell's attempt count, which is what the
    unresolved draws are counted against; it defaults to the number of draws that did
    resolve, i.e. to a cell that reports no loss.  ``z`` lets a caller that already holds
    ``Z`` supply it instead of having it re-derived from ``log_r``.
    """
    ok = np.isfinite(log_r) & np.all(np.isfinite(n_hat), axis=1)
    lr, nh = log_r[ok], n_hat[ok]
    tests = f5_test_names(q0s)
    out = {"n_analyzed": int(lr.size), "n_attempted": int(
        n_attempted if n_attempted is not None else lr.size), "honest_floor_rate": None,
        "tail_records": None, "tail_records_missing": None}
    for t in tests:
        out[f"p_{t}"] = None
        out[f"stat_{t}"] = None
    for q0 in q0s:
        for f in TAIL_EXTRA_FIELDS:
            out[f"tail_{tail_slug(q0)}_{f}"] = None
    if lr.size < 100:
        return out
    cos_theta = nh[:, 2]
    phi = np.arctan2(nh[:, 1], nh[:, 0])
    capped = cap is not None and np.isfinite(cap)

    if not capped:
        g = S.gof_uniform(cos_theta, -1.0, 1.0, "cos_theta")
        if g:
            out["p_direction_uniformity"] = F.simes_global([x.pvalue for x in g])
            out["stat_direction_uniformity"] = g[0].statistic
        ind = S.binned_independence_test(lr, cos_theta)
        out["p_independence"] = ind.pvalue
        out["stat_independence"] = ind.statistic
        gp = S.gof_uniform(phi, -math.pi, math.pi, "phi")
        if gp:
            out["p_frame_invariance"] = F.simes_global([x.pvalue for x in gp])
            out["stat_frame_invariance"] = gp[0].statistic
        if q0s:
            a = float(kappa) - 0.5
            ubv = ub if ub is not None else [0.0, 0.0, 1.0]
            # Z of the draws the loader returned.
            z_ret = (np.asarray(z)[ok] if z is not None
                     else S.z_from_log_w(S.log_w_from_log_r(lr), a))
            # Z of every attempt the probe recorded, overflowed ones included.  `log_r_ref`
            # is the probe's record of the radius each attempt carried; a caller that holds
            # Z itself passes `z_intended`.  Where neither exists -- an injected control
            # whose sample is the returned one by construction -- the returned sample stands
            # in, and the count member then reads exactly the conditioned sample, which is
            # what makes it able to see the conditioning.
            if z_intended is not None:
                z_int = np.asarray(z_intended, dtype=float)
            elif log_r_ref is not None:
                lr_all = np.asarray(log_r_ref, dtype=float)
                z_int = S.z_from_log_w(S.log_w_from_log_r(lr_all[np.isfinite(lr_all)]), a)
            else:
                z_int = z_ret
            c_ret = CEN.censoring_threshold_z(nh, kappa, 1.0, float(theta_ratio), ubv,
                                              precision)
            out.update(loader_tail_tests(z_int, z_ret, c_ret, out["n_attempted"], q0s,
                                         label=label))

    # Anisotropy: after dividing by the declared (sqrt(kappa) theta_perp, ..., sqrt(kappa)
    # theta_par) scale the three components of an isotropic direction are exchangeable, so a
    # mis-applied theta ratio shows up as |n_3| and |n_1| no longer having the same law.  A
    # two-sample rank test, with no moment anywhere: for kappa <= 3/2 the loaded population
    # has no finite variance.
    ks2 = stats.ks_2samp(np.abs(nh[:, 2]), np.abs(nh[:, 0]))
    out["p_anisotropy"] = float(ks2.pvalue)
    out["stat_anisotropy"] = float(ks2.statistic)

    if capped:
        # Conditional transform.  Given the direction, the accepted radius is the target
        # truncated at log R_max = log(lambda) - log(sqrt(kappa)) - log(max_j |n_j|), so
        # U = F(R)/F(R_max) is exactly uniform on (0,1) under the declared bounded law.
        max_abs = np.max(np.abs(nh), axis=1)
        log_rmax = math.log(float(cap)) - 0.5 * math.log(kappa) - np.log(max_abs)
        with np.errstate(divide="ignore", invalid="ignore"):
            f_r = -np.expm1(S.log_q_from_log_w(S.log_w_from_log_r(lr), kappa - 0.5))
            f_max = -np.expm1(S.log_q_from_log_w(S.log_w_from_log_r(log_rmax), kappa - 0.5))
        u = np.where(f_max > 0, f_r / f_max, np.nan)
        u = u[np.isfinite(u) & (u >= 0.0) & (u <= 1.0)]
        if u.size >= 100:
            ks = stats.kstest(u, "uniform")
            out["p_cap_law"] = float(ks.pvalue)
            out["stat_cap_law"] = float(ks.statistic)
    return out


def analyse_p4(ctx: Context) -> dict:
    proto = ctx.proto
    rows = ctx.load_phase("p4")
    cases = list(proto.get("matrix", "loader_cases"))
    seen_cases = sorted({r.get("case") for r in rows if r.get("case")})
    missing = sorted(set(cases) - set(seen_cases))
    if missing:
        raise AnalysisError(f"phase p4: cases {missing} of PROTOCOL.md 4.1 produced no row")

    alpha_F5 = proto.alpha("F5_loader_battery")
    # Amendment 1.3.0: the family gains an exceedance statistic on each loader cell's
    # recovered radial law, at the thresholds `config/protocol.json -> F5_tail.q0`.  It joins
    # the same joint Holm as the other members -- one more way for the candidate to be
    # rejected, not a replacement for anything.
    q0s = [float(q) for q in proto.get("F5_tail", "q0")]
    f5_tests_names = f5_test_names(q0s)
    columns = ["tag", "stdlib", "arch", "execution", "case", "role", "layer", "scope",
               "method", "precision", "kappa", "shape_a", "theta_ratio", "ub_x", "ub_y",
               "ub_z", "cap", "seed", "attempts", "n_returned", "n_analyzed",
               "nonfinite_attempt", "nonfinite_returned", "cap_reject",
               "cap_reject_unrepresentable", "cap_exhausted", "max_attempt_run",
               "engine_calls", "seconds", "conditional", "loss_fraction", "loss_kind",
               "honest_floor_rate", "honest_floor_log10", "in_family", "protocol_sha256"]
    for t in f5_tests_names:
        columns += [f"p_{t}", f"stat_{t}", f"applies_{t}", f"f5_holm_rejected_{t}"]
    for q0 in q0s:
        columns += [f"tail_{tail_slug(q0)}_{f}" for f in TAIL_EXTRA_FIELDS]
    columns += ["tail_records", "tail_records_missing",
                "seed_homogeneity_chi2", "seed_homogeneity_p"]
    for pref in ("nonfinite_attempt_rate",):
        columns += [f"{pref}_count", f"{pref}_n", f"{pref}_rate", f"{pref}_ci_lo",
                    f"{pref}_ci_hi", f"{pref}_interval_kind", f"{pref}_is_upper_bound"]

    measured = [r for r in rows if r.get("layer") in ("uncapped", "capped")]
    out_rows: list[dict] = []
    pooled: dict = defaultdict(lambda: {"log_r": [], "n_hat": [], "log_r_ref": [],
                                        "rows": []})
    unlabelled_conditional = 0

    for r in sorted(measured, key=lambda r: (r.get("tag"), r.get("case"), r.get("method"),
                                             jint(r, "seed"))):
        log_r, n_hat, _status, log_r_ref = loader_sample(ctx, r)
        cap = r.get("cap")
        cap = float(cap) if cap is not None else None
        kappa = jnum(r, "kappa")
        attempts = jint(r, "attempts", 0)
        n_returned = jint(r, "n_returned", 0)
        nf_attempt = jint(r, "nonfinite_attempt", 0)
        nf_returned = jint(r, "nonfinite_returned", 0)
        # An uncapped cell attempts a fixed number of draws and returns one vector per
        # attempt, so its exceedance denominator is `attempts`; a capped cell's unit is the
        # return, and the tail members do not apply to it anyway.
        n_att = attempts if cap is None else n_returned
        tests = loader_tests(log_r, n_hat, kappa, cap, n_attempted=n_att, q0s=q0s,
                             precision=r.get("precision", "double"),
                             theta_ratio=jnum(r, "theta_ratio", 1.0),
                             label=f"{r.get('case')}|{r.get('method')}|seed "
                                   f"{jint(r, 'seed', 0)}",
                             log_r_ref=(log_r_ref if cap is None else None),
                             ub=r.get("ub"))
        n_analyzed = tests["n_analyzed"]
        # Anything that removes a returned draw before a statistic sees it is correlated with
        # the radius here: a non-finite component is exactly what a large radius produces.
        lost = max(0, n_returned - n_analyzed)
        loss_fraction = ((nf_attempt / attempts) if attempts else float("nan")) if cap is None \
            else ((lost / n_returned) if n_returned else float("nan"))
        conditional = bool(lost > 0 or nf_attempt > 0 or nf_returned > 0)
        ub = r.get("ub") or [0.0, 0.0, 1.0]
        floor = honest_floor(kappa, r.get("precision"))
        row = {"tag": r.get("tag"), "stdlib": r.get("stdlib"), "arch": r.get("arch"),
               "execution": r.get("execution"), "case": r.get("case"), "role": r.get("role"),
               "layer": r.get("layer"), "scope": "seed", "method": r.get("method"),
               "precision": r.get("precision"), "kappa": kappa,
               "shape_a": jnum(r, "shape_a"), "theta_ratio": jnum(r, "theta_ratio"),
               "ub_x": float(ub[0]), "ub_y": float(ub[1]), "ub_z": float(ub[2]),
               "cap": cap, "seed": jint(r, "seed"), "attempts": attempts,
               "n_returned": n_returned, "n_analyzed": n_analyzed,
               "nonfinite_attempt": nf_attempt, "nonfinite_returned": nf_returned,
               "cap_reject": jint(r, "cap_reject", 0),
               "cap_reject_unrepresentable": jint(r, "cap_reject_unrepresentable", 0),
               "cap_exhausted": jint(r, "cap_exhausted", 0),
               "max_attempt_run": jint(r, "max_attempt_run"),
               "engine_calls": jint(r, "engine_calls"), "seconds": jnum(r, "seconds"),
               "conditional": conditional,
               "loss_fraction": loss_fraction,
               "loss_kind": ("non-finite attempts of an uncapped run, over attempts"
                             if cap is None else
                             "returned draws the analysis could not recover, over returns"),
               "honest_floor_rate": floor["floor"],
               "honest_floor_log10": floor["log10_floor"],
               "in_family": False, "protocol_sha256": proto.sha256}
        for t in f5_tests_names:
            row[f"p_{t}"] = tests[f"p_{t}"]
            row[f"stat_{t}"] = tests[f"stat_{t}"]
            row[f"applies_{t}"] = tests[f"p_{t}"] is not None
            row[f"f5_holm_rejected_{t}"] = False
        for q0 in q0s:
            for f in TAIL_EXTRA_FIELDS:
                row[f"tail_{tail_slug(q0)}_{f}"] = tests[f"tail_{tail_slug(q0)}_{f}"]
        row["tail_records"] = tests["tail_records"]
        row["tail_records_missing"] = tests["tail_records_missing"]
        row.update(rate_columns(nf_attempt, attempts, "nonfinite_attempt_rate",
                                ctx.interval_conf))
        if conditional and not np.isfinite(loss_fraction):
            unlabelled_conditional += 1
        out_rows.append(row)

        key = (r.get("tag"), r.get("case"), r.get("method"))
        pooled[key]["log_r"].append(log_r)
        pooled[key]["n_hat"].append(n_hat)
        pooled[key]["log_r_ref"].append(log_r_ref)
        pooled[key]["rows"].append(r)

    f5_tests: list[dict] = []
    pooled_rows: dict = {}
    for key in sorted(pooled, key=str):
        tag, case, method = key
        blob = pooled[key]
        r0 = blob["rows"][0]
        cap = r0.get("cap")
        cap = float(cap) if cap is not None else None
        kappa = jnum(r0, "kappa")
        log_r = np.concatenate(blob["log_r"])
        n_hat = np.concatenate(blob["n_hat"])
        log_r_ref = np.concatenate(blob["log_r_ref"])
        attempts = sum(jint(r, "attempts", 0) for r in blob["rows"])
        n_returned = sum(jint(r, "n_returned", 0) for r in blob["rows"])
        tests = loader_tests(log_r, n_hat, kappa, cap,
                             n_attempted=attempts if cap is None else n_returned, q0s=q0s,
                             precision=r0.get("precision", "double"),
                             theta_ratio=jnum(r0, "theta_ratio", 1.0),
                             label=f"{case}|{method}",
                             log_r_ref=(log_r_ref if cap is None else None),
                             ub=r0.get("ub"))
        nf_attempt = sum(jint(r, "nonfinite_attempt", 0) for r in blob["rows"])
        nf_returned = sum(jint(r, "nonfinite_returned", 0) for r in blob["rows"])
        n_analyzed = tests["n_analyzed"]
        lost = max(0, n_returned - n_analyzed)
        loss_fraction = ((nf_attempt / attempts) if attempts else float("nan")) if cap is None \
            else ((lost / n_returned) if n_returned else float("nan"))
        conditional = bool(lost > 0 or nf_attempt > 0 or nf_returned > 0)
        chi2, p_hom = F.seed_homogeneity(
            [jint(r, "nonfinite_attempt", 0) for r in blob["rows"]],
            [jint(r, "attempts", 0) for r in blob["rows"]])
        ub = r0.get("ub") or [0.0, 0.0, 1.0]
        floor = honest_floor(kappa, r0.get("precision"))
        in_family = bool(tag == PRIMARY_TAG and method == "CANDIDATE")
        row = {"tag": tag, "stdlib": r0.get("stdlib"), "arch": r0.get("arch"),
               "execution": r0.get("execution"), "case": case, "role": r0.get("role"),
               "layer": r0.get("layer"), "scope": "pooled", "method": method,
               "precision": r0.get("precision"), "kappa": kappa,
               "shape_a": jnum(r0, "shape_a"), "theta_ratio": jnum(r0, "theta_ratio"),
               "ub_x": float(ub[0]), "ub_y": float(ub[1]), "ub_z": float(ub[2]),
               "cap": cap, "seed": None, "attempts": attempts, "n_returned": n_returned,
               "n_analyzed": n_analyzed, "nonfinite_attempt": nf_attempt,
               "nonfinite_returned": nf_returned,
               "cap_reject": sum(jint(r, "cap_reject", 0) for r in blob["rows"]),
               "cap_reject_unrepresentable": sum(jint(r, "cap_reject_unrepresentable", 0)
                                                 for r in blob["rows"]),
               "cap_exhausted": sum(jint(r, "cap_exhausted", 0) for r in blob["rows"]),
               "max_attempt_run": max(jint(r, "max_attempt_run", 0) for r in blob["rows"]),
               "engine_calls": sum(jint(r, "engine_calls", 0) for r in blob["rows"]),
               "seconds": sum(jnum(r, "seconds", 0.0) for r in blob["rows"]),
               "conditional": conditional, "loss_fraction": loss_fraction,
               "loss_kind": ("non-finite attempts of an uncapped run, over attempts"
                             if cap is None else
                             "returned draws the analysis could not recover, over returns"),
               "honest_floor_rate": floor["floor"],
               "honest_floor_log10": floor["log10_floor"],
               "seed_homogeneity_chi2": chi2, "seed_homogeneity_p": p_hom,
               "in_family": in_family, "protocol_sha256": proto.sha256}
        for t in f5_tests_names:
            row[f"p_{t}"] = tests[f"p_{t}"]
            row[f"stat_{t}"] = tests[f"stat_{t}"]
            row[f"applies_{t}"] = tests[f"p_{t}"] is not None
            row[f"f5_holm_rejected_{t}"] = False
        for q0 in q0s:
            for f in TAIL_EXTRA_FIELDS:
                row[f"tail_{tail_slug(q0)}_{f}"] = tests[f"tail_{tail_slug(q0)}_{f}"]
        row["tail_records"] = tests["tail_records"]
        row["tail_records_missing"] = tests["tail_records_missing"]
        row.update(rate_columns(nf_attempt, attempts, "nonfinite_attempt_rate",
                                ctx.interval_conf))
        if conditional and not np.isfinite(loss_fraction):
            unlabelled_conditional += 1
        pooled_rows[key] = row
        out_rows.append(row)
        if in_family:
            for t in f5_tests_names:
                if tests[f"p_{t}"] is not None:
                    f5_tests.append({"label": f"{case}|{method}|{t}", "p": tests[f"p_{t}"],
                                     "conditional": conditional,
                                     "loss_fraction": loss_fraction})

    # The Holm decision is made once, by the family, and written into the rows -- so that no
    # reader of this CSV and no figure has to re-derive a rejection threshold from an alpha
    # and a count of tests, which is a second place for the rule to live.
    family_F5 = F.family_F5(proto, f5_tests)
    rejected = set(family_F5.offenders)
    for key, row in pooled_rows.items():
        for t in f5_tests_names:
            row[f"f5_holm_rejected_{t}"] = bool(f"{key[1]}|{key[2]}|{t}" in rejected)
    write_csv(ctx.out("loader_validation.csv"), columns, out_rows)
    write_validation_matrix(ctx, pooled_rows, family_F5, alpha_F5, f5_tests_names)
    return {"rows": out_rows, "pooled": pooled_rows, "F5": family_F5,
            "unlabelled_conditional": unlabelled_conditional, "f5_tests": f5_tests,
            "f5_test_names": list(f5_tests_names), "q0s": q0s}


def write_validation_matrix(ctx: Context, pooled_rows: dict, family_F5, alpha: float,
                            f5_tests_names) -> None:
    rejected = set(family_F5.offenders)
    L = ["# Experiment 7 - complete-loader validation matrix", "",
         f"Protocol `config/protocol.json` SHA-256 `{ctx.proto.sha256}`.", "",
         "Rows are the curated configurations of PROTOCOL.md 4.1, pooled over the frozen "
         "seeds; columns are the tests of family F5 -- the five frozen with the protocol, "
         "then the upper-tail members at each threshold in "
         "`F5_tail.q0`. Decisions are Holm-"
         f"corrected jointly over **all** cells and tests at a familywise alpha of {alpha:g}, "
         "not within each cell.", "",
         "`n.a.` marks a test that is not a property of the cell's law rather than one that "
         "passed: under a cap the accepted set couples the radius to the direction and is "
         "not rotationally symmetric in the azimuth, so direction uniformity, independence "
         "and frame invariance do not apply there and the cap law is tested by its own "
         "conditional transform instead.", "",
         "**A conditional cell cannot support the fidelity claim on its own.** A cell is "
         "conditional when the statistics ran on a subsample that something correlated with "
         "the radius had already filtered -- a non-finite component is exactly what a large "
         "radius produces. Every such cell carries its loss fraction in the last column and "
         "in `loader_validation.csv`.", "",
         "The two upper-tail members are amendment 1.3.0's, with the nulls corrected in "
         "amendment 2.0.0. The **count** is the number of ATTEMPTS whose intended `Z` "
         "exceeds `z0 = -log q0`, taken from the radius the probe recorded for every "
         "attempt including the overflowed ones, against `Binomial(attempts, q0)` -- exact, "
         "and applicable at every threshold, which is why no cell is marked not-applicable "
         "for a count any more. The **excess** is a Kolmogorov-Smirnov test of the returned "
         "draws above `z0`, each transformed by its own representability threshold "
         "`C(n) = Z(log max() - log max_j|g_j(n)|)`; where nothing can overflow that "
         "transform is the identity on the frozen `Exp(1)` statistic. `loader_validation."
         "csv` publishes, per cell and threshold, the intended count, the returned count, "
         "the expected count, the span of the censoring interval, and how many transformed "
         "values had to be clipped to [0,1] -- which should be none.", "",
         "| case | method | tag | " + " | ".join(f5_display(t) for t in f5_tests_names)
         + " | n analysed | conditional |",
         "|---|---|---|" + "---|" * (len(f5_tests_names) + 2)]
    for key in sorted(pooled_rows, key=lambda k: (k[1], k[2], k[0])):
        tag, case, method = key
        r = pooled_rows[key]
        cells = []
        for t in f5_tests_names:
            p = r.get(f"p_{t}")
            if p is None:
                cells.append("n.a.")
            elif f"{case}|{method}|{t}" in rejected:
                cells.append(f"**FAIL** (p={p:.3g})")
            elif not r["in_family"]:
                cells.append(f"p={p:.3g} (not in F5)")
            else:
                cells.append(f"pass (p={p:.3g})")
        cond = ("conditional, loss fraction "
                f"{r['loss_fraction']:.3g}" if r["conditional"] else "unconditional")
        L.append(f"| {case} | {method} | {tag} | " + " | ".join(cells)
                 + f" | {r['n_analyzed']} | {cond} |")
    L += ["", "Only the CANDIDATE cells of the primary environment "
              f"(`{PRIMARY_TAG}`) enter the F5 decision: within the 2.x line the stream is a "
              "function of the engine alone, so the other environment's rows are the same "
              "draws and would enter Holm twice. The LEGACY rows are the comparator and are "
              "printed with their p-values but are not gated -- where the released 1.0.0 form "
              "discards a non-negligible share of the intended draws, the draws it keeps are "
              "no longer distributed as the target, and a test detecting that is the "
              "conditioning result rather than a defect in this analysis.", "",
          f"F5 global: {'passes' if family_F5.passed else 'FAILS'}. {family_F5.detail}", ""]
    with open(ctx.out("validation_matrix.md"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(L) + "\n")


# ---------------------------------------------------------------------------
# Negative controls, by measured power
# ---------------------------------------------------------------------------
def analyse_negative_controls(ctx: Context, p1: dict, p4: dict) -> dict:
    """Power of the frozen battery against each injected defect, not a pass/fail bit.

    Each control's detection rule is the frozen family rule itself, evaluated on the injected
    sample: a control counts as detected when the family that certifies the corresponding
    property *fails*.  That keeps the measured power a property of the battery rather than of
    a statistic chosen here.

    Every replicate is drawn so that the measurement's own null is exact.  NC1 permutes the
    measured cell's radii against its directions, which is the independence null with the
    cell's own marginals; NC3 simulates the radial law outright.  Neither resamples a
    measured sample WITH replacement: the empirical joint of a finite sample carries a real
    dependence of order its own chi-square, so a with-replacement resample inherits it as a
    non-centrality and rejects at close to certainty with nothing injected at all.  The
    ``level`` row published for NC1 is the measurement of exactly that quantity, so the power
    figures below can be read against a stated false-positive rate rather than assumed.
    """
    proto = ctx.proto
    required = float(proto.get("F6", "required_power"))
    effects = [float(x) for x in proto.get("F6", "injection_loss_fractions")]
    n_inj = int(proto.get("F6", "injections_per_effect"))
    q0s = [float(q) for q in proto.get("F5_tail", "q0")]
    f5_names = f5_test_names(q0s)
    nc3_rule = proto.get("corrections", "NC3_effect")

    rows: list[dict] = []
    controls: list[dict] = []

    def f5_cell(tests: dict, label: str, q: float, names=f5_names) -> list[dict]:
        return [{"label": f"{label}|{k}", "p": tests[f"p_{k}"], "conditional": True,
                 "loss_fraction": q} for k in names if tests.get(f"p_{k}") is not None]

    # --- NC1: radius-direction coupling at the certified loss fractions --------------
    src = None
    for key, r in sorted(p4["pooled"].items(), key=str):
        if key[0] == PRIMARY_TAG and key[2] == "CANDIDATE" and key[1] == "C1":
            src = r
    if src is None:
        rows.append({"control": "NC1_radius_direction_coupling", "note":
                     "case C1 produced no pooled CANDIDATE cell on the primary environment"})
    nc1_sample = None
    if src is not None:
        nc1_sample = _pooled_loader_sample(ctx, "C1", "CANDIDATE")
    kappa1 = jnum(src, "kappa") if src is not None else float("nan")
    nc1_ratio = jnum(src, "theta_ratio", 1.0) if src is not None else 1.0
    nc1_ub = (src.get("ub") or [0.0, 0.0, 1.0]) if src is not None else [0.0, 0.0, 1.0]
    nc1_precision = (src.get("precision") or "double") if src is not None else "double"

    # The effect list carries a leading 0: the same procedure with nothing removed, which
    # measures the rejection rate of the measurement itself.  A power figure quoted without
    # it is a number whose null is unstated.
    for q in [0.0] + effects:
        det_silent, det_reported, n_used = 0, 0, 0
        if nc1_sample is not None:
            lr, nh = nc1_sample
            n_used = lr.size
            for i in range(n_inj):
                rng = np.random.default_rng(stable_seed("NC1", f"{q:g}", i))
                # Permute the radii against the directions.  That is the independence null
                # carrying the cell's own two marginals, it uses every measured draw exactly
                # once, and it gives each replicate a fresh joint sample for the injection to
                # act on -- which injecting into one fixed sample 200 times would not.
                lr_p = lr[rng.permutation(lr.size)]
                lr_i, nh_i = _inject_radius_direction_coupling(lr_p, nh, q, rng)
                # The injected defect is a loader that silently conditions: the draws it
                # loses never reach the caller and nothing in its bookkeeping records them,
                # so the cell it hands the analysis reports no loss.  That is the case the
                # battery has to catch on its own.  A loss the loader DOES report arrives as
                # unresolved attempts, is counted toward every threshold by
                # `F5_tail.unresolved_Z_counts_toward_every_threshold`, and leaves
                # the cell labelled conditional with its loss fraction under PROTOCOL.md 5.2.
                # Both readings are measured; only the first is the control.
                t = loader_tests(lr_i, nh_i, kappa1, None,
                                 n_attempted=lr_i.size, q0s=q0s,
                                 precision=nc1_precision, theta_ratio=nc1_ratio,
                                 ub=nc1_ub)
                if not F.family_F5(proto, f5_cell(t, "NC1", q)).passed:
                    det_silent += 1
                # The same injection read as a loss the loader DOES report: the removed
                # draws are restored to the count member's denominator and to its intended
                # sample, which is what a loader that notices its own losses would supply.
                t_rep = dict(t)
                t_rep.update(loader_tests(lr_i, nh_i, kappa1, None, n_attempted=lr_p.size,
                                          q0s=q0s, precision=nc1_precision,
                                          theta_ratio=nc1_ratio, ub=nc1_ub,
                                          z_intended=S.z_from_log_w(
                                              S.log_w_from_log_r(lr_p), kappa1 - 0.5)))
                if not F.family_F5(proto, f5_cell(t_rep, "NC1", q)).passed:
                    det_reported += 1
        level = (q == 0.0)
        row = _power_row(
            ctx, "NC1_radius_direction_coupling", "level" if level else "loss_fraction", q,
            "F5 loader battery including the upper-tail members of amendment 1.3.0, Holm "
            "jointly over the cell's applicable tests",
            n_inj if nc1_sample is not None else 0, det_silent,
            None if level else required, not level,
            f"C1 CANDIDATE, {n_used} pooled returned samples",
            ("no draw is removed; this row is the rejection rate of the measurement itself, "
             "against which the rows below are read"
             if level else
             "the q fraction of draws with the largest R*max_j|n_j| is removed and the cell "
             "reports no loss, which is the radius-direction coupling a finite-precision "
             "loss induces in a loader that silently conditions"))
        rows.append(row)
        if not level:
            controls.append({"name": "NC1_radius_direction_coupling", "effect": q,
                             "injections": row["injections"], "detections": det_silent})
        rows.append(_power_row(
            ctx, "NC1_radius_direction_coupling",
            "level" if level else "loss_fraction", q,
            "the same F5 battery, with the removed draws declared as unresolved attempts",
            n_inj if nc1_sample is not None else 0, det_reported,
            None if level else required, False,
            f"C1 CANDIDATE, {n_used} pooled returned samples",
            "the same injection, read as a loss the loader reports: the removed draws are "
            "counted toward every threshold and the denominator stays at the pre-injection "
            "attempt count, so the exceedance count is restored wherever the loss fraction "
            "is at or below q0"))

    # --- NC2: the capped sample against the uncapped law, at the weakest cap ---------
    weak = None
    for key, r in sorted(p4["pooled"].items(), key=str):
        if key[0] == PRIMARY_TAG and key[2] == "CANDIDATE" and r.get("cap") is not None:
            if weak is None or float(r["cap"]) > float(weak["cap"]):
                weak = r
    det2, n2, cap2 = 0, 0, None
    if weak is not None:
        cap2 = float(weak["cap"])
        lr, nh = _pooled_loader_sample(ctx, weak["case"], "CANDIDATE")
        a = jnum(weak, "kappa") - 0.5
        z = S.z_from_log_w(S.log_w_from_log_r(lr), a)
        z = z[np.isfinite(z)]
        n2 = z.size
        for i in range(n_inj):
            rng = np.random.default_rng(stable_seed("NC2", weak["case"], i))
            zi = z[rng.integers(0, z.size, z.size)] if z.size else z
            a2, p_ad = F.anderson_darling_exp1(zi)
            cell = [{"label": "NC2", "pvalues": [
                p_ad, ks_pvalue(float(stats.kstest(zi, "expon").statistic), zi.size),
                cvm_pvalue(float(stats.cramervonmises(zi, "expon").statistic))]}]
            if not F.family_F1(proto, cell).passed:
                det2 += 1
    row2 = _power_row(ctx, "NC2_capped_vs_uncapped_weak_cap", "cap_lambda", cap2,
                      "F1 radial law, Simes over AD/KS/CvM on Z",
                      n_inj if weak is not None else 0, det2, required, True,
                      f"{weak['case'] if weak else 'n.a.'} CANDIDATE, {n2} returned samples",
                      "the bounded sample is tested against the UNCAPPED law at the weakest "
                      "declared cap, which is the smallest version of this effect the "
                      "protocol declares",
                      "replicates are drawn from the measured capped sample WITH "
                      "replacement, which inflates the null of the three F1 statistics; the "
                      "detection rate below is therefore an upper bound on the power at this "
                      "cap, and the effect is far larger than the inflation")
    rows.append(row2)
    controls.append({"name": "NC2_capped_vs_uncapped_weak_cap", "effect": cap2,
                     "injections": row2["injections"], "detections": det2})

    # --- NC3: survivor conditioning, measured against the representability boundary ---
    # corrections.NC3_effect picks the worst cell among those the fidelity claim rests on --
    # the uncapped loader cells.  What CHANGED in amendment 1.4.0 is what is injected into
    # it, and the reason is that the old injection was not an injection of a defect.
    #
    # The frozen NC3 removed every draw above a single cutoff `-log q` from a pure Exp(1)
    # sample and required the battery to detect it.  But that is, to within the direction
    # dependence, exactly what honest overflow does to a CORRECT loader: near kappa = 1/2 the
    # law puts probability outside the type, and a loader that declines to return it is
    # behaving correctly.  The frozen battery "detected" it because its null was the
    # untruncated law, i.e. because it rejected correct behaviour -- its measured NC3 power
    # was Type-I error wearing a power label, and the same defect failed the candidate on
    # cases C3 and C4.
    #
    # NC3 is therefore split into the two questions that were tangled together:
    #
    #   NC3a  honest censoring ALONE, at the cell's own direction-dependent boundary,
    #         correctly reported.  The battery must NOT reject.  This is a level, not a
    #         power, and it is published as such.
    #   NC3b  honest censoring PLUS an extra silent conditioning of the draws the loader
    #         could have returned.  That is the defect, and the battery must detect it with
    #         power >= F6.required_power at the pre-registered effect sizes in
    #         F6.nc3_excess_loss_fractions.
    #
    # Both are simulated from the exact law by its own generative construction -- the
    # shape-boosting identity of config/honest_floor.md -- so the null is exact rather than
    # approximated, and the direction is drawn uniformly and independently of the radius,
    # which is what a correct loader produces.
    worst = None
    for key, r in sorted(p4["pooled"].items(), key=str):
        if key[0] != PRIMARY_TAG or key[2] != "CANDIDATE" or r.get("cap") is not None:
            continue
        lf = r.get("loss_fraction")
        lf = float(lf) if lf is not None and np.isfinite(float(lf)) else 0.0
        if worst is None or lf > worst[0]:
            worst = (lf, r)
    worst_q = worst[0] if worst else 0.0
    worst_row = worst[1] if worst else None
    worst_n = jint(worst_row, "attempts", 0) if worst_row is not None else 0
    worst_label = (f"{worst_row['case']} CANDIDATE {worst_row['precision']} "
                   f"kappa={worst_row['kappa']:g}" if worst_row is not None else None)
    nc3_excess = [float(x) for x in proto.get("F6", "nc3_excess_loss_fractions")]

    def nc3_replicate(rng, n_attempt, kappa, ratio, ub, precision, extra):
        """One cell's worth of attempts: intended Z, which of them the type allows back,
        and each returned draw's own representability threshold."""
        a = float(kappa) - 0.5
        log_r = 0.5 * (np.log(rng.gamma(1.5, 1.0, size=n_attempt))
                       - np.log(rng.gamma(a + 1.0, 1.0, size=n_attempt))
                       + rng.exponential(size=n_attempt) / a)
        gv = rng.normal(size=(n_attempt, 3))
        n_hat = gv / np.sqrt(np.einsum("ij,ij->i", gv, gv))[:, None]
        log_m = np.log(CEN.max_component(n_hat, kappa, 1.0, ratio, ub))
        log_max = math.log(CEN.MAX_FINITE[precision])
        returned = (log_r + log_m) <= log_max
        recorded = np.ones(n_attempt, dtype=bool)
        if extra > 0.0:
            m = int(round(extra * n_attempt))
            elig = np.flatnonzero(returned)
            if m > 0 and elig.size > m:
                # Condition on the largest component, which is the quantity representability
                # acts on, among the draws the type would have allowed back.  The injected
                # loss is therefore strictly IN EXCESS of the honest floor rather than
                # overlapping it.
                key_ = (log_r + log_m)[elig]
                drop = elig[np.argpartition(key_, elig.size - m)[elig.size - m:]]
                recorded[drop] = False
                returned[drop] = False
        return log_r[recorded], returned[recorded], n_hat[recorded], n_attempt

    nc3_cases = [("level", 0.0)] + [("in_family", q) for q in nc3_excess]
    for kind, q in nc3_cases:
        det, n_used = 0, int(worst_n)
        if worst_row is not None and n_used > 0:
            kappa3 = jnum(worst_row, "kappa")
            ratio3 = jnum(worst_row, "theta_ratio", 1.0)
            ub3 = worst_row.get("ub") or [0.0, 0.0, 1.0]
            prec3 = worst_row.get("precision") or "double"
            a3 = kappa3 - 0.5
            for i in range(n_inj):
                rng = np.random.default_rng(stable_seed("NC3", f"{q:.12g}", i))
                lr_rec, ret, nh_rec, n_att = nc3_replicate(rng, n_used, kappa3, ratio3, ub3,
                                                           prec3, q)
                if int(np.sum(ret)) < 100:
                    continue
                z_int = S.z_from_log_w(S.log_w_from_log_r(lr_rec), a3)
                t = loader_tests(lr_rec[ret], nh_rec[ret], kappa3, None,
                                 n_attempted=n_att, q0s=q0s, precision=prec3,
                                 theta_ratio=ratio3, ub=ub3, z_intended=z_int)
                if not F.family_F5(proto, f5_cell(t, "NC3", q)).passed:
                    det += 1
        row = _power_row(
            ctx, "NC3_survivor_conditioning",
            "level" if kind == "level" else "excess_loss_fraction", q,
            "F5 loader battery including the upper-tail members of amendment 1.4.0, Holm "
            "jointly over the cell's applicable tests",
            n_inj if (worst_row is not None and n_used > 0) else 0, det,
            None if kind == "level" else required, kind == "in_family",
            f"{worst_label}; n = {n_used}; measured loss fraction {worst_q:.3g}",
            ("honest overflow alone, at the cell's own direction-dependent representability "
             "boundary, correctly reported: a CORRECT loader, so this row is the battery's "
             "false-rejection rate and carries no power threshold"
             if kind == "level" else
             "honest overflow, plus a further q fraction of the draws the type WOULD have "
             "allowed back, removed silently -- conditioning in excess of what the floor "
             "forces, which is the defect the battery exists to exclude"),
            nc3_rule if kind == "in_family" else None)
        rows.append(row)
        if kind == "in_family":
            controls.append({"name": "NC3_survivor_conditioning", "effect": q,
                             "injections": row["injections"], "detections": det})

    cols = ["control", "effect_kind", "effect_size", "detection_statistic", "injections",
            "detections", "power", "power_ci_lo", "power_ci_hi", "power_interval_kind",
            "power_is_upper_bound", "required_power", "passed", "in_family", "source_cell",
            "injection_definition", "note", "protocol_sha256"]
    write_csv(ctx.out("negative_controls.csv"), cols, rows)
    return {"rows": rows, "controls": controls,
            "F6": F.family_F6(proto, controls)}



def _power_row(ctx, name, effect_kind, effect, statistic, n, det, required, in_family,
               source, definition, note=None) -> dict:
    """One measured rejection rate with its interval.

    ``required is None`` marks a row whose quantity is not a power -- the level row, where
    nothing is injected.  It then carries no threshold and no verdict, because a rejection
    rate measured under the null is neither passed nor failed by the power rule.
    """
    ri = rate_columns(det, n, "power", ctx.interval_conf)
    return {"control": name, "effect_kind": effect_kind, "effect_size": effect,
            "detection_statistic": statistic, "injections": n, "detections": det,
            "power": ri["power_rate"], "power_ci_lo": ri["power_ci_lo"],
            "power_ci_hi": ri["power_ci_hi"],
            "power_interval_kind": ri["power_interval_kind"],
            "power_is_upper_bound": ri["power_is_upper_bound"],
            "required_power": required,
            "passed": (None if required is None
                       else bool(n > 0 and (ri["power_ci_lo"] or 0.0) >= required)),
            "in_family": bool(in_family), "source_cell": source,
            "injection_definition": definition, "note": note,
            "protocol_sha256": ctx.proto.sha256}


_LOADER_CACHE: dict = {}


def _pooled_loader_sample(ctx: Context, case: str, method: str):
    key = (case, method)
    if key in _LOADER_CACHE:
        return _LOADER_CACHE[key]
    rows = [r for r in ctx.load_phase("p4")
            if r.get("case") == case and r.get("method") == method
            and r.get("tag") == PRIMARY_TAG and r.get("layer") in ("uncapped", "capped")]
    lrs, nhs = [], []
    for r in sorted(rows, key=lambda r: jint(r, "seed")):
        lr, nh, _, _ = loader_sample(ctx, r)
        ok = np.isfinite(lr) & np.all(np.isfinite(nh), axis=1)
        lrs.append(lr[ok])
        nhs.append(nh[ok])
    out = (np.concatenate(lrs), np.concatenate(nhs)) if lrs else (np.zeros(0),
                                                                  np.zeros((0, 3)))
    _LOADER_CACHE[key] = out
    return out


def _inject_radius_direction_coupling(log_r, n_hat, q, rng):
    """Remove the ``q`` fraction of draws with the largest ``R max_j |n_j|``.

    That is the coupling a finite-precision loss actually creates: the draws a loader cannot
    return are the ones whose largest component overflows, which is a joint property of the
    radius and the direction and of neither alone.  A random tie-break keeps the injection
    from being a deterministic function of the sample.
    """
    n = log_r.size
    m = int(round(q * n))
    if m <= 0:
        return log_r, n_hat
    key = log_r + np.log(np.max(np.abs(n_hat), axis=1)) + rng.random(n) * 1e-12
    drop = np.argpartition(key, n - m)[n - m:]
    keep = np.ones(n, dtype=bool)
    keep[drop] = False
    return log_r[keep], n_hat[keep]


# ---------------------------------------------------------------------------
# P5 -- portability
# ---------------------------------------------------------------------------
PORT_COLUMNS = ["kind", "label", "arch", "os", "stdlib", "compiler", "tag", "execution",
                "available", "completed", "reason_not_run", "n_rows", "prediction",
                "n_common", "n_digests_compared", "n_differences", "identical",
                "log_ratio", "se", "margin", "p_value", "equivalent", "disagrees",
                "informative", "in_decision", "protocol_sha256"]


def analyse_p5(ctx: Context, p1: dict) -> dict:
    proto = ctx.proto
    rows = ctx.load_phase("p5")
    klass = [r for r in rows if r.get("layer") == "class" and r.get("method") == "CANDIDATE"]
    if not klass:
        raise AnalysisError("phase p5: no CANDIDATE `class` rows; the same-architecture F7 "
                            "statistic is the digest over those rows")
    check_ladder_coverage(ctx, "p5", klass, ("CANDIDATE",))

    def digests(src, method):
        out = defaultdict(dict)
        for r in src:
            if r.get("method") != method:
                continue
            out[r.get("tag")][(r.get("precision"), round(jnum(r, "kappa"), 12),
                               jint(r, "seed"))] = (r.get("digest_sha256"),
                                                    r.get("digest_fnv1a64"),
                                                    jint(r, "n_finite"),
                                                    jint(r, "nonfinite_output"),
                                                    jint(r, "n_attempts_reported"),
                                                    jint(r, "n_nonfinite_reported"))
        return out

    cand = digests(klass, "CANDIDATE")
    legacy = digests([r for r in p1["class_rows"]], "LEGACY")
    env_by_tag = {}
    for r in rows:
        env_by_tag.setdefault(r.get("tag"), r)

    out_rows: list[dict] = []
    stdlib_pairs: list[dict] = []
    tags = sorted(cand)
    for i in range(len(tags)):
        for j in range(i + 1, len(tags)):
            ta, tb = tags[i], tags[j]
            ea, eb = env_by_tag.get(ta, {}), env_by_tag.get(tb, {})
            if ea.get("arch") != eb.get("arch"):
                continue
            for method, table, prediction, in_dec in (
                    ("CANDIDATE", cand, "bitwise equality", True),
                    ("LEGACY", legacy, "difference expected: the 1.0.0 form draws its Gamma "
                                       "from the standard library", False)):
                a, b = table.get(ta, {}), table.get(tb, {})
                common = sorted(set(a) & set(b), key=str)
                diffs = sum(1 for k in common if a[k] != b[k])
                identical = bool(common and diffs == 0)
                label = f"{ea.get('arch')}: {ta} vs {tb} [{method}]"
                out_rows.append({
                    "kind": "cross_stdlib_bitwise", "label": label, "arch": ea.get("arch"),
                    "os": None, "stdlib": f"{ea.get('stdlib')} vs {eb.get('stdlib')}",
                    "compiler": f"{ea.get('compiler')} vs {eb.get('compiler')}",
                    "tag": f"{ta}|{tb}", "execution": ea.get("execution"),
                    "available": True, "completed": True, "reason_not_run": None,
                    "n_rows": len(common), "prediction": prediction,
                    "n_common": len(common), "n_digests_compared": len(common),
                    "n_differences": diffs, "identical": identical,
                    "in_decision": in_dec, "protocol_sha256": proto.sha256})
                if in_dec:
                    stdlib_pairs.append({"label": label, "identical": identical})

    environments = []
    for tag in sorted(env_by_tag):
        e = env_by_tag[tag]
        n_rows = sum(1 for r in klass if r.get("tag") == tag)
        rec = {"arch": e.get("arch"), "os": None, "compiler": e.get("compiler"),
               "stdlib": e.get("stdlib"), "execution": e.get("execution"),
               "available": True, "completed": bool(n_rows), "n_rows": n_rows,
               "reason_not_run": None if n_rows else "no P5 counter rows for this tag"}
        environments.append(rec)
        out_rows.append({"kind": "environment", "label": f"{e.get('arch')}/{tag}",
                         "arch": e.get("arch"), "os": None, "stdlib": e.get("stdlib"),
                         "compiler": e.get("compiler"), "tag": tag,
                         "execution": e.get("execution"), "available": True,
                         "completed": bool(n_rows), "reason_not_run": rec["reason_not_run"],
                         "n_rows": n_rows, "prediction": None, "in_decision":
                             e.get("execution") == "native",
                         "protocol_sha256": proto.sha256})

    # Environments the G4 matrix expects and this host cannot provide.  An absent environment
    # is recorded as absent, with the reason, and never as a zero or a pass.
    expected_path = os.path.join(ROOT, ".github", "ci", "expected_environments.json")
    if os.path.exists(expected_path):
        with open(expected_path, encoding="utf-8") as fh:
            expected = json.load(fh)
        have_arch = {e["arch"] for e in environments}
        for spec in expected:
            for etag in spec.get("tags", []):
                if spec.get("arch") in have_arch and etag in env_by_tag:
                    continue
                reason = ("this host provides neither that architecture nor that standard "
                          "library; the exact command is in results/portability_remote.md, "
                          "and the CI workflow .github/workflows/portability.yml runs it")
                if not spec.get("in_scope_protocol_2_0_0", True):
                    reason = ("out of scope under protocol 2.0.0, which limits acceptance "
                              "to the supported environment and withdraws the "
                              "cross-architecture claim; recorded so that what a broader "
                              "claim would require stays visible. " + reason)
                environments.append({"arch": spec.get("arch"), "os": spec.get("os"),
                                     "compiler": None, "stdlib": etag, "execution": None,
                                     "available": False, "completed": False, "n_rows": 0,
                                     "reason_not_run": reason})
                out_rows.append({
                    "kind": "environment", "label": f"{spec.get('env_id')}/{etag}",
                    "arch": spec.get("arch"), "os": spec.get("os"), "stdlib": etag,
                    "compiler": None, "tag": etag, "execution": None, "available": False,
                    "completed": False, "reason_not_run": reason, "n_rows": 0,
                    "prediction": None, "in_decision": False,
                    "protocol_sha256": proto.sha256})

    fam = F.family_F7(proto, stdlib_pairs, [])
    write_csv(ctx.out("portability.csv"), PORT_COLUMNS, out_rows)
    return {"rows": out_rows, "F7": fam, "environments": environments}


# ---------------------------------------------------------------------------
# P6 -- cost
# ---------------------------------------------------------------------------
def analyse_p6(ctx: Context) -> dict:
    proto = ctx.proto
    rows = ctx.load_phase("p6")
    timed = [r for r in rows if r.get("layer") in ("uncapped", "capped")]
    repro = [r for r in rows if r.get("layer") == "reproducibility"]
    if not timed:
        raise AnalysisError("phase p6: no timed block rows")
    if not repro:
        raise AnalysisError("phase p6: no reproducibility rows; G5 requires same-seed "
                            "reproducibility to be checked, not assumed")
    cases = list(proto.get("matrix", "benchmark_cases"))
    missing = sorted(set(cases) - {r.get("case") for r in timed})
    if missing:
        raise AnalysisError(f"phase p6: benchmark cases {missing} produced no timed block")

    cols = ["scope", "tag", "stdlib", "arch", "case", "method", "precision", "kappa", "cap",
            "block", "slot", "candidate_first", "seed", "n_returned", "n_finite",
            "n_nonfinite", "attempts_reported", "seconds", "seconds_per_returned",
            "ratio_candidate_over_legacy", "ratio_lo", "ratio_hi", "ratio_conf",
            "boot_kind", "boot_resamples", "boot_rng_seed", "boot_n_seeds",
            "pre_registered_bound", "repeat_identical", "protocol_sha256"]
    out_rows: list[dict] = []
    per_block: dict = defaultdict(dict)

    for r in sorted(timed, key=lambda r: (r.get("tag"), r.get("case"), jint(r, "block"),
                                          r.get("method"))):
        out_rows.append({
            "scope": "block", "tag": r.get("tag"), "stdlib": r.get("stdlib"),
            "arch": r.get("arch"), "case": r.get("case"), "method": r.get("method"),
            "precision": r.get("precision"), "kappa": jnum(r, "kappa"), "cap": r.get("cap"),
            "block": jint(r, "block"), "slot": jint(r, "slot"),
            "candidate_first": bool(r.get("candidate_first")), "seed": jint(r, "seed"),
            "n_returned": jint(r, "n_returned"), "n_finite": jint(r, "n_finite"),
            "n_nonfinite": jint(r, "n_nonfinite"),
            "attempts_reported": jint(r, "attempts_reported"),
            "seconds": jnum(r, "seconds"),
            "seconds_per_returned": jnum(r, "seconds_per_returned"),
            "protocol_sha256": proto.sha256})
        per_block[(r.get("tag"), r.get("case"), jint(r, "block"))][r.get("method")] = r

    # corrections.G5_interval: the interval on the ratio is the cluster bootstrap at
    # `statistics.bootstrap.conf`, which is 0.99.  PROTOCOL.md 6 says 95 per cent in prose;
    # the machine-readable value governs and is the conservative choice for a gate on an
    # upper limit.  `cluster_value_interval` reads it from there, and every row records the
    # level it used in `ratio_conf`.
    bound = float(proto.get("gates", "G5", "time_per_return_ratio_bound"))
    by_case: dict = defaultdict(lambda: defaultdict(list))
    overall: dict = defaultdict(list)
    for (tag, case, block), pair in sorted(per_block.items(), key=str):
        c, l = pair.get("CANDIDATE"), pair.get("LEGACY")
        if not c or not l:
            continue
        cs, ls = jnum(c, "seconds_per_returned"), jnum(l, "seconds_per_returned")
        if not (np.isfinite(cs) and np.isfinite(ls) and ls > 0):
            continue
        seed = jint(c, "seed")
        by_case[(tag, case)][seed].append(cs / ls)
        if tag == PRIMARY_TAG:
            overall[seed].append(cs / ls)

    for key in sorted(by_case, key=str):
        tag, case = key
        bi = cluster_value_interval(by_case[key], np.median, proto, "p6", tag, case)
        r0 = next(r for r in timed if r.get("tag") == tag and r.get("case") == case)
        out_rows.append({
            "scope": "case_ratio", "tag": tag, "stdlib": r0.get("stdlib"),
            "arch": r0.get("arch"), "case": case, "method": "CANDIDATE/LEGACY",
            "precision": r0.get("precision"), "kappa": jnum(r0, "kappa"),
            "cap": r0.get("cap"), "ratio_candidate_over_legacy": bi["boot_estimate"],
            "ratio_lo": bi["boot_lo"], "ratio_hi": bi["boot_hi"],
            "ratio_conf": bi["boot_conf"], "boot_kind": bi["boot_kind"],
            "boot_resamples": bi["boot_resamples"], "boot_rng_seed": bi["boot_rng_seed"],
            "boot_n_seeds": bi["boot_n_seeds"], "pre_registered_bound": bound,
            "protocol_sha256": proto.sha256})

    ov = cluster_value_interval(overall, np.median, proto, "p6", "overall", PRIMARY_TAG)
    out_rows.append({
        "scope": "overall_ratio", "tag": PRIMARY_TAG, "case": "B0-B4",
        "method": "CANDIDATE/LEGACY",
        "ratio_candidate_over_legacy": ov["boot_estimate"], "ratio_lo": ov["boot_lo"],
        "ratio_hi": ov["boot_hi"], "ratio_conf": ov["boot_conf"],
        "boot_kind": ov["boot_kind"], "boot_resamples": ov["boot_resamples"],
        "boot_rng_seed": ov["boot_rng_seed"], "boot_n_seeds": ov["boot_n_seeds"],
        "pre_registered_bound": bound, "protocol_sha256": proto.sha256})

    repeat_ok = all(bool(r.get("repeat_identical")) for r in repro)
    for r in sorted(repro, key=lambda r: (r.get("tag"), r.get("case"), r.get("method"))):
        out_rows.append({"scope": "reproducibility", "tag": r.get("tag"),
                         "stdlib": r.get("stdlib"), "arch": r.get("arch"),
                         "case": r.get("case"), "method": r.get("method"),
                         "precision": r.get("precision"), "kappa": jnum(r, "kappa"),
                         "cap": r.get("cap"), "seed": jint(r, "seed"),
                         "repeat_identical": bool(r.get("repeat_identical")),
                         "protocol_sha256": proto.sha256})

    write_csv(ctx.out("performance.csv"), cols, out_rows)
    return {"rows": out_rows, "ratio": ov["boot_estimate"], "ratio_hi": ov["boot_hi"],
            "ratio_lo": ov["boot_lo"], "repeat_identical": repeat_ok, "bound": bound}


# ---------------------------------------------------------------------------
# G0 and G6 evidence
# ---------------------------------------------------------------------------
def validate_log_q(ctx: Context) -> tuple[float, list[dict]]:
    shapes = sorted({round(float(k) - 0.5, 12) for k in ctx.proto.get("matrix",
                                                                      "kappa_ladder")})
    rows = S.validate_log_q_evaluator(shapes)
    cols = ["shape_a", "log_w_min", "log_w_max", "points", "max_abs_log_q_error",
            "worst_at_log_w", "tolerance", "pass", "oracle", "protocol_sha256"]
    for r in rows:
        r["oracle"] = "mpmath.betainc at 60 decimal digits"
        r["protocol_sha256"] = ctx.proto.sha256
    write_csv(ctx.out("log_q_evaluator_validation.csv"), cols, rows)
    worst = max((r["max_abs_log_q_error"] for r in rows), default=float("nan"))
    return float(worst), rows


def verify_checksum_manifest(path: str, base: str) -> int | None:
    """Recompute every hash a checksum manifest lists.

    The gate asks whether the manifests VERIFY, which is a computation over file contents;
    Experiment 6's G6 asked whether the file existed.
    """
    if not os.path.exists(path):
        return None
    bad = 0
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            digest, _, name = line.partition("  ")
            target = os.path.join(base, name)
            if not os.path.exists(target) or sha256_file(target) != digest.strip():
                bad += 1
    return 0 if bad == 0 else bad


def g6_evidence(ctx: Context, p1: dict, override_path: str | None) -> dict:
    ev: dict = {"protocol_sha256": ctx.proto.sha256}

    raw_manifest = os.path.join(HERE, "raw", "raw_checksums.sha256")
    top_manifest = os.path.join(HERE, "checksums.sha256")
    a = verify_checksum_manifest(raw_manifest, os.path.join(HERE, "raw"))
    b = verify_checksum_manifest(top_manifest, HERE)
    ev["make_verify_exit_code"] = None if (a is None or b is None) else int(
        0 if (a == 0 and b == 0) else 1)

    envjson = os.path.join(HERE, "raw", "environment.json")
    if os.path.exists(envjson):
        with open(envjson, encoding="utf-8") as fh:
            env = json.load(fh)
        ev["dependencies_clean"] = bool(env.get("git", {}).get("dependencies_clean"))
        ev["archive_identifier"] = env.get("archive_identifier")
    else:
        ev["dependencies_clean"] = None
        ev["archive_identifier"] = None

    # "The P1 baseline comparison resolved": every P1 configuration carries both the previous
    # release (LEGACY, 1.0.0) and the candidate, on the same seed and the same stream.
    have = {(r["tag"], r["method"], r["precision"], round(r["kappa"], 12), r["seed"])
            for r in p1["scalar_rows"]}
    unmatched = [k for k in have if k[1] == "CANDIDATE"
                 and (k[0], "LEGACY", k[2], k[3], k[4]) not in have]
    ev["baseline_comparison_resolved"] = bool(have) and not unmatched

    # `make reverify` runs after this script, so the one thing this script cannot observe is
    # whether its own output regenerated byte for byte.  It is read from a receipt the
    # release step writes; absent, the key stays None and G6 fails naming it, which is the
    # honest state of a run that has not been re-verified yet.
    ev["make_reverify_identical"] = None
    path = override_path or ctx.out("g6_evidence.json")
    if os.path.exists(path):
        with open(path, encoding="utf-8") as fh:
            supplied = json.load(fh)
        for k in ("make_verify_exit_code", "make_reverify_identical", "dependencies_clean",
                  "baseline_comparison_resolved", "archive_identifier"):
            if k in supplied:
                ev[k] = supplied[k]
    return ev


def rng_stream_break_documented() -> bool:
    """Computed, not asserted: the stream break is stated in the frozen protocol and in the
    released header's own version history."""
    needles = ("The random stream changes for every seed",
               "No seed-for-seed continuity")
    try:
        with open(os.path.join(HERE, "PROTOCOL.md"), encoding="utf-8") as fh:
            text = fh.read()
    except OSError:
        return False
    if not all(n in text for n in needles):
        return False
    header = os.path.join(ROOT, "cpp", "bi_kappa_distribution.H")
    try:
        with open(header, encoding="utf-8") as fh:
            htext = fh.read()
    except OSError:
        return False
    return "2.0.0" in htext


# ---------------------------------------------------------------------------
# Portability ingest (the CI contract)
# ---------------------------------------------------------------------------
def run_portability_ingest(bundle: str, proto: F.Protocol) -> int:
    """``analyze.py --portability-ingest <DIR>``.

    The contract is stated in .github/workflows/portability.yml, job ``g4-decision``: group by
    the environment record's architecture, require bitwise equality across standard libraries
    within an architecture, run the two one-sided tests on the log rate ratio across
    architectures, ignore every record whose execution is not native, write
    results/portability.csv and a machine-readable verdict, and exit non-zero unless G4
    closes.  None of those rules is reimplemented here: ``exp7_portability.ingest`` applies
    them and ``exp7_families.family_F7`` decides.
    """
    if not os.path.isdir(bundle):
        print(f"analyze.py: the portability bundle {bundle} is not a directory",
              file=sys.stderr)
        return 2
    res = PORT.ingest(bundle, proto)
    fam = res.pop("_family", None)
    out_dir = os.path.join(HERE, "results")
    os.makedirs(out_dir, exist_ok=True)

    rows = []
    for e in res.get("environments", []):
        rows.append({"kind": "environment",
                     "label": f"{e.get('arch')}/{e.get('stdlib')}", "arch": e.get("arch"),
                     "os": e.get("os"), "stdlib": e.get("stdlib"),
                     "compiler": e.get("compiler"), "tag": None,
                     "execution": e.get("execution"),
                     "available": bool(e.get("available")),
                     "completed": bool(e.get("completed")),
                     "reason_not_run": e.get("reason_not_run"), "n_rows": e.get("n_rows"),
                     "prediction": None,
                     "in_decision": e.get("execution") == "native" and bool(
                         e.get("completed")),
                     "protocol_sha256": proto.sha256})
    for p in res.get("stdlib_pairs", []):
        rows.append({"kind": "cross_stdlib_bitwise", "label": p.get("label"),
                     "available": True, "completed": True,
                     "reason_not_run": None,
                     "prediction": "bitwise equality", "n_common": p.get("n_common"),
                     "n_digests_compared": p.get("n_digests_compared"),
                     "n_differences": p.get("n_differences"),
                     "identical": bool(p.get("identical")), "in_decision": True,
                     "protocol_sha256": proto.sha256})
    tost_rows = {r.get("label"): r for r in (fam.rows if fam else []) or []
                 if r.get("kind") == "cross_arch_equivalence"}
    for p in res.get("arch_pairs", []):
        t = tost_rows.get(p["label"], {})
        rows.append({"kind": "cross_arch_equivalence", "label": p["label"],
                     "available": True, "completed": True, "reason_not_run": t.get("reason"),
                     "prediction": "equality of rates within the pre-registered margin",
                     "log_ratio": t.get("log_ratio"), "se": t.get("se"),
                     "margin": t.get("margin"), "p_value": t.get("p_value"),
                     "equivalent": t.get("equivalent"), "disagrees": t.get("disagrees"),
                     "informative": t.get("informative"), "in_decision": True,
                     "protocol_sha256": proto.sha256})
    write_csv(os.path.join(out_dir, "portability.csv"), PORT_COLUMNS, rows)

    archs = sorted({e.get("arch") for e in res.get("native", [])}, key=str)
    reasons = []
    if not (res.get("F7") or {}).get("passed"):
        reasons.append("F7 found a disagreement between environments")
    if res.get("pending"):
        reasons.append(f"{len(res['pending'])} expected environment(s) produced no rows")
    if len(archs) < 2:
        reasons.append(f"only {archs} has native results, so the cross-architecture claim "
                       "is untested")
    # family_F7 classifies a cross-architecture cell with too few events as UNDERPOWERED:
    # it neither establishes equivalence nor demonstrates a difference.  Neither
    # exp7_portability.main nor exp7_gates.gate_G4 acts on that classification, and this
    # entry point keeps parity with them rather than inventing a stricter rule -- but the
    # count is recorded here, because a gate that closed on cells that proved nothing is a
    # thing a reader has to be able to see.
    tost = [r for r in (fam.rows if fam else []) or []
            if r.get("kind") == "cross_arch_equivalence"]
    underpowered = [r.get("label") for r in tost if r.get("informative") is False]
    verdict = {"gate": "G4", "closed": not reasons, "reasons": reasons,
               "native_architectures": archs,
               "translated_excluded": len(res.get("translated", [])),
               "cross_arch_cells": len(tost),
               "cross_arch_cells_underpowered": len(underpowered),
               "underpowered_labels": sorted(underpowered),
               "protocol_sha256": proto.sha256, "F7": res.get("F7")}
    with open(os.path.join(out_dir, "g4_verdict.json"), "w", encoding="utf-8") as fh:
        json.dump(verdict, fh, indent=2, sort_keys=True, default=str)
        fh.write("\n")
    print(json.dumps(verdict, indent=2, sort_keys=True, default=str))
    if reasons:
        print("G4 is NOT closed: " + "; ".join(reasons), file=sys.stderr)
        print("The exact command for every environment is in "
              "results/portability_remote.md.", file=sys.stderr)
        return 1
    print("G4 closes on this evidence.", file=sys.stderr)
    return 0


# ---------------------------------------------------------------------------
# Provenance
# ---------------------------------------------------------------------------
def write_provenance(ctx: Context, oracle_audit: list[dict]) -> None:
    """The one file here that records wall-clock time and build identity.

    Every value in it is copied out of the recorded run -- the manifest's `run_utc`, the
    probe's own hash of itself, the compile line baked into the binary -- rather than read
    from the clock or the host at analysis time, so this file too regenerates byte for byte
    and `make reverify` stays possible.  No absolute path appears anywhere in it.
    """
    m = ctx.manifest
    L = ["# Experiment 7 - provenance of this analysis", "",
         "This is the only file under results/ that carries wall-clock times and build "
         "identity. Every value below is copied from raw/manifest.csv and the environment "
         "record written at run time; nothing is read from the clock or the host while the "
         "analysis runs, so this file regenerates byte for byte and `make reverify` remains "
         "a meaningful check.", "",
         f"- protocol: `config/protocol.json` version "
         f"{ctx.proto.get('protocol_version')}, SHA-256 `{ctx.proto.sha256}`",
         f"- mode: {'smoke' if ctx.smoke else 'production'}",
         f"- manifest: `{ctx.rel(ctx.manifest_path)}`, "
         f"{len({r.get('file') for r in m})} distinct output files recorded in {len(m)} "
         f"manifest rows", ""]
    runs = sorted({r.get("run_utc", "") for r in m})
    L += [f"- probe run times (UTC, from the writing processes): {runs[0]} .. {runs[-1]}"
          if runs else "- probe run times: none recorded", ""]
    L += ["## Builds that produced the data", "",
          "| tag | compiler | stdlib | arch | execution | probe SHA-256 | cxxflags |",
          "|---|---|---|---|---|---|---|"]
    seen = set()
    for r in sorted(m, key=lambda r: (r.get("env_tag", ""), r.get("probe_sha256", ""))):
        key = (r.get("env_tag"), r.get("probe_sha256"))
        if key in seen:
            continue
        seen.add(key)
        L.append(f"| {r.get('env_tag')} | {r.get('compiler')} | {r.get('stdlib')} | "
                 f"{r.get('arch')} | {r.get('execution')} | `{r.get('probe_sha256')}` | "
                 f"`{r.get('cxxflags')}` |")
    L += ["", "## Adjudication", ""]
    if oracle_audit:
        files = sorted({r.get("file", "") for r in oracle_audit})
        L += [f"- {len(files)} audit stream(s) adjudicated by "
              f"`{oracle_audit[0].get('oracle')}`:"]
        L += [f"  - `{f}`" for f in files]
    else:
        L += ["- none recorded"]
    L += ["", "## Analysis inputs", "",
          "| file | SHA-256 |", "|---|---|"]
    for name in ("PROTOCOL.md", "config/protocol.json", "config/power_study.json",
                 "config/honest_floor.md", "results/schema.md", "analyze.py",
                 "exp7_families.py", "exp7_gates.py", "exp7_stats.py", "exp7_io.py",
                 "exp7_portability.py"):
        p = os.path.join(HERE, name)
        if os.path.exists(p):
            L.append(f"| `{name}` | `{sha256_file(p)}` |")
    L.append("")
    with open(ctx.out("provenance.md"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(L))


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------
def to_jsonable(obj):
    if isinstance(obj, dict):
        return {str(k): to_jsonable(v) for k, v in sorted(obj.items(), key=lambda kv: str(
            kv[0]))}
    if isinstance(obj, (list, tuple)):
        return [to_jsonable(v) for v in obj]
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        v = float(obj)
        return v if math.isfinite(v) else str(v)
    if isinstance(obj, float):
        return obj if math.isfinite(obj) else str(obj)
    if isinstance(obj, (np.bool_,)):
        return bool(obj)
    if isinstance(obj, GATES.Gate):
        return obj.as_dict()
    if isinstance(obj, F.FamilyResult):
        return obj.as_dict()
    return obj


def run(ctx: Context, g6_override: str | None) -> int:
    proto = ctx.proto
    ctx.check_protocol_hash()
    print(f"  protocol {proto.get('protocol_version')} sha256 {proto.sha256[:12]} matches "
          f"the run manifest "
          f"({len({r.get('file') for r in ctx.manifest})} recorded files)")

    # A failed run must not leave a verdict standing beside partially regenerated source
    # data.  The two files that carry one are removed before anything is written and are
    # rewritten only at the end, so an analysis that stops on a missing phase leaves no
    # verdict anywhere rather than a stale one.
    for name in ("analysis_report.md", "exp7_results.json"):
        stale = os.path.join(ctx.out_dir, name)
        if os.path.exists(stale):
            os.remove(stale)

    oracle_audit, oracle_disag = ctx.oracle_rows()

    p1 = analyse_p1(ctx)
    print(f"  P1: {len(p1['scalar_rows'])} scalar rows")
    p2 = analyse_p2(ctx, p1)
    print(f"  P2: {len(p2['rows'])} failure-envelope rows")
    p3 = analyse_p3(ctx, p2)
    print(f"  P3: {len(p3['bin_rows'])} conditioning bins, {len(p3['tail_rows'])} tail rows")
    p4 = analyse_p4(ctx)
    print(f"  P4: {len(p4['rows'])} loader rows")
    nc = analyse_negative_controls(ctx, p1, p4)
    print(f"  negative controls: {len(nc['rows'])} measured power rows")
    p5 = analyse_p5(ctx, p1)
    print(f"  P5: {len(p5['rows'])} portability rows")
    p6 = analyse_p6(ctx)
    print(f"  P6: {len(p6['rows'])} performance rows")

    worst_log_q, _ = validate_log_q(ctx)

    families = dict(p1["families"])
    families["F5_loader_battery"] = p4["F5"]
    families["F6_negative_controls"] = nc["F6"]
    families["F7_portability"] = p5["F7"]

    disag_records = [r for r in oracle_disag if r.get("kind") == "disagreement"]
    headers = [r for r in oracle_disag if r.get("kind") == "header"]
    if not headers:
        raise AnalysisError(
            "raw/oracle_disagreements.jsonl carries no `header` record, so an empty "
            "adjudication cannot be told from an oracle that never ran; PROTOCOL.md 2.3 "
            "requires one per audited file")
    # The oracle reports a per-file disagreement count and also writes one record per
    # disagreement; the larger of the two is used, so neither a missing record nor a
    # missing count can hide one.
    oracle_disagreements = max(len(disag_records),
                               max((int(r.get("disagreements", 0)) for r in oracle_audit),
                                   default=0))
    conversion_failures = max((int(r.get("conversion_failures", 0)) for r in oracle_audit),
                              default=0)
    audited_records = sum(int(r.get("n_records", 0)) for r in oracle_audit)

    candidate_avoidable = (p1["candidate_avoidable_p1"] + p2["candidate_avoidable_p2"]
                           + p3["candidate_avoidable_p3"])
    benign_ok = (bool(p1["benign"]) and all(all(v) for v in p1["benign"].values())
                 and p1["benign_digest_agrees"])

    evidence = {
        "log_q_worst_abs_error": worst_log_q,
        "candidate_avoidable_total": int(candidate_avoidable),
        "benign_control_exercises_candidate": benign_ok,
        "accounting_failures": int(p2["accounting_failures"] + p3["accounting_failures"]),
        "oracle_disagreements": int(oracle_disagreements),
        "oracle_conversion_failures": int(conversion_failures),
        "oracle_audited_records": int(audited_records),
        "honest_minus_floor_unadjudicated": int(
            max(0, p2["honest_minus_ref_total"] + p3["honest_minus_ref_total"]
                - len(disag_records))),
        "audit_coverage": p3["audit_coverage"],
        "conditional_cells_without_loss_fraction": int(p4["unlabelled_conditional"]),
        "environments": p5["environments"],
        "time_per_return_ratio": p6["ratio"],
        "time_per_return_ratio_hi": p6["ratio_hi"],
        "same_seed_reproducible": bool(p6["repeat_identical"]),
        "rng_stream_break_documented": rng_stream_break_documented(),
    }
    evidence.update(g6_evidence(ctx, p1, g6_override))

    result = GATES.evaluate(proto, families, evidence)
    GATES.write_report(
        ctx.out("analysis_report.md"), proto, result, families, evidence,
        "deterministically from the recorded run; this report carries no wall-clock "
        "stamp, so `make reverify` can diff it. Run times and build identity are in "
        "`provenance.md`")

    payload = {
        "experiment": "exp7_confirmatory",
        "mode": "smoke" if ctx.smoke else "production",
        "protocol_version": proto.get("protocol_version"),
        "protocol_sha256": proto.sha256,
        "verdict": result["verdict"],
        "gates": {k: g.as_dict() for k, g in result["gates"].items()},
        "families": {k: v.as_dict() for k, v in families.items()},
        "evidence": evidence,
        "expected_false_failure_bound": proto.get("statistics",
                                                  "expected_false_failure_bound"),
        "F2": {"observed_misses": p1["f2_misses"], "resolved_intervals": p1["f2_resolved"],
               "informative_intervals": sum(1 for c in p1["f2_cells"] if c["informative"]),
               "total_intervals": len(p1["f2_cells"])},
        "F5_tests": p4["f5_tests"],
        "F5_test_names": p4["f5_test_names"],
        "F5_tail_q0": p4["q0s"],
        "F6_controls": nc["controls"],
        # The amendment-1.3.0 corrections, copied from the protocol into the result so that
        # a reader of this file can see which reading of each rule the numbers were produced
        # under without opening a second document.
        "corrections_applied": proto.get("corrections"),
        "amendments": [a.get("version") for a in proto.get("amendments")],
        "performance": {"ratio": p6["ratio"], "lo": p6["ratio_lo"], "hi": p6["ratio_hi"],
                        "bound": p6["bound"]},
        "cross_phase_native_disagreements": p2["cross_phase_disagreements"],
        "primary_environment_tag": PRIMARY_TAG,
    }
    with open(ctx.out("exp7_results.json"), "w", encoding="utf-8") as fh:
        json.dump(to_jsonable(payload), fh, indent=2, sort_keys=True)
        fh.write("\n")

    write_provenance(ctx, oracle_audit)

    print(f"\n{result['verdict']}")
    for k in sorted(result["gates"]):
        g = result["gates"][k]
        print(f"  {k} {GATES.GATE_NAMES[k]:<26} {g.status}")
    print(f"\n  wrote {ctx.rel(ctx.out_dir)}/ "
          f"(analysis_report.md, exp7_results.json and the source-data CSVs)")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--smoke", action="store_true",
                    help="analyse raw/smoke/ only and write results/smoke/ only")
    ap.add_argument("--portability-ingest", metavar="DIR", default=None,
                    help="ingest a CI portability bundle, decide F7/G4, and exit non-zero "
                         "unless the gate closes")
    ap.add_argument("--g6-evidence", metavar="PATH", default=None,
                    help="JSON receipt written by the release step, carrying the `make "
                         "verify` and `make reverify` outcomes this script cannot observe")
    args = ap.parse_args()
    proto = F.Protocol()
    if args.portability_ingest:
        return run_portability_ingest(args.portability_ingest, proto)
    ctx = Context(smoke=args.smoke, proto=proto)
    try:
        return run(ctx, args.g6_evidence)
    except AnalysisError as exc:
        print(f"\nanalyze.py: {exc}", file=sys.stderr)
        print("No verdict was computed and no result file carries one.", file=sys.stderr)
        return 3


if __name__ == "__main__":
    sys.exit(main())
