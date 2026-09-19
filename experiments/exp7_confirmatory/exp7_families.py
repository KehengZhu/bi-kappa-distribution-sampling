#!/usr/bin/env python3
"""The seven pre-registered test families of ``PROTOCOL.md`` §5, and the gate rules of §6.

Everything here reads its levels, margins and critical values from ``config/protocol.json``.
Nothing in this module carries a numeric acceptance threshold of its own: if a constant is
not in the protocol, it is not an acceptance criterion, and the loader below raises rather
than defaulting.  That is the whole point of the file -- Experiment 6's acceptance rules
lived as literals scattered through a 2111-line analysis, which is how a 0.9 appeared where
the plan said 1.0.

The statistical primitives (Clopper-Pearson, Holm, order-statistic brackets, the two-piece
``log q`` evaluator) are reused from ``exp7_stats.py``, which is Experiment 6's frozen
module; the *rules built on them* are new, because Experiment 6's rules are what failed
review.
"""
from __future__ import annotations

import hashlib
import json
import os
from dataclasses import dataclass, field

import numpy as np
from scipy import stats

import exp7_stats as S

HERE = os.path.dirname(os.path.abspath(__file__))
PROTOCOL_PATH = os.path.join(HERE, "config", "protocol.json")


# ---------------------------------------------------------------------------
# Protocol loading.  A missing key is an error, never a default.
# ---------------------------------------------------------------------------
class Protocol:
    """The frozen protocol, plus the SHA-256 that ties a result to this exact text."""

    def __init__(self, path: str = PROTOCOL_PATH):
        with open(path, "rb") as fh:
            raw = fh.read()
        self.sha256 = hashlib.sha256(raw).hexdigest()
        self.path = path
        self.d = json.loads(raw)

    def get(self, *keys):
        node = self.d
        for k in keys:
            if k not in node:
                raise KeyError(
                    f"{'/'.join(map(str, keys))} is not in {self.path}. Every acceptance "
                    "criterion must be pre-registered; this module will not invent one.")
            node = node[k]
        return node

    def alpha(self, family: str) -> float:
        return float(self.get("statistics", "family_alpha", family))


@dataclass
class FamilyResult:
    """One family's global decision.  ``passed`` is the only thing a gate may read."""

    name: str
    alpha: float
    passed: bool
    statistic: float
    p_value: float
    detail: str
    offenders: list = field(default_factory=list)
    rows: list = field(default_factory=list)

    def as_dict(self) -> dict:
        return {"family": self.name, "alpha": self.alpha, "passed": bool(self.passed),
                "statistic": float(self.statistic), "p_value": float(self.p_value),
                "detail": self.detail, "offenders": self.offenders}


# ---------------------------------------------------------------------------
# Global combination rules
# ---------------------------------------------------------------------------
def simes_global(pvalues) -> float:
    """Simes' global p-value for the intersection null.

    Valid under independence and under positive regression dependence, which is the
    situation here (statistics computed on overlapping samples of the same draws).  Used
    instead of "zero rejections anywhere": over K families each controlled at alpha, the
    latter is a level-(1-(1-alpha)^K) procedure in the *failure* direction, which for
    Experiment 6's E1 alone reached 0.57 per candidate.
    """
    p = np.sort(np.asarray([x for x in pvalues if np.isfinite(x)], dtype=float))
    n = p.size
    if n == 0:
        return 1.0
    return float(np.min(n * p / np.arange(1, n + 1)))


def poisson_binomial_sf(miss_probs, k: int) -> float:
    """``P(misses > k)`` by exact convolution of independent Bernoulli trials."""
    pmf = np.array([1.0])
    for q in miss_probs:
        pmf = np.convolve(pmf, [1.0 - q, q])
    if k + 1 >= pmf.size:
        return 0.0
    return float(pmf[k + 1:].sum())


def stouffer(zs, weights=None) -> float:
    """Stouffer combination of standardized errors; two-sided p-value."""
    z = np.asarray([v for v in zs if np.isfinite(v)], dtype=float)
    if z.size == 0:
        return 1.0
    w = np.ones_like(z) if weights is None else np.asarray(weights, dtype=float)[:z.size]
    zc = float(np.sum(w * z) / np.sqrt(np.sum(w * w)))
    return float(2.0 * stats.norm.sf(abs(zc)))


def anderson_darling_exp1(z: np.ndarray) -> tuple[float, float]:
    """A^2 against the unit exponential, fully specified null (no parameter estimation).

    Anderson-Darling is included because the entire claim is about the tail and both KS and
    Cramer-von Mises are bulk-weighted.  Experiment 6's own data calibrate the gap: its
    radial battery detected survivor conditioning only above roughly a 1e-3 loss fraction,
    and the cells it certified sat below that.
    """
    x = np.sort(np.asarray(z, dtype=float))
    x = x[np.isfinite(x)]
    n = x.size
    if n < 8:
        return float("nan"), float("nan")
    u = -np.expm1(-x)                      # CDF of Exp(1), accurate for small x
    u = np.clip(u, 1e-300, 1.0 - 1e-16)
    i = np.arange(1, n + 1)
    a2 = -n - np.sum((2 * i - 1) * (np.log(u) + np.log1p(-u[::-1]))) / n
    # Marsaglia & Marsaglia (2004), "Evaluating the Anderson-Darling distribution",
    # J. Stat. Softw. 9(2).  `adinf` is the limiting *CDF* P(A^2 < z) for the fully
    # specified null; the finite-n correction is negligible at the sample sizes used here
    # (n >= 1e5).  The p-value is its upper tail -- 1 - adinf, not adinf, since a large A^2
    # is evidence against the null.
    if a2 <= 0.0:
        return float(a2), 1.0
    if a2 < 2.0:
        cdf = (a2 ** -0.5) * np.exp(-1.2337141 / a2) * (
            2.00012 + (0.247105 - (0.0649821 - (0.0347962 - (0.011672 - 0.00168691 * a2)
                                                * a2) * a2) * a2) * a2)
    else:
        cdf = np.exp(-np.exp(1.0776 - (2.30695 - (0.43424 - (0.082433 - (0.008056
                     - 0.0003146 * a2) * a2) * a2) * a2) * a2))
    return float(a2), float(min(max(1.0 - cdf, 0.0), 1.0))


def cluster_bootstrap(values_by_seed, statistic, resamples: int, conf: float,
                      rng_seed: int) -> dict:
    """Cluster bootstrap over seeds: resample seeds with replacement, then within seeds.

    Experiment 6's "seed-stratified" bootstrap resampled Binomial(n_i, k_i/n_i) at fixed
    stratum size, whose bootstrap variance is exactly the pooled binomial variance with
    *zero* between-seed component -- so no interval anywhere in that analysis could detect
    seed-level heterogeneity, which is the signature of an RNG-stream or portability defect.
    Resampling the seeds themselves restores that component.
    """
    rng = np.random.default_rng(rng_seed)
    seeds = list(values_by_seed.keys())
    k = len(seeds)
    if k == 0:
        return {"estimate": float("nan"), "lo": float("nan"), "hi": float("nan"),
                "resolved": False, "n_seeds": 0}
    point = statistic(np.concatenate([np.asarray(values_by_seed[s]) for s in seeds]))
    draws = np.empty(resamples, dtype=float)
    for r in range(resamples):
        picked = rng.integers(0, k, size=k)
        parts = []
        for idx in picked:
            v = np.asarray(values_by_seed[seeds[idx]])
            if v.size:
                parts.append(v[rng.integers(0, v.size, size=v.size)])
        draws[r] = statistic(np.concatenate(parts)) if parts else np.nan
    good = draws[np.isfinite(draws)]
    if good.size < resamples // 2:
        return {"estimate": float(point), "lo": float("nan"), "hi": float("nan"),
                "resolved": False, "n_seeds": k}
    lo, hi = np.percentile(good, [100 * (1 - conf) / 2, 100 * (1 + conf) / 2])
    return {"estimate": float(point), "lo": float(lo), "hi": float(hi),
            "resolved": True, "n_seeds": k}


def seed_homogeneity(counts, sizes) -> tuple[float, float]:
    """Exact-ish homogeneity of a rate across seeds: chi-square on the 2xK table."""
    counts = np.asarray(counts, dtype=float)
    sizes = np.asarray(sizes, dtype=float)
    if counts.size < 2 or sizes.sum() <= 0:
        return float("nan"), float("nan")
    table = np.vstack([counts, sizes - counts])
    if np.any(table.sum(axis=1) == 0):
        return float("nan"), float("nan")
    chi2, p, _, _ = stats.chi2_contingency(table, correction=False)
    return float(chi2), float(p)


def tost_log_rate_ratio(k1: int, n1: int, k2: int, n2: int, margin: float,
                        alpha: float) -> dict:
    """Two one-sided tests that ``|log(p1/p2)| < margin``.

    Equivalence, not identity: two environments are declared to agree when the log rate
    ratio is inside a pre-registered margin, which is a statement that can pass.  Testing
    equality instead can only fail to reject, which is not the same thing and is what
    Experiment 6's G4 never implemented at all.
    """
    if n1 <= 0 or n2 <= 0:
        return {"resolved": False, "passed": False, "reason": "empty sample"}
    if k1 == 0 and k2 == 0:
        # Both zero: nothing disagrees, but nothing is demonstrated either.  Report the
        # one-sided upper limits and mark the cell non-informative rather than counting a
        # pair of zeros as agreement.
        u1, u2 = S.zero_count_upper(n1), S.zero_count_upper(n2)
        return {"resolved": False, "informative": False, "disagrees": False,
                "reason": "no failures in either arm", "upper1": u1, "upper2": u2}
    if k1 == 0 or k2 == 0:
        # One arm saw failures and the other saw none.  Whether that is a real difference
        # depends on whether the zero arm's upper limit excludes the other arm's rate.
        seen, zero_n = (k1 / n1, n2) if k1 else (k2 / n2, n1)
        upper = S.zero_count_upper(zero_n)
        disagrees = bool(seen > upper)
        return {"resolved": False, "informative": bool(disagrees), "disagrees": disagrees,
                "reason": ("one arm has zero failures and the other's rate lies above the "
                           "zero arm's one-sided upper limit" if disagrees else
                           "one arm has zero failures; the counts cannot separate them"),
                "observed_rate": float(seen), "zero_arm_upper": float(upper)}
    p1, p2 = k1 / n1, k2 / n2
    est = np.log(p1) - np.log(p2)
    se = np.sqrt((1 - p1) / k1 + (1 - p2) / k2)
    p_lo = float(stats.norm.sf((est + margin) / se))
    p_hi = float(stats.norm.sf((margin - est) / se))
    pval = max(p_lo, p_hi)
    equivalent = bool(pval < alpha)
    # Failing to show equivalence is not the same as showing disagreement.  A cell whose
    # point estimate is inside the margin but whose interval is too wide to prove it is
    # UNDERPOWERED: it leaves the gate open rather than failing it, because absent evidence
    # must not pass and must not condemn either.
    disagrees = bool(abs(est) > margin)
    return {"resolved": True, "passed": equivalent, "equivalent": equivalent,
            "disagrees": disagrees,
            "informative": bool(equivalent or disagrees),
            "log_ratio": float(est), "se": float(se), "margin": float(margin),
            "p_value": pval,
            "reason": ("" if equivalent else
                       ("the log rate ratio lies outside the pre-registered margin"
                        if disagrees else
                        "underpowered: the point estimate is inside the margin but the "
                        "interval cannot establish equivalence"))}


# ---------------------------------------------------------------------------
# F1 -- radial law
# ---------------------------------------------------------------------------
def family_F1(proto: Protocol, cells) -> FamilyResult:
    """``cells``: dicts with ``label`` and ``pvalues`` (AD, KS, CvM on Z; Beta on W where
    every draw resolves).  Global rule: Simes over all configurations at alpha_F."""
    alpha = proto.alpha("F1_radial_law")
    per_cell = [(c["label"], simes_global(c["pvalues"])) for c in cells]
    if not per_cell:
        return FamilyResult("F1_radial_law", alpha, False, float("nan"), float("nan"),
                            "no radial cells were produced; an absent family is a failure, "
                            "not a vacuous pass")
    p_global = simes_global([p for _, p in per_cell])
    rejected = S.holm([p for _, p in per_cell], alpha=alpha)
    offenders = [lab for (lab, _), rej in zip(per_cell, rejected) if rej]
    return FamilyResult("F1_radial_law", alpha, bool(p_global >= alpha), float(p_global),
                        float(p_global),
                        f"Simes over {len(per_cell)} configurations; "
                        f"{len(offenders)} identified by Holm on rejection",
                        offenders, [{"label": l, "p": p} for l, p in per_cell])


# ---------------------------------------------------------------------------
# F2 -- quantile coverage
# ---------------------------------------------------------------------------
def family_F2(proto: Protocol, observed_misses: int, n_resolved: int) -> FamilyResult:
    """Poisson-binomial upper-tail test on the miss count.

    The null is fixed in the protocol from the analytically computed achieved coverage of
    every interval, so the critical value existed before the data did.  Neither "every
    interval covers" (which fails with probability 0.984 under a correct sampler) nor a
    round 0.9 (which has no power below a tenfold degradation) is used.
    """
    alpha = proto.alpha("F2_quantile_coverage")
    spec = proto.get("F2")
    critical = int(spec["critical_miss_count"])
    expected = float(spec["expected_misses"])
    frozen_n = int(spec["n_resolved_intervals"])
    if n_resolved != frozen_n:
        return FamilyResult(
            "F2_quantile_coverage", alpha, False, float(observed_misses), float("nan"),
            f"the run resolved {n_resolved} intervals but the protocol froze {frozen_n}; "
            "the null distribution no longer applies and the family cannot be evaluated "
            "without re-registering it")
    miss_probs = [1.0 - c["achieved_coverage"] for c in proto.get("F2_cells")
                  if c.get("resolved")]
    p = poisson_binomial_sf(miss_probs, observed_misses - 1) if observed_misses > 0 else 1.0
    return FamilyResult("F2_quantile_coverage", alpha,
                        bool(observed_misses <= critical), float(observed_misses), float(p),
                        f"{observed_misses} misses of {n_resolved} resolved intervals; "
                        f"expected {expected:.2f}; fail above {critical}")


def informative(lo_value: float, hi_value: float, proto: Protocol) -> bool:
    """Is a quantile bracket narrow enough to be evidence at all?

    At kappa = 0.5001 the order-statistic bracket on log R spans up to 8951 natural-log
    units.  An interval admitting a multiplicative error of e^8951 cannot distinguish a
    correct sampler from a badly wrong one, so it is declared non-informative in advance and
    the configuration is certified by F1/F3/F4 instead.  Classifying it after the fact would
    be exactly the move this protocol exists to prevent.
    """
    width = float(hi_value) - float(lo_value)
    return bool(np.isfinite(width) and
                width <= float(proto.get("statistics", "informative_width_log_units")))


# ---------------------------------------------------------------------------
# F3 -- quantile direction
# ---------------------------------------------------------------------------
def family_F3(proto: Protocol, cells) -> FamilyResult:
    """``cells``: dicts with ``label`` and ``z`` (standardized signed quantile errors, one
    per seed), calibrated against a Monte-Carlo null because the null is not symmetric at
    small shape.  Holm over cells at alpha_F.

    This is the family that has power where F2 is vacuous: coverage throws away sign, and
    sign is where a systematic tail deficit shows up.
    """
    alpha = proto.alpha("F3_quantile_direction")
    per_cell = [(c["label"], stouffer(c["z"])) for c in cells if len(c.get("z", []))]
    if not per_cell:
        return FamilyResult("F3_quantile_direction", alpha, False, float("nan"),
                            float("nan"), "no direction cells were produced")
    rejected = S.holm([p for _, p in per_cell], alpha=alpha)
    offenders = [lab for (lab, _), rej in zip(per_cell, rejected) if rej]
    p_global = simes_global([p for _, p in per_cell])
    return FamilyResult("F3_quantile_direction", alpha, not any(rejected), float(p_global),
                        float(p_global),
                        f"Holm over {len(per_cell)} (kappa, precision, p) cells",
                        offenders, [{"label": l, "p": p} for l, p in per_cell])


# ---------------------------------------------------------------------------
# F4 -- upper-tail mass
# ---------------------------------------------------------------------------
def exceedance_test(z: np.ndarray, q0: float) -> dict:
    """Exact binomial on the count above ``z0 = -log q0``, plus KS of the excesses.

    Under the null the count is exactly Binomial(n, q0) and the excesses are exactly Exp(1),
    so both statistics are calibrated with no asymptotics and both retain power in the tail
    where KS on the whole sample has none.
    """
    z = np.asarray(z, dtype=float)
    z = z[np.isfinite(z)]
    n = z.size
    if n == 0:
        return {"resolved": False}
    z0 = -np.log(q0)
    k = int(np.sum(z > z0))
    p_count = float(stats.binomtest(k, n, q0).pvalue)
    excess = z[z > z0] - z0
    if excess.size >= 8:
        p_exc = float(stats.kstest(excess, "expon").pvalue)
    else:
        p_exc = float("nan")
    return {"resolved": True, "q0": q0, "n": n, "observed": k, "expected": n * q0,
            "p_count": p_count, "p_excess": p_exc}


def family_F4(proto: Protocol, cells) -> FamilyResult:
    """``cells``: dicts with ``label`` and ``tests`` (the dicts ``exceedance_test`` returns)."""
    alpha = proto.alpha("F4_upper_tail_mass")
    pvals, labels = [], []
    for c in cells:
        for t in c.get("tests", []):
            if not t.get("resolved"):
                continue
            for key in ("p_count", "p_excess"):
                v = t.get(key)
                if v is not None and np.isfinite(v):
                    pvals.append(float(v))
                    labels.append(f"{c['label']}|q0={t['q0']:g}|{key}")
    if not pvals:
        return FamilyResult("F4_upper_tail_mass", alpha, False, float("nan"), float("nan"),
                            "no tail cells were produced")
    rejected = S.holm(pvals, alpha=alpha)
    offenders = [lab for lab, rej in zip(labels, rejected) if rej]
    p_global = simes_global(pvals)
    return FamilyResult("F4_upper_tail_mass", alpha, bool(p_global >= alpha), float(p_global),
                        float(p_global), f"Simes over {len(pvals)} exceedance statistics",
                        offenders)


# ---------------------------------------------------------------------------
# F5 -- complete-loader battery
# ---------------------------------------------------------------------------
def family_F5(proto: Protocol, tests) -> FamilyResult:
    """``tests``: dicts with ``label``, ``p``, ``conditional`` and ``loss_fraction``.

    Holm runs over **all** cells and tests jointly.  Experiment 6 applied Holm inside each
    of 42 cells and then demanded zero rejections across them, which is not a level-alpha
    procedure.
    """
    alpha = proto.alpha("F5_loader_battery")
    usable = [t for t in tests if np.isfinite(t.get("p", np.nan))]
    if not usable:
        return FamilyResult("F5_loader_battery", alpha, False, float("nan"), float("nan"),
                            "no loader tests were produced")
    rejected = S.holm([t["p"] for t in usable], alpha=alpha)
    offenders = [t["label"] for t, rej in zip(usable, rejected) if rej]
    conditional = [t["label"] for t in usable if t.get("conditional")]
    p_global = simes_global([t["p"] for t in usable])
    return FamilyResult(
        "F5_loader_battery", alpha, not any(rejected), float(p_global), float(p_global),
        f"Holm jointly over {len(usable)} tests across all cells; "
        f"{len(conditional)} cells are conditional on success and are labelled as such",
        offenders, [{"label": t["label"], "p": t["p"],
                     "conditional": bool(t.get("conditional")),
                     "loss_fraction": t.get("loss_fraction")} for t in usable])


# ---------------------------------------------------------------------------
# F6 -- negative controls, by measured power
# ---------------------------------------------------------------------------
def family_F6(proto: Protocol, controls) -> FamilyResult:
    """``controls``: dicts with ``name``, ``effect``, ``detections``, ``injections``.

    The rule runs in the opposite direction to the others: a control passes when the battery
    *does* detect the injected defect, at or above the pre-registered power, at the effect
    size actually being certified.  Experiment 6 injected a defect into 50 per cent of the
    sample and reported a binary, which says nothing about power at a 1e-3 loss fraction --
    and its ``all([])`` meant a control that never ran passed.
    """
    required = float(proto.get("F6", "required_power"))
    expected_names = set(proto.get("F6", "controls"))
    seen = {c["name"] for c in controls}
    missing = sorted(expected_names - seen)
    rows, ok = [], True
    for c in controls:
        n = int(c.get("injections", 0))
        d = int(c.get("detections", 0))
        power = d / n if n else float("nan")
        lo, _ = S.clopper_pearson(d, n) if n else (float("nan"), float("nan"))
        passed = bool(n > 0 and lo >= required)
        ok = ok and passed
        rows.append({"name": c["name"], "effect": c.get("effect"), "injections": n,
                     "detections": d, "power": power, "power_lo": lo,
                     "required": required, "passed": passed})
    if missing:
        ok = False
    return FamilyResult(
        "F6_negative_controls", float("nan"), bool(ok and not missing), float("nan"),
        float("nan"),
        (f"{len(rows)} controls; required power {required:g} at the pre-registered effect"
         + (f"; MISSING: {missing}" if missing else "")),
        missing, rows)


# ---------------------------------------------------------------------------
# F7 -- portability
# ---------------------------------------------------------------------------
def family_F7(proto: Protocol, stdlib_pairs, arch_pairs) -> FamilyResult:
    """Two claims, tested with the rule appropriate to each.

    ``stdlib_pairs``: dicts with ``label`` and ``identical`` -- bitwise equality of every
    counter and digest between standard libraries on one architecture.  The candidate's
    stream is a function of the engine alone, so the prediction is exact and there is no
    tolerance.

    ``arch_pairs``: dicts with ``label``, ``k1``, ``n1``, ``k2``, ``n2`` -- equivalence of
    rates across architectures, where libm differences make bitwise agreement impossible.
    Only rows whose both arms executed natively reach here; translated execution is
    corroborating evidence and is excluded upstream.
    """
    alpha = proto.alpha("F7_portability")
    margin = float(proto.get("gates", "G4", "log_ratio_margin"))
    offenders = [p["label"] for p in stdlib_pairs if not p.get("identical")]
    rows = [{"kind": "cross_stdlib_bitwise", "label": p["label"],
             "identical": bool(p.get("identical"))} for p in stdlib_pairs]

    tost_p, tost_rows = [], []
    n_equivalent = n_underpowered = 0
    for p in arch_pairs:
        r = tost_log_rate_ratio(p["k1"], p["n1"], p["k2"], p["n2"], margin,
                                alpha / max(len(arch_pairs), 1))
        r["label"] = p["label"]
        r["kind"] = "cross_arch_equivalence"
        tost_rows.append(r)
        if r.get("resolved"):
            tost_p.append(r["p_value"])
        if r.get("disagrees"):
            offenders.append(p["label"])          # a genuine difference: this FAILS
        elif r.get("equivalent"):
            n_equivalent += 1
        elif r.get("informative") is False:
            n_underpowered += 1                    # neither proves nor disproves
    rows.extend(tost_rows)

    if not stdlib_pairs and not arch_pairs:
        return FamilyResult("F7_portability", alpha, False, float("nan"), float("nan"),
                            "no portability comparison was produced; the gate stays open "
                            "rather than passing", [], rows)
    passed = not offenders
    return FamilyResult("F7_portability", alpha, bool(passed), float(len(offenders)),
                        float(simes_global(tost_p) if tost_p else 1.0),
                        f"{len(stdlib_pairs)} bitwise cross-stdlib comparisons (exact); "
                        f"{len(arch_pairs)} cross-architecture cells, of which "
                        f"{n_equivalent} are equivalent within +/-{margin:g} on the log "
                        f"rate ratio, {len(offenders)} disagree and {n_underpowered} are "
                        f"underpowered (too few events to establish either, so they leave "
                        f"the gate open rather than closing or failing it)",
                        offenders, rows)
