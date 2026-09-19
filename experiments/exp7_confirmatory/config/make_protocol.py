#!/usr/bin/env python3
"""Emit ``config/protocol.json``, the machine-readable half of ``PROTOCOL.md``.

Everything the analysis is allowed to decide with is computed here and frozen, before any
Experiment 7 data exists.  In particular the per-cell *achieved* coverage of every quantile
interval is a closed-form function of ``(n, p, correction)``, so the null distribution of the
F2 miss count -- and therefore its critical value -- is known in advance rather than being
chosen once the misses have been counted.

Deterministic: same inputs, byte-identical output.  ``make protocol-check`` re-runs it and
diffs, so the committed file cannot drift from the code that generated it.
"""
from __future__ import annotations

import json
import os
import numpy as np
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))

# --- frozen matrix ---------------------------------------------------------------------
KAPPA_LADDER = [0.5001, 0.501, 0.505, 0.51, 0.55, 0.60, 0.75, 1.0, 1.25, 1.49, 1.5, 2.0, 5.0]
PRECISIONS = ["float", "double"]
SEEDS_PRODUCTION = [7001, 7002, 7003, 7004, 7005]
SEEDS_PERFORMANCE = [7006, 7007, 7008, 7009, 7010]
QUANTILE_LEVELS = [0.5, 0.9, 0.99, 0.999, 0.9999]
TAIL_Q0 = [1e-2, 1e-3, 1e-4]

N_SCALAR = 1_000_000
N_MECHANISM = 1_000_000
N_CONDITIONING = 1_000_000
N_LOADER = 100_000

# --- frozen statistics -----------------------------------------------------------------
ALPHA_TOTAL = 0.05
FAMILY_ALPHA = {
    "F1_radial_law": 0.010,
    "F2_quantile_coverage": 0.005,
    "F3_quantile_direction": 0.010,
    "F4_upper_tail_mass": 0.010,
    "F5_loader_battery": 0.010,
    "F7_portability": 0.005,
}
INTERVAL_CONF = 0.95
BOOTSTRAP_CONF = 0.99
BOOTSTRAP_RESAMPLES = 10_000
LOG_Q_SWITCH = -30.0
LOG_Q_ORACLE_TOL = 1e-10

# A quantile bracket wider than this many natural-log units cannot distinguish a correct
# sampler from a badly wrong one, so it is declared non-informative in advance.
INFORMATIVE_WIDTH_LOG_UNITS = 1.0

# A returned value can be finite and still be wrong.  Experiment 6 had no category for that:
# its classifier scored any finite return as a success even when the intended draw was not
# representable, which is exactly the regime where the legacy form is worst -- its own oracle
# measured up to 0.48 relative error on the surviving radius once the denominator went
# subnormal.  A draw is declared FINITE_BUT_WRONG when its returned radius has lost more than
# half of the significand of its type, i.e. relative error against the arbitrary-precision
# oracle above 2^-(digits/2): 1.5e-8 in double, 2.4e-4 in float.  Stating it in bits rather
# than as an absolute number makes it type-aware and leaves the log path a wide margin -- its
# own error is about eps*|log R|, which is 1.5e-13 at the double overflow threshold and
# 5.3e-6 at the float one.
ACCURACY_MAX_REL_ERROR_BITS_LOST_FRACTION = 0.5

# G5 acceptability bound on candidate/LEGACY time per returned sample.
PERFORMANCE_BOUND = 2.0
# G4 equivalence margin on the log rate ratio between architectures.
PORTABILITY_LOG_RATIO_MARGIN = 0.15
# F6 required power at the pre-registered minimum detectable effect.
NC_REQUIRED_POWER = 0.90
NC_INJECTION_LOSS_FRACTIONS = [1e-3, 1e-4]
NC_INJECTIONS = 200


def order_statistic_interval(n: int, p: float, conf: float, bonferroni: int):
    """0-based index bracket ``(lo, hi)``; ``None`` when the sample cannot resolve ``p``."""
    if n <= 0 or not (0.0 < p < 1.0):
        return None
    alpha = (1.0 - conf) / bonferroni
    lo = int(stats.binom.ppf(alpha / 2, n, p)) - 1
    hi = int(stats.binom.isf(alpha / 2, n, p))
    lo = max(lo, 0)
    hi = min(hi, n - 1)
    if hi <= lo or hi >= n:
        return None
    if hi == n - 1 and p > 1 - 3.0 / n:
        return None
    return (lo, hi)


def achieved_coverage(n: int, p: float, lo: int, hi: int) -> float:
    """Exact coverage of ``[x_(lo), x_(hi)]`` for the ``p`` quantile.

    The interval covers the quantile iff the number of sample points below it lies in
    ``[lo+1, hi]``, i.e. ``P(lo+1 <= B <= hi)`` for ``B ~ Binomial(n, p)``.
    """
    return float(stats.binom.cdf(hi, n, p) - stats.binom.cdf(lo, n, p))


def f2_cells() -> tuple[list[dict], dict]:
    """Every quantile interval in the run, with its achieved coverage frozen."""
    cells, miss_probs = [], []
    for precision in PRECISIONS:
        for kappa in KAPPA_LADDER:
            for seed in SEEDS_PRODUCTION:
                for p in QUANTILE_LEVELS:
                    br = order_statistic_interval(N_SCALAR, p, INTERVAL_CONF,
                                                  len(QUANTILE_LEVELS))
                    if br is None:
                        cells.append({"precision": precision, "kappa": kappa, "seed": seed,
                                      "p": p, "resolved": False, "achieved_coverage": None})
                        continue
                    lo, hi = br
                    cov = achieved_coverage(N_SCALAR, p, lo, hi)
                    cells.append({"precision": precision, "kappa": kappa, "seed": seed,
                                  "p": p, "resolved": True, "lo_index": lo, "hi_index": hi,
                                  "achieved_coverage": cov})
                    miss_probs.append(1.0 - cov)

    # Poisson-binomial upper tail of the miss count, by exact convolution.
    pmf = np.array([1.0])
    for q in miss_probs:
        pmf = np.convolve(pmf, [1.0 - q, q])
    alpha = FAMILY_ALPHA["F2_quantile_coverage"]
    surv = 1.0 - np.cumsum(pmf)                       # surv[k] = P(misses > k)
    critical = int(np.argmax(surv <= alpha))          # smallest k with P(misses > k) <= alpha
    return cells, {
        "n_resolved_intervals": len(miss_probs),
        "expected_misses": float(sum(miss_probs)),
        "alpha": alpha,
        "critical_miss_count": critical,
        "reject_if_misses_exceed": critical,
        "attained_level": float(surv[critical]),
        "rule": ("fail F2 iff the observed miss count exceeds critical_miss_count; the null "
                 "is the Poisson-binomial of the frozen per-cell achieved coverages"),
    }


def main() -> None:
    cells, f2 = f2_cells()
    protocol = {
        "protocol_version": "1.1.0",
        "amendments": [
            {"version": "1.1.0",
             "before_any_data": True,
             "reason": "An independent audit of Experiment 6's evidence, completed after "
                       "1.0.0 was written and before any Experiment 7 draw existed, found "
                       "three gaps that no decision rule in 1.0.0 covered: its terminal "
                       "classifier had no FINITE_BUT_WRONG category and scored a finite "
                       "return as a success even for an unrepresentable draw; its oracle "
                       "audit sampled bulk failures at 1 in 4957 where the plan says every "
                       "public-path failure; and its accuracy statistics pooled float and "
                       "double. This amendment adds the accuracy threshold, the audit "
                       "stratum rates, and the requirement that accuracy be reported per "
                       "precision. It adds criteria; it relaxes none, and no Experiment 7 "
                       "datum had been generated when it was made."},
        ],
        "experiment": "exp7_confirmatory",
        "frozen_before_any_data": True,
        "candidate": {
            "name": "CANDIDATE",
            "description": "released loader cpp/bi_kappa_distribution.H at version 2.0.0, "
                           "radius built in the log domain",
            "version": "2.0.0",
            "log_gamma_primitive": "LOG-ID",
            "log_gamma_citation": "Ahrens, J.H. and Dieter, U. (1974), Computing 12, 223-246",
            "boosted_gamma_primitive": "Marsaglia-Tsang",
            "boosted_gamma_citation": "Marsaglia, G. and Tsang, W.W. (2000), "
                                      "ACM TOMS 26(3), 363-372, doi:10.1145/358407.358414",
            "uniform_mapping": "(k + 1/2) / 2^(digits-1) from engine bits; open interval; "
                               "no endpoint redraw and no endpoint clamp",
        },
        "comparator": {
            "name": "LEGACY",
            "description": "the same header at version 1.0.0, vendored unmodified",
            "version": "1.0.0",
            "radius": "sqrt(X1)/sqrt(X2) with std::gamma_distribution",
        },
        "seeds": {
            "production": SEEDS_PRODUCTION,
            "performance": SEEDS_PERFORMANCE,
            "declared_not_derived": True,
            "disjoint_from": {"exp1": [1001, 1005], "exp2": [2001, 2005],
                              "exp3": [3001, 3003, 3101], "exp4_exp6": [4001, 4010]},
        },
        "matrix": {
            "kappa_ladder": KAPPA_LADDER,
            "precisions": PRECISIONS,
            "quantile_levels": QUANTILE_LEVELS,
            "tail_q0": TAIL_Q0,
            "n_scalar": N_SCALAR,
            "n_mechanism": N_MECHANISM,
            "n_conditioning": N_CONDITIONING,
            "n_loader": N_LOADER,
            "loader_cases": ["C0", "C1", "C2", "C3", "C4", "C5", "C6"],
            "benchmark_cases": ["B0", "B1", "B2", "B3", "B4"],
        },
        "statistics": {
            "alpha_total": ALPHA_TOTAL,
            "family_alpha": FAMILY_ALPHA,
            "alpha_sum_check": round(sum(FAMILY_ALPHA.values()), 12),
            "interval_conf": INTERVAL_CONF,
            "bootstrap": {"kind": "cluster over seeds", "conf": BOOTSTRAP_CONF,
                          "resamples": BOOTSTRAP_RESAMPLES},
            "log_q_switch": LOG_Q_SWITCH,
            "log_q_oracle_tolerance": LOG_Q_ORACLE_TOL,
            "informative_width_log_units": INFORMATIVE_WIDTH_LOG_UNITS,
            "expected_false_failure_bound": round(sum(FAMILY_ALPHA.values()), 12),
        },
        "F2": f2,
        "F2_cells": cells,
        "F6": {
            "required_power": NC_REQUIRED_POWER,
            "injection_loss_fractions": NC_INJECTION_LOSS_FRACTIONS,
            "injections_per_effect": NC_INJECTIONS,
            "controls": ["NC1_radius_direction_coupling", "NC2_capped_vs_uncapped_weak_cap",
                         "NC3_survivor_conditioning"],
            "missing_control_is_failure": True,
        },
        "accuracy": {
            "bits_lost_fraction": ACCURACY_MAX_REL_ERROR_BITS_LOST_FRACTION,
            "max_relative_error": {
                # 2^-(digits * fraction); digits = 53 (double), 24 (float)
                "double": 2.0 ** -(53 * ACCURACY_MAX_REL_ERROR_BITS_LOST_FRACTION),
                "float": 2.0 ** -(24 * ACCURACY_MAX_REL_ERROR_BITS_LOST_FRACTION),
            },
            "rule": "a returned draw whose radius differs from the arbitrary-precision "
                    "oracle by more than max_relative_error for its type is classified "
                    "FINITE_BUT_WRONG and counted as a loss, not as a success",
            "report_per_precision": True,
            "never_pool_precisions": True,
        },
        "audit": {
            "rule": "every public-path failure is adjudicated; the remaining strata are "
                    "sampled at the declared rates and BOTH the audited and the total count "
                    "of every stratum is reported, so coverage is visible rather than implied",
            "strata_rates": {
                "public_path_failure": 1.0,
                "near_limit_2_log_units": 1.0,
                "method_disagreement": 1.0,
                "finite_but_wrong_candidate": 1.0,
                "subnormal_denominator": 0.1,
                "uniform_sample": 1e-4,
            },
            "require_header_record_even_when_zero_disagreements": True,
            "require_audited_file_sha256": True,
        },
        "power_study": {
            "artifact": "config/power_study.json",
            "computed_before_any_data": True,
            "note": "Anderson-Darling, KS and Cramer-von Mises on Z are blind to survivor "
                    "conditioning below a 1e-3 loss fraction at the frozen n; the F4 "
                    "exceedance tests are not. P1 output must therefore carry exceedance "
                    "counts and excesses, not only an ECDF grid.",
        },
        "gates": {
            "G4": {"cross_stdlib_rule": "bitwise equality, no tolerance",
                   "cross_arch_rule": "two one-sided tests on the log rate ratio",
                   "log_ratio_margin": PORTABILITY_LOG_RATIO_MARGIN,
                   "translated_execution_excluded": True},
            "G5": {"time_per_return_ratio_bound": PERFORMANCE_BOUND},
            "G6": {"requires_make_verify_exit_zero": True,
                   "requires_byte_identical_regeneration": True,
                   "requires_fresh_extraction_check": True},
        },
        "attribution_checklist": [
            {"method": "Gamma-ratio / Beta-prime Kappa construction",
             "doi": "10.1063/5.0117628", "source": "Zenitani & Nakano (2022)",
             "equation": "Kappa velocity as a Gamma ratio", "domain": "kappa > 1/2",
             "output_scale": "linear", "limitation": "no finite-precision treatment"},
            {"method": "log-transformation Gamma generators",
             "doi": "10.18637/jss.v055.i04", "source": "Xi, Tan & Liu (2013)",
             "equation": "Algorithm 1, ratio-of-uniforms with log output below a < 0.01",
             "domain": "alpha > 0", "output_scale": "log below a hard-coded switch",
             "limitation": "piecewise-linear bounding tables fitted at efficiency R = 0.9"},
            {"method": "small-shape Gamma generation with log-scale output",
             "doi": "10.1007/s00180-016-0692-0", "source": "Liu, Martin & Syring (2017)",
             "equation": "Sect. 3 p. 1770 envelope; acceptance rate Eq. (2) p. 1771",
             "domain": "0 < alpha < 1", "output_scale": "log",
             "limitation": "acceptance rate has a pole at alpha = 1; undefined at alpha >= 1"},
            {"method": "Gamma generators for shape below one",
             "doi": "arXiv:2411.01415", "source": "Zenitani (2024)",
             "equation": "shape-boosting and rejection variants", "domain": "alpha < 1",
             "output_scale": "linear", "limitation": "linear output underflows at small shape"},
            {"method": "Kappa denominator-zero hazard",
             "doi": "10.1029/2025JA034669", "source": "Zenitani, Usami & Matsukiyo (2026)",
             "equation": "low-parameter piecewise-rejection alternative",
             "domain": "low kappa", "output_scale": "linear",
             "limitation": "changes the high-level construction"},
            {"method": "shape boosting identity (the selected primitive)",
             "doi": "10.1007/BF02293108", "source": "Ahrens & Dieter (1974), Computing 12",
             "equation": "X = Y U^(1/a), Y ~ Gamma(a+1,1)", "domain": "a > 0",
             "output_scale": "linear; written here on the log scale",
             "limitation": "log X carries about one ulp of |log U|/a"},
            {"method": "boosted Gamma draw",
             "doi": "10.1145/358407.358414", "source": "Marsaglia & Tsang (2000)",
             "equation": "d = shape - 1/3, c = 1/sqrt(9d), squeeze u < 1 - 0.0331 x^4",
             "domain": "shape >= 1", "output_scale": "linear",
             "limitation": "rejection-based; expected attempts below 1.01"},
        ],
        "acceptance_rate_predictor": {
            "note": "Named here before the run so that it is a test rather than a "
                    "description.  For any acceptance-rejection primitive the measured "
                    "acceptance rate converges to the ratio of areas integral(h)/integral(eta), "
                    "not to a rate quoted under a different normalization.  The selected "
                    "primitive is rejection-free at the boosting step, so this applies only "
                    "to the Marsaglia-Tsang inner loop, whose predicted attempts per draw is "
                    "computed from its own envelope.",
            "rule": "area_ratio",
        },
    }

    out = os.path.join(HERE, "protocol.json")
    with open(out, "w") as fh:
        json.dump(protocol, fh, indent=2, sort_keys=False)
        fh.write("\n")
    print(f"wrote {out}")
    print(f"  F2: {f2['n_resolved_intervals']} resolved intervals, "
          f"expected misses {f2['expected_misses']:.2f}, "
          f"fail if misses > {f2['critical_miss_count']} "
          f"(attained level {f2['attained_level']:.4f})")
    print(f"  alpha budget: {protocol['statistics']['alpha_sum_check']}")


if __name__ == "__main__":
    main()
