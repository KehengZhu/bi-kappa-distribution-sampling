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
        "protocol_version": "1.3.0",
        "amendments": [
            {"version": "1.3.0",
             "before_any_data": True,
             "reason": "The frozen F5 loader battery contained no upper-tail statistic, and "
                       "section 2.4's own power study had already established that nothing "
                       "else has power against survivor conditioning below a 1e-3 loss "
                       "fraction. Section 5.4 then required negative control NC1 to be "
                       "detected at 1e-3 AND 1e-4 with power 0.90. The protocol therefore "
                       "demanded a power the battery it froze could not deliver, and a "
                       "correct candidate would have failed F6, G3 and the verdict. "
                       "Measured at production scale, n = 5e5 over 40 injections: the F5 "
                       "independence test detects 0 of 40 at either effect size; an "
                       "upper-tail exceedance statistic detects 40 of 40 at both. The "
                       "arithmetic is not subtle -- removing 50 draws from 500000 cannot "
                       "move an 8x8 chi-square whose cells expect 7800, and removes "
                       "essentially every exceedance above q0 = 1e-4. F5 therefore gains an "
                       "exceedance statistic on each loader cell's recovered radial law, "
                       "mirroring F4. This ADDS a test to the family rather than relaxing "
                       "one, and it closes a real gap: a battery certifying a heavy-tailed "
                       "law had no statistic that looks at the tail. Four smaller "
                       "corrections travel with it, listed under 'corrections' below. No "
                       "Experiment 7 holdout datum existed when this was made."},
            {"version": "1.2.0",
             "before_any_data": True,
             "reason": "Amendment 1.1.0 required every public-path failure to be "
                       "adjudicated. Costing it showed that is infeasible and, worse, "
                       "unnecessary: the honest failures alone come to 17.5 million across "
                       "the frozen matrix, about 5.4 CPU-hours of oracle time, and almost "
                       "all of them are draws whose intended value lies thousands of log "
                       "units outside the type's range, where no arithmetic could have "
                       "returned them and the classification is not in doubt. Committing "
                       "to a rule that cannot be followed is how Experiment 6 ended up "
                       "sampling bulk failures at 1 in 4957 while its plan said 'every'. "
                       "The rule is therefore restated by MARGIN rather than by category: "
                       "everything whose classification could go either way is adjudicated "
                       "in full, and a draw beyond a 20-log-unit margin -- four thousand "
                       "times the worst float disagreement between the working-precision "
                       "reference and the 100-digit oracle -- is sampled, with the margin "
                       "asserted for every one of them. This narrows nothing that was "
                       "decision-relevant and makes the coverage claim one that is "
                       "actually kept."},
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
            "rule": "every DECISION-RELEVANT attempt is adjudicated in full; an attempt is "
                    "decision-relevant when its classification could go either way. An "
                    "attempt whose intended log-component lies more than "
                    "margin_log_units beyond a type limit is UNAMBIGUOUS -- no arithmetic "
                    "could have returned it -- and is sampled instead, with the margin "
                    "itself asserted for every such attempt in working precision. BOTH the "
                    "audited and the total count of every stratum is reported, so coverage "
                    "is visible rather than implied.",
            "margin_log_units": 20.0,
            "margin_justification":
                "The working-precision reference and the 100-digit oracle differ by at most "
                "about 7.3e-12 natural-log units in double and 4.2e-3 in float (Experiment "
                "6's own oracle record). A 20-log-unit margin is over four thousand times "
                "the worst float discrepancy, so beyond it a misclassification is "
                "impossible rather than merely improbable. That is an argument the sampling "
                "rate does not have to carry.",
            "strata_rates": {
                "within_margin_of_a_type_limit": 1.0,
                "method_disagreement": 1.0,
                "finite_but_wrong_candidate": 1.0,
                "avoidable_loss_candidate": 1.0,
                "subnormal_or_zero_denominator_representable_target": 1.0,
                "unambiguous_failure_beyond_margin": 1e-3,
                "uniform_sample": 1e-4,
            },
            "require_header_record_even_when_zero_disagreements": True,
            "require_audited_file_sha256": True,
            "require_margin_assertion_for_unambiguous": True,
        },
        "F5_tail": {
            "statistic": "exceedance count above z0 = -log q0 on the recovered radial Z of "
                         "each uncapped loader cell, exact binomial against "
                         "Binomial(n_attempted, q0), plus KS of the excesses against Exp(1)",
            "q0": TAIL_Q0,
            "rationale": "added by amendment 1.3.0; without it the loader battery has no "
                         "statistic with power against the effect it certifies absent",
            "capped_cells_excluded": True,
            "capped_cells_reason":
                "under a cap the accepted radius is truncated at a direction-dependent "
                "bound, so the count above z0 is not Binomial(n, q0)",
            "count_member_requires_q0_above_honest_floor": True,
            "count_member_applicability":
                "The count member applies at a threshold q0 only when q0 exceeds the cell's "
                "honest-overflow rate f. The reason is what is observable, not a "
                "convenience. When q0 > f the threshold z0 = -log q0 lies below the "
                "overflow threshold z_f = -log f, so every attempt above z0 is either a "
                "returned survivor above z0 or an overflowed attempt, and both are counted "
                "exactly. When q0 < f every attempt above z0 has overflowed, and separating "
                "those above z0 from those merely above z_f needs the intended value of a "
                "draw that has none -- the loader returns no number for it. Such a cell is "
                "reported not-applicable for that threshold, never passed and never failed "
                "on it, with its honest floor published on the row so a reader can see why. "
                "Without this, the member would silently require every uncapped cell's loss "
                "to fall below 1e-4, which case C4 exists precisely to violate: its floor "
                "is 8.25e-4, so honest overflow alone would give p ~ 1.7e-223 and fail a "
                "correct candidate.",
            "excess_member_applies_at_every_threshold": True,
            "unresolved_Z_counts_toward_every_threshold": True,
            "unresolved_Z_note":
                "'Unresolved' here means the Z transform itself underflowed, so the draw is "
                "further into the tail than any threshold. It does NOT mean the velocity "
                "overflowed: an overflowed draw has a perfectly ordinary Z and must be "
                "counted by it, not treated as exceeding everything.",
        },
        "corrections": {
            "F1_contents": "Anderson-Darling, KS and Cramer-von Mises on Z only. The Beta "
                           "test on W is dropped: the schema emits no W sample, Z is a "
                           "monotone transform of log W so the two are near-redundant, and "
                           "G0 already validates the log-q evaluator against an "
                           "arbitrary-precision incomplete beta to 1e-10, which is a "
                           "stronger check than an overlap comparison between the routes.",
            "NC3_effect": "the largest loss fraction among the UNCAPPED LOADER CELLS "
                          "(C0-C4), not the global maximum over the scalar ladder. The "
                          "earlier wording said 'whose loss is not dominated by honest "
                          "overflow', which is unusable: G1 requires the candidate's "
                          "avoidable loss to be exactly zero, so every loss it has is "
                          "honest overflow and that qualifier empties the set. What it was "
                          "there to exclude is the degenerate scalar configurations -- "
                          "float kappa = 0.5001 loses 98 per cent of its draws, where "
                          "detection is trivial and says nothing -- not any cell from among "
                          "C0-C4.",
            "G5_interval": "0.99, the cluster-bootstrap level in statistics.bootstrap.conf. "
                           "Section 6 said 95 per cent in prose; the machine-readable value "
                           "governs and is the conservative choice for a gate on an upper "
                           "limit.",
            "F2_informativeness": "misses are counted over every RESOLVED interval, which "
                                  "is what the frozen Poisson-binomial null was built over. "
                                  "The informativeness flag is published as a map per "
                                  "section 5.1 and does not change the count; excluding "
                                  "non-informative cells would make F2 unevaluable against "
                                  "its own frozen null.",
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
