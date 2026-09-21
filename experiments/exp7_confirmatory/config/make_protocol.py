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
# Seed block for the SECOND holdout, drawn by the rule stated in PROTOCOL.md section 3 and
# not by choice: take the highest seed declared anywhere in this repository (7505, the
# Experiment 7 selftest fixtures), round up to the next multiple of 1000 (8000), and take the
# next ten integers.  The rule admits exactly one answer, so no seed was selected after any
# result was seen, and the block is disjoint from every other by construction.  The first
# holdout's block, 7001-7010, is spent: PROTOCOL.md section 8 forbids reusing it, and its
# NO-GO result is preserved in commit 45d3ef8.
SEEDS_FIRST_HOLDOUT = [7001, 7002, 7003, 7004, 7005, 7006, 7007, 7008, 7009, 7010]
SEEDS_SECOND_HOLDOUT = [8001, 8002, 8003, 8004, 8005, 8006, 8007, 8008, 8009, 8010]
SEEDS_PRODUCTION = [9001, 9002, 9003, 9004, 9005]
SEEDS_PERFORMANCE = [9006, 9007, 9008, 9009, 9010]

# Protocol 3.0.0: P1 and P5 give every (precision, kappa) configuration of a replicate its
# own engine stream, so that the 26 configurations of a replicate are independent rather
# than driven by one shared mt19937.  The stream seed is
#
#     base + P1_STREAM_SEED_STRIDE * (13 * precision_index + kappa_index)
#
# with precision_index 0 for double and 1 for float, and kappa_index the 0-based position in
# KAPPA_LADDER.  The formula is declared here, in PROTOCOL.md section 3, and in
# src/exp7_common.H, it admits exactly one answer for each configuration, and every value it
# produces is 9001-9005 modulo 10000, so it cannot collide with a declared block anywhere in
# this repository.
P1_STREAM_SEED_STRIDE = 10_000
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
# The environments acceptance is limited to under protocol 2.0.0.  Both must be present,
# native and complete for G4 to close; see gates.G4.cross_arch_withdrawn_reason.
SUPPORTED_ENVIRONMENTS = [
    {"arch": "arm64", "os": "macOS", "stdlib": "libc++", "tag": "libcxx"},
    {"arch": "arm64", "os": "macOS", "stdlib": "libstdc++", "tag": "libstdcxx"},
]
# F6 required power at the pre-registered minimum detectable effect.
NC_REQUIRED_POWER = 0.90
NC_INJECTION_LOSS_FRACTIONS = [1e-3, 1e-4]
NC_INJECTIONS = 200
# NC3's effect under protocol 2.0.0 is conditioning IN EXCESS of the representability floor,
# measured as a fraction of attempts removed from among the draws the type would have allowed
# back.  These are pre-registered from the development calibration in
# docs/revision/experiments/f5_tail_calibration.md, on simulated cells only, before any
# second-holdout datum existed.
NC3_EXCESS_LOSS_FRACTIONS = [1e-3, 1e-4]


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


# Monte Carlo settings for the per-configuration miss-count law.  Fixed here so that
# `make protocol-check` regenerates config/protocol.json byte for byte.
F2_GROUP_MC_REPLICATES = 20_000_000
F2_GROUP_MC_SEED = 20260920


def group_miss_pmf() -> tuple[np.ndarray, np.ndarray]:
    """The miss-count law of ONE configuration: five dependent level indicators.

    The five brackets of a configuration are read off one sample of ``N_SCALAR`` draws, so
    their miss indicators are dependent however the streams are seeded.  That dependence is
    exactly characterised: an interval at level ``p`` covers the quantile iff the number of
    sample points at or below it lies in ``[lo + 1, hi]``, and under the null that count is
    the running sum of a multinomial over the six intervals the five levels cut.  Nothing
    about the sampler enters, so the law below is the null's, not the candidate's.
    """
    brackets = [order_statistic_interval(N_SCALAR, p, INTERVAL_CONF, len(QUANTILE_LEVELS))
                for p in QUANTILE_LEVELS]
    if any(b is None for b in brackets):
        raise SystemExit("F2: a quantile level resolves no bracket")
    pv = np.diff([0.0] + list(QUANTILE_LEVELS) + [1.0])
    rng = np.random.default_rng(F2_GROUP_MC_SEED)
    counts = np.zeros(len(QUANTILE_LEVELS) + 1, dtype=np.int64)
    done, chunk = 0, 500_000
    while done < F2_GROUP_MC_REPLICATES:
        c = min(chunk, F2_GROUP_MC_REPLICATES - done)
        k = np.cumsum(rng.multinomial(N_SCALAR, pv, size=c)[:, :len(QUANTILE_LEVELS)],
                      axis=1)
        miss = np.zeros(k.shape, dtype=bool)
        for j, (lo, hi) in enumerate(brackets):
            miss[:, j] = (k[:, j] < lo + 1) | (k[:, j] > hi)
        counts += np.bincount(miss.sum(axis=1), minlength=len(QUANTILE_LEVELS) + 1)
        done += c
    return counts / float(F2_GROUP_MC_REPLICATES), counts


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

    n_groups = len(PRECISIONS) * len(KAPPA_LADDER) * len(SEEDS_PRODUCTION)
    group_size = len(QUANTILE_LEVELS)
    if len(miss_probs) != n_groups * group_size:
        raise SystemExit("F2: not every level resolves a bracket; the grouping is wrong")

    pmf_c, mc = group_miss_pmf()
    alpha = FAMILY_ALPHA["F2_quantile_coverage"]
    pmf = np.array([1.0])
    for _ in range(n_groups):
        pmf = np.convolve(pmf, pmf_c)
    surv = 1.0 - np.cumsum(pmf)                       # surv[k] = P(misses > k)
    critical = int(np.argmax(surv <= alpha))          # smallest k with P(misses > k) <= alpha
    expected = float((pmf * np.arange(pmf.size)).sum())

    # Independence across the five levels of one configuration is NOT assumed; it is also
    # not far off, and the difference is published rather than left implicit.  This is the
    # critical value the withdrawn Poisson-binomial null gave.
    pmf_i = np.array([1.0])
    for q in miss_probs:
        pmf_i = np.convolve(pmf_i, [1.0 - q, q])
    surv_i = 1.0 - np.cumsum(pmf_i)
    critical_if_intervals_independent = int(np.argmax(surv_i <= alpha))

    # The far tail of the per-configuration law is the only part the Monte Carlo resolves
    # poorly, so the decision is shown not to depend on it.
    stability = []
    for label, scale in (("cells_4_5_zeroed", 0.0), ("cells_4_5_doubled", 2.0),
                         ("cells_3_4_5_doubled", None)):
        q = pmf_c.copy()
        if scale is None:
            q[3] *= 2.0
            q[4] *= 2.0
            q[5] *= 2.0
        else:
            q[4] *= scale
            q[5] *= scale
        q = q / q.sum()
        t = np.array([1.0])
        for _ in range(n_groups):
            t = np.convolve(t, q)
        st = 1.0 - np.cumsum(t)
        stability.append({"perturbation": label,
                          "critical_miss_count": int(np.argmax(st <= alpha))})

    return cells, {
        "n_resolved_intervals": len(miss_probs),
        "n_groups": n_groups,
        "group_size": group_size,
        "expected_misses": expected,
        "alpha": alpha,
        "critical_miss_count": critical,
        "reject_if_misses_exceed": critical,
        "attained_level": float(surv[critical]),
        "group_miss_pmf": [float(x) for x in pmf_c],
        "group_miss_counts": [int(x) for x in mc],
        "critical_if_intervals_independent": critical_if_intervals_independent,
        "far_tail_stability": stability,
        "per_level_miss_probability": [1.0 - achieved_coverage(
            N_SCALAR, p, *order_statistic_interval(N_SCALAR, p, INTERVAL_CONF,
                                                   len(QUANTILE_LEVELS)))
            for p in QUANTILE_LEVELS],
        "monte_carlo": {
            "replicates": F2_GROUP_MC_REPLICATES,
            "rng": "numpy.random.default_rng(%d).multinomial" % F2_GROUP_MC_SEED,
            "rng_seed": F2_GROUP_MC_SEED,
            "what_is_simulated":
                "one configuration's sample of n_scalar draws, as the multinomial cell "
                "counts it puts into the six intervals the five quantile levels cut (0, "
                "p1], (p1, p2], ..., (p5, 1). The running sums of those counts are the "
                "K_k = #{Z_i <= z_pk}, and the interval at level k misses exactly when "
                "K_k lies outside [lo_k + 1, hi_k]. Nothing about the sampler enters: the "
                "null is the exact law, so the transform of the draws is Uniform(0,1) by "
                "construction and only the order statistics matter.",
        },
        "rule": ("fail F2 iff the observed miss count exceeds critical_miss_count; the null "
                 "is the n_groups-fold convolution of the per-configuration miss-count "
                 "distribution, which is exact because PROTOCOL.md section 3 gives every "
                 "configuration its own engine stream"),
        "why_not_poisson_binomial":
            "Protocol 2.0.0 and earlier compared the miss count against a Poisson-binomial "
            "over INDEPENDENT intervals. The 650 intervals were not independent, for two "
            "reasons of different size. (a) Every configuration of a replicate was driven "
            "by one mt19937 seeded with the replicate's seed, so the 26 configurations saw "
            "the same uniforms and their quantile errors co-moved almost exactly; the "
            "effective number of independent units was about five, not 650. (b) The five "
            "levels of one configuration are order statistics of one sample and are "
            "dependent whatever the streams do. Protocol 3.0.0 removes (a) by construction "
            "-- independent streams per configuration -- and computes (b) exactly instead "
            "of assuming it away. The statistic, the family alpha and the set of intervals "
            "counted are unchanged.",
    }


def main() -> None:
    cells, f2 = f2_cells()
    protocol = {
        "protocol_version": "3.0.0",
        "amendments": [
            {"version": "3.0.0",
             "before_any_data": False,
             "governs": "the third confirmatory holdout, on seeds 9001-9010",
             "reason":
                 "The second holdout, on seeds 8001-8010 under protocol 2.0.0, returned "
                 "NO-GO on gate G1 and is preserved unmodified in commit e5c9837. Four of "
                 "the five gates that had failed or stood open on 7001-7010 passed, and "
                 "both defects amendment 2.0.0 identified are closed: G2 recorded zero "
                 "oracle disagreements in 14235020 adjudicated attempts, F5 passed jointly "
                 "over 52 tests, G4 closed bitwise, G5 came in at 0.854x. G1 failed for "
                 "two reasons, and this amendment answers both. PROTOCOL.md section 8 "
                 "permits exactly one path after a failure -- identify a concrete defect, "
                 "fix it, freeze a new implementation hash AND a new protocol document, "
                 "draw a further disjoint seed block, and rerun -- and this is that "
                 "document.\n"
                 "\n"
                 "(1) IMPLEMENTATION. One avoidable loss, at float kappa = 0.505, seed "
                 "8005, attempt 157855, seen identically in both standard-library streams. "
                 "Its largest intended component lay 3.12e-08 natural-log units below "
                 "log(FLT_MAX) -- about half an ulp -- so the correctly rounded float is "
                 "finite and the draw is returnable. Two working-precision quantities on "
                 "the path to it are coarser than that margin: the float log radius, whose "
                 "error is 3.5e-06 after the division by a = 0.005, and the float order-"
                 "unity vector g, worth another 3.9e-08. The loader exponentiated one ulp "
                 "high and returned a non-finite component. Release 2.2.0 draws the same "
                 "variates in the same order in float -- the sampled law does not change -- "
                 "and evaluates the deterministic map from those variates to the returned "
                 "velocity in double, rounding once, where the component is materialized. "
                 "That puts the resolution of the decision near 1e-14 against a 6e-08 ulp. "
                 "A double instantiation is unaffected: its accumulator is its working "
                 "type, and 2.2.0 returns the same bits as 2.1.0 for every double "
                 "configuration in the matrix, capped and uncapped. The losing draw and "
                 "34 magnitudes straddling the float representability boundary are frozen "
                 "as deterministic selftest fixtures, replayed from recorded bit patterns "
                 "rather than from an engine, so the regression survives any later change "
                 "of seeds, phase or ladder.\n"
                 "\n"
                 "(2) PROTOCOL. F2 compared its miss count against a Poisson-binomial over "
                 "INDEPENDENT intervals, and the intervals were not independent. Every "
                 "configuration of a replicate was driven by one mt19937 seeded with the "
                 "replicate's seed, so all 26 configurations saw the same uniform stream "
                 "and their quantile errors co-moved: on 8001-8010, seed 8004 gave 103 of "
                 "130 signed errors positive and seed 8005 gave 31 of 130, and 20 of the "
                 "23 misses came from those two replicates. The effective number of "
                 "independent units was about five, not 650. Protocol 3.0.0 removes the "
                 "coupling rather than modelling it: P1 and P5 give every configuration of "
                 "a replicate its own declared engine stream, so the 130 configurations "
                 "are independent by construction. The remaining dependence -- the five "
                 "levels of one configuration are order statistics of one sample -- is "
                 "real and is now computed exactly instead of being assumed away. The "
                 "statistic (the miss count over every resolved interval), the family "
                 "alpha (0.005) and the set of intervals counted are all unchanged. The "
                 "critical value moves from 13 to 14 of 650, so the test is not relaxed in "
                 "any material sense; on the second holdout's own count of 23 it would "
                 "still have rejected, at p = 4.3e-07. Nothing here was chosen to make a "
                 "past number pass.\n"
                 "\n"
                 "SCOPE is unchanged from 2.0.0: arm64 macOS under both standard "
                 "libraries, cross-architecture withdrawn. No threshold is relaxed, no "
                 "cell removed, no family alpha changed, no result excluded. Every change "
                 "was made and committed before any datum on seeds 9001-9010 existed."},
            {"version": "2.0.0",
             "before_any_data": False,
             "governs": "the second confirmatory holdout, on seeds 8001-8010 (superseded)",
             "reason":
                 "The first holdout, on seeds 7001-7010 under protocol 1.3.0, returned "
                 "NO-GO. It is preserved unmodified in commit 45d3ef8 and is not reopened. "
                 "PROTOCOL.md section 8 permits exactly one path after a failure -- "
                 "identify a concrete implementation defect, fix it, freeze a new "
                 "implementation hash AND a new protocol document, draw a further disjoint "
                 "seed block, and rerun -- and this is that document. Two defects were "
                 "identified, one in the implementation and one in this protocol.\n"
                 "\n"
                 "(1) IMPLEMENTATION. The released loader decided representability by "
                 "comparing a logarithm against log(max()). log(max()) is a rounded value "
                 "and exp of it need not be finite: in float, exp(log(FLT_MAX)) is exactly "
                 "+inf. A float kappa = 0.51 draw whose log R equalled log(FLT_MAX) bit for "
                 "bit therefore passed the test, overflowed when it was exponentiated, and "
                 "was returned as (-inf, +inf, +inf) without being counted -- while all "
                 "three components it should have produced were representable, the largest "
                 "at 1.92e38 against a limit of 3.40e38. One such draw appeared in "
                 "14233536 audited attempts and failed G1 and G2. The corrected loader "
                 "materializes each component and decides from the component, so the test "
                 "and the result agree by construction; it also counts the representability "
                 "of the velocity rather than of the normalized coordinate in capped mode. "
                 "The implementation is versioned 2.1.0 and 2.0.0 is superseded before "
                 "release, so the two samplers are never confusable by version string.\n"
                 "\n"
                 "(2) PROTOCOL. Amendment 1.3.0's two upper-tail members of family F5 were "
                 "stated against the untruncated law. That is not the law the data obey: "
                 "near kappa = 1/2 the target puts probability outside the floating-point "
                 "type, so honest overflow right-censors the returned tail, at a "
                 "direction-dependent point. Measured over 2000 replicates of a perfectly "
                 "correct loader at the frozen production size, the 1.3.0 excess member "
                 "rejected with probability 1.000 on cases C3 and C4 against a nominal "
                 "0.010. F5's failure in the first holdout was therefore a property of the "
                 "null, not of the candidate. Both members keep their statistics and their "
                 "thresholds and are given the nulls the data obey; see F5_tail. The count "
                 "member's side condition is withdrawn because the quantity it could not "
                 "observe turned out to be recorded already, so that member now applies "
                 "where it used to be declared not applicable -- strictly more tests. "
                 "Negative control NC3 is respecified for the same reason: what it injected "
                 "was honest overflow, which a correct loader is right to produce, so its "
                 "measured power was type-I error carrying a power label.\n"
                 "\n"
                 "(3) SCOPE. Acceptance is limited to the supported environment, arm64 "
                 "macOS under both standard libraries. The cross-architecture claim is "
                 "withdrawn from section 9 and published as a limitation rather than left "
                 "as a permanently open gate on evidence this project cannot obtain.\n"
                 "\n"
                 "Nothing here relaxes a threshold, removes a cell, changes a family's "
                 "alpha, or excludes a result. Every change was made and committed before "
                 "any datum on seeds 8001-8010 existed."},
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
            "description": "released loader cpp/bi_kappa_distribution.H at version 2.2.0, "
                           "radius built in the log domain and carried, with the rest of "
                           "the deterministic map, in the accumulator of "
                           "bikappa_detail::log_accumulator -- double for a float "
                           "instantiation, the working type otherwise -- with one rounding "
                           "at materialization",
            "version": "2.2.0",
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
            "derivation_rule":
                "The highest seed declared anywhere in this repository is 8010 (the second "
                "holdout's performance block). Round up to the next multiple of 1000, "
                "which is 9000, and take the next ten integers: 9001-9005 for production "
                "and 9006-9010 for the performance block. The rule admits exactly one "
                "answer, so the block is a consequence of the repository's state and not a "
                "choice made after seeing a result.",
            "p1_stream_seed_stride": P1_STREAM_SEED_STRIDE,
            "p1_stream_seed_rule":
                "Protocol 3.0.0 only. In phases P1 and P5 the engine of configuration "
                "(precision, kappa) under replicate seed `base` is seeded with "
                "base + 10000 * (13 * precision_index + kappa_index), precision_index 0 "
                "for double and 1 for float and kappa_index the 0-based position in the "
                "kappa ladder. The 26 configurations of a replicate therefore run on "
                "independent streams instead of sharing one, which is what family F2's "
                "global rule requires and what protocol 2.0.0 and earlier did not "
                "provide. The formula is declared, deterministic and admits exactly one "
                "answer per configuration; every value it produces is congruent to "
                "9001-9005 modulo 10000 and so collides with no declared block. It is "
                "recorded per row as `stream_seed` beside the replicate's `seed`, which "
                "remains the unit the analysis groups and clusters by. P2, P3, P4 and P6 "
                "are unchanged: their families decide by Holm, by Simes or by an exact "
                "per-attempt identity, all of which are valid under arbitrary dependence "
                "between configurations, so there is nothing for the change to fix there "
                "and no reason to disturb them.",
            "spent_blocks": {"first_holdout": SEEDS_FIRST_HOLDOUT,
                             "second_holdout": SEEDS_SECOND_HOLDOUT,
                             "spent_reason":
                                 "used by the NO-GO holdouts preserved in commits 45d3ef8 "
                                 "(7001-7010) and e5c9837 (8001-8010); PROTOCOL.md section "
                                 "8 forbids recomputing any result on them, and they may "
                                 "now serve only as preserved failure and development "
                                 "evidence"},
            "disjoint_from": {"exp1": [1001, 1005], "exp2": [2001, 2005],
                              "exp3": [3001, 3003, 3101], "exp4_exp6": [4001, 4010],
                              "exp7_first_holdout": [7001, 7010],
                              "exp7_second_holdout": [8001, 8010],
                              "exp7_selftest": [7501, 7505]},
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
            "nc3_excess_loss_fractions": NC3_EXCESS_LOSS_FRACTIONS,
            "nc3_effect_definition":
                "NC3's effect size is the fraction of ATTEMPTS removed from among the draws "
                "the floating-point type would have allowed the loader to return -- "
                "conditioning in excess of the representability floor, not the floor "
                "itself. Protocol 1.3.0 injected a single direction-independent cutoff into "
                "an uncensored sample, which is what honest overflow does to a CORRECT "
                "loader; the battery 'detected' it only because its null was the "
                "untruncated law, so the measured power was type-I error carrying a power "
                "label. A level row at zero excess is published beside the power rows and "
                "carries no threshold.",
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
            "statistic":
                "Two upper-tail members per uncapped loader cell per threshold q0. COUNT: "
                "the number of ATTEMPTS whose intended Z exceeds z0 = -log q0, against "
                "Binomial(n_attempted, q0), exactly. EXCESS: a Kolmogorov-Smirnov test of "
                "the per-draw conditional probability integral transform of the RETURNED "
                "draws above z0 against Uniform(0,1), where each draw's truncation point is "
                "its own representability threshold C(n) = Z(log max() - log max_j|g_j(n)|).",
            "q0": TAIL_Q0,
            "rationale":
                "Amendment 1.3.0 added an upper-tail statistic because the battery "
                "certifying a heavy-tailed law had nothing in it that looked at the tail. "
                "It stated both members against the untruncated law, which is not the law "
                "the data obey: near kappa = 1/2 the bi-Kappa law puts probability outside "
                "every finite floating-point range, so the loader cannot return the far "
                "tail and is right not to. Amendment 2.0.0 keeps both members and both "
                "statistics and corrects the nulls.",
            "capped_cells_excluded": True,
            "capped_cells_reason":
                "under a cap the accepted radius is truncated at a direction-dependent "
                "bound, so neither the attempt count nor the returned sample is the one "
                "these members are defined on; the cap-law member tests those cells by "
                "their own conditional transform instead",
            "count_member_applies_at_every_threshold": True,
            "count_member_basis":
                "The probe records log_r_ref -- the radius the attempt carried -- for every "
                "attempt, the overflowed ones included, and an uncapped cell runs its core "
                "mapping exactly once per attempt, so its record count equals its attempt "
                "count and the INTENDED Z of every attempt is on disk, uncensored. Against "
                "that sample the count above z0 is Binomial(n_attempted, q0) exactly. "
                "Protocol 1.3.0 computed this member from the radius RECOVERED FROM THE "
                "RETURNED VECTORS, which is censored by representability, and then patched "
                "around the censoring with a side condition "
                "(count_member_requires_q0_above_honest_floor) that declared the member not "
                "applicable wherever it would have misfired. The side condition is "
                "withdrawn: the member now applies at every threshold of every uncapped "
                "cell, which is strictly more tests and not fewer. Family F4 was already "
                "computed this way on the scalar phase, which is why F4 passed the first "
                "holdout while F5 did not.",
            "count_member_requires_q0_above_honest_floor": False,
            "count_member_withdrawn_side_condition":
                "count_member_requires_q0_above_honest_floor was true under 1.3.0 and is "
                "false under 2.0.0. It is withdrawn because the quantity it could not "
                "observe is observable, not because the cells it excluded became "
                "convenient: those cells (C3 and C4 at q0 = 1e-4) are now tested where they "
                "previously were not.",
            "excess_member_applies_at_every_threshold": True,
            "excess_member_null":
                "Z is independent of the direction, so conditional on its own direction a "
                "returned draw above z0 is Exp(1) truncated to (z0, C_i], with "
                "C_i = Z(log max() - log max_j |g_j(n_i)|) computed from that draw's own "
                "recovered direction and the cell's declared theta and ub. Its probability "
                "integral transform U_i = (1 - exp(-(Z_i - z0))) / (1 - exp(-(C_i - z0))) "
                "is i.i.d. Uniform(0,1) under the null. Where no attempt can overflow, C_i "
                "is effectively infinite and U_i reduces to 1 - exp(-(Z_i - z0)); a "
                "Kolmogorov-Smirnov statistic is invariant under a common monotone "
                "transform of the data and the null CDF, so the number returned is "
                "IDENTICAL to amendment 1.3.0's, not merely equivalent. Measured on case "
                "C0 the two agree to 1.1e-16. This is the same construction the protocol "
                "already uses for the capped cells.",
            "excess_member_why_the_frozen_null_was_invalid":
                "Honest overflow removes exactly the largest draws, so the surviving "
                "excesses are right-censored, and at a DIRECTION-DEPENDENT point: a draw is "
                "returned iff R max_j |g_j(n)| <= max(), and max_j |g_j| varies over the "
                "sphere by up to sqrt(3) theta_max / theta_min, so a cutoff inferred from "
                "the cell's total overflow rate is the average of that boundary rather than "
                "the boundary. Measured over 2000 replicates of a perfectly correct loader "
                "at the frozen production size, the 1.3.0 excess member's rejection rate "
                "was 1.000 on cases C3 and C4 against a nominal 0.010, and the 1.3.0 count "
                "member's was 1.000 on C4 at q0 = 1e-4. Those are the three rejections that "
                "failed the first holdout's gate G3.",
            "calibration":
                "Level and power are measured in "
                "docs/revision/experiments/f5_tail_calibration.md on simulated cells only, "
                "before any second-holdout datum existed, at the frozen production size of "
                "500000 attempts per cell.",
            "unresolved_Z_counts_toward_every_threshold": True,
            "unresolved_Z_note":
                "'Unresolved' here means the Z transform itself underflowed, so the draw is "
                "further into the tail than any threshold. It does NOT mean the velocity "
                "overflowed: under 2.0.0 an overflowed attempt carries an ordinary intended "
                "Z, taken from log_r_ref, and is counted by it.",
            "record_count_must_equal_attempt_count": True,
            "record_count_rule":
                "An uncapped cell must produce exactly one record per attempt. A shortfall "
                "means attempts went unrecorded, which is silent conditioning, and it is "
                "reported as a count rather than absorbed into a rate.",
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
                                  "is what the frozen null is built over. The "
                                  "informativeness flag is published as a map per section "
                                  "5.1 and does not change the count; excluding "
                                  "non-informative cells would make F2 unevaluable against "
                                  "its own frozen null.",
            "F2_null": "the Poisson-binomial over independent intervals is withdrawn and "
                       "replaced by the exact convolution over independent CONFIGURATIONS "
                       "of the per-configuration miss-count law, which is computed rather "
                       "than assumed. See F2.why_not_poisson_binomial and "
                       "seeds.p1_stream_seed_rule; the two go together, because the "
                       "independence the new null needs is supplied by the streams and not "
                       "by an assumption about them.",
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
                   "translated_execution_excluded": True,
                   # Protocol 2.0.0 narrows G4 to the environment this work supports.  The
                   # cross-architecture claim is WITHDRAWN, not assumed: it is removed from
                   # the claim boundary in section 9 and published as a limitation. What
                   # remains is the sharper of the two tests and the one that is decidable
                   # here -- bitwise equality of the two standard libraries on the supported
                   # architecture, which is a prediction with no tolerance at all.
                   "scope": "single_architecture",
                   "cross_arch_in_scope": False,
                   "cross_arch_withdrawn_reason":
                       "Native x86_64 hardware is not available to this project, and "
                       "PROTOCOL.md section 5.3 refuses emulation as closure. Rather than "
                       "leave the gate permanently OPEN on evidence that cannot be "
                       "obtained, acceptance is limited to the supported environment and "
                       "the cross-architecture claim is withdrawn from section 9. Contrary "
                       "cross-architecture evidence, if any is ever filed, still FAILS the "
                       "gate: narrowing the claim does not license ignoring a "
                       "disagreement.",
                   "supported_environments": SUPPORTED_ENVIRONMENTS,
                   "requires_every_supported_environment_native": True},
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
