# Experiment 7 source data

Every figure and every number quoted from this experiment is generated from the CSV files
described here. `make_figures.py` reads only these files; it never opens a raw binary or a
counter JSONL, so every plotted point is recoverable from a human-readable row.

This file carries no generation timestamp, and neither does any CSV, `analysis_report.md`,
`validation_matrix.md`, `exp7_results.json` or `figures/figure_manifest.json`. `make
reverify` regenerates every derived artifact and diffs it byte for byte, which a wall-clock
stamp would make impossible. The one file that does record wall-clock time and build
identity is `provenance.md`, and it copies those values out of `raw/manifest.csv` rather
than reading the clock, so it regenerates byte for byte too.

Every file below carries a `protocol_sha256` column holding the SHA-256 of the
`config/protocol.json` the run was judged by. `analyze.py` refuses to run at all unless that
hash matches the one every probe process recorded in `raw/manifest.csv` and in every counter
row: a result that cannot be tied to the rules it was judged by is not a confirmatory
result.

A `--smoke` analysis reads `raw/smoke/` only and writes `results/smoke/` only, with its
figures under `results/smoke/figures/`. The repository excludes
`experiments/*/results/smoke/` wholesale, so a rehearsal can reach neither a production
result file nor the tracked figure set, and the two modes refuse to read each other's data.

---

## 0. Conventions that hold in every file

### Missing values and non-finite values

These are **different statements** and have different codes. A blank cell is never written.

| code | meaning |
|---|---|
| `NA` | the quantity does not apply to this row, or was not measured. Never a zero |
| `nan` | the quantity was computed and came out not-a-number — for example a homogeneity test on a table with an empty margin |
| `inf`, `-inf` | the quantity was computed and is infinite — for example an order-statistic bracket whose upper end is an overflowed draw |
| `true`, `false` | booleans, lower case |

### Units

Every logarithm is a natural logarithm. Every time is in seconds. Every count is a count of
**attempts** unless the column name says `returned`. `kappa` and `shape_a = kappa - 1/2` are
dimensionless; `log_r`, `log_q` and every `log_*` column are in natural-log units;
probabilities and ratios are dimensionless. The cap `lambda` is in thermal speeds.

### Rates and intervals

Every rate in this bundle is written as a block of seven columns, `<name>_count`,
`<name>_n`, `<name>_rate`, `<name>_ci_lo`, `<name>_ci_hi`, `<name>_interval_kind` and
`<name>_is_upper_bound`:

| `<name>_interval_kind` | meaning |
|---|---|
| `two_sided_clopper_pearson_0.95` | a non-zero count: the exact two-sided 95 % Clopper-Pearson interval. `_is_upper_bound` is `false` |
| `one_sided_upper_0.95` | a **zero** count: `_rate` is 0, `_ci_lo` is 0 and `_ci_hi` is the one-sided 95 % upper limit `1 - 0.05^(1/N)`. `_is_upper_bound` is `true` |
| `none_no_trials` | `N = 0`: rate and both ends are `NA`, because nothing was observed either way |

A zero count is **never** reported as a zero-width interval. The two kinds are not the same
bound and the column says which one each entry is, so "no failures observed" can never be
read as "the rate is zero". Figures draw a zero count as a downward triangle at its
one-sided upper limit and never at zero.

### Cluster-bootstrap intervals

Where an interval pools over seeds it is a **cluster bootstrap over the seed block**:
resample the seeds with replacement, then resample within each seed, 10 000 resamples, 99 %
percentile interval, with a recorded RNG seed per interval. The columns are `boot_estimate`,
`boot_lo`, `boot_hi`, `boot_conf`, `boot_resamples`, `boot_rng_seed`, `boot_kind`
(`cluster_over_seeds`), `boot_n_seeds`, `boot_n_per_seed` and `boot_resolved`.

`exp7_stats.seed_stratified_bootstrap` is deliberately **not** used anywhere: its own
docstring records that its variance has no between-seed component, so no interval built on
it can see the seed-level heterogeneity an RNG-stream or portability defect produces.

`boot_n_per_seed` is the resolution of the **within-seed** half of the resample, capped at
2000 draws per seed. The between-seed half — the component that was missing before — is
always exact, because it resamples the five seeds themselves. Capping the inner draw widens
the interval relative to the full sample and never narrows it. The decisive interval for a
failure probability is the Clopper-Pearson one in the `*_ci_*` columns, as the protocol
requires; the bootstrap columns are the pooled-over-seeds companion.

### The analytic honest floor

Every failure rate is reported beside `honest_floor_rate`, the closed-form probability that
the intended vector has **no representation** in the run's floating-point type. It is
derived in `config/honest_floor.md` and implemented in `analyze.py`:

```
log floor = lgamma(3/2 + a) - lgamma(3/2) - lgamma(a + 1)
            - 2 a log(MAX) + a log(kappa) + log E[(max_j |n_j|)^(2a)]
```

with `a = kappa - 1/2`, `MAX` the largest finite value of the type, and the direction
expectation evaluated on a fixed 512x512 midpoint grid over the sphere, so the number is the
same on every run and every host. It is reported twice: `honest_floor_rate` as a
probability, which underflows to 0 below about 1e-308, and `honest_floor_log10` as its base-10
logarithm, which does not. A configuration with no observed failure means nothing without
this number beside it: at double `kappa = 0.55` the floor is 1.5e-31, and at float
`kappa = 0.55` it is 1.4e-4.

`observed_over_floor` is the observed failure rate divided by the floor. It is `NA` where the
floor has underflowed to zero.

### Common identity columns

`tag` is the build that produced the row (`libcxx`, `libstdcxx`); `stdlib`, `compiler` and
`arch` are that build's identity; `execution` is `native` or `translated`, and translated
execution is recorded but excluded from every decision. `precision` is `float` or `double`
and accuracy statistics are **never** pooled across the two. `scope` is `seed` for one
replicate or `pooled` over the seed block. `method` is `LEGACY` (the released 1.0.0
comparator), `CANDIDATE` (2.0.0, the log-domain loader under test) or `QF` (the
quotient-first formation, a diagnostic that carries no claim). `layer` distinguishes the
layers described per file below.

---

## 1. `scalar_validation.csv` — P1, the scalar radial law

One row per `tag` x `method` x `precision` x `kappa` x `seed`, from the P1 `native` layer —
the instrumented replica of the released class, whose stream `make selftest` checks against
the class bit for bit.

| column | meaning | unit |
|---|---|---|
| `phase`, `tag`, `stdlib`, `arch`, `execution`, `layer`, `method`, `precision`, `seed` | identity, as above | -- |
| `kappa`, `shape_a` | spectral index and the denominator shape `kappa - 1/2`, formed in the working precision | dimensionless |
| `n_attempted` | attempts; uncapped, so also intended draws | count |
| `cat_finite` … `cat_overflow_returned_finite` | the ten terminal categories of `results/schema.md`. **In P1 these are the natively observable split only**: with no paired reference a draw whose denominator materialized as zero cannot be separated into avoidable and honest. `failure_envelope.csv`'s `paired` rows are the authoritative decomposition | count |
| `n_finite`, `n_avoidable`, `n_honest` | the three category sums | count |
| `nonfinite_output` | returned vectors with a non-finite component — the failure count | count |
| `x2_zero`, `x2_subnormal` | denominator states; LEGACY only, `0` for CANDIDATE, which never forms `X2` | count |
| `gamma_variates`, `uniform_variates`, `engine_calls` | variates consumed and `mt19937` invocations | count |
| `seconds`, `seconds_per_attempt` | wall time for the configuration | s |
| `max_finite_log_r` | largest finite `log R` seen | log units |
| `n_resolved`, `n_nonresolved` | draws whose `log R` is finite, and those lost to overflow. Unresolved draws are **not dropped**: they are sorted to `+inf`, so the order statistics are those of the whole sample | count |
| `loss_fraction` | `n_nonresolved / n_attempted` | probability |
| `conditional` | `true` when `n_nonresolved > 0`. Such a cell cannot by itself support the fidelity claim | -- |
| `failure_*`, `avoidable_*`, `honest_*` | rate blocks (see §0) for `nonfinite_output`, `n_avoidable` and `n_honest` over `n_attempted` | probability |
| `honest_floor_rate`, `honest_floor_log10`, `observed_over_floor` | the analytic floor (see §0) | probability, log10, ratio |
| `ad_statistic`, `ks_statistic`, `cvm_statistic` | exact Anderson-Darling `A^2`, Kolmogorov-Smirnov `D` and Cramer-von Mises `W^2` of `Z = -log I_W(a, 3/2)` against `Exp(1)`, computed by the probe with the estimators `results/schema.md` fixes. **These three are the whole of family F1** — `corrections.F1_contents`. The Beta test on `W` that PROTOCOL.md §5's table lists beside them is withdrawn: the probe emits no `W` sample, `Z` is a monotone transform of `log W`, and G0 already checks the `log q` evaluator against an arbitrary-precision incomplete beta to `1e-10`, which is a stronger check on the transform than comparing the two routes to each other | dimensionless |
| `ad_pvalue`, `ks_pvalue`, `cvm_pvalue` | their upper tails under the fully specified null. AD uses the Marsaglia & Marsaglia (2004) limiting CDF, KS `scipy.stats.kstwo.sf(D, n_resolved)`, CvM the Csorgo & Faraway limiting series | probability |
| `f1_cell_simes_p` | the Simes combination of the three, the configuration's F1 statistic | probability |
| `tail_records`, `tail_file` | the companion `TailRecord` stream and its length. The path is relative to the experiment directory | count, path |
| `class_n_finite`, `class_replica_agrees` | the released class's own count for the same configuration, and whether the instrumented replica reproduced it | count, -- |

### Quantile blocks `q0p5_*`, `q0p9_*`, `q0p99_*`, `q0p999_*`, `q0p9999_*`

One block per protocol level `p`.

| suffix | meaning |
|---|---|
| `_lo_index`, `_hi_index` | the 0-based order-statistic bracket frozen in `config/protocol.json` for the production sample size |
| `_achieved_coverage` | the protocol's analytically computed coverage of that bracket (nominal 99 %: 95 % Bonferroni-corrected over the five levels) |
| `_lo_value`, `_hi_value` | the sample `log R` at those indices. `nan` when the sample cannot fill the bracket; may be `inf` |
| `_point_value` | the sample `log R` at `ceil(p n) - 1`, the empirical quantile |
| `_target` | the exact target quantile of `log R`, obtained in the log domain and never from `R^2` |
| `_error` | `_point_value - _target` |
| `_width` | `_hi_value - _lo_value`, in natural-log units |
| `_resolved` | `true` when both bracket ends are finite. A cell that did not resolve enters neither the miss count nor the informativeness map |
| `_informative` | `true` when `_width <= 1` natural-log unit, the threshold frozen in `config/protocol.json`. **A non-informative cell is not evidence**: at `kappa = 0.5001` the bracket runs to 8951 natural-log units, and an interval admitting a multiplicative error of `e^8951` cannot distinguish a correct sampler from a wrong one. Those shapes are certified by the families that retain power at any shape. The flag is published as a map; it is computed from the frozen width threshold and not from the data |
| `_covered` | `true` when `_lo_value <= _target <= _hi_value`. `NA` when the cell did not resolve. Family F2's miss count is taken over every **resolved** interval — `corrections.F2_informativeness` — because that is the set the frozen Poisson-binomial null was built over. The informativeness flag is published as a map and changes no count; dropping the non-informative cells from the count would leave F2 with no null to be tested against |

### Tail blocks `tail0p01_*`, `tail0p001_*`, `tail0p0001_*`

One block per protocol `q0`.

| suffix | meaning |
|---|---|
| `_z0` | `-log q0` |
| `_observed_resolved` | resolved draws with `Z > z0` |
| `_observed_nonresolved` | unresolved draws, which exceed every threshold; the same number in all three blocks |
| `_observed_total` | `_observed_resolved + _observed_nonresolved`, **the count family F4 tests** |
| `_n_test` | the denominator it is tested against, which is `n_attempted` |
| `_expected_resolved` | `n_resolved * q0` |
| `_expected_per_attempt` | `n_attempted * q0`, the null expectation of `_observed_total` |
| `_n_excess` | resolved draws the excess test was computed from |
| `_p_count` | exact binomial p-value of `_observed_total` against `Binomial(n_attempted, q0)` |
| `_p_count_resolved_only` | what the same test would have said on `_observed_resolved` against `Binomial(n_resolved, q0)`. Published bookkeeping; no decision reads it |
| `_p_excess` | Kolmogorov-Smirnov p-value of the resolved excesses `Z - z0` against `Exp(1)` |

**The exceedance count includes the unresolved draws, and its denominator is the attempts.**
`config/protocol.json -> F5_tail.unresolved_exceedances_count_toward_every_threshold` is
true and amendment 1.3.0 gives the reason: a draw whose `Z` never resolved lies above every
threshold by construction, so testing only the resolved subset conditions on resolvability,
which is itself monotone in the tail — the exact bias the statistic exists to detect. The
rule is written for family F5's new upper-tail members and is applied to F4 as well, because
F4 is the same statistic on the scalar phase.

Two consequences are worth stating rather than leaving to be discovered. Where a cell's loss
fraction is at or below `q0`, the count is exactly unbiased: the unresolved draws replace
precisely the resolved exceedances they displaced. Where the loss fraction **exceeds** `q0`,
it is not: the lost draws then lie above the cell's own overflow threshold but not
necessarily above `z0`, so the count runs high and the test reads on the loss fraction
itself. `_p_count_resolved_only` and `_observed_nonresolved` are published so that a reader
can tell the two situations apart on any row.

---

## 2. `scalar_ecdf.csv` — the `Z` ECDF on the frozen grid

One row per `tag` x `method` x `precision` x `kappa` x grid point, pooled over seeds.
501 grid points, step 0.05.

| column | meaning | unit |
|---|---|---|
| `n_seeds`, `n_resolved_pooled` | seeds and resolved draws behind the row | count |
| `grid_index`, `z` | grid index `i` and `z = 0.05 i` | --, log units |
| `null_cdf` | `1 - exp(-z)`, the unit-exponential CDF at `z` | probability |
| `ecdf` | cumulative fraction of resolved draws with `Z <= z` | probability |
| `ecdf_residual` | `ecdf - null_cdf` | dimensionless |
| `band_lo`, `band_hi`, `band_alpha` | simultaneous Kolmogorov band at the frozen F1 familywise level | dimensionless |

---

## 3. `honest_floor.csv` — the analytic floor across the frozen ladder

One row per `precision` x `kappa`. A reference curve, not a measurement: it has no interval.

| column | meaning | unit |
|---|---|---|
| `honest_floor_rate`, `honest_floor_log10` | the floor, and its base-10 logarithm for the values that underflow | probability, log10 |
| `type_max`, `log_type_max` | the largest finite value of the type, and its logarithm | -- |
| `configuration` | the configuration the floor is derived for: isotropic, unrotated, `theta_perp = theta_par = 1`, which is what P1 and P2 run | -- |

---

## 4. `failure_envelope.csv` — P2, the mechanism decomposition

`layer` distinguishes the two layers of the phase, and the distinction is load-bearing:

| `layer` | meaning |
|---|---|
| `paired` | all three formations receive the **same** declared primitives, so "was this draw recoverable?" is decided per draw. This is the authoritative decomposition and what G1 and G2 read |
| `native` | each method runs its own stream, as a user would get. Same seed and stream as P1's `native` row for the same configuration, so the two must agree exactly |

`scope` is `seed` for one replicate or `pooled` over the seed block. `diagnostic_only` is
`true` for the `QF` rows: no claim rests on them.

| column | meaning | unit |
|---|---|---|
| `cat_*`, `n_finite`, `n_avoidable`, `n_honest`, `nonfinite_output` | the terminal-category decomposition | count |
| `candidate_avoidable_count` | the candidate's avoidable-loss count for this configuration, `NA` on a row that is not the candidate's. G1 requires the total to be **exactly zero**, with no tolerance proportional to a data-dependent count | count |
| `accounting_ok` | the probe's own evaluation of `n_attempted = n_finite + n_avoidable + n_honest` | -- |
| `accounting_residual` | the same identity recomputed here: `n_attempted - (finite + avoidable + honest)`. Zero when the identity holds exactly. Both are reported so that the check is not taken on trust | count |
| `ref_nonrepresentable` | attempts whose intended vector has no representation in the type — the oracle floor of this configuration | count |
| `honest_minus_ref` | the candidate's honest count minus that floor. G2 requires it to be zero, or each discrepancy individually adjudicated and listed | count |
| `rotation_recoverable`, `near_limit` | attempts whose local vector overflowed while the rotated one is representable, and attempts within the audit margin of a type limit | count |
| `x2_zero`, `x2_subnormal` | denominator states of the shared primitives; identical across the three paired rows | count |
| `n_rel_err`, `mean_rel_err_radius`, `max_rel_err_radius` | the relative radius error against the reference, **per precision**. Never pool the two: the `float` and `double` worst cases differ by orders of magnitude that come from the type and not from the method | dimensionless |
| `max_rel_error_threshold` | the FINITE_BUT_WRONG threshold for the type, from the protocol: 1.05e-8 in `double`, 2.44e-4 in `float` | dimensionless |
| `max_log_r_ref`, `max_finite_log_component` | largest intended `log R`, and largest intended `max_j log|V_j|` among successes | log units |
| `failure_*`, `avoidable_*`, `honest_*` | rate blocks (see §0) | probability |
| `honest_floor_rate`, `honest_floor_log10`, `observed_over_floor` | the analytic floor (see §0) | -- |
| `audit_margin_log_units` | 20.0, the margin within which an attempt is decision-relevant | log units |
| `audit_margin_assertions`, `audit_margin_assertion_failures` | attempts declared unambiguous, each of which had its margin evaluated in working precision, and how many failed the inequality. **Must be zero** | count |
| `audit_min_unambiguous_margin` | the smallest margin any of them cleared; `inf` when there were none | log units |
| `audited` | records this configuration contributed to the audit stream | count |
| `audit_<stratum>_total`, `_audited`, `_rate` | attempts that fell in each of the seven strata, how many were written, and the protocol's declared rate. Assignment is by priority and the first match wins, so the counts partition the sample | count, count, probability |
| `seed_homogeneity_chi2`, `seed_homogeneity_p` | exact-table homogeneity of the failure rate across seeds, on `pooled` rows. A required diagnostic, reported with its p-value | dimensionless, probability |
| `cross_phase_native_agrees` | on a `native` row, whether P1's row for the same configuration, seed and stream reports the same counters. `NA` on a `paired` row | -- |
| `boot_*` | cluster-bootstrap interval for the pooled failure rate (see §0) | -- |

---

## 5. `conditioning_bins.csv` — P3, success against the intended state

One row per `tag` x `method` x `precision` x `kappa` x seed x bin, and a `pooled` row per bin.
The bins are fixed decades of the **intended** upper-tail probability of the radius,
`q = Pr(R > r) = I_W(a, 3/2)`, computed from the per-attempt `CondRecord` stream.

| column | meaning | unit |
|---|---|---|
| `bin_kind` | `target_upper_tail_q` | -- |
| `bin_index`, `bin_lo`, `bin_hi` | the bin and its edges in `q`; the last bin runs down to 0 | --, probability |
| `n_in_bin`, `n_success` | attempts in the bin, and those that delivered a vector | count |
| `success_definition` | what counts as a success: a three-vector reached the caller, i.e. terminal category `finite`, `finite_but_wrong` or `overflow_returned_finite`. A draw whose radius was inaccurate still entered the returned sample, and a tail-retention statement has to count it | -- |
| `success_rate`, `ci_lo`, `ci_hi`, `interval_kind`, `is_upper_bound` | the rate and its interval, per §0 | probability |

---

## 6. `tail_metrics.csv` — P3, the tail consequence

One `pooled` row per `tag` x `method` x `precision` x `kappa` x `p`, for the protocol tail
levels. **Every ratio here says explicitly whether it is per attempt or conditional on
success**, because they are different quantities and the difference is the result.

| column | meaning | unit |
|---|---|---|
| `p`, `log_r_target` | the target percentile of the intended radius and its exact `log R` quantile | --, log units |
| `n_attempts`, `n_success` | attempts, and those that delivered a vector | count |
| `success_definition` | as in §5 | -- |
| `n_above_target_all_attempts` | attempts whose **intended** radius exceeded the target, whether or not they succeeded | count |
| `n_above_target_and_returned` | of those, how many delivered a vector | count |
| `failure_atom_mass_per_attempt` | `1 - n_success / n_attempts`: the mass the per-attempt law puts on failure | probability |
| `per_attempt_law_has_failure_atom` | `true` when that mass is non-zero | -- |
| `per_attempt_quantile_reported` | always `false`, with the reason in the next column | -- |
| `per_attempt_quantile_reason` | a defective per-attempt law is a sub-probability law with an atom on failure, so its upper quantiles are undefined above `1 - atom`. The tail **mass** above the target quantile is reported instead of a quantile error | -- |
| `tail_mass_per_attempt` (+`_lo`, `_hi`) | `n_above_target_and_returned / n_attempts`: the returned tail mass **per attempt**, in which a failed attempt stays outside the returned mass | probability |
| `returned_tail_ratio_per_attempt` (+`_lo`, `_hi`) | the same divided by the target mass `1 - p`. 1.0 is fidelity | dimensionless |
| `tail_mass_conditional_on_success` (+`_lo`, `_hi`) | `n_above_target_and_returned / n_success`: the tail mass **conditional on success**, which is the law a caller actually receives. For independent attempts it is also the law produced by silently redrawing until success | probability |
| `returned_tail_ratio_conditional_on_success` (+`_lo`, `_hi`) | the same divided by `1 - p` | dimensionless |
| `boot_conf`, `boot_resamples`, `boot_kind`, `boot_n_seeds`, `boot_n_per_seed` | the interval definition behind every `_lo`/`_hi` above (see §0) | -- |

The word "returned" never appears in this file without one of the two qualifiers.

---

## 7. `audit_coverage.csv` — adjudication coverage, P2 and P3 together

One row per audit stratum.

| column | meaning | unit |
|---|---|---|
| `stratum` | one of the seven strata of `config/protocol.json`, in priority order | -- |
| `declared_rate` | the rate the protocol declares for it | probability |
| `attempts_in_stratum`, `attempts_audited` | the stratum's size and how many attempts were written to the audit stream | count |
| `achieved_rate` | `attempts_audited / attempts_in_stratum`; `nan` for an empty stratum, where coverage is undefined rather than zero | probability |
| `meets_declared_rate` | whether the achieved rate reaches the declared one; `NA` for an empty stratum | -- |
| `phases` | `p2+p3` | -- |

---

## 8. `loader_validation.csv` — P4, the complete three-dimensional loader

One row per `tag` x `case` x `method` x seed (`scope = seed`) and one pooled over the seed
block (`scope = pooled`). The statistics are computed from the returned vectors alone, by
recovering `(log R, unit direction)` through a field-aligned basis **re-derived independently
of the loader's**, so the frame test is a genuine check of the rotation rather than a
tautology.

| column | meaning | unit |
|---|---|---|
| `case`, `role` | `C0`…`C6` and the role PROTOCOL.md §4.1 gives it | -- |
| `theta_ratio`, `ub_x`, `ub_y`, `ub_z`, `cap` | `theta_par/theta_perp` with `theta_perp = 1`, the field direction, and the cap in thermal speeds. `cap` is `NA` when uncapped | -- |
| `layer` | `uncapped`, `capped` | -- |
| `attempts`, `n_returned`, `n_analyzed` | attempts made, samples returned, and samples from which a radius and direction could be recovered | count |
| `nonfinite_attempt`, `nonfinite_returned` | attempts whose vector had a non-finite component, and how many of the returned samples did | count |
| `cap_reject`, `cap_reject_unrepresentable` | attempts the box rejected, and how many of those were also unrepresentable — which tells "the box rejected a representable sample" from "the box rejected one that could not have been returned at all" | count |
| `cap_exhausted`, `max_attempt_run` | requests that hit the internal attempt limit, and the longest rejection run behind a single return | count |
| `engine_calls`, `seconds` | `mt19937` invocations and wall time | count, s |
| `conditional` | **`true` when the statistics ran on a subsample that something correlated with the radius had already filtered.** A non-finite component is exactly what a large radius produces, so any such cell is conditional. A conditional cell cannot support the fidelity claim on its own | -- |
| `loss_fraction` | the fraction lost, always present on a conditional cell | probability |
| `loss_kind` | which fraction it is: non-finite attempts over attempts for an uncapped cell, unrecoverable returns over returns for a capped one | -- |
| `nonfinite_attempt_rate_*` | the rate block (see §0) for `nonfinite_attempt` over `attempts` | probability |
| `honest_floor_rate`, `honest_floor_log10` | the analytic floor for the cell's `kappa` and precision | -- |
| `in_family` | `true` on the cells that enter the F5 decision: the CANDIDATE cells of the primary environment, pooled over seeds. Within version 2.0.0 the stream is a function of the engine alone, so the other environment's rows are the same draws and would enter Holm twice. LEGACY rows are reported with their p-values but are not gated | -- |
| `seed_homogeneity_chi2`, `seed_homogeneity_p` | homogeneity of the non-finite rate across seeds, on pooled rows | -- |

### The tests of family F5

For each test `T` there are four columns: `p_T`, `stat_T`, `applies_T` and
`f5_holm_rejected_T`. The last is the family's own Holm decision, taken once over every cell
and every test jointly and written into the row, so that no reader and no figure has to
re-derive a rejection threshold from an alpha and a count of tests. It is `false` on a
`scope = seed` row and on a cell that is not `in_family`, neither of which is gated: the
decision is taken on the pooled candidate cells of the primary environment.

The first five tests were frozen with the protocol. The upper-tail members were added by
amendment 1.3.0, one pair per threshold in `config/protocol.json -> F5_tail.q0`, and named
after it: `tail_1em4_count` is the exceedance count at `q0 = 1e-4`. They **add** tests to
the family rather than replacing any of it.

| test | statistic | applies to |
|---|---|---|
| `direction_uniformity` | KS and Cramer-von Mises of `cos theta` in the field-aligned frame against `U(-1, 1)`, Simes-combined | uncapped cells |
| `independence` | chi-square on equal-count bins of `log R` against equal-width bins of `cos theta` — rank-based in the radius, whose scale spans hundreds of log units | uncapped cells |
| `frame_invariance` | KS and CvM of the azimuth in the re-derived field-aligned frame against `U(-pi, pi)`, Simes-combined | uncapped cells |
| `anisotropy` | two-sample KS between `abs(n_3)` and `abs(n_1)` after dividing by the declared scale. After descaling, the three components of an isotropic direction are exchangeable, so a mis-applied `theta` ratio breaks the equality. Moment-free on purpose: for `kappa <= 3/2` the loaded population has no finite variance | every cell |
| `cap_law` | KS against `U(0,1)` of `F(R)/F(R_max(n))`, where `log R_max = log(lambda) - log(sqrt(kappa)) - log(max_j abs(n_j))`. Given the direction, the accepted radius is the target truncated at `R_max`, so the transform is exactly uniform under the declared bounded law | capped cells |
| `tail_<q0>_count` | exceedance count above `z0 = -log q0` on the cell's recovered radial `Z`, the unresolved draws included, against `Binomial(attempts, q0)`. `stat_T` is the count itself | uncapped cells |
| `tail_<q0>_excess` | KS of the resolved excesses `Z - z0` against `Exp(1)`. `stat_T` is how many excesses it was computed from; fewer than eight leaves the member `NA` | uncapped cells |

`applies_T = false` means the test is **not a property of that cell's law**, not that it
passed. Under a cap the accepted set couples the radius to the direction and is not
rotationally symmetric in the azimuth, so direction uniformity, independence and frame
invariance do not apply there; the accepted radius is truncated at a direction-dependent
bound, so the count above `z0` is not `Binomial(n, q0)` and the upper-tail members do not
apply there either; the accepted set is symmetric under permuting the three
components, so the anisotropy test survives the cap. `p_T` and `stat_T` are `NA` there.

### The parts of every exceedance count

One block per threshold, on the same convention as §1's tail blocks.

| column | meaning |
|---|---|
| `tail_<q0>_observed_resolved` | recovered draws with `Z > z0` |
| `tail_<q0>_unresolved` | draws the cell attempted and the analysis could not resolve, each of which is above every threshold |
| `tail_<q0>_expected` | `attempts * q0` |
| `tail_<q0>_n` | the denominator, which is `attempts` |
| `tail_<q0>_p_count_resolved_only` | what the count would have said on the resolved draws alone. Published bookkeeping; no decision reads it |

A cell whose loss fraction exceeds `q0` is read by the count member as an excess, because
its unresolved draws lie above its own overflow threshold and not necessarily above `z0`.
The two columns above make that visible on the row rather than leaving it to be inferred
from a p-value.

`validation_matrix.md` is the same information as a table, with the Holm decisions, the
`n.a.` cells and every conditional cell's loss fraction shown.

---

## 9. `negative_controls.csv` — measured power, not a pass/fail bit

One row per control per effect size. The detection rule of each control is the **frozen
family rule itself**, evaluated on the injected sample: a control counts as detected when the
family that certifies the corresponding property fails. That keeps the measured power a
property of the battery rather than of a statistic chosen after the fact.

| column | meaning | unit |
|---|---|---|
| `control` | `NC1_radius_direction_coupling`, `NC2_capped_vs_uncapped_weak_cap`, `NC3_survivor_conditioning` | -- |
| `effect_kind`, `effect_size` | `loss_fraction` and its value, `cap_lambda` and the cap, or `level` with effect `0`, which is the same procedure with nothing injected | -- |
| `detection_statistic` | the family rule used as the detector | -- |
| `injection_definition` | exactly what was injected | -- |
| `source_cell` | the measured sample the injection was applied to, and its size | -- |
| `injections`, `detections` | replicates run, and how many were detected. `injections = 0` means the control could not be run and the row is a failure, never a vacuous pass | count |
| `power`, `power_ci_lo`, `power_ci_hi`, `power_interval_kind`, `power_is_upper_bound` | the measured power as a rate block (see §0) | probability |
| `required_power`, `passed` | the pre-registered threshold, and whether the Clopper-Pearson **lower** limit clears it. Both are `NA` on a `level` row: a rejection rate measured with nothing injected is not a power and is neither passed nor failed by the power rule | -- |
| `in_family` | `true` on the rows that enter family F6's decision. `false` rows are published diagnostics: they state the detection floor rather than leaving it implied | -- |
| `note` | free text, `NA` when there is nothing to say | -- |

### How a replicate is drawn, and why it is not a bootstrap

A power is a rejection rate under an alternative, so the replicates have to be drawn in a
way that leaves the measurement's own null exact. **None of them resamples a measured sample
with replacement.** The empirical joint distribution of a finite sample carries a real
dependence of the order of its own chi-square, so a with-replacement resample inherits that
as a genuine non-centrality: the independence member of F5 then rejects at close to
certainty with nothing injected at all, at any sample size, and the "power" it reports is
the resampling scheme's and not the battery's. NC1 therefore permutes the measured cell's
radii against its directions, which is the independence null carrying the cell's own two
marginals and uses every draw exactly once; NC3 simulates the law outright. The `level` row
published for NC1 measures the rejection rate with nothing injected, so every power figure
here can be read against a stated false-positive rate rather than an assumed one.

Injections:

* **NC1** removes the `q` fraction of draws with the largest `R max_j abs(n_j)`. That is the
  radius-direction coupling a finite-precision loss actually creates: the draws a loader
  cannot return are the ones whose largest component overflows, which is a joint property of
  the radius and the direction and of neither alone. Two readings are measured and both are
  published. The in-family row is the defect the control exists for: a loader that
  **silently** conditions, so the cell it hands the analysis reports no loss and the
  battery has to catch it unaided. The row beside it declares the same removed draws as
  unresolved attempts — the reading that applies when the loader reports its loss, which
  then counts toward every threshold and leaves the cell labelled conditional with its loss
  fraction under PROTOCOL.md §5.2.
* **NC2** tests the bounded sample against the **uncapped** law at the weakest declared cap,
  which is the smallest version of this effect the protocol declares. Its replicates are
  still drawn from the measured capped sample with replacement, which inflates the null of
  the three F1 statistics; its `note` says so, and its detection rate is therefore an upper
  bound on the power at that cap.
* **NC3** simulates survivor conditioning exactly: the finite survivors of a loader that
  loses a fraction `q` are the target conditioned on `Z < -log q`, so the effect is drawn from
  the null rather than approximated. Directions are drawn uniformly on the sphere and
  independently of the radius, which is what a correct loader produces; the conditioning acts
  on the radius alone, so the direction members of F5 have no power against it by
  construction and only the upper-tail members can respond. The in-family effect size is
  fixed by `config/protocol.json -> corrections.NC3_effect`: the loss fraction of the worst
  cell among those the fidelity claim rests on — the **uncapped loader cells** — and not the
  global maximum over the scalar ladder, which is float `kappa = 0.5001`, where almost every
  draw is unrepresentable and detection is trivial. The correction adds "whose loss is not
  dominated by honest overflow"; G1 requires the candidate's avoidable loss to be exactly
  zero, so every loss it has *is* honest overflow and reading that clause as a filter would
  empty the set. It is read here as it is illustrated there: excluding the degenerate
  configurations outside the uncapped loader cases, not excluding a case from among them.
  Each cell's honest floor sits beside its loss fraction in `loader_validation.csv`, so the
  reading can be checked. The diagnostic rows use the pre-registered injection fractions.

---

## 10. `portability.csv` — P5 and the G4 matrix

Three kinds of row, in the `kind` column. **An environment that did not run is a row with
`available = false`, never an absent row and never a zero**, so a missing platform cannot be
read as a passing one.

| column | meaning |
|---|---|
| `kind` | `environment`, `cross_stdlib_bitwise`, or `cross_arch_equivalence` |
| `label` | what the row compares, or the environment it describes |
| `arch`, `os`, `stdlib`, `compiler`, `tag` | the environment's identity; `NA` where unknown because the environment did not run |
| `execution` | `native` or `translated`. Translated execution — Rosetta 2, or Docker with the Rosetta binfmt handler — is recorded as corroborating evidence and excluded from the decision. `NA` for an environment that did not run |
| `available` | whether this host or bundle provides the environment at all |
| `completed` | whether it produced counter rows |
| `reason_not_run` | why it did not, in words, `NA` when it did. For an absent environment this names `results/portability_remote.md`, which carries the exact command |
| `n_rows` | counter rows behind the row | 
| `prediction` | what the protocol predicts for the comparison: bitwise equality for the candidate across standard libraries within one architecture, a difference for the comparator, equality of rates across architectures |
| `n_common`, `n_digests_compared`, `n_differences`, `identical` | the bitwise comparison: configurations present in both arms, digests compared, and differences found. There is no tolerance, because the claim is exact |
| `log_ratio`, `se`, `margin`, `p_value`, `equivalent`, `disagrees`, `informative` | the cross-architecture two one-sided tests on the log rate ratio against the pre-registered margin. `equivalent` and `disagrees` are different outcomes, and a cell that is neither is **underpowered**: it leaves the gate open rather than closing or failing it, which `informative = false` records |
| `in_decision` | whether the row enters the F7/G4 decision |

`analyze.py --portability-ingest <DIR>` writes this file and `g4_verdict.json` from a CI
bundle, and exits non-zero unless G4 closes.

---

## 11. `performance.csv` — P6, cost

`scope` distinguishes four kinds of row, and **the block-level timings are preserved beside
the summary** so that the summary can be recomputed from them.

| `scope` | meaning |
|---|---|
| `block` | one timed block of one method on one case: `block`, `slot` (the order the method ran within the block), `candidate_first` (the realized randomization, recorded so it can be checked rather than asserted), `seed`, `n_returned`, `n_finite`, `n_nonfinite`, `attempts_reported`, `seconds`, `seconds_per_returned` |
| `case_ratio` | the paired candidate/LEGACY ratio of time per returned sample for one case: `ratio_candidate_over_legacy`, `ratio_lo`, `ratio_hi`, `ratio_conf`, `boot_kind`, `boot_resamples`, `boot_rng_seed`, `boot_n_seeds`, `pre_registered_bound` |
| `overall_ratio` | the same pooled over the benchmark cases; this is the G5 statistic |
| `reproducibility` | `repeat_identical`: whether two same-seed runs produced identical counters |

The ratio is paired **within a block**, so it is not affected by drift between blocks; the
interval is a cluster bootstrap over the performance seed block (see §0). Its level is
**0.99**, the value in `config/protocol.json -> statistics.bootstrap.conf`, and every row
records the level it used in `ratio_conf`: PROTOCOL.md §6 says 95 per cent in prose, the
machine-readable value governs — `corrections.G5_interval` — and it is the conservative
choice for a gate on an upper limit. `attempts_reported`
is `0` for LEGACY, which has no attempt counter — a genuine zero, not a missing value.

---

## 12. `log_q_evaluator_validation.csv` — G0

The frozen two-piece `log q` evaluator against an arbitrary-precision incomplete beta, on
both sides of the switch, for every shape on the ladder.

| column | meaning | unit |
|---|---|---|
| `shape_a` | `kappa - 1/2` | dimensionless |
| `log_w_min`, `log_w_max`, `points` | the declared grid | log units, count |
| `max_abs_log_q_error`, `worst_at_log_w` | the worst absolute disagreement and where it occurred | log units |
| `tolerance`, `pass` | the pre-registered tolerance from `config/protocol.json`, and whether the shape cleared it | -- |
| `oracle` | `mpmath.betainc at 60 decimal digits` | -- |

---

## 13. The non-CSV outputs

| file | contents |
|---|---|
| `analysis_report.md` | the verdict on the first line, then the gate table, the family table, what the verdict means, the unresolved items, and the decisive source data. No timestamp |
| `validation_matrix.md` | the F5 matrix (see §8) |
| `exp7_results.json` | the machine-readable verdict: gates, families, the evidence dictionary the gates were evaluated on, the F2 miss count and informativeness totals, the F5 test list, the F6 controls and the performance ratio. No timestamp |
| `provenance.md` | the only file here that records wall-clock time and build identity, all of it copied from `raw/manifest.csv` and the environment record rather than read from the clock, so it regenerates byte for byte. Builds, compile lines, probe hashes, the adjudicated audit streams, and the SHA-256 of every analysis input |
| `g6_evidence.json` | **an input, not an output.** `make verify` and `make reverify` run *after* `analyze.py`, so the one thing the analysis cannot observe is whether its own output re-verified. It reads that from this receipt, written by the release step, with the keys `make_verify_exit_code`, `make_reverify_identical`, `dependencies_clean`, `baseline_comparison_resolved` and `archive_identifier`. Absent, those keys stay unset and G6 **fails**, naming what is missing -- which is the honest state of a run that has not been re-verified yet. The production sequence is `make analyze` -> `make checksums verify reverify` -> write the receipt -> `make analyze` again. `analyze.py` verifies the two checksum manifests itself by recomputing every hash, so `make_verify_exit_code` is a computation over file contents even when no receipt is supplied, not a check that a file exists |
| `g4_verdict.json` | written by `--portability-ingest`: whether G4 closes, and if not, each reason |
| `figures/captions.md`, `figures/figure_manifest.json` | the captions, and one manifest object per figure: panels, core conclusion, source CSVs with their SHA-256s, the script's SHA-256, the size in millimetres and the interval definition. No generation timestamp |

---

## 14. What cannot be read out of these files

* **A rate of zero.** There is none. A configuration with no observed failure is "no failures
  observed in N draws under the tested configuration", with its one-sided upper limit and its
  analytic floor beside it.
* **A pass from an absent measurement.** Every gate states the evidence it requires and fails
  when that evidence is missing. A phase that did not run is an error naming the phase, and
  no verdict is written at all.
* **A tail statement without its conditioning.** Every tail ratio says whether it is per
  attempt or conditional on success, and every loader cell says whether its statistics ran on
  a radius-filtered subsample.
* **An accuracy figure pooled across precisions.** Every accuracy column carries its
  precision, and the two are never combined.
