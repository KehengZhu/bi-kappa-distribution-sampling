# Experiment 6 source data

Generated 2026-09-19T17:46:08+00:00.

Every figure and table in this experiment is generated from the CSV files described here, and every plotted point is recoverable from a row in one of them. `make_figures.py` reads only these files; it never touches the bulk binaries under `raw/`.

Missing values are the literal string `NA`. A blank cell never appears, so a blank can never be read as a zero.

Rates are per attempt unless a column says otherwise. Two-sided intervals are 95% Clopper-Pearson; a zero count is reported as the one-sided 95% upper limit `1 - 0.05^(1/N)` and flagged with an `is_upper_bound` column, never as a zero probability. Tail ratios and method differences carry 99% seed-stratified bootstrap intervals from 10000 resamples. Quantile intervals are exact binomial order-statistic brackets, Bonferroni-corrected within each configuration.

## `failure_envelope.csv`

| column | meaning | unit |
|---|---|---|
| `phase` | experiment phase that produced the row (E2 or E5) | -- |
| `scope` | `pooled` over seeds, or `seed` for one replicate | -- |
| `layer` | `paired` (common declared primitives) or `native` (each method's own primitive) | -- |
| `method` | QF, SPLIT or LOG | -- |
| `precision` | float or double; a real pipeline in that type, never drawn in double and cast | -- |
| `stdlib` | standard library whose Gamma generator was used | -- |
| `kappa` | spectral index | dimensionless |
| `shape_a` | denominator Gamma shape kappa-1/2, formed in the working precision | dimensionless |
| `n_attempted` | attempts; for uncapped runs this equals intended draws | count |
| `n_finite` | attempts returning a finite three-vector | count |
| `n_avoidable` | attempts lost to arithmetic whose intended vector was representable | count |
| `n_honest_overflow` | attempts whose intended vector is not representable in the type | count |
| `failure_rate / avoidable_rate / honest_rate` | per-attempt probabilities | probability |
| `*_ci_lo / *_ci_hi` | two-sided 95% Clopper-Pearson interval | probability |
| `*_is_upper_bound` | true when the count was zero and the interval is the one-sided 95% rule-of-three bound | -- |
| `accounting_ok` | true when finite + avoidable + honest equals attempts exactly | -- |
| `rotation_recoverable` | draws whose pre-rotation local vector overflows while every returned component is representable | count |
| `audit_*_total / audit_*_audited` | size of each audit stratum and how much of it the 100-digit oracle adjudicated | count |

## `conditioning_bins.csv`

| column | meaning | unit |
|---|---|---|
| `bin_kind` | `target_upper_tail_q` or `log_radius` | -- |
| `bin_hi / bin_lo` | bin edges; for q bins these are upper-tail probabilities | probability / log units |
| `n_in_bin / n_success` | attempts in the bin and those returning a finite vector | count |
| `success_rate` | Pr(success | bin) | probability |
| `ci_lo / ci_hi` | two-sided 95% binomial interval | probability |
| `merged` | true when the bin holds fewer than 20 attempts and is not interpreted | -- |

## `tail_metrics.csv`

| column | meaning | unit |
|---|---|---|
| `p` | target percentile | probability |
| `target_log_r_quantile` | exact log R quantile of the uncapped target at this p | log units |
| `rho_attempt` | Pr(R > r_p and success) / (1-p); failures stay outside the returned mass | ratio |
| `rho_cond` | Pr(R > r_p | success) / (1-p); the law of finite survivors, which is also the law produced by silent redraw-to-success | ratio |
| `rho_*_lo / rho_*_hi` | 99% seed-stratified paired bootstrap percentile interval, 10000 resamples | ratio |
| `conditional_log_r_quantile` | empirical quantile among survivors | log units |
| `conditional_log_r_quantile_error` | survivor quantile minus target quantile | log units |
| `quantile_resolved` | false when the sample cannot resolve the level; the estimate is then not interpreted | -- |
| `unconditional_failure_probability` | per-attempt failure probability, reported before any conditional statistic | probability |

## `loader_validation.csv`

| column | meaning | unit |
|---|---|---|
| `family` | multiplicity family: radial, direction, independence, frame, cap_law, negative_control, or diagnostic | -- |
| `test` | statistic name | -- |
| `statistic / pvalue` | raw test statistic and p-value, preserved for source data | -- |
| `corrected_decision` | Holm-corrected decision at familywise alpha 0.01; `detected`/`MISSED` for negative controls | -- |
| `n_total / n_finite` | returned records and those with three finite components | count |

## `performance.csv`

| column | meaning | unit |
|---|---|---|
| `scope` | `block` (one timed block), `summary` (median and IQR over blocks), `ratio` (paired blockwise LOG/SPLIT), `reproducibility` | -- |
| `seconds_per_attempt / seconds_per_return` | wall time per attempt and per returned sample | s |
| `attempts_per_return` | cap rejection cost | ratio |
| `engine_calls_per_return` | mt19937 outputs consumed per returned sample | count |
| `median_time_per_return_ratio` | paired blockwise LOG/SPLIT ratio; 95% bootstrap interval in ratio_ci_lo/hi | ratio |
| `repeat_identical` | true when the same build, method and seed reproduced identical counters | -- |

## `portability.csv`

| column | meaning | unit |
|---|---|---|
| `arch` | hardware architecture | -- |
| `available / completed` | whether the environment exists on this host and whether the job ran | -- |
| `reason_not_run` | exact reason and the command to run the job elsewhere; never blank | -- |

## `scalar_validation.csv`

| column | meaning | unit |
|---|---|---|
| `method` | LOG-ID or LOG-PUB | -- |
| `measured_acceptance` | accepted log-Gamma draws divided by proposals | probability |
| `published_acceptance_eq2` | Liu, Martin and Syring (2017) Eq. (2) as printed | probability |
| `published_acceptance_area_ratio` | Gamma(a+1)/(1+w), the ratio of areas implied by the same paper's target and envelope | probability |
| `z_ks_stat / z_ks_p / z_cvm_stat / z_cvm_p` | KS and Cramer-von Mises against Exp(1) for Z = -log I_W(a,3/2) | -- |
| `w_ks_* / w_cvm_*` | the same against Beta(a,3/2) for W, computed only where W is numerically resolved | -- |
| `log_r_q*_estimate / _target / _error / _lo / _hi / _resolved / _covers_target` | log R quantiles with exact order-statistic intervals, Bonferroni-corrected over the five levels | log units |

## Bulk files under `raw/`

Not tracked in git, regenerable with the commands in `raw/environment.json`, and covered by `raw/raw_checksums.sha256`. Each begins with a 96-byte header naming its schema version, record size, configuration and record count; the reader refuses a file whose declared record count disagrees with its length.

| file pattern | contents |
|---|---|
| `pilot/pilot_*.bin` | per-draw `log R` and `log W` from the E1 scalar pilot |
| `conditioning/cond_*.bin` | per-attempt intended state and per-method terminal category |
| `loader/loader_*.bin` | returned three-vectors, intended `log R`, status and attempt count |
| `*/audit_*.bin` | primitives of audited attempts, handed to the 100-digit oracle |
| `*/*.jsonl` | per-configuration counters, one JSON object per line |
