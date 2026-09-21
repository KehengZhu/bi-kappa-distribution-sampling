# Experiment 7 — output schema

Every column and field the probe and the oracle emit, per phase, with units and
missing-value codes. The analysis is written against this document; the probe is the only
producer, and `raw/manifest.csv` names every file it wrote.

Nothing here is a decision rule. Levels, margins, critical values and acceptance thresholds
live in `config/protocol.json` and are read by `exp7_families.py`; what follows is only what
the run measured.

---

## 0. Conventions

| convention | meaning |
|---|---|
| **units** | every logarithm is a natural logarithm; every time is in seconds; every count is a count of attempts unless the name says `returned` |
| **`log R`** | the radius on the log scale, `log R = (log X1 − log X2)/2` |
| **`a`** | the denominator shape `kappa − 1/2`, formed **in the working precision** and then widened — `shape_a` in every row |
| **`Z`** | `Z = −log I_W(a, 3/2)` with `W = X2/(X1+X2)`; standard exponential under the target law at every shape |
| **attempt** | one execution of the core mapping. In uncapped mode it is also one returned sample; in capped mode a return may consume many |
| **JSON non-finite** | `±inf` and NaN are emitted as the **strings** `"inf"`, `"-inf"`, `"nan"`, never as bare tokens. Any numeric field may carry one |
| **JSON `cap`** | `null` means uncapped (`+inf`); a number is the cap in thermal speeds |
| **missing** | a field that does not apply to a row is **absent**, never zero. `case` is `"-"` and `kappa` is `nan` in manifest rows that span configurations; `seed` is `-1` there |
| **binary** | little-endian, IEEE-754, natural struct padding as produced by the two builds (both verified to agree on every record size) |

Every JSONL row carries the run's identity:

| field | type | meaning |
|---|---|---|
| `phase` | string | `p1` … `p6` |
| `layer` | string | see each phase |
| `method` | string | `LEGACY`, `CANDIDATE`, `QF`, or `all` in a manifest row |
| `tag` | string | environment tag: `libcxx` or `libstdcxx` |
| `execution` | string | `native` or `translated` (Rosetta/emulation; excluded from the F7 decision) |
| `protocol_sha256` | string | SHA-256 of the `config/protocol.json` the probe was built against |
| `compiler`, `stdlib`, `arch` | string | e.g. `clang 21.0.0 …`, `libc++`, `arm64` |
| `precision` | string | `float` or `double` |
| `rng` | string | always `mt19937` |
| `kappa` | float | spectral index |
| `shape_a` | float | `kappa − 1/2` in the working precision, widened to double |
| `seed` | int | the replicate: one of 9001–9005 (production) or 9006–9010 (P6) |
| `stream_seed` | int | the engine this configuration actually ran on. Present on P1 and P5 rows and on P2's `native` rows; `seed + 10000 × (13 × precision_index + kappa_index)` by PROTOCOL.md §3, which gives every configuration of a replicate an independent stream. Absent elsewhere, where the engine seed is the replicate |
| `n_attempted` | int | attempts the row was asked for. In P4 it is the attempts actually made |

### Terminal categories

Exactly one per method per attempt. Columns are named `cat_<name>`.

| name | meaning | accounting home |
|---|---|---|
| `finite` | a finite vector, with an accurate radius | finite |
| `denominator_zero` | the denominator Gamma materialized as exactly zero | avoidable |
| `quotient_first_loss` | the intermediate `X1/X2` overflowed (QF diagnostic only) | avoidable |
| `legacy_form_loss` | the split form lost a draw whose intended vector was representable | avoidable |
| `log_primitive_failure` | the log-domain primitive produced a non-finite log | avoidable |
| `finite_but_wrong` | a finite vector was returned for a representable draw, but its radius lost more than half the significand of its type | avoidable |
| `honest_overflow` | the intended vector is not representable and a non-finite vector was returned | honest |
| `overflow_returned_finite` | the intended vector is not representable and a finite vector was returned anyway | honest |
| `cap_reject` | the intended sample lies outside the declared box (capped rows only) | neither |
| `cap_exhausted` | capped mode hit the internal attempt limit | neither |

`n_finite`, `n_avoidable`, `n_honest` are the three sums, and
`n_attempted = n_finite + n_avoidable + n_honest` in every uncapped row. `accounting_ok`
reports that identity as a boolean, evaluated by the probe.

The FINITE_BUT_WRONG threshold is `accuracy.max_relative_error` from the protocol:
`1.0536712127723509e-08` for `double`, `0.000244140625` for `float`. It is echoed in every
paired row as `max_rel_error_threshold`.

---

## 1. `raw/manifest.csv` — one row per output file, written at run time

Written by the probe itself, from inside the process that produced the file, **after** the
atomic rename and after re-opening the file and checking its header, record count and byte
length against what was intended. A smoke run writes `raw/smoke/manifest.csv` instead and
never touches this one.

| column | type | meaning / missing value |
|---|---|---|
| `manifest_schema` | string | `exp7-manifest-1` |
| `run_utc` | string | ISO-8601 UTC, seconds, of the process that wrote the file |
| `phase` | string | `p1` … `p6` |
| `layer` | string | `summary` (the phase's JSONL), `tail`, `paired`, `audit`, `class`, `capped`, `uncapped` |
| `method` | string | `LEGACY`, `CANDIDATE`, `both`, `paired`, `all` |
| `precision` | string | `float`, `double`, `mixed`, `all` |
| `case` | string | `C0`…`C6`, or `-` |
| `kappa` | float | `nan` for files spanning configurations |
| `seed` | int | `-1` for files spanning seeds |
| `file` | string | path relative to the experiment directory |
| `format` | string | `bin` or `jsonl` |
| `n_records` | int | records (bin) or lines (jsonl), as verified after the rename |
| `record_size` | int | `sizeof(record)` in bytes; `-1` for jsonl |
| `schema` | int | record schema version; `-1` for jsonl |
| `record_kind` | int | 2 conditioning, 3 loader, 4 audit, 5 tail; `-1` for jsonl |
| `bytes` | int | file length in bytes |
| `sha256` | string | SHA-256 of the file as written |
| `probe` | string | `argv[0]` of the writing process |
| `probe_sha256` | string | SHA-256 of that executable, read by the process itself |
| `cxxflags` | string | the flags the binary was **compiled with**, baked in at build time |
| `compiler`, `stdlib`, `arch` | string | build identity |
| `env_tag` | string | `--tag` |
| `execution` | string | `native` or `translated` |
| `protocol_sha256` | string | protocol the build was generated against |
| `smoke` | bool | `true` only under `raw/smoke/` |

---

## 2. Binary record layouts

Every `.bin` opens with this 96-byte header:

```c
struct FileHeader {
  char     magic[8];            // "EXP7REC\0"
  uint32_t schema;              // record schema version
  uint32_t record_kind;         // 2 cond, 3 loader, 4 audit, 5 tail
  uint32_t record_size;         // sizeof(record)
  uint32_t reserved;            // 0
  double   kappa;               // 0 for files spanning configurations
  double   theta_perp, theta_par;
  double   ub[3];
  double   cap;                 // +inf when uncapped
  uint32_t seed;                // 0 for files spanning seeds
  uint32_t precision_is_float;  // 1 float, 0 double
  int64_t  n_records;           // written on close
};
```

`CondRecord` (kind 2, schema 1, **32 bytes**), one per P3 attempt:

| field | type | units | meaning |
|---|---|---|---|
| `log_r_ref` | double | log units | intended `log R`, from the recorded primitives in double |
| `log_speed_ref` | double | log units | intended `log |V|` |
| `max_log_component_ref` | double | log units | intended `max_j log |V_j|` |
| `flags` | uint32 | bitfield | see below |
| `cat_legacy`, `cat_candidate`, `cat_qf` | uint8 | enum | terminal category index |
| `reserved` | uint8 | — | 0 |

`log W` is **not** stored: it is exactly `−logaddexp(0, 2·log_r_ref)`.

`LoaderRecord` (kind 3, schema 1, **40 bytes**), one per P4 *returned* sample:

| field | type | units | meaning |
|---|---|---|---|
| `v[3]` | double | thermal speeds × θ | returned components, global frame; may be `±inf`/NaN |
| `log_r_ref` | double | log units | the radius this attempt carried or materialized; NaN in `class` dumps |
| `status` | uint32 | enum | terminal category index |
| `attempts` | uint32 | count | attempts this return consumed (1 when uncapped) |

`AuditRecord` (kind 4, schema 1, **72 bytes**), the oracle's input. It carries only
primitives, so the oracle re-derives every classification independently:

| field | type | meaning |
|---|---|---|
| `x1` | double | `Ga(3/2,1)` variate, exactly as produced in the working precision |
| `y` | double | `Ga(a+1,1)` variate of the boosting identity |
| `u` | double | the open-interval uniform of the boosting identity |
| `cos_theta`, `phi` | double | the paired layer's two **declared** direction primitives — see below |
| `log_x2_working` | double | what the working-precision log path computed |
| `kappa` | double | carried per record: one stream spans many configurations |
| `flags` | uint32 | bitfield, with `kFlagAudited` set |
| `precision_is_float` | uint32 | 1 float, 0 double |
| `cat_legacy`, `cat_candidate`, `cat_qf` | uint8 | the probe's verdicts, for the oracle to contradict |
| `reserved`, `pad` | uint8, uint32 | 0 |

`cos_theta` and `phi` are the paired diagnostic layer's own declared direction primitive.
They are **not** a record of what either released sampler drew: release 1.0.0 draws its two
angles through `std::uniform_real_distribution`, and release 2.0.0 draws no angle at all —
its direction comes from a rejection method on the unit disc that calls no transcendental.
What the audit stream carries, and what the oracle adjudicates, is the radius formations on
the layer's declared draw. All five declared primitives are exact functions of the engine,
so the stream is identical across standard libraries by construction; see §4.

`TailRecord` (kind 5, schema 1, **8 bytes**), one per P1 draw whose `Z` exceeded the lowest
protocol threshold `z0 = −log(1e−2)`:

| field | type | meaning |
|---|---|---|
| `z` | double | the exceeding `Z`, finite by construction |

Draws whose `Z` is not finite are **not** written: they exceed every threshold, and writing
them would turn a 10⁴-record file into a 10⁶-record one. They are carried as
`tail[].observed_nonresolved` in the summary row.

### `flags` bitfield

| bit | name | meaning |
|---:|---|---|
| 0 | `x2_zero` | the product form `y·u^(1/a)` materialized as exactly zero |
| 1 | `x2_subnormal` | it landed in the subnormal range |
| 2 | `qf_nonfinite` | the QF radius is not finite |
| 3 | `legacy_nonfinite` | the LEGACY radius is not finite |
| 4 | `candidate_nonfinite` | the carried `log R` is not finite |
| 5 | `honest_overflow_ref` | the intended vector is not representable in the run type |
| 6 | `cap_accept` | reserved; not set by the phases in this run |
| 7 | `audited` | the attempt was written to the audit stream |
| 8 | `near_limit` | within `audit.margin_log_units` of a type limit |
| 9 | `rotation_recoverable` | the local vector overflowed while the rotated one is representable |
| 10 | `cap_exhausted` | reserved; not set by the phases in this run |
| 11 | `legacy_success` | LEGACY's terminal category is `finite` |
| 12 | `candidate_success` | CANDIDATE's terminal category is `finite` |

---

## 3. P1 — `raw/p1/p1_<tag>.jsonl`

Four rows per (`kappa`, `precision`, `seed`): a `class` row for each method, then a `native`
row for each. 13 × 2 × 5 × 4 = 520 rows per build.

### `layer = "class"` — the released class, run end to end

This is what a user receives, and what the portability digests are taken over.

| field | type | meaning |
|---|---|---|
| `n_finite` | int | returned vectors with three finite components |
| `nonfinite_output` | int | returned vectors with any non-finite component |
| `thrown` | int | draws that threw (0 uncapped) |
| `n_attempts_reported` | int | `n_attempts()`; **CANDIDATE only**, 0 for LEGACY, which has no counter |
| `n_nonfinite_reported` | int | `n_nonfinite()`; CANDIDATE only |
| `seconds` | float | wall time for the configuration |
| `digest_sha256` | string | SHA-256 over the IEEE bit patterns of every returned component, in order, then over `n_finite`, `n_nonfinite`, `n_thrown` and (CANDIDATE) the two reported counters |
| `digest_fnv1a64` | string | FNV-1a 64 over the identical byte sequence, as a cross-check on the hash itself |
| `digest_bytes` | int | bytes fed to both digests |

The digest is the F7 same-architecture statistic: for CANDIDATE it is predicted to be
**equal** across `tag` values on one architecture, and for LEGACY it is predicted to differ.
Both predictions are tested on these rows.

### `layer = "native"` — the verified replica, with the radius observable

Same stream as the `class` row of the same method, seed and configuration — both run on this
configuration's `stream_seed`, not on the replicate's `seed`; `make selftest` checks the two
bit for bit. Everything F1–F4 needs comes from here.

| field | type | meaning |
|---|---|---|
| `cat_*`, `n_finite`, `n_avoidable`, `n_honest`, `nonfinite_output` | int | terminal categories. **In P1 these are the natively observable split only**: with no paired reference, a draw whose denominator materialized as zero cannot be separated into avoidable and honest. P2's `paired` rows are the authoritative decomposition |
| `x2_zero`, `x2_subnormal` | int | LEGACY only; 0 for CANDIDATE, which never forms `X2` |
| `gamma_variates`, `uniform_variates` | int | top-level variates consumed. `gamma_variates` is `2 × n_attempted` for both methods. `uniform_variates` is `2 × n_attempted` for LEGACY, one per angle; for CANDIDATE it is **not** a fixed multiple of `n_attempted` — one uniform per attempt for the shape-boosting identity, plus two per try of the direction's disc rejection, which averages `4/π = 1.273` tries |
| `engine_calls` | int | `mt19937` invocations |
| `seconds` | float | wall time |
| `max_finite_log_r` | float | largest finite `log R` seen; `"-inf"` if none |
| `n_resolved` | int | draws whose `log R` is finite |
| `n_nonresolved` | int | draws lost to overflow. **They are not dropped**: they are sorted to `+inf`, so the order statistics below are the order statistics of the whole sample |
| `loss_fraction` | float | `n_nonresolved / n_attempted` |
| `conditional` | bool | `true` when `n_nonresolved > 0`; such a cell cannot by itself support the fidelity claim |
| `quantiles` | array | five objects, one per protocol level |
| `ad_statistic`, `ks_statistic`, `cvm_statistic` | float | exact Anderson–Darling `A²`, Kolmogorov–Smirnov `D` and Cramér–von Mises `W²` of `Z` against `Exp(1)`, on the `n_resolved` sample. `"nan"` when `n_resolved < 8` |
| `tail` | array | three objects, one per protocol `q0` |
| `z_grid_points`, `z_grid_step` | int, float | 501 and 0.05 |
| `z_ecdf_counts` | array[501] | cumulative count of resolved draws with `Z ≤ 0.05·i`, `i = 0…500` |
| `tail_records`, `tail_file` | int, string | the companion `TailRecord` file and its length |
| `digest_*` | — | as for the class row, over the replica's returned vectors and counters |
| `n_rel_err`, `mean_rel_err_radius`, `max_rel_err_radius` | int, float | 0 and 0 in P1: the native layer has no reference to measure accuracy against. Accuracy is a P2/P3 quantity |

`quantiles[i]`:

| field | type | meaning |
|---|---|---|
| `p` | float | 0.5, 0.9, 0.99, 0.999, 0.9999 |
| `lo_index`, `hi_index` | int | the 0-based order-statistic bracket frozen in `config/protocol.json` for `n = 10⁶` |
| `achieved_coverage` | float | the protocol's analytically computed coverage of that bracket |
| `lo_value`, `hi_value` | float | the sample `log R` at those indices — the F2 interval. May be `"inf"` |
| `point_value` | float | the sample `log R` at `ceil(p·n) − 1`, the empirical quantile F3 uses |

F2 informativeness is `hi_value − lo_value ≤ statistics.informative_width_log_units`; the
probe does not classify it.

`tail[i]`:

| field | type | meaning |
|---|---|---|
| `q0` | float | 1e−2, 1e−3, 1e−4 |
| `z0` | float | `−log q0` |
| `observed_resolved` | int | resolved draws with `Z > z0` |
| `observed_nonresolved` | int | unresolved draws, which exceed every threshold; the same number in all three objects |
| `expected` | float | `n_attempted · q0` under the null |

The exceedance count F4 tests is `observed_resolved` against `Binomial(n_resolved, q0)` for a
cell with no loss, and the analysis decides how to treat `observed_nonresolved` in a
conditional cell; the probe reports the two separately rather than choosing.

The `A²`, `D` and `W²` above are computed with the estimators the analysis module uses, so
they agree to the last bit: `u_i = −expm1(−z_(i))` clipped to `[1e−300, 1 − 1e−16]`,
`A² = −n − Σ(2i−1)(log u_i + log1p(−u_(n+1−i)))/n`, `W² = 1/(12n) + Σ(u_i − (2i−1)/2n)²`,
`D = max_i max(i/n − u_i, u_i − (i−1)/n)`.

### `raw/p1/p1_tail_<tag>_<method>_k<kappa>_s<seed>_<precision>.bin`

`TailRecord` stream, kind 5. Records are written in **descending** `Z`. Their count is
`tail_records` in the summary row and `n_records` in the header.

---

## 4. P2 — `raw/p2/p2_<tag>.jsonl` and `raw/p2/audit_p2_<tag>.bin`

Five rows per (`kappa`, `precision`, `seed`): three `paired` (LEGACY, CANDIDATE, QF) and two
`native`. 13 × 2 × 5 × 5 = 650 rows per build.

### `layer = "paired"` — the mechanism decomposition

All three formations receive the *same* declared primitives `(x1, y, u, cosθ, φ)`, so "was
this draw recoverable?" is decided per draw. What the layer compares is the three **radius**
formations on one mathematical draw; the direction is a shared nuisance parameter and every
column is handed the identical one.

Each of those five primitives is an exact function of the engine — the Gamma, normal and
uniform generators live in the header rather than in `<random>` — so the **audit stream** of
the two builds has the same SHA-256 by construction rather than by result. That guarantee
covers the recorded primitives. It does not extend to the layer's derived counters: the
declared direction reaches the vector through `sin(φ)` and `cos(φ)` of one argument, which
clang fuses into Darwin's `__sincos_stret` while GCC does not, so a build can differ in the
last place of a direction component. A terminal category can only turn on that where a
component sits within an ulp of a type limit, which falls in the fully audited within-margin
stratum. The portability claim of record is P5's, taken over the released class.

| field | type | meaning |
|---|---|---|
| `diagnostic_only` | bool | `true` for the `QF` row; no claim rests on it |
| `cat_*`, `n_finite`, `n_avoidable`, `n_honest` | int | the authoritative decomposition |
| `accounting_ok` | bool | `n_attempted = n_finite + n_avoidable + n_honest` |
| `x2_zero`, `x2_subnormal` | int | of the shared primitives; identical across the three rows |
| `ref_nonrepresentable` | int | attempts whose intended vector has no representation in the type |
| `rotation_recoverable` | int | local vector overflowed, rotated one representable |
| `near_limit` | int | within the margin of a type limit |
| `max_log_r_ref` | float | largest intended `log R` |
| `max_rel_error_threshold` | float | the type's FINITE_BUT_WRONG threshold |
| `n_rel_err` | int | attempts on which the relative radius error is defined |
| `mean_rel_err_radius`, `max_rel_err_radius` | float | that error, **per precision** — never pool the two |
| `max_finite_log_component` | float | largest intended `max_j log|V_j|` among successes |
| `audited` | int | records this configuration contributed to the audit stream |
| `audit_<stratum>_total` | int | attempts that fell in the stratum |
| `audit_<stratum>_audited` | int | of those, how many were written |
| `audit_<stratum>_rate` | float | the protocol's declared rate for the stratum |
| `audit_margin_log_units` | float | `audit.margin_log_units`, 20.0 |
| `audit_margin_assertions` | int | attempts declared unambiguous, each of which had its margin evaluated |
| `audit_margin_assertion_failures` | int | of those, how many failed the inequality. **Must be zero** |
| `audit_min_unambiguous_margin` | float | the smallest margin any of them cleared; `"inf"` when there were none |

The seven strata, in the priority order an attempt is assigned to them — the first match
wins, so the counts partition the sample:

1. `within_margin_of_a_type_limit` (rate 1)
2. `method_disagreement` (1)
3. `finite_but_wrong_candidate` (1)
4. `avoidable_loss_candidate` (1)
5. `subnormal_or_zero_denominator_representable_target` (1)
6. `unambiguous_failure_beyond_margin` (1e−3)
7. `uniform_sample` (1e−4)

Because assignment is by priority, stratum 7 is an unbiased sample of *the attempts in no
other stratum*, not of all attempts. Both counts are reported for each, so the coverage is
exact rather than implied.

### `layer = "native"`

As P1's native rows, minus the scalar summary and the tail file, plus `accounting_ok`. Same
`stream_seed` and stream as P1's native row for the same configuration, so the two must agree
exactly; that cross-phase equality is a free consistency check. P2's `paired` rows keep the
replicate seed, as PROTOCOL.md §3 states — nothing compares them against P1.

### `raw/p2/audit_p2_<tag>.bin`

One `AuditRecord` stream for the whole phase, spanning every configuration. `kappa` and
`precision_is_float` are per record; the file header's `kappa` and `seed` are 0.

---

## 5. P3 — `raw/p3/p3_<tag>.jsonl`, `p3_cond_*.bin`, `audit_p3_<tag>.bin`

Three `paired` rows per (`kappa`, `precision`, `seed`) over the frozen cells — double
`kappa ∈ {0.501, 0.505, 0.51, 0.75}`, float `kappa ∈ {0.55, 0.60, 0.75}` — so 7 × 5 × 3 = 105
rows per build.

Fields are P2's `paired` fields, plus:

| field | type | meaning |
|---|---|---|
| `raw_file` | string | the `CondRecord` file for this configuration |

`raw/p3/p3_cond_<tag>_k<kappa>_s<seed>_<precision>.bin` holds one `CondRecord` per attempt,
in attempt order. This is the only phase that keeps a per-attempt record of a 10⁶-attempt
run, because conditioning is a statement about the joint distribution of the intended state
and success, which no summary preserves.

---

## 6. P4 — `raw/p4/p4_<tag>.jsonl` and the per-case record files

Per case and seed: one `uncapped`/`capped` row per method (the verified replica, with the
bookkeeping), plus — for the five uncapped cases — one `class` row per method dumping the
released class itself. 7 × 2 × 5 + 5 × 2 × 5 = 120 rows per build.

`layer` is `uncapped`, `capped` or `class`.

| field | type | meaning |
|---|---|---|
| `case`, `role` | string | `C0`…`C6` and its role from PROTOCOL.md §4.1 |
| `theta_ratio` | float | `θ∥/θ⊥`, with `θ⊥ = 1` |
| `ub` | array[3] | the field direction |
| `cap` | float or null | `null` when uncapped |
| `n_returned` | int | returned samples |
| `attempts` | int | attempts made. **Uncapped cases stop at exactly 10⁵ attempts; capped cases stop at 10⁵ returned samples** |
| `nonfinite_attempt` | int | attempts whose vector had a non-finite component |
| `nonfinite_returned` | int | of the returned samples, how many |
| `cap_reject` | int | attempts the box rejected |
| `cap_reject_unrepresentable` | int | of those, how many were also unrepresentable — the count that tells "the box rejected a representable sample" from "the box rejected one that could not have been returned at all" |
| `cap_exhausted` | int | requests that hit the internal attempt limit |
| `max_attempt_run` | int | longest rejection run behind a single return |
| `x2_zero` | int | LEGACY only |
| `engine_calls` | int | `mt19937` invocations |
| `seconds` | float | wall time |
| `raw_file` | string | the `LoaderRecord` file |
| `digest_*` | — | over the returned vectors and the four counters above |

`class` rows carry `n_finite`, `nonfinite_output`, `thrown`, `n_attempts_reported`,
`n_nonfinite_reported`, `seconds`, `digest_*` and `raw_file`, exactly as P1's class rows do.

---

## 7. P5 — `raw/p5/p5_<tag>.jsonl`

CANDIDATE only, the full P1 ladder, one `class` row per (`kappa`, `precision`, `seed`):
13 × 2 × 5 = 130 rows per build. Fields are identical to P1's `class` rows.

The same-architecture F7 test is equality of `digest_sha256` across `tag` values at equal
(`kappa`, `precision`, `seed`). The LEGACY half of the contrast — where equality is predicted
**not** to hold — comes from P1's class rows, which cover the same ladder and seeds.

There are two tags and no third. An earlier revision carried a diagnostic build,
`libcxx-nosincos`, to locate why this equality failed; the cause — clang fusing the
`sin(φ)`/`cos(φ)` pair of the released sampler into `__sincos_stret` — was then removed from
the released header, which now draws its direction by disc rejection and calls no
transcendental. The diagnostic build and its tag no longer exist, and no row carries them.

---

## 8. P6 — `raw/p6/p6_<tag>.jsonl`

Per case: 10 blocks × 2 methods timed rows, then 2 reproducibility rows. 5 × 22 = 110 rows
per build. `layer` is `uncapped`, `capped` or `reproducibility`.

Timed rows:

| field | type | meaning |
|---|---|---|
| `case` | string | `B0`…`B4` |
| `cap` | float or null | |
| `block` | int | 0–9 |
| `slot` | int | 0 or 1: the order this method ran **within** the block |
| `candidate_first` | bool | the realized randomization for that block, recorded so it can be checked rather than asserted |
| `n_returned` | int | samples drawn in the block |
| `n_finite`, `n_nonfinite` | int | of those |
| `attempts_reported` | int | `n_attempts()`; CANDIDATE only, 0 for LEGACY |
| `seconds` | float | wall time for the block |
| `seconds_per_returned` | float | the G5 statistic's per-block input |
| `seed` | int | 9006–9010, indexed by `block mod 5` from the declared vector |

Block size is calibrated once per case and method by an **untimed warm-up**, which is never
reported, so that every timed block clears both floors (≥ 2 s and ≥ 10⁶ attempts).

Reproducibility rows:

| field | type | meaning |
|---|---|---|
| `repeat_identical` | bool | two same-seed runs produced identical counters |
| `finite_1`, `finite_2`, `nonfinite_1`, `nonfinite_2`, `attempts_1`, `attempts_2` | int | the two runs, so the boolean can be audited |

---

## 9. Oracle — `raw/oracle_audit.jsonl` and `raw/oracle_disagreements.jsonl`

`exp7_oracle.exe` reads only the primitives of each audited attempt and re-derives every
classification at 100 decimal digits with its own code.

`raw/oracle_audit.jsonl`: one row per (audit file, precision) that had at least one record,
plus one `n_records: 0` row for a file with none.

| field | type | meaning |
|---|---|---|
| `tool`, `oracle` | string | `exp7_oracle`, `boost::multiprecision::cpp_dec_float_100` |
| `file` | string | the audit file |
| `file_sha256` | string | SHA-256 of that file, so the adjudication is tied to the bytes it adjudicated |
| `protocol_sha256` | string | |
| `precision` | string | `float` or `double`; **statistics are never pooled across the two** |
| `n_records` | int | records of this precision |
| `n_records_total` | int | records in the file, both precisions |
| `n_classified` | int | records whose categories were compared |
| `disagreements` | int | probe-versus-oracle category mismatches in the whole file |
| `conversion_failures` | int | records whose primitives did not convert exactly to 100 digits |
| `near_limit_records` | int | records within the margin of a type limit |
| `worst_overflow_margin_log_units` | float | largest `\|max_log_component − log_ovf\|` seen |
| `n_subnormal_denominator` | int | records whose denominator went subnormal |
| `max_abs_log_x2_error`, `mean_abs_log_x2_error` | float | working-precision `log X2` against the exact one, in log units |
| `n_radius_compared`, `n_radius_compared_subnormal` | int | records where both radii resolved |
| `max_rel_error_threshold` | float | the FINITE_BUT_WRONG threshold for this precision |
| `returned_finite_legacy`, `returned_finite_candidate` | int | records each method returned a finite vector for |
| `finite_but_wrong_legacy`, `finite_but_wrong_candidate` | int | of those, how many were inaccurate. The candidate's count is the one G1 cares about |
| `max_rel_err_legacy`, `mean_rel_err_legacy` | float | relative radius error against the exact value |
| `max_rel_err_candidate`, `mean_rel_err_candidate` | float | the same for the candidate |
| `max_rel_err_legacy_subnormal_denominator` | float | worst case restricted to subnormal denominators — where the split form is worst on exactly the draws it survives |
| `max_rel_err_candidate_subnormal_denominator` | float | the same for the candidate |

`raw/oracle_disagreements.jsonl` always carries a `header` and a `footer` record per audit
file, whether or not anything disagreed, so an empty adjudication is distinguishable from an
oracle that never ran.

| `kind` | fields |
|---|---|
| `header` | `tool`, `file`, `file_sha256`, `n_records_declared`, `protocol_sha256`, `max_rel_error_double`, `max_rel_error_float` |
| `disagreement` | `file`, `file_sha256`, `kappa`, `precision`, `x1`, `y`, `u`, `cos_theta`, `phi`, `probe` (3 category indices), `oracle` (3 indices), `margin_log_units` |
| `footer` | `file`, `file_sha256`, `n_records`, `disagreements`, `conversion_failures` |

`exp7_oracle.exe` exits non-zero on any disagreement or conversion failure, and `make oracle`
propagates that exit status.

---

## 10. `raw/environment.json`

Written by `preflight.py`, not by the probe. `./exp7_probe_<tag>.exe env` prints the
floating-point environment as one JSON object:

| field | meaning |
|---|---|
| `compiler`, `stdlib`, `arch`, `rng` | build identity |
| `bi_kappa_version` | `BI_KAPPA_VERSION_STRING` of the header compiled in |
| `protocol_version`, `protocol_sha256` | the protocol the build was generated against |
| `cxxflags` | the flags baked into the binary |
| `rounding_mode` | `FE_TONEAREST` etc. |
| `flt_eval_method` | `FLT_EVAL_METHOD` |
| `double_subnormals_flushed`, `float_subnormals_flushed` | whether FTZ/DAZ or the FPCR FZ bit is on — it would turn the small-shape denominator into a zero earlier than IEEE arithmetic would, which is the mechanism under study |
| `sizeof_long_double` | bytes |
| `fp_contract` | whether an FMA is available, and that the build disables contraction |
| `log_overflow_threshold_double/_float` | `log((2 − 2^−p)·2^emax)`, the exact threshold at which a real rounds to infinity |
| `log_max_double/_float` | `log(max())`, the predicate the released sampler itself applies |
| `max_rel_error_double/_float` | the FINITE_BUT_WRONG thresholds |
| `mt19937_min`, `mt19937_max` | engine range |
