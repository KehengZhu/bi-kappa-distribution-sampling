# Experiment 3 — per-sample timing against two other bi-Kappa samplers

This experiment times the uncapped sampler of this package against two published methods for
drawing isotropic Kappa velocities, after first checking that all three sample the intended
distribution. It is not used in a figure or table of the paper.

## Result

Both published methods were faster than this package's sampler on the test system.

- Zenitani's (2025) Pareto-envelope rejection sampler took 0.38 to 0.56 times as long per
  sample over κ ∈ [1.5, 50], with acceptance rates between 0.73 and 0.81.
- The normal-triple scale mixture of Abdul & Mace (2015) was faster by a factor of 1.47 to 1.89
  at every κ tested.

This package's sampler took 88 to 130 ns per sample (7.7 to 11.4 million samples per second),
and the cost is not monotone in κ. The fastest point is κ = 1.5, where the shape `κ − 1/2` of
the second Gamma variate equals 1. Standard-library Gamma generators commonly switch algorithm
at shape 1, which is a plausible cause, but the library was not instrumented to confirm it.
Anisotropy (θ∥ = 2θ⊥) and an oblique field direction each changed the time by at most 2.5 ns
per sample. A cap of λ = 20 added at most 3 ns per sample for κ ≥ 1.5, 5% at κ = 1 and 36% at
κ = 0.75, where it rejects 22% of attempts (Experiment 2).

## Methods compared

| tag | method | notes |
|---|---|---|
| `gamma_ratio_spherical` | this package, `cpp/bi_kappa_distribution.H` with `no_cap()` | `R = √X₁/√X₂` with `X₁ ~ Γ(3/2)`, `X₂ ~ Γ(κ − 1/2)`, and a uniform direction; equivalent to Zenitani & Nakano (2022), Alg. 1-1 |
| `scale_mixture_normals` | Abdul & Mace (2015), *Phys. Plasmas* **22**, 102107, Eq. (22) with Eqs. (19)–(20) | `v_i = θ√κ Z_i/√(χ²_ν)`, `ν = 2κ − 1`, three standard normals `Z_i` sharing one χ² variate |
| `pareto_rejection` | Zenitani (2025), *Res. Notes AAS* **9**, 299, Sec. 2 | rejection under a Pareto envelope with index `n = κ/2`, using uniform variates only |

The scale mixture uses the same construction as this package. Since `|Z|² ~ χ²₃ = 2 Γ(3/2)` and
`χ²_ν = 2 Γ(κ − 1/2)`, its radius has the same Gamma-ratio form. It differs only in drawing the
direction from three normal variates instead of two uniform variates. Abdul & Mace do not say how
the χ² variate with non-integer `ν` is generated. Here it is drawn as `2 Γ(ν/2)` with
`std::gamma_distribution`.

The envelope index `n = κ/2` recommended by Zenitani requires `0 < n < κ − 1/2`, that is κ > 1.
At κ = 0.75 and κ = 1 the rejection method is therefore recorded as not applicable. The measured
acceptance rates agree with Zenitani's published values to three digits (0.8060, 0.7856 and
0.7505 at κ = 1.5, 2 and 5, against 0.806, 0.785 and 0.750).

## Design

| | |
|---|---|
| validation | κ ∈ {0.75, 1, 1.5, 2, 5, 10}, seeds 3001–3003, 2 × 10⁵ draws per method and seed |
| timing | κ ∈ {0.75, 1, 1.5, 2, 5, 10, 50}, seed 3101, 10 timed batches of 10⁶ samples after one untimed warm-up batch |
| random-number generator | `std::mt19937` for every method, seeded identically |
| thermal speeds | θ⊥ = θ∥ = 1 for the comparison, because the two published methods are isotropic |

**Validation.** Before timing, each method is tested against the target distribution at every
applicable κ and seed. The radius is tested with Kolmogorov–Smirnov and Cramér–von Mises tests
of `W = 1/(1 + T) ~ Beta(κ − 1/2, 3/2)`, where `T = |v|²/(κθ²)`. The direction is tested with a
KS test of `cos θ` against U(−1, 1). The pass/fail decision applies a Holm–Bonferroni
correction at level α = 0.01 to the KS p-values of each kind (radius, direction) across all
methods, κ and seeds; the Cramér–von Mises p-values are recorded alongside. A method that fails
has its timings marked `usable: false`. All three methods pass. Without the correction one
test falls below 0.01 (direction, `gamma_ratio_spherical`, κ = 5, seed 3003, p = 0.0044). The
JSON output records the uncorrected results as well.

**Timing.** Every timing loop adds the returned components to a checksum, so the compiler cannot
remove the work being timed. The same checksum is applied to every method. Each configuration
reports the median, range and interquartile range of the ten batch times. This package's sampler
is additionally timed in three variants of its own: `aniso` (θ∥ = 2), `rotated` (θ∥ = 2 and field
direction `(0.3, −0.5, 0.8)`) and `capped20` (λ = 20).

## Rerunning

Run from this directory.

```bash
make run                                               # build, validate, time; write raw/
uv run --project ../../python python exp3_analyze.py   # read raw/, write results/
make verify                                            # check raw/ against the committed raw/checksums.sha256
```

`make run` rebuilds `exp3_bench.exe` against the current `cpp/bi_kappa_distribution.H`,
overwrites the committed `raw/validate.jsonl` and `raw/timing.jsonl` (`git checkout raw/`
restores them), writes the validation dumps `raw/val_*.bin`, and writes the checksums of the new
files to `raw/checksums_rerun.txt` (not committed). It never overwrites the committed
`raw/checksums.sha256`. `make verify` fails if any file listed there is missing or differs, so it
passes only on the original files. A rerun does not pass it: the timings differ from run to run,
and the validation dumps of the package's sampler differ because the current header draws
different variates from the version measured (see below). With version 3.0.0 only
`raw/validate.jsonl`, which records the rejection method alone, reproduces. The analysis reads
the validation dumps, which are not committed, so it needs a completed `make run`. `make clean`
removes the executable and `make distclean` also removes the validation dumps; the committed
files in `raw/` stay.

`exp3_bench.exe` can also be called directly:

```
exp3_bench.exe validate <out.bin> <kappa> <theta> <seed> <n>
exp3_bench.exe time <method> <variant> <kappa> <theta> <seed> <n> <repeats>
    method  = gamma_ratio_spherical | scale_mixture_normals | pareto_rejection
    variant = iso | aniso | rotated | capped20
```

**Runtime and disk use.** The timing phase draws about 4 × 10⁸ samples, which takes on the order
of a minute at the measured rates. The validation dumps take about 220 MB.

## Output files

| file | contents |
|---|---|
| `raw/validate.jsonl` | one line per (κ, seed): acceptance counts of the rejection method |
| `raw/timing.jsonl` | one line per (method, variant, κ): the ten batch times and the build identity |
| `raw/val_k<κ>_s<seed>.bin` | per draw and method: method id, `log\|v\|`, `cos θ` (not committed) |
| `raw/checksums.sha256` | SHA-256 of the two `.jsonl` files and every `val_*.bin` |
| `results/exp3_results.json` | validation results per method, κ and seed; timing summaries; build environment |
| `results/exp3_table.md` | the same results as Markdown tables |

## Header version and scope

The timings and validation dumps in `raw/` were recorded at commit 0139426 (2026-08-17), a
development version between release 1.0.0 and version 2.2.1, and measured
`cpp/bi_kappa_distribution.H` as committed there (SHA-256 beginning `6b138af5`). `raw/` does
not record the header's hash itself. That version computes the radius directly as `√X₁/√X₂`,
draws both Gamma variates with `std::gamma_distribution`, and draws the direction from `cos θ`
and `φ`. Since version 2.2.1 the header computes the radius on a logarithmic scale, generates
the Gamma variates itself, and draws the direction by Marsaglia's (1972) rejection method.
This benchmark has not been rerun on those versions, so its timings do not describe them.
Experiment 4 times the current radius calculation against the direct one.

`results/exp3_results.json` was regenerated later from the same `raw/` files. Its
`environment` block describes the machine and working tree at that time, not at the time of
measurement: `sampler_header_sha256` (beginning `8a3cfe3c`, version 2.2.1) is the header that
was on disk when the analysis ran, `git_commit` is b9a4609 with `git_dirty` true, and `cxx` is
the compiler then installed. The compiler that built the timed binary is recorded in each line
of `raw/timing.jsonl` (clang 21.0.0, clang-2100.1.1.101). Since the measurement
`exp3_bench.cpp` has changed only in its comments.

The measurements come from one machine and one toolchain: Apple clang 21 with libc++ on arm64
(Apple silicon), macOS 26.6.1, `-O2`. The measured version draws its Gamma variates from the
standard library, so its cost, and the ordering of the methods, may differ with another standard
library or processor.
