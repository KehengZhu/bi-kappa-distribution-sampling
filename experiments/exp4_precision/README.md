# Experiment 4 — finite-precision / low-κ audit of the released sampler

Answers **R1.3** (low-κ range), **R2.B** (κ domain and convention) and the **C5** decision,
and supplies the supported-range statement the Software section needs.

The question is *not* "can we invent a new low-κ algorithm". It is:

> What is the actual supported numerical range of the released implementation as a function
> of κ, precision and standard-library implementation, and can any **avoidable** failure be
> removed without changing the target law?

The object under test is the released C++ header `cpp/bi_kappa_distribution.H`. NumPy is not
used as a reference implementation anywhere in this experiment.

## Reproducing

```bash
make run                                              # build both toolchains + all raw output + checksums
uv run --project ../../python python exp4_analyze.py  # -> results/
```

`make run` regenerates everything from scratch. Nothing in `results/` depends on state that is
not in the committed source or the command line above. `make verify` re-checks a regenerated
`raw/` against `raw/checksums.sha256` — needed because the `logr_*.bin` bulk draws are
gitignored, so the checksums are the only way to confirm a regenerated set matches the one the
committed results were computed from. The environment block in `results/exp4_results.json`
pins the released header under test by content via `sampler_header_sha256`.

## Environment

Built twice from one source, because the answer could in principle depend on the standard
library's Gamma generator:

| Tag | Compiler | Standard library |
|---|---|---|
| `libcxx` | Apple clang 21.0.0 | libc++ |
| `libstdcxx` | Homebrew GCC 15.2.0 | libstdc++ |

arm64 (Apple M4 Max), `-O2 -std=c++11`, RNG `std::mt19937`, seeds 4001–4005,
10⁶ draws per (κ, precision, stdlib, seed) — 5×10⁶ per configuration, 1.8×10⁸ draws total for
the uncapped sweeps.

κ ladder: 0.5001, 0.501, 0.505, 0.51, 0.55, 0.60, 0.75, 1.0, 1.5. Precisions: `float` and
`double`, both as **real** pipelines — `bi_kappa_distribution<float>` instantiates
`std::gamma_distribution<float>` and does every operation in single precision. Nothing is
drawn in double and cast afterwards.

## The four modes

| Mode | What it does |
|---|---|
| `released` | Draws from `bi_kappa_distribution<T>` in `no_cap()` mode and counts non-finite outputs. Ground truth for "what the shipped code does". The released header is not instrumented. |
| `variates` | The mechanism decomposition. Drives a reference small-shape Gamma generator so that **identical variates** feed three radius formations, making the comparison exact draw-by-draw rather than distributional. |
| `capped` | Q5 — does the component-wise box cap prevent the invalid computation or merely resample it away? |
| `logdump` | Dumps `log R` from the log-domain path for distributional validation. |

### Why the mechanism decomposition is exact

For shape `a < 1` the standard boost identity is `X = Y·U^(1/a)` with `Y ~ Ga(a+1,1)`,
`U ~ U(0,1)` — which is what Marsaglia–Tsang style generators, and hence both libc++ and
libstdc++, use to reach `a < 1`. Shape `a+1 ∈ (1,2)` is numerically benign, so **all** of the
small-shape difficulty is concentrated in the single factor `U^(1/a)`:

```
X     = Y · pow(U, 1/a)        underflows to exactly 0 for small a
log X = log(Y) + log(U)/a      exact, and cannot underflow
```

Computing both from one `(Y, U)` pair means the same draw can be followed through every
formation at once, so "was this draw recoverable?" is decided per draw and not by comparing
histograms. `log R = (log x1 − log x2)/2` never overflows or underflows and is therefore the
arbiter.

The reference generator's failure statistics agree with the stdlib's own
`std::gamma_distribution` (`released` vs `variates` columns) to within seed noise, which is
what licenses using it as a stand-in.

## The three loss categories

Mutually exclusive, and they answer three different questions:

| Category | Meaning | Fixable? |
|---|---|---|
| **honest overflow** | the mathematical radius itself exceeds the largest representable number at this precision | **No.** No reformulation can eliminate this. |
| **spurious, ratio** | radius representable, but `sqrt(x1/x2)` lost it because the *quotient* overflowed before the square root | Yes — by `sqrt(x1)/sqrt(x2)`, already landed |
| **spurious, split** | radius representable, but `sqrt(x1)/sqrt(x2)` still lost it because the denominator Gamma variate underflowed to exactly 0 | Only by a log-domain construction — **not landed** |

## Headline results

Full tables in `results/exp4_table.md`; machine-readable in `results/exp4_results.json`.

**libc++ and libstdc++ agree to within seed noise at every configuration.** There is no
standard-library dependence to report.

Supported range of the released implementation as it currently stands:

| Precision | Zero observed failures | Degraded | Unusable |
|---|---|---|---|
| `double` | κ ≥ 0.55 | κ = 0.51 (5.7×10⁻⁴) | κ ≤ 0.505 |
| `float` | κ ≥ 0.75 | 0.55 ≤ κ ≤ 0.60 (5.7×10⁻³ … 2.9×10⁻⁵) | κ ≤ 0.51 |

This is consistent with Experiment 1, which independently saw non-finite draws only at
κ = 0.51 in double, at 5.8×10⁻⁴.

### Q2 — the `sqrt(x1)/sqrt(x2)` fix, on identical variates

| Precision | κ | loss, `sqrt(x1/x2)` | loss, `sqrt(x1)/sqrt(x2)` | recovered |
|---|---|---|---|---|
| float | 0.55 | 1.22×10⁻² | 5.7×10⁻³ | 6.5×10⁻³ |
| float | 0.60 | 1.5×10⁻⁴ | 2.9×10⁻⁵ | ≈5× |
| double | 0.51 | 8.2×10⁻⁴ | 5.7×10⁻⁴ | 2.5×10⁻⁴ |

Algebraically identical, and verified to agree to **2 ulp** (max relative difference
4.34×10⁻¹⁶ = 1.95·ε) at κ = 2 where neither formation can overflow. Regression tests
`test_radius_formation` A1–A4 in `cpp/test_suite.H`.

### Q3 — the log-domain construction

A log-domain small-shape Gamma removes the remaining *spurious* loss entirely, leaving only
honest overflow. What it would buy, in finite fraction:

| Precision | κ | current | log-domain bound |
|---|---|---|---|
| float | 0.55 | 0.99430 | 0.99986 |
| double | 0.505 | 0.97589 | 0.99918 |
| double | 0.501 | 0.52488 | 0.75810 |

Validated distributionally at κ = 0.55, 0.75, 2, 5 (3 seeds × 2×10⁵ draws each): 12/12 pass
KS **and** Cramér–von Mises against `W = 1/(1+T) ~ Beta(κ−1/2, 3/2)`, plus log-R quantile
probes at p = 0.5, 0.9, 0.99, 0.999. α = 0.01, fixed before the runs.

**This is measured but NOT landed in the released header.** Adopting it would replace
`std::gamma_distribution` with a hand-rolled generator, changing the RNG stream and breaking
seed-for-seed reproducibility against every previous version. That is a decision for the
authors, not a bug fix.

### Q5 — the cap masks the failure, it does not solve it

`non-finite returned` is **0 in every single configuration**, including those where 93% of
internal attempts produce a non-finite draw. A non-finite value can never satisfy the box
predicate, so the rejection loop silently redraws it.

| Precision | κ | λ | acceptance | non-finite per attempt | non-finite returned |
|---|---|---|---|---|---|
| double | 0.5001 | 5 | 3.7×10⁻⁴ | 0.928 | 0 |
| float | 0.51 | 5 | 3.6×10⁻² | 0.356 | 0 |
| double | 0.55 | 5 | 0.165 | 0 | 0 |

A user running capped mode at κ = 0.51 in single precision would see 36% of internal attempts
produce non-finite values and never learn of it. This is failure **hidden by truncation**, and
it must not be reported as failure solved. It also compounds the acceptance-rate cost: at
κ = 0.5001, λ = 5 the loop needs ≈2700 attempts per accepted sample.

## Numerical conventions (not optional)

* Never form `T = R²` at low κ — it overflows for `R > 1.3×10¹⁵⁴` even where `R` is
  representable, discarding exactly the heavy-tail draws the diagnostic exists to test.
* Never use `Y = T/(1+T)`. For small `κ−1/2` its mass piles up against 1.0 where no relative
  resolution remains, and a KS test on it measures rounding rather than the sampler.
* Use `W = 1/(1+T) ~ Beta(κ−1/2, 3/2)`, computed here as `expit(−2 log R)` directly from the
  log-domain radius, so `T` is never materialized at all.
* Use `hypot`-style norms, never `np.linalg.norm` or a naive sum of squares.

## Artifacts

| Path | Canonical? |
|---|---|
| `exp4_probe.cpp`, `exp4_analyze.py`, `GNUmakefile`, `README.md` | **canonical** (committed source) |
| `results/exp4_results.json`, `results/exp4_table.md` | **canonical** (compact summaries, committed) |
| `raw/*.jsonl` | **canonical**, tracked — per-configuration counters; small, and the direct input to every reported number |
| `raw/checksums.sha256` | **canonical**, tracked — makes a regenerated `raw/` verifiable |
| `raw/logr_*.bin` | regenerable (`make run`) — 18 MB of bulk log-R draws, gitignored |
| `*.exe` | build output, gitignored |
