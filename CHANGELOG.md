# Changelog

This project follows [Semantic Versioning](https://semver.org/). For a sampler, "breaking"
includes changing which random numbers a given seed produces, even when the API is untouched
and the sampled law is unchanged — a stored stream is part of what a user depends on.

## 2.0.0 — unreleased

### Changed — the bi-Kappa radius is built in the logarithmic domain

`bi_kappa_distribution` draws its radius as `R = sqrt(X1)/sqrt(X2)` with
`X1 ~ Gamma(3/2, 1)` and `X2 ~ Gamma(kappa - 1/2, 1)`. As `kappa` approaches `1/2` the
denominator's shape approaches zero, `X2` underflows to exactly zero, and the draw is lost —
a loss no rearrangement of the quotient can survive, because the information is gone before
the division happens.

`X2` is therefore no longer formed. The sampler carries `log X2`, obtained from the
Ahrens–Dieter shape-boosting identity, propagates `log R = (log X1 - log X2)/2`, builds the
order-unity vector `g` with `V = R g` before exponentiating anything, and decides whether each
returned component is representable from `log|V_j| = log R + log|g_j|`. The velocity cap is
tested in the same domain, so a sample far outside the box is rejected rather than overflowed
first and rejected afterwards.

What remains lost is velocity that has no representation in the floating-point type, which is
a property of the type and not of the sampler. Measured over 2×10⁶ draws, isotropic and
unrotated:

| precision | kappa | 1.0.0 | 2.0.0 | analytic representability floor |
|---|---:|---:|---:|---:|
| double | 0.5001 | 9.285e-01 | 8.673e-01 | 8.676e-01 |
| double | 0.505 | 2.396e-02 | 8.350e-04 | 8.250e-04 |
| double | 0.51 | 5.690e-04 | 1.000e-06 | 6.807e-07 |
| float | 0.55 | 5.620e-03 | 1.420e-04 | 1.376e-04 |
| float | 0.75 | 0 | 0 | 5.3e-20 |

The floor is derived in
`experiments/exp7_confirmatory/config/honest_floor.md`; every 2.0.0 value agrees with it
within 1.2 standard errors.

### Changed — the random stream, and its portability

**A given seed now produces a different sequence than it did under 1.x.** This is the
intended consequence of drawing different variates, not a side effect, and no seed-for-seed
continuity is offered. Pin 1.0.0 if a stored stream must be reproduced.

In exchange, the stream no longer depends on which standard library the code is built
against. `std::gamma_distribution` uses Ahrens–Dieter GS in libc++ and Marsaglia–Tsang in
libstdc++; `std::normal_distribution` returns the two variates of its generated pair in
opposite orders; and `std::uniform_real_distribution<float>` returns exactly `1.0f` under
libc++ at a measured rate of 3.0e-8 while libstdc++ clamps the same case to the largest value
below one, placing an atom there. Since the log path takes a logarithm of that uniform, the
two builds were running different discrete samplers.

`bi_kappa_distribution` now draws from primitives defined in its own header — a documented
mapping from engine bits to the open interval `(0,1)`, a Marsaglia–Tsang Gamma, and a polar
normal — so its output is a function of the engine alone. libc++ and libstdc++ builds agree
bit for bit. Across architectures they do not, and cannot: `log`, `exp` and `pow` are not
correctly rounded and their implementations differ.

The other samplers in this repository still call `<random>` distributions directly and do not
have this property.

### Fixed — `seed(int)` did not restart the stream

`bi_kappa_distribution::seed(int)` reseeded the engine without calling `reset()`. Under
libstdc++ the variate cached inside `std::gamma_distribution`'s member normal survived the
reseed, so `seed(s)` followed by an *odd* number of draws did not reproduce. Under libc++ it
always did, which is why this went unnoticed. `test_seed_api` passed only because it happened
to compare after an even number of draws.

`seed()` now calls `reset()`, and `reset()` clears every cached variate. A regression test
checks reproducibility for one through eight draws rather than for a single convenient count.

### Added — observability of non-finite draws

```cpp
unsigned long long n_attempts()  const;
unsigned long long n_nonfinite() const;
```

Both are cleared by `reset()`, by every `define()` overload, and by `seed()`.

`n_nonfinite()` counts attempts whose intended velocity had no representation in `RealType`.
Untruncated mode returns those with non-finite components, exactly as before, and now counts
them. **Capped mode rejects them**, as it rejects anything outside the box, and still counts
them — and that is the case which was previously invisible, because every returned sample is
finite and the acceptance rate alone cannot distinguish "the box rejected a representable
sample" from "the box rejected one that could never have been returned at all".

### Added — version macros

`BI_KAPPA_VERSION_MAJOR`, `_MINOR`, `_PATCH` and `_STRING`, so code that vendors the header
can tell which sampler it has.

### Unchanged

The public API is otherwise identical: the same constructors, `define()` overloads, setters,
getters, `operator()`, `no_cap()` and `param_type`. The sampled law is unchanged in both
untruncated and capped modes. `rotate_from_fieldAligned_frame` is untouched.

## 1.0.0

Initial release. Header-only C++11 samplers for bi-Kappa and bi-Maxwellian velocity
distributions, the uncapped law as the default target, an opt-in component-wise velocity cap
selecting a distinct conditional law, and a validation battery for an N×3 velocity sample.
