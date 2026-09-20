# Changelog

This project follows [Semantic Versioning](https://semver.org/). For a sampler, "breaking"
includes changing which random numbers a given seed produces, even when the API is untouched
and the sampled law is unchanged — a stored stream is part of what a user depends on.

## 2.1.0 — unreleased

### Fixed — representability is decided from the component, not from its logarithm

`bi_kappa_distribution` decided whether a velocity component could be returned by comparing a
logarithm — `log R`, or `log R + log|g_j|` — against `log(max())`. That comparison cannot be
made to agree with the arithmetic it predicts. `log(max())` is itself a rounded value and can
land on either side of the true logarithm of the largest finite number; on the libm this was
tested against it lands *above* in `float`, where `exp(log(FLT_MAX))` is exactly `+inf`, and
below in `double`. So the defect's visibility was a property of a rounding direction rather
than of the type.

A `float` draw at `kappa = 0.51` whose `log R` equalled `log(FLT_MAX)` bit for bit therefore
passed the test, overflowed when it was exponentiated, and was returned as
`(-inf, +inf, +inf)` — while all three components it should have produced were representable,
the largest at 1.92e38 against a limit of 3.40e38, because every `|g_j|` is below one and
pulls the product back under the limit. Worse, `n_nonfinite()` was **not** incremented: the
branch that produced the infinities never looked at what it had produced, so the loss was
invisible to the caller. One such draw appeared in 14 233 536 audited attempts.

The decision now lives in one place, `bikappa_detail::materializeComponents`, which
materializes each component and compares the component with `max()`. The test is applied to
the number that will be returned, so this invariant holds by construction and in both modes:

> the returned vector has a non-finite component **if and only if** `n_nonfinite()` counted it.

Two consequences for anyone who had pinned 2.0.0:

- **Returned values change, for draws 2.0.0 got wrong.** Where `exp(log R)` overflows but the
  components do not, 2.0.0 returned infinities and 2.1.0 returns the numbers. The random
  stream — which variates are drawn, and in what order — is unchanged, so every other draw is
  bit-identical.
- **`n_nonfinite()` changes under a cap.** It counted the representability of the *normalized*
  coordinate the cap predicate is written in; it now counts the representability of the
  velocity, which is what its documentation always said. The two differ by `theta` and by the
  rotation — up to 2.33× on the geometry this was measured on.

2.0.0 is superseded before release and should not be used: it and 2.1.0 are two different
samplers, and giving them one version string would make them indistinguishable to a consumer
that vendors the header.

### Changed — the header refuses to compile under `-ffast-math`

The representability test asks whether a materialized value is finite. `-ffast-math`, and
`-ffinite-math-only` on its own, license a compiler to assume no infinity or NaN ever occurs
and so to fold that test to `true` — which would return an overflowed component as if it were
a number *and* leave the counters reading zero, reinstating exactly the defect above with no
symptom. There is nothing to check at run time once the check has been compiled away, so the
header stops the build instead. Define `BI_KAPPA_ALLOW_FAST_MATH` to override; the overflow
counters then mean nothing.

### Changed — parameters are validated against NaN

`kappa`, both `theta`s and `ub` are now rejected when they are NaN or non-finite, and `ub` is
checked when the parameters are set rather than when the rotation first runs. The old
spellings (`kappa <= 0.5`, `theta <= 0`) are false for a NaN, so a NaN parameter used to
propagate into `log R`, where it could no longer be told from an unrepresentable velocity: the
sampler would report a non-finite draw and count it as honest overflow, which it is not.

### Changed — floating-point environment

The sampler now evaluates `exp(log R)` on every draw, so a draw whose radius genuinely
overflows raises `FE_OVERFLOW` and sets `errno` to `ERANGE` where 2.0.0 suppressed both by
testing a logarithm first. Near `kappa = 1/2` that is most draws — at double `kappa = 0.5001`
the representability floor is 0.868 — and it is the honest signal, not noise: those draws
really do have no representation. A caller that inspects `fetestexcept` or `errno` around the
sampler will see state it did not see before.

## 2.0.0 — unreleased, superseded by 2.1.0 before release

### Changed — the bi-Kappa radius is built in the logarithmic domain

`bi_kappa_distribution` draws its radius as `R = sqrt(X1)/sqrt(X2)` with
`X1 ~ Gamma(3/2, 1)` and `X2 ~ Gamma(kappa - 1/2, 1)`. As `kappa` approaches `1/2` the
denominator's shape approaches zero, `X2` underflows to exactly zero, and the draw is lost —
a loss no rearrangement of the quotient can survive, because the information is gone before
the division happens.

`X2` is therefore no longer formed. The sampler carries `log X2`, obtained from the
Ahrens–Dieter shape-boosting identity, propagates `log R = (log X1 - log X2)/2`, builds the
order-unity vector `g` with `V = R g` before exponentiating anything, and decides whether each
returned component is representable — by a logarithmic test that 2.1.0 replaces, for the
reasons given above. The velocity cap is
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
normal.

That was not by itself enough. The direction was drawn as `cos θ` uniform on `(−1,1)` and
`φ` uniform on `(0, 2π)`, which needs `sin(φ)` and `cos(φ)` of the same argument — and a
compiler is free to fuse that pair into one routine. Clang does, on Darwin, into
`__sincos_stret`, whose sine differs from the standalone `sin` by one unit in the last place
for roughly one argument in a thousand; GCC calls the two separately. Two builds of the same
header therefore returned different numbers for the same seed, and the difference belonged to
the compiler rather than to the sampler.

The direction is now drawn by the rejection method of Marsaglia (1972): a point uniform in
the unit disc, lifted to the sphere. It uses only `sqrt`, which IEEE-754 requires to be
correctly rounded, so it calls no library transcendental at all. It is also faster here —
12.1 ns against 16.3 under clang, 10.2 against 18.9 under gcc — despite averaging `4/π`
attempts.

With that, the sampled sequence is a function of the engine alone: full-loader digests over
six `kappa` values, capped and uncapped, rotated and anisotropic, agree bit for bit between
libc++ and libstdc++ in both precisions. Two things it still depends on, both stated rather
than assumed: whether the compiler contracts a multiply and an add (pass `-ffp-contract=off`
to pin the stream to what the source text says), and the architecture — `log` and `exp` are
not correctly rounded and their implementations differ, so cross-architecture agreement is
statistical, not bitwise.

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
