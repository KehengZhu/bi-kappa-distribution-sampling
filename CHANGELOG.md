# Changelog

This project follows [Semantic Versioning](https://semver.org/). A change in the random numbers
that a given seed produces counts as a breaking change.

## 2.2.1 — 2026-09-24

Changes since 1.0.0, the previous release. Versions 2.0.0 to 2.2.0 were never released.

### Breaking changes

- **bi-Kappa samples differ from 1.0.0 for every seed.** The radius is now formed on a
  logarithmic scale, so draws near `kappa = 1/2` are no longer lost when the
  `Gamma(kappa - 1/2)` variate underflows. The direction is drawn by Marsaglia's method. The
  sampled distribution is unchanged. `bi_maxwellian_distribution` returns the same values as
  in 1.0.0.
- **`general_velocity_generator` takes a speed-squared density.** The 1.0.0 interface
  `define(f, E_min, E_max, particle_mass, ...)`, with an energy density `f(E)`, is replaced by
  `define(g, speed_sq_min, speed_sq_max, probe_points = 1024, max_reject_tries = 100000)`, with
  `g(w)` and `w = |v|^2`. An energy density `f(E)` becomes `g(w) = f(m w / 2)`. An old four-argument call
  passes the mass as `probe_points`, and throws if the mass is below 64. The Python
  `GeneralVelocityGenerator` changed in the same way.
- A random engine passed to `bi_kappa_distribution::operator()(gen)` must produce a
  power-of-two number of distinct values, as `std::mt19937` and `std::mt19937_64` do.
  `std::minstd_rand` does not.

### Behavior changes

- `seed()` now also calls `reset()` on every sampler, so `seed(s)` followed by `n` draws gives
  the same values however many draws came before it. In 1.0.0, reseeding an object that had
  already drawn could reuse a normal variate cached from before. A fresh object, and
  `define(..., seed)`, give the same bi-Maxwellian values as in 1.0.0.
- `bi_kappa_distribution::define()` rejects NaN parameters, infinite thermal speeds,
  non-finite `ub` components and a zero `ub` when it is called.
- A bi-Kappa draw that overflows the floating-point range sets `FE_OVERFLOW` and
  `errno = ERANGE`.
- For bi-Kappa, the same seed gives the same values under libc++ and libstdc++ on the same
  architecture (compile with `-ffp-contract=off` to rule out differences from fused
  multiply-add). With `float`, each component is computed in `double` and rounded once.
- `general_velocity_generator` and `field_aligned_velocity_generator` require
  `probe_points >= 64`. If `speed_sq_min == speed_sq_max`, every sample has speed
  `sqrt(speed_sq_min)`.
- Compiling `bi_kappa_distribution.H` with `-ffast-math` prints a note, because `no_cap()` and
  the overflow counters are unreliable in that mode. Define `BI_KAPPA_ALLOW_FAST_MATH` to
  silence it.

### Added

- `no_cap()` on `bi_kappa_distribution` and `bi_maxwellian_distribution` samples the
  distribution without a velocity cap; `param_type::capped()` reports whether a finite cap is
  in force.
- `n_nonfinite()` and `n_attempts()` on `bi_kappa_distribution` count draws that overflow and
  all attempts, including draws rejected by the cap.
- `field_aligned_velocity_generator`: parallel speed from a user density `g(v_par^2)` with a
  fixed sign along `ub`, and Maxwellian perpendicular components capped at
  `max_normalized_velocity * theta_perp` (default 20).
- Version macros `BI_KAPPA_VERSION_MAJOR`, `_MINOR`, `_PATCH` and `_STRING`.
- `python/bikappa_validate.py`: tests an `(N, 3)` velocity sample from any bi-Kappa loader
  against the target distribution without using moments.
- `FieldAlignedVelocityGenerator` in `python/general_generators.py`, the Python version of
  `field_aligned_velocity_generator`.

### Unchanged

- The velocity cap `max_normalized_velocity` still defaults to 20 thermal speeds on
  `bi_kappa_distribution` and `bi_maxwellian_distribution`, as in 1.0.0.

## 1.0.0

Initial release. Header-only C++11 samplers for bi-Kappa and bi-Maxwellian velocity
distributions about an arbitrary field direction, with a component-wise velocity cap of 20
thermal speeds by default; a general velocity generator for a user-defined energy density; a
general position generator; Python versions of the two general generators; and notebooks that
plot samples.
