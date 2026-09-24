# bi-kappa-distribution-sampling

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20617011.svg)](https://doi.org/10.5281/zenodo.20617011)

Header-only C++11 samplers for loading particle velocities and positions in kinetic plasma
simulations (particle-in-cell and hybrid codes). The package samples bi-Kappa and
bi-Maxwellian velocity distributions about an arbitrary magnetic-field direction, and draws
velocities or positions from densities you supply. A Python script tests a bi-Kappa velocity
sample from any loader against the target distribution.

- `cpp/`: the five sampler headers, a regression suite, and a demo program.
- `python/`: the validator `bikappa_validate.py`, Python versions of the general samplers,
  and a notebook that plots samples.
- `experiments/`: the studies behind the accompanying paper.

## Quick start (C++)

Copy the headers you need from `cpp/` into your project, or add `cpp/` to the include path.
Any C++11 compiler works.

```cpp
#include "bi_kappa_distribution.H"
#include <cstdio>

int main()
{
    using BK = bi_kappa_distribution<double>;
    BK dist;

    // kappa, theta_perp, theta_par, field direction ub (need not be unit length).
    // The velocity cap is 20 thermal speeds unless you pass one.
    dist.define(2.0, 1.0, 2.0, {0.0, 0.0, 1.0});
    dist.seed(12345);
    BK::point_type v = dist();          // {vx, vy, vz} in the global frame

    // The same distribution without the cap:
    dist.define(2.0, 1.0, 2.0, {0.0, 0.0, 1.0}, BK::no_cap());
    v = dist();

    std::printf("%g %g %g\n", v[0], v[1], v[2]);
    return 0;
}
```

```bash
g++ -std=c++11 -O2 -I cpp example.cpp -o example
```

To build and run the code shipped in `cpp/`:

```bash
make -C cpp test                  # build and run the regression suite (exit status 0 on success)
cd cpp && make && ./main.exe      # short worked examples of all five samplers
make -C cpp CXX=clang++ test      # use a different compiler
```

### The velocity cap

`bi_kappa_distribution` and `bi_maxwellian_distribution` take a cap
`max_normalized_velocity`, in thermal speeds. It is applied to each component in
field-aligned coordinates: a draw with `|v_i| / theta_i > cap` for any component is redrawn.
The samples therefore follow the distribution conditioned on that box, not the distribution
itself.

- The default cap is 20. Release 1.0.0 used the same default, so code written against 1.0.0
  samples the same distribution.
- Pass `no_cap()` (which is `+infinity`) to sample the bi-Kappa or bi-Maxwellian distribution
  itself. `dist.param().capped()` reports whether a finite cap is in force.
- For bi-Kappa with `kappa <= 3/2`, a cap of 20 removes an appreciable part of the high-energy
  tail; the paper quantifies this. For bi-Maxwellian the default cap has no practical effect.

Details, including the frame in which the box is applied, are on the
[parameters page](https://kehengzhu.github.io/bi-kappa-distribution-sampling/api/parameters.html).

## The five samplers

| Header / class | Samples | `define(...)` arguments |
|---|---|---|
| `bi_kappa_distribution` | bi-Kappa velocity | `kappa, theta_perp, theta_par, ub[, cap[, seed]]` |
| `bi_maxwellian_distribution` | bi-Maxwellian velocity | `theta_perp, theta_par, ub[, cap[, seed]]` |
| `general_velocity_generator` | speed from your `g(w)`, `w = \|v\|^2`; isotropic direction | `g, w_min, w_max[, probe_points, max_reject_tries]` |
| `field_aligned_velocity_generator` | parallel speed from your `g(w)`, `w = v_par^2`, with a fixed sign along `ub`; Maxwellian perpendicular components | `g, w_min, w_max, theta_perp, ub, sign[, probe_points, max_reject_tries, cap]` |
| `general_position_generator` | position in 1-3 dimensions from your density `rho(x)` | `dimension, lower, upper, rho[, probe_points, max_reject_tries]` |

Common behavior:

- Construct, call `define(...)`, then call `dist()` for one sample. `define(...)` throws
  `std::invalid_argument` on an invalid parameter.
- `theta_perp` and `theta_par` are thermal-speed parameters in units of velocity. For the
  bi-Maxwellian, `theta = sqrt(2 k_B T / m)` and each component is `N(0, theta^2 / 2)`. For the
  bi-Kappa distribution they are the scales in its density (Eq. 2 of the paper) and are not
  tied to a temperature, which does not exist for `kappa <= 3/2`.
- `ub` is the magnetic-field direction. Samples are returned in the global frame.
- `seed(s)` with `s >= 0` makes the stream reproducible; `s < 0` seeds from
  `std::random_device`. Without a seed, runs are not reproducible.
- `dist(gen)` draws from your own engine, for example a shared `std::mt19937`. For
  `bi_kappa_distribution` the engine must produce a power-of-two number of distinct values,
  as `std::mt19937` and `std::mt19937_64` do (`std::minstd_rand` does not).
- `general_velocity_generator` and `field_aligned_velocity_generator` take no particle mass.
  For a distribution written in energy, use `w = 2E/m`. If `w_min == w_max`, every sample has
  speed `sqrt(w_min)`. `probe_points` must be at least 64.
- The `field_aligned_velocity_generator` caps each perpendicular component at
  `cap * theta_perp` (default 20); `no_cap()` removes it.
- Near `kappa = 1/2`, some bi-Kappa draws exceed the largest representable floating-point
  number. Without a cap they are returned with infinite components; with a cap they are
  redrawn. `n_nonfinite()` and `n_attempts()` count them in both cases.
- Compiling `bi_kappa_distribution.H` with `-ffast-math` prints a note, because in that mode
  `no_cap()` and the overflow counters are unreliable. A finite cap still works. Define `BI_KAPPA_ALLOW_FAST_MATH` to
  silence the note.

## Validate a sample (Python)

```bash
cd python
uv sync          # or: pip install numpy scipy   (the notebook and figure scripts also
                 #     need matplotlib, pandas, mpmath and ipykernel)
```

`bikappa_validate.py` tests an `(N, 3)` velocity sample, with `N >= 1000`, from any bi-Kappa
loader against the bi-Kappa distribution. It uses no moments, so it also applies at
`kappa <= 3/2`. It tests the radius, the direction, their independence, the counts in cells
of equal probability (with the radius resolved into the tail as far as the sample size
allows), the anisotropy, and the number of non-finite values.

```bash
uv run python bikappa_validate.py samples.txt --kappa 2 --theta-perp 1 --theta-par 2
```

| Option | Meaning |
|---|---|
| `sample` | `.npy`, text (`.txt`, `.csv`), or raw float64 (`.bin`) file with three columns |
| `--kappa`, `--theta-perp`, `--theta-par` | parameters of the target distribution (required) |
| `--bhat BX BY BZ` | field direction; omit if the third column is already the parallel component |
| `--attempts N` | total draws attempted, if the loader redrew samples outside a cap; adds the rejected fraction to the report |
| `--alpha A` | significance level, default 0.01 |
| `--binary` | read the file as raw float64 whatever its extension |
| `--json` | print the full report as JSON |

The exit status is 0 if every test passes and 1 otherwise. The tests compare the sample with
the uncapped bi-Kappa distribution, so a sample from `bi_kappa_distribution` should be drawn
with `no_cap()`; the default cap of 20 changes the distribution.

From Python:

```python
from bikappa_validate import validate_sample, format_report
report = validate_sample(v, kappa=2.0, theta_perp=1.0, theta_par=2.0)   # v: (N, 3) array
print(format_report(report)); assert report["passed"]
```

Other Python files:

- `general_generators.py`: `GeneralVelocityGenerator`, `FieldAlignedVelocityGenerator` and
  `GeneralPositionGenerator`, with the same rules as the C++ classes. Call a generator with a
  `numpy.random.Generator`: `v = gen(rng)`. Running the file writes example samples.
- `visualize_samples.ipynb`: plots sample files written by the C++ or Python examples.
- `test_bikappa_validate.py`: regression test of the validator
  (`uv run python test_bikappa_validate.py`).

## Reproducing the paper

Each experiment directory has a README with its commands and outputs. Run the commands from
that directory.

| Directory | In the paper | Command |
|---|---|---|
| `experiments/exp1_radial_directional` | Tables III and IV, Figs. 4 and 5 | `make cells tail moments marginals` |
| `experiments/exp2_cap_characterization` | Table II, Fig. 3 | `make run`, then `uv run --project ../../python python exp2_analyze.py` |
| `experiments/exp4_finite_precision` | Fig. 2 | the run sequence in its README, ending with `make analyze figures` |
| `experiments/exp3_benchmark` | not shown (speed comparison with two other samplers) | `make run`, then `uv run --project ../../python python exp3_analyze.py` |

`paper/figures/make_manuscript_assets.py` regenerates the paper's figures and tables from the
committed experiment results.

## Documentation

- API reference: [kehengzhu.github.io/bi-kappa-distribution-sampling](https://kehengzhu.github.io/bi-kappa-distribution-sampling/),
  or open `docs/api/index.html` locally. Regenerate it with `doxygen Doxyfile` from the
  repository root.
- [Parameter reference](https://kehengzhu.github.io/bi-kappa-distribution-sampling/api/parameters.html):
  every argument, default and valid range.
- [Usage examples](https://kehengzhu.github.io/bi-kappa-distribution-sampling/api/usage.html) and
  [choosing a sampler](https://kehengzhu.github.io/bi-kappa-distribution-sampling/api/choosing.html).
- [CHANGELOG.md](https://github.com/KehengZhu/bi-kappa-distribution-sampling/blob/main/CHANGELOG.md): changes between releases, including what differs from 1.0.0.

## How to cite

Please cite the paper and the software:

- K. Zhu and Y. A. Omelchenko, "Verifying the High-Energy Tail in Bi-Kappa Particle Loading:
  Three-Dimensional Tests, Velocity Bounds, and Finite Precision" (2026).
- K. Zhu and Y. A. Omelchenko, bi-kappa-distribution-sampling, Zenodo,
  [doi:10.5281/zenodo.20617011](https://doi.org/10.5281/zenodo.20617011). This concept DOI
  resolves to the latest release; Zenodo lists the DOI of each version.

GitHub's "Cite this repository" button reads [`CITATION.cff`](CITATION.cff).

## License

MIT. See [`LICENSE`](LICENSE).
