# Experiment 2 — effect of the component-wise velocity cap

`bi_kappa_distribution` has two modes, and they sample two different distributions.

| mode | `max_normalized_velocity` | distribution sampled |
|---|---|---|
| uncapped | `no_cap()` (= `+inf`) | the bi-Kappa distribution |
| capped | finite `λ` | the bi-Kappa distribution conditioned on the box `\|v_x\|/θ⊥ ≤ λ`, `\|v_y\|/θ⊥ ≤ λ`, `\|v_z\|/θ∥ ≤ λ` |

The capped distribution has bounded support, a truncated tail and different moments. The box
is a cube in normalized components, so the capped distribution is neither isotropic nor
axisymmetric about **B**. This experiment measures how far the capped distribution lies from
the uncapped one.

## What it computes

- **Rejected fraction.** Conditioning on an event of probability `p` multiplies the density by
  `1/p` inside the box, so the total-variation (TV) distance between the capped and uncapped
  distributions equals the rejected fraction `1 − p`. The experiment measures it empirically and
  also evaluates a closed form (below). It also records the distribution of attempts per
  accepted draw and compares it with `Geometric(p)`.
- **Tail fidelity.** The ratio of capped to uncapped speed quantiles at the 50th, 90th, 99th and
  99.9th percentiles, at matched seeds, and the largest gap between the two empirical CDFs of the
  speed.
- **Exact conditioning.** For every capped run, the capped output is compared with the draws of
  the uncapped run at the same seed that fall inside the box. The two sequences are bitwise
  identical in all 240 (case, λ, seed) pairs. This check works because each iteration of the
  sampler's internal loop consumes the same random variates whether or not a cap is set.
- **Independence of θ.** The θ values cancel in the normalized cap test, so the accept/reject
  decision is predicted to be the same for (θ⊥, θ∥) = (1, 1) and (1, 2) at the same seed. It
  agrees on all 6 × 10⁶ attempts compared.
- **Moments.** The capped-to-uncapped variance ratio for κ = 2, 5, 10. For κ ≤ 3/2 the
  uncapped second moment does not exist, so no ratio is reported. The analysis instead shows
  that the variance of the capped sample grows without bound as λ increases.
- **Angular structure.** The four-fold azimuthal Fourier coefficient `a4 = 2⟨cos 4φ⟩` and a KS
  test of `cos θ`, which measure how far the capped distribution departs from axisymmetry.

The closed form for the acceptance probability uses the fact that the θ values cancel. The box
is then a cube of half-side `c = λ/√κ` in isotropic coordinates, and with `M = maxᵢ|nᵢ|` for a
unit vector `n` uniform on the sphere,

```
P(accept) = E_M[ I_z(3/2, κ−1/2) ],   z = c²/(M² + c²),
f_M(m) = 3 − (12/π) arcsin( √((1−2m²)/(1−m²)) )   for 1/√3 ≤ m ≤ 1/√2,
f_M(m) = 3                                        for 1/√2 ≤ m ≤ 1.
```

The analysis checks `f_M` against a spherical Monte Carlo sample of 2 × 10⁶ points. It computes
the rejected fraction directly, through `1 − I_z(a, b) = I_{1−z}(b, a)`, rather than as
`1 − P(accept)`, because at large κ and λ that subtraction would lose every significant digit.

## Configuration

The field is along `ẑ` in every run. The cap is applied before the rotation into the global
frame, so the acceptance probability does not depend on the field direction.

| | anisotropic runs | isotropic control |
|---|---|---|
| `block` column in `raw/manifest.csv` | `A` | `C` |
| κ | 0.75, 1, 1.5, 2, 5, 10 | 0.75, 2 |
| (θ⊥, θ∥) | (1, 2) | (1, 1) |
| λ | `no_cap()`, 3, 5, 10, 20, 50, 100 | `no_cap()`, 3, 5, 10, 20, 50, 100 |
| seeds | 2001–2005 | 2001–2005 |
| draws per run | 10⁵ | 10⁵ |
| runs | 210 | 70 |

In total there are 280 runs and 2.8 × 10⁷ draws, all in double precision. λ = 20 is the default
cap of release 1.0.0.

## Paper

`paper/figures/make_manuscript_assets.py` draws Fig. 3 and writes Table II from
`results/exp2_results.json`. These show the rejected fraction (closed form) and the ratio of
the 99.9th-percentile speeds for the anisotropic runs.

## Rerunning

Run from this directory.

```bash
make run                                               # build, write raw/, write raw/checksums.sha256
uv run --project ../../python python exp2_analyze.py   # read raw/, write results/
make verify                                            # check raw/ against raw/checksums.sha256
```

`exp2_analyze.py` also accepts the input and output directories as arguments,
`exp2_analyze.py [raw_dir] [results_dir]`. `make clean` removes the executable and
`make distclean` also removes `raw/`.

**Runtime and disk use.** On an Apple M4 Max, `exp2_sample.exe` takes about 4 s and the analysis
about 20 s. `raw/` holds 280 binary files of 2.4 MB each, about 670 MB in total. These files are
not committed. `raw/manifest.csv` and `raw/checksums.sha256` are committed, so a regenerated
`raw/` can be checked against them.

## Output files

| file | contents |
|---|---|
| `raw/run_NNNN.bin` | draws of one run, 3 doubles per draw (not committed) |
| `raw/manifest.csv` | one row per run: κ, θ, field direction, mode, λ, seed, draw count, non-finite count, file name |
| `raw/checksums.sha256` | SHA-256 of `manifest.csv` and every `run_*.bin` |
| `results/exp2_results.json` | all measured quantities per configuration, the conditioning and θ checks, and the build environment |
| `results/exp2_table.md` | the same results as Markdown tables |

## Committed results

The committed results were produced with release 2.2.1 of `cpp/bi_kappa_distribution.H` on
macOS 26.7 on arm64 (Apple silicon), compiled by Apple clang 21 with libc++ and
`-Wall -Wextra -std=c++11 -O2`.
`results/exp2_results.json` records the SHA-256 of the header (`sampler_header_sha256`,
beginning `79842dc4`), the compiler, the git commit and the Python package versions.

The main results are these:

1. The rejected fraction decays as `λ^{−(2κ−1)}`. The measured slope between λ = 50 and λ = 100
   agrees with this exponent to three digits. At κ = 0.75 a cap of λ = 100 still rejects 9.7% of
   attempts, which is also the TV distance from the bi-Kappa distribution.
2. A small TV distance does not imply a preserved tail. At κ = 1.5 and λ = 50 the TV distance is
   6.3 × 10⁻⁴, but the 99.9th-percentile speed of the capped distribution is 20% below that of
   the uncapped one.
3. Take a cap as negligible when TV < 10⁻³ and the 99.9th-percentile speed is within 1% of the
   uncapped value. No λ in the ladder meets both conditions for κ ≤ 3/2. The smallest λ that
   does is 50 at κ = 2, 10 at κ = 5 and 5 at κ = 10.
4. The capped distribution has a four-fold azimuthal modulation about **B**. At κ = 0.75 and
   λ = 3 the coefficient `a4` is 11.6 standard deviations from zero.
5. The number of attempts per accepted draw follows `Geometric(p)`. The sampler's limit of 10⁶
   consecutive rejections is never approached: in the worst configuration the base-10 logarithm
   of the probability of reaching it is about −2.6 × 10⁵.

## Numerical notes

- The speed `|v|` is computed with nested `np.hypot`. `np.linalg.norm` squares internally and
  overflows in the κ ≤ 1 tails.
- `R²` is never formed. The closed form uses `z = c²/(m² + c²)`.
- The analysis evaluates the cap predicate `in_box` on the uncapped draws. It is a copy of the
  header's `withinNormalizedVelocityCap`. The header itself is used unmodified.
