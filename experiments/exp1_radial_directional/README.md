# Experiment 1 — three-dimensional validation of the uncapped sampler

This experiment tests the full three-dimensional output of `bi_kappa_distribution` in its
uncapped mode (`no_cap()`), where the sampled distribution is the bi-Kappa distribution itself.
The sampler is built from a radius `R` with `R² ~ BetaPrime(3/2, κ − 1/2)`, a uniform direction
independent of `R`, and a scaling and rotation into the field-aligned frame. The tests check the
radius, the direction, their independence and the frame mapping directly, so they apply for
`1/2 < κ ≤ 3/2`, where the second moment does not exist.

All tests use the same draws. At each of nine values κ ∈ {0.51, 0.55, 0.75, 1, 1.25, 1.5, 2, 5,
10} the sampler produces 100 seeds × 5 × 10⁵ = 5 × 10⁷ velocities in double precision, with
(θ⊥, θ∥) = (1, 2) and the field along `ẑ`. The κ at index `i` of that list uses seeds
`100001 + 1000 i` to `100100 + 1000 i`.

## What each script computes

| script | make target | computes | paper |
|---|---|---|---|
| `exp1_cells.py` | `make cells` | Pearson χ² over 400 cells of equal probability, 10 radius shells × 40 direction cells, at all nine κ | Table III (400-cell column) |
| `exp1_tail.py` | `make tail` | Pearson χ² over 14 radius shells at all nine κ | Table III (radius column), Fig. 4 |
| `exp1_moments.py` | `make moments` | variance of each velocity component at κ = 2, 5, 10 | Table IV |
| `exp1_frame.py` | `make frame` | invariance under rescaling and rotation at κ = 0.55, 2, 10 | Sec. VI C |
| `exp1_marginals.py` | `make marginals` | histograms of the three velocity components at κ = 1, 2, 10 | Fig. 5 |

Two further source files support these scripts. `exp1_draw.cpp` writes the velocities of one
seed at one κ to a binary file; every script calls it once per seed. `exp1_common.py` defines
the radius shells, the direction cells, and a Python copy of the sampler's field-aligned basis.

**Cell test (`exp1_cells.py`).** The ten radius shells are bounded by the exact deciles of `R`.
The 40 direction cells cross five equal intervals of `cos Θ` with eight equal intervals of `Φ`,
so they have equal solid angle. Because the radius and the direction are independent under the
target distribution, each of the 400 cells has probability 1/400, and 125 000 draws are expected
in each. The script reports two p-values: the ten shells alone, and all 400 cells together, which
tests the radius, the direction and their independence at once. It also compares the scatter of
the per-seed χ² values with the multinomial prediction. At κ = 0.51, 39 draws have a component
larger than the largest double (34 are expected from the exact tail). They lie beyond every shell
edge, so the radius test counts them in the outermost shell and the 400-cell test omits them.

**Radius test resolved into the tail (`exp1_tail.py`).** The outermost of the ten shells,
`R > r₉`, holds the fastest 10% of draws. The script divides it further at the 99th, 99.9th,
99.99th and 99.999th percentiles of `R`, which gives 14 shells and 13 degrees of freedom, with
500 draws expected beyond the last edge. Near κ = 1/2 these edges lie outside the double range
of `W = 1/(1 + R²)` (the 99.999th percentile is `R ≈ 10²⁵⁰` at κ = 0.51). The script therefore
computes them in log space from the leading term of the incomplete-Beta series, and checks that
form against the Beta quantile wherever both can be evaluated. It also recounts the 400-cell
counts and stops if they differ from `results/exp1_cells.json`, which confirms that it reads the
same draws.

**Variances (`exp1_moments.py`).** For each of the 100 seeds the script computes the sample
variance of each component. It reports the mean of the 100 values and its standard error, and
compares them with `θ² κ/(2κ − 3)`. The variance exists only for κ > 3/2, so the test is limited
to κ = 2, 5 and 10.

**Frame invariance (`exp1_frame.py`).** On the same seeds the script redraws each sample with
(θ⊥, θ∥) ∈ {(1, 1), (1, 2)} and three field directions, `ẑ`, `(1, 1, 1)/√3` and
`(0.3, −0.5, 0.8)/‖·‖`. Each redrawn sample is mapped back to normalized field-aligned
coordinates with the basis rebuilt in `exp1_common.field_basis`. It is then compared draw by
draw with the (1, 2), `ẑ` sample at the same seed. The script reports the largest relative
difference over 7.5 × 10⁸ comparisons.

**Component histograms (`exp1_marginals.py`).** Histograms of `v_i/θ_i` on 48 bins of width
0.25 over [−6, 6].

The moments and marginals scripts also recount the 400-cell counts and record whether they match
`results/exp1_cells.json` (`same_draws_as_cell_test` in the JSON output). In the frame script the
reference sample, (1, 2) with the field along `ẑ`, is the cell test's sample itself.

## Rerunning

Run from this directory. Python dependencies come from the project in `../../python` through
`uv`. Each target builds `exp1_draw.exe` first if it is missing or out of date.

```bash
make cells        # must run first: the other four targets read results/exp1_cells.json
make tail
make moments
make frame
make marginals
```

Each target runs `uv run --project ../../python python <script>.py`. To check the 14-shell
test at selected κ without writing any results, pass the κ values to the tail script:

```bash
uv run --project ../../python python exp1_tail.py 0.51 0.55
```

`make clean` removes `exp1_draw.exe`.

**Runtime and disk use.** One seed (5 × 10⁵ draws, written and counted) takes about 0.1 s on an
Apple M4 Max, so each target should take a few minutes. `frame` makes six draws per seed and runs
in parallel on half the available cores. The scripts write each seed's draws (12 MB) to a
temporary directory and delete them after counting, so the only files kept are those in
`results/`.

## Output files

| file | contents |
|---|---|
| `results/exp1_cells.json`, `results/exp1_cells.md` | per-κ p-values, overflow counts, per-seed scatter, shell deviations, and the pooled 10 × 40 cell counts |
| `results/exp1_tail.json`, `results/exp1_tail.md` | 14-shell edges, probabilities, counts, deviations and p-values |
| `results/exp1_moments.json`, `results/exp1_moments.md` | per-seed variances, their mean and standard error, and the deviation from theory |
| `results/exp1_frame.json`, `results/exp1_frame.md` | largest relative difference for each (κ, θ, field direction) combination |
| `results/exp1_marginals.json` | histogram counts of the three components |

`paper/figures/make_manuscript_assets.py` draws Figs. 4 and 5 and writes Tables III and IV from
these JSON files.

## Committed results

In the committed results every p-value in Table III is at least 0.062. The largest shell
deviation among the 126 counts of the 14-shell test is 2.5 standard deviations. All nine
variances lie within 1.6 standard errors of theory. The largest relative difference in the frame
test is 9.97 × 10⁻¹⁶.

The committed results were produced with release 2.2.1 of `cpp/bi_kappa_distribution.H`, on
macOS 26.7 on arm64 (Apple silicon). Each JSON file records the SHA-256 of the header it compiled
(`sampler_header_sha256`, beginning `b6559a9d`), the git commit, and the Python, NumPy and SciPy
versions.

## Numerical notes

- `R² = T` is never formed. `T` overflows for `R > 1.3 × 10¹⁵⁴` while `R` is still
  representable, so any statistic formed from `T` would discard valid draws in the heavy tail.
  The scripts compute `R` with nested `np.hypot` rather than `np.linalg.norm`, which squares
  internally, and compare `log R` with log-space shell edges.
- The shell edges are computed from quantiles of `W = 1/(1 + T) ~ Beta(κ − 1/2, 3/2)`, not of
  `Y = T/(1 + T) = 1 − W`. The two are equivalent in exact arithmetic. In the upper tail,
  however, `W` is small and keeps its relative precision, while `Y` rounds to 1.
