# Experiment 1 results

Uncapped sampler (`no_cap()`), 5 replicate seeds x 10^5 samples per configuration,
double precision, Apple clang / libc++, `std::mt19937`.

`sqrt(n)D` is the KS statistic scaled by sqrt(n). It is O(1) -- around 0.87 on average --
when the law is correct, and does not shrink with n, unlike a p-value, which collapses to
0 for any large enough sample regardless of how small the discrepancy is.
Entries are mean +/- sd over the 5 replicates.

## Block A -- radial law and directional uniformity (theta_perp:theta_par = 1:2, B || z)

| kappa | non-finite / 5x10^5 | radial sqrt(n)D | CvM | cos(theta) sqrt(n)D | phi sqrt(n)D | indep. chi2 p | direction cells chi2 p | max |dlog R| |
|---|---|---|---|---|---|---|---|---|
| 0.51 | 0 | 0.900 +/- 0.160 | 0.158 | 0.964 +/- 0.129 | 0.778 +/- 0.072 | 0.708 | 0.632 | 8.690 |
| 0.55 | 0 | 0.889 +/- 0.356 | 0.160 | 1.067 +/- 0.282 | 0.706 +/- 0.128 | 0.668 | 0.529 | 1.390 |
| 0.75 | 0 | 0.792 +/- 0.473 | 0.151 | 0.877 +/- 0.253 | 0.639 +/- 0.096 | 0.624 | 0.583 | 0.336 |
| 1 | 0 | 0.867 +/- 0.353 | 0.201 | 0.900 +/- 0.226 | 0.593 +/- 0.069 | 0.334 | 0.485 | 0.179 |
| 1.25 | 0 | 0.771 +/- 0.251 | 0.133 | 0.878 +/- 0.287 | 0.610 +/- 0.048 | 0.306 | 0.533 | 0.134 |
| 1.5 | 0 | 0.652 +/- 0.269 | 0.094 | 0.994 +/- 0.160 | 0.709 +/- 0.159 | 0.301 | 0.591 | 0.083 |
| 2 | 0 | 0.728 +/- 0.228 | 0.089 | 0.901 +/- 0.213 | 0.655 +/- 0.120 | 0.386 | 0.628 | 0.046 |
| 5 | 0 | 1.043 +/- 0.521 | 0.314 | 0.843 +/- 0.169 | 0.788 +/- 0.139 | 0.326 | 0.595 | 0.019 |
| 10 | 0 | 1.106 +/- 0.477 | 0.299 | 0.837 +/- 0.199 | 0.841 +/- 0.196 | 0.491 | 0.582 | 0.013 |

## Cell-count test (the manuscript's test), replicates pooled

Pearson chi^2 against equal occupancy of 10 equal-probability radius shells x 40 equal-solid-angle
direction cells (5 cos(theta) x 8 phi intervals). `radius` sums over direction and tests the radial law alone;
`all cells` tests radius, direction and their independence together.

| block | kappa | theta_par | B direction | draws | expected per cell | min cell | radius p | all cells p |
|---|---|---|---|---|---|---|---|---|
| A | 0.51 | 2 | z | 500000 | 1250 | 1132 | 0.353 | 0.674 |
| A | 0.55 | 2 | z | 500000 | 1250 | 1122 | 0.019 | 0.043 |
| A | 0.75 | 2 | z | 500000 | 1250 | 1144 | 0.040 | 0.117 |
| A | 1 | 2 | z | 500000 | 1250 | 1141 | 0.456 | 0.312 |
| A | 1.25 | 2 | z | 500000 | 1250 | 1120 | 0.983 | 0.036 |
| A | 1.5 | 2 | z | 500000 | 1250 | 1154 | 0.664 | 0.566 |
| A | 2 | 2 | z | 500000 | 1250 | 1147 | 0.806 | 0.707 |
| A | 5 | 2 | z | 500000 | 1250 | 1173 | 0.504 | 0.823 |
| A | 10 | 2 | z | 500000 | 1250 | 1152 | 0.126 | 0.475 |
| B | 0.55 | 1 | diag111 | 500000 | 1250 | 1122 | 0.019 | 0.043 |
| B | 0.55 | 1 | oblique | 500000 | 1250 | 1122 | 0.019 | 0.043 |
| B | 0.55 | 1 | z | 500000 | 1250 | 1122 | 0.019 | 0.043 |
| B | 0.55 | 2 | diag111 | 500000 | 1250 | 1122 | 0.019 | 0.043 |
| B | 0.55 | 2 | oblique | 500000 | 1250 | 1122 | 0.019 | 0.043 |
| B | 0.55 | 2 | z | 500000 | 1250 | 1122 | 0.019 | 0.043 |
| B | 2 | 1 | diag111 | 500000 | 1250 | 1147 | 0.806 | 0.707 |
| B | 2 | 1 | oblique | 500000 | 1250 | 1147 | 0.806 | 0.707 |
| B | 2 | 1 | z | 500000 | 1250 | 1147 | 0.806 | 0.707 |
| B | 2 | 2 | diag111 | 500000 | 1250 | 1147 | 0.806 | 0.707 |
| B | 2 | 2 | oblique | 500000 | 1250 | 1147 | 0.806 | 0.707 |
| B | 2 | 2 | z | 500000 | 1250 | 1147 | 0.806 | 0.707 |
| B | 10 | 1 | diag111 | 500000 | 1250 | 1152 | 0.126 | 0.475 |
| B | 10 | 1 | oblique | 500000 | 1250 | 1152 | 0.126 | 0.475 |
| B | 10 | 1 | z | 500000 | 1250 | 1152 | 0.126 | 0.475 |
| B | 10 | 2 | diag111 | 500000 | 1250 | 1152 | 0.126 | 0.475 |
| B | 10 | 2 | oblique | 500000 | 1250 | 1152 | 0.126 | 0.475 |
| B | 10 | 2 | z | 500000 | 1250 | 1152 | 0.126 | 0.475 |

## Orientation of the bounded radial diagnostic

Both W = 1/(1+T) ~ Beta(kappa-1/2, 3/2) and Y = T/(1+T) ~ Beta(3/2, kappa-1/2) are exact
bijections of T and carry identical information in exact arithmetic. In doubles they do
not. The table above uses W; the same test on Y is shown here for comparison.

| kappa | radial sqrt(n)D using W | radial sqrt(n)D using Y | Y values rounded to exactly 1 | W values underflowed to 0 |
|---|---|---|---|---|
| 0.51 | 0.900 | 220.437 | 348541 / 500000 (69.71%) | 293 / 500000 |
| 0.55 | 0.889 | 51.937 | 82120 / 500000 (16.42%) | 0 / 500000 |
| 0.75 | 0.792 | 0.792 | 55 / 500000 (0.01%) | 0 / 500000 |
| 1 | 0.867 | 0.867 | 0 / 500000 (0.00%) | 0 / 500000 |
| 1.25 | 0.771 | 0.771 | 0 / 500000 (0.00%) | 0 / 500000 |
| 1.5 | 0.652 | 0.652 | 0 / 500000 (0.00%) | 0 / 500000 |
| 2 | 0.728 | 0.728 | 0 / 500000 (0.00%) | 0 / 500000 |
| 5 | 1.043 | 1.043 | 0 / 500000 (0.00%) | 0 / 500000 |
| 10 | 1.106 | 1.106 | 0 / 500000 (0.00%) | 0 / 500000 |

## Block B -- anisotropy and arbitrary field frame

| kappa | theta_par | B direction | radial sqrt(n)D | MAD ratio (expect theta_par/theta_perp) | basis orthonormality err | max off-diag corr |
|---|---|---|---|---|---|---|
| 0.55 | 1 | diag111 | 0.889 +/- 0.356 | 1.0096 (expect 1) | 4.4e-16 | n/a |
| 0.55 | 1 | oblique | 0.889 +/- 0.356 | 1.0096 (expect 1) | 2.2e-16 | n/a |
| 0.55 | 1 | z | 0.889 +/- 0.356 | 1.0096 (expect 1) | 0.0e+00 | n/a |
| 0.55 | 2 | diag111 | 0.889 +/- 0.356 | 2.0192 (expect 2) | 4.4e-16 | n/a |
| 0.55 | 2 | oblique | 0.889 +/- 0.356 | 2.0192 (expect 2) | 2.2e-16 | n/a |
| 0.55 | 2 | z | 0.889 +/- 0.356 | 2.0192 (expect 2) | 0.0e+00 | n/a |
| 2 | 1 | diag111 | 0.728 +/- 0.228 | 0.9996 (expect 1) | 4.4e-16 | 3.94e-02 |
| 2 | 1 | oblique | 0.728 +/- 0.228 | 0.9996 (expect 1) | 2.2e-16 | 3.94e-02 |
| 2 | 1 | z | 0.728 +/- 0.228 | 0.9996 (expect 1) | 0.0e+00 | 3.94e-02 |
| 2 | 2 | diag111 | 0.728 +/- 0.228 | 1.9992 (expect 2) | 4.4e-16 | 3.94e-02 |
| 2 | 2 | oblique | 0.728 +/- 0.228 | 1.9992 (expect 2) | 2.2e-16 | 3.94e-02 |
| 2 | 2 | z | 0.728 +/- 0.228 | 1.9992 (expect 2) | 0.0e+00 | 3.94e-02 |
| 10 | 1 | diag111 | 1.106 +/- 0.477 | 0.9985 (expect 1) | 4.4e-16 | 6.12e-03 |
| 10 | 1 | oblique | 1.106 +/- 0.477 | 0.9985 (expect 1) | 2.2e-16 | 6.12e-03 |
| 10 | 1 | z | 1.106 +/- 0.477 | 0.9985 (expect 1) | 0.0e+00 | 6.12e-03 |
| 10 | 2 | diag111 | 1.106 +/- 0.477 | 1.9969 (expect 2) | 4.4e-16 | 6.12e-03 |
| 10 | 2 | oblique | 1.106 +/- 0.477 | 1.9969 (expect 2) | 2.2e-16 | 6.12e-03 |
| 10 | 2 | z | 1.106 +/- 0.477 | 1.9969 (expect 2) | 0.0e+00 | 6.12e-03 |

## Frame invariance (direct comparison)

For matched (kappa, theta, seed), the normalized radius recovered after rotating into an
arbitrary field direction is compared against the axis-aligned run draw by draw.

- comparisons: 60 run pairs, 6000000 draws
- bitwise identical radii: 1374640 / 6000000 (22.9107%)
- largest relative difference: 1.16e-15 (kappa=2, diag111, seed=1005)

Whole normalized velocity: every block-B run against the (theta_par = 2, B || z) run
with the same kappa and seed, after undoing each run's own scaling and rotation.

- comparisons: 75 run pairs, 7500000 draws
- largest relative difference |u - u_ref| / |u_ref|: 8.97e-16

Notes.

- Independence: chi^2 contingency of 4 radial quartiles x 40 equal-solid-angle direction
  cells (5 cos(theta) x 8 phi intervals). The direction-cell column tests the same cell
  counts, summed over radius, against equal occupancy: the joint direction law.
- Variance-based checks are reported only for kappa > 3/2, where the second moment
  exists. Below that the MAD ratio and the radial/quantile diagnostics carry the test.
- `max off-diag corr` is the largest off-diagonal correlation of the sample covariance
  in the recovered field-aligned frame; a misaligned rotation would inflate it. It is
  itself noisy at kappa = 2, where the fourth moment does not exist and the correlation
  estimator therefore has no finite variance.
