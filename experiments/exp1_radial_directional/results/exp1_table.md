# Experiment 1 results

Uncapped sampler (`no_cap()`), 5 replicate seeds x 10^5 samples per configuration,
double precision, Apple clang / libc++, `std::mt19937`.

`sqrt(n)D` is the KS statistic scaled by sqrt(n). It is O(1) -- around 0.87 on average --
when the law is correct, and does not shrink with n, unlike a p-value, which collapses to
0 for any large enough sample regardless of how small the discrepancy is.
Entries are mean +/- sd over the 5 replicates.

## Block A -- radial law and directional uniformity (theta_perp:theta_par = 1:2, B || z)

| kappa | non-finite / 5x10^5 | radial sqrt(n)D | CvM | cos(theta) sqrt(n)D | phi sqrt(n)D | indep. chi2 p | max |dlog R| |
|---|---|---|---|---|---|---|---|
| 0.51 | 292 | 0.815 +/- 0.078 | 0.136 | 0.835 +/- 0.238 | 1.071 +/- 0.162 | 0.589 | 29.309 |
| 0.55 | 0 | 0.751 +/- 0.035 | 0.098 | 0.921 +/- 0.189 | 0.976 +/- 0.263 | 0.590 | 1.918 |
| 0.75 | 0 | 0.930 +/- 0.264 | 0.167 | 0.867 +/- 0.244 | 0.979 +/- 0.273 | 0.459 | 0.256 |
| 1 | 0 | 0.910 +/- 0.090 | 0.158 | 0.840 +/- 0.215 | 0.975 +/- 0.442 | 0.486 | 0.082 |
| 1.25 | 0 | 0.841 +/- 0.150 | 0.126 | 0.896 +/- 0.381 | 0.994 +/- 0.295 | 0.501 | 0.084 |
| 1.5 | 0 | 0.970 +/- 0.283 | 0.211 | 0.938 +/- 0.270 | 0.763 +/- 0.157 | 0.413 | 0.033 |
| 2 | 0 | 0.890 +/- 0.185 | 0.143 | 0.813 +/- 0.222 | 0.773 +/- 0.285 | 0.592 | 0.058 |
| 5 | 0 | 0.873 +/- 0.243 | 0.145 | 0.780 +/- 0.184 | 0.940 +/- 0.301 | 0.563 | 0.020 |
| 10 | 0 | 0.891 +/- 0.209 | 0.164 | 0.805 +/- 0.149 | 0.843 +/- 0.298 | 0.703 | 0.016 |

## Orientation of the bounded radial diagnostic

Both W = 1/(1+T) ~ Beta(kappa-1/2, 3/2) and Y = T/(1+T) ~ Beta(3/2, kappa-1/2) are exact
bijections of T and carry identical information in exact arithmetic. In doubles they do
not. The table above uses W; the same test on Y is shown here for comparison.

| kappa | radial sqrt(n)D using W | radial sqrt(n)D using Y | Y values rounded to exactly 1 | W values underflowed to 0 |
|---|---|---|---|---|
| 0.51 | 0.815 | 220.379 | 348348 / 500000 (69.67%) | 2 / 500000 |
| 0.55 | 0.751 | 51.762 | 81843 / 500000 (16.37%) | 0 / 500000 |
| 0.75 | 0.930 | 0.930 | 62 / 500000 (0.01%) | 0 / 500000 |
| 1 | 0.910 | 0.910 | 0 / 500000 (0.00%) | 0 / 500000 |
| 1.25 | 0.841 | 0.841 | 0 / 500000 (0.00%) | 0 / 500000 |
| 1.5 | 0.970 | 0.970 | 0 / 500000 (0.00%) | 0 / 500000 |
| 2 | 0.890 | 0.890 | 0 / 500000 (0.00%) | 0 / 500000 |
| 5 | 0.873 | 0.873 | 0 / 500000 (0.00%) | 0 / 500000 |
| 10 | 0.891 | 0.891 | 0 / 500000 (0.00%) | 0 / 500000 |

## Block B -- anisotropy and arbitrary field frame

| kappa | theta_par | B direction | radial sqrt(n)D | MAD ratio (expect theta_par/theta_perp) | basis orthonormality err | max off-diag corr |
|---|---|---|---|---|---|---|
| 0.55 | 1 | diag111 | 0.751 +/- 0.035 | 1.0038 (expect 1) | 4.4e-16 | n/a |
| 0.55 | 1 | oblique | 0.751 +/- 0.035 | 1.0038 (expect 1) | 2.2e-16 | n/a |
| 0.55 | 1 | z | 0.751 +/- 0.035 | 1.0038 (expect 1) | 0.0e+00 | n/a |
| 0.55 | 2 | diag111 | 0.751 +/- 0.035 | 2.0077 (expect 2) | 4.4e-16 | n/a |
| 0.55 | 2 | oblique | 0.751 +/- 0.035 | 2.0077 (expect 2) | 2.2e-16 | n/a |
| 0.55 | 2 | z | 0.751 +/- 0.035 | 2.0077 (expect 2) | 0.0e+00 | n/a |
| 2 | 1 | diag111 | 0.890 +/- 0.185 | 0.9980 (expect 1) | 4.4e-16 | 1.72e-02 |
| 2 | 1 | oblique | 0.890 +/- 0.185 | 0.9980 (expect 1) | 2.2e-16 | 1.72e-02 |
| 2 | 1 | z | 0.890 +/- 0.185 | 0.9980 (expect 1) | 0.0e+00 | 1.72e-02 |
| 2 | 2 | diag111 | 0.890 +/- 0.185 | 1.9959 (expect 2) | 4.4e-16 | 1.72e-02 |
| 2 | 2 | oblique | 0.890 +/- 0.185 | 1.9959 (expect 2) | 2.2e-16 | 1.72e-02 |
| 2 | 2 | z | 0.890 +/- 0.185 | 1.9959 (expect 2) | 0.0e+00 | 1.72e-02 |
| 10 | 1 | diag111 | 0.891 +/- 0.209 | 0.9978 (expect 1) | 4.4e-16 | 5.27e-03 |
| 10 | 1 | oblique | 0.891 +/- 0.209 | 0.9978 (expect 1) | 2.2e-16 | 5.27e-03 |
| 10 | 1 | z | 0.891 +/- 0.209 | 0.9978 (expect 1) | 0.0e+00 | 5.27e-03 |
| 10 | 2 | diag111 | 0.891 +/- 0.209 | 1.9956 (expect 2) | 4.4e-16 | 5.27e-03 |
| 10 | 2 | oblique | 0.891 +/- 0.209 | 1.9956 (expect 2) | 2.2e-16 | 5.27e-03 |
| 10 | 2 | z | 0.891 +/- 0.209 | 1.9956 (expect 2) | 0.0e+00 | 5.27e-03 |

## Frame invariance (direct comparison)

For matched (kappa, theta, seed), the normalized radius recovered after rotating into an
arbitrary field direction is compared against the axis-aligned run draw by draw.

- comparisons: 60 run pairs, 6000000 draws
- bitwise identical radii: 1361180 / 6000000 (22.6863%)
- largest relative difference: 1.08e-15 (kappa=10, diag111, seed=1005)

Notes.

- Variance-based checks are reported only for kappa > 3/2, where the second moment
  exists. Below that the MAD ratio and the radial/quantile diagnostics carry the test.
- `max off-diag corr` is the largest off-diagonal correlation of the sample covariance
  in the recovered field-aligned frame; a misaligned rotation would inflate it. It is
  itself noisy at kappa = 2, where the fourth moment does not exist and the correlation
  estimator therefore has no finite variance.
