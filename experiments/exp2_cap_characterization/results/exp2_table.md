# Experiment 2 results -- the component-wise velocity cap

**The capped and uncapped modes sample two different probability laws.**
`no_cap()` samples the untruncated bi-Kappa distribution. A finite
`max_normalized_velocity = lambda` samples that distribution *conditioned on* the
component-wise box `|v_x|/theta_perp <= lambda AND |v_y|/theta_perp <= lambda AND
|v_z|/theta_par <= lambda`. The capped mode is a truncated/conditional target, not a
numerical approximation to the uncapped one, and must never be reported as the same
distribution.

5 replicate seeds (2001-2005) x 10^5 samples per configuration, double precision,
Apple clang / libc++, `std::mt19937`, `theta_perp:theta_par = 1:2` unless stated.

## How the rejected fraction was measured

`operator()` loops internally and reports no attempt count, and
`cpp/bi_kappa_distribution.H` was **not** modified. The shipped predicate
`withinNormalizedVelocityCap` was instead transcribed into `exp2_analyze.py:in_box`
and evaluated on the draws of the **uncapped** run at the same seed. Every loop
iteration of `operator()` consumes `x1, x2, cosTheta, phi` in the same order whether
or not a cap is in force, so the uncapped run *is* the capped run's attempt stream;
the mean of the predicate over it is `P(accept)` directly.

That correspondence is verified, not assumed: for all 240 (case, lambda, seed) pairs the capped run's
output is **bitwise identical** to the uncapped run's draws restricted to the box (240/240 pairs).

## 1. Headline -- rejected (redrawn) fraction, kappa x lambda

Fraction of attempts thrown away and redrawn, `1 - P(accept)`. Empirical values are
the mean over 5 seeds of the box-predicate rate on 10^5 uncapped attempts each;
`analytic` is the exact value from the closed-form derivation below.

**This same number is also the exact total-variation distance between the capped law
and the untruncated target**, because conditioning on an event of probability `p`
gives density ratio `1_box/p` and hence `TV = 1 - p`. So the rejection column is not
merely a cost: it is the distortion.

| kappa | lambda=3 | lambda=5 | lambda=10 | lambda=20 | lambda=50 | lambda=100 |
|---|---|---|---|---|---|---|
| 0.75 | 0.54881 (0.54834) | 0.43023 (0.43014) | 0.30585 (0.30583) | 0.21616 (0.21656) | 0.13643 (0.13702) | 0.09679 (0.09689) |
| 1 | 0.33549 (0.33525) | 0.20833 (0.20770) | 0.10507 (0.10533) | 0.05309 (0.05285) | 0.02108 (0.02116) | 0.01058 (0.01058) |
| 1.5 | 0.15199 (0.15175) | 0.05992 (0.05976) | 0.01555 (0.01555) | 0.00402 (0.00393) | 0.00064 (0.00063) | 0.00018 (0.00016) |
| 2 | 0.08124 (0.08079) | 0.02082 (0.02071) | 0.00289 (0.00280) | 0.00040 (0.00036) | 0.00002 (0.00002) | 0.00001 (2.88e-06) |
| 5 | 0.00854 (0.00857) | 0.00025 (0.00026) | 0.00000 (8.69e-07) | 0.00000 (1.97e-09) | 0.00000 (5.40e-13) | 0.00000 (1.06e-15) |
| 10 | 0.00161 (0.00167) | 0.00001 (4.26e-06) | 0.00000 (7.23e-11) | 0.00000 (2.61e-16) | 0.00000 (8.66e-24) | 0.00000 (1.70e-29) |

Format: empirical (analytic). The empirical column is a mean over 5x10^5 attempts,
so its own resolution is ~2e-6 and it reads 0.00000 wherever the analytic value is
below that; the analytic column is the one to quote in those cells.

### Mean attempts per accepted draw (analytic)

| kappa | lambda=3 | lambda=5 | lambda=10 | lambda=20 | lambda=50 | lambda=100 |
|---|---|---|---|---|---|---|
| 0.75 | 2.2140 | 1.7548 | 1.4406 | 1.2764 | 1.1588 | 1.1073 |
| 1 | 1.5043 | 1.2621 | 1.1177 | 1.0558 | 1.0216 | 1.0107 |
| 1.5 | 1.1789 | 1.0636 | 1.0158 | 1.0039 | 1.0006 | 1.0002 |
| 2 | 1.0879 | 1.0211 | 1.0028 | 1.0004 | 1.0000 | 1.0000 |
| 5 | 1.0086 | 1.0003 | 1.0000 | 1.0000 | 1.0000 | 1.0000 |
| 10 | 1.0017 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 |

### Retry-count distribution and the internal attempt limit

Attempts consumed per accepted draw, read off the gaps between accepted indices in
the uncapped attempt stream. If attempts are i.i.d. this is `Geometric(p)` on
`{1, 2, ...}`; the geometric reference is shown alongside. `operator()` throws after
`kMaxCapRejectTries = 1,000,000` consecutive rejections
(`bi_kappa_distribution.H:247`); `log10 P(throw)` per draw is `1e6 log10(1-p)`.

| kappa | lambda | mean attempts | (geometric) | p99 attempts | (geometric) | max observed | log10 P(hit try limit) |
|---|---|---|---|---|---|---|---|
| 0.75 | 3 | 2.2164 | 2.2164 | 8.0 | 8.0 | 22 | -2.61e+05 |
| 0.75 | 5 | 1.7551 | 1.7551 | 6.0 | 6.0 | 15 | -3.66e+05 |
| 0.75 | 10 | 1.4406 | 1.4406 | 4.0 | 4.0 | 11 | -5.15e+05 |
| 0.75 | 20 | 1.2758 | 1.2758 | 3.6 | 3.8 | 9 | -6.64e+05 |
| 0.75 | 50 | 1.1580 | 1.1580 | 3.0 | 3.0 | 6 | -8.63e+05 |
| 0.75 | 100 | 1.1072 | 1.1072 | 2.2 | 2.0 | 6 | -1.01e+06 |
| 1 | 3 | 1.5049 | 1.5049 | 5.0 | 5.0 | 11 | -4.75e+05 |
| 1 | 5 | 1.2632 | 1.2632 | 3.0 | 3.0 | 9 | -6.83e+05 |
| 1 | 10 | 1.1174 | 1.1174 | 3.0 | 3.0 | 6 | -9.77e+05 |
| 1 | 20 | 1.0561 | 1.0561 | 2.0 | 2.0 | 5 | -1.28e+06 |
| 1 | 50 | 1.0215 | 1.0215 | 2.0 | 2.0 | 4 | -1.67e+06 |
| 1 | 100 | 1.0107 | 1.0107 | 2.0 | 2.0 | 4 | -1.98e+06 |
| 1.5 | 3 | 1.1792 | 1.1792 | 3.0 | 3.0 | 8 | -8.19e+05 |
| 1.5 | 5 | 1.0637 | 1.0637 | 2.0 | 2.0 | 5 | -1.22e+06 |
| 1.5 | 10 | 1.0158 | 1.0158 | 2.0 | 2.0 | 4 | -1.81e+06 |
| 1.5 | 20 | 1.0040 | 1.0040 | 1.0 | 1.0 | 3 | -2.41e+06 |
| 1.5 | 50 | 1.0006 | 1.0006 | 1.0 | 1.0 | 2 | -3.2e+06 |
| 1.5 | 100 | 1.0002 | 1.0002 | 1.0 | 1.0 | 2 | -3.8e+06 |
| 2 | 3 | 1.0884 | 1.0884 | 2.0 | 2.0 | 6 | -1.09e+06 |
| 2 | 5 | 1.0213 | 1.0213 | 2.0 | 2.0 | 4 | -1.68e+06 |
| 2 | 10 | 1.0029 | 1.0029 | 1.0 | 1.0 | 3 | -2.55e+06 |
| 2 | 20 | 1.0004 | 1.0004 | 1.0 | 1.0 | 2 | -3.45e+06 |
| 2 | 50 | 1.0000 | 1.0000 | 1.0 | 1.0 | 2 | -4.64e+06 |
| 2 | 100 | 1.0000 | 1.0000 | 1.0 | 1.0 | 2 | -5.54e+06 |
| 5 | 3 | 1.0086 | 1.0086 | 1.0 | 1.0 | 3 | -2.07e+06 |
| 5 | 5 | 1.0003 | 1.0003 | 1.0 | 1.0 | 2 | -3.59e+06 |
| 5 | 10 | 1.0000 | 1.0000 | 1.0 | 1.0 | 1 | -6.06e+06 |
| 5 | 20 | 1.0000 | 1.0000 | 1.0 | 1.0 | 1 | -8.7e+06 |
| 5 | 50 | 1.0000 | 1.0000 | 1.0 | 1.0 | 1 | -1.23e+07 |
| 5 | 100 | 1.0000 | 1.0000 | 1.0 | 1.0 | 1 | -1.5e+07 |
| 10 | 3 | 1.0016 | 1.0016 | 1.0 | 1.0 | 3 | -2.78e+06 |
| 10 | 5 | 1.0000 | 1.0000 | 1.0 | 1.0 | 2 | -5.37e+06 |
| 10 | 10 | 1.0000 | 1.0000 | 1.0 | 1.0 | 1 | -1.01e+07 |
| 10 | 20 | 1.0000 | 1.0000 | 1.0 | 1.0 | 1 | -1.56e+07 |
| 10 | 50 | 1.0000 | 1.0000 | 1.0 | 1.0 | 1 | -2.31e+07 |
| 10 | 100 | 1.0000 | 1.0000 | 1.0 | 1.0 | 1 | -2.88e+07 |

The retry counts follow the geometric reference throughout, and the internal
attempt limit is unreachable: even at the worst configuration in the sweep the
log-probability of a single draw exhausting it is about -260 000. The attempt limit
is therefore not a practical failure mode, and the cost of the cap is entirely the
mean-attempt overhead plus -- far more importantly -- the change of target law.


### Closed form for the acceptance probability

In the field-aligned local frame the sampler emits
`v = sqrt(kappa) (theta_perp u_x, theta_perp u_y, theta_par u_z)`, so the shipped
predicate reduces to `sqrt(kappa) max_i |u_i| <= lambda`: **the theta's cancel**, and
the box event is a cube of half-side `c = lambda/sqrt(kappa)` in the isotropic
u-coordinates. Hence, with `u = R n`, `n` uniform on S^2, `M = max_i |n_i|`,

```
P(accept) = E_M[ I_z(3/2, kappa-1/2) ],   z = c^2 / (M^2 + c^2),
f_M(m) = 3 - (12/pi) arcsin( sqrt((1-2m^2)/(1-m^2)) )   for 1/sqrt3 <= m <= 1/sqrt2,
f_M(m) = 3                                              for 1/sqrt2 <= m <= 1.
```

`f_M` integrates to 0.999999999800 and its CDF agrees with a 2e+06-point spherical Monte Carlo to 3.36e-04.

Three consequences:

1. `P(accept)` depends on `(kappa, lambda)` only -- **not** on `theta_perp`,
   `theta_par` or the field direction (the cap is tested before the frame rotation).
2. The rejected fraction decays only *algebraically* in lambda, like
   `lambda^{-(2 kappa - 1)}`, so at heavy tails it is stubborn: see the table above.
3. The cap is a **cube** in the isotropic coordinates, so the capped law is not
   isotropic in direction and is not axisymmetric about **B**; see section 4.

### Theta-independence of acceptance -- isotropic (1,1) control vs anisotropic (1,2)

| kappa | lambda | reject frac (1,1) | reject frac (1,2) | abs diff | analytic |
|---|---|---|---|---|---|
| 0.75 | 3 | 0.54881 | 0.54881 | 0.00e+00 | 0.54834 |
| 0.75 | 5 | 0.43023 | 0.43023 | 0.00e+00 | 0.43014 |
| 0.75 | 10 | 0.30585 | 0.30585 | 0.00e+00 | 0.30583 |
| 0.75 | 20 | 0.21616 | 0.21616 | 0.00e+00 | 0.21656 |
| 0.75 | 50 | 0.13643 | 0.13643 | 0.00e+00 | 0.13702 |
| 0.75 | 100 | 0.09679 | 0.09679 | 0.00e+00 | 0.09689 |
| 2 | 3 | 0.08124 | 0.08124 | 0.00e+00 | 0.08079 |
| 2 | 5 | 0.02082 | 0.02082 | 0.00e+00 | 0.02071 |
| 2 | 10 | 0.00289 | 0.00289 | 0.00e+00 | 0.00280 |
| 2 | 20 | 0.00040 | 0.00040 | 0.00e+00 | 0.00036 |
| 2 | 50 | 0.00002 | 0.00002 | 0.00e+00 | 0.00002 |
| 2 | 100 | 0.00001 | 0.00001 | 0.00e+00 | 0.00000 |

These are **not** two independent estimates that happen to agree. The two runs
share the seed and therefore the RNG stream, so the prediction is the stronger one
that the accept/reject decision agrees on every individual attempt -- and it does:
0 disagreements in 6000000 attempts
across 60 (kappa, seed, lambda) comparisons. The cancellation of
`theta_perp` and `theta_par` in the cap predicate is exact, not statistical.

### Tail-exponent scaling of the rejected fraction

`P(R > r) ~ r^{-(2 kappa - 1)}` for the untruncated radial law, so the rejected
fraction must decay as `lambda^{-(2 kappa - 1)}`. Local log-log slope of the exact
rejection curve between the two widest caps:

| kappa | lambda pair | measured slope | predicted -(2 kappa - 1) |
|---|---|---|---|
| 0.75 | 50 -> 100 | -0.4999 | -0.5000 |
| 1 | 50 -> 100 | -0.9998 | -1.0000 |
| 1.5 | 50 -> 100 | -1.9994 | -2.0000 |
| 2 | 50 -> 100 | -2.9988 | -3.0000 |
| 5 | 50 -> 100 | -8.9911 | -9.0000 |
| 10 | 50 -> 100 | -18.9609 | -19.0000 |

This is the whole problem in one line. The cost of the cap is not exponentially
small in lambda, it is a power law, and the power is weakest exactly where kappa
distributions are physically interesting.

## 2. Distortion of the capped law relative to the uncapped target

Robust quantile ratios `q_capped / q_uncapped` at matched seeds, the largest
empirical-CDF gap on the speed `|v|`, and the exact total-variation distance
(analytic, since the empirical estimate saturates at ~2e-6). A ratio of 1 means no
distortion.

**Read the last two columns against the p99.9 column, not instead of it.** `TV` is an
upper bound on how much any *probability* can move, and the sup ECDF gap is a lower
bound on it -- but neither bounds a *quantile* ratio. A cap can amputate the entire
far tail while moving no probability by more than 1e-3, because the amputated region
carries almost no probability and enormous velocity. At `kappa = 1.5, lambda = 50`
the total-variation distance is 6.3e-04, and the p99.9 speed is still
20% too small.

| kappa | lambda | |v| p50 | |v| p90 | |v| p99 | |v| p99.9 | |v_z| p99.9 | sup ECDF gap on |v| | TV (analytic) |
|---|---|---|---|---|---|---|---|---|
| 0.75 | 3 | 0.3518 | 0.0292 | 4.02e-04 | 4.63e-06 | 6.03e-06 | 0.48706 | 5.483e-01 |
| 0.75 | 5 | 0.4396 | 0.0429 | 6.53e-04 | 7.62e-06 | 1.00e-05 | 0.38677 | 4.301e-01 |
| 0.75 | 10 | 0.5523 | 0.0714 | 0.0013 | 1.50e-05 | 1.99e-05 | 0.28228 | 3.058e-01 |
| 0.75 | 20 | 0.6519 | 0.1156 | 0.0024 | 2.95e-05 | 3.96e-05 | 0.20409 | 2.166e-01 |
| 0.75 | 50 | 0.7598 | 0.2027 | 0.0054 | 7.23e-05 | 9.76e-05 | 0.13141 | 1.370e-01 |
| 0.75 | 100 | 0.8200 | 0.2872 | 0.0096 | 1.41e-04 | 1.91e-04 | 0.09427 | 9.689e-02 |
| 1 | 3 | 0.6533 | 0.2445 | 0.0342 | 0.0039 | 0.0047 | 0.28360 | 3.353e-01 |
| 1 | 5 | 0.7666 | 0.3410 | 0.0546 | 0.0064 | 0.0078 | 0.18527 | 2.077e-01 |
| 1 | 10 | 0.8726 | 0.5101 | 0.1000 | 0.0124 | 0.0153 | 0.09839 | 1.053e-01 |
| 1 | 20 | 0.9337 | 0.6760 | 0.1721 | 0.0240 | 0.0298 | 0.05147 | 5.285e-02 |
| 1 | 50 | 0.9735 | 0.8394 | 0.3242 | 0.0560 | 0.0690 | 0.02081 | 2.116e-02 |
| 1 | 100 | 0.9867 | 0.9123 | 0.4871 | 0.1026 | 0.1257 | 0.01051 | 1.058e-02 |
| 1.5 | 3 | 0.8658 | 0.6183 | 0.2759 | 0.0995 | 0.1097 | 0.12703 | 1.518e-01 |
| 1.5 | 5 | 0.9449 | 0.7813 | 0.4160 | 0.1602 | 0.1800 | 0.05486 | 5.976e-02 |
| 1.5 | 10 | 0.9858 | 0.9291 | 0.6405 | 0.3004 | 0.3367 | 0.01513 | 1.555e-02 |
| 1.5 | 20 | 0.9963 | 0.9801 | 0.8430 | 0.5099 | 0.5727 | 0.00400 | 3.929e-03 |
| 1.5 | 50 | 0.9994 | 0.9970 | 0.9699 | 0.7987 | 0.8331 | 0.00064 | 6.304e-04 |
| 1.5 | 100 | 0.9999 | 0.9992 | 0.9926 | 0.9160 | 0.9295 | 0.00018 | 1.577e-04 |
| 2 | 3 | 0.9336 | 0.7904 | 0.5110 | 0.2645 | 0.2794 | 0.06881 | 8.079e-02 |
| 2 | 5 | 0.9827 | 0.9272 | 0.7133 | 0.4157 | 0.4500 | 0.01957 | 2.071e-02 |
| 2 | 10 | 0.9976 | 0.9885 | 0.9122 | 0.7046 | 0.7524 | 0.00288 | 2.800e-03 |
| 2 | 20 | 0.9996 | 0.9986 | 0.9841 | 0.8882 | 0.9157 | 0.00040 | 3.571e-04 |
| 2 | 50 | 1.0000 | 0.9999 | 0.9989 | 0.9917 | 0.9884 | 0.00002 | 2.299e-05 |
| 2 | 100 | 1.0000 | 1.0000 | 0.9998 | 0.9985 | 0.9987 | 0.00001 | 2.876e-06 |
| 5 | 3 | 0.9940 | 0.9792 | 0.9364 | 0.7953 | 0.8052 | 0.00776 | 8.566e-03 |
| 5 | 5 | 0.9999 | 0.9995 | 0.9956 | 0.9832 | 0.9839 | 0.00025 | 2.556e-04 |
| 5 | 10 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 0.00000 | 8.689e-07 |
| 5 | 20 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 0.00000 | 1.973e-09 |
| 5 | 50 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 0.00000 | 5.399e-13 |
| 5 | 100 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 0.00000 | 1.061e-15 |
| 10 | 3 | 0.9990 | 0.9965 | 0.9882 | 0.9510 | 0.9510 | 0.00152 | 1.673e-03 |
| 10 | 5 | 1.0000 | 1.0000 | 0.9999 | 0.9998 | 0.9998 | 0.00001 | 4.262e-06 |
| 10 | 10 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 0.00000 | 7.229e-11 |
| 10 | 20 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 0.00000 | 2.613e-16 |
| 10 | 50 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 0.00000 | 8.663e-24 |
| 10 | 100 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 0.00000 | 1.698e-29 |

### Absolute scale of the quantiles being compared (uncapped target, |v|)

Given so the ratios above can be read: at `kappa = 0.75` the untruncated p99.9 speed
is 1.47e+06 thermal units, while a capped draw obeys `|v| <= lambda sqrt(theta_perp^2 + theta_perp^2 + theta_par^2)` by construction. The
ratio is ~1e-5 and the tail is simply gone -- not compressed, gone.

| kappa | p50 | p90 | p99 | p99.9 |
|---|---|---|---|---|
| 0.75 | 6.003 | 153 | 1.521e+04 | 1.47e+06 |
| 1 | 3.075 | 17.43 | 176.1 | 1723 |
| 1.5 | 2.161 | 6.366 | 21.25 | 65.25 |
| 2 | 1.906 | 4.659 | 11.14 | 24.22 |
| 5 | 1.592 | 3.121 | 5.176 | 7.475 |
| 10 | 1.512 | 2.82 | 4.369 | 5.786 |

### Where the cap becomes negligible

Criterion: `TV < 0.001` **and** p99.9 speed quantile
distorted by `< 1%`. Both are required precisely
because they fail in different places.

| block | kappa | (theta_perp, theta_par) | smallest lambda in the ladder that qualifies |
|---|---|---|---|
| A | 0.75 | (1, 2) | **none** of 3, 5, 10, 20, 50, 100 |
| A | 1 | (1, 2) | **none** of 3, 5, 10, 20, 50, 100 |
| A | 1.5 | (1, 2) | **none** of 3, 5, 10, 20, 50, 100 |
| A | 2 | (1, 2) | 50 |
| A | 5 | (1, 2) | 10 |
| A | 10 | (1, 2) | 5 |
| C | 0.75 | (1, 1) | **none** of 3, 5, 10, 20, 50, 100 |
| C | 2 | (1, 1) | 50 |

## 3. Moments -- what could and could not legitimately be compared

The untruncated bi-Kappa second moment is `theta^2 kappa/(2 kappa - 3)`. It exists
**only for kappa > 3/2** and diverges at `kappa = 3/2`. For `kappa <= 3/2` there is
no untruncated reference value, so a capped-vs-uncapped variance comparison would be
comparing a number to a divergent integral. Those comparisons are **refused**, not
silently omitted.

| kappa | second moment of the untruncated target | comparison |
|---|---|---|
| 0.75 | **does not exist** | **REFUSED -- there is no reference value to compare against** |
| 1 | **does not exist** | **REFUSED -- there is no reference value to compare against** |
| 1.5 | **does not exist** | **REFUSED -- there is no reference value to compare against** |
| 2 | exists, = 2 theta^2 | reported below |
| 5 | exists, = 0.7143 theta^2 | reported below |
| 10 | exists, = 0.5882 theta^2 | reported below |

### Variance ratio to the untruncated theory, kappa > 3/2 only

`sample var / (theta^2 kappa/(2 kappa - 3))`, perpendicular then parallel component.

| kappa | mode | lambda | var_x / theory | var_z / theory |
|---|---|---|---|---|
| 2 | uncapped | - | 1.0432 | 0.9891 |
| 2 | capped | 3 | 0.4629 | 0.4607 |
| 2 | capped | 5 | 0.6388 | 0.6364 |
| 2 | capped | 10 | 0.8093 | 0.8074 |
| 2 | capped | 20 | 0.9071 | 0.8997 |
| 2 | capped | 50 | 0.9825 | 0.9732 |
| 2 | capped | 100 | 0.9991 | 0.9880 |
| 5 | uncapped | - | 1.0048 | 0.9999 |
| 5 | capped | 3 | 0.9452 | 0.9428 |
| 5 | capped | 5 | 0.9993 | 0.9949 |
| 5 | capped | 10 | 1.0048 | 0.9999 |
| 5 | capped | 20 | 1.0048 | 0.9999 |
| 5 | capped | 50 | 1.0048 | 0.9999 |
| 5 | capped | 100 | 1.0048 | 0.9999 |
| 10 | uncapped | - | 1.0010 | 0.9989 |
| 10 | capped | 3 | 0.9896 | 0.9882 |
| 10 | capped | 5 | 1.0008 | 0.9988 |
| 10 | capped | 10 | 1.0010 | 0.9989 |
| 10 | capped | 20 | 1.0010 | 0.9989 |
| 10 | capped | 50 | 1.0010 | 0.9989 |
| 10 | capped | 100 | 1.0010 | 0.9989 |

### The cap manufactures a finite variance where none exists

For `kappa <= 3/2` the capped sample has a perfectly finite variance -- it has
bounded support -- and that number is reported by any naive diagnostic. It is an
artifact of `lambda`, not a property of the plasma. Since the untruncated radial
density behaves as `f_R(r) ~ r^{-2 kappa}`, the second moment truncated at `~lambda`
grows as `lambda^{3 - 2 kappa}` for `kappa < 3/2` and as `log lambda` at
`kappa = 3/2` exactly -- without bound in both cases. Reporting it as *the* variance
of a kappa distribution with `kappa <= 3/2` would be wrong.

| kappa | var_x, lambda=3 | var_x, lambda=5 | var_x, lambda=10 | var_x, lambda=20 | var_x, lambda=50 | var_x, lambda=100 | measured growth (20 -> 50) | predicted |
|---|---|---|---|---|---|---|---|---|
| 0.75 | 1.274 | 2.611 | 6.78 | 17.83 | 65.32 | 175.7 | power-law slope 1.428 | `lambda^{1.5}` |
| 1 | 1.178 | 2.181 | 4.771 | 9.892 | 25.85 | 51.82 | power-law slope 1.003 | `lambda^{1}` |
| 1.5 | 1.027 | 1.616 | 2.557 | 3.543 | 4.939 | 5.953 | log-slope B = 1.46 (B over all four pairs: 1.15, 1.36, 1.42, 1.52, 1.46) | `A + B log lambda`, B constant |

The uncapped runs have no entry here on purpose: their sample variance is a finite
number produced by a divergent population moment and means nothing. The capped runs
do have a well-defined population variance -- it is just a variance of the box, not
of the plasma.

Two caveats on the fitted rates. At `kappa = 3/2` the log-slope `B` rises
from 1.15 to 1.47 across the ladder rather than sitting at a constant: the growth is
unmistakably slower than any power of lambda and consistent with `log lambda`, but
lambda = 50 is not yet deep enough in the asymptotic regime to pin `B` down. At
`kappa = 0.75` the measured power-law slope 1.41 likewise approaches the asymptotic
1.5 from below, because an O(1) contribution from the distribution core is still
present; fitting `var = A + B lambda^{1.5}` to the two widest caps gives A = 1.86,
B = 0.179, which reproduces the lambda = 10 point to 11%. Neither caveat touches the
conclusion, which is that the number diverges as the cap is relaxed and therefore
is not a property of the distribution being sampled.

## 4. Angular structure -- the cap breaks axisymmetry about B

The box is a **cube** in the isotropic u-coordinates, so a direction pointing at a
cube corner has `sqrt(3)` times more radial room than a direction along an axis. The
conditioned law therefore has a four-fold azimuthal modulation about **B** and a
non-uniform polar distribution -- structure the physical bi-Kappa distribution does
not have and which no choice of `theta` can absorb.

`a4 = 2<cos 4 phi>` is the leading azimuthal Fourier coefficient (0 for the target);
`z4 = a4 sqrt(n/2)` is it in units of its own sampling s.d., so |z4| > 3 is a
detection. `cos(theta) sqrt(n) D` is the KS statistic against U(-1,1), which is
O(1) (about 0.87) when the law is correct.

| kappa | lambda | a4 | z4 | cos(theta) sqrt(n)D |
|---|---|---|---|---|
| 0.75 | uncapped | -0.00117 | -0.3 | 0.712 |
| 0.75 | 3 | -0.05183 | -11.6 | 2.392 |
| 0.75 | 5 | -0.03539 | -7.9 | 1.887 |
| 0.75 | 10 | -0.02157 | -4.8 | 1.214 |
| 0.75 | 20 | -0.01284 | -2.9 | 0.948 |
| 0.75 | 50 | -0.00776 | -1.7 | 0.826 |
| 0.75 | 100 | -0.00469 | -1.0 | 0.769 |
| 1 | uncapped | -0.00222 | -0.5 | 0.863 |
| 1 | 3 | -0.04172 | -9.3 | 1.969 |
| 1 | 5 | -0.02338 | -5.2 | 1.275 |
| 1 | 10 | -0.01160 | -2.6 | 1.017 |
| 1 | 20 | -0.00572 | -1.3 | 0.931 |
| 1 | 50 | -0.00351 | -0.8 | 0.852 |
| 1 | 100 | -0.00306 | -0.7 | 0.834 |
| 1.5 | uncapped | -0.00095 | -0.2 | 0.656 |
| 1.5 | 3 | -0.02804 | -6.3 | 1.569 |
| 1.5 | 5 | -0.00983 | -2.2 | 0.910 |
| 1.5 | 10 | -0.00358 | -0.8 | 0.639 |
| 1.5 | 20 | -0.00171 | -0.4 | 0.659 |
| 1.5 | 50 | -0.00102 | -0.2 | 0.659 |
| 1.5 | 100 | -0.00095 | -0.2 | 0.655 |
| 2 | uncapped | -0.00130 | -0.3 | 0.768 |
| 2 | 3 | -0.01965 | -4.4 | 1.307 |
| 2 | 5 | -0.00608 | -1.4 | 0.786 |
| 2 | 10 | -0.00209 | -0.5 | 0.753 |
| 2 | 20 | -0.00136 | -0.3 | 0.758 |
| 2 | 50 | -0.00131 | -0.3 | 0.768 |
| 2 | 100 | -0.00130 | -0.3 | 0.768 |
| 5 | uncapped | -0.00036 | -0.1 | 0.750 |
| 5 | 3 | -0.00432 | -1.0 | 0.777 |
| 5 | 5 | -0.00052 | -0.1 | 0.745 |
| 5 | 10 | -0.00036 | -0.1 | 0.750 |
| 5 | 20 | -0.00036 | -0.1 | 0.750 |
| 5 | 50 | -0.00036 | -0.1 | 0.750 |
| 5 | 100 | -0.00036 | -0.1 | 0.750 |
| 10 | uncapped | +0.00056 | +0.1 | 0.701 |
| 10 | 3 | -0.00053 | -0.1 | 0.702 |
| 10 | 5 | +0.00056 | +0.1 | 0.700 |
| 10 | 10 | +0.00056 | +0.1 | 0.701 |
| 10 | 20 | +0.00056 | +0.1 | 0.701 |
| 10 | 50 | +0.00056 | +0.1 | 0.701 |
| 10 | 100 | +0.00056 | +0.1 | 0.701 |

## 5. Summary

1. `no_cap()` samples the bi-Kappa distribution itself; validation runs use it.
2. A finite `max_normalized_velocity` samples the bi-Kappa distribution conditioned on
   a component-wise box. It is a cube in normalized velocity components: it is not
   isotropic, it is not axisymmetric about **B** (section 4), and its shape depends
   on lambda.
3. The rejected fraction equals the total-variation distance from the target,
   `1 - P(accept)`, with the closed form above.
4. A variance from a capped run at `kappa <= 3/2` is a variance of the box, not of
   the bi-Kappa distribution.
