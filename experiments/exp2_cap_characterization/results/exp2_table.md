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
| 0.75 | 0.54810 (0.54834) | 0.43076 (0.43014) | 0.30662 (0.30583) | 0.21669 (0.21656) | 0.13732 (0.13702) | 0.09689 (0.09689) |
| 1 | 0.33619 (0.33525) | 0.20823 (0.20770) | 0.10514 (0.10533) | 0.05261 (0.05285) | 0.02101 (0.02116) | 0.01055 (0.01058) |
| 1.5 | 0.15184 (0.15175) | 0.05919 (0.05976) | 0.01522 (0.01555) | 0.00387 (0.00393) | 0.00066 (0.00063) | 0.00015 (0.00016) |
| 2 | 0.08081 (0.08079) | 0.02073 (0.02071) | 0.00275 (0.00280) | 0.00036 (0.00036) | 0.00002 (0.00002) | 0.00000 (2.88e-06) |
| 5 | 0.00854 (0.00857) | 0.00021 (0.00026) | 0.00000 (8.69e-07) | 0.00000 (1.97e-09) | 0.00000 (5.40e-13) | 0.00000 (1.06e-15) |
| 10 | 0.00167 (0.00167) | 0.00000 (4.26e-06) | 0.00000 (7.23e-11) | 0.00000 (2.61e-16) | 0.00000 (8.66e-24) | 0.00000 (1.70e-29) |

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
| 0.75 | 3 | 2.2128 | 2.2129 | 8.0 | 8.0 | 22 | -2.61e+05 |
| 0.75 | 5 | 1.7567 | 1.7567 | 6.0 | 6.0 | 16 | -3.66e+05 |
| 0.75 | 10 | 1.4422 | 1.4422 | 4.0 | 4.0 | 11 | -5.15e+05 |
| 0.75 | 20 | 1.2766 | 1.2766 | 3.6 | 3.8 | 8 | -6.64e+05 |
| 0.75 | 50 | 1.1592 | 1.1592 | 3.0 | 3.0 | 7 | -8.63e+05 |
| 0.75 | 100 | 1.1073 | 1.1073 | 2.0 | 2.0 | 6 | -1.01e+06 |
| 1 | 3 | 1.5065 | 1.5065 | 5.0 | 5.0 | 12 | -4.75e+05 |
| 1 | 5 | 1.2630 | 1.2630 | 3.0 | 3.0 | 9 | -6.83e+05 |
| 1 | 10 | 1.1175 | 1.1175 | 3.0 | 3.0 | 6 | -9.77e+05 |
| 1 | 20 | 1.0555 | 1.0555 | 2.0 | 2.0 | 5 | -1.28e+06 |
| 1 | 50 | 1.0215 | 1.0215 | 2.0 | 2.0 | 5 | -1.67e+06 |
| 1 | 100 | 1.0107 | 1.0107 | 2.0 | 2.0 | 4 | -1.98e+06 |
| 1.5 | 3 | 1.1790 | 1.1790 | 3.0 | 3.0 | 8 | -8.19e+05 |
| 1.5 | 5 | 1.0629 | 1.0629 | 2.0 | 2.0 | 5 | -1.22e+06 |
| 1.5 | 10 | 1.0155 | 1.0155 | 2.0 | 2.0 | 4 | -1.81e+06 |
| 1.5 | 20 | 1.0039 | 1.0039 | 1.0 | 1.0 | 3 | -2.41e+06 |
| 1.5 | 50 | 1.0007 | 1.0007 | 1.0 | 1.0 | 2 | -3.2e+06 |
| 1.5 | 100 | 1.0001 | 1.0001 | 1.0 | 1.0 | 2 | -3.8e+06 |
| 2 | 3 | 1.0879 | 1.0879 | 2.0 | 2.0 | 6 | -1.09e+06 |
| 2 | 5 | 1.0212 | 1.0212 | 2.0 | 2.0 | 4 | -1.68e+06 |
| 2 | 10 | 1.0028 | 1.0028 | 1.0 | 1.0 | 3 | -2.55e+06 |
| 2 | 20 | 1.0004 | 1.0004 | 1.0 | 1.0 | 2 | -3.45e+06 |
| 2 | 50 | 1.0000 | 1.0000 | 1.0 | 1.0 | 2 | -4.64e+06 |
| 2 | 100 | 1.0000 | 1.0000 | 1.0 | 1.0 | 2 | -5.54e+06 |
| 5 | 3 | 1.0086 | 1.0086 | 1.0 | 1.0 | 4 | -2.07e+06 |
| 5 | 5 | 1.0002 | 1.0002 | 1.0 | 1.0 | 2 | -3.59e+06 |
| 5 | 10 | 1.0000 | 1.0000 | 1.0 | 1.0 | 1 | -6.06e+06 |
| 5 | 20 | 1.0000 | 1.0000 | 1.0 | 1.0 | 1 | -8.7e+06 |
| 5 | 50 | 1.0000 | 1.0000 | 1.0 | 1.0 | 1 | -1.23e+07 |
| 5 | 100 | 1.0000 | 1.0000 | 1.0 | 1.0 | 1 | -1.5e+07 |
| 10 | 3 | 1.0017 | 1.0017 | 1.0 | 1.0 | 2 | -2.78e+06 |
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

Three consequences worth stating in the manuscript:

1. `P(accept)` depends on `(kappa, lambda)` only -- **not** on `theta_perp`,
   `theta_par` or the field direction (the cap is tested before the frame rotation).
2. The rejected fraction decays only *algebraically* in lambda, like
   `lambda^{-(2 kappa - 1)}`, so at heavy tails it is stubborn: see the table above.
3. The cap is a **cube** in the isotropic coordinates, so the capped law is not
   isotropic in direction and is not axisymmetric about **B**; see section 4.

### Theta-independence of acceptance -- isotropic (1,1) control vs anisotropic (1,2)

| kappa | lambda | reject frac (1,1) | reject frac (1,2) | abs diff | analytic |
|---|---|---|---|---|---|
| 0.75 | 3 | 0.54810 | 0.54810 | 0.00e+00 | 0.54834 |
| 0.75 | 5 | 0.43076 | 0.43076 | 0.00e+00 | 0.43014 |
| 0.75 | 10 | 0.30662 | 0.30662 | 0.00e+00 | 0.30583 |
| 0.75 | 20 | 0.21669 | 0.21669 | 0.00e+00 | 0.21656 |
| 0.75 | 50 | 0.13732 | 0.13732 | 0.00e+00 | 0.13702 |
| 0.75 | 100 | 0.09689 | 0.09689 | 0.00e+00 | 0.09689 |
| 2 | 3 | 0.08081 | 0.08081 | 0.00e+00 | 0.08079 |
| 2 | 5 | 0.02073 | 0.02073 | 0.00e+00 | 0.02071 |
| 2 | 10 | 0.00275 | 0.00275 | 0.00e+00 | 0.00280 |
| 2 | 20 | 0.00036 | 0.00036 | 0.00e+00 | 0.00036 |
| 2 | 50 | 0.00002 | 0.00002 | 0.00e+00 | 0.00002 |
| 2 | 100 | 0.00000 | 0.00000 | 0.00e+00 | 0.00000 |

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
the total-variation distance is 6.6e-4 -- by any probability-based measure the two
laws are indistinguishable -- and the p99.9 speed is still 24% too small.

| kappa | lambda | |v| p50 | |v| p90 | |v| p99 | |v| p99.9 | |v_z| p99.9 | sup ECDF gap on |v| | TV (analytic) |
|---|---|---|---|---|---|---|---|---|
| 0.75 | 3 | 0.3516 | 0.0288 | 3.94e-04 | 3.90e-06 | 5.60e-06 | 0.48699 | 5.483e-01 |
| 0.75 | 5 | 0.4385 | 0.0423 | 6.41e-04 | 6.44e-06 | 9.30e-06 | 0.38777 | 4.301e-01 |
| 0.75 | 10 | 0.5502 | 0.0707 | 0.0012 | 1.27e-05 | 1.85e-05 | 0.28278 | 3.058e-01 |
| 0.75 | 20 | 0.6493 | 0.1145 | 0.0023 | 2.48e-05 | 3.67e-05 | 0.20454 | 2.166e-01 |
| 0.75 | 50 | 0.7566 | 0.1997 | 0.0052 | 6.05e-05 | 9.02e-05 | 0.13231 | 1.370e-01 |
| 0.75 | 100 | 0.8195 | 0.2841 | 0.0096 | 1.19e-04 | 1.77e-04 | 0.09434 | 9.689e-02 |
| 1 | 3 | 0.6523 | 0.2437 | 0.0344 | 0.0036 | 0.0047 | 0.28452 | 3.353e-01 |
| 1 | 5 | 0.7660 | 0.3410 | 0.0549 | 0.0060 | 0.0077 | 0.18509 | 2.077e-01 |
| 1 | 10 | 0.8736 | 0.5069 | 0.1006 | 0.0116 | 0.0152 | 0.09842 | 1.053e-01 |
| 1 | 20 | 0.9344 | 0.6740 | 0.1724 | 0.0224 | 0.0297 | 0.05085 | 5.285e-02 |
| 1 | 50 | 0.9729 | 0.8389 | 0.3241 | 0.0521 | 0.0684 | 0.02076 | 2.116e-02 |
| 1 | 100 | 0.9862 | 0.9126 | 0.4870 | 0.0961 | 0.1246 | 0.01047 | 1.058e-02 |
| 1.5 | 3 | 0.8662 | 0.6191 | 0.2814 | 0.0947 | 0.1095 | 0.12696 | 1.518e-01 |
| 1.5 | 5 | 0.9457 | 0.7854 | 0.4258 | 0.1527 | 0.1800 | 0.05423 | 5.976e-02 |
| 1.5 | 10 | 0.9856 | 0.9321 | 0.6465 | 0.2858 | 0.3386 | 0.01483 | 1.555e-02 |
| 1.5 | 20 | 0.9963 | 0.9803 | 0.8498 | 0.4814 | 0.5588 | 0.00385 | 3.929e-03 |
| 1.5 | 50 | 0.9994 | 0.9967 | 0.9700 | 0.7611 | 0.8383 | 0.00065 | 6.304e-04 |
| 1.5 | 100 | 0.9999 | 0.9993 | 0.9927 | 0.9251 | 0.9480 | 0.00015 | 1.577e-04 |
| 2 | 3 | 0.9339 | 0.7908 | 0.5133 | 0.2600 | 0.2820 | 0.06802 | 8.079e-02 |
| 2 | 5 | 0.9829 | 0.9276 | 0.7199 | 0.4118 | 0.4556 | 0.01946 | 2.071e-02 |
| 2 | 10 | 0.9978 | 0.9890 | 0.9186 | 0.6963 | 0.7556 | 0.00273 | 2.800e-03 |
| 2 | 20 | 0.9998 | 0.9986 | 0.9884 | 0.8868 | 0.9323 | 0.00036 | 3.571e-04 |
| 2 | 50 | 1.0000 | 1.0000 | 0.9996 | 0.9945 | 0.9940 | 0.00002 | 2.299e-05 |
| 2 | 100 | 1.0000 | 1.0000 | 0.9999 | 0.9977 | 0.9988 | 0.00000 | 2.876e-06 |
| 5 | 3 | 0.9943 | 0.9792 | 0.9341 | 0.7900 | 0.8026 | 0.00770 | 8.566e-03 |
| 5 | 5 | 0.9998 | 0.9994 | 0.9962 | 0.9798 | 0.9870 | 0.00021 | 2.556e-04 |
| 5 | 10 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 0.00000 | 8.689e-07 |
| 5 | 20 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 0.00000 | 1.973e-09 |
| 5 | 50 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 0.00000 | 5.399e-13 |
| 5 | 100 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 0.00000 | 1.061e-15 |
| 10 | 3 | 0.9989 | 0.9961 | 0.9887 | 0.9452 | 0.9463 | 0.00154 | 1.673e-03 |
| 10 | 5 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 0.00000 | 4.262e-06 |
| 10 | 10 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 0.00000 | 7.229e-11 |
| 10 | 20 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 0.00000 | 2.613e-16 |
| 10 | 50 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 0.00000 | 8.663e-24 |
| 10 | 100 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 0.00000 | 1.698e-29 |

### Absolute scale of the quantiles being compared (uncapped target, |v|)

Given so the ratios above can be read: at `kappa = 0.75` the untruncated p99.9 speed
is 1.76e+06 thermal units, while a capped draw obeys `|v| <= lambda sqrt(theta_perp^2 + theta_perp^2 + theta_par^2)` by construction. The
ratio is ~1e-5 and the tail is simply gone -- not compressed, gone.

| kappa | p50 | p90 | p99 | p99.9 |
|---|---|---|---|---|
| 0.75 | 6.009 | 154.8 | 1.553e+04 | 1.765e+06 |
| 1 | 3.08 | 17.49 | 175.2 | 1826 |
| 1.5 | 2.158 | 6.35 | 20.85 | 68.92 |
| 2 | 1.903 | 4.663 | 11.09 | 24.51 |
| 5 | 1.59 | 3.119 | 5.179 | 7.512 |
| 10 | 1.513 | 2.821 | 4.368 | 5.831 |

### Where the cap becomes negligible

Criterion fixed in advance: `TV < 0.001` **and** p99.9 speed quantile
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
| 2 | uncapped | - | 1.0355 | 0.9792 |
| 2 | capped | 3 | 0.4590 | 0.4621 |
| 2 | capped | 5 | 0.6333 | 0.6388 |
| 2 | capped | 10 | 0.8075 | 0.8116 |
| 2 | capped | 20 | 0.9039 | 0.9049 |
| 2 | capped | 50 | 0.9589 | 0.9661 |
| 2 | capped | 100 | 0.9899 | 0.9783 |
| 5 | uncapped | - | 0.9975 | 0.9988 |
| 5 | capped | 3 | 0.9404 | 0.9428 |
| 5 | capped | 5 | 0.9935 | 0.9952 |
| 5 | capped | 10 | 0.9975 | 0.9988 |
| 5 | capped | 20 | 0.9975 | 0.9988 |
| 5 | capped | 50 | 0.9975 | 0.9988 |
| 5 | capped | 100 | 0.9975 | 0.9988 |
| 10 | uncapped | - | 0.9990 | 1.0011 |
| 10 | capped | 3 | 0.9875 | 0.9899 |
| 10 | capped | 5 | 0.9989 | 1.0011 |
| 10 | capped | 10 | 0.9990 | 1.0011 |
| 10 | capped | 20 | 0.9990 | 1.0011 |
| 10 | capped | 50 | 0.9990 | 1.0011 |
| 10 | capped | 100 | 0.9990 | 1.0011 |

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
| 0.75 | 1.277 | 2.607 | 6.77 | 17.88 | 65.18 | 178.3 | power-law slope 1.452 | `lambda^{1.5}` |
| 1 | 1.177 | 2.193 | 4.789 | 9.979 | 25.6 | 51.21 | power-law slope 1.000 | `lambda^{1}` |
| 1.5 | 1.026 | 1.615 | 2.541 | 3.521 | 4.869 | 6.005 | log-slope B = 1.64 (B over all four pairs: 1.15, 1.34, 1.41, 1.47, 1.64) | `A + B log lambda`, B constant |

The uncapped runs have no entry here on purpose: their sample variance is a finite
number produced by a divergent population moment and means nothing. The capped runs
do have a well-defined population variance -- it is just a variance of the box, not
of the plasma.

Two honest caveats on the fitted rates. At `kappa = 3/2` the log-slope `B` rises
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
| 0.75 | uncapped | -0.00317 | -0.7 | 0.843 |
| 0.75 | 3 | -0.05170 | -11.6 | 2.387 |
| 0.75 | 5 | -0.03547 | -7.9 | 1.767 |
| 0.75 | 10 | -0.02178 | -4.9 | 1.334 |
| 0.75 | 20 | -0.01430 | -3.2 | 1.055 |
| 0.75 | 50 | -0.01011 | -2.3 | 0.953 |
| 0.75 | 100 | -0.00729 | -1.6 | 0.913 |
| 1 | uncapped | -0.00348 | -0.8 | 0.817 |
| 1 | 3 | -0.04186 | -9.4 | 2.198 |
| 1 | 5 | -0.02479 | -5.5 | 1.402 |
| 1 | 10 | -0.01298 | -2.9 | 1.105 |
| 1 | 20 | -0.00803 | -1.8 | 0.999 |
| 1 | 50 | -0.00507 | -1.1 | 0.896 |
| 1 | 100 | -0.00373 | -0.8 | 0.863 |
| 1.5 | uncapped | -0.00008 | -0.0 | 0.962 |
| 1.5 | 3 | -0.02851 | -6.4 | 1.684 |
| 1.5 | 5 | -0.01089 | -2.4 | 1.149 |
| 1.5 | 10 | -0.00199 | -0.4 | 1.020 |
| 1.5 | 20 | -0.00054 | -0.1 | 0.971 |
| 1.5 | 50 | -0.00017 | -0.0 | 0.968 |
| 1.5 | 100 | -0.00012 | -0.0 | 0.964 |
| 2 | uncapped | +0.00159 | +0.4 | 0.811 |
| 2 | 3 | -0.01695 | -3.8 | 1.376 |
| 2 | 5 | -0.00372 | -0.8 | 0.934 |
| 2 | 10 | +0.00073 | +0.2 | 0.819 |
| 2 | 20 | +0.00142 | +0.3 | 0.812 |
| 2 | 50 | +0.00158 | +0.4 | 0.810 |
| 2 | 100 | +0.00159 | +0.4 | 0.811 |
| 5 | uncapped | -0.00137 | -0.3 | 0.744 |
| 5 | 3 | -0.00578 | -1.3 | 0.796 |
| 5 | 5 | -0.00143 | -0.3 | 0.746 |
| 5 | 10 | -0.00137 | -0.3 | 0.744 |
| 5 | 20 | -0.00137 | -0.3 | 0.744 |
| 5 | 50 | -0.00137 | -0.3 | 0.744 |
| 5 | 100 | -0.00137 | -0.3 | 0.744 |
| 10 | uncapped | -0.00157 | -0.4 | 0.722 |
| 10 | 3 | -0.00242 | -0.5 | 0.735 |
| 10 | 5 | -0.00157 | -0.4 | 0.722 |
| 10 | 10 | -0.00157 | -0.4 | 0.722 |
| 10 | 20 | -0.00157 | -0.4 | 0.722 |
| 10 | 50 | -0.00157 | -0.4 | 0.722 |
| 10 | 100 | -0.00157 | -0.4 | 0.722 |

## 5. Recommendation for the manuscript

1. Run every validation and every physics result with the cap **off**
   (`no_cap()`). That is the mode whose target law is the bi-Kappa distribution.
2. Document the finite `max_normalized_velocity` as an **optional pragmatic
   finite-velocity-box conditional target**: the bi-Kappa distribution conditioned on
   a component-wise box, useful only when a code needs bounded particle speeds.
3. Do **not** present it as a regularized or physically motivated kappa model. It is
   a cube in normalized velocity components: it is not isotropic, it is not
   axisymmetric about **B** (section 4), and its shape depends on the arbitrary
   choice of lambda.
4. State the cost and the distortion together, since they are the same number:
   rejected fraction = total-variation distance from the target = `1 - P(accept)`,
   with the closed form above.
5. Never quote a variance for `kappa <= 3/2` obtained from a capped run.
