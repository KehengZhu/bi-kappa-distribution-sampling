# Experiment 4 — supported numerical range of the released bi-Kappa sampler

All fractions are pooled over 5 seeds × 1,000,000 draws per configuration, uncapped.

## Q1 + Q2 — failure envelope and the effect of `sqrt(x1)/sqrt(x2)`

`ratio` is the pre-fix formation `sqrt(x1/x2)`; `split` is the current formation `sqrt(x1)/sqrt(x2)`. Both columns are evaluated on **identical variates**, so the difference is the fix and nothing else.

### libc++, float

| κ | shape(x2) | lost, `ratio` | lost, `split` | recovered by fix | honest overflow | still spurious |
|---|---|---|---|---|---|---|
| 0.5001 | 0.0001 | 0.99117 | 0.98966 | 0.0015 | 0.98242 | 0.00725 |
| 0.501 | 0.001 | 0.91559 | 0.90197 | 0.01363 | 0.8379 | 0.06407 |
| 0.505 | 0.005 | 0.64368 | 0.59706 | 0.04663 | 0.41316 | 0.1839 |
| 0.51 | 0.01 | 0.41451 | 0.35604 | 0.05847 | 0.17066 | 0.18538 |
| 0.55 | 0.05 | 0.01223 | 0.0057 | 0.00653 | 0.00014 | 0.00556 |
| 0.6 | 0.1 | 0.00015 | 2.88e-05 | 0.00012 | 0 | 2.88e-05 |
| 0.75 | 0.25 | 0 | 0 | 0 | 0 | 0 |
| 1 | 0.5 | 0 | 0 | 0 | 0 | 0 |
| 1.5 | 1 | 0 | 0 | 0 | 0 | 0 |

### libc++, double

| κ | shape(x2) | lost, `ratio` | lost, `split` | recovered by fix | honest overflow | still spurious |
|---|---|---|---|---|---|---|
| 0.5001 | 0.0001 | 0.93176 | 0.92849 | 0.00327 | 0.86799 | 0.0605 |
| 0.501 | 0.001 | 0.49217 | 0.47512 | 0.01704 | 0.2419 | 0.23323 |
| 0.505 | 0.005 | 0.02879 | 0.02411 | 0.00469 | 0.00082 | 0.02329 |
| 0.51 | 0.01 | 0.00082 | 0.00057 | 0.00025 | 1.00e-06 | 0.00057 |
| 0.55 | 0.05 | 0 | 0 | 0 | 0 | 0 |
| 0.6 | 0.1 | 0 | 0 | 0 | 0 | 0 |
| 0.75 | 0.25 | 0 | 0 | 0 | 0 | 0 |
| 1 | 0.5 | 0 | 0 | 0 | 0 | 0 |
| 1.5 | 1 | 0 | 0 | 0 | 0 | 0 |

### libstdc++, float

| κ | shape(x2) | lost, `ratio` | lost, `split` | recovered by fix | honest overflow | still spurious |
|---|---|---|---|---|---|---|
| 0.5001 | 0.0001 | 0.9912 | 0.9897 | 0.0015 | 0.98243 | 0.00727 |
| 0.501 | 0.001 | 0.91564 | 0.90177 | 0.01387 | 0.8379 | 0.06387 |
| 0.505 | 0.005 | 0.64387 | 0.59702 | 0.04685 | 0.41309 | 0.18394 |
| 0.51 | 0.01 | 0.41431 | 0.35619 | 0.05812 | 0.17063 | 0.18555 |
| 0.55 | 0.05 | 0.01219 | 0.00574 | 0.00645 | 0.00015 | 0.0056 |
| 0.6 | 0.1 | 0.00015 | 3.12e-05 | 0.00012 | 0 | 3.12e-05 |
| 0.75 | 0.25 | 0 | 0 | 0 | 0 | 0 |
| 1 | 0.5 | 0 | 0 | 0 | 0 | 0 |
| 1.5 | 1 | 0 | 0 | 0 | 0 | 0 |

### libstdc++, double

| κ | shape(x2) | lost, `ratio` | lost, `split` | recovered by fix | honest overflow | still spurious |
|---|---|---|---|---|---|---|
| 0.5001 | 0.0001 | 0.93145 | 0.92815 | 0.0033 | 0.86779 | 0.06036 |
| 0.501 | 0.001 | 0.49204 | 0.475 | 0.01704 | 0.24185 | 0.23314 |
| 0.505 | 0.005 | 0.02877 | 0.02416 | 0.00461 | 0.00084 | 0.02332 |
| 0.51 | 0.01 | 0.00084 | 0.00058 | 0.00026 | 6.00e-07 | 0.00058 |
| 0.55 | 0.05 | 0 | 0 | 0 | 0 | 0 |
| 0.6 | 0.1 | 0 | 0 | 0 | 0 | 0 |
| 0.75 | 0.25 | 0 | 0 | 0 | 0 | 0 |
| 1 | 0.5 | 0 | 0 | 0 | 0 | 0 |
| 1.5 | 1 | 0 | 0 | 0 | 0 | 0 |

## Released class, end to end (`bi_kappa_distribution<T>`, `no_cap()`)

| stdlib | precision | κ | non-finite fraction | exceptions | max finite log R |
|---|---|---|---|---|---|
| libc++ | double | 0.5001 | 0.9283 | 0 | 372.8 |
| libc++ | double | 0.501 | 0.47511 | 0 | 373.1 |
| libc++ | double | 0.505 | 0.02422 | 0 | 372.9 |
| libc++ | double | 0.51 | 0.00059 | 0 | 372.9 |
| libc++ | double | 0.55 | 0 | 0 | 177.4 |
| libc++ | double | 0.6 | 0 | 0 | 88.81 |
| libc++ | double | 0.75 | 0 | 0 | 35.76 |
| libc++ | double | 1 | 0 | 0 | 18.19 |
| libc++ | double | 1.5 | 0 | 0 | 9.529 |
| libc++ | float | 0.5001 | 0.98966 | 0 | 52.31 |
| libc++ | float | 0.501 | 0.90177 | 0 | 52.37 |
| libc++ | float | 0.505 | 0.59637 | 0 | 52.48 |
| libc++ | float | 0.51 | 0.35537 | 0 | 52.51 |
| libc++ | float | 0.55 | 0.00564 | 0 | 52.36 |
| libc++ | float | 0.6 | 3.26e-05 | 0 | 51.98 |
| libc++ | float | 0.75 | 0 | 0 | 29.31 |
| libc++ | float | 1 | 0 | 0 | 14.75 |
| libc++ | float | 1.5 | 0 | 0 | 8.077 |
| libstdc++ | double | 0.5001 | 0.92818 | 0 | 372.9 |
| libstdc++ | double | 0.501 | 0.47493 | 0 | 372.9 |
| libstdc++ | double | 0.505 | 0.02427 | 0 | 373 |
| libstdc++ | double | 0.51 | 0.00056 | 0 | 372.6 |
| libstdc++ | double | 0.55 | 0 | 0 | 160.6 |
| libstdc++ | double | 0.6 | 0 | 0 | 81.22 |
| libstdc++ | double | 0.75 | 0 | 0 | 33.43 |
| libstdc++ | double | 1 | 0 | 0 | 17.36 |
| libstdc++ | double | 1.5 | 0 | 0 | 7.781 |
| libstdc++ | float | 0.5001 | 0.98967 | 0 | 52.3 |
| libstdc++ | float | 0.501 | 0.9018 | 0 | 52.4 |
| libstdc++ | float | 0.505 | 0.59681 | 0 | 52.51 |
| libstdc++ | float | 0.51 | 0.35576 | 0 | 52.57 |
| libstdc++ | float | 0.55 | 0.00573 | 0 | 52.38 |
| libstdc++ | float | 0.6 | 3.30e-05 | 0 | 52.24 |
| libstdc++ | float | 0.75 | 0 | 0 | 29.92 |
| libstdc++ | float | 1 | 0 | 0 | 15.78 |
| libstdc++ | float | 1.5 | 0 | 0 | 8.032 |

## Q4 — supported range

| stdlib | precision | κ | finite fraction | spurious | honest | status |
|---|---|---|---|---|---|---|
| libc++ | double | 0.5001 | 0.07151 | 0.0605 | 0.86799 | unsupported |
| libc++ | double | 0.501 | 0.52488 | 0.23323 | 0.2419 | unsupported |
| libc++ | double | 0.505 | 0.97589 | 0.02329 | 0.00082 | unsupported |
| libc++ | double | 0.51 | 0.99943 | 0.00057 | 1.00e-06 | degraded |
| libc++ | double | 0.55 | 1 | 0 | 0 | supported |
| libc++ | double | 0.6 | 1 | 0 | 0 | supported |
| libc++ | double | 0.75 | 1 | 0 | 0 | supported |
| libc++ | double | 1 | 1 | 0 | 0 | supported |
| libc++ | double | 1.5 | 1 | 0 | 0 | supported |
| libc++ | float | 0.5001 | 0.01034 | 0.00725 | 0.98242 | unsupported |
| libc++ | float | 0.501 | 0.09803 | 0.06407 | 0.8379 | unsupported |
| libc++ | float | 0.505 | 0.40294 | 0.1839 | 0.41316 | unsupported |
| libc++ | float | 0.51 | 0.64396 | 0.18538 | 0.17066 | unsupported |
| libc++ | float | 0.55 | 0.9943 | 0.00556 | 0.00014 | degraded |
| libc++ | float | 0.6 | 0.99997 | 2.88e-05 | 0 | degraded |
| libc++ | float | 0.75 | 1 | 0 | 0 | supported |
| libc++ | float | 1 | 1 | 0 | 0 | supported |
| libc++ | float | 1.5 | 1 | 0 | 0 | supported |
| libstdc++ | double | 0.5001 | 0.07185 | 0.06036 | 0.86779 | unsupported |
| libstdc++ | double | 0.501 | 0.525 | 0.23314 | 0.24185 | unsupported |
| libstdc++ | double | 0.505 | 0.97584 | 0.02332 | 0.00084 | unsupported |
| libstdc++ | double | 0.51 | 0.99942 | 0.00058 | 6.00e-07 | degraded |
| libstdc++ | double | 0.55 | 1 | 0 | 0 | supported |
| libstdc++ | double | 0.6 | 1 | 0 | 0 | supported |
| libstdc++ | double | 0.75 | 1 | 0 | 0 | supported |
| libstdc++ | double | 1 | 1 | 0 | 0 | supported |
| libstdc++ | double | 1.5 | 1 | 0 | 0 | supported |
| libstdc++ | float | 0.5001 | 0.0103 | 0.00727 | 0.98243 | unsupported |
| libstdc++ | float | 0.501 | 0.09823 | 0.06387 | 0.8379 | unsupported |
| libstdc++ | float | 0.505 | 0.40298 | 0.18394 | 0.41309 | unsupported |
| libstdc++ | float | 0.51 | 0.64381 | 0.18555 | 0.17063 | unsupported |
| libstdc++ | float | 0.55 | 0.99426 | 0.0056 | 0.00015 | degraded |
| libstdc++ | float | 0.6 | 0.99997 | 3.12e-05 | 0 | degraded |
| libstdc++ | float | 0.75 | 1 | 0 | 0 | supported |
| libstdc++ | float | 1 | 1 | 0 | 0 | supported |
| libstdc++ | float | 1.5 | 1 | 0 | 0 | supported |

## Q5 — interaction with the component-wise cap

`non-finite per attempt` counts draws generated inside the rejection loop that were not finite. `non-finite returned` counts those that reached the caller. A non-finite draw can never satisfy the box predicate, so the cap **redraws** it.

| stdlib | precision | κ | λ | acceptance | non-finite per attempt | non-finite returned | exhausted |
|---|---|---|---|---|---|---|---|
| libc++ | double | 0.5001 | 5 | 0.00037 | 0.92823 | 0 | 0 |
| libc++ | double | 0.5001 | 20 | 0.00065 | 0.92817 | 0 | 0 |
| libc++ | double | 0.51 | 5 | 0.03602 | 0.00055 | 0 | 0 |
| libc++ | double | 0.51 | 20 | 0.06211 | 0.00058 | 0 | 0 |
| libc++ | double | 0.55 | 5 | 0.16463 | 0 | 0 | 0 |
| libc++ | double | 0.55 | 20 | 0.26978 | 0 | 0 | 0 |
| libc++ | double | 0.75 | 5 | 0.56789 | 0 | 0 | 0 |
| libc++ | double | 0.75 | 20 | 0.7845 | 0 | 0 | 0 |
| libc++ | double | 1.5 | 5 | 0.94206 | 0 | 0 | 0 |
| libc++ | double | 1.5 | 20 | 0.99671 | 0 | 0 | 0 |
| libc++ | float | 0.5001 | 5 | 0.00037 | 0.98971 | 0 | 0 |
| libc++ | float | 0.5001 | 20 | 0.00064 | 0.98973 | 0 | 0 |
| libc++ | float | 0.51 | 5 | 0.03632 | 0.35554 | 0 | 0 |
| libc++ | float | 0.51 | 20 | 0.06232 | 0.35551 | 0 | 0 |
| libc++ | float | 0.55 | 5 | 0.167 | 0.00576 | 0 | 0 |
| libc++ | float | 0.55 | 20 | 0.27762 | 0.00611 | 0 | 0 |
| libc++ | float | 0.75 | 5 | 0.57389 | 0 | 0 | 0 |
| libc++ | float | 0.75 | 20 | 0.78321 | 0 | 0 | 0 |
| libc++ | float | 1.5 | 5 | 0.94447 | 0 | 0 | 0 |
| libc++ | float | 1.5 | 20 | 0.99651 | 0 | 0 | 0 |
| libstdc++ | double | 0.5001 | 5 | 0.00037 | 0.92824 | 0 | 0 |
| libstdc++ | double | 0.5001 | 20 | 0.00065 | 0.92821 | 0 | 0 |
| libstdc++ | double | 0.51 | 5 | 0.03603 | 0.00061 | 0 | 0 |
| libstdc++ | double | 0.51 | 20 | 0.06269 | 0.00051 | 0 | 0 |
| libstdc++ | double | 0.55 | 5 | 0.16576 | 0 | 0 | 0 |
| libstdc++ | double | 0.55 | 20 | 0.27109 | 0 | 0 | 0 |
| libstdc++ | double | 0.75 | 5 | 0.57205 | 0 | 0 | 0 |
| libstdc++ | double | 0.75 | 20 | 0.78468 | 0 | 0 | 0 |
| libstdc++ | double | 1.5 | 5 | 0.93967 | 0 | 0 | 0 |
| libstdc++ | double | 1.5 | 20 | 0.99671 | 0 | 0 | 0 |
| libstdc++ | float | 0.5001 | 5 | 0.00036 | 0.98969 | 0 | 0 |
| libstdc++ | float | 0.5001 | 20 | 0.00063 | 0.98969 | 0 | 0 |
| libstdc++ | float | 0.51 | 5 | 0.03689 | 0.35512 | 0 | 0 |
| libstdc++ | float | 0.51 | 20 | 0.06322 | 0.35355 | 0 | 0 |
| libstdc++ | float | 0.55 | 5 | 0.16776 | 0.00562 | 0 | 0 |
| libstdc++ | float | 0.55 | 20 | 0.27519 | 0.00539 | 0 | 0 |
| libstdc++ | float | 0.75 | 5 | 0.56561 | 0 | 0 | 0 |
| libstdc++ | float | 0.75 | 20 | 0.78425 | 0 | 0 | 0 |
| libstdc++ | float | 1.5 | 5 | 0.93659 | 0 | 0 | 0 |
| libstdc++ | float | 1.5 | 20 | 0.99592 | 0 | 0 | 0 |

## Q3 — log-domain construction, distributional validation

`W = 1/(1+T) ~ Beta(κ−1/2, 3/2)`, obtained as `expit(−2 log R)` so that `T` is never formed. Quantiles are compared in `log R`. Threshold α = 0.01 fixed in advance.

| κ | seed | n | √n·D | KS p | CvM p | max |Δ log R| over p=0.5…0.999 | verdict |
|---|---|---|---|---|---|---|---|
| 0.55 | 4001 | 200,000 | 0.745 | 0.635 | 0.692 | 0.5219 | PASS |
| 0.55 | 4002 | 200,000 | 0.424 | 0.994 | 0.996 | 1.0936 | PASS |
| 0.55 | 4003 | 200,000 | 0.898 | 0.395 | 0.354 | 0.0758 | PASS |
| 0.75 | 4001 | 200,000 | 0.914 | 0.374 | 0.298 | 0.0222 | PASS |
| 0.75 | 4002 | 200,000 | 0.669 | 0.762 | 0.521 | 0.0524 | PASS |
| 0.75 | 4003 | 200,000 | 0.519 | 0.950 | 0.971 | 0.2086 | PASS |
| 2 | 4001 | 200,000 | 1.214 | 0.105 | 0.062 | 0.0262 | PASS |
| 2 | 4002 | 200,000 | 0.687 | 0.732 | 0.851 | 0.0300 | PASS |
| 2 | 4003 | 200,000 | 1.249 | 0.088 | 0.133 | 0.0097 | PASS |
| 5 | 4001 | 200,000 | 0.850 | 0.465 | 0.399 | 0.0123 | PASS |
| 5 | 4002 | 200,000 | 0.661 | 0.773 | 0.719 | 0.0098 | PASS |
| 5 | 4003 | 200,000 | 0.981 | 0.290 | 0.187 | 0.0048 | PASS |

