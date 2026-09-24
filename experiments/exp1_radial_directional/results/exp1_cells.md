# Experiment 1 large-sample cell test -- nine kappa, 5 x 10^7 draws each

100 seeds x 500000 draws per kappa (kappa index i uses seeds 100001 + 1000 i onward); (theta_perp, theta_par) = (1, 2), B || z, uncapped, double precision. Expected per cell: 125000. Overflow = draws with a component beyond the largest double (counted in the outermost shell of the radius test, excluded from the 400-cell test).

| kappa | radius p (10 shells) | all cells p (400) | overflow | NaN | seed chi^2/dof | seed KS p |
|---|---|---|---|---|---|---|
| 0.51 | 0.995 | 0.353 | 39 | 0 | 1.040 +/- 0.047 | 0.187 |
| 0.55 | 0.070 | 0.062 | 0 | 0 | 1.005 +/- 0.050 | 0.565 |
| 0.75 | 0.139 | 0.125 | 0 | 0 | 1.006 +/- 0.051 | 0.969 |
| 1 | 0.795 | 0.464 | 0 | 0 | 1.023 +/- 0.041 | 0.051 |
| 1.25 | 0.248 | 0.968 | 0 | 0 | 1.000 +/- 0.042 | 0.580 |
| 1.5 | 0.554 | 0.992 | 0 | 0 | 0.990 +/- 0.045 | 0.941 |
| 2 | 0.756 | 0.147 | 0 | 0 | 1.008 +/- 0.051 | 0.851 |
| 5 | 0.403 | 0.672 | 0 | 0 | 0.986 +/- 0.047 | 0.817 |
| 10 | 0.406 | 0.719 | 0 | 0 | 1.029 +/- 0.049 | 0.217 |

Shell deviations in binomial standard deviations:

- kappa = 0.51: +0.2 +0.2 +0.1 -0.0 -0.4 -0.5 -0.3 +0.8 +0.6 -0.6
- kappa = 0.55: +1.3 -2.5 +1.0 -0.1 +1.3 -0.7 -1.1 -1.6 +1.2 +1.1
- kappa = 0.75: +1.1 -0.2 +1.8 -0.2 -1.3 -1.5 +1.9 -0.5 +0.3 -1.6
- kappa = 1: +0.6 -1.1 -0.1 -1.5 +0.4 +1.1 +0.5 -0.6 +0.4 +0.2
- kappa = 1.25: +0.8 +0.9 -1.8 +1.2 +0.6 +0.8 +0.8 -1.0 -0.2 -2.0
- kappa = 1.5: -0.9 +0.6 -0.7 +1.4 -0.6 -0.4 +2.0 -0.4 -0.2 -0.8
- kappa = 2: -0.1 -0.9 +1.0 +0.1 +0.8 +0.5 -0.2 +0.1 -1.8 +0.6
- kappa = 5: +1.9 +1.1 -0.9 -0.1 -0.0 -1.5 +1.0 -0.9 -0.7 +0.2
- kappa = 10: -0.4 +0.2 -0.7 +1.0 -1.1 +0.4 +0.7 +1.9 -1.8 -0.3
