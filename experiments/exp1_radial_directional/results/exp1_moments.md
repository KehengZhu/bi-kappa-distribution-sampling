# Experiment 1 large-sample variances -- kappa = 2, 5, 10

100 runs x 500000 draws per kappa, the exp1_cells.py seeds; (theta_perp, theta_par) = (1, 2), B || z, uncapped, double precision. Mean of the run variances +/- standard error (sd / sqrt(100)); sd = spread of single runs.

| kappa | component | expected | mean +/- SE | sd of runs | diff % | diff / SE |
|---|---|---|---|---|---|---|
| 2 | v_perp1 | 2.0000 | 1.9916 +/- 0.0067 | 0.0666 | -0.42 | -1.26 |
| 2 | v_perp2 | 2.0000 | 1.9925 +/- 0.0078 | 0.0783 | -0.37 | -0.96 |
| 2 | v_par | 8.0000 | 8.1369 +/- 0.1606 | 1.6059 | +1.71 | +0.85 |
| 5 | v_perp1 | 0.7143 | 0.7146 +/- 0.0002 | 0.0017 | +0.04 | +1.53 |
| 5 | v_perp2 | 0.7143 | 0.7143 +/- 0.0002 | 0.0017 | -0.00 | -0.14 |
| 5 | v_par | 2.8571 | 2.8563 +/- 0.0008 | 0.0079 | -0.03 | -1.04 |
| 10 | v_perp1 | 0.5882 | 0.5882 +/- 0.0001 | 0.0012 | -0.01 | -0.33 |
| 10 | v_perp2 | 0.5882 | 0.5882 +/- 0.0001 | 0.0011 | -0.00 | -0.16 |
| 10 | v_par | 2.3529 | 2.3531 +/- 0.0004 | 0.0043 | +0.01 | +0.28 |

Same draws as the cell test: kappa = 2: True, kappa = 5: True, kappa = 10: True
