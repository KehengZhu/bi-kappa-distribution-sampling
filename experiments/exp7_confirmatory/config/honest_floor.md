# The honest representability floor

Frozen with the protocol, before any Experiment 7 data existed.

A loader that never throws a representable draw away still cannot return a velocity that has
no representation. Near `kappa = 1/2` the bi-Kappa law puts non-zero probability outside every
finite floating-point range, so the right question is not "does the sampler lose draws?" but
"does it lose more than the type forces?" That requires knowing the floor in closed form
rather than asserting that a measured rate is small.

## Derivation

With `a = kappa - 1/2`, `X1 ~ Gamma(3/2, 1)`, `Y ~ Gamma(a+1, 1)`, `U ~ Uniform(0,1)` and
`E = -log U ~ Exp(1)`, the shape-boosting identity gives `log X2 = log Y - E/a`, so

    log R = (log X1 - log X2) / 2 = (log X1 - log Y + E/a) / 2.

Hence, for any threshold `L`,

    P(log R > L) = P(E > 2aL - a(log X1 - log Y))
                 = E[(X1/Y)^a] · e^(-2aL),

because `E` is exponential and independent of `X1` and `Y`. The prefactor is exact:

    E[X1^a]  = Gamma(3/2 + a) / Gamma(3/2)
    E[Y^(-a)] = Gamma(1) / Gamma(a+1) = 1 / Gamma(a+1)

A returned component is `V_j = R g_j`, and in the isotropic, unrotated, `theta = 1`
configuration `g_j = sqrt(kappa) n_j` with `n` uniform on the sphere. The draw overflows iff
`R · sqrt(kappa) · max_j |n_j| > MAX`, so

    P(overflow) = [ Gamma(3/2 + a) / (Gamma(3/2) Gamma(a+1)) ]
                  · MAX^(-2a)
                  · E[ (sqrt(kappa) max_j |n_j|)^(2a) ]

The last factor is an expectation over the direction alone, evaluated once to high accuracy by
quadrature or Monte Carlo; it is within a few per cent of one for every shape of interest.

## Check against the candidate

Measured over 2x10^6 draws with the released 2.0.0 loader, seed 4001, isotropic and unrotated:

| precision | kappa | analytic floor | measured | z |
|---|---:|---:|---:|---:|
| double | 0.5001 | 8.676e-01 | 8.673e-01 | −1.19 |
| double | 0.505 | 8.250e-04 | 8.350e-04 | +0.49 |
| double | 0.51 | 6.807e-07 | 1.000e-06 | +0.55 |
| float | 0.55 | 1.376e-04 | 1.420e-04 | +0.54 |
| float | 0.75 | 5.252e-20 | 0 observed | — |

Every value agrees within 1.2 standard errors. The 1.0.0 loader over the same draws gives
9.285e-01, 2.396e-02, 5.690e-04 and 5.620e-03 — between 1.07 and 29 times the floor, and the
excess is the avoidable loss this work removes.

## The floor across the frozen ladder

`double`:

| kappa | 0.5001 | 0.501 | 0.505 | 0.51 | 0.55 | 0.60 | 0.75 | 1.0 |
|---|---|---|---|---|---|---|---|---|
| floor | 8.7e-01 | 2.4e-01 | 8.3e-04 | 6.8e-07 | 1.5e-31 | 2.2e-62 | 7.2e-155 | 5.9e-309 |

`float`:

| kappa | 0.5001 | 0.501 | 0.505 | 0.51 | 0.55 | 0.60 | 0.75 | 1.0 |
|---|---|---|---|---|---|---|---|---|
| floor | 9.8e-01 | 8.4e-01 | 4.1e-01 | 1.7e-01 | 1.4e-04 | 1.9e-08 | 5.2e-20 | 3.1e-39 |

## How the gates use this

G1 requires the candidate's **avoidable** loss — observed loss minus the floor, adjudicated
per draw by the arbitrary-precision oracle rather than by subtracting two rates — to be
exactly zero. G2 requires the candidate's honest count to equal the oracle floor exactly, or
each discrepancy to be individually adjudicated and listed.

Reporting rests on the same quantity. A configuration with no observed failure is reported as
"no failures observed in N draws" with its one-sided upper bound *and* its analytic floor, so
that a reader can see whether zero observations mean the floor is negligible (double
`kappa >= 0.55`, where it is below 1e-31) or merely that N was too small to reach it.
