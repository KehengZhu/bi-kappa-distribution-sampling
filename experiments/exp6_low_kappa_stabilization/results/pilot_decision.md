# E1 - pilot decision

Generated 2026-09-19T17:43:58+00:00.

## Frozen evaluator check, run before any sampling result was inspected

The two-piece `log q` evaluator (switch at `log W = -30`) was compared
against an `mpmath` incomplete beta at 60 decimal digits on a 41-point grid spanning
10 natural-log units on each side of the switch, for every pilot shape.

**PASS** - worst absolute `log q` disagreement 2.13e-14, tolerance 1e-10.
Source data: `log_q_evaluator_validation.csv`.

## Candidates

| candidate | configs | in-domain | log failures | configs with a rejection | quantile intervals covering target | median s per log-Gamma attempt | defined at 0.5001 | passes |
|---|---:|---:|---:|---:|---|---:|---|---|
| LOG-ID | 84 | 84 | 0 | 0 | 419/420 | 9.69e-08 | yes | yes |
| LOG-PUB | 84 | 72 | 0 | 0 | 417/420 | 6.63e-08 | yes | yes |

## Transcription check on LOG-PUB

Liu, Martin and Syring (2017) state the acceptance rate of their scheme as
`r(alpha) = {1 + w(alpha)}^-1` (Eq. 2, p. 1771).  Their target `h_alpha` is
normalized (`c^-1 = Gamma(alpha+1)`) while the envelope `eta_alpha` integrates to
`c(1+w)`, so the ratio of areas - which is what an acceptance-rejection scheme's
measured rate converges to - is `Gamma(alpha+1)/{1+w(alpha)}`.  The two agree to
leading order in the small-shape regime the paper addresses, and the paper's own
prose identifies `1/(1+w)` as the probability of *proposing* from the Exp(1) branch.

Measured against the area ratio: worst standardized deviation +2.24.
Measured against Eq. (2) as published: worst standardized deviation -127.25.

The sampled law is unaffected either way; the implementation is used as transcribed.


## Decision

**Selected primitive: LOG-PUB.**

both candidates passed; LOG-PUB has the lower median time per accepted log-Gamma attempt (6.63e-08 s vs 9.69e-08 s). LOG-ID is additionally rejection-free and defined for every shape, while LOG-PUB carries a published acceptance-rejection envelope whose stated domain is 0 < a < 1; the cost tie-break is recorded rather than overridden.

The selection was made from the pilot alone.  No later phase was read before this
file was written, and `config/frozen.json` records the choice together with the
SHA-256 of the primitive's source, which every later phase checks.
