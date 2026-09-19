PARTIAL

# Experiment 6 - analysis report

Generated 2026-09-19T17:46:11+00:00.

## Gates

| gate | name | verdict | evidence |
|---|---|---|---|
| G0 | primary-source audit | PASS | config/pilot.json records the exact equations, stated parameter domain, output scale and limitations of both candidates with their DOIs; the frozen log-q evaluator agrees with an arbitrary-precision incomplete beta to better than 1e-10; the LOG-PUB transcription reproduces the acceptance rate implied by the paper's own target and envelope and the sampled law satisfies E[log X] = psi(alpha). |
| G1 | scalar correctness | PASS | selected primitive LOG-PUB: 0 log-domain failures, 0 of 84 pilot configurations with a Holm-corrected rejection. |
| G2 | mechanism closure | PASS | 0 accounting-identity failures across 132 configurations; every method's honest count agrees with the reference floor to within the audited near-threshold band (0 exceptions); the log path loses 0 representable draws in double across the whole ladder; the independent cpp_dec_float_100 oracle adjudicated 589964 audited attempts with 0 disagreements. |
| G3 | complete-loader fidelity | PASS | 42 validation cells; 0 of the candidate's cells contain a failure; negative controls detected; 84 tail-metric rows quantify the conditioning. 8 released-path cells reject the radial law on their finite survivors, which is the conditioning result rather than a defect. |
| G4 | portability | OPEN | completed: [('arm64', 'libc++'), ('arm64', 'libstdc++')]; not available on this host: [('x86_64', 'libc++'), ('x86_64', 'libstdc++')]. An unavailable architecture is recorded with the exact command to run it elsewhere; emulation is not accepted as a substitute, so this gate stays open. |
| G5 | operational viability | PASS | worst paired LOG/SPLIT time per returned sample: 1.42x; same-seed reproducibility True. Whether that cost is acceptable for a header-only C++11 reference implementation is an author decision, not a statistical one. |
| G6 | reproducible artifact | PASS | dependency set clean: True; commit 5b5393e74bcf; exploratory flag: False. |

## What the verdict means

The gates that a manuscript mitigation claim requires (G0-G3, G6) passed; at least one gate that production adoption requires did not.

The log path stays an unreleased prototype. The measured operating envelope, the mechanism decomposition and the conditioning result stand on their own and are reportable within the boundary those gates establish. The released default in `cpp/bi_kappa_distribution.H` is not changed by this experiment, and no public API is altered.

## Headline numbers

- Selected log primitive: **LOG-PUB** - both candidates passed; LOG-PUB has the lower median time per accepted log-Gamma attempt (6.63e-08 s vs 9.69e-08 s). LOG-ID is additionally rejection-free and defined for every shape, while LOG-PUB carries a published acceptance-rejection envelope whose stated domain is 0 < a < 1; the cost tie-break is recorded rather than overridden.
- Accounting-identity failures: 0 across 132 paired configurations.
- Independent 100-digit oracle: 589964 audited attempts, 0 disagreements, 0 conversion failures.
- Complete-loader cells: 42, of which 0 of the candidate's contain a failure; negative controls detected.
- Released-path cells whose radial law is rejected on the finite survivors: 8. That is the conditioning result: discarding the failures changes the law of what is kept, and the log path passes the same test on the same configuration.

## Accuracy of the surviving draws

The oracle also re-derives each formation's radius exactly, which is how this experiment can say something the counters cannot. Surviving a draw is not the same as getting it right. On the 33,326 audited draws whose denominator Gamma landed in the subnormal range - the draws the split form is closest to losing altogether - the split form's returned radius carries up to 0.484 relative error, against 7.76e-06 for the log path. A subnormal denominator has already lost most of its significand, and the square root halves the exponent without restoring it; the log path never materializes the denominator at all. This is conditional on the denominator being subnormal, a stratum the audit covers at a fixed sampling rate, and it says nothing about draws where the denominator is normal - there the two agree to rounding.

## Claim boundary

Set by the independent novelty audit, not by this result. The Gamma-ratio Kappa construction, log-domain Gamma generation, the small-shape Gamma underflow and the low-parameter denominator hazard are all established prior art and are attributed as such.

Supportable on the evidence here:

- Forming the intermediate quotient `X1/X2` can overflow even where its square root is representable; evaluating the algebraically identical `sqrt(X1)/sqrt(X2)` removes that mode, and in paired tests the two formulations receive identical variates.
- That change does not remedy zeros produced by the small-shape Gamma primitive and does not extend the mathematical support of the target law.
- Working in the log domain is established for small-shape Gamma variates; the contribution here is its integration and validation in this loader.
- Failure rates are reported as observed counts with intervals, and a configuration with no observed failure is reported as `no failures observed in N draws under the tested configuration`, with its one-sided upper bound.

Not supportable, and not claimed:

- discovery of the small-shape Gamma underflow problem, or of division by zero in a low-parameter Kappa loader;
- a first or novel log-domain Gamma generator, or a novel Gamma-ratio, Beta-prime, Student-t or rejection Kappa sampler;
- that the split expression solves the low-kappa finite-precision problem;
- that the log representation is exact - it is a stable representation carrying about one ulp of `(log U)/a`;
- that the log path removes all avoidable loss, beyond what the final-vector oracle has closed;
- any universal mathematical lower bound on kappa, or `reliable down to`, without a sample size, platform, failure criterion and confidence bound;
- that both standard libraries use an identical `Y U^(1/a)` implementation;
- complete-loader correctness inferred from tests run after deleting non-finite outputs. Conditional results are labelled conditional throughout.

## Unresolved exceptions

- G4 (portability): completed: [('arm64', 'libc++'), ('arm64', 'libstdc++')]; not available on this host: [('x86_64', 'libc++'), ('x86_64', 'libstdc++')]. An unavailable architecture is recorded with the exact command to run it elsewhere; emulation is not accepted as a substitute, so this gate stays open.

## Decisive source data

- `baseline_reproduction.md`, `baseline_envelope.csv` (E0)
- `pilot_decision.md`, `scalar_validation.csv`, `log_q_evaluator_validation.csv` (E1)
- `failure_envelope.csv` (E2), `oracle_audit.jsonl`, `oracle_disagreements.jsonl`
- `conditioning_bins.csv`, `tail_metrics.csv` (E3)
- `loader_validation.csv`, `validation_matrix.md` (E4)
- `portability.csv` (E5), `performance.csv` (E6)
- `source_data_README.md` defines every column.
