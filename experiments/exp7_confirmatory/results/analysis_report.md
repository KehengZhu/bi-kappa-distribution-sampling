NO-GO

# Experiment 7 - confirmatory analysis report

Generated deterministically from the recorded run; this report carries no wall-clock stamp, so `make reverify` can diff it. Run times and build identity are in `provenance.md`.

Protocol `config/protocol.json` version 2.0.0, SHA-256 `914c5d1d2b6900d872402a545f46a99819038987416d1b882909eb4ca2bafa48`.

## Gates

| gate | name | verdict | evidence |
|---|---|---|---|
| G0 | primary-source audit | PASS | log-q evaluator agrees with the 60-digit oracle to 2.84e-14 (tolerance 1e-10); attribution checklist complete for 7/7 cited methods; acceptance-rate predictor pre-registered as 'area_ratio' |
| G1 | scalar correctness | FAIL | 3/4 scalar families pass (failed: ['F2_quantile_coverage']); candidate avoidable loss 2 (required exactly 0); benign control exercises the candidate: True |
| G2 | mechanism closure | PASS | 0 accounting-identity failures; 0 oracle disagreements in 14235020 adjudicated attempts; 0 honest-count discrepancies left unadjudicated |
| G3 | complete-loader fidelity | PASS | F5 passes jointly over its cells (Holm jointly over 52 tests across all cells; 18 cells are conditional on success and are labelled as such); F6 passes (5 controls; required power 0.9 at the pre-registered effect); 0 conditional cells lack a loss fraction |
| G4 | portability | PASS | 2 native environments completed ([('arm64', 'libc++'), ('arm64', 'libstdc++')]); 0 translated results recorded as corroborating only and excluded from the decision; 1 bitwise cross-standard-library comparisons; acceptance is limited to the supported environment; the cross-architecture claim is withdrawn from the claim boundary and published as a limitation, not assumed; 1 bitwise cross-stdlib comparisons (exact); 0 cross-architecture cells, of which 0 are equivalent within +/-0.15 on the log rate ratio, 0 disagree and 0 are underpowered (too few events to establish either, so they leave the gate open rather than closing or failing it) |
| G5 | operational viability | PASS | candidate/LEGACY time per returned sample 0.854 (cluster-bootstrap upper limit 0.906; pre-registered bound 2); same-seed reproducible: True; RNG-stream break documented: True |
| G6 | reproducible artifact | PASS | make verify exit 0; derived artifacts regenerate byte-identically: True; dependency set clean: True; baseline comparison resolved: True; protocol hash matches the frozen document: True; archive identifier bi-kappa 2.1.0 @ e5c98378d4e729395f9dc4c5a12c23ca15d50698 |

## Families

| family | alpha | passed | statistic | p | detail |
|---|---:|---|---:|---:|---|
| F1_radial_law | 0.01 | yes | 0.3586 | 0.3586 | Simes over 26 configurations; 0 identified by Holm on rejection |
| F2_quantile_coverage | 0.005 | NO | 23 | 1.464e-07 | 23 misses of 650 resolved intervals; expected 6.21; fail above 13 |
| F3_quantile_direction | 0.01 | yes | 0.8652 | 0.8652 | Holm over 130 (kappa, precision, p) cells |
| F4_upper_tail_mass | 0.01 | yes | 0.1938 | 0.1938 | Simes over 780 exceedance statistics |
| F5_loader_battery | 0.01 | yes | 0.8918 | 0.8918 | Holm jointly over 52 tests across all cells; 18 cells are conditional on success and are labelled as such |
| F6_negative_controls | - | yes | - | - | 5 controls; required power 0.9 at the pre-registered effect |
| F7_portability | 0.005 | yes | 0 | 1 | 1 bitwise cross-stdlib comparisons (exact); 0 cross-architecture cells, of which 0 are equivalent within +/-0.15 on the log rate ratio, 0 disagree and 0 are underpowered (too few events to establish either, so they leave the gate open rather than closing or failing it) |

## What the verdict means

At least one gate a manuscript claim requires did not pass. The failure is preserved. Per PROTOCOL.md §8 the only legitimate next step is to identify a concrete implementation defect, fix it, freeze a new implementation hash and a new protocol, draw a further disjoint seed block, and rerun. Repeating the run to obtain a more favourable result is not permitted.

Expected false-failure rate of this design, computed from the frozen family structure before the run: at most 0.05.

## Unresolved

- **G1 (scalar correctness) - FAIL.** 3/4 scalar families pass (failed: ['F2_quantile_coverage']); candidate avoidable loss 2 (required exactly 0); benign control exercises the candidate: True

## Decisive source data

- `PROTOCOL.md`, `config/protocol.json`, `config/power_study.json`, `config/honest_floor.md` - the frozen rules and what they can detect
- `results/scalar_validation.csv`, `results/tail_metrics.csv` (F1-F4)
- `results/failure_envelope.csv`, `results/oracle_audit.jsonl` (G2)
- `results/loader_validation.csv`, `results/negative_controls.csv` (F5, F6)
- `results/portability.csv`, `results/portability_remote.md` (F7, G4)
- `results/performance.csv` (G5)
- `results/source_data_README.md` defines every column.
