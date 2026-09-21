# Holdout 2 — seeds 8001–8010, protocol 2.0.0 — **NO-GO**

This directory is a copy, made for convenience, of the decisive artifacts of the second
confirmatory holdout. The authoritative record is commit `e5c9837`, which holds the complete
tree; nothing here may be regenerated, corrected or re-run.

**Verdict: NO-GO.** G0 PASS, G1 FAIL, G2 PASS, G3 PASS, G4 PASS, G5 PASS, G6 PASS.

Protocol 2.0.0's two corrections held. Four of the five gates that had failed or stood open
on 7001–7010 passed:

- **G2** — 0 oracle disagreements in 14 235 020 adjudicated attempts, against 2 before. The
  representability decision and the 100-digit adjudicator now agree everywhere.
- **G3** — F5 passes jointly over 52 tests. The three cells that failed under 1.3.0's
  untruncated nulls (C3's excess at `q0 = 1e-3`, C4's at `1e-2` and `1e-3`) pass at
  p = 0.258, 0.885 and 0.667, and the two count members 1.3.0 had to declare not applicable
  now apply and pass. NC3a, honest censoring alone, is rejected 2 times in 200; NC3b detects
  excess conditioning 200/200 at both pre-registered effect sizes.
- **G4** — both supported environments native, cross-standard-library agreement bitwise.
- **G5** — 0.854× time per returned sample, 99 % cluster-bootstrap upper limit 0.906,
  against a bound of 2.0.

**G1 failed, for two reasons, and neither was argued away at the time.**

- **One avoidable loss — a real implementation defect.** `float kappa = 0.505`, seed 8005,
  attempt 157855, one draw seen in both standard-library streams. Reproduced at 120 digits:
  the largest intended component lies `3.12e-08` natural-log units below `log(FLT_MAX)`,
  about half an ulp, so the correctly rounded `float` is finite and the draw is returnable.
  The `float` log radius carries `3.5e-06` of error, inherited from `logf(u)` and amplified
  200-fold by the division by `a = 0.005`, and the `float` order-unity vector carries another
  `3.9e-08` against a `float` ulp of `6.0e-08` there. Neither resolves the margin, so the
  loader exponentiated one ulp high and returned a non-finite component. The oracle, working
  exactly, agrees with the loader's *classification*, which is why G2 is clean and G1 is not.
  Release 2.2.0 carries the whole deterministic map in a wider accumulator and rounds once;
  see PROTOCOL.md §2.7.1. The draw is now a deterministic fixture in
  `cpp/test_suite.H::test_representability_boundary` (S5) and in `make selftest`.
- **F2 quantile coverage — not a defect of the candidate, and not of the statistic either.**
  23 misses of 650 resolved intervals against 6.21 expected, `p = 1.5e-07` under the frozen
  Poisson-binomial. That null is built over *independent* intervals and the intervals were
  not independent: every configuration of a replicate ran on one `mt19937` seeded with the
  replicate's seed, so all 26 saw the same uniforms and their quantile errors co-moved. Seed
  8004 gave 103 of 130 signed errors positive, seed 8005 gave 31 of 130, and 20 of the 23
  misses came from those two replicates; the first holdout's 0 misses of 650 is the same
  over-dispersion from the other side. F3, which tests the signed standardized quantile error
  against a Monte-Carlo-calibrated null, passed at p = 0.865. **The diagnosis was recorded and
  not acted on**: re-specifying a predicand after seeing it fail is what §8 forbids. Protocol
  3.0.0 acts on it before any new datum exists, by giving each configuration its own stream
  and computing the remaining within-configuration dependence exactly; see PROTOCOL.md
  §2.7.2. On this run's own count of 23 the replacement rule would still have rejected, at
  `p = 4.3e-07`.

**Seeds 8001–8010 are spent.** PROTOCOL.md §8 forbids recomputing any result on them. They
may serve only as preserved failure evidence and as material for diagnosing a defect, which
is what §2.7 used them for. `make selftest` asserts they appear in no block in use.
