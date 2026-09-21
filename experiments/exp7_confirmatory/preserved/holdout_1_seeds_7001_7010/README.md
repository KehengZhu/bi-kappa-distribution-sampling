# Holdout 1 — seeds 7001–7010, protocol 1.3.0 — **NO-GO**

This directory is a copy, made for convenience, of the decisive artifacts of the first
confirmatory holdout. The authoritative record is commit `45d3ef8`, which holds the complete
tree; nothing here may be regenerated, corrected or re-run.

**Verdict: NO-GO.** G0 PASS, G1 FAIL, G2 FAIL, G3 FAIL, G4 OPEN, G5 PASS, G6 FAIL.

Two of those failures were real and one was not, and the distinction is the whole content of
protocol 2.0.0 §2.6:

- **G1 and G2 — a real implementation defect.** One draw in 14 233 536 audited attempts, at
  `float kappa = 0.51`, whose `log R` equalled `log(FLT_MAX)` bit for bit. It passed the
  loader's `log R <= log(max())` test, overflowed when it was exponentiated, and was returned
  as `(-inf, +inf, +inf)` without being counted — while all three components it should have
  produced were representable. The two recorded avoidable losses and the two oracle
  disagreements are that one draw seen in the libc++ and libstdc++ streams, not four events.
  Its primitives are in `oracle_disagreements.jsonl` and are now compiled into the regression
  fixtures of `cpp/test_suite.H::test_representability_boundary`.
- **G3 — not a defect of the candidate.** The three F5 rejections in `validation_matrix.md`
  (C3's excess at `q0 = 1e-3`, C4's excess at `q0 = 1e-2` and `1e-3`) came from testing
  representability-censored survivors against an uncensored `Exp(1)` null. Measured afterwards
  over 2000 replicates of a perfectly correct loader, that null rejects with probability
  **1.000** on those cells. See `docs/revision/experiments/f5_tail_calibration.md`.
- **G4 OPEN** on absent native x86_64 evidence, and **G6 FAIL** on `make verify` never having
  been run.

**Seeds 7001–7010 are spent.** PROTOCOL.md §8 forbids recomputing any result on them. They may
serve only as preserved failure evidence and as development material, which is what the
calibration above used them for. `make selftest` asserts they appear in no block in use.
