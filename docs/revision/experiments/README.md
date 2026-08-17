# Cross-experiment notes

The experiments themselves do **not** live here. Each one is a self-contained, runnable
directory at the repository root, with its own `README.md`, committed scripts, fixed seeds,
`raw/manifest.csv` and `results/`:

| Experiment | Directory | Answers | Status |
|---|---|---|---|
| 1 — radial, directional, anisotropic, frame | `experiments/exp1_radial_directional/` | R1.3 (primary); R1.5, R2.A2, R2.A3, R2.C1, R2.C3 | **complete** |
| 2 — capped vs uncapped characterization | `experiments/exp2_cap_characterization/` | R1.1, R1.2, R1.3 item 5 | **complete** |
| 3 — performance benchmark | *not created* | R1.4, R2.A2 | not started — prerequisite below |
| 4 — finite-precision / low-κ audit | `experiments/exp4_precision/` | R1.3 (low-κ range), R2.B, the C5 decision | **complete** |

This directory holds only notes that span more than one experiment.

## Standing conventions

Every experiment in this programme must satisfy these. They are not stylistic.

1. **The released C++ header is the object under test.** NumPy is a comparison or reference
   only, never the thing being validated.
2. **Committed scripts, fixed seeds, one documented command.** No notebook-only hidden state,
   no hand-edited numerical output.
3. **Machine-readable manifest** in `raw/manifest.csv` or equivalent, plus an `environment`
   block in the results JSON recording compiler, version, standard library, architecture,
   RNG, library versions, git commit and working-tree dirty flag.
4. **Raw and summarized results stay distinguishable.** Bulk `raw/*.bin` is regenerable and
   gitignored; manifests, compact summaries and manuscript-quality figures are tracked. See the
   enumerated `.gitignore` exceptions at the bottom of `/.gitignore`.
5. **Thresholds are fixed before the run and not adjusted afterwards.** If a test fails,
   diagnose the cause before touching the threshold. If the test measured the wrong quantity,
   replace the measurement and say so — that is not the same as tuning.

## Numerical rules, learned the hard way

These came out of Experiments 1 and 4 and are binding on all later work:

- **Never form `T = R²`.** It overflows for `R > 1.3×10¹⁵⁴` even where `R` is representable,
  silently discarding exactly the heavy-tail draws a low-κ diagnostic exists to test.
- **Never use `Y = T/(1+T) ~ Beta(3/2, κ−1/2)` at low κ.** Its mass piles up against 1.0 where
  no relative resolution remains. At κ = 0.55, 16.4% of values round to exactly 1.0 and KS
  reports a spurious `√n·D = 51.8`.
- **Use `W = 1/(1+T) ~ Beta(κ−1/2, 3/2)`**, ideally computed as `expit(−2 log R)` so `T` is
  never materialized. Same data, same test: `√n·D = 0.751`.
- **Never use `np.linalg.norm`** or a naive sum of squares for the radius — both square
  internally. Use `hypot` chains.
- **Report moments only where the untruncated moment exists** (κ > 3/2). Below that, the
  comparison has no reference; say so rather than silently omitting it.

## Experiment 3 prerequisite

Not started, deliberately. The benchmark may only run once **every compared method has passed
its own distributional validation** — timing an incorrect sampler is worse than not timing at
all. Current state:

- Our Gamma-ratio implementation: **validated** (Experiment 1).
- Abdul & Mace (2015) Student-*t* route: now **fully specified from the primary source**
  (Eq. 22, `X = μ + σ√(ν/χ²_ν)·Z`, ν = 2κ−1, σ² = κθ²/(2κ−1)) and implementable, but not yet
  implemented or validated.
- Zenitani piecewise-rejection route: **not yet transcribed** from ZUM2026.

Note the published constraint the benchmark would be testing: Zenitani (2025) §3 states that
his method needs ≈4.5–4.7 uniform variates per particle versus the standard method's three
normals plus one gamma, and therefore *"would be computationally less expensive than the
standard method."* **If rejection wins, we report that.**
