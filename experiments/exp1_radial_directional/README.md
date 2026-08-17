# Experiment 1 — radial, directional, anisotropic and frame validation

Answers **R1.3** (primary) and supports R1.5, R2.A2, R2.A3, R2.C1, R2.C3.

This is the experiment that tests the construction the sampler is actually built on —
`T = R² ~ BetaPrime(3/2, κ−1/2)`, a uniform and independent direction, and the anisotropic
and field-frame mappings — rather than Cartesian marginals alone. It also covers
`1/2 < κ ≤ 3/2`, where the second moment does not exist and the manuscript's existing
variance-based validation cannot reach.

## Reproducing

```bash
make run                                                   # build + generate raw samples
uv run --project ../../python python exp1_analyze.py       # analyze -> results/
uv run --project ../../python python exp1_figure.py        # figure -> results/
```

`make run` regenerates every sample from scratch; nothing in `results/` depends on state
that is not recorded in `raw/manifest.csv`.

## What is run

**Mode: uncapped throughout** (`bi_kappa_distribution<double>::no_cap()`). The sampled law
is the full bi-Kappa distribution, and the core mapping executes exactly once per draw with
no acceptance–rejection loop. The capped mode is Experiment 2's subject, not this one.

| | Block A — radial core | Block B — anisotropy × frame |
|---|---|---|
| κ | 0.51, 0.55, 0.75, 1.0, 1.25, 1.5, 2, 5, 10 | 0.55, 2, 10 |
| (θ⊥, θ∥) | (1, 2) | (1, 1) and (1, 2) |
| **B** direction | `ẑ` | `ẑ`, `(1,1,1)/√3`, `(0.3,−0.5,0.8)/‖·‖` |
| seeds | 1001–1005 | 1001–1005 |
| N per run | 100 000 | 100 000 |
| runs | 45 | 90 |

Total 135 runs, 1.35×10⁷ draws. Environment (compiler, stdlib, numpy/scipy versions, RNG)
is recorded in `results/exp1_results.json`.

## Diagnostics

- **Radial law** — KS and Cramér–von Mises against the exact law, plus quantile probes at
  p = 0.5, 0.9, 0.99, 0.999 compared in `log R`.
- **Directional uniformity** — KS of `cos Θ` against U(−1,1) and of `Φ` against U(−π,π).
- **Radial–direction independence** — χ² contingency over radial quartiles × direction-cosine
  deciles, plus a two-sample KS of `cos Θ` between the innermost and outermost radial quartiles.
- **Anisotropy alignment** — ratio of median absolute deviations `MAD(v∥)/MAD(v⊥)`, which needs
  no moments and therefore works at every κ.
- **Frame invariance** — the recovered normalized radius is compared draw by draw against the
  axis-aligned run at the same seed.
- **Moments** — reported only for κ > 3/2, where they exist.

## Two implementation notes that are results in their own right

**Never form `T = R²`.** `T` overflows for `R > 1.3×10¹⁵⁴` even though `R` is perfectly
representable, so computing the diagnostic through `T` silently discards valid heavy-tail
draws exactly where the diagnostic matters. This is the same trap that the sampler itself
used to fall into at `cpp/bi_kappa_distribution.H:273`. The analysis works in `R` throughout
and uses `np.hypot`, not `np.linalg.norm` (which squares internally and overflows).

**The orientation of the bounded transform is not cosmetic.** `W = 1/(1+T) ~ Beta(κ−1/2, 3/2)`
and `Y = T/(1+T) ~ Beta(3/2, κ−1/2)` are exact bijections of `T` and carry identical
information in exact arithmetic. In doubles they do not: for small `κ−1/2` the mass of `Y`
piles up against 1, where no relative resolution remains, and 16% of `Y` values at κ = 0.55
round to exactly 1.0. A KS test on `Y` then measures rounding rather than the sampler, and
reports a spurious `√n·D = 51.8`. The same data tested on `W` gives 0.751 — a clean pass.
`W` is the correct bounded diagnostic; see `results/exp1_table.md`.

This supersedes the recommendation in `docs/reviewer_response_matrix.md` R1.3 item 2 and in
`docs/introduction_rewrite_proposal.md` §3a, both of which name `Y`.
