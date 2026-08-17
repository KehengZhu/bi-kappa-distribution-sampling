# Experiment 2 — characterization of the component-wise velocity cap

Answers **R1.1** and **R1.2** (primary) and R1.3 item 5; supports R1.5 and R1.6.

Quantifies what the finite `max_normalized_velocity` in `cpp/bi_kappa_distribution.H`
actually costs and actually changes. Experiment 1 validated the **uncapped** sampler
against the bi-Kappa law; this experiment is about the other mode.

## The target-law distinction, stated once and not blurred again

The released sampler has two modes and **they sample two different probability laws**.

| mode | `max_normalized_velocity` | law being sampled |
|---|---|---|
| uncapped | `no_cap()` (= `+inf`) | the untruncated bi-Kappa distribution — the intended target |
| capped | finite `lambda` | that distribution **conditioned on** the component-wise box `\|v_x\|/θ⊥ ≤ λ` **and** `\|v_y\|/θ⊥ ≤ λ` **and** `\|v_z\|/θ∥ ≤ λ` |

The capped mode is a **truncated/conditional target**, not a numerical approximation to
the uncapped one. It has bounded support, different tails, different moments, and a
different angular symmetry (it is a *cube* in normalized components, so it is neither
isotropic nor axisymmetric about **B**). It must never be described, plotted, or cited
as "the bi-Kappa distribution".

Everything below is a measurement of the gap between those two laws.

## Reproducing — one command each

```bash
make run                                                   # build + generate raw samples + checksums
uv run --project ../../python python exp2_analyze.py       # analyze -> results/
```

`make run` regenerates all 240 sample files from scratch; nothing in `results/` depends
on state that is not recorded in `raw/manifest.csv` and the two committed sources.
`make verify` re-checks a regenerated `raw/` against `raw/checksums.sha256`.

## Configuration matrix

Field direction is **ẑ** for every run. The cap is tested *before* the field-frame
rotation (`bi_kappa_distribution.H` step 4 vs step 5), so the box event — and hence the
acceptance probability and the whole capped law in the local frame — cannot depend on
`ub`. Experiment 1 already validates the rotation itself.

| | Block A — anisotropic ladder | Block C — isotropic control |
|---|---|---|
| κ | 0.75, 1.0, 1.5, 2.0, 5.0, 10.0 | 0.75, 2.0 |
| (θ⊥, θ∥) | (1, 2) | (1, 1) |
| λ | `no_cap()`, 3, 5, 10, 20, 50, 100 | `no_cap()`, 3, 5, 10, 20, 50, 100 |
| seeds | 2001–2005 | 2001–2005 |
| N per run | 100 000 | 100 000 |
| runs | 210 | 70 |

The λ ladder deliberately keeps **20** (the library default) and **100** (the value in
the manuscript's example) visible, as `docs/revision/planning/reviewer_response_matrix.md`
§"Experiment 2" asks.

280 runs, 2.8×10⁷ draws. `κ = 0.75` and `κ = 1.0` are the heavy-tailed cases where the
untruncated second moment does not exist; `κ = 1.5` is the boundary where it diverges.
Environment (compiler, target triple, flags, stdlib, numpy/scipy, git commit and
working-tree dirty flag) is recorded in `results/exp2_results.json`. The working tree
was dirty when these numbers were produced, for reasons outside this experiment, so the
environment block additionally pins the exact sampler by content:
`sampler_header_sha256` is the SHA-256 of `cpp/bi_kappa_distribution.H` as compiled.

## How the rejected fraction was measured

`operator()` loops internally and reports no attempt count, and
**`cpp/bi_kappa_distribution.H` was not modified.** Instead:

1. The shipped predicate `withinNormalizedVelocityCap` is transcribed verbatim into
   `exp2_analyze.py:in_box`.
2. It is evaluated on the draws of the **uncapped** run at the *same seed*. Every loop
   iteration of `operator()` consumes `x1, x2, cosTheta, phi` in the same order whether
   or not a cap is in force, so the uncapped run *is* the capped run's attempt stream.
   The mean of the predicate over it is `P(accept)` directly.
3. That correspondence is **verified, not assumed**: for all 240 (case, λ, seed) pairs
   the capped run's output is bitwise identical to the uncapped run's draws restricted
   to the box. This turns "capped = uncapped conditioned on the box" from a claim about
   the code into a measured fact about the shipped binary.

The same accept mask also yields the **retry-count distribution** (gaps between accepted
indices), which is compared against the `Geometric(p)` reference, and the probability of
hitting the sampler's internal `kMaxCapRejectTries = 10⁶` limit.

An exact closed form for `P(accept)` is also derived and used as an independent check
(and as the only usable value where the rejection rate falls below the 2×10⁻⁶ resolution
of 5×10⁵ Monte-Carlo attempts). The θ's cancel in the normalized predicate, leaving a
cube of half-side `c = λ/√κ` in the isotropic coordinates, so with `M = maxᵢ|nᵢ|` for
`n` uniform on S²,

```
P(accept) = E_M[ I_z(3/2, κ−1/2) ],   z = c²/(M² + c²),
f_M(m) = 3 − (12/π) arcsin( √((1−2m²)/(1−m²)) )   for 1/√3 ≤ m ≤ 1/√2,
f_M(m) = 3                                        for 1/√2 ≤ m ≤ 1.
```

`f_M` is validated against a 2×10⁶-point spherical Monte Carlo in the results JSON.

## Moments: what is compared and what is refused

The untruncated bi-Kappa second moment is `θ² κ/(2κ−3)`. It **exists only for κ > 3/2**
and diverges at κ = 3/2. Variance comparisons are therefore reported only for κ = 2, 5, 10
and are **explicitly refused** for κ = 0.75, 1.0, 1.5 — refused in the output, not
silently omitted, because a capped sample at those κ does have a perfectly finite
variance and a naive diagnostic will happily print it. That number is a property of λ,
not of the plasma: `results/exp2_table.md` §3 shows it growing without bound as the cap
is relaxed.

## Numerical safety rules observed

- `|v|` via chained `np.hypot`. `np.linalg.norm` squares internally and overflows on the
  κ ≤ 1 tails.
- `T = R²` is never formed. The analytic acceptance uses `z = c²/(m² + c²)`.
- The bounded radial variable `Y = T/(1+T)` is not used anywhere; it rounds to exactly 1
  at low κ (see Experiment 1). Nothing here needs a bounded radial diagnostic, but the
  rule is honoured.
- The rejection probability is integrated *directly* rather than as `1 − P(accept)`:
  at κ = 10, λ = 50 the acceptance is `1 − 9×10⁻²⁴` and the subtraction would return
  exactly zero (or a negative quadrature residue). The complement is taken inside the
  incomplete beta via `1 − I_z(a,b) = I_{1−z}(b,a)`, where `1 − z = m²/(m² + c²)` is exact.
- Thresholds (`TV < 10⁻³`, p99.9 quantile within 1%) were fixed in the script before any
  result was looked at and were not retuned.

## Artifacts: canonical vs regenerable

| path | status |
|---|---|
| `exp2_sample.cpp`, `GNUmakefile`, `exp2_analyze.py`, `README.md` | **canonical**, tracked |
| `raw/manifest.csv` | **canonical**, tracked — one row per run, full provenance |
| `raw/checksums.sha256` | **canonical**, tracked — makes a regenerated `raw/` verifiable |
| `results/exp2_results.json` | **canonical**, tracked — machine-readable summary + environment |
| `results/exp2_table.md` | **canonical**, tracked — headline and distortion tables |
| `raw/*.bin` | regenerable (`make run`), ~640 MB, gitignored on purpose |
| `exp2_sample.exe` | regenerable (`make`), gitignored |

## Headline findings

1. **The rejected fraction is the distortion.** Conditioning on an event of probability
   `p` gives density ratio `1_box/p`, so the total-variation distance between the capped
   law and the untruncated target is *exactly* `1 − p`. Cost and error are one number.
2. **It decays only as a power law**, `λ^{−(2κ−1)}`, verified to 3 digits. At κ = 0.75 a
   cap as wide as λ = 100 still throws away 9.7% of attempts and sits 0.097 in total
   variation from the intended law. There is no cap value that makes the heavy-tail cases
   clean.
3. **A small total-variation distance does not mean a small distortion of the tail.** At
   κ = 1.5, λ = 50 the TV distance is 6.3×10⁻⁴ — indistinguishable by any probability-based
   measure — while the p99.9 speed is still 24% too small. TV bounds probabilities, not
   quantiles. Under the two-part criterion above, **no λ in the ladder is negligible for
   κ ≤ 3/2**; λ = 50 suffices at κ = 2, λ = 10 at κ = 5, λ = 5 at κ = 10.
4. **The cap breaks axisymmetry about B.** The box is a cube in normalized coordinates, so
   corner directions get √3 more radial room than axis directions. The capped law carries
   a four-fold azimuthal modulation, detected at up to −11.6σ.
5. **Acceptance is independent of θ⊥ and θ∥** — exactly, not statistically: the isotropic
   control and the anisotropic runs agree on the accept/reject decision for all 6×10⁶
   individual attempts compared.
6. **The internal attempt limit is not a practical failure mode.** Retry counts match
   `Geometric(p)`; the worst configuration in the sweep has `log₁₀ P(hit the 10⁶ limit)`
   ≈ −2.6×10⁵ per draw. The problem with the cap is the target law, not robustness.

The recommendation for the manuscript is in `results/exp2_table.md` §5: validate with the
cap **off**, and document the finite cap as an optional pragmatic finite-velocity-box
conditional target — never as a physically regularized kappa model.
