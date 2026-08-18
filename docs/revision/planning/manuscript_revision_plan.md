# JTJ1001 — manuscript revision plan, section by section

**Status:** evidence frozen 2026-08-17. All four experiments complete; literature audit closed
except the blockers in `reviewer_response_matrix.md` §9.5. **This is a plan, not the rewrite.**
No manuscript text has been changed.

**Sources of truth.** `paper/overleaf/main.tex` (current: Abstract + 7 sections, Eqs. 1–36,
Figs. 1–6, Tables I–III). Evidence: `experiments/exp{1,2,3,4}_*/results/`. Claim verdicts:
`../literature/step1_claim_audit.md`. Reviewer analysis: `reviewer_response_matrix.md`.

**Reading rule for whoever writes the prose:** every quantitative sentence carries its
parameter range, or it does not ship. Every deletion listed below is a deletion, not a
softening — the matrix records which claims were retired and why.

---

## 0. The three findings that drive this revision

Whoever writes the prose should internalize these first, because they overturn the
manuscript's existing narrative rather than refine it.

1. **Negligible rejection does not imply negligible bias.** κ=1.5, λ=50: TV distance
   6.3×10⁻⁴ while the p99.9 speed is 24% too small (Exp 2). The Abstract's inference is a
   non-sequitur. **Delete it.**
2. **We are not faster; the rejection method is.** Zenitani (2025) runs 1.8×–2.6× faster than
   our implementation where κ > 1, and per-sample cost is not constant in κ (Exp 3). **Delete
   every speed claim.**
3. **The construction is prior art, including the anisotropic extension.** ZUM 2026
   Algorithms 3.1/3.2 write out the θ∥ ≠ θ⊥ loader explicitly. **Claim no construction
   novelty, and do not claim the anisotropic case is undocumented.**

---

## 1. Abstract — **REWRITE**

| Action | Detail | Answers |
|---|---|---|
| **DELETE** | *"negligible rejection rate"* and any inference from it to negligible bias or preserved moments. Not softened — the inference is invalid. | R1.1, R1.2 |
| **DELETE** | *"runs in constant time per sample."* Refuted by measurement (88–130 ns/sample, non-monotonic in κ). | R1.4 |
| **DELETE** | Any implication of discovering the Gamma-ratio construction, and any "first"/"novel" framing. | R2.A1, R2.A2 |
| **REWRITE** | *"exact sampling algorithm"* → distributionally exact **for the uncapped target in exact arithmetic**, with prior-art equivalence stated in the same breath. | R1.1, R2.A2 |
| **REWRITE** | *"core mapping requires no rejection step"* → "the uncapped construction uses a fixed sequence of high-level variates with no outer acceptance–rejection loop", noting library Gamma generators may reject internally. | R1.1 |
| **ADD** | One clause naming the actual contribution: open-source simulation-oriented implementation, arbitrary field-frame loading, explicit uncapped/capped target-law semantics, direct 3-D validation including 1/2 < κ ≤ 3/2. | R2.A2 |
| **ADD** | The supported numerical range (double: κ ≥ 0.55; float: κ ≥ 0.75). | R1.3 |

**Do not** state a throughput number in the Abstract; it is hardware-specific. Keep it in the
performance subsection with its environment.

## 2. §I Introduction — **SUBSTANTIALLY REWRITE**

| Action | Detail | Answers |
|---|---|---|
| **DELETE** | *"prohibitively low acceptance rates"* — refuted: a dedicated rejection sampler at 0.73–0.81 acceptance beats us. | R1.4, R1.7, R2.A2 |
| **DELETE** | *"resolve these computational bottlenecks"*, *"fast and accurate"* as a contribution claim. | R1.4, R2.A2 |
| **REWRITE** | The difficulty paragraph. Lead with the difficulties that are real and now evidenced: moment non-existence for κ ≤ 3/2; heavy-tail diagnostics that fail numerically (Exp 1 `Y` vs `W`); truncation that is not a free fix (Exp 2); finite precision as κ → 1/2⁺ (Exp 4). Demote "no closed-form inverse CDF" — it is the weakest difficulty and must stop being decisive. | R2.A2 |
| **REWRITE** | Chronology, in this order: Abdul & Mace 2014 (**1-D**, Bailey polar transform — not multivariate); Abdul & Mace 2015 (isotropic trivariate scale mixture, Eq. 22); Zenitani & Nakano 2022 (**isotropic Gamma-ratio — the earliest prior art for our core**); Abdul, Matthews & Mace 2021 (2D3V bi-Kappa PIC, 6.7×10⁷ particles/species); An 2022 (Chebyshev inversion); Zenitani 2025 (Pareto rejection); ZUM 2026 (generalized (r,q), **anisotropic loader written out**, reduces to bi-Kappa at r=0, q=κ+1). | R2.A1, R1.7 |
| **ADD** | Explicit statement that this paper does **not** introduce the Gamma-ratio construction, and that independent development is chronology, not priority. | R2.A1, R2.A2 |
| **ADD** | Foundational citations R1.7 asks for: Hellinger et al. 2006 and/or Matteini et al. 2007 (anisotropy bounded by instability thresholds); Chandran et al. 2010 (turbulent perpendicular heating); Pierrard et al. 2016 for anisotropy ⊗ Kappa — **electrons only, do not generalize to ions**. Kinetic-code initialization: Birdsall & Langdon, Hockney & Eastwood. Random variate generation: Devroye 1986. | R1.7 |
| **ADD** | A bounded contribution list matching §9.3 of the matrix, with nothing above Tier 1/2. | R2.A2 |

⚠ **Blocked here:** Marsch et al. 1982 is abstract-level only — cite for the qualitative
persistence-of-anisotropy proposition, quote no numeric detail. Prefer Marsch 2018 or Hellinger
2006 if a verifiable citation is needed.

## 3. New §"Background and Related Work" — **ADD DEDICATED SECTION**

Placed after the Introduction, before the mathematical framework. Absorbs scattered §II
material so literature is discussed once.

- Standard isotropic Kappa and bi-Kappa densities, their relationship, normalization domain
  κ > 1/2, and moment-existence thresholds (second moment requires κ > 3/2).
- Parameter definitions **with physical meaning and units**: κ, θ⊥, θ∥, **B** direction,
  density/normalization, and the cap λ. State the thermal-speed convention explicitly and do
  **not** call θ a temperature.
- The isotropic limit θ⊥ = θ∥ and its relation to the standard Kappa definition.
- Sampling-method survey by class: inversion (An 2022); rejection (Zenitani 2025); scale
  mixture / Student-*t* / Beta-prime (A&M 2015, Zenitani & Nakano 2022, ZUM 2026); with the
  Beta-prime/Gamma-ratio route placed in that literature as **prior art**.
- A short paragraph distinguishing the three target laws used later: bi-Kappa, capped
  bi-Kappa, and regularized Kappa (Scherer et al. — a *different physical model*, not ours).

**Answers:** R2.A1, R2.B, R2.A3, and part of R1.7.

## 4. §II Notations and probabilistic background — **RETAIN, REORGANIZE, CONDENSE**

| Current | Action |
|---|---|
| §II.A Notations / Table I | Keep; complete the physical meanings and units (R2.B). |
| §II.B expectation/variance | Keep, condense. **ADD** the PDF definition R1.7 asks for: non-negativity + unit integral, with the measure and variables named. |
| §II.C Marginal distributions | **MOVE to Validation.** State on arrival that a projected 1-D sample is computed *from* the final 3-D vectors and is **not** an input to the sampler — this is R2.C1's whole point. |
| §II.D Transformation of random variables (Eq. 11) | Keep and **EXPAND**. Say explicitly that Eq. (11) is a derivational density rule, not an executable step (R2.C2). |
| §II.E Gamma, §II.F Beta-prime | Keep. Cite prior art at point of use. |
| §II.G Bi-Kappa | Keep; align with the new Background section's convention. |
| §II.H Q–Q plots | **MOVE to Validation.** |

## 5. §III Sampling method — **REWRITE AND SPLIT**

### 5a. Uncapped core (currently §III.A–D) — rewrite, make self-contained

- **ADD the full Jacobian derivation** R2.C3 asks for: write `V = A U`,
  `A = diag(√κ θ⊥, √κ θ⊥, √κ θ∥)`; substitute into Eq. (23) showing the quadratic form becomes
  `u_x²+u_y²+u_z²`; apply `p_U(u) = p_V(Au)|det A|` with `|det A| = κ^{3/2} θ⊥² θ∥`, showing the
  normalization factors cancel.
- **ADD the geometric clarification**, tactfully: constant-density surfaces are **ellipsoids in
  physical V-space** and **spheres in normalized U-space**. The referee's wording has these
  swapped; clarify the coordinate labels rather than reproduce the inverted statement or call
  the referee wrong (R2.C3 item 5).
- **ADD** the spherical Jacobian `r² sin θ` → `4πr²` for `p_R`, and `dr/dt = 1/(2√t)` for
  `T = R²` (R2.C2 items 3–4).
- **ADD** a compact algorithm box / pseudocode for traceability.
- **SCOPE** "exact" and "rejection-free" to this subsection only, every occurrence.
- **CITE** Zenitani & Nakano 2022 and ZUM 2026 at the point the Gamma-ratio radius is derived.
- **NOTE** the implementation forms `sqrt(x1)/sqrt(x2)`, not `sqrt(x1/x2)` — with the Exp-4
  reason (the quotient overflows where the radius is representable). This is a small,
  defensible implementation contribution.

### 5b. §III.E Truncation — **SPLIT INTO ITS OWN SUBSECTION**

This is the most-rewritten part of the paper. Six things must appear, and the current text
asserts the opposite of two of them.

1. **Define the capped target as a separate law**, once:
   `f_cap(v) = f_κ(v)·1{v ∈ B_λ} / P(v ∈ B_λ)`, with
   `B_λ = {v : |v_x|/θ⊥ ≤ λ, |v_y|/θ⊥ ≤ λ, |v_z|/θ∥ ≤ λ}`.
   **Never** call it Eq. (20), the untruncated law, or a regularized Kappa.
2. **The rejected fraction *is* the distortion.** TV(f_cap, f_κ) = 1 − P(accept) exactly.
   Quote cost and error as one number. Give the closed form for P(accept) and its
   `λ^{−(2κ−1)}` decay (verified to ~0.5% in the log-log slope).
3. **Rejection fraction does not measure tail fidelity.** State the κ=1.5, λ=50 case
   (TV 6.3×10⁻⁴, p99.9 24% low) and the pre-registered two-part criterion. Report where it is
   met — κ=10 at λ≥5, κ=5 at λ≥10, κ=2 at λ≥50 — and that **no λ in the ladder qualifies for
   κ ≤ 3/2**. At κ=0.75, λ=100 still rejects 9.7% and TV < 10⁻³ would need λ ~ 10⁹. Do not
   imply a large-enough λ exists.
4. **Finite capped variance at κ ≤ 3/2 is a property of the box, not the plasma** — 1.28 → 65.2
   at κ=0.75 as λ goes 3 → 100, growing as `λ^{3−2κ}`. One sentence warning that a naive
   diagnostic will happily print it.
5. **θ⊥ and θ∥ cancel** from the normalized accept/reject decision — give the *algebraic*
   argument (the predicate reduces to `√κ·maxᵢ|uᵢ| ≤ λ`). Mention the 6×10⁶-decision empirical
   control as redundant confirmation only; do not present it as an independent physical result.
6. **The cap breaks axisymmetry about B.** The box is a cube in normalized coordinates, so the
   conditioned law carries a four-fold azimuthal modulation (up to −11.6σ). For PIC
   initialization this is a physics defect. **This is new material and belongs in the paper.**
- **ADD** a sentence that the released API now defaults to **uncapped**, and that a finite cap
  is opt-in.
- **ADD** the Exp-4 result that in capped mode a non-finite draw can never satisfy the box
  predicate, so the loop silently redraws it — truncation **hides** the failure rather than
  solving it.
- **DELETE** the forward reference asserting negligible bias in the §IV variances.

**Answers:** R1.1, R1.2, R1.3 item 5, R1.5.
**New figure:** rejection fraction / TV vs λ, stratified by κ, with the two-part criterion
marked. **New table:** κ, λ, attempts, accepted, rejected fraction, TV, p99.9 ratio, mode.

## 6. §IV Sampling and verification → **Validation — SUBSTANTIALLY EXPAND AND RESTRUCTURE**

Replace the current "histogram + variance ⇒ correctness" structure with an explicit hierarchy.
The existing structure cannot reach κ ≤ 3/2 and does not test the central construction.

| Order | Test | Evidence |
|---|---|---|
| 1 | **Central radial law** — KS + Cramér–von Mises against the exact law, plus log-R quantile probes at p = 0.5, 0.9, 0.99, 0.999 | Exp 1, κ = 0.51 … 10 |
| 2 | **Directional uniformity** — cos Θ ~ U(−1,1), Φ ~ U(−π,π) | Exp 1 |
| 3 | **Radial–direction independence** — χ² contingency + two-sample KS between inner/outer radial quartiles | Exp 1 |
| 4 | **Anisotropy scaling** — `MAD(v∥)/MAD(v⊥)`, which needs no moments and therefore works at every κ | Exp 1, ratio 2 recovered to ≈0.4% |
| 5 | **Field-frame invariance** — draw-by-draw against the axis-aligned run at the same seed, 3 **B** directions | Exp 1, ≈1×10⁻¹⁵ relative |
| 6 | **Cartesian marginals / bi-Maxwellian limit** — kept, but demoted to a secondary check | current Figs. 2–4 |
| 7 | **Moments** — reported **only for κ > 3/2** | current Table II, regenerated |

- **STATE** mode (uncapped), λ (n/a), N, seeds, RNG, κ, θ for every dataset — this is R1.2's
  core request and it applies to the regenerated figures too.
- **ADD** the methodological result that `W = 1/(1+T) ~ Beta(κ−1/2, 3/2)` must be used and
  `Y = T/(1+T)` must not: at κ=0.55, 16.4% of `Y` rounds to exactly 1.0 and KS falsely reports
  `√n·D = 51.8`, where `W` gives 0.751. Also that `T = R²` must never be formed (overflows for
  R > 1.3×10¹⁵⁴). This is short, concrete, and genuinely useful to readers.
- **ADD** the supported-range statement from Exp 4 (double κ ≥ 0.55; float κ ≥ 0.75), and state
  that validation claims do not extend beyond it.
- **DO NOT** use variance as a correctness metric where the second moment does not exist.
- **REPLACE** "confirms the algorithm's correctness" with a bounded statement naming the tested
  distributions, parameters, diagnostics, and precision.

**Answers:** R1.3 (all six items), R1.5, R2.C1, R2.A2.

### 6a. Figures 2–4 and Tables II–III — **REGENERATE AND RELABEL**

Every R1.6 item, plus one inference the referee did not make:

- Captions gain θ⊥, θ∥, N, κ, sampling mode/cap, and unambiguous curve identities
  (histogram / KDE / analytic bi-Kappa / bi-Maxwellian reference).
- Ordinate label "Density" → **"Probability density"**, with inverse-velocity units; if
  nondimensionalized, say "dimensionless" and define the normalization.
- Abscissa gains velocity units or a normalized-axis definition.
- State that histograms are normalized probability-density estimates matched to the KDE and
  reference PDFs.
- State that the figure sequence probes convergence to the bi-Maxwellian limit as κ increases.
- **Author inference:** the plotting notebook clips histograms to a central percentile envelope
  while the KDE uses the full sample. **Either remove the clipping or disclose it** and apply
  ranges consistently — otherwise the comparison is visually biased.
- **R1.5:** state the two largest κ=2 relative errors explicitly (3.66% `vx`, 3.77% `v⊥`) in the
  validation text, attribute them to heavier-tailed finite-sample variability, and — since the
  table is being regenerated anyway — **report replicate variability / CIs** so statistical
  error is distinguishable from cutoff bias. Also check the `v⊥` definition matches the variance
  formula used.

**Answers:** R1.6 (all six), R1.5, R1.2.

## 7. New §"Performance" — **ADD, AND IT REPORTS AGAINST US**

| Action | Detail |
|---|---|
| **ADD** | Methods, environment (compiler, stdlib, flags, CPU, RNG), protocol (correctness gate before timing; 10 batches × 10⁶; warm-up discarded), and results as median with spread. |
| **ADD** | The comparative table: ours 88–130 ns/sample; A&M 2015 variant 1.4×–1.6× faster; Zenitani 2025 rejection **1.8×–2.6× faster** where κ > 1. |
| **ADD** | That per-sample cost is **not constant in κ**, and that this refutes the previous constant-time claim directly. |
| **ADD** | The two narrow favourable results: anisotropy and arbitrary-**B** rotation cost within 1–2 ns; the cap's overhead is confined to low κ (+36% at κ=0.75). |
| **ADD** | The one genuine tradeoff: Zenitani's recommended envelope index `n = κ/2` requires κ > 1, so it is undefined for `1/2 < κ ≤ 1`, where our route works. |
| **STATE** | Single machine / single toolchain; adequate to refute constant-time and fix a ≈2× ordering, **not** a cross-platform characterization. |
| **FORBID** | "fast", "resolves bottlenecks", "outperforms", "constant time" — anywhere in the paper. |

**Answers:** R1.4, R2.A2 (performance half).

## 8. §V Magnetic field-aligned frame transform — **RETAIN, EXPAND VALIDATION**

- Keep Eq. (35). **ADD** that the transform is orthogonal and therefore preserves the
  normalized radial law — which is why Exp 1's frame-invariance test is a meaningful check.
- **ADD** the Exp 1 Block B evidence (3 **B** directions, ≈1×10⁻¹⁵ relative agreement) and the
  Exp 3 result that the rotation costs nothing measurable.
- **ADD** documentation of the edge-case basis construction (the `maxcomp` branch).
- **DO NOT** present field rotation as sampler novelty. Present it as the integration
  capability prior implementations do not ship — ZUM 2026 returns field-aligned components
  `v∥, v⊥1, v⊥2` and does not rotate into an arbitrary global frame. **This is now the
  strongest single differentiator** and should be stated precisely, not expansively.

**Answers:** R2.A3, R2.C3 item 6, part of R2.A2.

## 9. §VI Code availability — **REWRITE AND EXPAND**

- Enumerate what is actually released: header-only C++ samplers, test suite, the four
  experiment directories with `make run` / `make verify`. **Python is notebook-level** — say so;
  do not claim a packaged module.
- **CORRECT** the API description: the cap is now **off by default**; a finite cap is opt-in;
  `no_cap()` and `param().capped()` report/select the target law. The previous text described
  an uncapped mode the released API did not have (R1.1 item 5) — say plainly that the API was
  fixed.
- State supported κ range per precision, the retry limit, and failure behaviour.
- State license, version, Zenodo DOI **matching the tested release**, RNG and seeding
  semantics, and the exact reproduction commands for every figure and table.
- **DO NOT** claim cross-platform bitwise reproducibility. Exp 4 shows agreement in *failure
  statistics* between libc++ and libstdc++, which is a different and weaker statement.

**Answers:** R1.1, R1.2, R2.A2.

## 10. New §"Discussion / Limitations" — **ADD**

Short, and it does most of the good-faith work in the revision:

- Independent development chronology, stated factually and **not** as priority.
- Zenitani & Nakano 2022 (isotropic) and ZUM 2026 (anisotropic, including the θ∥ ≠ θ⊥ loader)
  as prior art; the core construction is not novel here.
- Parameterization/convention differences across the Kappa literature.
- Finite-precision limits and the honest-overflow vs spurious-loss distinction; the log-domain
  mitigation is **measured and validated but not landed**, because it would change the RNG
  stream and break seed-for-seed reproducibility — an authors' decision, disclosed.
- Moment divergence for κ ≤ 3/2 as mathematics, not a bug.
- Cap tradeoffs, including the broken axisymmetry.
- Benchmark limitations (one machine, one toolchain) and the scope of validation.

**Answers:** R2.A1, R2.A2, R1.4.

## 11. §VII Summary — **REWRITE**

- **DELETE** "exact, rejection-free algorithm" as a paper-level contribution claim; "confirm the
  algorithm's correctness"; "to within a few percent" as a global characterization; constant-time
  and any speed claim.
- **REWRITE** to summarize only what is demonstrated: an open-source, tested, simulation-oriented
  bi-Kappa loader equivalent to published constructions; arbitrary field-frame support; explicit
  uncapped/capped target-law semantics with the cap quantified; direct 3-D validation including
  1/2 < κ ≤ 3/2; a documented supported numerical range; and measured performance that does not
  claim superiority.
- Bound every statement to the tested κ, θ, N, precision, and mode.

## 12. Title — **RETAIN**

Neutral as it stands. Do not add "novel"/"new". A software/implementation qualifier is optional
and only if it improves accuracy.

---

## 13. Cross-cutting audit passes (do these last, mechanically)

1. **Terminology sweep** — every occurrence of "exact", "rejection-free", "constant time",
   "fast", "prohibitive", "outperform", "confirms correctness". Each either gets scoped to the
   uncapped high-level core with its parameter range, or is deleted.
2. **Mode sweep** — every figure, table, and numerical statement declares uncapped vs capped.
3. **Moment sweep** — no variance/second-moment claim survives at κ ≤ 3/2.
4. **Citation sweep** — every new reference verified against the primary source in
   `paper/reference/README.md`; A&M 2014 never described as multivariate; the 2026 EPS paper
   cited as **Zenitani & Umeda**, not "Zenitani et al."; no claim that existing loaders were
   undemonstrated in multidimensional PIC (A&M+Matthews 2021 is 2D3V).
5. **Renumbering pass** — equations, figures, and tables will move substantially. Do the
   point-by-point response's section/line references **only after** the manuscript is stable.

## 14. What is deliberately NOT in this plan

- The rebuttal letter and the list of changes — after the manuscript is stable.
- Final figure numbering.
- Landing the log-domain Gamma mitigation — authors' decision, RNG-stream-breaking.
- Making `bi_maxwellian_distribution` symmetric with the new bi-Kappa default — an authors'
  consistency choice; the Exp-2 argument does not transfer to a Gaussian (truncation at 20
  thermal speeds removes ~10⁻⁸⁸ of the mass).
- Any new simulation application. Out of scope for a method/software paper.
