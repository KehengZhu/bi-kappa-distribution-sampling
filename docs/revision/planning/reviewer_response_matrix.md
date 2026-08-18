# JTJ1001 major-revision reviewer-response matrix

**Manuscript:** JTJ1001, “Sampling the Bi-Kappa Distribution”<br>
**Journal:** APS Open Science<br>
**Purpose:** Authoritative analysis and planning record for the major revision. This is not a rebuttal draft and does not authorize manuscript, code, or experiment changes.

# 1. Revision context

## 1.1 Source hierarchy and repository baseline

The following sources were inspected. Reviewer and editor sources take precedence over author/coauthor strategy and over this matrix's recommendations.

| ID | Source | Role and authority |
|---|---|---|
| D1 | APS Open Science email, subject `Your_manuscript JTJ1001 Zhu`, from `apsos@aps.org`, July 29, 2026 | Original major-revision decision. Contains the complete Second Referee report in the email body and the editor's resubmission instructions. |
| R1 | Attachment `jtj1001_report_1_1.pdf` to D1, 2 pages | Complete First Referee report. It contains exactly seven numbered comments. |
| E1 | Author's July 30, 2026 email to the editor | Disclosure of Zenitani, Usami, and Matsukiyo (2026), the mathematical overlap, and the proposed software/implementation/validation repositioning. This is an author proposal, not an editorial requirement. |
| E2 | Editor's August 4, 2026 reply | Confirms that the revised materials will return to the previous referees and says the authors may want to give them the July 30 context. |
| Y1 | Yuri Omelchenko correspondence, July 29–August 4, 2026 | Coauthor strategy. Useful for planning, but not a referee or editor requirement. |
| M1 | `paper/overleaf/main.tex`, `paper/overleaf/refs.bib`, and `paper/overleaf/latex-build/main.pdf` | Canonical current manuscript source and compiled manuscript in this repository. The nested manuscript Git checkout is clean; the PDF is 11 pages. |
| M0 | `APSOS__Sampling_the_Bi_Kappa_Distribution.pdf` | Slightly older 11-page submission copy. It is not used for current section/equation locations where M1 differs. |
| C1 | Current C++, Python/notebook, README, tests, and example driver in this repository | Evidence for what the companion software actually implements and validates. |

No earlier formal response matrix, rebuttal draft, or major-revision plan was found. `DISTRIBUTION_NOTES.md` is software-release planning, and `Re Kappa paper.pdf` is a pre-submission email print; neither is reviewer authority.

## 1.2 Decision status and editor instructions

- The decision is **major revision**, based on reports from **two referees**.
- The July 29 decision characterizes the concerns as substantial.
- A resubmission must include a complete point-by-point response, a marked-up manuscript PDF, and a list of changes. The decision says only the manuscript PDF is needed for the manuscript upload, not source files or separate figure files.
- On August 4 the editor wrote: “The revised materials will be sent back to the previous referees for further evaluation based on their reports. You may want to give the referees the new context you described in your correspondence of July 30.”
- Therefore, the revision and eventual response must tell the prior referees about the Zenitani overlap explicitly. It must not present the overlap merely as private editorial background.

## 1.3 Zenitani discovery and required positioning change

During revision preparation, the authors found Zenitani, Usami, and Matsukiyo (2026), *Loading Non-Maxwellian Velocity Distributions in Particle Simulations*, published in *JGR: Space Physics* in March 2026. Their generalized `(r,q)` construction reduces to the standard bi-Kappa case for `r = 0` and `q = κ + 1`. Its Beta-prime/Gamma-ratio radial draw, uniform spherical direction, and parallel/perpendicular scaling are mathematically equivalent to the manuscript's core construction.

The July 30 disclosure records independent earlier development by the present authors, with development history dating to November 2025. That chronology is relevant context, but it does **not** restore priority or mathematical novelty after publication of the Zenitani paper. The revised manuscript must not base novelty on discovering the rejection-free Gamma-ratio construction, nor imply that rejection sampling is the only existing route to Kappa samples.

The defensible working direction is narrower and software-centered:

- an open-source, simulation-oriented implementation of bi-Kappa loading;
- C++ integration and a properly supported Python path;
- arbitrary magnetic-field-frame support;
- an explicit distinction between the untruncated target distribution and an optional capped/truncated distribution;
- rigorous radial, marginal, anisotropic, and frame-transform validation;
- documented low-κ floating-point behavior;
- reproducible figures and tests; and
- evidence-based performance characterization.

This direction itself must not be overstated. Open source, magnetic-field support, numerical robustness, or performance can be claimed as contributions only to the extent that the released code and new evidence support them.

## 1.4 Current manuscript and implementation baseline

The current manuscript consists of Abstract; I Introduction; II Notations and probabilistic background; III Sampling method; IV Sampling and verification; V Magnetic field-aligned frame transform; VI Code availability; and VII Summary. It currently contains Eqs. (1)–(35), Figs. 1–6, and Tables I–III.

Important current-state facts are:

- The Abstract, Sec. VI, and Sec. VII describe the method/package as exact and rejection-free; the Abstract and Summary also claim constant time and broad correctness confirmation.
- Sec. III.E implements a component-wise box cutoff, not a radial speed cutoff, and asserts negligible rejection and variance bias without reporting `λ`, the rejected fraction, or evidence.
- Sec. IV reports `θ⊥ = 1`, `θ∥ = 2`, `κ ∈ {2,5,10}`, and `N = 100,000`, but does not identify whether truncation was used. It validates Cartesian marginals and moments, not the central radial law directly.
- The current C++ class always requires a positive `max_normalized_velocity` (default 20), rejects outside the cap, and has no untruncated/off mode. Thus the released C++ sampler is not presently rejection-free end to end, despite the manuscript's “optional” language.
- The high-level core uses a fixed sequence of Gamma and uniform draws, but `std::gamma_distribution` may use a variable-cost internal algorithm. A strict worst-case “constant time” statement is not justified.
- The C++ implementation does provide header-only samplers, deterministic/external RNG support, and arbitrary magnetic-field rotation. Existing tests cover API behavior, caps, seeding, and frame-transform algebra, but not the radial law, full bi-Kappa distribution, low-κ robustness, moments, or speed.
- The Python bi-Kappa code is currently notebook-level rather than a packaged module. The committed C++ example uses `κ = 2`, `N = 200,000`, and cap 100, which does not reproduce the manuscript's full `κ` sweep or stated `N`. The provenance of Tables II–III and Figs. 2–6 is not self-documenting from one committed driver.

## 1.5 Coauthor strategy, kept separate from reviewer requirements

Yuri's August 4 recommendation is to keep the paper focused on method/software rather than add a new simulation application; retain theory as necessary explanation; treat existing and new tests as software validation; explain independent development and the recent Zenitani publication; consider implementing and comparing the Zenitani piecewise-rejection approach; and make the software openly usable. July 29–30 messages similarly suggest C++/Python examples, magnetic-field support, and investigation of Zenitani's “zero division” discussion.

These are **author/coauthor strategy choices**, not demands made by either referee. Earlier author-side statements that no standard-Kappa samplers existed were superseded by the Zenitani discovery and must not be used as factual premises.

# 2. First Referee — detailed response matrix

The attachment contains exactly seven numbered comments. The wording below is retained from `jtj1001_report_1_1.pdf`.

### R1.1 — Separate exact untruncated sampling from cutoff sampling

**Reviewer comment**

> The abstract and Summary describe an exact, rejection-free sampler of the Bi-Kappa distribution, while Sec. III.E adds an optional velocity cutoff that throws away extreme speeds and redraws. That produces a truncated distribution, not the full formula in Eq. (20). Sec. VI calls the core code “exact and rejection-free” based on Sec. III, but Sec. III includes both the core method and truncation. The authors should clearly separate these two cases: with truncation off, samples follow Eq. (20); with truncation on, samples follow a cut-off version of Bi-Kappa. The term “rejection-free” should apply only to the core steps in Sec. III.A–D, not to the optional cutoff step.

**What the referee is actually asking**

This is a **correctness and terminology issue**, not cosmetic wording. Two different target probability laws are being conflated: the full Eq. (20) distribution and the conditional distribution obtained after rejecting samples outside a bounded region. “Exact” and “rejection-free” must be scoped to the untruncated high-level construction only.

**Relevant current manuscript locations**

- Abstract: “exact sampling algorithm whose core mapping requires no rejection step,” followed by optional truncation.
- Sec. III.A–D, especially Eqs. (23)–(29): untruncated construction.
- Sec. III.E: component-wise box cutoff in `|vx|/θ⊥`, `|vy|/θ⊥`, and `|v∥|/θ∥`.
- Sec. VI: package described as implementing the “exact, rejection-free algorithm of Sec. III.”
- Sec. VII: “exact, rejection-free algorithm” and “no rejection step in the bi-Kappa core.”
- Current C++ API/README: positive cap is mandatory and defaults to 20; there is no off mode.

**Required manuscript changes**

1. Split the method presentation into an explicitly untruncated core and a separately named optional capped sampler.
2. State each target law. For the cap, identify the support geometry as a component-wise box and define the retained law as Eq. (20) conditioned on the cap event, with its own normalization.
3. Restrict “distributionally exact” to the untruncated target in exact arithmetic. Restrict “rejection-free” to the manuscript-level Gamma-ratio/direction/scaling construction, and acknowledge that library Gamma generators may have internal variable-cost/rejection behavior.
4. Make Abstract, Sec. III, Sec. VI, and Summary use the same terminology.
5. Before the revised software can be described as offering both modes, add and document an actual no-cap/off mode in the C++ API. Otherwise describe the released C++ behavior honestly as capped.
6. Do not say truncation “enforces finite moments” without making clear that it changes the target distribution.

**New analysis / experiment required**

**New analysis required.** Formalize the two target laws and reconcile the paper with the actual public API. Quantitative cutoff evidence is consolidated under Experiment 2 in Sec. 5.

**Evidence now available**

**Experiment 2 answers this comment and is complete.** The referee was right on both counts:
the capped mode is a different target law, and it is not rejection-free.

| Referee's point | Evidence | Where |
|---|---|---|
| Capped sampling is a *different* law, not Eq. (20) | The capped output is **bitwise identical** to the uncapped draws restricted to the box, 240/240 (case, λ, seed) pairs. "Capped = uncapped conditioned on the box" is a measured fact about the shipped binary, not a claim about the code. | `subsequence_check` |
| The cap is not rejection-free | Rejected fraction tabulated over κ × λ with a matching closed form; retry counts match `Geometric(p)` | `exp2_table.md` §1–2 |
| Truncation changes the target, and cannot be waved off | TV distance from the target is *exactly* `1 − P(accept)` — cost and distortion are one number, decaying only as `λ^−(2κ−1)` | `tail_exponent_scaling` |
| Do not say truncation "enforces finite moments" innocently | Capped variance at κ ≤ 3/2 is finite but is a property of λ, not the plasma: 1.28 → 65.2 at κ=0.75 as λ 3 → 100 | `exp2_table.md` §3 |

Two results strengthen the referee's position beyond what was asked:

- **The cap breaks axisymmetry about B.** The box is a cube in normalized coordinates, so the
  conditioned law carries a four-fold azimuthal modulation, detected at up to −11.6σ. For PIC
  initialization a non-gyrotropic initial condition is a physics defect, not merely a
  statistical one.
- **The θ's cancel exactly** from the accept/reject decision, so the capped law in the local
  frame depends on `(κ, λ)` only. This is an algebraic identity, confirmed by 0 disagreements
  in 6×10⁶ attempts — the empirical control is redundant confirmation, not an independent
  physical result, and must be presented that way.

**API status:** `no_cap()` and `capped()` now exist, with tests `test_no_cap_semantics` B1–B4,
so item 5 ("add and document an actual no-cap/off mode") is **discharged in code**. The
remaining engineering question — whether uncapped should become the *default* — is tracked in
§9.6.

**Remaining blocker**

**None for evidence.** Items 1–4 and 6 are manuscript prose; item 5 is landed in code.

**Recommended response strategy**

Agree fully. Thank the referee for identifying a real scope inconsistency, revise terminology globally, and explain the code/API correction. Do not defend the current wording. Concede specifically that the released C++ API had no uncapped mode when the manuscript claimed one, and state that it now does.

**Dependency / overlap**

Directly overlaps R1.2 and R1.3; also affects R1.4, R2.A2, the Zenitani matrix, and the software/reproducibility revision.

**Priority:** BLOCKER<br>
**Status:** **EVIDENCE COMPLETE** (was: REQUIRES NEW WORK). Prose not yet written.

### R1.2 — Disclose truncation settings and validation provenance

**Reviewer comment**

> Sec. III.E states that a large cutoff λ makes the discarded fraction, and thus any bias in the Sec. IV variances, negligible, while Table II and the Summary compare sample variances with the untruncated expressions in Eqs. (31)–(33). Sec. IV does not state whether truncation was used in those runs, nor does it give λ or the discarded fraction. Without that information, it is unclear whether Table II tests Eq. (20) or the truncated distribution of Sec. III.E, and whether the “negligible bias” claim is shown or only asserted. The authors should report these settings in Sec. IV (and in Table II if truncation was used). If truncation was off, they should say so and revise the forward reference in Sec. III.E accordingly.

**What the referee is actually asking**

This is a **reproducibility and unsupported-claim issue**. The validation target and run configuration cannot be reconstructed, so the reported agreement with untruncated theory is not interpretable.

**Relevant current manuscript locations**

- Sec. III.E: unquantified “discarded fraction” and “negligible” variance bias.
- Sec. IV opening: `θ⊥ = 1`, `θ∥ = 2`, `κ ∈ {2,5,10}`, `N = 10^5`, but no cutoff state, `λ`, seed, or rejection count.
- Eqs. (31)–(33) and Table II: untruncated moments.
- Sec. VII: variance agreement claim.
- Current repository: C++ driver now says `N = 200,000`, `κ = 2`, cap 100; this does not establish the provenance of the manuscript's multi-κ results.

**Required manuscript changes**

1. State, for every validation dataset, whether the untruncated or capped sampler generated it.
2. If capped, report `λ`, the cap geometry, attempted and accepted draws, rejected fraction, seed/RNG, sample size, and whether theoretical values are for the capped or uncapped law.
3. Do not compare capped samples to Eqs. (31)–(33) without either deriving/estimating capped moments and uncertainty or quantitatively showing that cap bias is below the reported statistical error.
4. Put the settings in Sec. IV and in the Table II caption or a compact run-configuration table.
5. If uncapped, state that explicitly and remove/rewrite the Sec. III.E forward claim about bias in those results.
6. Commit one reproducible driver/configuration that regenerates the table and associated figures.

**New analysis / experiment required**

**New numerical experiment required.** Existing provenance is insufficient to determine the original cutoff state reliably; regenerate the validation with explicit mode/settings and measure rejection when capped.

**Evidence now available**

**Experiment 2 answers this comment and is complete — and the answer is stronger than the
referee asked for.** The referee asked us to *report* λ and the discarded fraction and to
*support* the negligible-bias claim. The first is now a full table; the second turned out to be
unsupportable in principle.

| Referee's request | Evidence | Where |
|---|---|---|
| Report `λ` and the discarded fraction | Full κ × λ table, λ ladder `{3, 5, 10, 20, 50, 100}` plus uncapped, with a closed-form `P(accept)` cross-check | `exp2_table.md` §1 |
| State mode for every dataset | `raw/manifest.csv` carries `mode` and `lambda` per run; every reported number is keyed to it | `raw/manifest.csv` |
| Report `N`, seeds, RNG, attempts/accepted | 280 runs, 5 seeds × 10⁵, `std::mt19937` seeded per run, all in the manifest | `raw/manifest.csv` |
| Support or drop "negligible bias" | **The inference is invalid.** At κ=1.5, λ=50: TV = 6.3×10⁻⁴ yet p99.9 speed is **24% too small** | `negligibility` |
| Commit one reproducible driver | `make run` + `make verify`; bit-for-bit deterministic against committed checksums | `GNUmakefile` |

**The decisive finding.** A small discarded fraction bounds *probabilities*; it says nothing
about *quantiles*, which is precisely what a Kappa distribution exists to represent. The
Abstract's negligible-rejection argument is therefore a **non-sequitur and must be deleted, not
softened**. On the pre-registered two-part criterion (`TV < 10⁻³` **and** p99.9 within 1%),
negligibility holds at κ=10 for λ≥5, κ=5 for λ≥10, κ=2 for λ≥50, and **for no λ in the ladder
at κ ≤ 3/2**. At κ=0.75 the decay exponent is `2κ−1 = 0.5`, so λ=100 still rejects 9.7% and
`TV < 10⁻³` would need λ ~ 10⁹ — so the manuscript must not imply that a large-enough λ exists.

Note for the response: the library default λ=20 and the committed example's λ=100 are **both**
non-negligible everywhere in κ ≤ 3/2. See §9.6.

**Remaining blocker**

**None for evidence.** What remains is manuscript prose plus regenerating Table II / Figs. 2–4
from an auditable driver with the mode stated (shared with R1.6).

**Recommended response strategy**

Agree. State that the original manuscript omitted essential settings. Supply auditable regenerated results rather than infer the old settings from today's example code. Then go further than the referee asked: report that the negligibility *argument* is invalid and has been removed, not merely quantified.

**Dependency / overlap**

Depends on R1.1's two-mode definition and is largely answered by Experiment 2. It also overlaps R1.3, R1.5, R1.6, and reproducibility requirements.

**Priority:** HIGH<br>
**Status:** **EVIDENCE COMPLETE** (was: REQUIRES NEW WORK). Prose not yet written.

### R1.3 — Expand validation to the radial law, low κ, and one truncated case

**Reviewer comment**

> Sec. IV tests only ((\theta_\perp, \theta_\parallel) = (1, 2)) and (\kappa \in \{2, 5, 10\}), using marginals, variances, and Q–Q plots. The Summary’s claim that this “confirms the algorithm’s correctness” overstates those results. No test is given for (1/2 < \kappa \le 3/2) or for truncation with a reported rejection rate, although Sec. III.E discusses both. The method is based on (T = R^2 \sim \beta'(3/2, \kappa - 1/2)), but this radial law is not checked directly. The authors should limit the Summary to the cases actually shown, either verify (T) (or (R)) against Eq. (27) or avoid stronger joint claims, and report one truncated example with (\lambda) and the rejection fraction if truncation is part of the method.

**What the referee is actually asking**

This is a **validation gap and overclaim**. Matching one-dimensional Cartesian marginals does not directly test the central Beta-prime radial construction or fully establish the three-dimensional joint law. The manuscript also discusses a domain (`1/2 < κ ≤ 3/2`) that its moment-based validation cannot cover.

**Relevant current manuscript locations**

- Sec. II.H says `T` or `R` can be checked, but the reported Q–Q plots test only Cartesian components.
- Sec. III.B, Eqs. (26)–(28): `T = R² ~ BetaPrime(3/2, κ − 1/2)`.
- Sec. III.E: low-κ moment divergence and truncation.
- Sec. IV.A–C: one anisotropy pair, `κ = 2,5,10`; no direct radial or low-κ test.
- Sec. VII: “confirm the algorithm's correctness.”

**Required manuscript changes**

1. Add direct untruncated radial-law validation over a range that includes `1/2 < κ ≤ 3/2`.
2. Test the radial law through a numerically bounded transform so extreme tails do not make the diagnostic itself unstable — but use `W = 1/(1+T) ~ Beta(κ−1/2, 3/2)`, **not** `Y = T/(1+T) ~ Beta(3/2, κ−1/2)`. Experiment 1 measured both: at κ = 0.55, 16.4% of `Y` values round to exactly 1.0 and the KS test falsely reports `√n·D = 51.8`, while `W` on the same data gives 0.751. Never form `T = R²` either; it overflows for `R > 1.3e154` where `R` is still representable.
3. Test directional uniformity and radial–direction independence, or otherwise moderate any claim about the full joint law.
4. Do not use variance as a correctness metric where the theoretical variance does not exist (`κ ≤ 3/2`). Use distributional distances and quantiles there.
5. Include at least one fully specified capped example with `λ` and rejection fraction if the capped mode remains part of the paper.
6. Replace global “confirms correctness” language with a bounded statement naming the tested distributions, parameters, diagnostics, and numerical precision.

**New analysis / experiment required**

**New numerical experiment required.** See Experiments 1, 2, and 4 in Sec. 5.

**Evidence now available**

**Experiment 1 answers this comment directly and is complete.**
`experiments/exp1_radial_directional/` — 135 runs, 1.35×10⁷ draws, released C++ sampler,
uncapped, κ from 0.51 to 10, 5 seeds × 10⁵ per configuration.

| Referee's request | Evidence | Where |
|---|---|---|
| Verify `T` (or `R`) against Eq. (27) directly | Radial-law KS + Cramér–von Mises + log-R quantile probes at p = 0.5, 0.9, 0.99, 0.999 — **PASS** | `results/exp1_table.md`, `exp1_results.json` |
| Test `1/2 < κ ≤ 3/2` | κ = 0.51, 0.55, 0.75, 1.0, 1.25, 1.5 all covered and passing | Block A |
| Do not use variance where it does not exist | Anisotropy measured by `MAD(v∥)/MAD(v⊥)`, which needs no moments; moments reported only for κ > 3/2 | Block B |
| Joint law, not just marginals | Directional uniformity + radial–direction independence (χ² contingency and two-sample KS) — **PASS** | Block A |
| Frame handling | Arbitrary **B**-frame invariance, draw-by-draw against the axis-aligned run at the same seed — **PASS** | Block B |
| Truncated example with `λ` and rejection fraction | **Experiment 2 — COMPLETE.** Full κ × λ rejection table with a matching closed form; e.g. κ=2, λ=20 → 3.6×10⁻⁴ rejected | `experiments/exp2_cap_characterization/` |

Experiment 4 additionally bounds the numerical range over which any of this can be claimed:
double is clean for κ ≥ 0.55, float for κ ≥ 0.75. The Summary's scope statement must not
exceed those bounds.

**Remaining blocker**

**None.** All six items are discharged. Item 5 was closed by Experiment 2.

**Recommended response strategy**

Agree and provide new evidence. Explicitly acknowledge that Cartesian marginals alone were insufficient for the central construction and that moments cannot validate the low-κ regime.

**Dependency / overlap**

Overlaps R1.1–R1.2, R1.5, R2.A2, R2.C1, and the low-κ robustness work.

**Priority:** BLOCKER<br>
**Status:** **EVIDENCE COMPLETE** (was: REQUIRES NEW WORK). Prose not yet written.

### R1.4 — Support or moderate speed and constant-cost claims

**Reviewer comment**

> The Introduction states that naive rejection sampling has prohibitively low acceptance rates and that the present method is a fast algorithm that resolves this bottleneck; the abstract adds that sampling runs in constant time per sample. Constant cost is plausible from the core algorithm, but the manuscript reports no timings, acceptance rates, or comparisons to support those claims. The authors should either add minimal performance data or moderate the Introduction so that speed and bottleneck claims are not made without evidence.

**What the referee is actually asking**

This is an **unsupported performance claim**. Algorithmic intuition is not evidence of throughput, acceptance, scaling, or comparative advantage. “Constant time” is also ambiguous because both the capped outer loop and the library Gamma generator can take variable work.

**Relevant current manuscript locations**

- Abstract: “runs in constant time per sample.”
- Introduction: “prohibitively low acceptance rates,” “fast and accurate,” and “resolve these computational bottlenecks.”
- Sec. III.E: rejection loop for capped samples.
- Sec. VI and VII: software-performance implications and repeated constant-time claim.
- No current performance figure, table, benchmark script, or acceptance-rate result exists.

**Required manuscript changes**

1. Remove “prohibitively,” “resolves the bottleneck,” “outperforms,” and strict “constant time” unless quantitative evidence supports each statement.
2. If performance remains a contribution, add a reproducible benchmark against a relevant alternative. The equivalent Zenitani Gamma-ratio construction is prior art, not an independent mathematical competitor; a fair comparison may instead evaluate implementation variants and Zenitani's piecewise-rejection approach.
3. Report absolute throughput/latency, κ dependence, cap dependence, acceptance/retry behavior, failures, implementation language, RNG, compiler, hardware, and uncertainty across repeats.
4. Describe the core conservatively as using a fixed number of *high-level* random-variate calls; do not imply fixed worst-case primitive cost.
5. Report unfavorable or neutral outcomes without marketing reinterpretation.

**New analysis / experiment required**

**Benchmark required** if any speed, bottleneck-resolution, or outperformance claim is retained. Otherwise those claims must be retired.

**Evidence now available**

**Experiment 3 is complete, and it settles this comment against us.**
`experiments/exp3_benchmark/` — correctness gate then timing, 10 batches × 10⁶ draws per
configuration, `std::mt19937` for every method, Apple clang 21 / libc++ / arm64.

**The released implementation is the slowest of the three methods tested.**

| κ | released Gamma-ratio | A&M 2015 normal-triple | Zenitani 2025 Pareto rejection |
|---|---|---|---|
| 1.5 | 88.1 ns | 46.6 ns | **48.9 ns** (acc. 0.806) |
| 2 | 129.9 ns | 88.3 ns | **49.9 ns** (acc. 0.786) |
| 5 | 121.0 ns | 80.5 ns | **53.1 ns** (acc. 0.751) |
| 10 | 124.2 ns | 82.3 ns | **53.7 ns** (acc. 0.740) |
| 50 | 117.4 ns | 74.8 ns | **54.1 ns** (acc. 0.733) |

- Zenitani (2025) rejection costs **0.38×–0.56×** ours over κ ∈ [1.5, 50] — i.e. **1.8×–2.6×
  faster**, *despite* being a rejection method. The literature's implied expectation was
  correct and ours was wrong.
- A&M 2015's scale mixture is **≈1.4×–1.6× faster** than ours — and it is the *same
  construction*, differing only in buying the direction from three normals rather than two
  uniforms. So part of our cost is our own implementation choice, not the algorithm.
- **"Constant time per sample" is contradicted by direct measurement:** the baseline ranges
  88–130 ns/sample non-monotonically in κ, fastest at κ=1.5 and slowest at κ=2. (`shape =
  κ−1/2` crosses 1 at κ=1.5, where library Gamma generators switch algorithm — plausible
  mechanism, not instrumented.)

**Two favourable results, both narrow:** anisotropy and arbitrary-**B** rotation cost within
≈1–2 ns of the isotropic baseline — C2 is free. And the cap's overhead is confined to low κ
(+36% at κ=0.75, nil for κ ≥ 1.5), tracking Exp 2's rejection rates exactly.

**Faithfulness of the transcription is itself evidenced:** measured Pareto acceptance
reproduces Zenitani's published 0.806 / 0.785 / 0.750 at κ = 1.5 / 2 / 5 to three digits, on
a quantity we did not fit.

**Bounded scope, stated so it is not overread:** `n = κ/2` requires κ > 1, so the rejection
method is **inapplicable** at κ = 0.75 and 1.0 — recorded as a result, not skipped. Single
machine, single toolchain; adequate to refute constant-time and establish a ≈2× ordering,
**not** a cross-platform characterization.

**Remaining blocker**

**None.** The comment is answerable now — by retirement, not support.

**Recommended response strategy**

Agree, and **retire the claims rather than defend them.** Tell the referee plainly that we benchmarked, that the dedicated rejection method of Zenitani (2025) is roughly twice as fast as ours where it applies, that our own per-sample cost is not constant in κ, and that every speed and bottleneck-resolution claim has been deleted accordingly. Retain only the measured absolute throughput and the two narrow favourable results. Reporting an unfavourable benchmark is stronger evidence of good faith than omitting it, and the referee explicitly offered "or moderate the Introduction" as an acceptable outcome.

**Dependency / overlap**

Overlaps R2.A2 and the Zenitani novelty-impact matrix; depends on R1.1's mode definitions.

**Priority:** HIGH<br>
**Status:** **EVIDENCE COMPLETE** (was: REQUIRES NEW WORK). Prose not yet written; the required prose is deletion plus a short measured-performance paragraph.

### R1.5 — Discuss the larger κ = 2 variance errors

**Reviewer comment**

> Table II shows relative errors of about 3.7% for some (\kappa = 2) cases, and much smaller errors for (\kappa = 5) and (10). The Summary describes all results as agreeing “to within a few percent,” which is broadly consistent with the table, but the larger errors at (\kappa = 2) should be noted briefly in Sec. IVB.

**What the referee is actually asking**

This is a **presentation and statistical-interpretation issue**. The largest reported discrepancies (3.66% for `vx` and 3.77% for `v⊥`) should not be hidden by a broad summary phrase. Heavy tails can increase finite-sample variance-estimator variability, but that explanation must be demonstrated or stated cautiously rather than assumed.

**Relevant current manuscript locations**

- Sec. IV.B, Eqs. (31)–(33), and Table II.
- Sec. VII: “to within a few percent.”

**Required manuscript changes**

1. State the two largest `κ = 2` relative errors explicitly in Sec. IV.B.
2. Explain that lower κ gives heavier tails and potentially larger finite-sample variability, without labeling the observed deviations as mere noise unless replicated uncertainty supports that conclusion.
3. If Table II is regenerated, report replicate variability or confidence intervals and distinguish statistical error from cutoff bias.
4. Check the definition of `v⊥` used in Table II; its theoretical/sample quantity must be unambiguous and consistent with the stated variance formula.

**New analysis / experiment required**

**Textual clarification only** for the referee's minimum request. Replicated uncertainty becomes part of Experiment 1 if stronger interpretation is retained.

**Recommended response strategy**

Agree and add a brief, numerically specific discussion. Do not claim the discrepancy proves either bias or correctness without uncertainty analysis.

**Dependency / overlap**

Overlaps R1.2–R1.3 and the regenerated validation table.

**Priority:** MEDIUM<br>
**Status:** UNADDRESSED

### R1.6 — Make Figs. 2–4 self-contained and quantitatively labeled

**Reviewer comment**

> Figures 2–4 compare normalized histograms, KDE fits, and bi-Maxwellian reference curves [Eq. (30)] for Vx, Vy, and Vz, but the captions list only κ. Please add information on (θ⊥,θ∥) and N, a brief identification of each curve, and note that the plots test convergence to the bi-Maxwellian limit as κ increases. Replace the y-axis label “Density” with “Probability density,” give units on both axes, and state that the histograms are normalized to match the KDE and reference PDFs.

**What the referee is actually asking**

This is a **figure documentation and units issue**. The plots are not self-contained enough to reproduce or interpret, and “Density” does not identify a probability density or its units.

**Relevant current manuscript locations**

- Sec. IV.A and Eq. (30).
- Figs. 2–4 and their captions.
- Plotting notebook: histogram display is trimmed to a central percentile envelope while KDE calculations use the full sample; axes currently say “Density” and omit units.

**Required manuscript changes**

1. Add `θ⊥`, `θ∥`, `N`, κ, sampling mode/cap, and curve identities to every caption.
2. State that the histograms are normalized probability-density estimates and identify histogram, KDE, analytic bi-Kappa (if shown), and bi-Maxwellian curves unambiguously.
3. Label the ordinate “Probability density” and give inverse-velocity units. If variables are nondimensionalized, say “dimensionless” and define the normalization rather than invent physical units.
4. Give velocity units or normalized-axis definitions on each abscissa.
5. State that the figure sequence probes convergence toward the bi-Maxwellian limit as κ increases.
6. **Author inference beyond the explicit comment:** either remove the notebook's percentile clipping or disclose it and apply plotting ranges consistently so histogram/KDE comparisons are not visually biased.

**New analysis / experiment required**

**Existing results can be reused** only after provenance and cutoff status are established. Figure regeneration is required, but no distinct new scientific experiment is needed solely for this comment.

**Recommended response strategy**

Agree and regenerate the figures/captions from a committed, auditable plotting path.

**Dependency / overlap**

Overlaps R1.2's run settings and the reproducibility work; does not replace R1.3's direct radial validation.

**Priority:** MEDIUM<br>
**Status:** UNADDRESSED

### R1.7 — Add foundational, sampling, and PDF references/definitions

**Reviewer comment**

> The Introduction would benefit from adding a few more original references. The second paragraph cites [4, 5] for solar-wind temperature anisotropy and turbulence-driven perpendicular heating; please add one or two foundational references on these topics. The third paragraph states that Maxwellian generators are widely available and that naive rejection sampling can have very low acceptance rates, but gives no source. Please cite representative kinetic-code initialization papers and a standard reference on random variate generation/rejection sampling. Please add a brief definition of a PDF (including normalization) or cite a standard reference at the start of the section.

**What the referee is actually asking**

This is a **citation, literature, and accessibility issue**. Several broad historical and computational statements lack primary or standard sources, and the probability-density notation assumes background not provided to general readers.

**Relevant current manuscript locations**

- Introduction paragraphs 2–3.
- Start of Sec. II and Table I.
- `refs.bib` already contains some potentially relevant kinetic-code and random-variate sources (e.g., Birdsall–Langdon, Hockney–Eastwood, Verboncoeur, Devroye), but they are not cited in the relevant text.
- Existing citations on anisotropy/heating do not by themselves fulfill the request for one or two foundational original references.

**Required manuscript changes**

1. Add and verify one or two foundational primary references for solar-wind temperature anisotropy and turbulence/cyclotron-related perpendicular heating.
2. Cite representative particle/kinetic-code initialization literature for commonly used Maxwellian loaders.
3. Cite a standard random-variate/rejection-sampling reference; use it to state when acceptance becomes poor rather than making an unbounded “prohibitively low” assertion.
4. At the start of the probability background, define a PDF by non-negativity and unit integral (or cite a standard source), and state the measure/variables used.
5. Verify every added bibliographic field and DOI against a primary source before submission.

**New analysis / experiment required**

**Textual clarification only.** This requires literature verification, not a numerical experiment.

**Recommended response strategy**

Agree. Add the requested primary/standard citations and a concise PDF definition. Do not treat currently uncited BibTeX entries as though the comment were already addressed.

**Dependency / overlap**

Overlaps R2.A1–A2 and R2.B; the new Related Work/Background section should prevent duplicate literature discussion.

**Priority:** HIGH<br>
**Status:** UNADDRESSED

# 3. Second Referee — detailed response matrix

The report appears in the July 29 decision email. Its preamble is:

> This paper proposes an exact algorithm for sampling the bi-kappa distribution, which is highly relevant to space science. However, the manuscript requires significant improvements before a final decision can be made. The authors need address the following concerns and questions:

The referee then gives the following headings and seven items. There are no additional Second Referee items in the retrieved decision correspondence.

> A. Introduction and Motivation. While the introduction establishes why the bi-kappa distribution is important in space plasma physics, the following points need to be addressed:

### R2.A1 — Review the state of Kappa-sampling literature

**Reviewer comment**

> What is the current state of progress in kappa distribution sampling? Given its importance in space science, a comprehensive overview of existing literature is needed.

**What the referee is actually asking**

This is a **state-of-the-literature and completeness issue**. The manuscript needs a current, structured account of what has already been sampled, by which mathematical constructions and software implementations, and where the present package fits. The Zenitani discovery makes this a blocker rather than a routine citation addition.

**Relevant current manuscript locations**

- Introduction: frames direct Kappa sampling as a largely unresolved computational challenge.
- Sec. II: mathematical background, but no review of sampling methods.
- Sec. III: presents the Beta-prime/Gamma-ratio construction without current prior-art context.
- `refs.bib` contains Abdul & Mace and several Zenitani entries, but the manuscript does not cite them; no Student-t literature is currently present.

**Required manuscript changes**

1. Add a dedicated Related Work/Background subsection covering isotropic Kappa sampling, the equivalence to multivariate Student-t/scale-mixture constructions where applicable, rejection/inversion approaches, and simulation-loading literature.
2. Cite and explain Zenitani et al. (2026), including the `r = 0`, `q = κ + 1` reduction and mathematical equivalence to the present core.
3. Distinguish mathematical methods from available, simulation-ready implementations; do not equate lack of a convenient package with lack of an algorithm.
4. Identify which capabilities are pre-existing, which are independently reimplemented, and which repository features are actually additional.
5. Verify primary sources and avoid claims of comprehensiveness until the search is complete.

**New analysis / experiment required**

**New analysis required.** A systematic literature analysis and equivalence discussion are required, but no numerical experiment is necessary for A1 itself.

**Recommended response strategy**

Agree fully. Explain that the revision now includes literature discovered after submission and corrects the original framing. Independence chronology can be disclosed, but not used to erase prior publication.

**Dependency / overlap**

Overlaps R1.7, R2.A2–A3, R2.B, and every Zenitani-affected claim.

**Priority:** BLOCKER<br>
**Status:** UNADDRESSED

### R2.A2 — Explain difficulty, solved problems, novelty, and evidence of performance

**Reviewer comment**

> Using language accessible to general readers who may lack specialized knowledge in computational algorithms and the kappa sampling problem, please explain why kappa sampling is difficult. What specific mathematical or physical properties cause these challenges? Which of these problems have already been solved by existing methods, and which ones does this paper aim to solve? Furthermore, please clarify what is novel about the proposed method and explain how it outperforms current approaches.

**What the referee is actually asking**

This combines **accessibility, literature, novelty, scope, and comparative-evidence issues**. The revision must separate real mathematical/numerical difficulties from problems already solved by standard distributions or prior algorithms. It must then define a truthful residual contribution. Any “outperforms” claim requires a fair benchmark; if no advantage is demonstrated, the claim must be removed rather than argued rhetorically.

**Relevant current manuscript locations**

- Abstract: no inverse CDF, exact/rejection-free, constant-time, and correctness claims.
- Introduction paragraph 3: “intractable” inversion, “prohibitively low” rejection acceptance, and “fast and accurate” solution.
- Sec. III: Gamma-ratio derivation now known to be equivalent to Zenitani et al. (2026).
- Sec. III.E: heavy tails, finite-moment thresholds, and cutoff.
- Secs. VI–VII: implementation and contribution claims.
- No current comparative benchmark or low-κ numerical study exists.

**Required manuscript changes**

1. Explain in non-specialist terms that heavy tails, κ-dependent moment existence, anisotropic parameterization, field-aligned coordinates, and finite-precision random-variate behavior can complicate practical loading.
2. Also explain that absence of an elementary inverse CDF does not imply absence of exact sampling: the radial law is a standard Beta-prime/Student-t-related construction, and existing work already exploits equivalent representations.
3. State explicitly that the paper does not newly discover the Gamma-ratio sampler.
4. Recast contributions around verified implementation, API/integration, frame handling, reproducibility, validation, cutoff semantics, robustness, and measured performance—only where the revised repository supports them.
5. Define the comparison baseline and benchmark it if claiming outperformance. If results show parity, tradeoffs, or inferiority, report that result and remove superiority language.
6. Avoid treating a lack of a closed-form inverse CDF as the sole or decisive reason the problem is “difficult.”

**New analysis / experiment required**

**Benchmark required** for the explicit outperformance part of the comment. Literature/novelty rewriting is required regardless; Experiments 1–4 support any retained implementation claims.

**Evidence now available**

The comment has two halves and they now have different answers.

*"How does it outperform current approaches?"* — **It does not, and Experiment 3 measures that.**
The released implementation is the slowest of the three methods benchmarked: Zenitani (2025)
rejection is 1.8×–2.6× faster where κ > 1, and Abdul & Mace (2015)'s variant of our own
construction is 1.4×–1.6× faster. The answer to the referee is a direct, measured "it does
not outperform them", not a rhetorical reframe. This is the strongest available response
because it is falsifiable and unfavourable.

*"What is novel?"* — the construction is not. Prior art is Zenitani & Nakano (2022) for the
isotropic Gamma-ratio route, A&M (2015) Eq. (22) for the trivariate scale mixture, and
ZUM (2026) Algorithms 3.1/3.2 for the anisotropic loader **including** the θ∥ ≠ θ⊥ scaling
(see the §9.4 withdrawal — do not claim the anisotropic extension is undocumented). What
survives is implementation, arbitrary-**B**-frame loading, target-law semantics, validation
depth, reproducibility, and finite-precision characterization — the ledger in §9.3.

*"Why is Kappa sampling difficult?"* — the honest accessible answer is now evidence-backed:
heavy tails mean moments may not exist (κ ≤ 3/2); a naive bounded diagnostic silently fails
(Exp 1's `Y` vs `W`); truncation is not a free fix (Exp 2); and finite precision bites as
κ → 1/2⁺ (Exp 4). None of these is "no closed-form inverse CDF", which is the weakest of the
difficulties and must stop being presented as the decisive one.

**Remaining blocker**

**None for evidence.** Items 1–6 are prose, drawing on Exps 1–4 and the literature audit.

**Recommended response strategy**

Agree and substantially reframe. Tell the referee directly that the newly identified Zenitani work changes the novelty claim, and that the benchmark shows we do **not** outperform — quoting the numbers. Do not attempt to recover superiority on a different axis; state the tradeoff (our route covers κ ≤ 1, where the recommended rejection envelope is undefined) and stop there.

**Dependency / overlap**

Overlaps R1.3–R1.4, R1.7, R2.A1/A3/B, and the entire novelty-impact matrix.

**Priority:** BLOCKER<br>
**Status:** **EVIDENCE COMPLETE** (was: REQUIRES NEW WORK). Prose not yet written; this is the largest single writing task in the revision.

### R2.A3 — Distinguish isotropic-Kappa and bi-Kappa sampling challenges

**Reviewer comment**

> What unique challenges does bi-kappa sampling present compared to isotropic kappa sampling?

**What the referee is actually asking**

This is a **motivation and scope-clarity issue**. The manuscript must explain whether anisotropy changes the stochastic core or mainly adds parameterization, coordinate transforms, implementation, and validation obligations.

**Relevant current manuscript locations**

- Introduction: motivates bi-Kappa physically but does not compare sampling difficulty with isotropic Kappa.
- Sec. II.F, Eq. (20): anisotropic bi-Kappa definition.
- Sec. III.A/D, Eqs. (23)–(25) and (29): anisotropic scaling to/from a spherical core.
- Sec. V, Eq. (35): arbitrary magnetic-field frame.

**Required manuscript changes**

1. State that, for this parameterization, bi-Kappa does not require a new radial probability law: an affine scaling reduces it to the same isotropic dimensionless core.
2. Identify genuine added implementation concerns: separate parallel/perpendicular scales, consistent convention/units, mapping from a field-aligned frame to an arbitrary Cartesian frame, anisotropic cutoff geometry, and component/joint validation.
3. Explain that these are important simulation-integration issues but should not be inflated into a new stochastic construction.
4. Relate the isotropic limit (`θ⊥ = θ∥`) explicitly to the standard Kappa definition adopted in the new background section.

**New analysis / experiment required**

**Textual clarification only.** Existing transform tests can be cited after they are documented; broader validation is consolidated in Experiment 1.

**Recommended response strategy**

Agree and clarify. A strong answer should say plainly that the radial sampling problem is not uniquely harder; the added value lies in anisotropic parameter handling and simulation coordinates.

**Dependency / overlap**

Overlaps R2.B/C3 and the magnetic-field/software positioning.

**Priority:** HIGH<br>
**Status:** UNADDRESSED

### R2.B — Add standard-Kappa/bi-Kappa background, methods, and parameter meanings

**Reviewer comment**

> Background and FrameworkPlease add a dedicated section briefly introducing the mathematical formulations of both the standard kappa and bi-kappa functions, along with current sampling methodologies. This section should explicitly include the definitions and physical meanings of all parameters.

The missing space after “Framework” is retained from the source email.

**What the referee is actually asking**

This is a **missing background, notation, physical-meaning, and literature issue**. The current probability preliminaries do not substitute for a coherent domain-specific definition of isotropic Kappa, bi-Kappa, their parameter conventions, physical interpretation, and established sampling routes.

**Relevant current manuscript locations**

- Sec. II.A/Table I: symbols are listed, but physical meanings and conventions are incomplete.
- Sec. II.F, Eqs. (20)–(21): bi-Kappa and bi-Maxwellian only.
- Sec. III: present sampling method.
- No dedicated standard isotropic-Kappa formula or current-methods comparison exists.

**Required manuscript changes**

1. Add a dedicated “Background and Related Work” section before the new method/implementation section.
2. Give the adopted normalized isotropic Kappa and bi-Kappa formulas, their relationship, normalization domain `κ > 1/2`, and finite-moment thresholds.
3. Define mathematical and physical meanings, domains, and units for `κ`, `θ⊥`, `θ∥`, velocity components, magnetic-field direction, density/normalization, and any cutoff `λ`.
4. Explain the thermal-speed/temperature convention; Kappa literature uses multiple conventions, so the manuscript must not call `θ` a temperature if it is a speed parameter without qualification.
5. Summarize current sampling methodologies and place the Beta-prime/Gamma-ratio route within that literature.
6. State the isotropic limit and distinguish standard Kappa, bi-Kappa, and capped bi-Kappa targets.

**New analysis / experiment required**

**Textual clarification only**, supported by the new literature analysis for R2.A1.

**Recommended response strategy**

Agree. Reorganize the scattered current Sec. II material rather than merely append another formula.

**Dependency / overlap**

Overlaps R1.7, R2.A1–A3, and R2.C3.

**Priority:** HIGH<br>
**Status:** PARTIALLY ADDRESSED IN CURRENT MANUSCRIPT

> C. Methodology and Mathematical Clarifications

### R2.C1 — Explain the projected one-dimensional sample's purpose

**Reviewer comment**

> (Page 2, Section C): Please clearly explain the utility of the projected one-dimensional sample and explicitly describe how this 1D sample relates to the final three-dimensional (3D) sample.

**What the referee is actually asking**

This is a **method-versus-validation ambiguity**. A reader may think the one-dimensional projected sample is an intermediate used to construct the 3D draw. In the current manuscript it is a marginal diagnostic derived from generated 3D vectors, not a generation step.

**Relevant current manuscript locations**

- Current Sec. II.C, especially Eq. (9): projected 1D marginal/sample discussion.
- Eq. (10): KDE diagnostic.
- Sec. II.H and Sec. IV.C/Figs. 5–6: one-dimensional Q–Q diagnostics.
- Sec. III.A–D: actual 3D generator, which does not use a projected sample.

**Required manuscript changes**

1. State immediately that a projected/marginal 1D sample is computed *from* the final 3D vectors for analysis and validation; it is not input to the sampler.
2. Define the projection or selected component mathematically and show how the empirical 1D sample relates to the analytic marginal obtained by integrating the 3D density over the other components.
3. Explain what it can test (component marginals) and what it cannot prove alone (radial law, directional uniformity, or full dependence structure).
4. Move generic validation material out of the core sampling derivation if that ordering continues to create ambiguity.

**New analysis / experiment required**

**Textual clarification only.** Direct joint/radial evidence is a separate requirement under R1.3.

**Recommended response strategy**

Clarify the misunderstanding and accept responsibility for the ordering/wording that caused it. Do not claim that three correct marginals alone prove the full joint law.

**Dependency / overlap**

Overlaps R1.3 and the planned separation of Mathematical Framework, Sampling Method, and Validation.

**Priority:** MEDIUM<br>
**Status:** PARTIALLY ADDRESSED IN CURRENT MANUSCRIPT

### R2.C2 — Trace the role of the transformation/Jacobian equation

**Reviewer comment**

> What is the specific role of Equation (11) in the overall bi-kappa sampling algorithm?

**What the referee is actually asking**

This is a **derivation and algorithm-flow ambiguity**. Eq. (11) states a general change-of-variables rule, but the manuscript does not clearly trace where its Jacobians enter the isotropization and radial-law derivations or distinguish derivation from executable steps.

**Relevant current manuscript locations**

- Current Sec. II.D, Eqs. (11)–(12): change-of-variables formulas.
- Sec. III.A, Eqs. (23)–(25): `V → U` scaling, with no explicit determinant derivation.
- Sec. III.B, Eqs. (26)–(28): Cartesian-to-radius and `R → T` transformations; Eq. (11) is cited.

**Required manuscript changes**

1. Say explicitly that Eq. (11) is a derivational probability-density rule, not a random draw in the implementation.
2. Apply it step by step to `V = A U`, with `|det A| = κ^{3/2} θ⊥² θ∥`, showing cancellation of anisotropic scale factors in `pU`.
3. Apply the spherical-coordinate Jacobian `r² sin θ` and angular integration `4πr²` to derive `pR`.
4. Apply the scalar Jacobian `dr/dt = 1/(2√t)` for `T = R²`.
5. Add a compact algorithm/derivation flow so readers can see how these transformations justify the sampling steps.

**New analysis / experiment required**

**Textual clarification only.** The derivation must be checked algebraically, but no numerical experiment is needed.

**Recommended response strategy**

Agree and make the two distinct uses of the Jacobian rule explicit.

**Dependency / overlap**

Overlaps R2.C3 and R1.3's radial-law validation.

**Priority:** HIGH<br>
**Status:** PARTIALLY ADDRESSED IN CURRENT MANUSCRIPT

### R2.C3 — Derive isotropization and correct ellipsoid/sphere language

**Reviewer comment**

> The transition from Equation (23) to Equation (25) is central to the proposed method. Please provide more detailed mathematical derivations to show how Equation (25) is obtained from Equation (23). Additionally, explain why the u-space forms an ellipsoid when (\theta_\perp \neq \theta_\parallel), and clarify how Equation (25) achieves spherical symmetry.

**What the referee is actually asking**

This is a **central derivation and geometric-clarity issue**. The current manuscript states the scaled density without showing the Jacobian. The wording also reveals a likely notation-induced misunderstanding: for the manuscript's Eq. (24) definition, constant-density surfaces are ellipsoids in physical `V`-space and spheres in normalized `U`-space—not ellipsoids in `U`-space.

**Relevant current manuscript locations**

- Sec. III.A, Eqs. (23)–(25).
- Sec. II.D, Eq. (11).
- Sec. III.D, Eq. (29), inverse scaling.
- Sec. V, Eq. (35), subsequent orthogonal frame rotation.

**Required manuscript changes**

1. Write the scaling in matrix form `V = A U`, `A = diag(√κ θ⊥, √κ θ⊥, √κ θ∥)`.
2. Substitute `V = A U` into Eq. (23), showing that the quadratic form becomes `ux² + uy² + uz²`.
3. Use `pU(u) = pV(Au)|det A|` and show the normalization-factor cancellation explicitly.
4. Explain geometrically that `vx²/(κθ⊥²) + vy²/(κθ⊥²) + vz²/(κθ∥²) = constant` is an ellipsoid in `V`-space when `θ⊥ ≠ θ∥`; dividing by its semi-axis scales maps it to a sphere in `U`-space.
5. Address the referee respectfully by clarifying the coordinate labels rather than reproducing a false statement that the defined `U`-space is ellipsoidal.
6. State that the later magnetic-frame transform is orthogonal and therefore preserves the normalized radial law.

**New analysis / experiment required**

**New analysis required.** Supply and independently verify the full Jacobian/geometric derivation; no new numerical experiment is required for the derivation itself.

**Recommended response strategy**

Agree that the derivation is too compressed, clarify the likely `V`/`U` notation misunderstanding, and provide the missing algebra without framing the referee as incorrect.

**Dependency / overlap**

Overlaps R2.A3, R2.B, R2.C2, and R1.3.

**Priority:** HIGH<br>
**Status:** PARTIALLY ADDRESSED IN CURRENT MANUSCRIPT

# 4. Zenitani / novelty-impact matrix

The central rule is: mathematical correctness can remain even when novelty does not. “Exact” and “rejection-free” may describe the untruncated high-level distributional construction if carefully scoped, but they cannot be used to imply discovery or comparative superiority.

| Current wording/location | Still defensible? | Why / impact of Zenitani et al. (2026) | Recommended replacement positioning |
|---|---|---|---|
| Title: “Sampling the Bi-Kappa Distribution” | Yes | The title does not itself claim priority. | Retain unless the revised journal framing needs an implementation qualifier; do not add “novel” or “new.” |
| Abstract: “no closed-form inverse cumulative distribution function is available” and direct sampling is “non-trivial” | Partly | The inverse-CDF statement may be true, but it does not establish absence of exact samplers. Standard-distribution/Student-t/Beta-prime constructions bypass inversion. | State the practical loading problem and available construction neutrally; cite existing methods. |
| Abstract: “We present an exact sampling algorithm” | Only with qualification | Distributional exactness of an untruncated construction can be true, but the Gamma-ratio construction is not newly discovered here. Current C++ defaults to truncation. | “We implement and validate a simulation-oriented bi-Kappa sampler based on the Beta-prime/Gamma-ratio construction, independently developed and equivalent to Zenitani et al. (2026).” Scope exactness to uncapped mode. |
| Abstract: “core mapping requires no rejection step” | Yes, narrowly | It is a property of the high-level uncapped construction, not novel. Gamma RNG internals may use rejection, and the current C++ wrapper always rejects against a cap. | “The uncapped high-level construction uses two Gamma variates and a uniform direction without an outer acceptance-rejection loop.” |
| Abstract: optional truncation has “negligible rejection rate” | **No — and the inference behind it is invalid** | Originally: no `λ`, rejected fraction, or bias reported. **Experiment 2 now shows the argument itself is a non-sequitur**: at κ=1.5, λ=50, TV = 6.3×10⁻⁴ yet the p99.9 speed is 24% low. A negligible discarded fraction does not bound tail fidelity, which is the whole point of a Kappa distribution. On a two-part criterion, **no λ in the ladder is negligible for κ ≤ 3/2**. | **Delete the negligibility argument, do not soften it.** Report `λ`, rejected fraction (= TV distance, closed form available), and a tail-quantile ratio, per κ. State that no practical λ suffices for κ ≤ 1. |
| Abstract/Summary: “runs in constant time per sample” / “produces one sample per call in constant time” | **No — now refuted, not merely unsupported** | **Experiment 3 measures 88–130 ns/sample, non-monotonic in κ** (fastest κ=1.5, slowest κ=2), consistent with the library Gamma generator switching algorithm as `shape = κ−1/2` crosses 1. Cap retries are variable-cost too. | **Delete the constant-time claim.** Report measured throughput with its hardware/toolchain, and say only that the uncapped algorithm uses a fixed number of *high-level* variate calls — which is a code-structure statement, not a cost statement. |
| Introduction: reliable Kappa generation presents a computational challenge because inversion is intractable and naive rejection is prohibitively inefficient | **No — the rejection half is now refuted by our own measurement** | It omits established non-inversion methods, and Experiment 3 measures a *dedicated* rejection sampler at 0.73–0.81 acceptance running **1.8×–2.6× faster than our own non-rejection method**. "Prohibitively inefficient" is not merely unsourced, it is wrong for this problem. An 2022 also undercuts "inversion is intractable". | Give a balanced method survey. State that a well-chosen envelope makes rejection *competitive and in our measurements faster*, and that generic rejection is inefficient only when the envelope is poor. |
| Introduction: “To resolve these computational bottlenecks, this paper introduces a fast and accurate algorithm” | No, as a novelty/superiority claim | The core is equivalent to published prior work, and no benchmark supports “fast” or bottleneck resolution. | Reframe around implementation, validation, frame support, reproducibility, and measured tradeoffs. |
| Implied claim that rejection sampling is the only available alternative | No | Zenitani's generalized Gamma-ratio method and related standard-distribution constructions are non-rejection alternatives at the high-level target-law stage. | Enumerate method classes and state explicitly that equivalent rejection-free constructions predate the revised paper's publication. |
| Sec. II.E / III.B: Beta-prime variable as a ratio of Gamma draws, culminating in `T ~ BetaPrime(3/2, κ−1/2)` | Correct but not novel | This is the central mathematical overlap for `r=0`, `q=κ+1`. | Retain as self-contained derivation needed for implementation; cite Zenitani and relevant statistical literature at the point of use. |
| Sec. III.A/D: uniform sphere plus parallel/perpendicular scaling | Correct but not novel as a construction | It is part of the same equivalent method. | Present as the implementation pipeline and connect it to anisotropic/frame APIs, not as a newly discovered mapping. |
| “Exact” terminology throughout | Conditionally | “Exact” may mean exact target-law sampling in ideal arithmetic, not exact floating-point values, finite moments, or novelty. It is false for a capped sample if the target named is Eq. (20). | Use “distributionally exact uncapped target in exact arithmetic” sparingly; use “capped conditional distribution” separately. |
| “Rejection-free” terminology throughout | Conditionally | It can describe absence of an outer rejection loop in the uncapped Gamma-ratio construction. It cannot describe the current default C++ package or capped mode, and says nothing about Gamma implementation internals. | Qualify every occurrence with “uncapped high-level core”; remove it from blanket package descriptions. |
| Sec. VI: C++ samplers implement the “exact, rejection-free algorithm” | No, for the current release | The C++ API mandates a cap and redraw loop. This is an implementation/documentation mismatch independent of priority. | Add a true uncapped API and document both targets, or call the released code capped. Emphasize tested API, integration, frame transform, licensing, and reproducibility. |
| Speed/performance claims | **No — measured and refuted** | Experiment 3: the released implementation is the **slowest** of three methods. Zenitani (2025) rejection is 1.8×–2.6× faster where κ > 1; A&M 2015's variant of our own construction is 1.4×–1.6× faster. | **Remove every superiority claim.** Report measured absolute throughput (≈8–11 M samples/s, stated with hardware) plus the two narrow favourable results: anisotropy and arbitrary-**B** rotation are free, and the cap's cost is confined to low κ. State the one genuine tradeoff — our route is defined for κ ≤ 1, where Zenitani's recommended envelope index is not. |
| Sec. VII: “We presented an exact, rejection-free algorithm” | Mathematically narrow but not novel | As a paper-level contribution claim it implies discovery and ignores the current capped implementation. | “We provide an open-source implementation and validation of a bi-Kappa loading construction equivalent to the generalized method of Zenitani et al. (2026), with explicit uncapped/capped modes and field-frame support.” Only use this if implementation catches up. |
| Sec. VII: “Numerical experiments … confirm the algorithm's correctness” | No, too broad | Existing tests cover only selected marginals/moments and not the radial law, low κ, truncation, or numerical failures. Zenitani makes independent rigorous validation a more plausible contribution, so the evidence must be stronger. | Enumerate exactly what was tested and with what diagnostics; reserve broad correctness language for evidence that includes Experiments 1, 2, and 4. |
| Summary/conclusion claim of a self-contained toolkit suitable for PIC initialization | Partly | The C++ library and frame support are real, but Python packaging, uncapped mode, provenance, and distributional tests are incomplete. | Describe concrete released components and tested use cases; do not imply untested production robustness. |

# 5. Proposed new validation / benchmark work

The work below is consolidated into **four** computational packages. Caption edits, literature review, derivations, API documentation, and regeneration from auditable scripts are revision tasks, but are not counted as separate scientific experiments.

## Experiment 1 — Direct untruncated radial, directional, and anisotropic validation

**Reviewer comments answered:** R1.3 directly; supports R1.5, R2.A2, R2.A3, R2.C1–C3.

**Scientific question**

Does the uncapped implementation produce the central target law `T = R² ~ BetaPrime(3/2, κ−1/2)`, a uniform independent direction, and the intended anisotropic/frame-transformed 3D distribution across low-to-high κ?

**Parameter sweep**

- `κ = {0.51, 0.55, 0.75, 1.0, 1.25, 1.5, 2, 5, 10}`. The first six cover the requested low-κ regime; do not use variance diagnostics at or below 1.5.
- At minimum, isotropic scales and `θ⊥:θ∥ = 1:2`; add a more extreme but physically defensible anisotropy only if it represents a supported use case.
- Field-aligned `B = ẑ` plus at least two non-axis-aligned, normalized magnetic-field directions.
- Multiple independent fixed seeds at a predeclared `N` (e.g. at least five replicates of `10^5`); add a larger stress sample only if memory/time is documented.

**Quantities to measure**

- ECDF/Q–Q comparison for the radius against the exact law, worked in `log R`. **Never form
  `T = R²`** — it overflows for `R > 1.3×10¹⁵⁴` even where `R` is representable, silently
  discarding exactly the heavy-tail draws the diagnostic exists to test.
- Use the bounded complement `W = 1/(1+T) ~ Beta(κ−1/2, 3/2)` for KS/Cramér–von Mises
  distances and tail quantiles. ⚠ **This supersedes the earlier recommendation of
  `Y = T/(1+T) ~ Beta(3/2, κ−1/2)`, which is withdrawn.** The two are exact bijections in
  exact arithmetic but not in floating point: for small `κ−1/2` the mass of `Y` piles up
  against 1.0 where no relative resolution remains. Experiment 1 measured 16.4% of `Y` values
  rounding to exactly 1.0 at κ = 0.55, producing a spurious KS failure at `√n·D = 51.8`, while
  `W` on the same data gave 0.751 — a clean pass. Experiment 4 computes `W` directly as
  `expit(−2 log R)`, so `T` is never materialized at all.
- Uniformity of `cos Θ` and `Φ`, and empirical radial–direction independence.
- Recovery of normalized radial invariants after anisotropic scaling and after arbitrary orthogonal field rotation.
- Cartesian marginal diagnostics as secondary checks.
- Replicate variability for moments only where those moments exist.

**Expected output**

- One compact multi-panel radial/directional validation figure.
- One table of κ, `N`, seeds/replicates, distributional distances, tail-quantile errors, and numerical failure counts.
- Optional appendix/supplement table for anisotropy/frame invariants.

**Necessity**

**Necessary.** It answers the First Referee's explicit central-law and low-κ request and supplies the most defensible validation contribution after the novelty change.

**STATUS: COMPLETE.** `experiments/exp1_radial_directional/`. 135 runs, 1.35×10⁷ draws,
released C++ sampler (not a Python reimplementation), uncapped throughout, 5 seeds × 10⁵ per
configuration, κ from 0.51 to 10.

| Check | Result |
|---|---|
| Radial law vs. exact | PASS |
| Directional uniformity | PASS |
| Radial–direction independence | PASS |
| Arbitrary **B**-frame invariance | PASS |
| Anisotropy ratio θ∥/θ⊥ = 2 | recovered to ≈0.4% |
| Non-finite draws | only at κ = 0.51, ≈5.8×10⁻⁴ |

Two methodological findings are results in their own right and are recorded above: the `W`
versus `Y` orientation, and the prohibition on forming `T = R²` or using `np.linalg.norm`.

## Experiment 2 — Untruncated versus component-wise cutoff characterization

**Reviewer comments answered:** R1.1, R1.2, and R1.3; supports R1.5 and R1.6.

**Scientific question**

How does the current component-wise cap change the target distribution, rejection rate, tails, and finite moments relative to the exact uncapped law?

**Parameter sweep**

- Modes: uncapped/off and capped.
- Representative `κ = {0.55, 1.0, 1.5, 2, 5, 10}`.
- `λ = {5, 10, 20, 50, 100}`, subject to pruning after analytic/short pilot estimates; keep the current default 20 and the example's 100 visible.
- `θ⊥:θ∥ = 1:2`; because the cap is normalized component-wise, confirm whether acceptance is invariant to these scales in exact arithmetic.
- Fixed `N` accepted samples plus total attempted draws, across multiple seeds.

**Quantities to measure**

- Attempted/accepted draws, rejection fraction, retry-count distribution, and max-retry failures.
- Empirical and, where feasible, numerically integrated mass outside the component-wise box.
- ECDF/quantile changes for radial and component distributions.
- Bias in second moments only for `κ > 3/2`; for `κ ≤ 3/2`, report robust quantiles/tail probabilities rather than “bias” relative to a divergent untruncated moment.
- Clarify that the cap guarantees finite moments by defining a different bounded distribution; it does not regularize Eq. (20) while leaving it unchanged.

**Expected output**

- Rejection-fraction-versus-`λ` plot stratified by κ.
- Table containing `κ`, `λ`, attempts, accepted draws, rejected fraction, moment/quantile effects, and mode.
- One explicitly reported capped example in the main manuscript; fuller sweep may be supplemental.

**Necessity**

**Necessary if truncation remains in the manuscript or released API.** Removing the capped mode from the paper would reduce the required scope, but it would not resolve the current code/paper mismatch by itself.

**STATUS: COMPLETE.** `experiments/exp2_cap_characterization/`. 280 runs, 2.8×10⁷ draws,
5 seeds × 10⁵, double, κ ∈ {0.75, 1, 1.5, 2, 5, 10} × λ ∈ {3, 5, 10, 20, 50, 100} plus uncapped
reference. λ = 100 is in the ladder because this section asks to keep the default 20 and the
example's 100 visible. Regeneration is bit-for-bit deterministic, re-checkable via `make verify`
against committed checksums.

**Acceptance was measured without instrumenting the released header.** The shipped predicate was
transcribed and evaluated on the *uncapped* run at the same seed; because `operator()` consumes
`x1, x2, cosΘ, φ` in the same order whether or not a cap is in force, the uncapped run **is** the
capped run's attempt stream. Verified, not assumed: the capped output is **bitwise identical** to
the uncapped draws restricted to the box in **240/240** (case, λ, seed) pairs
(`results/exp2_results.json:subsequence_check`).

**Five results, each of which changes what the manuscript may say.**

0. **⚠ "Negligible rejection ⇒ negligible bias" is a NON-SEQUITUR, and the manuscript makes it.**
   At κ = 1.5, λ = 50 the total-variation distance is 6.3×10⁻⁴ — the two laws are
   indistinguishable by any probability-based measure — yet the **p99.9 speed is still 24% too
   small**. TV bounds probabilities, not quantiles. The Abstract's "negligible rejection rate"
   argument must be **deleted, not softened**: a small discarded fraction does not license any
   statement about tail fidelity, which is precisely what a Kappa distribution is for. Any
   negligibility claim needs a two-part criterion (TV **and** a tail quantile), which is what
   Experiment 2 pre-registered.

1. **The rejected fraction IS the distortion.** Conditioning on an event of probability `p` gives
   density ratio `1_box/p`, so the total-variation distance from the target is exactly `1 − p`.
   Cost and bias are the same number; they must be quoted together.
2. **The θ's cancel exactly.** The predicate reduces to `√κ·maxᵢ|uᵢ| ≤ λ`, so `P(accept)` depends
   on `(κ, λ)` only — not on `θ⊥`, `θ∥`, or the field direction (the cap is tested *before* the
   frame rotation). Not a statistical coincidence: **0 disagreements in 6×10⁶ attempts** across
   60 comparisons on the shared RNG stream
   (`results/exp2_results.json:theta_independence_drawwise`). A closed form for `P(accept)` is
   derived and matches empirically.
3. **The cost decays only algebraically, as `λ^−(2κ−1)`.** Measured log-log slopes match the
   prediction to ~0.5%: −0.4996 vs −0.5 at κ=0.75, −18.80 vs −19 at κ=10. **The cap is not
   exponentially cheap in λ, and it is weakest exactly where Kappa distributions are physically
   interesting.** At κ = 0.75, λ = 50 still rejects 13.7%.
4. **⚠ The cap breaks axisymmetry about B.** The box is a *cube* in the isotropic
   u-coordinates, so a direction pointing at a cube corner has `√3` times more radial room than
   one along an axis. The conditioned law acquires a four-fold azimuthal modulation that the
   physical bi-Kappa does not have and **that no choice of θ can absorb**: `a₄ = 2⟨cos 4φ⟩`
   reaches −0.052 at κ=0.75, λ=3, i.e. **11.6 sampling standard deviations**. For a PIC
   initialization this is the most consequential finding in the experiment — a non-gyrotropic
   initial condition is a physics defect, not just a statistical one.

**Moments:** comparisons for κ ≤ 3/2 were **refused and labelled as refused**, not silently
omitted — the untruncated second moment `θ²κ/(2κ−3)` does not exist there, so a capped-vs-uncapped
variance ratio would compare a number against a divergent integral. Reported for κ = 2, 5, 10
only, where at λ=3 the cap removes over half the variance (ratio 0.459 vs 1.036 uncapped).

**Where the cap actually becomes negligible**, on the pre-registered two-part criterion
(`TV < 10⁻³` **and** p99.9 within 1%): κ=10 at λ=5, κ=5 at λ=10, κ=2 at λ=50, and **no λ in the
ladder qualifies for κ ≤ 3/2**. The heavy-tail cases are not merely expensive, they are
hopeless: at κ=0.75 the decay exponent is `2κ−1 = 0.5`, so λ=100 still rejects 9.7% and reaching
`TV < 10⁻³` would need λ ~ 10⁹. **Say that, rather than implying a large-enough λ exists.**

**Manuscript treatment (recommended):** validation and physics with cap **OFF**; document the
finite cap as an optional pragmatic finite-velocity-box **conditional target**; never present it
as a regularized or physically motivated Kappa model; quote cost and TV distortion as one number;
never quote a variance from a capped run at κ ≤ 3/2. **Consider making `no_cap()` the documented
default in the manuscript's examples** — the current default λ=20 is not negligible at any
κ ≤ 3/2, and neither is the committed example's λ=100.

**A trap for the unwary, worth a sentence in the Truncation section:** the capped sample has a
perfectly finite variance at κ ≤ 3/2 that any naive diagnostic will happily print — 1.28 → 65.2
at κ=0.75 as λ goes 3 → 100, growing as `λ^(3−2κ)` (measured slope 1.41 vs asymptotic 1.5), and
as `log λ` at κ=3/2. **It is a variance of the box, not of the plasma.**

## Experiment 3 — Reproducible absolute and comparative performance benchmark

**Reviewer comments answered:** R1.4 and the performance/outperformance part of R2.A2.

**Scientific question**

What are the measured cost, acceptance behavior, κ dependence, and failure modes of the implementations, and how do they compare with a relevant Zenitani piecewise-rejection implementation without pretending that two equivalent Gamma-ratio formulas are different algorithms?

**Parameter sweep**

- `κ = {0.55, 0.75, 1, 1.5, 2, 5, 10, 50}`.
- Uncapped and representative caps (at least default 20; include 100 if retained as an example).
- C++ scalar/header implementation and Python/NumPy vectorized implementation reported separately.
- Methods: present Gamma-ratio implementation; an independently verified implementation of the relevant Zenitani piecewise-rejection method; any equivalent Zenitani Gamma-ratio path labeled as an implementation cross-check, not a distinct algorithm.
- Repeated warm and cold batches over predeclared sample sizes; pin compiler flags, hardware, OS, RNG, and library versions.

**Quantities to measure**

- Samples/s and time/sample with median and spread over repeats.
- Startup/allocation versus steady-state time where relevant.
- Acceptance/rejection and retry counts for every rejection-based layer.
- Inf/NaN/overflow/zero-denominator/max-attempt failures.
- Output-distribution checks to ensure timing does not compare incorrect samplers.

**Expected output**

- One throughput/latency figure versus κ.
- One benchmark table with environment, methods, modes, acceptance, and failures.
- Committed benchmark script/config and machine-readable output.

**Necessity**

**Necessary if any speed, bottleneck-resolution, or outperformance claim is retained.** If all comparative claims are removed, the comparative part is useful rather than mandatory, but absolute performance characterization remains strongly recommended for the intended software-centered paper.

**STATUS: COMPLETE.** `experiments/exp3_benchmark/`. Correctness gate before timing, enforced
by the harness. Three methods, each transcribed from its primary source: the released
Gamma-ratio header; Abdul & Mace (2015) Eq. (22) — labelled an **implementation variant**, not
a rival algorithm, because it provably reduces to the same construction; and Zenitani (2025)
Pareto-envelope rejection, a genuinely distinct algorithm. 10 batches × 10⁶ draws per
configuration, one RNG type for all, κ ∈ {0.75, 1, 1.5, 2, 5, 10, 50}.

**Headline, unfavourable and reported as such: the released implementation is the slowest of
the three.** Zenitani (2025) rejection is 1.8×–2.6× faster where it applies (κ > 1); A&M 2015's
normal-triple variant is 1.4×–1.6× faster. Per-sample cost is **not constant in κ** (88–130 ns,
non-monotonic), so the Abstract's constant-time claim is refuted by measurement, not merely
unsupported. Anisotropy and arbitrary-**B** rotation cost ≈nothing (within 1–2 ns). The cap's
overhead is confined to low κ (+36% at κ=0.75), tracking Exp 2's rejection rate.

Measured Pareto acceptance reproduces Zenitani's published 0.806/0.785/0.750 at κ=1.5/2/5 to
three digits — independent evidence the transcription is faithful, which is the precondition
this section sets before benchmarking a published method.

Scope limits recorded on purpose: `n = κ/2` requires κ > 1, so the rejection method is
inapplicable at κ ≤ 1 (recorded, not skipped); and this is one machine and one toolchain —
enough to refute constant-time and fix a ≈2× ordering, not a cross-platform characterization.

## Experiment 4 — Low-κ floating-point and RNG robustness

**Reviewer comments answered:** R1.3; supports R2.A2 and the low-κ software positioning disclosed to the editor.

**Scientific question**

As `κ → 1/2+`, when does the Gamma denominator with shape `κ−1/2` become numerically zero/subnormal or produce overflow in the ratio, and can a stable implementation preserve the target distribution over a documented range?

**Parameter sweep**

- Moderate low κ from Experiment 1 plus `κ = 1/2 + {10^-1, 10^-2, 10^-3, 10^-4, 10^-6}` where representable and meaningful.
- C++ `std::gamma_distribution` and NumPy Gamma generation, with exact compiler/library versions.
- `float` and `double`; include a stable log-ratio or bounded-Beta route as a candidate mitigation, not as a presumed winner.
- Uncapped first; then identify how a cap changes failures/retries rather than hiding them.

**Quantities to measure**

- Frequency of zero/subnormal Gamma denominators.
- Inf, NaN, overflow, underflow, and max-attempt counts.
- Agreement of the bounded `W = 1/(1+T) ~ Beta(κ−1/2, 3/2)` with its Beta law, plus log-R
  quantile error. ⚠ **Not `Y = T/(1+T)`** — see the Experiment 1 note above; `Y` is unusable as
  a low-κ diagnostic and that earlier recommendation is withdrawn.
- Runtime and acceptance side effects of any mitigation.

**Required conceptual distinctions**

1. **Mathematical normalization boundary:** the distribution is normalizable only for `κ > 1/2`.
2. **Floating-point/RNG issue:** theoretically positive Gamma draws can underflow or be returned as zero in finite precision, making `X1/X2` overflow or divide by zero.
3. **Moment divergence:** the second moment diverges for `κ ≤ 3/2`; this is not the same as a divide-by-zero bug and must not be “fixed” by numerical tricks.

**Expected output**

- One robustness table by κ, precision, library, failure mode, and mitigation.
- A documented supported κ/precision range and explicit behavior outside it.
- A small diagnostic plot for `Y` or log-radius near the boundary.

**Necessity**

**Necessary** if the revised paper claims low-κ robustness, as proposed to the editor. At minimum, moderate low-κ distributional validation is required by R1.3 even if extreme-boundary robustness is not claimed.

**STATUS: COMPLETE.** `experiments/exp4_precision/`. Released C++ header under test; NumPy is
not used as a reference anywhere. Built twice from one source — Apple clang 21 / libc++ and
GCC 15.2 / libstdc++, arm64 — with real `float` and real `double` pipelines (nothing drawn in
double and cast). κ ladder 0.5001 … 1.5, seeds 4001–4005, 10⁶ draws each, 1.8×10⁸ total.

**Supported range of the released implementation as it stands:**

| Precision | Zero observed failures | Degraded | Unusable |
|---|---|---|---|
| `double` | κ ≥ 0.55 | κ = 0.51 (5.7×10⁻⁴) | κ ≤ 0.505 |
| `float` | κ ≥ 0.75 | 0.55 ≤ κ ≤ 0.60 | κ ≤ 0.51 |

Independently consistent with Experiment 1, which saw non-finite draws only at κ = 0.51 in
double at 5.8×10⁻⁴.

**libc++ and libstdc++ agree to within seed noise at every configuration.** There is no
standard-library dependence to report — a negative result, and worth stating as one.

**Q2 — `sqrt(x1)/sqrt(x2)` landed.** Evaluated on *identical variates* so the difference is
the fix and nothing else: float κ=0.55 loss 1.22×10⁻² → 5.7×10⁻³; float κ=0.60 1.5×10⁻⁴ →
2.9×10⁻⁵; double κ=0.51 8.2×10⁻⁴ → 5.7×10⁻⁴. Verified to agree to 2 ulp at κ = 2 where neither
formation can overflow. Regression tests `test_radius_formation` A1–A4 in `cpp/test_suite.H`.

**Q3 — log-domain construction measured, NOT landed.** Removes the remaining *spurious* loss
entirely, leaving only honest overflow (float κ=0.55: 0.99430 → 0.99986 finite; double
κ=0.505: 0.97589 → 0.99918). Validated 12/12 on KS *and* Cramér–von Mises against
`W ~ Beta(κ−1/2, 3/2)` at κ = 0.55, 0.75, 2, 5, α = 0.01 fixed in advance. It is not landed
because it would replace `std::gamma_distribution` with a hand-rolled generator, changing the
RNG stream and breaking seed-for-seed reproducibility against every prior version. **That is
an authors' decision, not a bug fix.**

**Q5 — the cap MASKS the failure.** Non-finite draws returned to the caller: **0 in every
configuration**, including where 93% of internal attempts are non-finite. A non-finite value
can never satisfy the box predicate, so the rejection loop silently redraws it. At κ = 0.51 in
float, 36% of internal attempts are non-finite and the user never learns of it. This must be
reported as failure *hidden* by truncation, never as failure solved.

**C5 verdict: PROMOTE, narrowly — as characterization plus a landed fix, not as novelty.** The
three conditions are met (a real failure in a plausibly-used released configuration; a
mitigation that provably preserves the target law; a measurable, reproducible improvement).
But ZUM2026 §4 already says "caution is needed" on `1/2 < κ ≤ 3/2`, so the honest framing is
that we *quantified an acknowledged caution* and *removed one avoidable cause of it*. No
discovery claim. See the contribution ledger in §9.

# 6. Manuscript-level revision map

| Current/proposed section | Action | Planned scope |
|---|---|---|
| Title | Retain provisionally | Current title is neutral. Change only if a software/implementation qualifier improves accuracy; do not add novelty language. |
| Abstract | **Rewrite** | Remove discovery implication, unsupported constant-time/negligible-rejection/correctness claims, and blanket exact/rejection-free language. State prior-art equivalence, concrete software contribution, explicit modes, and only completed validation/performance results. |
| I Introduction | **Substantially rewrite** | Replace bottleneck rhetoric with accessible motivation, balanced difficulty statement, existing solutions, remaining implementation gaps, and a bounded contribution list. Add requested foundational and computational citations. |
| New Background and Related Work | **Add dedicated section** | Standard isotropic Kappa and bi-Kappa formulas; conventions and parameter physics; normalization/moment conditions; Student-t/Beta-prime relation; prior sampling methods; Zenitani equivalence. This absorbs relevant current Sec. II material and answers R2.A1/B. |
| Mathematical Framework | **Retain but reorganize/condense** | Keep only probability tools needed later. Define PDF/normalization. Move projected marginals and Q–Q/KDE details to Validation. Expand change-of-variable/Jacobian derivation. |
| Sampling Method — uncapped core | **Rewrite and make self-contained** | Present `V→U`, radial Beta-prime/Gamma ratio, direction, and inverse scaling. Cite prior art. Use a flow diagram/pseudocode only if it improves traceability. Scope exact/rejection-free terms precisely. |
| Truncation / bounded mode | **Split into its own subsection** | Define component-wise cap geometry and conditional target, API off/on behavior, finite-moment implications, retry handling, and evidence from Experiment 2. Do not call it Eq. (20). |
| Validation | **Substantially expand and split by question** | Add radial/directional/joint checks, low-κ diagnostics, anisotropy/frame invariants, marginal/Maxwellian-limit plots, and moment tests only where defined. Report seeds, `N`, modes, κ, θ, cap, and statistical uncertainty. |
| Existing Figs. 2–4 / Tables II–III / Figs. 5–6 | **Regenerate and relabel** | Establish provenance, fix units/captions/normalization, remove or disclose percentile clipping, discuss κ=2 errors, and avoid treating marginals as proof of the full joint law. Consolidate figures if the expanded validation would otherwise become repetitive. |
| Performance | **Add subsection if claims retained** | Report Experiment 3 methods, environment, timing, acceptance, failures, and honest comparison. If no benchmark is added, remove performance superiority claims everywhere. |
| Magnetic-field transform | **Retain and expand validation/documentation** | Keep Eq. (35) concept, state orthogonality/invariance, document edge-case basis construction, and connect to existing transform tests and Experiment 1. Do not present field rotation as mathematical-sampler novelty. |
| Software / Code Availability | **Rewrite and expand** | Enumerate actual C++ and Python artifacts, license/DOI/version, RNG and seeding, supported modes/κ range, cap defaults, build/test instructions, and exact figure/benchmark reproduction commands. Reconcile claims with API reality. |
| Discussion / Limitations | **Add** | Independent-development chronology; Zenitani equivalence; non-novel core; parameterization conventions; finite-precision limits; moment divergence; cap tradeoffs; benchmark limitations; scope of validation. |
| Summary / Conclusions | **Rewrite** | Retire discovery, blanket correctness, constant-time, and unsupported speed claims. Summarize only demonstrated software, validation, robustness, frame, and reproducibility contributions. |

This map implies substantial reorganization. Equation/figure numbers will likely change; the final rebuttal should cite revised section labels and line numbers only after the revised manuscript is stable.

# 7. Cross-reviewer consolidation

Each reviewer comment still requires its own explicit response. The table identifies shared revision work so that evidence is produced once and cited in multiple responses.

| Consolidated revision package | Comments served | Shared deliverable | Traceability guardrail |
|---|---|---|---|
| Literature, prior art, and novelty correction | R1.7; R2.A1; R2.A2; R2.B | New Background/Related Work, verified references, explicit Zenitani equivalence, revised contribution list | R1.7 still gets its requested foundational/PDF citations; R2.A1 gets the full sampling survey; R2.A2 gets the difficulty/novelty answer. |
| Uncapped versus capped target-law separation | R1.1; R1.2; R1.3 | Two-mode definitions, API/documentation reconciliation, terminology audit | R1.1 gets the conceptual split; R1.2 gets exact run settings; R1.3 gets a reported capped example. |
| Direct radial/joint/low-κ validation | R1.3; supports R1.5; R2.A2; R2.C1–C3 | Experiment 1 plus the distributional part of Experiment 4 | Do not substitute marginal plots for the direct radial test or use divergent moments at low κ. |
| Truncation characterization | R1.1–R1.3; supports R1.5/R1.6 | Experiment 2 and capped-mode table | Report box geometry, `λ`, rejection, and distributional/moment effects separately. |
| Performance evidence or claim retirement | R1.4; R2.A2 | Experiment 3 and global performance-language audit | A response that only says the draw count is fixed does not answer either comment. |
| Isotropic versus bi-Kappa framework and parameter meanings | R2.A3; R2.B; supports R1.7 | Standard/bi-Kappa definitions, conventions, physical meanings, isotropic limit | State that anisotropy does not create a new radial law; list real integration challenges. |
| Transformation/Jacobian derivation | R2.C2; R2.C3 | Explicit `V→U`, spherical, and `R→T` Jacobians plus geometry explanation | Answer C2's algorithmic role and C3's algebra/ellipsoid language separately. |
| One-dimensional projections versus 3D validation | R2.C1; R1.3 | Reordered validation section and bounded claims | Explicitly say projections are outputs/diagnostics, not generator inputs and not proof of the joint law. |
| Figure/table provenance and statistical reporting | R1.2; R1.5; R1.6 | Reproducible run manifest, regenerated figures/tables, units/captions, uncertainty | Preserve R1.5's explicit κ=2 discussion and every R1.6 caption/axis request. |
| Software-centered contribution | R1.1; R1.4; R2.A2/A3/B | Honest API inventory, tested frame support, reproducibility, supported range, benchmark | Do not convert Yuri's preferred strategy into a referee demand or claim unsupported uniqueness. |

# 8. Major-revision acceptance checklist

## Source and traceability

- [ ] Use `paper/overleaf/main.tex` and `refs.bib` as the revision source of truth.
- [ ] Preserve an immutable copy/reference of the July 29 decision, First Referee PDF, July 30 disclosure, August 4 editor reply, and relevant Yuri correspondence.
- [ ] Keep editor instructions, reviewer requirements, coauthor strategy, and author inference labeled separately.
- [ ] Verify final revised section/equation/figure/table locations after renumbering; do not reuse obsolete page/line references.

## Point-by-point coverage

- [ ] R1.1 answered: uncapped and capped laws separated; terminology/API reconciled.
- [ ] R1.2 answered: mode, `λ`, rejection, `N`, seed, and provenance reported.
- [ ] R1.3 answered: radial law, low κ, bounded correctness claim, and capped example addressed.
- [ ] R1.4 answered: performance data supplied or claims removed.
- [ ] R1.5 answered: larger κ=2 errors discussed.
- [ ] R1.6 answered: all requested caption, curve, normalization, label, and unit changes made.
- [ ] R1.7 answered: foundational, kinetic-code, random-variate, and PDF sources/definition added.
- [ ] R2.A1 answered: current sampling literature reviewed comprehensively and accurately.
- [ ] R2.A2 answered: accessible difficulty explanation, solved/unsolved split, corrected novelty, and evidence-based comparison supplied.
- [ ] R2.A3 answered: isotropic versus bi-Kappa challenges stated without manufacturing novelty.
- [ ] R2.B answered: dedicated formulas/methods/parameter-meaning background added.
- [ ] R2.C1 answered: projected 1D sample identified as a validation diagnostic and related explicitly to 3D samples.
- [ ] R2.C2 answered: Eq. (11)'s derivational role and Jacobian uses traced.
- [ ] R2.C3 answered: Eq. (23)→(25) derivation and ellipsoid-to-sphere geometry supplied.

## Novelty and claim discipline

- [ ] Zenitani et al. (2026) explicitly cited and discussed in the manuscript and point-by-point response.
- [ ] `r = 0`, `q = κ + 1` equivalence stated accurately and checked against the published notation.
- [ ] Mathematical novelty claim for the Gamma-ratio construction removed.
- [ ] Independent-development chronology stated factually, without converting it into a priority claim.
- [ ] Student-t/Beta-prime and other relevant prior literature verified and cited.
- [ ] No statement implies rejection sampling is the only pre-existing alternative.
- [ ] “Exact” and “rejection-free” terminology is consistent and restricted to the uncapped high-level construction.
- [ ] “Constant time,” “fast,” “prohibitively low,” “outperforms,” and “confirms correctness” appear only if defined and supported; otherwise removed.
- [ ] Manuscript and rebuttal claims are mutually consistent.

## Software and computational evidence

- [ ] C++ API provides/documented uncapped behavior if the manuscript claims it.
- [ ] Capped behavior, default `λ`, support geometry, retry limit, and failure behavior documented.
- [ ] Python deliverable described exactly as released; no packaged-module claim if it remains notebook-only.
- [ ] Experiment 1 radial/directional/anisotropic validation completed and committed reproducibly.
- [ ] Experiment 2 cutoff characterization completed if truncation remains.
- [ ] Experiment 3 benchmark completed if performance/comparative claims remain.
- [ ] Experiment 4 low-κ robustness completed if low-κ numerical robustness is claimed.
- [ ] Mathematical singularity, floating-point/RNG failure, and moment divergence are distinguished everywhere.
- [ ] Experiments are reproducible from committed code, fixed configs/seeds, environment metadata, and machine-readable outputs.
- [ ] Figures and tables are regenerated from committed code; no undocumented manual notebook state is required.
- [ ] Benchmark claims are supported by data and include unfavorable/neutral results where observed.
- [ ] Distributional validation is performed before timing any implementation used in a comparison.
- [ ] Existing frame-transform, seeding, and API tests still pass after later code changes.

## Manuscript quality and submission package

- [ ] Standard and bi-Kappa definitions use one explicit parameter convention and correct units/physical meanings.
- [ ] PDF definition and normalization are stated or cited.
- [ ] All figures have self-contained captions, normalized-histogram statements, curve identities, parameters, modes, units, and sample sizes.
- [ ] All tables identify mode/configuration and report appropriate uncertainty.
- [ ] All new references and DOIs are verified against primary sources.
- [ ] Code/archive version and Zenodo DOI cited in the manuscript match the tested release.
- [ ] Clean manuscript PDF produced from the canonical source.
- [ ] Marked-up manuscript PDF produced.
- [ ] Complete point-by-point response produced, with one explicit response per R1/R2 item.
- [ ] List of changes produced.
- [ ] July 30 novelty context provided to both returning referees as the editor suggested.
- [ ] Final response avoids claiming that Yuri's strategy choices were referee requirements.

## Stop gate

- [ ] The matrix has been reviewed by the authors before any manuscript revision, rebuttal drafting, code change, benchmark, or new scientific run begins.

# 9. Alignment

This section exists so that this matrix, `introduction_rewrite_proposal.md`,
`../literature/step1_claim_audit.md`, `paper/reference/README.md`, the experiment READMEs and
`refs.bib` can be re-synchronized against evidence as it arrives, without anyone having to
re-derive which document is currently stale. **Update this section every time an experiment
finishes or a primary source is obtained.**

## 9.1 Alignment protocol

When new evidence lands, walk these in order. The order matters: a claim cannot be written
before the evidence that licenses it exists, and a document cannot be trusted while a document
it depends on is stale.

1. **Record the raw result** in the experiment's own `README.md` + `results/`. Nothing is
   summarized before it is recorded.
2. **Update this matrix** — the relevant `R*` item's *Evidence now available* block, its
   **Status**, and the Experiment status in §5.
3. **Update the contribution ledger** (§9.3). Move items between tiers only on evidence.
4. **Update `../literature/step1_claim_audit.md`** if the evidence touches a literature claim;
   re-issue that claim's PASS / MODIFY / DROP / BLOCKED verdict.
5. **Update `introduction_rewrite_proposal.md`** last, because its prose depends on all of the
   above.
6. **Touch `refs.bib` only** when bibliographic data are primary-source verified.
7. **Log it** in §9.2 with the date and what changed.

**Superseded-recommendation rule.** When a recommendation in any of these documents is
withdrawn, do not delete it silently. Mark it withdrawn *in place*, state what replaced it and
what evidence caused the change. A returning referee who read the earlier reasoning will
otherwise see an unexplained reversal. Current live example: `Y = T/(1+T)` → `W = 1/(1+T)`.

**Claim-discipline rule.** Every quantitative statement carries its scope. `0.73–0.8`
acceptance is Zenitani (2025) *over κ ≥ 3/2 only*. `κ > 3/2` is the domain of the published
Gamma-ratio loaders under their thermal-speed convention, not ours. Numbers travel with their
parameter range or they do not travel.

## 9.2 Alignment log

| Date | Trigger | Documents updated | Net effect |
|---|---|---|---|
| 2026-08-17 | Step-1 literature audit (10 primary sources read) | `step1_claim_audit.md`, `paper/reference/README.md` | Prior-art concession moved back four years, 2026 → **2022** (Zenitani & Nakano Table I, Algorithm 1-1). Three citation errors caught. |
| 2026-08-17 | **Experiment 1 complete** | exp1 `README.md`, this matrix §5 + R1.3 | R1.3 discharged except the capped example. `Y` recommendation **withdrawn** in favour of `W`; `T = R²` and `np.linalg.norm` prohibited. |
| 2026-08-17 | **Abdul & Mace 2015 + Abdul, Matthews & Mace 2021 obtained** | `step1_claim_audit.md`, `paper/reference/README.md`, §9.3 | P4-b unblocked and **inverted** — see §9.4. A&M 2015 is a genuine trivariate loader and 2021 is 2D3V at 6.7×10⁷ particles/species. Kills any "not demonstrated in multidimensional PIC" framing. |
| 2026-08-17 | **Solar-wind anisotropy sources obtained** | `paper/reference/README.md`, `step1_claim_audit.md` | P2-a unblocked via Hellinger 2006 + Matteini 2007 + Pierrard 2016. Marsch 1982 full text still **BLOCKED** (abstract-level only). Ion anisotropy⊗Kappa observation **BLOCKED** — no primary source found. |
| 2026-08-17 | **Experiment 4 complete** | exp4 `README.md`, this matrix §5 + §9.3, `cpp/test_suite.H` | Supported range established. `sqrt(x1)/sqrt(x2)` fix quantified and regression-tested. Cap shown to **mask** failure. C5 → PROMOTE narrowly, as characterization + fix. |
| 2026-08-17 | **Experiment 2 complete** | exp2 `README.md`, this matrix §5 + R1.3 + §4 + §9.3 | R1.3 fully discharged. Cap is a conditional target whose rejected fraction **is** its TV distortion, decaying only as `λ^−(2κ−1)`; it **breaks axisymmetry about B** (a₄ at 11.6σ); and **"negligible rejection ⇒ negligible bias" is shown to be a non-sequitur** (TV 6.3×10⁻⁴ with a 24%-low p99.9). Manuscript validation must use cap OFF, and the Abstract's negligibility argument must be deleted rather than softened. |
| 2026-08-17 | **Provenance freeze across all three experiments** | exp1/exp4 `GNUmakefile` + `*_analyze.py` + `README.md`, exp2 `README.md`, this matrix §5 | exp1 had committed `checksums.sha256` but **no committed command to generate or verify it**; exp1 and exp4 recorded no git revision or sampler hash. All three now emit checksums, expose `make verify`, and pin the sampler by content. Re-running the (deterministic) analyses changed **only** the `environment` block. Two wrong numbers corrected in §5 against `exp2_results.json`: 200/200 → **240/240** pairs, 5×10⁶ → **6×10⁶** attempts. |
| 2026-08-17 | **Cap default flipped to uncapped** | `cpp/bi_kappa_distribution.H`, `cpp/test_suite.H`, `cpp/main.cpp`, `README.md`, `usage.dox`, this matrix R1.1 + §9.6 | The default target law is now the bi-Kappa distribution; a finite cap is opt-in. Verified behaviour-neutral for every reported number: all three experiments set the cap explicitly and **regenerate bit-for-bit identically** against the modified header. Tests B5–B8 added. |
| 2026-08-17 | **Experiment 3 complete** | exp3 `README.md` + `results/`, this matrix §5 + R1.4 + §9.3, §9.5 blocker 2 cleared | R1.4 answerable. **Unfavourable headline: our implementation is the slowest of three tested.** Zenitani (2025) Pareto rejection 1.8×–2.6× faster where κ > 1; A&M 2015 normal-triple variant 1.4×–1.6× faster; per-sample cost **not constant in κ** (88–130 ns). Every speed/superiority claim retired by measurement. Anisotropy + **B**-rotation shown free. Transcription faithfulness evidenced by reproducing Zenitani's published acceptance to 3 digits. |
| 2026-08-17 | **§9.4 differentiator corrected** | this matrix §9.4 + §9.3 | The claim that the anisotropic bi-Kappa loader is "never written down" in prior literature is **withdrawn** — ZUM2026 Algorithms 3.1/3.2 write it out with distinct θ∥, θ⊥, as this repository's own audit (P4-g, §3) already recorded. Tier-2 candidate retired. Surviving differentiators restated around arbitrary-**B**-frame loading, implementation, validation and target-law semantics. |

## 9.3 Contribution ledger, by evidence

**Tier 1 — established, evidence in hand.**

| Contribution | Evidence |
|---|---|
| Simulation-oriented open-source C++ implementation of bi-Kappa loading | `cpp/`, released header, test suite |
| Arbitrary magnetic-field-frame loading | Exp 1 Block B — frame invariance PASS on 3 **B** directions |
| Explicit uncapped/capped target-law semantics in the API | `no_cap()` sentinel + `capped()`; tests `test_no_cap_semantics` B1–B4 |
| Quantitative characterization of what the cap does to the target law | Exp 2 — closed-form `P(accept)`, `λ^−(2κ−1)` cost law, TV ≡ rejected fraction, and the broken-axisymmetry result |
| Direct 3-D validation of the central radial law, incl. `1/2 < κ ≤ 3/2` | Exp 1 — 1.35×10⁷ draws, radial/directional/independence/anisotropy all PASS |
| Reproducibility: committed scripts, fixed seeds, manifests, environment capture | Exp 1 + Exp 4 provenance; `.gitignore` whitelist |

**Tier 2 — candidate, evidence-dependent.**

| Contribution | Status |
|---|---|
| Low-κ finite-precision characterization + the `sqrt(x1)/sqrt(x2)` fix | **EARNED** by Exp 4, narrowly. Frame as quantifying an acknowledged caution (ZUM2026 §4), never as discovery. |
| Log-domain small-shape Gamma mitigation | **MEASURED AND VALIDATED, NOT LANDED.** Authors' call — it changes the RNG stream. |
| Performance characterization | **EARNED as characterization, and it is a NEGATIVE result for us.** Exp 3: our implementation is the slowest of three; Zenitani (2025) rejection is 1.8×–2.6× faster; cost is not constant in κ. Publishable as honest measurement + claim retirement. The only favourable parts are narrow: arbitrary-**B** rotation and anisotropy are free, and the cap's cost is confined to low κ. **No superiority claim of any kind survives.** |
| ~~The anisotropic bi-Kappa scale-mixture written down explicitly~~ | **RETIRED 2026-08-17.** ZUM2026 Algorithms 3.1/3.2 write the anisotropic loader down explicitly, with distinct θ∥, θ⊥. See the withdrawal notice in §9.4. Do not resurrect. |

**FROZEN 2026-08-17.** All four experiments are complete and the literature audit is closed
except the §9.5 blockers, so the ledger below is the final evidence-backed contribution
statement. The section-by-section prose plan that follows from it is
`manuscript_revision_plan.md`. Move an item between tiers only on new evidence, and log it
in §9.2.

**Rejected / retired — do not resurrect.**

| Claim | Why |
|---|---|
| Discovery of the Gamma-ratio / Beta-prime construction | Zenitani & Nakano 2022 Table I Algorithm 1-1; A&M 2015 Eq. (22); ZUM2026 Eq. (6) |
| The anisotropic bi-Kappa loader is undocumented in the literature | **Retired 2026-08-17.** ZUM2026 Algorithms 3.1/3.2 write it out with distinct θ∥, θ⊥; see §9.4. |
| Any speed, throughput, or efficiency advantage | **Retired 2026-08-17 by measurement.** Exp 3: ours is the slowest of three methods tested. |
| "Multivariate" Student-*t* equivalence attributed to A&M **2014** | 2014 is strictly 1-D (Bailey polar transform, one deviate at a time) |
| Existing Kappa loaders not demonstrated in multidimensional PIC | A&M+Matthews 2021 is 2D3V, 6.7×10⁷ bi-Kappa particles/species |
| Zenitani et al. mischaracterize Abdul & Mace 2015 | All seven characterizations checked verbatim; every one accurate |
| `Y = T/(1+T)` as a low-κ diagnostic | Exp 1 — 16.4% round to exactly 1.0 at κ = 0.55; spurious `√n·D = 51.8` |
| "Fast" / "resolves bottleneck" / "prohibitively low acceptance" / "constant time" | No benchmark exists. Zenitani (2025) publishes the *opposite* expectation. |
| Cross-platform bitwise reproducibility | Never demonstrated. Exp 4 shows agreement in *failure statistics*, which is not the same thing. |
| The cap as a physical regularized-Kappa model | It is a pragmatic component-wise box → a conditional target. Regularized Kappa is a separately defined smooth distribution (Scherer et al.). |

## 9.4 The differentiator — **corrected 2026-08-17, and it is narrower than previously written**

> ⚠ **Withdrawn in place.** This section previously claimed that "the anisotropic bi-Kappa scale
> mixture is *used but never written down* anywhere in the prior literature", listing ZUM2026
> with only the quote *"We can easily extend it for θ∥ ≠ θ⊥."* **That claim is not supportable
> and is withdrawn.** It also contradicted this repository's own literature audit, which had the
> correct reading all along (`../literature/step1_claim_audit.md` P4-g and §3): ZUM2026's
> Algorithm 3.1 box carries the anisotropic scaling **explicitly**, in executable form:
>
> ```
> v∥   ← θ∥ R x (2U₃ − 1)
> v⊥1  ← 2 θ⊥ R x √(U₃(1−U₃)) cos(2πU₄)
> v⊥2  ← 2 θ⊥ R x √(U₃(1−U₃)) sin(2πU₄)
> return v∥, v⊥1, v⊥2
> ```
>
> and the surrounding text says *"Trivial modifications for θ∥ ≠ θ⊥ are also included in the
> procedure."* Re-verified against the primary PDF during this session. **Do not write any
> sentence claiming the anisotropic loader is undocumented in the literature** — a returning
> referee holding ZUM2026 would open the algorithm box and find it.

**What the prior art does contain.** The anisotropic loader is written down (ZUM2026 Algorithms
3.1 and 3.2, with distinct θ∥, θ⊥). The isotropic scale mixture is written down repeatedly:
Abdul & Mace (2015) §III B carries a general covariance **R** in Eq. (17) and then sets
**R** = σ²**I**; Abdul (2018) PhD §2.3 is the fullest written statement and is isotropic only;
Abdul, Matthews & Mace (2021) runs a genuine bi-Kappa in 2D3V citing the 2015 loader.

**What actually survives as ours**, stated so that nothing depends on the retracted claim:

| Surviving | Why the prior art does not cover it |
|---|---|
| **Arbitrary magnetic-field-frame loading** (C2) | ZUM2026 returns `v∥, v⊥1, v⊥2` — components **in the field-aligned frame**. It does not rotate into an arbitrary **B** direction in a simulation's global frame. Our Eq. (35) transform does, and Exp 1 Block B validates it on 3 **B** directions to ~1×10⁻¹⁵. |
| **Open-source, tested, simulation-oriented C++ implementation** (C1) | ZUM2026's release is a Zenodo Jupyter notebook (audit §1-7). No C++ library, no field-frame handling. |
| **Explicit uncapped/capped target-law semantics** (C3) | No prior source separates the two target laws, and none quantifies the cap's distortion. |
| **Direct 3-D validation incl. 1/2 < κ ≤ 3/2** (C4) | Prior loaders are presented with far lighter validation; none tests radial law + directional uniformity + independence + frame invariance at these κ. |
| **Finite-precision characterization** (C5) | ZUM2026 §4 flags the zero-division caution qualitatively; Exp 4 quantifies it and removes one avoidable cause. |

The honest framing throughout: **the construction is prior art; the engineering, the validation,
the frame handling, and the target-law discipline are the contribution.** A novel construction
is not claimed, and neither is a novel anisotropic extension.

## 9.5 Open blockers

| # | Blocker | Blocks | Status |
|---|---|---|---|
| ~~1~~ | ~~Experiment 2~~ | ~~R1.1, R1.2, R1.3 item 5~~ | **CLEARED 2026-08-17** |
| ~~2~~ | ~~Experiment 3 prerequisite: a faithful, independently distributionally validated implementation of at least one competing method~~ | ~~R1.4, R2.A2~~ | **CLEARED 2026-08-17.** Both A&M 2015 Eq. (22) and Zenitani (2025) §2 transcribed from primary source and implemented; all three methods pass the distributional gate before any timing was believed. Zenitani's published acceptance reproduced to 3 digits. |
| 3 | Marsch et al. (1982) full text | any page/figure-level claim from it | Not open access; abstract verified via Crossref publisher deposit. Cite for abstract-level propositions only. |
| 4 | A primary ion observation fitted with an anisotropic bi-Kappa | the anisotropy⊗κ claim *for ions* | **BLOCKED — none found.** Pierrard et al. (2016) supplies it for **electrons**. Do not silently generalize to ions. |
| 5 | Scherer, Fichtner & Lazar (2017) EPL 120, 50002 | regularized-Kappa citation at its origin | Optional; the 2019 anisotropic paper serves our comparison better. |
| 6 | Published versions of Bale 2009 and Pierrard 2016 | quoted detail from either | On-disk copies are arXiv preprints. Check before any quotation ships. |

## 9.6 The cap/default decision — settled 2026-08-17

Recorded here because R1.1 item 5 and R1.2 both depend on what the released API actually does
by default, and because it is a deliberate behaviour change that a returning referee may notice.

**Previous behaviour.** `max_normalized_velocity` defaulted to `20.0` in `param_type`, in the
value constructor, and in `define(...)`. A caller who omitted the argument silently sampled the
**conditional box law**, not bi-Kappa. `no_cap()` existed but had to be asked for. The committed
example used `100.0`; `usage.dox` and the README quick-start both showed `20.0`.

**Final behaviour.** The default is `no_cap()` on every construction path. A finite cap is
opt-in. `param_type::capped()` still reports which law is active. Nothing else changed: the
core mapping, the RNG consumption order, and the capped-mode predicate are untouched.

**Why, on evidence — not preference.**

| Evidence | Consequence for a defaulted caller |
|---|---|
| Exp 2: TV distance to bi-Kappa is *exactly* the rejected fraction, decaying only as `λ^−(2κ−1)` | At λ=20 the default was non-negligible **everywhere in κ ≤ 3/2** on the pre-registered criterion |
| Exp 2: small TV does not bound tail quantiles (κ=1.5, λ=50 → TV 6.3×10⁻⁴, p99.9 24% low) | The old default could not be defended as "close enough" even where TV looked tiny |
| Exp 2: the box is a cube in normalized coordinates → four-fold azimuthal modulation at 11.6σ | A defaulted PIC initialization was **not gyrotropic** — a physics defect, silently |
| Exp 4: a non-finite draw can never satisfy the box predicate, so the loop redraws it | Capped mode **hid** generation failure; at κ=0.51 float, 36% of internal attempts were non-finite and 0 were reported |

The manuscript claims the sampler targets Eq. (20). The default now does.

**Compatibility consequences.** Callers who passed a cap explicitly are unaffected — that is
every experiment in this repository, and both existing test suites. Callers who relied on the
default now get the untruncated law: a behaviour change, and the intended one, since the old
behaviour was the defect R1.1 identified. Verified neutral for the scientific record: exp1,
exp2 and exp4 were rebuilt against the modified header and **regenerate bit-for-bit identical
raw output**, checked against the committed `checksums.sha256`. Only the recorded
`sampler_header_sha256` changed in `results/`.

**One residual footgun, documented rather than papered over.** The seed is the sixth `define(...)`
argument, after the cap, so a caller who omits the cap *cannot* pass a seed positionally —
`define(κ, θ⊥, θ∥, ub, 12345)` sets a **cap** of 12345, not a seed. Adding an `int`-taking
overload would silently reinterpret existing calls that pass an integer cap, so it was not
added. README and `usage.dox` both call this out, and `dist.seed(s)` remains the unambiguous
route.

**Not changed: `bi_maxwellian_distribution`.** Its default is still `20.0` and it has no
`no_cap()`. The Exp-2 argument does not transfer — for a Gaussian, truncation at 20 thermal
speeds removes ~10⁻⁸⁸ of the mass, which is negligible under any criterion including the
two-part one. Making the two classes symmetric is an **authors' consistency choice**, not an
evidence-driven fix, and is deliberately left open.

**Test coverage.** `test_no_cap_semantics` B1–B4 (pre-existing) plus:

| Test | Guards |
|---|---|
| B5 | the default is uncapped on *all five* construction paths, including `param_type` itself |
| B6 | a defaulted sampler behaviourally produces out-of-box draws, not merely a flag saying so |
| B7 | omitting the cap is **bitwise** identical to naming `no_cap()` at the same seed |
| B8 | the two laws remain measurably different — guards against a future change making the cap a no-op |

`test_radius_formation` A1–A4 continue to cover the `sqrt(x1)/sqrt(x2)` formation, which was
already landed before this session (`bi_kappa_distribution.H:273`); A4 pins the two formations
to 2 ulp at κ=2. Full suite passes, `main.exe` exits 0.

**The capped target, for the manuscript.** Define it once, as

```
f_cap(v) = f_kappa(v) * 1{v in B_lambda} / P(v in B_lambda),
B_lambda = { v : |v_x|/theta_perp <= lambda, |v_y|/theta_perp <= lambda, |v_z|/theta_par <= lambda }
```

and never as Eq. (20), the untruncated law, a regularized Kappa, or a physically motivated tail
regularization. `P(v in B_lambda)` has the closed form in the Exp-2 README, depends on
`(κ, λ)` only, and equals `1 − TV(f_cap, f_kappa)`.
