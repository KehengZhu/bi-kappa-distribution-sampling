# Experiment 7 — pre-registered confirmatory protocol

**Protocol 3.0.0. Status: frozen before any datum on seeds 9001–9010 existed.** This document
and `config/protocol.json` are committed in a source-only commit. The analysis refuses to run
unless the SHA-256 of both files matches the values recorded in the run manifest, and every
decision rule below is evaluated by code that reads `config/protocol.json` rather than by a
constant written into the analysis.

**This is the third confirmatory holdout.** The first ran on seeds 7001–7010 under protocol
1.3.0 and returned **NO-GO**; the second ran on 8001–8010 under protocol 2.0.0 and also
returned **NO-GO**, on gate G1 alone. Both are preserved unmodified, in commits `45d3ef8` and
`e5c9837`, neither is reopened, and both seed blocks are spent. §8 permits exactly one path
after a failure — identify a concrete defect, fix it, freeze a new implementation hash *and a
new protocol document*, draw a further disjoint seed block, and rerun — and §2.7 records the
two defects the second holdout found, one in the implementation and one in this protocol.
Nothing in 3.0.0 relaxes a threshold, removes a cell, changes a family's α, or excludes a
result.

---

## 1. Why a second experiment rather than a rerun of Experiment 6

Experiment 6 measured the right things and its raw data stand. What it cannot do is
*confirm* anything, for three reasons that are properties of its design rather than of its
numbers:

1. **Its primitive-selection rule was not the frozen one.** The plan required that every
   pre-declared log-quantile lie inside its acceptance interval. The analysis instead
   required 90 % of them to (`exp6/analyze.py:369`). Under a perfectly correct sampler the
   frozen rule fails with probability 0.984 — the expected miss count over the 420 intervals
   is 4.10 — so it is not a conservative rule but one with a 98.4 % type-I error rate, and
   the 90 % replacement has no power against anything short of a tenfold degradation. Either
   way, the rule that selected the primitive is not the rule that was pre-declared.
2. **Its acceptance-rate check changed its predictand after seeing the result.** The
   transcribed method's measured acceptance rate deviated from the published Eq. (2) by −127
   standardized deviations; the analysis then compared it against Γ(α+1)/(1 + w) instead
   (+2.2 σ). That re-derivation is mathematically correct, and this protocol adopts it — but
   adopting it *after* the measurement means the check stopped being a test.
3. **Several gates cannot fail.** G3 evaluates `all([])` when no negative-control cell was
   produced, so it passes on absent evidence. G4 is a hard-coded constant on this host. G5
   passes when a ratio is finite. G6 checks that a checksum file exists rather than that it
   verifies. The E0 baseline verdict is computed and never read.

Experiment 6 is therefore retained, unmodified, as **exploratory** evidence: its mechanism
decomposition, its oracle adjudication and its conditioning result are descriptive
measurements that do not depend on the acceptance rules above. Experiment 7 is the
confirmatory test, run on a disjoint seed block, against the rules fixed here, with a
different implementation.

## 2. What is under test

The candidate is the released loader `cpp/bi_kappa_distribution.H` at version **2.2.0**, whose
radius is built in the log domain and carried, with the rest of the deterministic map from the
drawn variates to the returned velocity, in the accumulator of §2.7.1 — `double` for a `float`
instantiation, the working type otherwise — with one rounding where the component is
materialized.  2.0.0 carried the same construction and the representability defect of §2.6.1;
2.1.0 corrected that and carried the resolution defect of §2.7.1.  Both are superseded before
release and neither is the candidate. The comparator is the same file at version 1.0.0, vendored
unmodified as `src/legacy/bi_kappa_distribution_v1.H` and referred to as **LEGACY**.

Neither the Gamma-ratio construction, log-domain Gamma generation, the small-shape underflow,
nor the low-parameter denominator hazard is claimed as new here; all are prior art and are
attributed in `README.md`. What is under test is the integration of an established primitive
into a complete anisotropic bi-Kappa loader and the validation of that loader's
three-dimensional output, tail, cap conditioning, portability and cost.

### 2.1 The selected primitive, and why it is not the one Experiment 6 froze

Experiment 6 piloted two log-scale primitives and selected **LOG-PUB** (Liu, Martin & Syring
2017) over **LOG-ID** (Ahrens–Dieter shape boosting on the log scale) on a median-cost
tie-break. The plan's own selection rule makes that tie-break conditional — the faster method
wins *unless* it has "materially greater complexity or weaker documented bounds" — and on the
evidence the condition is met:

- **LOG-PUB's cost diverges at the top of its domain.** Its acceptance rate is
  `Γ(α+1)/(1 + w)` with `w = α/(e(1−α))`, which has a pole at `α = 1`, i.e. at `kappa = 3/2`.
  Measured: 1.24 attempts per draw at `α = 0.25`, 2.29 at `α = 0.75`, 4.49 at `α = 0.9`, and
  **37.6 at `α = 0.99`**, where it costs 31× the released path. Experiment 6's kappa ladder
  steps from 1.0 to 1.5 and so passes over this pole without sampling it.
- **It is undefined at and above `α = 1`**, where it falls back to `std::gamma_distribution`.
  At its own benign control, `kappa = 2`, *every one of its 200 000 pilot draws* took that
  fallback, so the control exercised the standard library rather than the candidate.
- **LOG-ID has neither property**: it is valid for every `α > 0`, is rejection-free at the
  boosting step, costs one Gamma and one uniform per call independently of `kappa`, and is
  roughly ten lines.

Experiment 7 therefore uses the boosting identity. Two further changes follow from defects
that the audit of Experiment 6 identified in how that identity was implemented, and both are
corrections rather than preferences:

- **The uniform is built from engine bits, not from `std::uniform_real_distribution`.** The
  plan requires "a documented mapping from engine bits to an open-interval `U ∈ (0,1)`
  wherever a logarithm is taken". Delegating to the standard library does not provide one:
  `uniform_real_distribution<float>` returns exactly `1.0f` under libc++ at a measured rate of
  3.0e-8 — whose logarithm is zero, silently turning a `Gamma(a)` draw into a `Gamma(a+1)`
  one — while libstdc++ clamps the same case to the largest value below 1, placing an atom
  there. The two builds were therefore running different discrete samplers. The replacement
  maps `b = digits − 1` engine bits to `(k + 1/2)/2^b`, which is exactly representable, is
  strictly inside `(0,1)` for every `k`, and needs neither a redraw nor a clamp.
- **The Gamma and normal primitives are implemented in the header.** `std::gamma_distribution`
  uses Ahrens–Dieter GS in libc++ and Marsaglia–Tsang in libstdc++, so the stream, the variate
  count and the failure envelope were properties of the library. With Marsaglia–Tsang and a
  polar normal in the header, the sampled sequence is a function of the engine alone. This is
  what lets G4 below be an equality rather than an interval overlap.

### 2.2 Declared changes in numerical behaviour

- The random stream changes for every seed. No seed-for-seed continuity with 1.x is claimed
  or possible; the stream change is the point, not a side effect.
- Within the 2.x line the stream is a function of the engine and the parameters alone, so
  libc++ and libstdc++ builds agree bit for bit at a fixed multiply-add contraction setting.
  Reaching that took two changes, not one: moving the Gamma, normal and uniform primitives
  into the header, and then replacing the `cos θ, φ` direction with a rejection method that
  uses no trigonometry. The second was necessary because clang fuses a `sin`/`cos` pair on
  one argument into `__sincos_stret`, whose sine differs from the standalone `sin` by an ulp
  on about one argument in a thousand, while GCC does not — so the first change alone left
  the two builds disagreeing, on the compiler's account rather than the library's.
- Across architectures nothing is bitwise: `log` and `exp` are not correctly rounded and the
  arm64 and x86_64 implementations differ, so cross-architecture agreement is statistical.
  G4 tests each of the two claims with the rule appropriate to it.
- `seed(int)` now calls `reset()`. Before 2.0.0 it did not, and under libstdc++ the cached
  second normal variate of `std::gamma_distribution` survived a reseed, so `seed(s)` followed
  by an odd number of draws did not reproduce. This is a defect fix; it changes behaviour.

### 2.3 Amendment 1.1.0 — three criteria added, none relaxed

Made **before any Experiment 7 draw existed**, in response to an independent audit of
Experiment 6's evidence that completed after §1–§9 were first written. Recorded here and in
`config/protocol.json` under `amendments`, with the reason, so that the addition is auditable
rather than silent. It adds acceptance criteria; it weakens none.

1. **A returned value can be finite and still be wrong.** Experiment 6's terminal classifier
   read `return method_finite ? kCatFinite : kCatHonestOverflow`, so a method that returned a
   finite vector for a mathematically unrepresentable draw was scored as a *success*. That is
   exactly the regime where the legacy form is worst: Experiment 6's own oracle measured up
   to 0.48 relative error on the surviving radius once the denominator went subnormal, and no
   gate looked at it. A draw is now classified **FINITE_BUT_WRONG**, and counted as a loss,
   when its returned radius differs from the arbitrary-precision oracle by more than
   `2^-(digits/2)` — losing more than half the significand of its type, i.e. 1.1e-8 in
   `double` and 2.4e-4 in `float`. Stating the threshold in bits makes it type-aware, and it
   leaves the candidate a wide margin: its own error is about `eps·|log R|`, which is 1.5e-13
   at the `double` overflow threshold and 5.3e-6 at the `float` one.
2. **Audit coverage is stated by margin, and the margin is kept.** Experiment 6's stratified
   audit sampled bulk failures at 1 in 4957 while its plan said "send every public-path
   failure"; with zero disagreements in 7039 of them, the one-sided 95 % bound still admits
   about 14 900 misclassified failures among 34.9 million.

   Amendment 1.1.0 first responded by requiring *every* public-path failure. Costing that
   showed it is infeasible and, more to the point, unnecessary: the honest failures alone
   come to 17.5 million across the frozen matrix — about 5.4 CPU-hours of oracle time — and
   almost all of them are draws whose intended value lies thousands of log units outside the
   type's range, where the classification is not in doubt because no arithmetic could have
   returned them. Writing down a rule that cannot be followed is exactly how Experiment 6
   arrived at 1 in 4957.

   **Amendment 1.2.0** therefore states the rule by margin. An attempt is *decision-relevant*
   — and adjudicated in full — when its intended log-component lies within **20 natural-log
   units** of a type limit, when the two methods reach different terminal categories, when the
   candidate returned a FINITE_BUT_WRONG value, when a subnormal or zero denominator
   coincided with a representable target (the avoidable-loss candidates), or when it falls in
   the uniform 1-in-10⁴ sample. Everything else is *unambiguous*: it is sampled at 1 in 10³,
   **and the margin is asserted in working precision for every one of them**, so the claim
   rests on an inequality rather than on a sampling rate.

   Twenty log units is four thousand times the worst disagreement ever measured between the
   working-precision reference and the 100-digit oracle (7.3e-12 in `double`, 4.2e-3 in
   `float`). Beyond it a misclassification is impossible, not merely improbable. Both the
   audited and the total count of every stratum are reported.
3. **Accuracy is never pooled across precisions.** Experiment 6's headline contrast — 0.48
   against 7.8e-6 — put a `float` worst case next to a `float` figure without labelling
   either, while the `double` worst cases are 0.46 and 1.1e-13. Six of the orders of magnitude
   in that contrast came from the type, not from the method. Every accuracy statistic carries
   its precision.

### 2.4 What the frozen battery can detect, computed before the run

`config/power_study.py` measures the power of the frozen families against survivor
conditioning — the effect the battery is used to certify *absent* — at the frozen sample size,
from simulated draws only. The finite survivors of a loader that loses a fraction `q` of its
draws are the target conditioned on `Z < -log q`, so the question is answerable exactly.

Power at `n = 10^6` (200 replicates per cell), against a loss fraction of:

| statistic | 1e-2 | 1e-3 | 1e-4 | 1e-5 |
|---|---:|---:|---:|---:|
| Anderson–Darling on `Z` | 1.00 | 0.99 | 0.01 | 0.01 |
| Kolmogorov–Smirnov | 1.00 | 0.14 | 0.00 | 0.01 |
| Cramér–von Mises | 1.00 | 0.21 | 0.01 | 0.01 |
| **F1 global** | 1.00 | 0.93 | 0.01 | 0.01 |
| **F4 exceedance** | 1.00 | 1.00 | 1.00 | 0.07 |

Minimum loss fraction detectable at power ≥ 0.90: **1e-3 for F1, 1e-4 for F4**. Neither
detects 1e-5, and the protocol says so here rather than discovering it afterwards.

Level, checked over 2000 replicates at `n = 5×10^4`: 0.0095 (AD), 0.0085 (KS), 0.011 (CvM),
0.0085 (F1), 0.0055 (F4), against a nominal 0.01; the Anderson–Darling p-value is uniform
under the null (KS against `U(0,1)`, p = 0.25). The statistics hold their size, so the power
figures above are comparable.

This is why F4 exists. Experiment 6 had no equivalent, and its radial battery's detection
floor was about a 0.3 % truncation at `n = 5×10^5` — while the candidate's residual loss in
the cells it certified was 0.076 % and 0.015 %, below that floor by factors of 4 and 18. "The
log path passes the radial test" therefore meant "its residual truncation is under 0.3 %",
not "it is correct". The P1 output must carry what F4 needs — exceedance counts above
`z0 = -log q0` and the excesses themselves — or the confirmatory run cannot detect the effect
it exists to rule out.

### 2.5 Amendment 1.3.0 — the loader battery gains a tail statistic

Also made **before any holdout datum existed**, and it repairs an inconsistency between
§2.4 and §5.4 that was mine.

§2.4 establishes that nothing bulk-weighted has power against survivor conditioning below a
1e-3 loss fraction. §5.4 then required negative control NC1 to be *detected* at 1e-3 **and**
1e-4 with power ≥ 0.90. But the frozen F5 loader battery — direction uniformity,
independence, frame invariance, cap law, anisotropy — contains no upper-tail statistic at
all. The protocol therefore demanded a power the battery it froze could not deliver, and a
perfectly correct candidate would have failed F6, then G3, then the verdict.

Measured at production scale, `n = 5×10^5`, 40 injections per effect:

| effect | F5 independence | upper-tail exceedance |
|---|---:|---:|
| loss fraction 1e-3 | 0 / 40 | 40 / 40 |
| loss fraction 1e-4 | 0 / 40 | 40 / 40 |

The arithmetic is not subtle. Removing 50 draws from 500 000 cannot move an 8×8 χ² whose
cells expect 7 800 apiece; the same 50 draws are essentially every exceedance above
`q0 = 1e-4`.

**F5 therefore gains an exceedance statistic on each uncapped loader cell's recovered radial
law**, mirroring F4: the count above `z0 = -log q0` for `q0 ∈ {1e-2, 1e-3, 1e-4}` against
`Binomial(n_attempted, q0)`, plus a KS of the excesses against `Exp(1)`. Capped cells are
excluded, because under a cap the accepted radius is truncated at a direction-dependent bound
and the count is not binomial.

Two conditions are part of the rule, not caveats on it.

> **Superseded in part by §2.6.2.** The two conditions below were 1.3.0's answer to the
> censoring problem, and the second of them — the count member's applicability rule — is
> **withdrawn**, because the quantity it could not observe turns out to be recorded. The
> members themselves, their statistics and their thresholds are unchanged. Read §2.6.2 for
> what governs.

**A draw whose `Z` transform underflowed counts toward every threshold** — it is further into
the tail than any of them — and the denominator stays `n_attempted`. Testing only the
resolved subset would condition on resolvability, which is monotone in the tail, i.e. exactly
the bias these statistics exist to detect. *This is not the same as a draw whose velocity
overflowed.* An overflowed draw has a perfectly ordinary `Z` and must be counted by it.

**The count member applies at a threshold `q0` only where `q0` exceeds the cell's
honest-overflow rate `f`.** The reason is what is observable. When `q0 > f`, the threshold
`z0 = -log q0` lies below the overflow threshold `z_f = -log f`, so every attempt above `z0`
is either a returned survivor above `z0` or an overflowed attempt, and both are counted
exactly. When `q0 < f`, every attempt above `z0` has overflowed, and separating those above
`z0` from those merely above `z_f` needs the intended value of a draw that has none — the
loader returned no number for it. Such a cell is reported **not applicable** at that
threshold: never passed on it, never failed on it, with its honest floor printed on the row.

Without that condition the member would silently demand that every uncapped cell lose less
than 1e-4, which case C4 exists precisely to violate — its floor is 8.25e-4, so honest
overflow alone would give `p ≈ 1.7e-223` and fail a correct candidate on the physics of the
floating-point type. The excess member is unaffected and applies at every threshold.

This **adds** a test to the family rather than relaxing one — it gives the candidate one more
way to be rejected — and it closes a genuine gap: a battery certifying a heavy-tailed law had
nothing in it that looks at the tail.

Four smaller corrections travel with it, recorded under `corrections` in
`config/protocol.json`:

- **F1 is Anderson–Darling, KS and Cramér–von Mises on `Z` only.** The Beta test on `W` is
  dropped. The probe emits no `W` sample; `Z` is a monotone transform of `log W`, so the two
  routes are near-redundant; and G0 already validates the `log q` evaluator against an
  arbitrary-precision incomplete beta to `1e-10`, which is a stronger check on the transform
  than comparing the two routes to each other would be.
- **NC3's effect size** is the loss fraction of the worst cell *among those the fidelity claim
  rests on* — the uncapped loader cells whose loss is not dominated by honest overflow — not
  the global maximum, which is float `kappa = 0.5001`, where almost every draw is
  unrepresentable and detection is trivial.
- **G5's interval is 99 %**, the cluster-bootstrap level in `statistics.bootstrap.conf`. §6
  said 95 % in prose; the machine-readable value governs, and it is the conservative choice
  for a gate on an upper limit.
- **F2 counts misses over every *resolved* interval**, which is what its frozen
  Poisson-binomial null was built over. The informativeness flag of §5.1 is published as a map
  and does not change the count — excluding non-informative cells would make F2 unevaluable
  against its own null.

## 2.6 Amendment 2.0.0 — what the first holdout found, and what changes

Made **after** the first holdout was read, which is what makes it a different kind of
amendment from the three above and why it carries a new major version. §8 allows it only on
the terms §8 states, and those terms are met: two concrete defects were identified, both are
corrected, the implementation hash is new, this document is new, and the seed block is new.

### 2.6.1 The implementation defect

The released loader decided representability by comparing a logarithm — `log R`, or
`log R + log|g_j|` — against `log(max())`. That comparison cannot be made to agree with the
arithmetic it predicts, because `log(max())` is itself a rounded value and can land on either
side of the true logarithm of the largest finite number. On the tested libm it lands *above*
in `float`: `exp(log(FLT_MAX))` is exactly `+inf`, while in `double` `exp(log(DBL_MAX))` is
finite. So the defect's visibility is a property of a rounding direction, not of the type.

One draw in 14 233 536 audited attempts hit it, at `float kappa = 0.51`. Its `log R` was
`0x1.62e43p+6`, bit for bit `log(FLT_MAX)`; it passed `log R <= logMax`; `exp(log R)`
overflowed; and the sampler returned `(-inf, +inf, +inf)` — and, the overflow never having
been detected, counted nothing. All three components it should have produced were
representable, the largest at 1.92e38 against a limit of 3.40e38, a margin of 0.571 natural-log
units, because every `|g_j|` is below one and pulls the product back under the limit. The draw
was scored `log_primitive_failure`, which is an avoidable loss, and it failed **G1** (avoidable
loss must be exactly zero) and **G2** (oracle disagreements must be zero). The two recorded
losses and the two oracle disagreements are the same draw seen in the libc++ and libstdc++
streams, not four events.

The correction, in `bikappa_detail::materializeComponents`, materializes each component and
decides from the component: the test is applied to the number that will be returned, so the
test and the result agree by construction, whichever way the library rounds. A second defect
of the same family is corrected with it: in capped mode the non-representability counter was
testing the *normalized* coordinate the cap predicate is written in rather than the velocity,
and the two differ by `theta` and by the rotation — up to 2.33× on the C5/C6 geometry. The
counter now answers the question `n_nonfinite()` documents.

The implementation under test at that point was **2.1.0**. It is superseded in turn by 2.2.0
under §2.7.1; 2.0.0 and 2.1.0 are both superseded before release, so no two samplers share a
version string.

### 2.6.2 The protocol defect

Amendment 1.3.0 gave F5 two upper-tail members and stated both against the **untruncated**
law. That is not the law the data obey. Near `kappa = 1/2` the target puts non-zero
probability outside every finite floating-point range, so the loader cannot return the far
tail — and is right not to. Honest overflow therefore right-censors the returned sample, and
at a **direction-dependent** point: a draw comes back iff `R max_j |g_j(n)| <= max()`, and
`max_j |g_j|` varies over the sphere by up to `sqrt(3) theta_max / theta_min`. A cutoff
inferred from a cell's total overflow rate is the average of that boundary, not the boundary.

Measured over 2000 replicates of a **perfectly correct** loader at the frozen production size:

| member | C4 (double κ=0.505) | C3 (float κ=0.55) | C1 (double κ=0.51) | C0 (double κ=2) |
|---|---:|---:|---:|---:|
| 1.3.0 excess, q₀=1e−2 | **1.000** | 0.115 | 0.009 | 0.009 |
| 1.3.0 excess, q₀=1e−3 | **1.000** | **1.000** | 0.009 | 0.012 |
| 1.3.0 count, q₀=1e−4 | **1.000** | 0.584 | 0.010 | 0.007 |

against a nominal 0.010. The three F5 rejections that failed **G3** in the first holdout were
those cells. The failure was a property of the null, not of the candidate.

Both members keep their statistic, their threshold and their place in the joint Holm. What
changes is the null each is compared against.

**The count member** is computed from the **intended** `Z` of every attempt instead of from
the radius recovered from the returned vectors. The probe records `log_r_ref` — the radius the
attempt carried — for every attempt including the overflowed ones, and an uncapped cell runs
its core mapping exactly once per attempt, so its record count equals its attempt count and
that sample is complete and uncensored. It was on disk for the first holdout too; the analysis
simply was not reading it. Against it the count above `z₀` is `Binomial(n_attempted, q₀)`
**exactly**. Family F4 was already computed this way on the scalar phase, which is why F4
passed the first holdout while F5 did not.

That also **withdraws the side condition** 1.3.0 needed. `count_member_requires_q0_above_honest_floor`
declared the member "not applicable" wherever it would have misfired; the quantity it could not
observe turns out to be observable, so the member now applies at every threshold of every
uncapped cell. C3 and C4 at `q₀ = 1e−4` are now tested where they previously were not. This
adds tests; it removes none.

**The excess member** stays on the returned draws — the fidelity claim is about the population
the loader hands back — and is given the null those draws obey. `Z` is independent of the
direction, so conditional on its own direction a returned draw above `z₀` is `Exp(1)` truncated
to `(z₀, C_i]`, with

    C_i = Z( log max() − log max_j |g_j(n_i)| )

computed from that draw's own recovered direction and the cell's declared `theta` and `ub`.
The member is a Kolmogorov–Smirnov test of the per-draw probability integral transform

    U_i = (1 − exp(−(Z_i − z₀))) / (1 − exp(−(C_i − z₀)))

against `Uniform(0,1)`. Where no attempt can overflow, `C_i` is effectively infinite, `U_i`
reduces to `1 − exp(−(Z_i − z₀))`, and because a Kolmogorov–Smirnov statistic is invariant
under a common monotone transform of the data and the null CDF, **the number returned is
identical to 1.3.0's** — measured on C0, the two agree to 1.1e−16. The replacement equals the
frozen test wherever the frozen test was valid and is defined where it was not. It is also the
construction this protocol already used for the capped cells, whose accepted radius is
truncated at a direction-dependent bound in exactly the same way.

Neither member needs quadrature, a grid, or numerical integration: both are per-draw arithmetic
in the frozen `Z` evaluator that G0 validates against an arbitrary-precision incomplete beta.

Level and power are measured in `docs/revision/experiments/f5_tail_calibration.md`, on
simulated cells and on the spent block, before any datum on 8001–8010 existed.

### 2.6.3 Negative control NC3, respecified for the same reason

1.3.0's NC3 removed every draw above a single direction-independent cutoff from an uncensored
sample and required the battery to detect it. But that is, up to the direction dependence,
exactly what honest overflow does to a **correct** loader. The battery "detected" it only
because its null was the untruncated law — the measured NC3 power was type-I error wearing a
power label, and it was the same defect that failed the candidate on C3 and C4.

NC3 is split into the two questions that were tangled together:

- **NC3a**, honest censoring alone at the cell's own direction-dependent boundary, correctly
  reported. The battery must **not** reject. This is a level, is published as one, and carries
  no power threshold.
- **NC3b**, honest censoring **plus** a further `q` fraction of the draws the type *would* have
  allowed back, removed silently. That is the defect. Required power ≥ 0.90 at the
  pre-registered `F6.nc3_excess_loss_fractions`.

NC1 and NC2 are unchanged in construction; NC1 now computes its censoring boundary from the
source cell's real `theta` and `ub` rather than from an isotropic stand-in.

### 2.6.4 Scope

Acceptance is limited to the **supported environment**: arm64 macOS under both standard
libraries. The cross-architecture claim is **withdrawn** — removed from §9 and published as a
limitation — rather than left as a permanently open gate on evidence this project cannot
obtain, since §5.3 refuses emulation as closure and no native x86_64 host is available. What
remains in G4 is the sharper of the two predictions and the one that is decidable here: within
2.x the stream is a function of the engine alone, so libc++ and libstdc++ must agree **bitwise**,
with no tolerance at all. Narrowing the claim does not license ignoring contrary evidence: a
cross-architecture comparison filed anyway that *disagrees* still fails G4.

## 2.7 Amendment 3.0.0 — what the second holdout found, and what changes

Made **after** the second holdout was read, on the terms §8 states and no others: two concrete
defects were identified, both are corrected, the implementation hash is new, this document is
new, and the seed block is new. The second holdout's own record is in commit `e5c9837` and is
not edited, reinterpreted or reopened.

What it established stands. Four of the five gates that had failed or stood open on 7001–7010
passed: G2 recorded **0 oracle disagreements in 14 235 020 adjudicated attempts**, against 2
before; F5 passed jointly over 52 tests, including the three cells whose 1.3.0 nulls had been
invalid; G4 closed with bitwise cross-standard-library agreement; G5 came in at 0.854× against
a bound of 2. G1 failed, for two reasons.

### 2.7.1 The implementation defect

One avoidable loss, at **float κ = 0.505, seed 8005, attempt 157855**, seen identically in the
libc++ and libstdc++ streams — one draw, not two events. Its largest intended component lies
`3.12e-08` natural-log units below `log(FLT_MAX)`. That is about half an ulp, so the correctly
rounded `float` is finite and the draw is one the type can hold; the loader returned a
non-finite component and counted it as a loss.

The cause is resolution, not a wrong comparison. 2.1.0 already decided representability from
the materialized component rather than from a logarithm, and G2 is clean because the
classification the loader reached agrees with the oracle's. What it could not do is *resolve*
the margin. Two working-precision quantities on the path to that component are coarser than
`3.12e-08`:

- **The `float` log radius.** `log R = (log X1 − log X2)/2` with
  `log X2 = log Y + (log U)/a`. The division by `a = κ − 1/2 = 0.005` multiplies the error of
  `log U` by 200, and the result carries about `3.5e-06` — a hundred times the margin.
- **The `float` order-unity vector.** `g_j = sqrt(κ) θ_j n_j`, formed in `float` from a
  `float` direction, carries about `3.9e-08` of relative error. That alone is larger than the
  margin: an ulp of `float` at the overflow threshold is `6.0e-08` relative, so rounding `g`
  before multiplying by `R` decides the question by a coin toss.

**Release 2.2.0 separates the variates from the map.** The variates — `X1`, `Y`, `U`, and the
two disc coordinates — are still drawn in `RealType`, from the same engine bits, through the
same acceptance tests, so the law being sampled does not change. The deterministic function
that turns those variates into a returned velocity — the logarithms, the sum, `sqrt(κ)`, the
direction lift, the rotation and the product `R g_j` — is evaluated in the accumulator
`bikappa_detail::log_accumulator<RealType>`, and the value is rounded to `RealType` exactly
once, where the component is materialized. For `float` the accumulator is `double`: the
resolution of the representability decision moves from about `5.3e-06` to about `1.6e-14`,
against an ulp of `6.0e-08`. For `double` the accumulator *is* `double`, so the arithmetic is
unchanged and 2.2.0 returns **the same bits as 2.1.0 for every `double` configuration in the
matrix**, capped and uncapped, with identical attempt and non-finite counts.

One `float` acceptance test moves with it, and it is a correction rather than a side effect:
the disc rejection `d1² + d2² < 1` needs 49 bits to be exact for two `float` coordinates, so
the `float` form could accept a point outside the disc or reject one inside it. It is now
evaluated in the accumulator, where it is exact.

**Two deterministic regressions are frozen with the fix**, in `make selftest`, replayed from
recorded bit patterns rather than from an engine, so that they survive any later change of
seeds, of phase or of the ladder:

- the losing draw itself — its five `float` variates as the run recorded them — which must now
  be returned, with its margin reported; and
- the boundary it sits on: 32 magnitudes stepping in eighths of an ulp from −2 to +2 ulp of
  `FLT_MAX`, plus the two points `1e-12` either side of the rounding midpoint, each decided as
  the hardware's own `float` conversion decides it. The midpoint itself is excluded and the
  exclusion is part of the rule: it is the one magnitude whose answer comes from
  round-half-to-even rather than from an inequality, no finite-precision reconstruction of it
  can be relied on to land on the right side, and an exact tie has probability zero under a
  continuous law. The same sweep in `double` is stated at `1e-10` relative, which is a
  thousand times what `double`'s own `eps |log R|` can resolve there and a hundred times finer
  than the accuracy the protocol requires of the returned radius.

**One consequence is a narrowing and is recorded as one.** For a `float` run the stabilized
candidate's radius is now the `double` evaluation of the same primitives, which is exactly what
the paired diagnostic layer's *reference* computes. The paired layer therefore cannot judge the
`float` candidate's accuracy — the two agree by construction — and no `float` FINITE_BUT_WRONG
verdict from that layer is evidence of anything. Accuracy at `float` is adjudicated by the
100-digit oracle alone, which shares no arithmetic with either, and amendment 1.2.0's margin
rule already sends every attempt within 20 natural-log units of a type limit to it in full. The
`double` paired layer is unaffected, because there the reference still does not materialize
`X2` and the candidate still does not share its formation.

### 2.7.2 The protocol defect

F2 compared its miss count against a Poisson–binomial over **independent** intervals. The
intervals were not independent, for two reasons of very different size.

**The one that mattered is a property of the experiment, not of the statistic.** Every
configuration of a replicate was driven by one `std::mt19937` seeded with the replicate's seed,
so all 26 of them — thirteen κ values in two precisions — consumed the same uniform stream, and
`X1` in particular is a function of the engine alone and is *identical* across the κ ladder at
a fixed precision. Their quantile errors then co-move almost exactly. On 8001–8010 this is
visible directly: seed 8004 gave 103 of 130 signed quantile errors positive
(`p = 1.2e-11`), seed 8005 gave 31 of 130 (`p = 1.8e-09`), and 20 of the 23 misses came from
those two replicates. The effective number of independent units was about five, not 650, and a
count statistic compared against a 650-trial null is over-dispersed by roughly that factor.
The first holdout's 0 misses of 650 — probability about `2e-03` under the same null — is the
same over-dispersion seen from the other side.

**The second is irreducible and was never modelled at all.** The five levels of one
configuration are order statistics of one sample, so their miss indicators are dependent
whatever the streams do.

3.0.0 answers the first by **removing the coupling** rather than by modelling it, and the
second by **computing it**.

- **P1 and P5 give every configuration of a replicate its own engine stream**, by the declared
  formula in §3. The 130 configurations are then independent by construction, which is what
  F2's global rule needs and what no amount of re-specification could have supplied from
  coupled data. P2, P3, P4 and P6 are untouched: their families decide by Holm, by Simes or by
  an exact per-attempt identity, every one of which is valid under arbitrary dependence between
  configurations, so there is nothing there for this to fix and no reason to disturb them.
  `make selftest` checks that the 130 derived streams are distinct and that none of them is a
  spent, performance or fixture seed.
- **The null becomes the 130-fold convolution of the per-configuration miss-count law.** That
  law is computed, not assumed: an interval at level `p` covers its quantile iff the number of
  sample points at or below it lies in `[lo+1, hi]`, and under the null the vector of those
  five counts is the running sum of a multinomial over the six intervals the five levels cut.
  Nothing about the sampler enters — the null's own probability integral transform is
  `Uniform(0,1)` by definition — so the law is a property of `(n, the five levels, the five
  brackets)` and is identical for every configuration. `config/make_protocol.py` evaluates it
  by Monte Carlo at 2×10⁷ replicates under a recorded seed, publishes the resulting
  distribution in `config/protocol.json`, and the analysis re-convolves it and refuses to run
  if the frozen critical value does not come back out.

**What does not change.** The statistic is the same: the number of order-statistic brackets
that fail to cover their quantile, over every resolved interval. The family α is the same,
0.005. The set of intervals counted is the same, 650, and the informativeness map still changes
no count. The per-level miss probabilities are the same analytically computed achieved
coverages.

**What the change costs the candidate.** The critical value moves from **13 to 14** of 650, and
the attained level from 0.00464 to 0.00250. The Monte Carlo's marginal miss rates reproduce the
analytic achieved coverages to `6e-05`, and the critical value is unchanged when the two
worst-resolved cells of the per-configuration law are zeroed, doubled, or doubled together with
the cell below them; all four checks are published in `config/protocol.json` under
`F2.far_tail_stability`. **On the second holdout's own count of 23 the new rule would still
have rejected, at `p = 4.3e-07`.** The replacement is therefore not a relaxation dressed as a
correction: it is a design change that supplies the independence the frozen statistic always
required, and it leaves that statistic able to fail.

### 2.7.3 Scope

Unchanged from §2.6.4: arm64 macOS under both standard libraries, cross-architecture withdrawn
and published as a limitation.

## 3. Seeds

**Production seeds: 9001–9005. Performance-block seeds: 9006–9010.** Both blocks are declared
here, in `config/protocol.json`, and in `src/exp7_common.H`, and are disjoint from every seed
used anywhere else in this repository (exp1 1001–1005, exp2 2001–2005, exp3 3001–3003 and
3101, exp4 and exp6 4001–4010, the first holdout 7001–7010, the second holdout 8001–8010, the
Experiment 7 selftest fixtures 7501–7505). Experiment 6 derived five of its performance seeds
implicitly as `4001 + block`, which is why 4006–4010 appear in its manifest and in no
declaration; the second block above exists so that no replicate in Experiment 7 is derived
rather than declared.

**How this block was drawn.** Not by choice. The rule is: take the highest seed declared
anywhere in this repository — 8010, the second holdout's performance block — round up to the
next multiple of 1000, which is 9000, and take the next ten integers. The rule admits exactly
one answer, so the block is a consequence of the repository's state rather than a selection
made after a result was seen, and it remains disjoint from anything a future experiment adds
below it. `make selftest` asserts that 7001–7010 and 8001–8010 appear in no block in use, so a
rerun on a spent block fails before it writes a byte.

**Configuration streams in P1 and P5.** Under 3.0.0 the replicate seed is not the engine seed.
The configuration `(precision, κ)` of replicate `base` runs on

    stream_seed = base + 10000 × (13 × precision_index + kappa_index)

with `precision_index` 0 for `double` and 1 for `float`, and `kappa_index` the 0-based position
in the ladder of §4. §2.7.2 is why: the 26 configurations of a replicate shared one stream, and
family F2's global rule is a count statistic whose null is built over independent units. The
formula is declared, is a formula rather than a table, and admits exactly one answer per
configuration; every value it produces is congruent to a production seed modulo 10000 and so
can collide with no block declared anywhere in this repository. Each P1 and P5 row records both
numbers — `seed`, the replicate, which remains the unit the analysis groups and cluster-
bootstraps by, and `stream_seed`, the engine the configuration actually ran on — so the
independence the null assumes is auditable from the raw rows rather than asserted here.

**P2's native replica layer follows P1's**, with the same configuration stream and the same
`stream_seed` recorded, because the schema predicts that the two reproduce each other exactly
and that prediction is a live check (`cross_phase_native_agrees`) rather than a description.
**P2's paired layer, P3, P4 and P6 keep the replicate seed**, for the reason in §2.7.2: their
families decide by Holm, by Simes or by an exact per-attempt identity, all valid under
arbitrary dependence between configurations.

**7001–7010 and 8001–8010 are spent.** They carry the two preserved NO-GO holdouts. No result
may be recomputed on either; they may now serve only as preserved failure evidence and as
material for diagnosing a defect, which is what §2.6 and §2.7 used them for. If this holdout
fails too, the recovery path is in §8.

## 4. Test matrix

Frozen sizes. `N` is intended attempts, never returned samples, except where stated.

| phase | methods | configurations | replication |
|---|---|---|---|
| **P1 scalar** | LEGACY, CANDIDATE | kappa ∈ {0.5001, 0.501, 0.505, 0.51, 0.55, 0.60, 0.75, 1.0, 1.25, 1.49, 1.5, 2, 5}; float and double | 10^6 × 5 seeds |
| **P2 mechanism** | LEGACY, CANDIDATE, oracle | same ladder, double and float, isotropic, ẑ, uncapped | 10^6 × 5 seeds |
| **P3 conditioning** | LEGACY, CANDIDATE, oracle | double kappa ∈ {0.501, 0.505, 0.51, 0.75}; float kappa ∈ {0.55, 0.60, 0.75} | 10^6 × 5 seeds |
| **P4 loader** | LEGACY, CANDIDATE | C0–C6 of §4.1 | uncapped 10^5 intended × 5 seeds; capped 10^5 **returned** × 5 seeds, all attempts retained |
| **P5 portability** | CANDIDATE | full P1 ladder on every available environment | 10^6 × 5 seeds |
| **P6 performance** | LEGACY, CANDIDATE | B0–B4 of §4.2 | 10 timed blocks ≥ 2 s and ≥ 10^6 attempts each, seeds 9006–9010 |

The ladder adds **1.25 and 1.49** to Experiment 6's. They exist because that is where the
rejected primitive's cost diverges, and a ladder that steps over the only region where a
candidate misbehaves is not a test of it. 1.49 and 1.5 straddle `α = 1`.

### 4.1 Complete-loader cases

| ID | precision | kappa | θ∥/θ⊥ | b̂ | cap | role |
|---|---|---:|---:|---|---|---|
| C0 | double | 2 | 2 | ẑ | none | benign reference |
| C1 | double | 0.51 | 2 | (1,2,3)/√14 | none | primary rescued case |
| C2 | double | 0.51 | 1/2 | (−2,1,2)/3 | none | inverse anisotropy + rotation |
| C3 | float | 0.55 | 2 | (1,2,3)/√14 | none | single-precision boundary |
| C4 | double | 0.505 | 1 | ẑ | none | honest-overflow semantics |
| C5 | double | 0.51 | 2 | (1,2,3)/√14 | 5 | strong bounded target |
| C6 | double | 0.51 | 2 | (1,2,3)/√14 | 20 | weak bounded target |

### 4.2 Benchmark cases

B0 double kappa = 0.51 uncapped; B1 float kappa = 0.55 uncapped; B2 double kappa = 0.75
uncapped; B3 double kappa = 0.51 cap 5; B4 double kappa = 0.51 cap 20.

## 5. Statistical protocol

Total type-I budget **α = 0.05**, partitioned across seven families. Each family has one
global statistic; a family is the unit of decision. "Zero rejections anywhere" is not used as
a rule: over K families at level α it is a level-(1 − (1 − α)^K) procedure in the failure
direction, which for Experiment 6's E1 alone was 0.57 per candidate.

| family | contents | α_F | global rule |
|---|---|---:|---|
| **F1** radial law | Anderson–Darling, KS and Cramér–von Mises on `Z = −log I_W(a, 3/2) ~ Exp(1)`, per configuration (see §2.5: the Beta route on `W` is withdrawn) | 0.010 | Simes global test over all configurations; on rejection, Holm within the family names the configuration |
| **F2** quantile coverage | every order-statistic interval, all kappa × precision × seed × level | 0.005 | Upper-tail test on the miss count against the 130-fold convolution of the per-configuration miss-count law of §2.7.2, whose per-level marginals are the **analytically computed** achieved coverages. Report count, μ and p; never a ratio |
| **F3** quantile direction | standardized signed quantile error per (kappa, precision, p), Stouffer-combined over seeds, Holm over cells | 0.010 | null calibrated by parametric Monte Carlo from the exact law at the same n and a, ≥ 2000 replicates, because the null is not symmetric at small a |
| **F4** upper-tail mass | exceedance counts above `z₀ = −log q₀` for `q₀ ∈ {1e−2, 1e−3, 1e−4}` (exact binomial), and KS of the excesses against Exp(1) | 0.010 | Holm within the family, Simes globally |
| **F5** loader battery | direction uniformity, rank-based independence, frame invariance, cap law, anisotropy, and — per §2.5, with the nulls corrected in §2.6 — an upper-tail count and an upper-tail excess test on each uncapped cell at each `q₀` | 0.010 | Holm over **all** cells and tests jointly, not per cell |
| **F6** negative controls | measured power curves | — | reversed: require power ≥ 0.90 at the pre-registered minimum detectable effect |
| **F7** portability | pairwise environment comparison | 0.005 | equivalence, not identity — see §5.3 |

### 5.1 Why F2 replaces "every interval covers", and what its null is

Each interval is an exact binomial order-statistic bracket at 95 % Bonferroni-corrected over
the five levels, i.e. nominal 99 %. Its achieved coverage is a closed-form function of `n`,
`p` and the correction, so the marginal miss probability of every interval is known before any
data exist.

The **joint** law is the part 2.0.0 got wrong, and §2.7.2 says how. The 650 intervals are 130
configurations of five levels. Under §3 the configurations run on independent streams, so the
130 group counts are independent; the five levels inside a configuration are order statistics
of one sample and are not. The null of the total is therefore the 130-fold convolution of the
per-configuration miss-count law, and that law is computed rather than assumed: an interval at
level `p` covers iff the count of sample points at or below the quantile lies in `[lo+1, hi]`,
and under the null the vector of those five counts is the running sum of a multinomial over
the six intervals the five levels cut. `config/protocol.json` stores the per-cell achieved
coverages, the per-configuration law, the resulting critical value and the stability checks;
the analysis reads them, re-convolves, and refuses to run if the critical value does not come
back out.

**F2 is also declared non-decisive at small shape.** The brackets are on `log R` and their
width scales as `1/a`: at `kappa = 0.5001` the median width runs from 57 natural-log units at
`p = 0.5` to 8951 at `p = 0.9999`. An interval that admits a multiplicative error of e^8951
is not evidence. Every quantile cell is therefore classified in advance by its width:
a cell wider than **1 natural-log unit is marked non-informative**, counts neither for nor
against, and is published in an informativeness map. The kappa values those cells cover are
certified by F1, F3 and F4, which retain power at any shape. This is the honest statement of
where this experiment can certify anything, and it is fixed before the run rather than
discovered after it.

### 5.2 Intervals and pooling

- Failure probabilities: two-sided 95 % Clopper–Pearson; a zero count is reported as the
  one-sided 95 % upper limit `1 − 0.05^(1/N)`, and the column records which of the two each
  entry is, since they are not the same kind of bound.
- Seed-level results are preserved and a per-configuration seed-homogeneity test (exact
  multinomial) is a required diagnostic with its own reported p-value. Experiment 6's
  "seed-stratified" bootstrap had zero between-seed variance by construction and so could not
  detect stream or portability defects; Experiment 7 uses a **cluster bootstrap over seeds**
  (resample seeds with replacement, then within seeds), 10 000 resamples, 99 % percentile
  intervals, with a distinct recorded RNG seed per interval.
- No statistic is computed on a subsample filtered by anything correlated with the radius
  unless the cell is labelled **conditional** and reported with its loss fraction. A
  conditional cell cannot by itself support the fidelity claim.

### 5.3 F7: equivalence, not identity

Two claims, tested differently:

- **Same architecture, different standard library** — the candidate's stream is a function of
  the engine alone, so the prediction is **bitwise equality** of every counter and of the
  returned vectors. Any difference fails F7 outright; there is no tolerance, because the claim
  is exact. This is a sharper test than Experiment 6 could run.
- **Different architecture** — `log`/`exp`/`pow` differ between libm implementations, so the
  prediction is equality of the *rates*. Pre-registered equivalence margin: the log rate ratio
  lies within **±0.15** (two one-sided tests at α_F/2 each, per configuration, Holm-corrected
  across configurations). Zero-count cells are compared by an exact conditional test against
  the same margin.

Emulation is corroborating evidence and never closure. Rosetta 2 and Docker/Rosetta results
are recorded with `execution = translated` and are excluded from the F7 decision. An
environment that is unavailable leaves G4 **open** with the exact command required to run it
elsewhere; it is never recorded as a pass, and never replaced by a zero.

### 5.4 Negative controls, calibrated at the effect that matters

Experiment 6 injected a defect into 50 % of the sample and tested a cap-5 sample against the
uncapped law — effects far larger than anything at issue — and reported a binary. Experiment 7
measures power:

- **NC1 radius–direction coupling**, injected at the loss fractions actually observed in the
  cells being certified (1e−3 and 1e−4), 200 injections each, requiring power ≥ 0.90.
- **NC2 capped against uncapped** at the *weakest* cap, λ = 20, not the strongest.
- **NC3 survivor conditioning**, respecified in §2.6.3. **NC3a** injects the cell's own honest,
  direction-dependent censoring, correctly reported — a *correct* loader — and the battery must
  not reject it; the row is a level and carries no threshold. **NC3b** injects that plus a
  further `q` fraction of the draws the type would have allowed back, removed silently, and
  requires power ≥ 0.90 at the pre-registered `F6.nc3_excess_loss_fractions`. 1.3.0 injected
  honest overflow alone and called detecting it power; that is what §2.6.3 corrects. Experiment
  6 had no such control at all.

A missing control is a gate failure, not a vacuous pass.

## 6. Gate decision rules

Every rule below is computed. None is satisfied by a prose string, by the existence of a file,
or by a number being finite.

| gate | passes iff |
|---|---|
| **G0** primary-source audit | the `log q` evaluator agrees with an arbitrary-precision incomplete beta to < 1e−10 on the declared grid; the attribution checklist in `config/protocol.json` — equation, stated domain, output scale, limitation, DOI — is machine-checked present for every cited method; and the pre-registered acceptance-rate predictor for any rejection-based primitive is named **in this document before the run** |
| **G1** scalar correctness | F1, F2, F3 and F4 all pass on in-domain configurations including a benign control that exercises the candidate itself; and the candidate's avoidable loss is **exactly zero** at every configuration, with no tolerance proportional to a data-dependent count |
| **G2** mechanism closure | the accounting identity `N = finite + avoidable + honest` holds exactly per seed and configuration; oracle disagreements and conversion failures are zero; the candidate's honest count equals the oracle floor exactly, or each discrepancy is individually adjudicated and listed |
| **G3** complete-loader fidelity | F5 passes jointly; F6 power thresholds are met and the expected number of control cells exists; every conditional cell is labelled with its loss fraction |
| **G4** portability | within the scope §2.6.4 declares: every environment in `gates.G4.supported_environments` completed **natively**, F7 passed, and at least one cross-standard-library comparison was actually made and was **bitwise** equal, with no tolerance. A supported environment that did not run FAILS the gate. A cross-architecture comparison is not required — the claim is withdrawn — but one that is filed and **disagrees** still fails |
| **G5** operational viability | the cluster-bootstrap median of the candidate/LEGACY time per returned sample is below the pre-registered bound of **2.0×** with its 95 % interval, same-seed reproducibility holds, and the RNG-stream break is documented |
| **G6** reproducible artifact | `make verify` exits zero on both manifests **from a fresh extraction**; the production run's dependency set is clean and the manifest records it consistently; `make reverify` reproduces every derived artifact byte for byte; the P1 baseline comparison resolved; and an exact-version archive identifier is present |

**Verdict.** `GO` iff G0–G6 all pass. `PARTIAL` iff G0–G3 and G6 pass but some production gate
does not. `NO-GO` otherwise. The verdict is the first line of `results/analysis_report.md`.

## 7. Expected false-failure rate

Computed from the frozen family structure before the run, and published with the result so
that a failure can be read against a stated rate rather than argued about afterwards. With
α = 0.05 partitioned as in §5, the probability that a fully correct candidate fails at least
one family is **at most 0.05** by construction, since each family's global statistic is a
level-α_F test and the α_F sum to 0.05. `config/protocol.json` records each α_F and the
resulting bound; the analysis recomputes and prints it.

## 8. What counts as a holdout failure, and what may follow

Any pre-registered family or gate rule failing at its pre-registered level on seeds 9001–9010
is a holdout failure. No re-pooling, no post-hoc exclusion of out-of-domain or near-limit
cells, no change to the seeds, the levels, the family membership, or the informativeness
threshold after the holdout is read.

The only legitimate path after a failure is: identify a concrete implementation defect, fix
it, freeze a new implementation hash **and a new protocol document**, draw a further disjoint
seed block, and rerun. The failed run is preserved in the repository and reported. A rerun
undertaken to obtain a more favourable result, without an identified defect, is not permitted
and would invalidate everything above.

## 9. Claim boundary

Fixed by the independent novelty audit in `docs/revision/literature/`, not by this result.
Supportable if every gate passes:

> Attributed integration of established log-scale Gamma generation into a complete anisotropic
> bi-Kappa Gamma-ratio loader, with mechanism-resolved finite-precision, tail, bounded-law,
> three-dimensional and performance validation, and with cross-standard-library reproducibility
> established bitwise on the one architecture tested.

Not supportable, and not to be written anywhere: discovery of the small-shape Gamma underflow
or of the low-parameter denominator hazard; a new or first log-domain Gamma generator; a new
Gamma-ratio, Beta-prime, Student-t or rejection Kappa sampler; that the log representation is
exact; any universal mathematical lower bound on kappa; **and, under 2.0.0, any claim about
behaviour on an architecture other than the one tested.** Cross-architecture agreement is
untested here, is withdrawn from the claim, and must be stated as a limitation wherever the
result is reported — not as an expectation, and not silently omitted. Rates are reported as observed counts
with intervals, and a configuration with no observed failure is reported as "no failures
observed in N draws under the tested configuration" with its one-sided upper bound — never as
"reliable down to".
