# Experiment 7 — pre-registered confirmatory protocol

**Status: frozen before any Experiment 7 data existed.** This document and
`config/protocol.json` are committed in a source-only commit. The analysis refuses to run
unless the SHA-256 of both files matches the values recorded in the run manifest, and every
decision rule below is evaluated by code that reads `config/protocol.json` rather than by a
constant written into the analysis.

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

The candidate is the released loader `cpp/bi_kappa_distribution.H` at version 2.0.0, whose
radius is built in the log domain. The comparator is the same file at version 1.0.0, vendored
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
- Within version 2.0.0 the stream is a function of the engine and the parameters alone, so
  libc++ and libstdc++ builds agree bit for bit. Across architectures they do not: `log`,
  `exp` and `pow` are not correctly rounded and the arm64 and x86_64 libm implementations
  differ, so cross-architecture agreement is statistical, not bitwise. G4 tests each of these
  two claims with the rule appropriate to it.
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
2. **Audit coverage is a declared rate, not an implementation detail.** Experiment 6's
   stratified audit sampled bulk failures at 1 in 4957 where the plan says "send every
   public-path failure"; with zero disagreements in 7039 of them, the one-sided 95 % bound
   still admits about 14 900 misclassified failures among 34.9 million. Every public-path
   failure, every near-limit draw, every method disagreement and every FINITE_BUT_WRONG
   candidate draw is now adjudicated in full; the remaining strata are sampled at rates fixed
   in `config/protocol.json`, and **both the audited and the total count of every stratum are
   reported**.
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

## 3. Seeds

**Production seeds: 7001–7005. Performance-block seeds: 7006–7010.** Both blocks are declared
here, in `config/protocol.json`, and in `src/exp7_common.H`, and are disjoint from every seed
used anywhere else in this repository (exp1 1001–1005, exp2 2001–2005, exp3 3001–3003 and
3101, exp4 and exp6 4001–4010). Experiment 6 derived five of its performance seeds implicitly
as `4001 + block`, which is why 4006–4010 appear in its manifest and in no declaration; the
second block above exists so that no seed in Experiment 7 is derived rather than declared.

No result may be recomputed on a different seed block. If the holdout fails, the recovery path
is in §8.

## 4. Test matrix

Frozen sizes. `N` is intended attempts, never returned samples, except where stated.

| phase | methods | configurations | replication |
|---|---|---|---|
| **P1 scalar** | LEGACY, CANDIDATE | kappa ∈ {0.5001, 0.501, 0.505, 0.51, 0.55, 0.60, 0.75, 1.0, 1.25, 1.49, 1.5, 2, 5}; float and double | 10^6 × 5 seeds |
| **P2 mechanism** | LEGACY, CANDIDATE, oracle | same ladder, double and float, isotropic, ẑ, uncapped | 10^6 × 5 seeds |
| **P3 conditioning** | LEGACY, CANDIDATE, oracle | double kappa ∈ {0.501, 0.505, 0.51, 0.75}; float kappa ∈ {0.55, 0.60, 0.75} | 10^6 × 5 seeds |
| **P4 loader** | LEGACY, CANDIDATE | C0–C6 of §4.1 | uncapped 10^5 intended × 5 seeds; capped 10^5 **returned** × 5 seeds, all attempts retained |
| **P5 portability** | CANDIDATE | full P1 ladder on every available environment | 10^6 × 5 seeds |
| **P6 performance** | LEGACY, CANDIDATE | B0–B4 of §4.2 | 10 timed blocks ≥ 2 s and ≥ 10^6 attempts each, seeds 7006–7010 |

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
| **F1** radial law | Anderson–Darling, KS and Cramér–von Mises on `Z = −log I_W(a, 3/2) ~ Exp(1)`, per configuration; the Beta test on `W` additionally where *every* draw resolves | 0.010 | Simes global test over all configurations; on rejection, Holm within the family names the configuration |
| **F2** quantile coverage | every order-statistic interval, all kappa × precision × seed × level | 0.005 | Poisson–binomial upper-tail test on the miss count against Σ(1 − c_i), where each `c_i` is the **analytically computed** achieved coverage. Report count, μ and p; never a ratio |
| **F3** quantile direction | standardized signed quantile error per (kappa, precision, p), Stouffer-combined over seeds, Holm over cells | 0.010 | null calibrated by parametric Monte Carlo from the exact law at the same n and a, ≥ 2000 replicates, because the null is not symmetric at small a |
| **F4** upper-tail mass | exceedance counts above `z₀ = −log q₀` for `q₀ ∈ {1e−2, 1e−3, 1e−4}` (exact binomial), and KS of the excesses against Exp(1) | 0.010 | Holm within the family, Simes globally |
| **F5** loader battery | direction uniformity, rank-based independence, frame invariance, cap law, anisotropy | 0.010 | Holm over **all** cells and tests jointly, not per cell |
| **F6** negative controls | measured power curves | — | reversed: require power ≥ 0.90 at the pre-registered minimum detectable effect |
| **F7** portability | pairwise environment comparison | 0.005 | equivalence, not identity — see §5.3 |

### 5.1 Why F2 replaces "every interval covers"

Each interval is an exact binomial order-statistic bracket at 95 % Bonferroni-corrected over
the five levels, i.e. nominal 99 %. Its achieved coverage is a closed-form function of `n`,
`p` and the correction, so the null distribution of the miss count is known before any data
exist. `config/protocol.json` stores the per-cell achieved coverages and the resulting
critical value; the analysis reads them.

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
- **NC3 survivor conditioning**: an uncapped sample with its non-finite attempts removed, at
  the loss fraction of the candidate's worst cell. This measures directly whether the battery
  can see the effect it is being used to certify absent. Experiment 6 had no such control.

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
| **G4** portability | F7 passes on every **native** environment; bitwise equality holds across standard libraries on each architecture; unavailable environments leave the gate open with an exact command |
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

Any pre-registered family or gate rule failing at its pre-registered level on seeds 7001–7010
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
> three-dimensional, portability and performance validation.

Not supportable, and not to be written anywhere: discovery of the small-shape Gamma underflow
or of the low-parameter denominator hazard; a new or first log-domain Gamma generator; a new
Gamma-ratio, Beta-prime, Student-t or rejection Kappa sampler; that the log representation is
exact; any universal mathematical lower bound on kappa. Rates are reported as observed counts
with intervals, and a configuration with no observed failure is reported as "no failures
observed in N draws under the tested configuration" with its one-sided upper bound — never as
"reliable down to".
