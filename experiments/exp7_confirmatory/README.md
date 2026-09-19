# Experiment 7 — pre-registered confirmatory test of the log-domain bi-Kappa loader

The released loader builds its radius from a Gamma ratio,

```
R = sqrt(X1) / sqrt(X2),    X1 ~ Ga(3/2, 1),    X2 ~ Ga(a, 1),    a = κ − 1/2,
```

so as κ → 1/2 the denominator shape goes to zero, `X2` underflows, and draws are lost.
Version 2.0.0 never forms `X2`: it carries `log X2 = log Y + log(U)/a` with `Y ~ Ga(a+1,1)`,
propagates `log R = (log X1 − log X2)/2`, builds the order-unity vector
`g = Q(b̂) diag(√κ θ) n` **first**, and decides each returned component's representability
from `log|V_j| = log R + log|g_j|` before exponentiating anything. Its direction `n` is drawn
by the rejection method of Marsaglia (1972) — two open-interval uniforms per try, retried
until the pair falls in the unit disc, then lifted onto the sphere — so the sampler calls no
library transcendental and consumes a variable number of uniforms per attempt. Experiment 6
measured
that change; it could not confirm it, for reasons that are properties of its design and are
set out in `PROTOCOL.md` §1. This experiment is the confirmatory test: a disjoint seed
block, decision rules frozen before any datum existed, and a different implementation.

`PROTOCOL.md` and `config/protocol.json` were committed in a source-only commit and are the
authority for every level, margin and critical value. Nothing in this directory carries an
acceptance threshold of its own.

## What is under test

| name | radius | what runs |
|---|---|---|
| **LEGACY** | `sqrt(X1)/sqrt(X2)`, `X2` materialized by `std::gamma_distribution` | `src/legacy/bi_kappa_distribution_v1.H` — the released header exactly as shipped at 1.0.0, vendored unmodified |
| **CANDIDATE** | `log R = (log X1 − log X2)/2`, carried to the final component test | `cpp/bi_kappa_distribution.H` at 2.0.0 |

`QF` (`sqrt(X1/X2)`, whose intermediate quotient overflows even where its square root is
representable) survives only as a **third diagnostic column of the paired layer**. It is not
a method, it never runs natively, and no claim rests on it; it costs one square root per
attempt and it is the only way to show that the 1.0.0 split form was itself already a repair.

Nothing here is a new sampling method. Log-domain Gamma generation, shape boosting, the
Gamma-ratio Kappa construction, the small-shape underflow and the low-parameter denominator
hazard are all prior art, and are attributed in `config/protocol.json`'s
`attribution_checklist` with equation, stated domain, output scale, limitation and DOI for
each. What is under test is the integration of an established primitive into a complete
anisotropic bi-Kappa loader, and the validation of that loader's three-dimensional output,
tail, cap conditioning, portability and cost.

`make check-legacy` re-derives the original 1.0.0 header from the vendored comparator — by
undoing the include-guard rename and the `namespace bikappa_v1` wrap, and nothing else — and
diffs it against `git show 0fc2c95:cpp/bi_kappa_distribution.H`. A comparator that has
drifted is not a comparator.

## Two layers, never pooled

* **Paired diagnostic layer.** LEGACY, CANDIDATE and QF receive the *same* declared
  primitives `(x1, y, u, cosθ, φ)`, so one mathematical attempt can be followed through
  every formation and "was this draw recoverable?" is decided per draw instead of by
  comparing histograms. What it compares is the three **radius** formations on one
  mathematical draw; the direction is a shared nuisance parameter, handed identically to
  every column. **This is a mechanism experiment. It is not a draw-by-draw reconstruction
  of any standard library's Gamma generator, nor of either release's direction sampling**:
  it declares its own primitives and says so — see "The paired layer's declared direction"
  below. Every one of the five is an exact function of the engine, so the two builds' audit
  streams have the same SHA-256 by construction rather than by result.

* **Native layer.** Each method runs its own real implementation end to end. This is what a
  user receives. It cannot separate avoidable from honest loss on a draw whose denominator
  materialized as zero, because the intended value is then unknowable from the output —
  which is precisely what the paired layer is for. Where an intermediate has to be
  observable the native layer uses a replica, and `make selftest` checks each replica
  against the released class it stands for, component by component, in both precisions and
  in both a rotated and an unrotated configuration.

### The paired layer's declared direction

The layer's direction primitive is `(cosθ, φ)`, as it was in Experiment 6, and it stays that
way now that release 2.0.0 has stopped drawing angles at all. The alternative — giving the
layer a Marsaglia disc-rejection direction so that it matches 2.0.0's native path — was
considered and rejected, for three reasons.

* **It is a nuisance parameter, not a treatment.** The layer exists to push one mathematical
  draw through three radius formations. All that is required of the direction is that every
  column receive the identical one, which either choice satisfies.
* **Neither choice makes the layer a reconstruction of both native paths.** Release 1.0.0
  draws `cosθ` and `φ` through `std::uniform_real_distribution`; release 2.0.0 draws no angle
  at all. A shared primitive must be declared either way, and `(cosθ, φ)` is the one that
  LEGACY's replica, its pre-rotation cap test and the 1.0.0 comparator are already defined in
  terms of.
* **It is the representation three readers are written against.** `AuditRecord`, the
  100-digit oracle and the analysis reader all take the direction as an angle pair. Switching
  it would change a binary record layout, force a schema-version bump and invalidate every
  existing audit stream — for a quantity no claim is about.

The cost is stated rather than hidden. `(cosθ, φ)` reaches the vector through `sin(φ)` and
`cos(φ)` of one argument, and clang fuses that pair into Darwin's `__sincos_stret` while GCC
calls the two separately. Measured directly on this host: clang's unit direction differs in
the last place of the `y` component on **0.20 %** of two million draws from the value a
standalone `sin` gives, which is GCC's value on every one of them. So the layer's *derived*
numbers are a function of the compiler as well as of the source text, while its *recorded*
primitives are not — the audit stream is exact arithmetic on the engine and hashes the same
either way.

This is not repairable by de-fusing the shared helper. LEGACY's replica has to reproduce
1.0.0 exactly, fusion included, and that replica is the experiment's negative control; a
helper that suppressed the fusion would break the thing it exists to reproduce.

The residual is bounded. A terminal category can turn on an ulp of direction only where a
component sits within an ulp of a type limit, which falls in the within-margin stratum —
audited at rate 1 and adjudicated at 100 digits from the primitives themselves. In the
reduced run of P2 taken here, all 78 `paired` rows agree across the two builds, and the
oracle reports no disagreement on either build's audit stream. The portability claim of
record is P5's, taken over the released class, and rests on nothing in this layer.

## The phases

| phase | what it measures | output |
|---|---|---|
| **P1** scalar | the radial law, both methods, the full ladder, float and double, 10⁶ × 5 seeds | order statistics at the frozen F2 brackets, exact AD/KS/CvM of `Z` against Exp(1), the `Z` ECDF on a fixed grid, exceedance counts and the excesses themselves, terminal categories, variate counts, timing, and a digest over every returned bit pattern |
| **P2** mechanism | the avoidable/honest decomposition, paired and native, uncapped, isotropic, ẑ | per-configuration counters, the accounting identity per seed, and the stratified audit stream |
| **P3** conditioning | whether losing draws matters: the joint law of the intended state and success | one 32-byte record per attempt |
| **P4** loader | the complete anisotropic, rotated, capped 3-D loader, cases C0–C6 | one record per returned sample, plus every attempt and rejection counter |
| **P5** portability | CANDIDATE on the full ladder, one tagged environment per build | a digest per configuration, to be compared for **exact equality** |
| **P6** performance | cases B0–B4, 10 timed blocks each ≥ 2 s and ≥ 10⁶ attempts, randomized method order, one untimed warm-up, seeds 7006–7010 | per-block time per returned sample, and same-seed reproducibility |

The ladder is κ ∈ {0.5001, 0.501, 0.505, 0.51, 0.55, 0.60, 0.75, 1.0, 1.25, 1.49, 1.5, 2, 5}.
**1.25 and 1.49 are new** relative to Experiment 6 and are load-bearing: they sample the
region where the rejected primitive's cost diverges, and 1.49 and 1.5 straddle `a = 1`. A
ladder that steps over the only region where a candidate misbehaves is not a test of it.

## Seeds

Production **7001–7005**, performance **7006–7010**. Both blocks are written out as explicit
vectors in `src/exp7_common.H`, and the probe compares them — and the kappa ladder — against
`config/protocol.json` at start-up and refuses to run on any difference. **No seed anywhere
in this experiment is derived by arithmetic from another.** Experiment 6 obtained five of
its performance seeds as `4001 + block`, which is why 4006–4010 appear in its manifest and
in no declaration.

P6 indexes the declared performance vector by `block mod 5`. The randomization of the
*method order within a block* is a separate, declared constant that never reaches a variate,
and the realized order is recorded in every row so that it can be checked rather than
asserted.

`make selftest` runs on a third declared vector, **7501–7505**, disjoint from both. Its
fixtures write no data, but several of its checks are gate predicates — "CANDIDATE has zero
avoidable loss at κ = 0.501 in double" is gate G1 — and running those on a production seed
before the holdout is read would be reading the holdout.

## What a returned value has to be to count as a success

Amendment 1.1.0 to the protocol adds a category Experiment 6 did not have. Its terminal
classifier read `method_finite ? finite : honest_overflow`, so a method that returned a
number for a draw that has no representation was scored as a success — and so was a method
whose surviving radius had lost half its significand, which its own oracle measured at up to
0.48 relative error once the denominator went subnormal.

A finite return is therefore a success here only when it is also accurate. A returned radius
that differs from the arbitrary-precision value by more than `2^-(digits/2)` — 1.05e−8 in
`double`, 2.44e−4 in `float` — is **FINITE_BUT_WRONG** and counts as a loss. The two ways a
loss can be dressed as a success are kept apart:

* `finite_but_wrong` — the intended draw *was* representable and the arithmetic degraded it.
  Avoidable: a formation that carries the radius in the log domain returns it accurately.
* `overflow_returned_finite` — the intended draw was **not** representable and a number came
  back anyway. Honest: no formation can rescue a draw that has no representation, so this is
  not evidence that anything was thrown away.

With those two, `N = finite + avoidable + honest` still holds exactly, and `accounting_ok`
reports it per seed and configuration.

Every accuracy statistic carries its precision and none is pooled across the two. Experiment
6's headline contrast — 0.48 against 7.8e−6 — put a `float` worst case next to a `float`
figure without labelling either, while its `double` worst cases were 0.46 and 1.1e−13; six
of the orders of magnitude in that contrast came from the type rather than from the method.

## The oracle

`src/exp7_oracle.cpp` is a separate program, built as C++14 against
`boost::multiprecision::cpp_dec_float_100`. It exists to disagree with the probe: it reads
only the *primitives* of each audited attempt and re-derives every classification from
scratch at 100 decimal digits with its own code, sharing no arithmetic with
`src/exp7_loaders.H` — not even the `Method` enum, which it restates. Keeping it separate is
also what lets the probe and the sampling headers stay strictly C++11 (`make cxx11-check`
proves it under both compilers with `-pedantic-errors`), since Boost.Multiprecision requires
C++14.

Two defects of Experiment 6's oracle are fixed. Its `make oracle` ended in `|| echo`, so a
disagreement printed a line and the target still succeeded; here the program's non-zero exit
propagates and fails the target. And it recorded the audit file's path but not its content
hash, so nothing tied an adjudication to the bytes adjudicated; every summary line, and the
header and footer of the disagreement stream, now carry the file's SHA-256. The disagreement
stream gets that header and footer whether or not anything disagreed, so an empty
adjudication is distinguishable from an oracle that never ran.

### Audit coverage is stated by margin, not by rate

Experiment 6 sampled its bulk failures at 1 in 4957 while its plan said "send every
public-path failure"; with zero disagreements in 7039 of them, the one-sided 95 % bound still
admitted about 14 900 misclassified failures among 34.9 million. Requiring *every*
public-path failure instead is both infeasible and beside the point, because almost all of
them are thousands of log units outside the type's range, where no arithmetic could have
returned the value.

Amendment 1.2.0 therefore states the rule by margin, and the rates arrive in the C++ from
`config/protocol.json` through a generated header rather than as literals:

| stratum | audited | why |
|---|---|---|
| within **20 natural-log units** of a type limit | all | the classification could go either way |
| the formations reached different terminal categories | all | the decomposition is exactly what is in question |
| CANDIDATE returned a finite but inaccurate radius | all | this count must be zero |
| CANDIDATE suffered avoidable loss | all | this count must be zero |
| subnormal or zero denominator with a representable target | all | the avoidable-loss candidates |
| a failure more than the margin beyond a type limit | 1 in 10³ | unambiguous — **and the margin is asserted in working precision for every one of them** |
| everything else | 1 in 10⁴ | an unbiased sample of the attempts in no other stratum |

Twenty log units is about four thousand times the worst disagreement ever measured between
the working-precision reference and the 100-digit oracle. Beyond it a misclassification is
impossible rather than merely improbable, and that is an argument a sampling rate cannot
make. Both the audited and the total count of every stratum are reported, and
`audit_margin_assertion_failures` — which must be zero — records the inequality being
evaluated rather than assumed. A single disagreement in any stratum fails gate G2.

Strata are assigned by priority, first match wins, so the counts partition the sample and
the last stratum is an unbiased sample of the remainder rather than of everything.

## Sample size against disk

Experiment 6's `raw/` reached 4.1 GB, almost all of it per-draw `.bin` dumps that no
statistic read twice. Here per-attempt records are kept only where the analysis genuinely
needs per-draw data:

| phase | what is kept | projected size, both builds |
|---|---|---|
| P1 | summaries, plus the excesses above `z0 = −log(1e−2)` in an 8-byte-per-record companion file | ≈ 40 MB |
| P2 | per-configuration counters, plus the audit stream | ≈ 0.8 GB |
| P3 | one 32-byte record per attempt — conditioning is a statement about a joint law, which no summary preserves | ≈ 2.2 GB |
| P4 | one 40-byte record per returned sample, plus all attempt and rejection counters | ≈ 1.0 GB |
| P5, P6 | summaries only | ≈ 1 MB |

**Projected total ≈ 4.1 GB**, and projected wall clock, measured by extrapolating a reduced
run on this host: about 8 minutes per build for P1–P6, and about **4.5 hours** for
`make oracle` over both builds' audit streams at roughly 1.1 ms per record. No phase
approaches 24 wall-clock hours or 100 GB.

For P1 the summary is sufficient for families F1–F4 and nothing more is needed:

* **F2 and F3** read order statistics, and the probe emits the sample value at exactly the
  index brackets `config/protocol.json` froze, plus the empirical quantile at each level.
* **F1** reads Anderson–Darling, Kolmogorov–Smirnov and Cramér–von Mises on `Z`. The probe
  holds the sorted sample in memory and computes all three exactly, with the same estimators
  and the same CDF clipping the analysis module uses, so the numbers agree to the last bit
  rather than approximately. The `Z` ECDF on a 501-point grid is emitted alongside as an
  independent check and for the figures.
* **F4** reads exceedance counts *and the excesses*, because the power study in
  `PROTOCOL.md` §2.4 makes F4 the only family with power against a 10⁻⁴ loss fraction. The
  excesses above the lowest protocol threshold nest the other two, so a single companion
  file of about 10⁴ records per configuration answers all three exactly.

A draw the method could not resolve is sorted to `+inf` rather than dropped. Failure here is
overflow, so an unresolved draw *is* the largest draw; dropping it would silently turn the
sample into its own survivors, which is the effect the experiment exists to detect. Every
cell reports `n_resolved`, `n_nonresolved` and `loss_fraction`, and is flagged `conditional`
when anything was lost.

## The manifest is written at run time

Experiment 6 rebuilt `raw/manifest.csv` from whatever was on disk when the analysis ran, so
`probe_sha256` was the hash of the binary present at analysis time and 165 rows carried the
wrong `cxxflags`. Here the probe writes one row per output file, from inside the process
that produced it, **after** the atomic `.part`-and-rename and after re-opening the file and
checking its header, record count and byte length against what was intended. A file that
fails its own contract stops the run rather than being recorded. The flags come from a
`-D` the compiler baked into the binary, and the probe hashes its own executable.

The atomic writer is Experiment 6's, unchanged: the file is built under a `.part` name, its
header is rewritten with the final record count, and only then is it renamed, so a truncated
job never leaves a file that looks complete.

## Reproducing

```bash
make cxx11-check                 # strict C++11 under both compilers, -pedantic-errors
make check-legacy                # the comparator is still the released 1.0.0 header
make protocol-check              # config/protocol.json is what make_protocol.py produces
make selftest                    # accounting identity, replica-vs-class equality, digests
make preflight                   # environment + provenance; refuses a dirty production run
make p1 p2 p3 p4 p5 p6
make oracle                      # adjudicate every audit stream at 100 decimal digits
make analyze && make figures
make checksums && make verify && make reverify
```

Every simulation target depends on `preflight`, so a dirty dependency set stops the run
before it produces data rather than after. In Experiment 6 `preflight` was a prerequisite of
nothing, and a dirty tree ran to completion with only the manifest saying so.
`ALLOW_DIRTY_DEV=1 make preflight` permits a development run against an uncommitted
dependency; the manifest is then marked exploratory and its output must not feed a
manuscript figure or table.

`make preflight` is the opening step of a run, not a status command: it rewrites
`raw/environment.json` with the commit and timestamp the run starts from, and that file is
covered by `checksums.sha256`. Running it again after `make checksums` invalidates the
bundle.

Every simulation target accepts `SMOKE=1`, which shrinks the run to 2000 attempts on one seed
and diverts **everything** it writes — data files and the manifest alike — under `raw/smoke/`.
`make SMOKE=1 oracle` adjudicates only the smoke streams. `make checksums` prunes
`raw/smoke` outright. In Experiment 6 smoke output leaked into `make oracle`, `make
checksums` and the production manifest; here nothing outside `raw/smoke/` is touched by a
smoke run, and no production target looks inside it. `SMOKE` is forwarded to `analyze` and
`figures` as `--smoke`.

## Environment

| tag | compiler | standard library |
|---|---|---|
| `libcxx` | Apple clang | libc++ |
| `libstdcxx` | Homebrew GCC | libstdc++ |

arm64 (Apple silicon), `-Wall -Wextra -std=c++11 -O2 -ffp-contract=off`, RNG `std::mt19937`.
**Floating-point contraction is disabled explicitly.** Fused multiply-add changes the
rounding of exactly the products whose overflow is under study, and the compiler applies it
inconsistently across inlining contexts; with contraction on, the instrumented replica stops
reproducing the released class bit for bit. The rounding mode, subnormal handling and exact
overflow thresholds of each build are recorded by `./exp7_probe_*.exe env` and stored in
`raw/environment.json`.

No x86_64 host is available here. `PROTOCOL.md` §5.3 requires that platform to be recorded
with the exact command to run it elsewhere, not replaced by a zero and not substituted by
emulation; a Rosetta or Docker/Rosetta run is recorded with `execution = translated` and is
excluded from the F7 decision. Pass `EXECUTION=translated` to tag such a run.

### Portability: what the digests are for, and what they already show

The candidate's *stream* is a function of the engine alone — the Gamma, normal and uniform
primitives live in the header rather than in `<random>` — so `PROTOCOL.md` §5.3 predicts
**bitwise equality** of every counter and of every returned vector between the two builds on
one architecture. Testing that with interval overlap would be a weaker claim than the one
being made, so each configuration carries a SHA-256 over the IEEE bit patterns of every
returned component followed by a canonical serialization of every counter, and an FNV-1a
alongside it over the identical byte sequence as a check on the hash itself. The same
digests are taken for LEGACY, where equality is predicted **not** to hold; that contrast is
part of the result.

That prediction did not hold when it was first run here, and the cause was in the sampler
rather than in the libraries: the direction was drawn from `cosθ` and `φ`, and clang fuses
the `sin(φ)`/`cos(φ)` pair into Darwin's `__sincos_stret`, whose sine differs from the
standalone `sin` by one ulp on roughly one argument in a thousand, while GCC calls the two
separately. Two builds of one header then returned different numbers for the same seed, and
the difference belonged to the compiler rather than to the sampler.

The released header now draws its direction by the rejection method of Marsaglia (1972),
which uses only `sqrt` — correctly rounded by IEEE-754 — and so calls no transcendental at
all. It averages `4/π = 1.273` tries per direction, and the probe reports the realized
uniform count rather than assuming one: a CANDIDATE attempt consumes `3.53` uniforms on
average against LEGACY's fixed 2, and P6 is where the cost of that is measured.
The cause was removed rather than documented: every full-loader digest over the κ ladder,
capped and uncapped, rotated and anisotropic, in both precisions, is now identical under
clang/libc++ and g++-15/libstdc++, where the two previously disagreed. A reduced P5 run on
this host compares 26 CANDIDATE configurations and **none** differs. LEGACY, which still
draws angles, still differs — that half of the contrast is unchanged and is what makes the
equality a result rather than a tautology.

## Frozen statistics

Fixed in `config/protocol.json` and `exp7_families.py` before any production number existed.
The full statement is `PROTOCOL.md` §5; what the probe contributes is:

* `Z = −log I_W(a, 3/2) ~ Exp(1)` as the radial diagnostic of record, computed from `log W`
  by a two-piece evaluator whose switch point is the protocol's `log_q_switch`. `R²` is
  never formed: it overflows for `R > 1.3e154` while `R` is representable, which would
  silently discard exactly the heavy-tail draws the diagnostic exists to test.
* Failure probabilities are reported as counts with their `n`, never as bare rates; the
  interval construction is the analysis's.
* Seed-level results are preserved in full — every row is one seed — so the cluster
  bootstrap over seeds has something to resample. Experiment 6's "seed-stratified" bootstrap
  had zero between-seed variance by construction and so could not detect stream or
  portability defects.

## Artifacts

| path | tracked | note |
|---|---|---|
| `src/`, `GNUmakefile`, `PROTOCOL.md`, `config/`, `preflight.py`, `analyze.py`, `make_figures.py`, `exp7_*.py` | yes | the complete source of every number |
| `src/exp7_protocol.H` | generated | written from `config/protocol.json` by `src/gen_protocol_header.py`; never edited |
| `results/*.csv`, `*.md`, `*.json` | yes | source data, tables and reports |
| `figures/*.pdf`, `*.svg` | yes | canonical figures, editable text |
| `figures/*.png` | no | previews, regenerable by `make figures` |
| `raw/manifest.csv`, `raw/environment.json`, `raw/raw_checksums.sha256`, `checksums.sha256` | yes | provenance |
| `raw/**/*.jsonl` | yes | per-configuration counters; every rate in `results/` derives from these |
| `raw/**/*.bin` | no | bulk per-draw output, several GB, regenerable; covered by `raw/raw_checksums.sha256` |
| `raw/smoke/**` | no | never production data |
| `*.exe` | no | build output |

## Gates

`results/analysis_report.md` opens with exactly one of `GO`, `PARTIAL` or `NO-GO`, followed
by the G0–G6 table of `PROTOCOL.md` §6 and links to the decisive source data.

## Deviations from the protocol

1. **None that change a decision rule, a seed, a level or a matrix entry.** The kappa ladder,
   both seed blocks, the phase sizes, the loader cases, the benchmark cases, the quantile
   levels, the tail thresholds, the accuracy thresholds and the audit strata rates are all
   read from `config/protocol.json`, and the probe aborts at start-up if the copies declared
   in `src/exp7_common.H` disagree with it.

2. **`QF` is emitted as a third paired-layer column.** `PROTOCOL.md` §4 names two methods.
   QF is not a third one: it never runs natively, every row that carries it is marked
   `diagnostic_only`, and no gate reads it. It is kept because it is nearly free and it is
   what shows that the 1.0.0 split form was already a repair of an earlier formation.

3. **The probe evaluates F1's three statistics itself.** They are functions of the sorted
   sample, which the probe already holds and the analysis would otherwise need 8 MB per
   configuration to reconstruct. The estimators and the CDF clipping are the analysis
   module's, written out in `results/schema.md`, and the `Z` ECDF grid is emitted so the
   analysis can check them rather than take them.

4. **The uniform audit sample is of the remainder.** Strata are assigned by priority, so the
   1-in-10⁴ stratum samples the attempts that fell in no other stratum rather than all
   attempts. Both counts are reported for every stratum, so the coverage is exact.

5. **`make selftest` uses its own declared seed block, 7501–7505.** `PROTOCOL.md` §3
   declares two blocks and requires that no seed be derived; it does not name a fixture
   block, and running the selftest's gate-predicate checks on 7001–7005 would read the
   holdout before the run. The third vector is declared in `src/exp7_common.H` alongside the
   other two, and the probe checks all three for pairwise disjointness at start-up.

6. **x86_64 is not covered here**; gate G4 stays explicitly open with the command to run it
   elsewhere, and emulation is not accepted as a substitute.
