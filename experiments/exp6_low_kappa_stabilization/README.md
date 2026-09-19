# Experiment 6 — finite-precision validity and log-domain stabilization near κ = 1/2

The released bi-Kappa loader builds its radius from a Gamma ratio,

```
R = sqrt(X1) / sqrt(X2),    X1 ~ Ga(3/2, 1),    X2 ~ Ga(a, 1),    a = κ − 1/2,
```

so as κ → 1/2 the denominator shape goes to zero, `X2` underflows, and draws are lost. This
experiment asks three separate questions about that and keeps the answers separate:

1. **How much of the loss is avoidable?** A draw is *avoidably* lost when the vector the
   loader intended to return is perfectly representable and the arithmetic threw it away.
   It is *honestly* lost when the intended vector is not representable at all — no
   formation can rescue that, and pretending otherwise would be a claim about the
   floating-point type, not about the sampler.
2. **Does losing those draws matter?** Failure is concentrated in the tail the loader
   exists to reproduce, so discarding failures and testing the survivors tests
   `Pr(V ∈ A | success)` rather than the target `Pr(V ∈ A)`. The experiment measures the
   conditioning rather than assuming it away.
3. **Can an established log-domain Gamma method remove the avoidable part inside a complete
   anisotropic loader, and at what cost?**

## Scope

This experiment does **not** change the released default in `cpp/bi_kappa_distribution.H`,
does not change the public API, does not publish a release, and does not overwrite
Experiment 4. The log path is a research prototype confined to `src/`. Experiment 4 remains
the historical comparator; nothing here writes to it.

## What is under test

| name | radius | role |
|---|---|---|
| `QF` | `sqrt(x1/x2)` | the pre-fix formation; the intermediate quotient can overflow even when its square root is representable |
| `SPLIT` | `sqrt(x1)/sqrt(x2)` | **the released formation**, and the baseline in every comparison |
| `LOG` | `log R = (log X1 − log X2)/2`, carried in the log domain to the final component test | the candidate |

`LOG` needs a small-shape Gamma variate on the logarithmic scale. Two attributed primitives
were implemented from primary sources and compared in the pilot:

* **LOG-ID** — `log X = log Y + log(U)/a` with `Y ~ Ga(a+1,1)`, `U ~ U(0,1)`: the classical
  shape-boosting identity of Ahrens and Dieter (1974), which is how Marsaglia–Tsang style
  generators reach shapes below one, written on the log scale. Valid for every `a > 0`.
* **LOG-PUB** — the small-shape acceptance–rejection generator of Liu, Martin and Syring
  (2017, *Comput. Stat.* **32**, 1767–1775, doi:10.1007/s00180-016-0692-0), Sect. 3, p. 1770,
  which returns `log Y` by design. Stated domain `0 < a < 1`; outside it the implementation
  falls back to `std::gamma_distribution` and flags the draw `out_of_domain`, so nothing is
  ever attributed to the published method outside the range its authors state.

Xi, Tan and Liu (2013, *J. Stat. Softw.* **55**(4), doi:10.18637/jss.v055.i04) was read in
full and recorded in `config/pilot.json` as the considered alternative: its Algorithm 1 is
valid for all `a > 0` and delivers log-scale output below a hard-coded `a < 0.01` switch,
using piecewise-linear bounding tables fitted at efficiency `R = 0.9`.

**Nothing here is a new sampling method.** Log-domain Gamma generation, the Gamma-ratio
Kappa construction, and the low-parameter denominator hazard are all established. What is
under test is the integration of an established primitive into a complete anisotropic
bi-Kappa loader, and the validation of that loader's three-dimensional output, tail, cap
conditioning, portability and cost.

### A note on the published acceptance rate

Liu, Martin and Syring give the acceptance rate of their scheme as `r(α) = {1 + w(α)}⁻¹`
(Eq. 2, p. 1771). Their target `h_α` is normalized (`c⁻¹ = Γ(α+1)`) while their envelope
`η_α` integrates to `c(1+w)`, so the ratio of areas — which is what a measured acceptance
rate converges to — is `Γ(α+1)/{1+w(α)}`. The two agree to leading order as `α → 0`, the
regime the paper addresses, and the paper's own prose identifies `1/(1+w)` as the
probability of *proposing* from the `Exp(1)` branch. The measurement follows the area ratio
to within sampling error at every shape tested, and the sampled law is correct either way
(`E[log X] = ψ(α)` is verified independently in `make selftest`). The implementation is used
exactly as transcribed; both rates are recorded in `results/scalar_validation.csv`.

## Two layers, never pooled

* **Paired diagnostic layer.** `QF`, `SPLIT` and the log-domain propagation receive the
  *same* declared primitives `(x1, y, u, cosθ, φ)`, so one mathematical attempt can be
  followed through every formation and "was this draw recoverable?" is decided per draw
  instead of by comparing histograms. This is a mechanism experiment. It is **not** a
  draw-by-draw reconstruction of any standard library's Gamma generator: it declares its own
  primitives. Because the layer requires shared primitives, its `LOG` column is always
  LOG-ID arithmetic, and the rows say so, even when LOG-PUB is the selected native primitive.
* **Native layer.** The released loader runs against its own `std::gamma_distribution`, and
  the log loader against its own primitive. This is what a user receives. It cannot separate
  avoidable from honest loss on a draw whose denominator materialized as zero, because the
  intended value is then unknowable from the output — which is precisely what the paired
  layer is for.

## Why the log path can rescue a rotated draw

The released loader forms the local field-aligned vector and then rotates it. An orthogonal
rotation can spread one over-large local component across three representable global ones,
so a draw whose *returned* components are all representable can still be lost in the
intermediate. Writing `V = R g` with `g = Q(b̂) diag(√κ θ) n` order-unity, and deciding
representability from `log|V_j| = log R + log|g_j|` before exponentiating, removes that mode
as well. It is rare — a few draws in 10⁶ at κ = 0.505 — and it is counted separately as
`rotation_recoverable` rather than folded into the headline number.

## The oracle

`src/exp6_oracle.cpp` is a separate program, built as C++14 against
`boost::multiprecision::cpp_dec_float_100`. It exists to disagree with the probe: it reads
only the *primitives* of each audited attempt and re-derives every classification from
scratch at 100 decimal digits with its own code, sharing no arithmetic with
`src/exp6_loaders.H`. Keeping it separate is also what lets the probe and the sampling
headers stay strictly C++11 (`make cxx11-check` proves it), since Boost.Multiprecision
requires C++14.

The oracle is not decorative. During development it rejected 31 of 266 000 audited
attempts, all in single precision, and the cause was a defect in the probe rather than in
the oracle: the paired layer was forming the reference shape `kappa - 1/2` in double
instead of in the working precision. At float `kappa = 0.5001` that shifts `log X2` by
about 0.03, which is enough to move draws across the representability threshold and
mis-split avoidable from honest loss. The defect was fixed and the affected phases rerun.

The oracle also re-derives each formation's radius exactly, which is how the experiment can
say something the counters cannot: on draws whose denominator landed in the subnormal range
— the draws the split form is closest to losing — the split form's *surviving* radius is
already degraded, while the log path never materializes the denominator at all.

### Stratified audit sampling

The oracle costs about 1.1 ms per record. Auditing every public-path failure at the frozen
sample sizes would be several CPU-hours spent almost entirely on one trivially classified
mode: the denominator materialized as exactly zero, thousands of log units from any
threshold. The audit is therefore stratified, with the strata and rates fixed in
`src/exp6_common.H`, and **both the audited and the total count of every stratum are
reported** in `failure_envelope.csv` so the coverage is visible rather than implied:

| stratum | audited |
|---|---|
| within 2 natural-log units of the type's overflow threshold | all |
| the formations reached different terminal categories | 1 in 100 |
| the denominator landed in the subnormal range | 1 in 10 |
| every formation failed, far from any threshold | 1 in 10 000 |
| unbiased sample of all attempts | 1 in 10 000 |

Every decision-relevant draw — anything where the classification could genuinely go either
way — is audited in full. A single disagreement in any stratum fails gate G2.

## Reproducing

```bash
make selftest                       # accounting identities, replica-vs-class equality, fixtures
make preflight                      # environment + provenance; refuses a dirty production run
make pilot && make analyze          # E1: select and freeze the log primitive
make baseline mechanism conditioning loader portability performance
make oracle                         # adjudicate the audit streams at 100 decimal digits
make analyze && make figures
make checksums && make verify
```

Every simulation target accepts `SMOKE=1`, which shrinks the run to 1000 attempts on one
seed and diverts its output to `raw/smoke/`, where `analyze.py` ignores it unless `--smoke`
is passed. `ALLOW_DIRTY_DEV=1 make preflight` permits a development run against an
uncommitted dependency; the manifest is then marked `exploratory` and its output must not
feed a manuscript figure or table.

`make preflight` defines "clean" narrowly and on purpose: what must be unmodified relative
to `HEAD` is the dependency set — the released header, the experiment sources, the
GNUmakefile, the configs and the analysis scripts. Unrelated edits elsewhere in the working
tree do not change what this experiment computes, and the plan requires them to be
preserved. The overall repository state is recorded either way in `git.repo_dirty`.

## Environment

| tag | compiler | standard library |
|---|---|---|
| `libcxx` | Apple clang | libc++ |
| `libstdcxx` | Homebrew GCC | libstdc++ |

arm64 (Apple silicon), `-Wall -Wextra -std=c++11 -O2 -ffp-contract=off`, RNG
`std::mt19937`. **Floating-point contraction is disabled explicitly.** Fused multiply-add
changes the rounding of exactly the products whose overflow is under study, and the compiler
applies it inconsistently across inlining contexts; with contraction on, the instrumented
replica stops reproducing the released class bit for bit. The rounding mode, subnormal
handling and exact overflow thresholds of each build are recorded by `./exp6_probe_*.exe env`
and stored in `raw/environment.json`.

No x86_64 host was available. That platform is recorded in `results/portability.csv` with
`available=false` and the exact command to run it elsewhere; emulation is not accepted as a
substitute, so gate G4 stays open rather than being quietly passed.

## Frozen statistics

Fixed in `exp6_stats.py` before any production number was looked at, and not retuned:

* Failure probabilities: two-sided 95 % Clopper–Pearson; a zero count is reported as the
  one-sided 95 % upper limit `1 − 0.05^(1/N)`, never as "zero probability". Seed-level
  results are preserved alongside the pooled ones.
* Radial law: `W = X₂/(X₁+X₂) ~ Beta(a, 3/2)`, computed as `expit(−2 log R)` and carried as
  `log W = −logaddexp(0, 2 log R)`. `R²` is never formed: it overflows for `R > 1.3e154`
  while `R` is representable, which would silently discard exactly the heavy-tail draws the
  diagnostic exists to test.
* At the smallest shapes `W` itself rounds to zero, so the diagnostic of record is
  `Z = −log I_W(a, 3/2) ~ Exp(1)`, computed from `log W`. The ordinary Beta test on `W` is
  run **only where every draw resolves** — not merely most of them. Testing the survivors of
  a partial underflow is testing a truncation, and it rejects the truncation rather than the
  sampler. The criterion is numerical and is applied before any statistic is computed.
* The two-piece `log q` evaluator (switch at `log W = −30`) is validated against an
  `mpmath` incomplete beta at 60 digits over ten natural-log units on each side of the
  switch, for every pilot shape, before any sampling result is inspected
  (`results/log_q_evaluator_validation.csv`).
* Familywise type-I error controlled at α = 0.01 with Holm correction over the predeclared
  radial, direction, independence, frame and cap-law tests of each configuration.
* Tail ratios and method differences: paired seed-stratified bootstrap, 10 000 resamples,
  99 % percentile intervals. Quantile intervals: exact binomial order-statistic brackets,
  Bonferroni-corrected within each configuration; a level the sample cannot resolve is
  marked unresolved rather than extrapolated.
* Passing means consistency at the tested sample size, not proof of exactness.

## Negative controls

The sample size is not trusted until the diagnostics are shown to detect a defect that is
actually there.

* **NC1** couples radius to direction inside a random half of the sample by sorting each and
  re-pairing, which leaves both one-dimensional marginals bit-identical and changes only the
  joint. A diagnostic that cannot see this is not powerful enough to certify independence.
* **NC2** tests the capped sample against the uncapped target.

Both rows pass only when the injected defect **is** detected.

## Artifacts

| path | tracked | note |
|---|---|---|
| `src/`, `GNUmakefile`, `analyze.py`, `make_figures.py`, `preflight.py`, `exp6_*.py`, `config/` | yes | the complete source of every number |
| `results/*.csv`, `*.md`, `*.json`, `*.tex` | yes | source data, tables and reports |
| `figures/*.pdf`, `*.svg` | yes | canonical figures, editable text |
| `figures/*.png` | no | 600-dpi previews, regenerable by `make figures` |
| `raw/manifest.csv`, `raw/environment.json`, `raw/raw_checksums.sha256`, `checksums.sha256` | yes | provenance |
| `raw/**/*.jsonl` | yes | per-configuration counters, ~2 MB; every rate in `results/` is derived from these |
| `raw/**/*.bin` | no | bulk per-draw output, several GB, regenerable; covered by `raw/raw_checksums.sha256` |
| `*.exe` | no | build output |

## Gates

`results/analysis_report.md` opens with exactly one of `GO`, `PARTIAL` or `NO-GO`, followed
by the G0–G6 table and links to the decisive source data. G0–G3 and G6 are what a manuscript
mitigation claim requires; all seven are what production adoption requires.

## Deviations from the plan, and why

1. **The audit is stratified, not exhaustive.** Rates and per-stratum totals are reported;
   see above.
2. **The oracle is a separate C++14 program** rather than an inline call inside the probe.
   Boost.Multiprecision requires C++14, and the probe must stay strictly C++11; the
   separation also makes the oracle genuinely independent rather than a shared code path.
3. **`-ffp-contract=off`** was added to the build, for the reason given under Environment.
4. **Seeds 4001–4005** are reused from Experiment 4 as the plan freezes them, rather than
   taking a fresh per-experiment block as the other experiments here do. That is deliberate:
   E0's job is to reproduce Experiment 4, and identical seeds make the comparison exact
   rather than statistical.
5. **Extra source files** beyond the plan's listing — `src/exp6_loaders.H`, `exp6_stats.py`,
   `exp6_io.py`, `preflight.py` — separate the formations, the frozen statistics, the record
   readers and the provenance check from the phase drivers.
6. **x86_64 is not covered**; G4 is left explicitly open.
