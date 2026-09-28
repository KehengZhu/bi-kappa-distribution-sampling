# Experiment 4 — finite-precision losses of the radius calculation

The sampler forms the radius from two Gamma variates,

```
R = √X₁ / √X₂,    X₁ ~ Γ(3/2),    X₂ ~ Γ(α),    α = κ − 1/2.
```

As κ → 1/2 the shape α goes to zero. A floating-point `X₂` can then underflow to zero, and the
attempt returns a non-finite velocity even when every final velocity component would be
representable. This is an avoidable loss. Some attempts also have a final component larger than
the largest finite value of the output type. Such a loss is unavoidable, since no calculation
can return that component. Both kinds of failure occur at large radii, so they remove draws from
the high-energy tail.

This experiment measures both kinds of loss for two radius calculations, in single and double
precision:

- **Direct calculation**, which forms `X₂` and evaluates `√X₁/√X₂`, with both Gamma variates
  drawn by `std::gamma_distribution`. It is `cpp/bi_kappa_distribution.H` as committed at
  53dbd30 (2026-08-18), the last version of the header before the stabilized calculation
  replaced it. The data files label it `LEGACY`.
- **Stabilized calculation**, which never forms `X₂`. It obtains `log X₂ = log G + (log ξ)/α`
  with `G ~ Γ(α + 1)` and `ξ ~ U(0, 1)`, evaluates `log R = (log X₁ − log X₂)/2`, and keeps
  `log R` until the final components are formed. This is `cpp/bi_kappa_distribution.H` of
  version 3.0.0. The data files label it `CANDIDATE`.

A third column, `QF` (`√(X₁/X₂)`, whose intermediate quotient can overflow when its square root
would not), appears only in the paired comparison described below. It is a diagnostic, not a
sampler under test. Release 1.0.0 formed the radius in this quotient-first way.

The study was numbered 7 when it was designed, so its source and data files carry the prefix
`exp7`.

## Reproducing Fig. 2

Fig. 2 of the paper is `figures/fp1_failure_envelope.pdf`. `make_figures.py` draws it from three
CSV files in `results/`, and `paper/figures/make_manuscript_assets.py` copies it into the paper.

The release archive (the GitHub and Zenodo downloads) contains this part of the experiment
only. From this directory,

```bash
make figures
```

runs `uv run --project ../../python python make_figures.py`, which needs `uv` and the Python
project in `../../python`. It redraws Fig. 2 as `figures/fp1_failure_envelope.pdf` and `.svg`,
with a PNG preview that is not part of the release, and rewrites `figures/captions.md` and the
Fig. 2 entry of `figures/figure_manifest.json`. With the shipped CSV files all four shipped
figure files are reproduced byte for byte. The other four figures are drawn from results that
are in the repository only; in the archive they are not drawn, and `figure_manifest.json` keeps
their entries.

| file in the release archive | contents |
|---|---|
| `results/failure_envelope.csv` | the non-finite fraction of every setting, with its 95% Clopper–Pearson interval, split into avoidable and unavoidable losses; the markers and the dotted line of Fig. 2 |
| `results/honest_floor_curve.csv` | the closed-form fraction of draws with a component too large for the type, on a fine grid of κ; the grey region of Fig. 2 |
| `results/direct_zero_denominator_curve.csv` | the closed-form probability that the direct calculation's `X₂` rounds to zero; the dashed line of Fig. 2 |
| `results/source_data_README.md` | definitions of every column |
| `figures/fp1_failure_envelope.pdf`, `.svg` | Fig. 2 |
| `figures/captions.md` | captions of the five figures of this experiment |
| `figures/figure_manifest.json` | for each figure, its source CSV files with their SHA-256, and the SHA-256 of `make_figures.py` |
| `make_figures.py`, `GNUmakefile` | the plotting script and the make targets |

Everything else described below — the protocol, the C++ probe, the 100-digit recomputation, the
analysis code and the run records — is in the GitHub repository.

## Design

| | |
|---|---|
| κ | 0.5001, 0.501, 0.505, 0.51, 0.55, 0.6, 0.75, 1, 1.25, 1.49, 1.5, 2, 5 |
| precision | `float` and `double` |
| attempts | 10⁶ per seed and setting, seeds 11001–11005 (5 × 10⁶ per setting); `p2` 10⁷ per seed (5 × 10⁷ per setting); timing on seeds 11006–11010 |
| builds | Apple clang with libc++, and GCC 15 with libstdc++, both on arm64 (Apple silicon) |
| flags | `-Wall -Wextra -std=c++11 -O2 -ffp-contract=off`; RNG `std::mt19937` |

Floating-point contraction is disabled because a fused multiply-add changes the rounding of
exactly the products whose overflow is studied here.

Each attempt is classified as returned, avoidably lost, or unavoidably lost. A finite result
also counts as a loss when its radius differs from the exact value by more than `2^−(p/2)`,
where `p` is the number of significand bits: 1.05 × 10⁻⁸ in `double` and 2.44 × 10⁻⁴ in
`float`. The counts satisfy `N = returned + avoidable + unavoidable` in every setting and seed,
and the analysis checks this.

The comparison is made in two ways.

- **Paired.** All three calculations receive the same random variates for each attempt, so the
  outcome of one attempt can be compared across calculations. Fig. 2 uses this comparison.
- **Native.** Each implementation runs as released, drawing its own variates. This is what a
  user receives.

The paired comparison uses 5 × 10⁷ attempts per setting so that a setting with no failure has
a one-sided 95% upper bound of 6.0 × 10⁻⁸ on its failure fraction; with 5 × 10⁶ the bound would
be 6.0 × 10⁻⁷, above the fractions observed at neighbouring κ.

**Recomputation at 100 digits.** The probe records the random variates of selected attempts.
`src/exp7_oracle.cpp`, a separate program built on `boost::multiprecision::mpfr_float_100` (MPFR)
that shares no arithmetic with the probe, recomputes each recorded attempt from its variates and
classifies it again. The recorded attempts are:

- every attempt within 20 natural-log units of a limit of the output type;
- every attempt on which the calculations reach different outcomes;
- every inaccurate or avoidable loss of the stabilized calculation;
- every attempt whose `X₂` is subnormal or zero while the exact result is representable;
- 1 in 10³ of the failures further than 20 log units beyond a type limit;
- 1 in 10⁴ of the remaining attempts.

For every recorded attempt that a calculation lost avoidably, the recomputation also reports
whether the exact largest component lies within the rounding band of the overflow threshold
(see Results).

The experiment runs in six phases, each a make target.

| target | measures |
|---|---|
| `p1` | the radial distribution of both calculations over the full κ list: Anderson–Darling, KS and Cramér–von Mises tests of `Z = −log I_W(α, 3/2) ~ Exp(1)`, quantiles of `log R`, and exceedance counts in the upper tail |
| `p2` | returned, avoidable and unavoidable counts, paired and native; uncapped, isotropic, field along `ẑ` |
| `p3` | the probability that an attempt succeeds as a function of the radius it was meant to produce, and the resulting loss of tail mass |
| `p4` | the complete loader with anisotropy, rotation and optional cap, in seven configurations C0–C6 |
| `p5` | a SHA-256 digest of every returned value and counter, compared between the two builds |
| `p6` | time per returned sample, stabilized against direct, in five benchmark cases |

`PROTOCOL.md` specifies the design, every statistical test and every threshold, fixed before the
run. `config/protocol.json` holds the same settings in machine-readable form; it is generated by
`config/make_protocol.py`, and the C++ probe reads it through the generated header
`src/exp7_protocol.H`.

## Results

The committed results were produced with `cpp/bi_kappa_distribution.H` as committed at 0e7783f,
which computes in the working precision in every instantiation: in `float` all arithmetic is
done in `float`. `raw/environment.json` records its SHA-256 (beginning `0573ea1c`) and the
SHA-256 of every other source file, and `results/provenance.md` records the compilers and run
times. The header released as 3.0.0 differs from that file only in its comments, so its SHA-256
differs but it compiles to the same code.

- In double precision the stabilized calculation had no avoidable loss in any setting. Every
  non-finite output it returned had a component too large for `double`.
- In single precision it returned one draw as non-finite whose largest component could have
  been stored, at κ = 0.501, out of 5 × 10⁷ attempts; the same draw occurred in both builds.
  Its exact largest component lay 2.7 × 10⁻⁶ (relative) below the value at which a `float`
  overflows. The single-precision radius carries a relative error of about 10⁻⁵ at that size,
  so a component that close to the limit can round either way. The analysis counts such a draw
  as a rounding-band loss rather than an avoidable one (see below). Every other non-finite
  output had a component too large for `float`.
- The direct calculation lost additional draws near κ = 1/2. In double precision at κ = 0.505,
  for example, it failed on 2.4% of attempts, against 0.083% for the stabilized calculation. In
  single precision at κ = 0.55 the two rates were 0.57% and 0.013%.
- Unavoidable losses dominate at the smallest κ. At κ = 0.5001 they affect 86.8% of attempts in
  double precision and 98.2% in single precision.
- The 100-digit recomputation agreed with the probe's classification on all 115 021 896
  recorded attempts.
- The tests of the radial distribution (`p1`) and of the complete loader (`p4`) found no
  departure from the target distribution at the significance levels set in `PROTOCOL.md`.
- The two builds returned bitwise identical output for the stabilized calculation in all 130
  configurations compared. The direct calculation, which draws its Gamma variates from the
  standard library, differed between the builds in all 130.
- With libc++ the stabilized calculation took 0.71 to 0.92 times as long per returned sample as
  the direct one, depending on the benchmark case. With libstdc++ it took 1.15 to 1.40 times as
  long. The single-precision case is the fastest of the five.

Only one processor architecture (64-bit ARM) was tested.

**Rounding-band losses.** A calculation in the working type cannot tell on which side of the
overflow threshold a component falls when the exact component lies within its rounding error
of the threshold. A failure at the step that forms the output component, whose exact largest
component `V` lies below the threshold by a relative distance of at most `4 eps max(1, |log V|)`
(4.2 × 10⁻⁵ in `float`, 6.3 × 10⁻¹³ in `double`), is therefore counted separately, in the
`rounding_band_*` columns of `results/failure_envelope.csv`, and not as an avoidable loss. In
this run that applied to the one single-precision draw above. The direct calculation had 28
failures in single precision whose exact largest component also lay within the band (9 at
κ = 0.501, 8 at κ = 0.505, 11 at κ = 0.51). Their denominator `X₂` had underflowed to zero, so
the failure did not occur at the final step, and they are counted as avoidable losses. The
stabilized calculation never forms `X₂`, so every one of its failures is a failure of the final
step.

## The full study (repository only)

The files below are in the GitHub repository and not in the release archive.

**Rerunning.** Run from this directory. The builds need `clang++` with libc++, optionally a
Homebrew `g++-15`, `g++-14` or `g++-13` for the second build, the Boost headers (set `BOOST_INC`,
default `/opt/homebrew/include`) and MPFR (set `MPFR_PREFIX`, default `/opt/homebrew`). Python
dependencies come from `../../python` through `uv`. `make check-legacy`, which every phase runs
through `make preflight`, reads commit 53dbd30 with git, so it needs a clone with history.

```bash
make selftest                 # build; consistency checks on seeds 7501-7505, which write no data
make preflight                # record the environment in raw/environment.json
make p1 p2 p3 p4 p5 p6        # simulation phases, both builds
make oracle                   # recompute the recorded attempts at 100 digits
make analyze                  # write results/
make figures                  # write figures/
make checksums                # rewrite raw/raw_checksums.sha256 and checksums.sha256
make verify                   # check both checksum files
make reverify                 # regenerate results/ and figures/ and compare them with the committed files
```

Three further targets check the sources: `make cxx11-check` compiles the probe as strict
C++11 under both compilers, `make check-legacy` confirms that the vendored comparator
`src/legacy/bi_kappa_distribution_v1.H` matches `cpp/bi_kappa_distribution.H` at commit 53dbd30
apart from its include guard, namespace and comment banner, and `make protocol-check` confirms
that `config/protocol.json` is what `config/make_protocol.py` produces.

`make reverify` and `make protocol-check` compare every word, integer and row exactly and
every floating-point number to a relative tolerance of 1e-9, so that a different numpy or
scipy, which moves the last digits of some p-values, does not fail them. The protocol's one
Monte Carlo estimate, the F2 miss-count distribution, is checked by re-deriving the protocol
from its committed counts and by testing a fresh run's counts for consistency with them,
because numpy does not keep that random stream fixed across versions. `STRICT=1` requires
byte-for-byte equality instead; the committed files meet it under the package versions
recorded in `raw/environment.json`. `compare_regenerated.py` does both comparisons.

- Every simulation target runs `preflight` first. `preflight` stops if any source file the
  experiment depends on has uncommitted changes. `ALLOW_DIRTY_DEV=1 make preflight` overrides
  this and marks the run as exploratory.
- `preflight` rewrites `raw/environment.json`, and `make checksums` rewrites both checksum
  lists. These files are part of the committed record, so running either changes it. After an
  edit to a committed source, README or result that leaves the run data unchanged,
  `make bundle-checksums` rewrites `checksums.sha256` alone.
- Adding `SMOKE=1` to any target runs 2000 attempts on one seed and writes its data under
  `raw/smoke/` and `results/smoke/`. The `preflight` it runs first still rewrites
  `raw/environment.json`; `git checkout raw/environment.json` restores the committed record.
- `ORACLE_JOBS` sets the number of processes that share the 100-digit recomputation. On macOS
  the default is the number of performance cores.
- `make clean` removes the executables and caches. `make distclean` also removes the smoke
  runs, the PNG previews and the per-attempt binary files. Neither removes a committed file.

**Runtime and disk use.** On the test machine (Apple silicon, 12 performance cores) the six
phases took 34 minutes for both builds together, 21 of them in `p2`. The 100-digit
recomputation takes about 29 µs per recorded attempt, about 55 minutes of single-core time for
the 115 million recorded attempts, and took 4 minutes with `ORACLE_JOBS=16`. The analysis takes
about 7 minutes. The phases write about 11 GB of per-attempt binary files under `raw/`, most of
it the two `p2` audit streams. These `.bin` files are not committed.

**Checking the committed files.** `shasum -a 256 -c checksums.sha256` checks the committed
sources, settings, results and PDF and SVG figures; it lists committed files only. `make
verify` also checks `raw/raw_checksums.sha256`, which lists the uncommitted `.bin` files, so it
passes only after the phases have been rerun.

**Acceptance criteria.** `PROTOCOL.md` sets seven acceptance criteria, G0 to G6, before the
run. `results/analysis_report.md` records that G0 to G5 are met: the checks of the evaluation
code, the tests of the radial distribution and of the complete loader, the accounting of every
loss, the bitwise agreement of the two builds, and the cost bound. G6 requires, among other
provenance checks, a record that `make verify` and `make reverify` succeeded on the finished
results, supplied to the analysis as `g6_evidence.json`. No such record was supplied, so G6 is reported as failed and the report's
first line, its overall verdict, reads NO-GO. The failure of G6 concerns only that missing
record, not a measured result.

**Portability on other machines.** `results/portability_remote.md` gives the commands for
running the phase `p5` comparison on other machines, including x86_64, and
`.github/workflows/portability.yml` runs it on hosted runners when started by hand.
`config/protocol.json` limits the portability criterion to the tested architecture.

**Notes on the committed records.** These files record the run as it was made and are not
edited afterwards.

- `PROTOCOL.md` is titled "Experiment 7" and refers to an earlier study, "Experiment 6", whose
  files are not in this repository. It also uses the terms holdout, for a production run on
  a block of seeds not used before, and battery, for a fixed set of statistical tests.
  `config/protocol.json` lists, under `seeds.spent_blocks`, the seed blocks of earlier runs of
  this study, which are not included here; this run uses none of them.
- `PROTOCOL.md`, `config/protocol.json` and the comment banner of
  `src/legacy/bi_kappa_distribution_v1.H` describe the comparator as the header of release 1.0.0,
  and the banner names commit `0fc2c95`, which is not in this repository's history. That
  description is inaccurate. The
  vendored file reproduces the header at commit 53dbd30, which forms the radius as `√X₁/√X₂`;
  release 1.0.0 formed it as `√(X₁/X₂)`. `make check-legacy` checks against 53dbd30.
- `PROTOCOL.md` §2.9.2 states that the rounding band applies to every calculation alike. The
  analysis applies it only to failures at the step that forms the output component, as
  described under Results. The two rules differ only for the 28 single-precision underflow
  losses of the direct calculation, and the same draws of the quotient-first diagnostic, which
  the analysis counts as avoidable. The counts of the stabilized calculation are the same under
  both rules.
- Comments in the probe, the recomputation, the analysis and the protocol call the stabilized
  header "release 2.0.0", and some name versions 2.1.0 to 2.3.0. Those were development
  version numbers that were never released (see `CHANGELOG.md`); the committed run used
  version 3.0.0. The C++ sources and the protocol are left as they were run, since
  `raw/environment.json` records their hashes. Several Python scripts (`preflight.py`,
  `analyze.py`, `make_figures.py`, `exp7_portability.py`, `src/check_legacy.py`) have been
  corrected since the run without changing any committed number; their hashes in
  `raw/environment.json` are therefore those of the versions that ran.
- `raw/environment_open.json` is the environment record of an earlier run of this study
  (2026-09-20, under an earlier protocol and header, commit `8c36438`, which is not in this
  repository's history). No file of the committed run refers to it.
- `results/portability_remote.md` was written before the production run; a note at the top of
  its status table says what has run since.

## Output files

| path | contents |
|---|---|
| `raw/p1/` … `raw/p6/` | per-configuration counters as JSON lines (committed) and per-attempt binary files (not committed) |
| `raw/manifest.csv` | one row per output file, written by the process that produced it: phase, method, precision, κ, seed, size, SHA-256, build |
| `raw/environment.json` | commit, source hashes, compilers and floating-point environment of each build |
| `raw/oracle_audit.jsonl`, `raw/oracle_disagreements.jsonl`, `raw/oracle_rounding_band.jsonl` | 100-digit recomputation results, any disagreements, and the rounding-band records |
| `results/*.csv` | source data for every figure and number; `results/source_data_README.md` defines every column |
| `results/exp7_results.json` | the same results in one JSON file |
| `results/analysis_report.md`, `results/validation_matrix.md` | outcome of every statistical test and consistency check |
| `results/provenance.md` | run times and build identity |
| `results/schema.md` | binary record formats of the probe's output |
| `figures/*.pdf`, `figures/*.svg` | figures; `figures/captions.md` has their captions and `figures/figure_manifest.json` their source data |
| `checksums.sha256` | SHA-256 of the committed sources, settings, results and figures |
| `raw/raw_checksums.sha256` | SHA-256 of every file under `raw/p*/`, the manifest and the recomputation output |
