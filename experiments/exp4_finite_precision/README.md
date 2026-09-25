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

- **Direct calculation**, which forms `X₂` and evaluates `√X₁/√X₂`. It is the calculation of the
  header before release 2.0.0, vendored as `src/legacy/bi_kappa_distribution_v1.H` with only its
  include guard renamed and a namespace added. The data files label it `LEGACY`.
- **Stabilized calculation**, which never forms `X₂`. It obtains `log X₂ = log G + (log ξ)/α`
  with `G ~ Γ(α + 1)` and `ξ ~ U(0, 1)`, evaluates `log R = (log X₁ − log X₂)/2`, and keeps
  `log R` until the final components are formed. This is `cpp/bi_kappa_distribution.H`. The data
  files label it `CANDIDATE`.

A third column, `QF` (`√(X₁/X₂)`, whose intermediate quotient can overflow when its square root
would not), appears only in the paired comparison described below. It is a diagnostic, not a
sampler under test.

Source and output files carry the prefix `exp7`.

## Paper

Fig. 2 of the paper is `figures/fp1_failure_envelope.pdf`. It is drawn by `make_figures.py` from
`results/failure_envelope.csv` and copied into the paper by
`paper/figures/make_manuscript_assets.py`. The four other figures in `figures/` are not used in
the paper.

The release archive (the GitHub and Zenodo downloads) contains only what reproduces Fig. 2 from
the committed results: this README, `GNUmakefile`, `make_figures.py`,
`results/failure_envelope.csv` with its column description `results/source_data_README.md`, and
the figure files. There, `make figures` redraws Fig. 2. The protocol, the analysis code, the
sampler probe, the 100-digit recomputation and the run records described below are in the
GitHub repository.

## Design

| | |
|---|---|
| κ | 0.5001, 0.501, 0.505, 0.51, 0.55, 0.6, 0.75, 1, 1.25, 1.49, 1.5, 2, 5 |
| precision | `float` and `double` |
| attempts | 10⁶ per seed and setting, seeds 10001–10005 (5 × 10⁶ per setting); `p2` 10⁷ per seed (5 × 10⁷ per setting); timing on seeds 10006–10010 |
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

The experiment runs in six phases, each a make target.

| target | measures |
|---|---|
| `p1` | the radial distribution of both calculations over the full κ list: Anderson–Darling, KS and Cramér–von Mises tests of `Z = −log I_W(α, 3/2) ~ Exp(1)`, quantiles of `log R`, and exceedance counts in the upper tail |
| `p2` | returned, avoidable and unavoidable counts, paired and native; uncapped, isotropic, field along `ẑ` |
| `p3` | the probability that an attempt succeeds as a function of the radius it was meant to produce, and the resulting loss of tail mass |
| `p4` | the complete loader with anisotropy, rotation and optional cap, in seven configurations C0–C6 |
| `p5` | a SHA-256 digest of every returned value and counter, compared between the two builds |
| `p6` | time per returned sample, stabilized against direct, in five benchmark cases |

`PROTOCOL.md` specifies the design, every statistical test and every threshold in full.
`config/protocol.json` holds the same settings in machine-readable form; it is generated by
`config/make_protocol.py`, and the C++ probe reads it through the generated header
`src/exp7_protocol.H`.

## Results

The committed results were produced with release 2.3.0 of `cpp/bi_kappa_distribution.H`, which
computes in the working precision in every instantiation: in `float` all arithmetic is done in
`float`. `raw/environment.json` records its SHA-256 (beginning `4c7c95ea`) and the SHA-256 of
every other source file, and `results/provenance.md` records the compilers and run times.

- In double precision the stabilized calculation had no avoidable loss in any setting. Every
  non-finite output it returned had a component too large for `double`.
- In single precision it returned three draws as non-finite whose largest component could have
  been stored: two at κ = 0.505 and one at κ = 0.51, out of 5 × 10⁷ attempts per setting, the
  same three draws in both builds. Their exact largest components lay 4.3 × 10⁻⁶, 9.8 × 10⁻⁶
  and 5.2 × 10⁻⁷ (relative) below the value at which a `float` overflows. The single-precision
  radius carries a relative error of about 10⁻⁵ at that size, so such a component can round
  either way; one draw at κ = 0.51 rounded the other way and was returned finite although its
  exact value lay 6.2 × 10⁻⁶ above the limit. These three draws are avoidable losses under the
  protocol. Every other non-finite output had a component too large for `float`.
- The direct calculation lost additional draws near κ = 1/2. In double precision at κ = 0.505,
  for example, it failed on 2.4% of attempts, against 0.083% for the stabilized calculation. In
  single precision at κ = 0.55 the two rates were 0.57% and 0.014%.
- Unavoidable losses dominate at the smallest κ. At κ = 0.5001 they affect 86.8% of attempts in
  double precision and 98.2% in single precision.
- The 100-digit recomputation agreed with the probe's classification on all 115 029 238
  recorded attempts.
- The tests of the radial distribution (`p1`) and of the complete loader (`p4`) found no
  departure from the target distribution at the significance levels set in `PROTOCOL.md`.
- The pre-registered verdict in `results/analysis_report.md` is NO-GO, because gate G1 requires
  the stabilized calculation to have no avoidable loss at all, and the three single-precision
  draws above fail it. The gate counts them once per build, as six.
- The two builds returned bitwise identical output for the stabilized calculation in all 130
  configurations compared. The direct calculation, which draws its Gamma variates from the
  standard library, differed between the builds in all 130.
- With libc++ the stabilized calculation took 0.70 to 0.91 times as long per returned sample as
  the direct one, depending on the benchmark case. With libstdc++ it took 1.15 to 1.40 times as
  long. The single-precision case is the fastest of the five.

Only one processor architecture (64-bit ARM) was tested. `results/portability_remote.md` gives
the commands for running the comparison on x86_64.

**2026-09-25, protocol 4.0.0.** The results above replace those of protocol 3.0.0, which were
produced with release 2.2.0 on seeds 9001–9010. Release 2.2.0 computed a `float` sample in
`double` and rounded each component once; release 2.3.0 computes it in `float`, so that the
header works where only single precision is available. `double` output is unchanged. The
mechanism phase `p2` now runs 10⁷ attempts per seed instead of 10⁶. With 5 × 10⁶ attempts per
setting a setting with no failure could only be bounded at 6.0 × 10⁻⁷ (one-sided 95%), above the
nonzero fractions observed at neighbouring κ; with 5 × 10⁷ the bound is 6.0 × 10⁻⁸. The
recomputation now uses MPFR instead of `cpp_dec_float_100`. `PROTOCOL.md` §2.8 records the
amendment.

## Rerunning

Run from this directory. The builds need `clang++` with libc++, optionally a Homebrew `g++-15`,
`g++-14` or `g++-13` for the second build, and the Boost headers (set `BOOST_INC`, default
`/opt/homebrew/include`). Python dependencies come from `../../python` through `uv`.

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
C++11 under both compilers, `make check-legacy` confirms that the vendored comparator matches
`cpp/bi_kappa_distribution.H` at commit `0fc2c95` apart from its include guard and namespace, and
`make protocol-check` confirms that `config/protocol.json` is what `config/make_protocol.py`
produces.

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
- `preflight` rewrites `raw/environment.json`, and `make checksums` rewrites `checksums.sha256`.
  Both files are part of the committed record, so running either changes it.
- Adding `SMOKE=1` to any target runs 2000 attempts on one seed and writes everything under
  `raw/smoke/` and `results/smoke/`, leaving the committed data untouched.
- `ORACLE_JOBS` sets the number of processes that share the 100-digit recomputation. On macOS
  the default is the number of performance cores.

**Runtime and disk use.** On the test machine (Apple silicon, 12 performance cores) the six
phases took 34 minutes for both builds together, 21 of them in `p2`. The 100-digit
recomputation takes about 29 µs per recorded attempt, about 55 minutes of single-core time for
the 115 million recorded attempts, and took 4 minutes with `ORACLE_JOBS=16`. The analysis takes
about 7 minutes. The phases write about 11 GB of per-attempt binary files under `raw/`, most of
it the two `p2` audit streams. These `.bin` files are not committed.

**Checking the committed files.** `shasum -a 256 -c checksums.sha256` checks the committed
sources, settings, results and figures. `make verify` also checks `raw/raw_checksums.sha256`,
which lists the uncommitted `.bin` files, so it passes only after the phases have been rerun.

## Output files

| path | contents |
|---|---|
| `raw/p1/` … `raw/p6/` | per-configuration counters as JSON lines (committed) and per-attempt binary files (not committed) |
| `raw/manifest.csv` | one row per output file, written by the process that produced it: phase, method, precision, κ, seed, size, SHA-256, build |
| `raw/environment.json` | commit, source hashes, compilers and floating-point environment of each build |
| `raw/oracle_audit.jsonl`, `raw/oracle_disagreements.jsonl` | 100-digit recomputation results and any disagreements |
| `results/*.csv` | source data for every figure and number; `results/source_data_README.md` defines every column |
| `results/exp7_results.json` | the same results in one JSON file |
| `results/analysis_report.md`, `results/validation_matrix.md` | outcome of every statistical test and consistency check |
| `results/provenance.md` | run times and build identity |
| `results/schema.md` | binary record formats of the probe's output |
| `figures/*.pdf`, `figures/*.svg` | figures; `figures/captions.md` has their captions and `figures/figure_manifest.json` their source data |
| `checksums.sha256` | SHA-256 of the sources, settings, results and figures |
| `raw/raw_checksums.sha256` | SHA-256 of every file under `raw/p*/`, the manifest and the recomputation output |
