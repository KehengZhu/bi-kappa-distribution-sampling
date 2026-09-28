# CI helpers for the portability comparison of the finite-precision experiment

Four files, used by `.github/workflows/portability.yml` and by the hand-run recipe in
`experiments/exp4_finite_precision/results/portability_remote.md`. Between them they produce
the *inputs* to the experiment's portability criterion (`gates.G4` in
`experiments/exp4_finite_precision/config/protocol.json`) and check that the inputs are what
they claim to be. They decide nothing statistical: the tests of family F7 and the portability
decision belong to `experiments/exp4_finite_precision/analyze.py --portability-ingest`,
which reads `config/protocol.json`.

This directory and the workflow are in the repository only; the release archive leaves them
out. The workflow runs only when started by hand (`workflow_dispatch`).

| file | what it does |
|---|---|
| `env_probe.cpp` | a C++11 program that prints the arithmetic it was compiled into: rounding mode, subnormal flushing (measured, plus the MXCSR or FPCR word), `LDBL_MANT_DIG`, `FLT_EVAL_METHOD`, the exact overflow thresholds, the standard-library version macros, the engine's range, and the loader version it included |
| `make_env_record.py` | merges that with the toolchain's `--version` and `-dumpmachine`, the host's identity, and the compile line make actually used, into one `environment_<tag>.json`. Computes `execution` (native or translated) rather than taking it on trust, and fails if the frozen build flags are absent from the compile line |
| `check_evidence.py` | presence and shape of the uploaded evidence against `expected_environments.json`: every environment there, every run native, the frozen sizes (`matrix.n_scalar`) and seeds (`seeds.production`) of `config/protocol.json` used. Exits non-zero when a required environment is missing or incomplete, so that an incomplete run cannot end as a green check; a missing optional environment is reported only |
| `expected_environments.json` | the environments the workflow runs, in one place: for each, `env_id`, `arch`, `os`, `runner`, the standard-library builds it must provide (`tags`), whether it is `required`, and whether `config/protocol.json` scopes G4 to it (`in_protocol_scope`) |

## The environments

| `env_id` | runner | builds | required | in the protocol's G4 scope |
|---|---|---|---|---|
| `macos-arm64` | `macos-latest` | Apple clang with libc++, Homebrew GCC with libstdc++ | yes | yes |
| `linux-x86_64` | `ubuntu-latest` | clang with libc++, GCC with libstdc++ | yes | no |
| `linux-arm64` | `ubuntu-24.04-arm` | clang with libc++, GCC with libstdc++ | no | no |

`config/protocol.json` limits G4 to arm64 macOS with both standard libraries, and the
committed run met it on the development host (`experiments/exp4_finite_precision/raw/p5/`).
The macOS leg repeats that comparison on a hosted runner. The Linux legs supply the
cross-architecture comparison, which is outside the protocol's scope. `linux-x86_64` is
nevertheless required, because `analyze.py --portability-ingest` does not report the
criterion as met while only one architecture has native results. `linux-arm64` is optional
because the availability of its runner label depends on the account.

## The bundle layout

One directory per environment, named `portability-<env_id>`:

```
portability-linux-x86_64/
  environment_libcxx.json      environment_libstdcxx.json
  portability_libcxx.jsonl     portability_libstdcxx.jsonl
  run.log
  run_manifest.json
```

`portability_<tag>.jsonl` is the probe's `raw/p5/p5_<tag>.jsonl`, renamed. `run.log` is the
output of `make p5`, from which `make_env_record.py` takes the compile lines.
`analyze.py --portability-ingest` is given the directory that holds these, together with
the `coverage.json` that `check_evidence.py` writes beside them.

## What is not here

Per-attempt binaries. The portability decision uses counter rows and, for the
cross-standard-library comparison, the per-row digest `digest_sha256` of the returned
vectors; the bulk output stays on the runner.
