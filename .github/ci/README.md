# CI helpers for Experiment 7's portability gate

Four files, used by `.github/workflows/portability.yml` and by the hand-run recipe in
`experiments/exp4_finite_precision/results/portability_remote.md`. Between them they produce
the *inputs* to gate G4 and check that the inputs are what they claim to be. They decide
nothing statistical: F7 and the G4 verdict belong to
`experiments/exp4_finite_precision/analyze.py --portability-ingest`, which reads
`config/protocol.json`.

| file | what it does |
|---|---|
| `env_probe.cpp` | a C++11 program that prints the arithmetic it was compiled into: rounding mode, subnormal flushing (measured, plus the MXCSR or FPCR word), `LDBL_MANT_DIG`, `FLT_EVAL_METHOD`, the exact overflow thresholds, the standard-library version macros, the engine's range, and the loader version it included |
| `make_env_record.py` | merges that with the toolchain's `--version` and `-dumpmachine`, the host's identity, and the compile line make actually used, into one `environment_<tag>.json`. Computes `execution` (native or translated) rather than taking it on trust, and fails if the frozen build flags are absent from the compile line |
| `check_evidence.py` | presence and shape of the uploaded evidence against `expected_environments.json`: every required environment there, every run native, the frozen sizes and seeds used. Exits non-zero when an environment is missing, so that an incomplete run cannot end as a green check |
| `expected_environments.json` | the G4 matrix, in one place, read by both the workflow and the hand-run recipe |

## The bundle layout

One directory per environment, named `portability-<env_id>`:

```
portability-linux-x86_64/
  environment_libcxx.json      environment_libstdcxx.json
  portability_libcxx.jsonl     portability_libstdcxx.jsonl
  run.log
  run_manifest.json
```

`analyze.py --portability-ingest` is given the directory that holds these, together with
the `coverage.json` written beside them.

## What is not here

Per-draw binaries. G4 is decided on counter rows and, for the cross-standard-library
claim, a per-row digest of the returned vectors; the bulk output stays on the runner.
