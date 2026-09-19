# E0 - clean baseline reproduction

Generated 2026-09-19T17:43:33+00:00.

Experiment 4's saved result records a dirty working tree, and the clean git object
at the commit it names does not reproduce the sampler-header hash it recorded, so its
numbers are a historical baseline rather than reproducible evidence.  This phase reruns
the same protocol - same seeds 4001-4005, same 10^6 attempts per seed, same isotropic
unrotated uncapped configuration, same released class - from a recorded source state,
and states plainly whether the numbers come back.

**Verdict: reproduced, residual difference attributed to floating-point contraction.**

Cells compared: 36.  Bitwise-identical counts: 32.
Cells that are bitwise identical once the historical build's floating-point contraction setting is restored: 4.  Cells still unexplained: 0.

Experiment 4 was built without an explicit `-ffp-contract` flag and so inherited the
compiler default, which enables fused multiply-add; this experiment disables it, because
contraction changes the rounding of exactly the products whose overflow is under study.
That setting was not recorded in the Experiment 4 manifest. `make baseline-contract-check`
rebuilds the libstdc++ probe with the historical setting and reruns the single-precision
baseline; the rows it produces carry the tag `libstdcxx-contractfast` and are a provenance
control, not production data.
Cells that change a qualitative operating-range statement: 0.

Same seed and same class means the two runs should agree bit for bit unless the
released header changed between them.  A standardized difference is reported as well,
because agreement across standard libraries is statistical, not bitwise.

| precision | stdlib | kappa | Exp4 rate | Exp6 rate | difference | z | identical |
|---|---|---:|---:|---:|---:|---:|---|
| double | libc++ | 0.5001 | 0.928304 | 0.928304 | +0 | +0.00 | yes |
| double | libc++ | 0.501 | 0.475106 | 0.475106 | +0 | +0.00 | yes |
| double | libc++ | 0.505 | 0.0242196 | 0.0242196 | +0 | +0.00 | yes |
| double | libc++ | 0.51 | 0.0005864 | 0.0005864 | +0 | +0.00 | yes |
| double | libc++ | 0.55 | 0 | 0 | +0 | +0.00 | yes |
| double | libc++ | 0.6 | 0 | 0 | +0 | +0.00 | yes |
| double | libc++ | 0.75 | 0 | 0 | +0 | +0.00 | yes |
| double | libc++ | 1 | 0 | 0 | +0 | +0.00 | yes |
| double | libc++ | 1.5 | 0 | 0 | +0 | +0.00 | yes |
| double | libstdc++ | 0.5001 | 0.928183 | 0.928183 | +0 | +0.00 | yes |
| double | libstdc++ | 0.501 | 0.474933 | 0.474933 | +0 | +0.00 | yes |
| double | libstdc++ | 0.505 | 0.0242746 | 0.0242746 | +0 | +0.00 | yes |
| double | libstdc++ | 0.51 | 0.0005622 | 0.0005622 | +0 | +0.00 | yes |
| double | libstdc++ | 0.55 | 0 | 0 | +0 | +0.00 | yes |
| double | libstdc++ | 0.6 | 0 | 0 | +0 | +0.00 | yes |
| double | libstdc++ | 0.75 | 0 | 0 | +0 | +0.00 | yes |
| double | libstdc++ | 1 | 0 | 0 | +0 | +0.00 | yes |
| double | libstdc++ | 1.5 | 0 | 0 | +0 | +0.00 | yes |
| float | libc++ | 0.5001 | 0.989657 | 0.989657 | +0 | +0.00 | yes |
| float | libc++ | 0.501 | 0.90177 | 0.90177 | +0 | +0.00 | yes |
| float | libc++ | 0.505 | 0.59637 | 0.59637 | +0 | +0.00 | yes |
| float | libc++ | 0.51 | 0.355371 | 0.355371 | +0 | +0.00 | yes |
| float | libc++ | 0.55 | 0.0056358 | 0.0056358 | +0 | +0.00 | yes |
| float | libc++ | 0.6 | 3.26e-05 | 3.26e-05 | +0 | +0.00 | yes |
| float | libc++ | 0.75 | 0 | 0 | +0 | +0.00 | yes |
| float | libc++ | 1 | 0 | 0 | +0 | +0.00 | yes |
| float | libc++ | 1.5 | 0 | 0 | +0 | +0.00 | yes |
| float | libstdc++ | 0.5001 | 0.989671 | 0.989671 | +0 | +0.00 | yes |
| float | libstdc++ | 0.501 | 0.901802 | 0.901802 | +2e-07 | +0.00 | contraction |
| float | libstdc++ | 0.505 | 0.596809 | 0.596813 | +3.2e-06 | +0.01 | contraction |
| float | libstdc++ | 0.51 | 0.35576 | 0.355762 | +2.2e-06 | +0.01 | contraction |
| float | libstdc++ | 0.55 | 0.0057344 | 0.0057346 | +2e-07 | +0.00 | contraction |
| float | libstdc++ | 0.6 | 3.3e-05 | 3.3e-05 | +0 | +0.00 | yes |
| float | libstdc++ | 0.75 | 0 | 0 | +0 | +0.00 | yes |
| float | libstdc++ | 1 | 0 | 0 | +0 | +0.00 | yes |
| float | libstdc++ | 1.5 | 0 | 0 | +0 | +0.00 | yes |
