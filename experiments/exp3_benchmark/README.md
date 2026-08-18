# Experiment 3 — reproducible absolute and comparative performance benchmark

Answers **R1.4** (primary) and the "how does it outperform" part of **R2.A2**.

## The result, stated first because it is unfavourable to the manuscript

**The released implementation is not the fastest method tested. It is the slowest.**

Zenitani (2025)'s Pareto-envelope rejection sampler costs **0.38×–0.56×** what our
Gamma-ratio implementation costs, over κ ∈ [1.5, 50] — i.e. it is **1.8×–2.6× faster**,
despite being a rejection method with ~0.73–0.81 acceptance. Abdul & Mace (2015)'s
normal-triple scale mixture — which is the *same construction as ours*, differing only in
how the direction is bought — is also faster, by ≈1.4×–1.6×.

This is reported because it is what the measurement says. Every "fast", "resolves
computational bottlenecks", and "outperforms" claim in the manuscript must go, and R1.4
must be answered by **retiring** the performance claim, not by supporting it.

## Correctness before timing, enforced by the harness

No timing number here is meaningful unless the method being timed samples the intended
target law. `make run` therefore runs validation first, and `exp3_analyze.py` marks a
method's timings `usable: false` if it failed its own gate.

All three methods **pass**: radial law (KS **and** Cramér–von Mises against
`W = 1/(1+T) ~ Beta(κ−1/2, 3/2)`), directional uniformity, and zero non-finite draws,
at every κ tested, over 3 seeds × 2×10⁵ draws.

**A transcription check worth more than the gate.** The measured Pareto acceptance
reproduces Zenitani's published values to three digits:

| κ | measured here | Zenitani (2025) reports |
|---|---|---|
| 1.5 | 0.8060 | 0.806 |
| 2 | 0.7856 | 0.785 |
| 5 | 0.7505 | 0.750 |
| 50 | 0.7327 | → `√(πe)/4 ≈ 0.7314` asymptotically |

Agreement this close on a quantity we did not fit is strong evidence the algorithm was
transcribed faithfully rather than approximated — which is the condition the revision
plan sets before any published method may be benchmarked.

## The three methods, and which are actually distinct

| Tag | Source | Status |
|---|---|---|
| `gamma_ratio_spherical` | **The released header**, `cpp/bi_kappa_distribution.H` in `no_cap()` mode | baseline; equivalent to Zenitani & Nakano (2022) Alg. 1-1 and ZUM (2026) Alg. 3.1 |
| `scale_mixture_normals` | Abdul & Mace (2015), *Phys. Plasmas* **22**, 102107, Eq. (22) with Eqs. (19)–(20) | **implementation variant of the baseline, not a rival algorithm** |
| `pareto_rejection` | Zenitani (2025), *RNAAS* **9**, 299, §2 procedure, envelope index `n = κ/2` | **genuinely distinct algorithm**; uniform variates only |

`scale_mixture_normals` is labelled a variant because it provably is one. A&M 2015 Eq. (22)
with their Eqs. (19)–(20) reduces to `v_i = θ√κ · Z_i / √(χ²_ν)`, `ν = 2κ−1`; since
`|Z|² ~ χ²₃ = 2·Ga(3/2,1)` and `χ²_ν = 2·Ga(κ−1/2,1)`, the radius is *exactly* the
baseline's Gamma ratio. Inflating it into a competing algorithm is precisely what R2.A2
warns against, so the benchmark reports it as what it is: the same construction buying its
direction from three normals instead of two uniforms — and that choice is worth ≈1.5×.

> **Disclosure.** A&M 2015 never states how the non-integer-ν χ² deviate is generated.
> `χ²_ν = 2·Ga(ν/2,1)` via `std::gamma_distribution` is **our** choice, not theirs. Any cost
> difference attributable to a different χ² route is not attributable to Abdul & Mace.

### Where the rejection method does not apply

Zenitani's recommended envelope index `n = κ/2` requires `0 < n < κ − 1/2`, i.e. **κ > 1**.
At κ = 0.75 and κ = 1.0 it is **inapplicable**, and the harness records that as a result
rather than silently skipping it. His quoted efficiencies begin at κ = 1.5. So the honest
comparative statement is bounded: the rejection method wins *where it is defined*, and the
Gamma-ratio route covers a κ range it does not.

## Reproducing

```bash
make run                                              # validate, then time; writes checksums
uv run --project ../../python python exp3_analyze.py  # -> results/
make verify                                           # re-check a regenerated raw/
```

## Configuration

| | |
|---|---|
| Validation | κ ∈ {0.75, 1, 1.5, 2, 5, 10}, seeds 3001–3003, 2×10⁵ draws each |
| Timing | κ ∈ {0.75, 1, 1.5, 2, 5, 10, 50}, seed 3101, 10⁶ per batch, **10 batches** + 1 discarded warm-up |
| RNG | `std::mt19937` for **every** method, seeded identically — no method gets a cheaper generator |
| θ | isotropic (θ⊥ = θ∥ = 1) for the cross-method comparison, because M2 and M3 as published are isotropic |

Environment (compiler, target, flags, stdlib, arch, git revision, and the sampler header's
SHA-256) is recorded in `results/exp3_results.json`.

### Fairness measures, all deliberate

- The comparison is the **isotropic core only**. The released implementation's anisotropy
  and field-rotation paths are timed as *its own variants*, never charged against M2/M3.
- Timing loops accumulate a checksum, so the optimizer cannot delete the work being timed.
  The same checksum is applied to all three methods, so it cannot bias the comparison.
- One warm-up batch is discarded — first-call allocation is not steady-state per-sample cost.
- Every configuration reports the **distribution** over 10 batches (median, min–max, IQR),
  never a single timing.

## Secondary findings

**1. Per-sample cost is not constant in κ.** The baseline ranges 88–130 ns/sample
non-monotonically across the ladder. The fastest point is κ = 1.5 (88 ns), the slowest
κ = 2 (130 ns). Note `shape(x₂) = κ − 1/2` crosses 1 exactly at κ = 1.5, and library Gamma
generators switch algorithm at shape 1 — the plausible mechanism, though we did not
instrument the standard library to confirm it. Either way: **the Abstract's "constant time
per sample" is contradicted by direct measurement**, independently of the cause.

**2. Anisotropy and arbitrary-**B** rotation are essentially free.** `aniso` and `rotated`
sit within ≈1–2 ns of `iso` at every κ. C2 (arbitrary field-frame loading) costs nothing
measurable — a genuinely favourable result, and the only performance statement in this
experiment that flatters the implementation.

**3. The cap's cost tracks its rejection rate, as Experiment 2 predicts.** `capped20` is
indistinguishable from `iso` for κ ≥ 1.5 but costs +36% at κ = 0.75 (140.9 vs 103.9
ns/sample) — where Exp 2 measures 21.7% rejection at λ = 20. Cost and distortion move
together, exactly as the "rejected fraction *is* the TV distance" result implies.

## What the manuscript may and may not say

**May, with the parameter range attached:** absolute throughput ≈8–11 M samples/s for the
released implementation on the recorded hardware; anisotropy and frame rotation add no
measurable cost; the capped mode's overhead is confined to low κ.

**May not, and this is now settled by data rather than by caution:** "fast", "resolves
computational bottlenecks", "constant time per sample", "outperforms", or any implication
that rejection sampling is inefficient for this problem. The dedicated rejection method is
the fastest thing in the table.

## Artifacts

| Path | Canonical? |
|---|---|
| `exp3_bench.cpp`, `exp3_analyze.py`, `GNUmakefile`, `README.md` | **canonical** (committed source) |
| `results/exp3_results.json`, `results/exp3_table.md` | **canonical** (compact summaries, committed) |
| `raw/validate.jsonl`, `raw/timing.jsonl` | **canonical**, tracked — small, and the direct input to every reported number |
| `raw/checksums.sha256` | **canonical**, tracked |
| `raw/val_*.bin` | regenerable (`make run`) — bulk validation dumps, gitignored |
| `exp3_bench.exe` | build output, gitignored |

## Caveat on the timing numbers

These are single-machine, single-toolchain wall-clock measurements (Apple clang / libc++ /
arm64). They are reproducible on that machine and adequate to refute a constant-time claim
and to establish a ≈2× ordering, which is what R1.4 needs. They are **not** a
cross-platform performance characterization, and no claim of one is made. In particular the
Gamma-generator cost that dominates the baseline is a standard-library implementation
detail and may order differently elsewhere.
