# Experiment 3 -- performance benchmark

Answers **R1.4** and the performance part of **R2.A2**.


## 1. Correctness gate (phase 1)

No timing below is believed for a method that fails here. alpha = 0.01, fixed before the runs, with Holm-Bonferroni across the family of simultaneous tests.


> **On the multiplicity correction.** The gate runs one radial and one directional test per (method, kappa, seed) — 48 tests. At alpha = 0.01 the expected number of false rejections is ~0.5, so a rule of "every single test must clear alpha" fails ~38% of the time on a *correct* sampler. Exactly one uncorrected rejection occurred: `gamma_ratio_spherical` directional uniformity at kappa = 5, seed 3003, p = 0.0044 — while the other two seeds at that kappa give p = 0.59 and p = 0.85, and Experiment 1 validated directional uniformity for this same sampler over 1.35e7 draws. The per-test alpha is unchanged; only the multiplicity is accounted for. Uncorrected verdicts are retained in the JSON under `*_uncorrected`.


| method | kappa tested | radial law | direction | non-finite | verdict |
|---|---|---|---|---|---|
| `gamma_ratio_spherical` | 0.75, 1, 1.5, 2, 5, 10 | PASS | PASS | 0 | **USABLE** |
| `scale_mixture_normals` | 0.75, 1, 1.5, 2, 5, 10 | PASS | PASS | 0 | **USABLE** |
| `pareto_rejection` | 1.5, 2, 5, 10 | PASS | PASS | 0 | **USABLE** |

### What each method actually is

- **`gamma_ratio_spherical`** — released implementation, cpp/bi_kappa_distribution.H (no_cap()); equivalent to Zenitani & Nakano 2022 Alg. 1-1 and ZUM 2026 Alg. 3.1
- **`scale_mixture_normals`** — Abdul & Mace 2015, Phys. Plasmas 22, 102107, Eq. (22) with Eqs. (19)-(20); IMPLEMENTATION VARIANT of the same construction, not a distinct algorithm; chisq_nu = 2*Ga(nu/2,1) is our choice -- the paper does not specify it
- **`pareto_rejection`** — Zenitani 2025, Res. Notes AAS 9, 299, Section 2 procedure, envelope index n = kappa/2 (the author's recommendation); genuinely distinct algorithm, uniform variates only

## 2. Timing, isotropic core (phase 2)

`std::mt19937` for every method; median of independent batches, with the full spread. One warm-up batch discarded.


| kappa | method | ns/sample (median) | min–max | IQR | Msamples/s | acceptance |
|---|---|---|---|---|---|---|
| 0.75 | `gamma_ratio_spherical` | 103.9 | 102.8–106.7 | 0.5 | 9.62 | 1.0000 |
| 0.75 | `scale_mixture_normals` | 62.9 | 62.4–63.6 | 0.3 | 15.91 | 1.0000 |
| 0.75 | `pareto_rejection` | — | — | — | — | **inapplicable**: n = kappa/2 violates 0 < n < kappa - 1/2 |
| 1 | `gamma_ratio_spherical` | 110.7 | 108.2–111.8 | 0.7 | 9.03 | 1.0000 |
| 1 | `scale_mixture_normals` | 67.6 | 66.6–68.3 | 0.5 | 14.79 | 1.0000 |
| 1 | `pareto_rejection` | — | — | — | — | **inapplicable**: n = kappa/2 violates 0 < n < kappa - 1/2 |
| 1.5 | `gamma_ratio_spherical` | 88.1 | 87.5–89.4 | 0.8 | 11.35 | 1.0000 |
| 1.5 | `scale_mixture_normals` | 46.6 | 45.5–47.3 | 1.0 | 21.46 | 1.0000 |
| 1.5 | `pareto_rejection` | 48.9 | 48.0–51.1 | 0.6 | 20.44 | 0.8060 |
| 2 | `gamma_ratio_spherical` | 129.9 | 129.3–131.5 | 0.9 | 7.70 | 1.0000 |
| 2 | `scale_mixture_normals` | 88.3 | 87.4–90.4 | 0.2 | 11.32 | 1.0000 |
| 2 | `pareto_rejection` | 49.9 | 49.6–51.5 | 0.4 | 20.04 | 0.7856 |
| 5 | `gamma_ratio_spherical` | 121.0 | 120.7–121.7 | 0.5 | 8.26 | 1.0000 |
| 5 | `scale_mixture_normals` | 80.5 | 79.8–86.0 | 0.6 | 12.41 | 1.0000 |
| 5 | `pareto_rejection` | 53.1 | 52.6–54.3 | 0.7 | 18.82 | 0.7505 |
| 10 | `gamma_ratio_spherical` | 124.2 | 123.2–125.8 | 0.9 | 8.05 | 1.0000 |
| 10 | `scale_mixture_normals` | 82.3 | 81.1–83.5 | 1.0 | 12.15 | 1.0000 |
| 10 | `pareto_rejection` | 53.7 | 53.2–54.1 | 0.3 | 18.62 | 0.7403 |
| 50 | `gamma_ratio_spherical` | 117.4 | 116.9–118.6 | 1.2 | 8.52 | 1.0000 |
| 50 | `scale_mixture_normals` | 74.8 | 74.5–75.8 | 0.6 | 13.37 | 1.0000 |
| 50 | `pareto_rejection` | 54.1 | 53.9–54.7 | 0.3 | 18.49 | 0.7327 |

## 3. Cost of the released implementation's own features

Attributable rather than folded into one number. `iso` is the baseline; `aniso` adds theta_par != theta_perp; `rotated` adds the arbitrary-**B** frame rotation; `capped20` adds the component-wise cap at lambda = 20.


| kappa | iso | aniso | rotated | capped20 |
|---|---|---|---|---|
| 0.75 | 103.9 | 104.0 | 105.3 | 140.9 |
| 1 | 110.7 | 108.2 | 110.5 | 116.2 |
| 1.5 | 88.1 | 88.9 | 88.9 | 88.7 |
| 2 | 129.9 | 129.4 | 130.5 | 129.2 |
| 5 | 121.0 | 121.5 | 122.7 | 122.2 |
| 10 | 124.2 | 124.6 | 126.0 | 124.6 |
| 50 | 117.4 | 117.4 | 118.3 | 120.1 |

## 4. What this licenses the manuscript to say

- Zenitani (2025) Pareto rejection costs **0.38x to 0.56x** the released Gamma-ratio implementation over kappa in [1.5, 50] (ratio < 1 means the rejection method is FASTER).
- Measured acceptance for the Pareto envelope is reported per kappa above; compare against the 0.73–0.8 the author reports **for kappa >= 3/2 only**.
- The released implementation's per-sample cost is **not** constant in kappa; read column `iso` in section 3 before writing any constant-time claim.

**Wording that remains forbidden regardless of these numbers:** "fast", "resolves computational bottlenecks", "prohibitively low acceptance", "constant time per sample", "outperforms" — unless the specific sentence is tied to the specific measurement and parameter range above.
