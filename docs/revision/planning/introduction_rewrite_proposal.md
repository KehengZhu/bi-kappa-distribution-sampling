# Introduction + contribution statement — rewrite proposal

**Scope:** answers R2.A1, R2.A2, R2.A3 (partly), R1.4, R1.7 and the Zenitani novelty matrix
(§4 of `reviewer_response_matrix.md`). This is a drafting proposal, not an applied edit.
`main.tex` is untouched.

**Source of truth for the prior art:**
`paper/reference/ZenitaniUsamiMatsukiyo2026_JGR131_e2025JA034669.pdf`, read in full.
All statements attributed to Zenitani et al. below are verified against that PDF.
See `paper/reference/README.md` for the full primary-source inventory.

---

## 1. What Zenitani et al. (2026) actually says — the facts that constrain us

These are the load-bearing facts. Every one of them was checked in the PDF, not recalled.

1. **Their Introduction already surveys prior Kappa sampling.** Verbatim: *"despite its
   importance, numerical procedures for Kappa distributions have been largely unknown to the
   community, until Abdul and Mace (2014, 2015) used its relation to Student's t distributions.
   More recently, Zenitani and Nakano (2022, 2023) proposed numerical procedures for Kappa,
   relativistic, and loss-cone distributions."*
   → Our current Introduction's premise ("inversion intractable + naive rejection inefficient,
   therefore a gap") is contradicted **in print, by name**, in a paper the referees will now see.
   This is the single most urgent deletion.

2. **They also state the limitations of the two generic routes** (acceptance–rejection needs a good
   envelope or "becomes increasingly inefficient, especially in three dimensions"; inverse transform
   "relies on the inverse of the CDF, which is not always available analytically" and is "difficult
   to extend to multiple dimensions"). We may still make this point — but it is now *their* framing,
   and must be cited to them rather than presented as our own diagnosis.

3. **Their §3 Algorithm 3.1 ("beta-prime method") is our core construction.** For the (r,q)
   distribution they draw `X1 ~ Ga(3/(2(1+r)), 1)`, `X2 ~ Ga(q − 3/(2(1+r)), 1)`, form
   `x ← [(q−1)X1/X2]^{1/(2(1+r))}`, then scatter over the sphere with two uniforms and scale
   anisotropically:
   `v∥ ← θ∥ x (2U1 − 1)`, `v⊥1 ← 2θ⊥ x √(U1(1−U1)) cos(2πU2)`, `v⊥2 ← … sin(2πU2)`.

4. **Their Eq. (12) is the reduction we must concede.** `f_κ(v; κ, θ∥, θ⊥) = f_rq(v; 0, κ+1, θ∥, θ⊥)`
   — i.e. `r = 0`, `q = κ + 1`. At those values Algorithm 3.1 becomes `X1 ~ Ga(3/2,1)`,
   `X2 ~ Ga(κ−1/2,1)`, `x = √(κ X1/X2)`, uniform direction, `(θ⊥, θ⊥, θ∥)` scaling. That is
   our Eqs. (23)–(29), term for term, **including the anisotropic case** — the θ∥ ≠ θ⊥ extension is
   explicitly in their algorithm box, not left as an exercise. There is no daylight here and we
   should not go looking for any.

5. **They flag the low-κ failure and then decline to solve it.** Verbatim caution: when
   `q ≤ 1 + 3/(2(1+r))` — on the Kappa line, exactly `κ ≤ 3/2` — the denominator gamma has shape
   `< 1`, has non-zero density at 0, and *"we need to take care of the division by zero."* Their
   recommendation is to **switch methods**: *"Below the blue curve … we recommend the piecewise
   rejection method, because the beta-prime method may encounter the zero-division problem."*
   → This is an opening, but a much narrower one than it looks. The measured behaviour is in §3a;
   it does not support promoting low-κ robustness to the paper's headline claim.

6. **Their capped analogue is the regularized Kappa (Scherer et al. 2017), not a box cutoff.**
   Their §4 treats a *physically motivated* Gaussian-cutoff Kappa with two rejection samplers. Our
   Sec. III.E component-wise box is a different and cruder object. The Introduction should not let a
   reader conflate them; the Background section should name the regularized Kappa as the principled
   alternative we are *not* implementing.

7. **Their release is a Jupyter notebook archived on Zenodo.** Data Availability Statement, verified.
   Not a C++ library, not a PIC-integrable API, no field-frame handling. This is narrow but real,
   and it is where the software contribution lives — provided we do not overstate it.

**Net effect on positioning:** the Gamma-ratio/Beta-prime construction is prior art as of March 2026,
for the anisotropic case, from a named group, in a journal our referees read. Independent development
is a fact worth *stating* and worth nothing as a *claim*. The paper survives on software +
validation. The low-κ regime is a *candidate* addition on top of that, not part of the floor —
see §3, which measured it and found the first version of this argument wrong.

---

## 2. Literature map for the new survey paragraph

Found via Consensus; all links are the canonical paper pages. The four marked **cite** are the
minimum needed for R2.A1 to be answerable; the rest strengthen R1.7.

| Role in the argument | Paper | In `refs.bib`? |
|---|---|---|
| **cite** — first practical Kappa loader; Student-t route | [A method to generate kappa distributed random deviates for particle-in-cell simulations](https://consensus.app/papers/details/6e8e4a5a972e5eaf9d5745e684558263/) (Abdul & Mace, 2014, *Comput. Phys. Commun.*) | ✅ `abdul2014kappa` |
| **cite** — Beta-prime/Gamma route, relativistic Kappa | [Loading a relativistic Kappa distribution in particle simulations](https://consensus.app/papers/details/9a7990a4beb35f2788c5a55863c0d2ac/) (Zenitani & Nakano, 2022, *Phys. Plasmas*) | ✅ `zenitani2022relativistic` |
| **cite** — the conflicting paper | [Loading Non-Maxwellian Velocity Distributions in Particle Simulations](https://consensus.app/papers/details/437b2946348b50f0a7c8a7ed0a39ce35/) (Zenitani, Usami & Matsukiyo, 2026, *JGR Space Phys.*) | ✅ `zenitani2026nonmaxwellian` |
| **cite** — inversion *is* viable numerically; kills our "intractable" line | [Fast Inverse Transform Sampling of Non-Gaussian Distribution Functions in Space Plasmas](https://consensus.app/papers/details/862b53e9cef85b7f80cee540e24f46c3/) (An et al., 2022, *JGR Space Phys.*) | ❌ **add** |
| Pareto-envelope rejection, ~0.73–0.8 acceptance — the quantitative counter to "prohibitively low" | [A Simple Procedure for Generating a Kappa Distribution in PIC Simulation](https://consensus.app/papers/details/31bdb0eafc315fe2bfaf5a23de76f8b1/) (Zenitani, 2025, *RNAAS*) | ✅ `zenitani2025simple` |
| Fourth 2026 Kappa sampler — q-exponential inverse transform, GPU-oriented | [An approximate Kappa generator for particle simulations](https://consensus.app/papers/details/32fa4c98fe7a53c39d3251cb7aa1a7ae/) (Zenitani et al., 2026, *Earth Planets Space*) | ❌ **add** |
| Loss-cone / KLC transformation algorithms | [Loading Loss-Cone Distributions in Particle Simulations](https://consensus.app/papers/details/f368faa1537e59118ea32b23b82da87a/) (Zenitani & Nakano, 2023, *JGR Space Phys.*) | ✅ `zenitani2023losscone` |
| Downstream use of a bi-Kappa loader in PIC — shows the demand is real | [2D PIC simulations of the whistler instability with bi-kappa electrons](https://consensus.app/papers/details/9410a1ad6cfa5e4e864a9df1983262b4/) (Abdul et al., 2021, *Phys. Plasmas*) | ✅ `abdul2021whistler` |
| Numerical basis for An et al.'s Chebsampling | [Fast inverse transform sampling in one and two dimensions](https://consensus.app/papers/details/27814d0fabd751c5acef9672ba43b9d8/) (Olver & Townsend, 2013) | ❌ optional |
| **R1.7** — foundational turbulence-driven ⊥ ion heating | [Perpendicular ion heating by low-frequency Alfvén-wave turbulence in the solar wind](https://consensus.app/papers/details/29ea76f179b25935ae15294a2cfdd5ba/) (Chandran et al., 2010, *ApJ* 720, 503) | ❌ **add** |

> ⚠ **Two rows above are now known to be wrong or misleading.** Corrected in §2a; the table is
> left as drafted so the reversal is visible rather than silent.
>
> - Row `abdul2021whistler` — "shows the demand is real" **understates it, in a way that costs
>   us a claim.** It is a 2D3V GPU PIC run with 6.7×10⁷ bi-Kappa particles per species. It does
>   not show demand; it shows the problem was already **solved and deployed at scale**.
> - Row `zenitani et al., 2026, Earth Planets Space` — the authorship is wrong. It is
>   **Zenitani & Umeda**, two authors.

### 2a. Corrections and additions from the primary-source rounds

Everything below was read in the primary source. See
`../literature/step1_claim_audit.md` §7 and `paper/reference/README.md`.

| Role | Source | Status |
|---|---|---|
| **The multivariate Kappa loader.** Isotropic **3-D**, genuinely trivariate: a normal vector over a **single shared** χ²_ν, Eq. (22), ν = 2κ−1. "One-dimensional" in the title is 1D3V *simulation geometry*. | Abdul & Mace (2015), *Phys. Plasmas* **22**(10), 102107, doi:10.1063/1.4933005 | ❌ **add** — verified |
| **Downstream multidimensional deployment.** 2D3V, bi-Kappa (Summers & Thorne form), T⊥/T∥ = 3, 6.7×10⁷ particles/species, using the 2015 loader. | Abdul, Matthews & Mace (2021), *Phys. Plasmas* **28**(6), 062104, doi:10.1063/5.0047638 | ✅ present as `abdul2021whistler` — **verify the fields** |
| **R1.7 instability bound.** Supplies the four-threshold formula + coefficients everyone reuses. | Hellinger, Trávníček, Kasper & Lazarus (2006), *GRL* **33**, L09101, doi:10.1029/2006GL025925 | ❌ **add** — verified |
| **R1.7 radial evolution** — the source for the word "expanding". | Matteini et al. (2007), *GRL* **34**, L20105, doi:10.1029/2007GL030920 | ❌ **add** — verified |
| **R1.7 anisotropy ⊗ Kappa in one observational fit** — bi-Maxwellian core + **anisotropic bi-Kappa halo**, ~124 000 events, 0.3–3.95 AU. ⚠ **ELECTRONS.** | Pierrard, Lazar, Poedts, Štverák, Maksimovic & Trávníček (2016), *Solar Phys.* **291**(7), 2165–2179, doi:10.1007/s11207-016-0961-7 | ❌ **add** — verified (⚠ on-disk copy is the arXiv preprint) |
| Optional: instabilities *actively* regulating, not merely coinciding. | Bale et al. (2009), *PRL* **103**, 211101 | ❌ optional (⚠ arXiv preprint on disk) |
| **R1.7 foundational ion observation.** Abstract-level only. | Marsch et al. (1982), *JGR* **87**(A1), 52–72, doi:10.1029/JA087iA01p00052 | ❌ **add**, but ⚠ **full text NOT obtained** — quote no numeric detail; author order is Marsch, **Mühlhäuser, Schwenn, Rosenbauer**, Pilipp, Neubauer |

**Two propositions remain BLOCKED and must not be written:**

1. **Any ion observation fitted with an anisotropic bi-Kappa.** None found. Pierrard 2016 is
   electrons; `yoon2023bikappa` already in our bib is *theory*. Do not generalize.
2. **Any page/figure-level detail from Marsch et al. (1982).** Abstract only, and the
   publisher-deposited abstract carries typos.

**And one framing is now dead:** any sentence implying prior Kappa loaders were not demonstrated
in multidimensional PIC. See §2b.

### 2b. What the Abdul chain leaves us — the surviving differentiator

The anisotropic bi-Kappa scale mixture is **used but written down by nobody**:

- Abdul & Mace (2015) §III B carries a general covariance **R** in Eq. (17), then immediately
  sets **R** = σ²**I**. Isotropic only.
- Abdul (2018) PhD thesis §2.3 — the fullest written statement of the loader anywhere — is
  isotropic only.
- Abdul, Matthews & Mace (2021) runs a genuine bi-Kappa but cites 2015 without the extension.
- ZUM2026 Algorithm 3.1: *"We can easily extend it for θ∥ ≠ θ⊥."*

A documented, validated, field-aligned bi-Kappa implementation is defensible on exactly this
basis. A novel *construction* is not. **The ZUM2026 "easily extend" sentence must appear in our
text**, not be finessed — a returning referee will find it, and quoting it ourselves converts a
vulnerability into evidence of care.

One further usable point: Abdul & Mace (2015) requires a χ²_ν deviate with generally
non-integer ν — a true Gamma variate — and **never says how it is generated**. No appendix, no
gamma-generator citation. Zenitani (2025) makes the same portability point. Safe to state, and
scoped to *portability*, not to correctness or efficiency.

**Not verifiable from Consensus, must be checked against the primary source before use**
(the checklist item "verify every added bibliographic field and DOI" applies):

- Regularized Kappa: Scherer, Fichtner & Lazar (2017), *EPL* **120**(5), 50002 — needed for §1
  fact 6. (Not *A&A*; that is why the earlier search failed.)
- ~~Foundational solar-wind proton anisotropy~~ — **resolved in §2a**, except Marsch 1982 full text.
- Gamma variate generation: Marsaglia & Tsang (2000), *ACM TOMS* 26, 363 — the algorithm behind
  `std::gamma_distribution`, and the reason the "constant time" claim needs qualification.
- Kroese, Taimre & Botev (2011), *Handbook of Monte Carlo Methods* — standard reference for the
  Student-t / chi-squared construction, alongside `devroye1986nonuniform` (already present).

Note: `devroye1986nonuniform`, `birdsall1985plasma`, `hockney1981computer`, `verboncoeur2005review`
and `cartwright2000loading` are **already in `refs.bib` but never cited in the text**. R1.7 asks for
exactly these. The rewrite below uses all five.

---

## 3. Low-κ numerical-risk audit — measured, not assumed

Run `python3 python/lowkappa_risk_audit.py` to reproduce. `N = 2×10⁶` per κ, seed 20260815,
numpy 2.3.3. This audit was run **because** the first draft of this document asserted that the
bounded transform `Y = T/(1+T) ~ Beta(3/2, κ−1/2)` rescues the Gamma-ratio path at low κ.
It does not. The assertion was wrong and is withdrawn.

### 3a. The bounded transform is not neutral — it is the worst route tested

Non-finite radii per 2×10⁶ draws:

| κ | b = κ−1/2 | X2 == 0 | `sqrt(X1/X2)` | `sqrt(X1)/sqrt(X2)` | **bounded Y** |
|---|---|---|---|---|---|
| 1.50 | 1.000 | 0 | 0 | 0 | 0 |
| 0.75 | 0.250 | 0 | 0 | 0 | **222** |
| 0.60 | 0.100 | 0 | 0 | 0 | **52,205** |
| 0.55 | 0.050 | 0 | 0 | 0 | **322,318** |
| 0.51 | 0.010 | 1,193 | 1,670 | 1,193 | **1,388,113** |
| 0.505 | 0.005 | 48,036 | 57,337 | 48,036 | **1,666,853** |

At κ = 0.55 the direct Gamma ratio fails on **zero** of two million draws while the bounded
route fails on **16%** of them. The mechanism is rounding, not underflow: `Y = X1/(X1+X2)`
evaluates to exactly `1.0` as soon as `X2/X1 < ε/2 ≈ 1.1e-16`, whereupon `T = Y/(1−Y)` divides
by zero. `X2` itself remains representable down to `4.9e-324`. Closed form, confirmed against
the counts above to three digits:

```
P(bounded route fails)  ≈ (ε/2)^b / Γ(1+b)        threshold 1.1e-16
P(direct route fails)   ≈ (4.9e-324)^b / Γ(1+b)   threshold 4.9e-324
```

So the bounded variable fails roughly `10^(307.6 b)` times more often. It discards ~292 decades
of the dynamic range of `T`. Your diagnosis was right; the magnitude is worse than "relocates it."

**What survives — less than first written here.** The claim that `Y` is at least an excellent
*validation* statistic was itself too generous, and Experiment 1 measured it. `Y` is bounded and
cannot overflow, but boundedness is not the property that matters: for small `b = κ−1/2` the mass
of `Y` piles up against 1, where doubles have no relative resolution left. At κ = 0.55, **16.4% of
`Y` values round to exactly 1.0**, and a KS test against `Beta(3/2, b)` then reports `√n·D = 51.8`
— a catastrophic false failure that is measuring rounding, not the sampler.

The fix is the *complement*, `W = 1/(1+T) ~ Beta(κ−1/2, 3/2)`. It is the same bijection of `T`
carrying identical information in exact arithmetic, but it puts that mass near 0, where relative
resolution runs to ~1e-308. The same data tested on `W` gives `√n·D = 0.751`, a clean pass. See
`experiments/exp1_radial_directional/results/exp1_table.md`.

So the honest statement is: the orientation of the bounded transform is not cosmetic, and `Y` is
the wrong orientation. `W` is the correct bounded diagnostic. This supersedes R1.3 item 2 in the
reviewer matrix, which names `Y`.

### 3b. One real defect found in our own code — **FIXED, and now regression-tested**

`cpp/bi_kappa_distribution.H` computed `std::sqrt(x1 / x2)`, forming `T` before the root.
At κ = 0.51 this loses ~500 draws per 2×10⁶ to *intermediate* overflow of `T`, whose true radii
are `1.4e154 … 6.7e161` — comfortably representable. `sqrt(x1)/sqrt(x2)` recovers every one of
them at identical cost. Free, strict improvement, independent of anything else here.

**Status: landed** at `cpp/bi_kappa_distribution.H:273`, with Experiment 4 quantifying it on
*identical variates* (so the difference is the fix and nothing else) and regression tests
`test_radius_formation` A1–A4 in `cpp/test_suite.H` guarding it:

| Precision | κ | loss, `sqrt(x1/x2)` | loss, `sqrt(x1)/sqrt(x2)` |
|---|---|---|---|
| float | 0.55 | 1.22×10⁻² | 5.7×10⁻³ |
| float | 0.60 | 1.5×10⁻⁴ | 2.9×10⁻⁵ |
| double | 0.51 | 8.2×10⁻⁴ | 5.7×10⁻⁴ |

Verified to agree to **2 ulp** (max relative difference 1.95·ε) at κ = 2, where neither
formation can overflow — i.e. no measurable change of law at ordinary κ, checked by magnitude
rather than by how often the two round differently.

### 3c. Log-domain only helps if the generator provides it

`exp(0.5*(log X1 − log X2))` applied to an already-drawn `X2` buys nothing: `log(0) = −inf`.
Identical failure counts to `split`. But both numpy and libstdc++ realise a small-shape gamma as
the Marsaglia–Tsang boost `X = Ga(b+1)·U^(1/b)`, so

```
log X2 = log Ga(b+1) + (1/b)·log U
```

is exactly computable — `(1/b)·log U` is O(−10⁵), far inside double range — even where `X2`
flushes to zero. Implemented this way: **zero** non-finite results at every κ tested, and KS
against numpy's gamma gives p = 0.63 / 0.12 / 0.67 at b = 0.5 / 0.05 / 0.01.

This also lets failures be split into two categories that the current code conflates:

| κ | spurious (true R representable, lost to underflow) | honest (true R > 1.8e308) |
|---|---|---|
| 0.510 | 1,623 | 1 |
| 0.505 | 55,587 | 1,601 |
| 0.501 | 498,782 | 484,782 |
| 0.5001 | 128,197 | 1,734,894 |

The log-space route eliminates the *spurious* column. It cannot eliminate the *honest* one — no
reformulation can, because the sample genuinely exceeds double range. The honest claim is
therefore "converts silent NaN/Inf into a diagnosable, quantified, documentable condition,"
not "solves the zero-division problem."

### 3d. The findings that most affect positioning

**In double precision, the range R1.3 actually asks about is already clean.** Zero failures at
κ ≥ 0.55 across 2×10⁶ draws. The hazard Zenitani et al. flag bites only below κ ≈ 0.52, and
κ < 0.55 has no established physical motivation in space plasmas, where κ is typically 1.5–10.
So what `1/2 < κ ≤ 3/2` needs is **validation**, which we lack, not a **fix**, which isn't
needed. That is the measurement-based argument for your demotion, and it is stronger than the
argument from caution.

**In single precision it is a different story, and this cuts the other way.** The underflow
threshold moves from `4.9e-324` to `1.4e-45`, so `P(fail) ≈ (1.4e-45)^b`:

| κ | b | float32 failures / 2×10⁶ | predicted |
|---|---|---|---|
| 0.75 | 0.25 | 0 | 6.1e-12 |
| 0.60 | 0.10 | 64 | 3.3e-5 |
| 0.55 | 0.05 | 11,427 | 5.7e-3 |
| 0.51 | 0.01 | 710,993 | 3.6e-1 |

Prediction and measurement agree to two digits. At κ = 0.55 in single precision **one draw in
175 fails**, and single-precision particle data is common in PIC codes. So C5 may have more life
in it than §3a alone suggests — but in the `float` pipeline, not the `double` one, which is a
much more specific and more checkable claim than the one I originally wrote.

*Caveat:* this used numpy's double generator cast to `float32`, so the cast supplies the
underflow. A genuine single-precision Gamma would underflow during generation instead. The
threshold arithmetic is the same either way, but Experiment 4 Q1 must test a real `float`
pipeline before any of this is stated in the paper.

### 3e. Experiment 4 — the caveat in §3d is discharged, and the prediction holds

Experiment 4 tested a genuine `float` pipeline: `bi_kappa_distribution<float>` instantiating
`std::gamma_distribution<float>`, every operation in single precision, nothing drawn in double
and cast. 10⁶ draws × 5 seeds per configuration, both libc++ and libstdc++.

**The NumPy-cast predictions in §3d are confirmed to two digits by the real pipeline:**

| κ | §3d predicted | Exp 4 measured (libc++) | Exp 4 measured (libstdc++) |
|---|---|---|---|
| 0.75 | 0 | 0 | 0 |
| 0.60 | 3.3×10⁻⁵ | 2.9×10⁻⁵ | 3.1×10⁻⁵ |
| 0.55 | 5.7×10⁻³ | 5.7×10⁻³ | 5.7×10⁻³ |
| 0.51 | 3.6×10⁻¹ | 3.56×10⁻¹ | 3.56×10⁻¹ |

So the argument in §3d stands as written, and **"one draw in 175 fails at κ = 0.55 in single
precision" is now a measurement of the released code**, not an extrapolation from a cast.

Three further results change the positioning:

1. **No standard-library dependence.** libc++ and libstdc++ agree to within seed noise at every
   configuration. Both realise small-shape Gamma through the same boost identity, so there is
   nothing to report here — a negative result, and worth stating as one, because it is the
   thing a referee would otherwise ask about.
2. **§3c's log-domain route is confirmed and validated, but is NOT landed.** It removes the
   *spurious* column entirely, leaving only honest overflow (float κ=0.55: 0.99430 → 0.99986
   finite). Validated 12/12 on KS **and** Cramér–von Mises against `W ~ Beta(κ−1/2, 3/2)` at
   κ = 0.55, 0.75, 2, 5. It is not landed because it would replace `std::gamma_distribution`
   with a hand-rolled generator, changing the RNG stream and breaking seed-for-seed
   reproducibility against every prior release. **That is an authors' decision, not a bug fix.**
3. **The cap masks the failure rather than preventing it.** Non-finite draws returned to the
   caller: **0 in every configuration**, including where 93% of internal attempts are
   non-finite. A non-finite value can never satisfy the box predicate, so the rejection loop
   silently redraws it. At κ = 0.51 in float, 36% of internal attempts are non-finite and the
   user never learns of it. This must be written as failure *hidden* by truncation.

**Net effect on C5: PROMOTE, narrowly.** The three conditions of the decision rule are met. But
ZUM2026 §4 already says *"caution is needed"* on `1/2 < κ ≤ 3/2`, so the defensible framing is
that we **quantified an acknowledged caution and removed one avoidable cause of it** — not that
we discovered anything. Tier 2, characterization plus a landed fix. No novelty claim.

---

## 4. The contribution statement

This is the piece to agree on before any prose gets written, because everything else follows from it.

**Retire entirely:** discovery of the Gamma-ratio construction; "resolves computational
bottlenecks"; "fast"; "prohibitively low acceptance rates"; "constant time per sample";
any implication that rejection sampling was the only alternative.

**Tier 1 — establishable, and the paper stands on these alone:**

- **C1 — Simulation-integrable implementation.** Header-only C++ with externally supplied,
  seeded RNG state, plus a Python path, released open source under a citable DOI. Zenitani et al.
  release a Jupyter notebook; the gap is integration into compiled PIC/hybrid codes, not
  mathematics. *True today*, modulo the Python-packaging caveat in matrix §1.4. Wording
  constraint in the §5 drafting notes.
  **Reproducibility must be stated at two distinct strengths, and only one of them is currently
  earned.** (i) *Controlled seeded RNG state and a reproducible workflow* — same seed, same
  binary, same output; true now. (ii) *Cross-platform bitwise reproducibility* — same seed across
  compilers and standard libraries; **not claimable**, because `std::gamma_distribution` is not
  specified to a particular algorithm and libstdc++ / libc++ / numpy differ, especially in the
  small-shape branch that Experiment 4 Q1 is about. Claim (i); claim (ii) only if Q1 happens to
  demonstrate it, which is unlikely.
- **C2 — Arbitrary magnetic-field-frame loading.** Orthogonal rotation from the field-aligned
  `(⊥,∥)` frame to global Cartesian, with invariance of the normalized radial law stated and
  tested. Not in Zenitani et al. *API true today; becomes a validation claim after Experiment 1's
  frame-invariance panel.*
- **C3 — Explicit uncapped/capped target-law separation.** Two named modes, two named target laws,
  the capped one defined as Eq. (20) conditioned on the component-wise box event. *Blocked until
  the C++ no-cap mode exists* (R1.1 item 5) — a small, well-understood code change.
- **C4 — Validation and reproducibility framework.** Direct validation of the central radial law
  `T = R²`, directional uniformity, and radial–direction independence — not Cartesian marginals
  alone — plus one committed driver regenerating every figure and table from recorded seeds,
  modes, `λ`, and `N`. Requires Experiments 1 and 2. **This is the substantive scientific
  deliverable of the revision**, and it is the one the referees asked for most directly.

**Tier 2 — candidate, to be established or dropped by the revision:**

- **C5 — Low-κ numerical robustness.** *Candidate numerical contribution, not a novelty claim.*
  §3 shows the original framing was wrong twice over: the bounded transform makes matters
  dramatically worse, and the range R1.3 actually asks about needs no fix at all. What may
  survive is narrower — a log-space small-shape Gamma construction that removes all spurious
  failures below κ ≈ 0.52 and makes the irreducible overflow fraction measurable, plus the
  `sqrt(x1)/sqrt(x2)` correction of §3b. Whether that clears the bar for a stated contribution is
  for Experiment 4 to decide, not for this document.
  **Guard:** the relation `T ~ β'(a,b) ⟺ T/(1+T) ~ Beta(a,b)` is textbook, and small-shape Gamma
  handling has a long literature (Ahrens & Dieter 1974; Best 1983; Devroye 1986 ch. VII; Zenitani
  2024b). Nothing in this area may be called novel before that audit is done.

**Honesty ledger.** C1 and C2 are true today. C3 is false until the API changes. C4 is a promise
until Experiments 1–2 run. C5 may not survive at all. The referees return holding their original
reports; the contribution list must not get ahead of this ledger.

---

## 5. Draft Introduction

Five paragraphs replacing the current three, plus an explicit contribution list. Paragraphs 1–2
are lightly revised from the existing text (they were never the problem); 3–5 are new.
`\cite` keys marked `NEW:` do not exist in `refs.bib` yet.

```latex
\section{Introduction}\label{sec:intro}

% ---- P1: unchanged in substance; this paragraph was never the problem. ----
It is well known that without external forces, particle systems in thermal
equilibrium stabilize into a Maxwellian velocity distribution. In numerous
physical environments---ranging from the solar wind and planetary
magnetospheres to fusion plasmas---particle velocity populations, however,
frequently exhibit characteristic non-Maxwellian
tails~\cite{livadiotis2013toolbox}. This feature is empirically well fitted by
a Kappa ($\kappa$) or bi-Kappa (anisotropic)
distribution~\cite{vasyliunas1968survey}, which reduces to the Maxwellian as
$\kappa\to\infty$ but retains a power-law suprathermal tail at finite $\kappa$.
We use the Kappa form here as an empirical description of these populations.
Several theoretical routes to it have been proposed for collisionless plasmas,
in which collective interaction with long-range fields rather than binary
collisions governs the
dynamics~\cite{livadiotis2013toolbox,livadiotis2018universe}; the present work
does not depend on adopting any particular one.

% ---- P2: R1.7. Two foundational refs added; heating claim now sourced. ----
Bi-Kappa distributions are widely used to characterize space plasma
populations~\cite{verscharen2019multiscale}, especially protons and minor
(heavy) ions~\cite{yoon2023bikappa}. Temperature anisotropy is a persistent
feature of the expanding solar wind, and the observed distribution of
$T_\perp/T_\parallel$ is bounded by kinetic instability
thresholds~\cite{NEW:marsch1982helios,NEW:hellinger2006anisotropy}; the
perpendicular heating that sustains it is commonly attributed to the
dissipation of pre-existing turbulence at ion-gyroradius
scales~\cite{NEW:chandran2010perpendicular,lazar2022temperature}. Simulations
intended to represent such anisotropic suprathermal populations therefore
require an initial loading that carries both the anisotropy and the power-law
tail; where the suprathermal component is not the object of study, a
bi-Maxwellian loading remains appropriate.
Maxwellian loaders are standard in particle-in-cell (PIC) and hybrid
codes~\cite{birdsall1985plasma,hockney1981computer,verboncoeur2005review,cartwright2000loading},
but the bi-Kappa case requires dedicated machinery.

% ---- P3: R2.A2 + R2.A3. Says plainly what is hard AND what is not. ----
% This paragraph is where we stop overclaiming difficulty.
It is worth stating precisely which parts of this problem are difficult, since
the absence of a closed-form inverse cumulative distribution function is often
taken as the whole of it. It is not. Numerical inverse-transform methods are
available for broad classes of one- and low-dimensional distributions---
Chebyshev-based inverse transform sampling has been applied directly to plasma
distribution functions in one and two
dimensions~\cite{NEW:an2022chebsampling}---and generic acceptance--rejection is
available whenever a workable envelope exists, with the caveat that envelope
quality degrades in three dimensions~\cite{devroye1986nonuniform}. Nor does
anisotropy alter the stochastic core: for the standard gyrotropic bi-Kappa
parameterization adopted here, the affine rescaling
$\mathbf{V}=\mathrm{diag}(\sqrt{\kappa}\,\theta_\perp,\sqrt{\kappa}\,\theta_\perp,\sqrt{\kappa}\,\theta_\parallel)\,\mathbf{U}$
reduces the anisotropic problem to the same isotropic dimensionless radial law
as the standard Kappa distribution, so within this parameterization anisotropy
does not introduce an additional random-variate sampling problem. What does remain
genuinely awkward is practical: the second moment diverges for
$\kappa\le3/2$, so variance-based verification is unavailable exactly where the
tail is heaviest; finite-precision behavior degrades as $\kappa\to1/2^{+}$; and
a loader must be embedded in a simulation's coordinate conventions, magnetic
field geometry, and random-number infrastructure rather than existing as a
standalone formula.

% ---- P4: R2.A1. THE survey paragraph. Concedes priority explicitly. ----
Several routes to Kappa-distributed variates are already established.
Abdul and Mace~\cite{abdul2014kappa} identified the equivalence between the
Kappa distribution and the multivariate Student-$t$ distribution and adapted
the standard Student-$t$ generator accordingly, subsequently applying it to
multi-dimensional PIC loading~\cite{abdul2021whistler}. Zenitani and
Nakano~\cite{zenitani2022relativistic,zenitani2023losscone} developed
procedures for relativistic Kappa and loss-cone distributions built on the
generalized Beta-prime distribution, expressing the radial variate as a ratio
of Gamma variates. Purely rejection-based
recipes~\cite{zenitani2025simple} and approximate inverse-transform
recipes~\cite{NEW:zenitani2026approximate} have also been proposed, the former
reporting acceptance efficiencies of $0.73$--$0.8$ for the isotropic Kappa
distribution. Most directly relevant to the present work, Zenitani, Usami and
Matsukiyo~\cite{zenitani2026nonmaxwellian} give Monte Carlo procedures for the
generalized $(r,q)$ distribution, which reduces to the standard bi-Kappa
distribution at $r=0$, $q=\kappa+1$. Their Beta-prime algorithm at those
parameter values---a radial draw from a ratio of two Gamma variates, an
isotropic direction on the unit sphere, and a parallel/perpendicular
scaling---is mathematically equivalent to the construction we present in
Sec.~\ref{sec:method}. We developed that construction independently, but we
make no claim of novelty for it, and we cite their work at each point of use.
This paper does not introduce the Gamma-ratio construction; it implements and
validates an established one.

% ---- P5: the gap + contributions. Nothing here outruns the evidence. ----
What remains open is not the mapping but its delivery and its verification. The
implementations we surveyed do not provide the combination of a header-only
\texttt{C++} interface, externally controlled random-number state, loading in an
arbitrary magnetic-field frame, and validation against the three-dimensional law
being sampled rather than against its one-dimensional marginals alone.
Accordingly, this paper contributes:
\begin{enumerate}
  \item an open-source, header-only \texttt{C++} bi-Kappa loader with
        deterministic and externally supplied random-number generation,
        together with a Python implementation, released under a citable DOI
        (Sec.~\ref{sec:code});
  \item loading in an arbitrary magnetic-field frame via an orthogonal
        transform that preserves the normalized radial law
        (Sec.~\ref{sec:transform});
  \item an explicit separation of the untruncated target
        distribution from an optional component-wise capped variant, with each
        target law stated and each mode independently characterized
        (Secs.~\ref{sec:method} and~\ref{sec:truncation}); and
  \item direct validation of the central radial law $T=R^{2}$, of directional
        uniformity, and of radial--direction independence, in addition to the
        Cartesian marginal and moment checks, together with a reproducible
        driver that regenerates every figure and table reported here from
        recorded seeds, modes, and parameters
        (Sec.~\ref{sec:verification}).
\end{enumerate}
% PLACEHOLDER, do not draft until Experiment 4 reports.  If and only if the
% low-kappa work yields a defensible result, add a fifth item here and a
% correspondingly scoped sentence about the shape < 1 regime.  See
% docs/introduction_rewrite_proposal.md Sec. 3 and Sec. 8 -- the first draft of
% this paragraph asserted a bounded-transform fix that measurement refuted.
This work is released in alignment with Open Science and the Heliophysics Open
Modeling Environment (HOME) initiative~\cite{corti2026openscience}.
```

**Drafting notes**

- P4's last two sentences are the load-bearing concession. They should stay adjacent and stay
  blunt — a hedged version reads worse to a referee who has read Zenitani et al., not better.
  Resist the temptation to soften "we make no claim of novelty for it."
- **The chronology does not belong in the manuscript.** P4 says only "developed independently";
  the November 2025 development record and the "prior to publication" timeline go in the
  **response letter**, where they do real work — they explain why the original submission did not
  cite a paper published after it. In the body they add nothing scientific and read as a priority
  defence, which is the impression the whole revision is trying to avoid.
- P1 presents the Kappa form as an empirical description and does not commit to a derivation;
  P2 says simulations of *suprathermal anisotropic* populations need this loading, rather than
  implying bi-Maxwellian initialization is generally wrong. Anisotropy alone does not imply a
  power-law tail, and a referee would catch that inference.
- P5's numbered list is ordered by what is provable today: items 1–2 are true now, item 3 lands
  with the API change, item 4 with Experiments 1–2. Nothing about low κ appears, by design — see
  the placeholder comment. If any item slips, delete it; do not weaken its wording and keep it.
- P5 deliberately says "the implementations we surveyed do not provide the combination of …"
  rather than "no simulation-ready package exists." The first is verifiable from the survey in
  P4; the second is an exhaustive prior-art claim we cannot support and would be asked to defend.
  The same care applies to the Sec. VI software claims.
- P3 answers R2.A3 inline, scoped to the *gyrotropic* parameterization adopted here rather than
  to bi-Kappa sampling in general, so no semantic gap is left open. The Background section can
  cross-reference it rather than repeat it.
- The word "exact" does not appear in this draft, and neither does "rejection-free". Both should
  reappear only in Sec. III, tightly scoped, per R1.1.

---

## 6. Deletions required elsewhere (Introduction-adjacent)

Sentence-level, so the terminology audit can be checked off mechanically:

| Location | Current | Action |
|---|---|---|
| Abstract | "We present an exact sampling algorithm" | → "We present an open-source implementation and validation of a bi-Kappa loading algorithm" |
| Abstract | "runs in constant time per sample" | delete; replace with fixed *high-level* variate count, or measured throughput from Experiment 3 |
| Abstract | "at the cost of a negligible rejection rate" | delete "negligible" until Experiment 2 reports `λ` and the rejected fraction |
| Intro ¶3 | "Analytical inversion … is intractable" | → numerically tractable but not the only route; cite An et al. |
| Intro ¶3 | "prohibitively low acceptance rates" | delete; Zenitani (2025) reports 0.73–0.8 for the isotropic case |
| Intro ¶3 | "To resolve these computational bottlenecks … fast and accurate algorithm" | delete wholesale |
| Intro ¶3 | "ensures rigorous adherence to prescribed velocity limits … correctly reproduce the expected truncated suprathermal tails" | delete; this asserts the capped mode is correct for Eq. (20), which is the R1.1 error |
| Sec. VII | "We presented an exact, rejection-free algorithm" | → per matrix §4 replacement wording |

---

## 7. `refs.bib` additions

Four to add, all with DOIs to verify against the primary source before submission:

- `an2022chebsampling` — An, X., et al. (2022), *JGR Space Phys.* 127, e2021JA030031,
  doi:10.1029/2021JA030031
- `zenitani2026approximate` — Zenitani, S., et al. (2026), *Earth, Planets and Space*,
  "An approximate Kappa generator for particle simulations"
- `chandran2010perpendicular` — Chandran, B. D. G., et al. (2010), *ApJ* 720, 503
- `scherer2017regularized` — Scherer, Fichtner & Lazar (2017), *A&A* — for the Background
  section's treatment of the regularized Kappa as the principled alternative to our box cap

Plus, pending a primary-source check: `marsch1982helios` / `hellinger2006anisotropy`,
`marsaglia2000gamma`, `kroese2011handbook`.

If any low-κ claim survives Experiment 4, the small-shape Gamma literature must also be cited
at the point of use: Ahrens & Dieter (1974), Best (1983), Devroye (1986 ch. VII §2.6) — all
named in Zenitani et al. §2.2 — plus Zenitani (2024b).

---

## 8. Experiment 4 — design

Revised in light of §3. The original scope ("does the bounded transform rescue low κ?") is
answered: no. The remaining questions are narrower and mostly about our own implementation.

**Q1. Where does the current C++ actually break?** Sweep κ ∈ {0.5001, 0.501, 0.505, 0.51, 0.55,
0.75, 1.0, 1.5}, `float` and `double`, `std::gamma_distribution` vs numpy. Record exact-zero
denominators, Inf, NaN, and max-attempt failures. The Python audit predicts zero failures for
κ ≥ 0.55 in double; **confirm this holds for libstdc++ and libc++ too**, since the small-shape
Gamma algorithm is implementation-defined and the `U^(1/b)` boost is not mandated by the standard.
**Priority within Q1 is the `float` pipeline**, where §3d measures ~1 failure in 175 at κ = 0.55
— inside the physically relevant range — against zero in double. Use a genuine single-precision
generator, not a cast.

**Decision rule for C5 hinges mostly on this sub-question.** A real failure at κ = 0.55 in
`float` is a defensible robustness contribution; a failure only at κ = 0.505 in `double` is a
footnote.

**Q2. Is the §3b `sqrt(x1)/sqrt(x2)` fix worth landing on its own?** It is free and strictly
better. Land it regardless of what else survives; verify no distributional change at moderate κ.

**Q3. Does a log-space small-shape Gamma preserve the target law?** The Python KS check
(p = 0.63 / 0.12 / 0.67) is necessary but not sufficient — it can only compare on the range
numpy can represent. Add a moment-based and a quantile-based check on `log R` itself, where the
comparison is not truncated.

**Q4. What is the honest supported range?** Produce the table §3c started: for each κ and
precision, the spurious-failure count (removable) and the honest-overflow count (not removable).
The deliverable is a documented supported (κ, precision) envelope and a defined behaviour outside
it, not a claim of having eliminated the failure.

**Q5. Does the cap interact?** In capped mode the enormous radii are rejected anyway, so the
failure may be invisible there. Check whether capping masks the condition rather than handling
it — a silent mask is a defect worth reporting.

**Decision rule, fixed in advance.** C5 becomes a stated contribution only if Q1 confirms a real
failure in the released code within a range someone would plausibly use, *and* Q3–Q4 show the
mitigation preserves the target law with a measurable improvement in supported range. If κ ≥ 0.55
turns out clean in every implementation tested, the honest outcome is a documented operating
range and a footnote — and the paper rests on Tier 1. Write the result either way.

---

## 9. Sequencing

Revised: freeze the contribution statement *after* the experiments, not before. The first draft
of this document had it backwards and produced a claim that measurement then refuted.

1. ~~**Literature/novelty audit — read the primary sources.**~~ — **DONE 2026-08-17**, results in
   `docs/step1_claim_audit.md`. All ten PDFs in `paper/reference/` read. Three claims below are
   wrong as drafted (P4-a "multivariate"; P4-b mis-cited to 2014/2021 when the source is Abdul &
   Mace **2015**; P4-c understates Zenitani & Nakano 2022, whose Table I Algorithm 1-1 **is our
   construction** — so the concession moves from 2026 back to **2022**). The EPS paper is
   Zenitani & **Umeda**. Eight sources remain unobtained and block four claims. The sub-questions
   originally listed here are answered there; retained below for traceability.

   Original checklist — every *substantive* sentence in the related-work paragraph has to be
   checked against the paper it rests on, because these are the claims a returning referee is
   most likely to test. Verifying
   DOIs is the smaller half. Every *substantive* sentence in the related-work paragraph has to be
   checked against the paper it rests on, because these are the claims a returning referee is
   most likely to test:
   - An et al. (2022): how broad is the class Chebsampling actually covers? P3 currently says
     "broad classes of one- and low-dimensional distributions" — confirm that is what they
     demonstrate, and in how many dimensions.
   - Zenitani (2025): the 0.73–0.8 acceptance figure — over which κ range and which settings?
     P4 quotes it as the counter to "prohibitively low," so its scope matters.
   - Zenitani & Nakano (2022, 2023): which distribution machinery does each actually use? P4
     attributes the Beta-prime/Gamma-ratio route to 2022 and loss-cone methods to 2023; confirm.
   - Abdul & Mace (2014): confirm the Student-*t* equivalence is stated as P4 describes, and
     check whether 2014 or the 2015 companion is the right citation for multi-dimensional loading.
   - Scherer et al. (2017): confirm the regularized Kappa is the object Zenitani et al. §4 treats.
   Also complete the small-shape-Gamma prior-art check that C5 would depend on
   (Ahrens & Dieter 1974; Best 1983; Devroye 1986 ch. VII §2.6; Zenitani 2024b).
2. ~~**Fix the untruncated/capped API** (R1.1 item 5)~~ — **DONE 2026-08-15.**
   `bi_kappa_distribution<T>::no_cap()` disables the cap; `param().capped()` queries the mode;
   the uncapped path runs the core mapping once with no rejection loop. The §3b
   `sqrt(x1)/sqrt(x2)` correction landed in the same change. Existing test suite passes;
   README documents the two target laws. **Open decision:** the default is still `20.0`
   (non-breaking). Whether the default should become `no_cap()` — so the library's out-of-the-box
   behaviour matches the law the paper presents as primary — is a release-policy call for the
   authors, and it affects how Sec. VI can be worded.
3. ~~**Experiment 1** — radial, directional, anisotropic, frame-invariance validation.~~
   — **DONE 2026-08-17.** `experiments/exp1_radial_directional/`. 135 runs, 1.35e7 draws,
   uncapped, 5 seeds x 1e5 per configuration. The radial law `T = R² ~ β'(3/2, κ−1/2)`,
   directional uniformity, radial–direction independence and frame invariance all pass at
   every κ in {0.51 … 10}, **including the whole `1/2 < κ ≤ 3/2` range R1.3 asks about**
   (`√n·D` ∈ [0.75, 0.97] against an expected ≈0.87). Frame invariance holds to 1.1e-15
   relative. Non-finite draws occur only at κ = 0.51 (5.8e-4), confirming §3d's
   double-precision prediction on libc++ as well as numpy. **C4 is now evidence, not a promise.**
4. **Experiment 4** — per §8 above.
5. **Experiment 2** — truncation characterization.
6. **Experiment 3** — performance benchmark, optional but recommended (see §11).
7. **Freeze the contribution statement against actual results.**
8. Rewrite Introduction, Abstract, Summary; apply the §6 deletions in one pass.
9. Background/Related Work section (R2.B) — elaborates P4, does not restate it.
    Freeze the title here too (§10).
10. Point-by-point rebuttal.

The P1–P4 prose in §5 is stable and independent of every experiment — it is pure literature and
framing. It can be drafted at step 1 if useful. **P5 cannot**, and the placeholder comment in the
draft marks exactly where the dependency lies.

---

## 10. Title — decide at step 7, not now

The reviewer matrix judged the current title neutral because it makes no priority claim. That is
true as far as it goes, but it understates the risk: a second-round referee who reads
*"Sampling the Bi-Kappa Distribution"* will reasonably expect the contribution to be the sampler.
If the frozen contribution statement is implementation + validation, the title should say so.

Candidates, to weigh once Experiments 1–4 report:

- *Implementation and Validation of Bi-Kappa Sampling for Particle Simulations* — shortest, and
  accurate under every outcome currently plausible.
- *Open-Source Implementation and Validation of Bi-Kappa Loading for Particle Simulations* —
  more explicit, foregrounds the software contribution; "loading" is the term the PIC literature
  and Zenitani et al. use.
- Retain *Sampling the Bi-Kappa Distribution* — only if C5 survives strongly enough that a
  methods framing is still defensible.

Freeze this together with the contribution statement at step 7. Do not add "novel" or "new" under
any outcome.

---

## 11. On the benchmark (Experiment 3)

Worth running even though no superiority claim will be made, because there are now four
comparable routes — Beta-prime/Gamma-ratio, piecewise rejection, Student-*t*, and approximate
inverse transform — and a parameter-dependent comparison is genuinely useful to a reader
choosing among them. Two conditions: validate distributional correctness before timing anything,
and report whatever the numbers say. A benchmark that has to come out favourable is not a
benchmark.
