# Step-1 claim audit — every planned literature claim vs. its primary source

**Scope:** sequencing step 1 of `introduction_rewrite_proposal.md` §9. Each claim planned for the
revised Introduction (P2–P5), for §1 "what Zenitani et al. actually says", and for `refs.bib`
gets one verdict: **PASS** (stands as drafted), **MODIFY** (true but mis-scoped/mis-attributed),
**DROP** (not supportable), or **BLOCKED** (source not in hand — no verdict possible).

**Method.** All ten PDFs in `paper/reference/` were converted with `markitdown` and read. Every
verdict below cites the passage it rests on. Nothing here is from memory or from metadata.

**Headline:** 3 claims are wrong as drafted and must change; 1 reference is mis-authored; the
prior-art concession must move back **four years**, from 2026 to 2022. Eight sources are still
unobtained and block four further claims.

---

## 0. Verdict summary

> **Update, 2026-08-17.** A second round obtained five further primary sources. Two blockers
> are resolved and one of them **inverts a planned claim**; two new blockers are recorded. See
> §7. The verdicts below are updated in place; the original reasoning is retained.

| # | Claim | Verdict |
|---|---|---|
| P2-a | Solar-wind anisotropy bounded by instability thresholds | ~~BLOCKED~~ → **MODIFY** (§7.2) |
| P2-b | Turbulent perpendicular ion heating at gyroradius scales | **MODIFY** |
| P2-c | Maxwellian loaders standard in PIC/hybrid codes | PASS (unverified textbooks) |
| P3-a | Chebyshev inverse transform, "broad classes / low-dimensional" | **MODIFY** |
| P3-b | Rejection envelope quality degrades in 3-D → Devroye | **MODIFY** |
| P3-c | Second moment diverges for κ ≤ 3/2 | PASS |
| P3-d | Affine rescaling reduces anisotropic → isotropic core | PASS |
| P3-e | Finite-precision degrades as κ → 1/2⁺ | PASS |
| P4-a | A&M 2014 = *multivariate* Student-*t* equivalence | **MODIFY — error** |
| P4-b | A&M applied to multi-dimensional PIC loading | ~~BLOCKED~~ → **DROP — inverted** (§7.1) |
| P4-c | Z&N 2022/2023 built on generalized Beta-prime | **MODIFY — understated** |
| P4-d | Zenitani 2025 acceptance 0.73–0.8 | **MODIFY — scope** |
| P4-e | "Zenitani et al. 2026" approximate generator | **MODIFY — authorship** |
| P4-f | (r,q) reduces to bi-Kappa at r=0, q=κ+1 | PASS |
| P4-g | Their Algorithm 3.1 ≡ our construction, incl. anisotropy | PASS |
| §1-1 | Our "gap" premise contradicted in print | PASS |
| §1-5 | Low-κ zero-division caution | **MODIFY — incomplete** |
| §1-6 | Regularized Kappa = Scherer | **MODIFY — wrong year** |
| §1-7 | Their release is a Zenodo Jupyter notebook | PASS |

---

## 1. The three claims that are wrong as drafted

### P4-a — "multivariate" is wrong. **MODIFY.**

Draft: *"Abdul and Mace [abdul2014kappa] identified the equivalence between the Kappa distribution
and the **multivariate** Student-t distribution and adapted the standard Student-t generator
accordingly."*

`AbdulMace2014_CPC185_2383.pdf` demonstrates the equivalence for the **one-dimensional** Kappa
distribution. Their Eq. (2) is obtained by *"Integrating Eq. (1), the three-dimensional kappa
distribution, over two velocity space coordinates"*, and the generator they adopt is Bailey's
**univariate** polar transform, `b = sqrt(ν(a₁^(−2/ν) − 1))·cos(2πa₂)`, of which they note
*"only one deviate is generated at a time"*.

Independently corroborated by Zenitani (2025) §1: *"There exist t-generators that only require
uniform variates (Bailey 1994; Abdul & Mace 2014), however, they are only available for **1D and
2D** distributions."*

**Action:** delete "multivariate". A&M 2014 is a 1-D method paper.

### P4-b — the multi-dimensional loader is Abdul & Mace **2015**, which we have not read. **BLOCKED.**

Draft: *"…subsequently applying it to multi-dimensional PIC loading [abdul2021whistler]."*

Two independent sources attribute the 3-D/multivariate Kappa loader to a paper that is neither
in `refs.bib` nor in `paper/reference/`:

- Zenitani (2025) §1: *"Recent works (Abdul & Mace **2015**; Zenitani & Nakano 2022) employ this
  strategy"* — the strategy being *"to divide a multivariate normal distribution by a
  gamma-distributed random number."*
- Zenitani, Usami & Matsukiyo (2026) §3: *"Along the red line, one should use dedicated methods
  for the Kappa distribution (Abdul & Mace, **2015**; Zenitani, 2025; Zenitani & Nakano, 2022)."*
  The red line is the Kappa line, Eq. (12).

Full citation, from ZUM2026's reference list: **Abdul, R. F., & Mace, R. L. (2015).
One-dimensional particle-in-cell simulations of electrostatic Bernstein waves in plasmas with
kappa velocity distributions. *Physics of Plasmas*, 22(10), 102107, doi:10.1063/1.4933005.**

Note the tension: the title says *one-dimensional* PIC simulations, yet two independent groups
cite it for the multivariate loader. That cannot be resolved without reading it.

**Action:** do not write this sentence until A&M 2015 is obtained. `abdul2021whistler` is also
unread — the doc's literature table describes it as "bi-kappa electrons" while the `refs.bib`
title says "Kappa Velocity Distributions"; that needs checking too.

### P4-c — Z&N 2022 is understated, and it moves the concession back to 2022. **MODIFY.**

Draft: *"Zenitani and Nakano [2022, 2023] developed procedures for relativistic Kappa and
loss-cone distributions built on the generalized Beta-prime distribution, expressing the radial
variate as a ratio of Gamma variates."*

Three problems, one of them consequential.

**(i) Their headline method is rejection-based, which the draft omits.** Abstract of
`ZenitaniNakano2022_PoP29_113904.pdf`: *"A procedure for loading particle velocities from a
relativistic kappa distribution … is presented. It is based on the **rejection method** and the
beta prime distribution."*

**(ii) 2023 is not a Beta-prime paper.** `ZenitaniNakano2023` covers the Ashour-Abdalla–Kennel
loss-cone (subtracted Maxwellian) via a summation algorithm, the Dory-type via a gamma variate,
and a kappa loss-cone (KLC) construction — plus *"two transformation algorithms"* alongside
acceptance-rejection. Lumping it under "generalized Beta-prime" is inaccurate. Split the cites.

**(iii) The consequential one — their Table I, Algorithm 1-1 is our construction, published 2022.**
As preparation for the relativistic case they set out the *non-relativistic* isotropic Kappa
algorithm in full:

```
generate X1, X2 ~ U(0,1)
generate X3   ~ Ga(3/2, 2)
generate v²_ν ~ Ga(κ − 1/2, 2)
v  ← sqrt( κθ² · X3 / v²_ν )
vx ← v(2X1 − 1)
vy ← 2v·sqrt(X1(1−X1))·cos(2πX2)
vz ← 2v·sqrt(X1(1−X1))·sin(2πX2)
```

That is our Eqs. (23)–(29) in the isotropic limit — same Gamma ratio, same shape parameters, same
two-uniform spherical scatter. The proposal's §1 dates the overlap to March 2026. **It is
actually November 2022**; ZUM2026 adds the (r,q) generalization and the anisotropic
`(θ⊥, θ⊥, θ∥)` scaling.

**Action:** the concession in P4 must name Zenitani & Nakano (2022) as the origin of the
Gamma-ratio Kappa construction, with ZUM2026 as the generalization that first covers the
anisotropic case in print. This makes the concession stronger, not weaker — and a returning
referee who checks will find it either way.

---

## 2. Scope and attribution corrections

### P4-d — the 0.73–0.8 figure is real but κ > 3/2 only. **MODIFY.**

`Zenitani2025_RNAAS.pdf` verified in full. Pareto envelope, uniform variates only, genuinely 3-D
(*"a random number generator for a 3D Kappa distribution that only requires uniform variates"*).

The efficiency depends on the envelope index `n`:
- near-optimum `n = 2κ/3 − 1/5` → *"≈0.8 regardless of κ"*;
- recommended `n = κ/2` → *"≈0.73–0.8. It starts from 0.806 at κ = 1.5, 0.785 at κ = 2, 0.750 at
  κ = 5, and then it is asymptotic to √(πe)/4 ≈ 0.731."*

**The scope caveat the draft misses:** the paper states *"the standard Kappa distribution is
defined for κ > 3/2"*. Every quoted number sits at κ ≥ 1.5. We must not let it read as covering
the `1/2 < κ ≤ 3/2` regime our own R1.3 work targets.

### P4-e — wrong authors. **MODIFY.**

Not "Zenitani et al." The paper is **Zenitani, S., & Umeda, T. (2026), *Earth, Planets and Space*
**78**:119, doi:10.1186/s40623-026-02465-0** — two authors, verified on the title page. Method as
described: q-exponential approximation to the CDF, inverse transform, *"practically accurate
results, in particular for κ < 4"*, GPU-oriented. The word "approximate" is theirs, so the
draft's "approximate inverse-transform recipes" is fair.

Fix the planned bib key `zenitani2026approximate` before it propagates.

### P3-a — state their actual scope. **MODIFY.**

`An2022` abstract: *"Chebsampling, to sample general distribution functions in **one and two
dimensions**."* The draft's "broad classes of one- and low-dimensional distributions" is vaguer
and slightly broader than what they demonstrate. Say "one and two dimensions".

The paper supports our point regardless — it states the classical limitation itself (*"The
application of inverse transform sampling is limited in practice, however, because it [requires
the inverse CDF]"*) and then shows it is numerically tractable anyway. That is exactly the
"inversion is not intractable" argument, sourced.

### P3-b — attribute the 3-D caveat to the paper that actually says it. **MODIFY.**

Draft attributes *"envelope quality degrades in three dimensions"* to `devroye1986nonuniform`,
which we have not obtained. ZUM2026 says it in print: acceptance-rejection *"becomes increasingly
inefficient, especially in three dimensions."* Cite ZUM2026 for the caveat, or obtain Devroye.

### P2-b — Chandran 2010 is right, the framing needs a hedge. **MODIFY.**

`Chandran2010_ApJ720_503.pdf` verified: *"Perpendicular Ion Heating by Low-Frequency Alfvén-Wave
Turbulence in the Solar Wind"*, ApJ 720, 503–515. Stochastic heating by AW/KAW turbulence at
gyroradius scales; *"The heating is anisotropic, increasing v⊥² much more than v∥²"* — precisely
the mechanism P2 wants.

Two scope facts the draft should respect: the result is for **β ≲ 1**, and the paper explicitly
notes that at β ≳ 1 Landau and transit-time damping instead drive **parallel** heating. So
"commonly attributed" overstates a single-mechanism consensus. Write it as one established
mechanism, with the β condition, rather than as the attribution.

### §1-5 — the low-κ caution is more permissive than the proposal records. **MODIFY.**

The quoted §3 sentence is verbatim correct: *"Below the blue curve … q ≤ 1 + 3/(2(1+r)), we
recommend the piecewise rejection method, because the beta-prime method may encounter the
zero-division problem."* On the Kappa line this is κ ≤ 3/2.

But the proposal omits their §4 summary, which is notably softer:

> *"The piecewise rejection method is the only choice for 0 < κ ≤ 1/2. **Either method can be
> used for 1/2 < κ ≤ 3/2, but caution is needed** for the post-rejection method. For κ > 3/2, we
> recommend the post-rejection method, which outperforms the piecewise rejection method."*

**This is good news for C5, and it reframes it.** ZUM2026 do not declare the Gamma-ratio route
unusable on `1/2 < κ ≤ 3/2`; they flag a qualitative hazard and advise caution. Our §3 audit
measured *where the hazard actually bites* — zero failures at κ ≥ 0.55 in double over 2×10⁶
draws, onset below κ ≈ 0.52, and ~1 failure in 175 at κ = 0.55 in single precision. That is a
**quantification of an acknowledged caution**, not a contradiction of it and not a new discovery.
It is a modest, defensible, and checkable contribution — and it is a better claim than the one
the proposal currently holds in reserve.

Their mechanism statement also matches ours: *"the gamma distribution has non-zero density at
x = 0. This may cause a zero-division problem"* (their §4, for shape κ − 1/2 ≤ 1).

### §1-6 — regularized Kappa is Scherer **2017**, not 2019. **MODIFY.**

ZUM2026 §4 opens: *"The regularized Kappa distribution (**Scherer et al., 2017**) is essentially
a Kappa distribution with a high-energy cutoff."* Their reference list gives **Scherer, K.,
Fichtner, H., & Lazar, M. (2017). Regularized κ-distributions with non-diverging moments.
*Europhysics Letters*, 120(5), 50002** — still not obtained.

The two Scherer 2019 papers we *do* have are related but are **not** ZUM's object:
- `SchererLazarHusidicFichtner2019_ApJ880_118` — *Moments of the Anisotropic Regularized
  κ-distributions*. Genuinely more apt for us, since our cap is a cutoff on a **bi**-Kappa, and it
  supplies the moments. Cite it for that, but never as "the distribution Zenitani et al. treat".
- `SchererFichtnerFahrLazar2019_ApJ881_93` — *On the Applicability of κ-distributions*. Useful
  and slightly awkward: for **κ < 2** the standard κ-distribution's pressure gets a significant
  *"unphysical contribution … from unrealistic, superluminal particles"*. This strengthens §3d's
  "κ < 0.55 has no established physical motivation" — the interpretive trouble starts as high as
  κ < 2 — but it also means we should not lean on physical motivation for low κ at all. Frame the
  low-κ work as numerical-robustness and validation-completeness, not physical relevance.

---

## 3. Claims that pass unchanged

- **P4-f, P4-g** — the core concession. Eq. (12) reduction at `r = 0, q = κ + 1` confirmed
  (the "red line" in their Figure 1b). Algorithm 3.1 confirmed verbatim, anisotropic scaling
  included in the algorithm box: `v∥ ← θ∥Rx(2U3−1)`, `v⊥1 ← 2θ⊥Rx√(U3(1−U3))cos(2πU4)`,
  `v⊥2 ← …sin(2πU4)`. There is no daylight, exactly as the proposal says.
- **§1-1** — verified verbatim at ZUM2026 line 55: *"numerical procedures for Kappa distributions
  have been largely unknown to the community, until Abdul and Mace (2014, 2015) used its relation
  to Student's t distributions. More recently, Zenitani and Nakano (2022, 2023) proposed numerical
  procedures…"* Our "gap" premise is contradicted in print, by name.
- **§1-7** — *"The Jupyter notebook for this article is archived in Zenodo (Zenitani et al.,
  2026)"*, doi:10.5281/zenodo.17148070. No C++ library, no field-frame handling. C1 survives.
- **P3-c** — corroborated by ZUM2026 §4: the regularized form allows *"the kappa index of
  0 < κ < 3/2 that are not accessible in the standard Kappa distribution."*
- **P3-d, P3-e** — our own mathematics and our own measurement; ZUM2026's Algorithm 3.1
  independently corroborates the affine-rescaling structure.

---

## 4. Findings that change the plan, beyond the claim list

**(a) A sourced, quantitative replacement for "prohibitively low acceptance rates."**
The proposal says delete the phrase. We can do better than deletion, with two papers now read:

- Abdul & Mace (2014) §2 measured `p(accept) ≈ 1/13 ≈ 0.077` using *"a simple uniform (rectangle)
  distribution"* as envelope, and found their transform method *"roughly 20 times faster than the
  accept–reject method … over the full range of N tested."*
- Zenitani (2025) reaches 0.73–0.8 with a **Pareto** envelope on the same target.

The honest, defensible sentence is therefore: acceptance efficiency is governed by envelope
quality — a naive rectangular envelope gives ≈8%, a well-matched Pareto envelope ≈73–80%. That
answers R1.7's "cite a standard reference and state when acceptance becomes poor rather than
making an unbounded assertion" with measured numbers from the primary literature.

**(b) A published claim that the rejection route may be *cheaper* than ours.**
Zenitani (2025) §3: the standard Kappa generator (A&M 2015; Z&N 2022) *"requires three normal
variates and one gamma variate per particle"*; via Marsaglia–Tsang that is *"four normal variates
and one uniform variate"*; his method calls *"≈4.5–4.7 uniform variates per particle. Therefore,
it would be computationally less expensive than the standard method."*

This directly constrains R1.4 and Experiment 3. We must not claim a speed advantage, and if we
benchmark, this is the specific published claim the benchmark is testing.

**(c) A real convention gap that is ours to occupy.**
Z&N 2022 Table I is captioned *"κ > 3/2 is required"*; Zenitani 2025 states *"the standard Kappa
distribution is defined for κ > 3/2"*. Both follow from the `θ² = 2[(κ−3/2)/κ](T/m)` thermal-speed
convention (Abdul & Mace 2014 Eq. 1). Our manuscript normalizes on `κ > 1/2`. That is a genuine
convention difference, it is exactly what R2.B demands be spelled out, and it is the cleanest
statement of where our validation adds something: the published Gamma-ratio Kappa loaders are
presented on `κ > 3/2`, and nobody has validated the construction below it.

**(d) Small-shape Gamma prior art, full citations captured** from ZUM2026's reference list, for
use only if C5 survives Experiment 4: Ahrens & Dieter (1974), *Computing* 12(3), 223–246,
doi:10.1007/BF02293108; Best (1983), *Computing* 30(2), 185–188, doi:10.1007/BF02280789;
Devroye (1986), Springer.

---

## 5. Still blocked — sources not in hand

No claim resting on these may be written yet.

| Priority | Source | Blocks |
|---|---|---|
| ~~1~~ | ~~Abdul & Mace (2015)~~ | **OBTAINED 2026-08-17 — see §7.1** |
| ~~2~~ | ~~Abdul, Matthews & Mace (2021)~~ | **OBTAINED 2026-08-17 — see §7.1** |
| ~~3~~ | ~~Hellinger et al. (2006) / Bale et al. (2009)~~ | **OBTAINED 2026-08-17 — see §7.2.** Marsch et al. (1982) full text still missing. |
| 4 | Scherer, Fichtner & Lazar (2017), *EPL* **120**(5), 50002 | §1-6 at its origin |
| 5 | Marsaglia & Tsang (2000), *ACM TOMS* **26**, 363 | the "constant time" qualification |
| 6 | Devroye (1986) | P3-b, if not re-attributed to ZUM2026 |
| 7 | Kroese, Taimre & Botev (2011) | Student-*t*/χ² construction reference |
| 8 | Birdsall & Langdon; Hockney & Eastwood; Verboncoeur (2005); Cartwright (2000) | P2-c — low risk, textbook-level statements |

Items 1–3 are the ones that block prose. Items 4–8 are verification of statements we are
confident about.

---

## 6. Net effect on the revision

1. **The concession moves from 2026 to 2022.** P4 must name Zenitani & Nakano (2022) Algorithm 1-1
   as the origin of the Gamma-ratio Kappa construction, with ZUM2026 as the anisotropic
   generalization. Rewrite P4 accordingly; the "we make no claim of novelty" sentence stays and
   gets, if anything, blunter.
2. **Three citation errors are fixed** before drafting: no "multivariate" on A&M 2014; A&M 2015
   obtained before any multi-dimensional-loading sentence; Zenitani & **Umeda** on the EPS paper.
3. **C5's framing improves.** ZUM2026 say "caution is needed" on `1/2 < κ ≤ 3/2`, not "unusable".
   Quantifying that caution — where it bites, in which precision — is a legitimate contribution
   and needs no novelty claim. Experiment 4's decision rule is unchanged.
4. **A new, cleaner statement of where we add value:** the published Kappa loaders are presented
   on `κ > 3/2` under a thermal-speed convention that requires it. Our validation targets
   `κ > 1/2`. That is factual, checkable, and does not claim anyone was wrong.
5. **Tier 1 (C1–C4) is untouched by everything above.** Nothing found in ten primary sources
   threatens the implementation, frame-transform, mode-separation, or validation contributions.
   The paper still stands on them.

---

# 7. Round 2 — 2026-08-17

Five further primary sources obtained. Two of the three prose-blocking items are resolved.
One resolution **inverts a claim we had planned to make**, which is the most important finding
in this round.

## 7.1 Abdul & Mace (2015) and Abdul, Matthews & Mace (2021) — OBTAINED

Both read in full. Also obtained: two open-access UKZN theses by the same author that document
the loader more fully than either journal paper.

| File | Citation |
|---|---|
| `AbdulMace2015_PoP22_102107.pdf` | Abdul & Mace (2015), *Phys. Plasmas* **22**(10), 102107, doi:10.1063/1.4933005 |
| `AbdulMatthewsMace2021_PoP28_062104.pdf` | Abdul, Matthews & Mace (2021), *Phys. Plasmas* **28**(6), 062104, doi:10.1063/5.0047638 |
| `Abdul2013_MScThesis_UKZN10413-12288.pdf` | Abdul (2013), MSc, UKZN, hdl:10413/12288 |
| `Abdul2018_PhDThesis_UKZN10413-22458.pdf` | Abdul (2018), PhD, UKZN, hdl:10413/22458 |

**The tension noted in §1 is resolved, and not in the direction we assumed.**

1. **What it samples:** the **isotropic three-dimensional** Kappa. Not a bi-Kappa, not a
   product of 1-D Kappas. Abstract: *"isotropic three-dimensional kappa velocity
   distributions"*; §III B is titled *"Generating random deviates following the 3D kappa
   velocity distribution"*.
2. **Dimensionality of the method: genuinely trivariate.** Eq. (22) is
   `X = μ + σ√(ν/χ²_ν)·Z` with `Z` a `p`-dimensional standard normal vector and **a single
   `χ²_ν` deviate shared across all three components** — the shared scale mixture is exactly
   what makes it a true trivariate Kappa. `ν = 2κ − 1`, `σ² = κθ²/(2κ−1)`.
3. **"One-dimensional" in the title is simulation geometry only.** §VI opens *"Results of 1D3V
   particle-in-cell simulations"*; the abstract says *"one-and-two-halves dimensional"*. One
   spatial dimension, three velocity components. Any statement equating that with 1-D
   velocity-space sampling is flatly wrong.
4. **Algorithm:** normal-over-chi-squared scale mixture, citing Kotz & Nadarajah (2004), Shaw &
   Lee (2008), Hofert (2013) as prior statistics literature. They explicitly *"avoid the
   computational waste produced through accept-reject sampling"*. They claim novelty only for
   *applying* it to PIC Kappa loading.
5. **Relation to A&M 2014:** 2014 is 1-D, uniform-variates-only (Bailey polar transform). 2015
   is the 3-D generalization, requires a Gamma/χ² generator, and **explicitly repudiates**
   applying the 1-D generator component-wise (p. 102107-5: a product of 1-D Kappas is
   *"innately anisotropic"*, *"not a monotonic function of kinetic energy"*, and *"can produce
   artificial instabilities when a stable plasma configuration is expected"*).
6. **Zenitani-side characterizations are all ACCURATE.** Seven separate claims across Zenitani
   (2025), ZUM2026, Zenitani & Umeda (2026) and Z&N (2022) were checked verbatim against
   AM2015. Every one checks out, including the operation count *"three normal variates and one
   gamma variate per particle"*. **Do not build any argument on a claimed misreading.**

**Abdul, Matthews & Mace (2021):** a **2D3V** electromagnetic GPU/CUDA PIC study running a
standard bi-Kappa (Summers & Thorne form, their Eq. 2), `T⊥/T∥ = 3.0`, **6.7×10⁷ particles per
species**, explicitly using the 2015 loader (§III). Tail quality validated to `|v| = 5 v_th,e`.

### Verdicts

| Planned statement | Verdict | Action |
|---|---|---|
| "A&M applied to multi-dimensional PIC loading" (P4-b) | **DROP as a gap claim; state as fact** | It is true, and it is *their* contribution, not a gap we fill. |
| "Existing Kappa loaders were not demonstrated in multidimensional PIC" | **DROP** | 2021 is 2D3V at 6.7×10⁷ particles/species. |
| "A&M's loader is one-dimensional in velocity space" | **DROP** | False for 2015. True only of 2014. |
| "A&M 2015 is a dedicated methods paper" | **MODIFY** | It is a Bernstein-wave PIC study whose §III B (≈2 pp.) introduces the loader. The dedicated methods paper is the 2014 CPC companion, which is 1-D only. |
| "A&M 2015 samples a bi-Kappa" | **MODIFY** | Isotropic 3-D Kappa. The anisotropic case is *used* in 2021 but written out nowhere. |
| "Prior loaders require a Gamma generator, hurting portability" | **PASS** | And stronger than expected: A&M 2015 requires `χ²_ν` with generally non-integer `ν`, and **never says how it is generated** — no appendix, no gamma-generator citation. Zenitani (2025) makes the same point. Scope it to *portability*, not correctness or efficiency. |
| "Zenitani et al. mischaracterize A&M 2015" | **DROP** | All seven characterizations verified accurate. |

### The differentiator this leaves

**The anisotropic bi-Kappa scale mixture is used but never written down.** A&M 2015 §III B
carries a general covariance matrix **R** in Eq. (17) and then sets **R** = σ²**I**; the 2018
thesis §2.3 is isotropic only; the 2021 paper cites 2015 without extension; ZUM2026 says
*"We can easily extend it for θ∥ ≠ θ⊥."* Everyone treats it as trivial; nobody documents it.

A documented, validated, field-aligned bi-Kappa implementation is defensible on this basis. A
novel *construction* is not, and the ZUM2026 "easily extend" sentence must be acknowledged in
the text rather than hidden.

## 7.2 Foundational solar-wind anisotropy — P2-a resolved, with two caveats

| File | Citation | Status |
|---|---|---|
| `HellingerTravnicekKasperLazarus2006_GRL33_L09101.pdf` | Hellinger, Trávníček, Kasper & Lazarus (2006), *GRL* **33**, L09101, doi:10.1029/2006GL025925 | Publisher PDF, fully verified |
| `Matteini2007_GRL34_L20105.pdf` | Matteini et al. (2007), *GRL* **34**, L20105, doi:10.1029/2007GL030920 | Publisher PDF, verified |
| `Pierrard2016_SolPhys291_2165_arXiv.pdf` | Pierrard, Lazar, Poedts, Štverák, Maksimovic & Trávníček (2016), *Solar Phys.* **291**(7), 2165–2179, doi:10.1007/s11207-016-0961-7 | ⚠ **arXiv preprint, not the published version** |
| `Bale2009_PRL103_211101_arXiv.pdf` | Bale, Kasper, Howes, Quataert, Salem & Sundkvist (2009), *PRL* **103**, 211101, doi:10.1103/PhysRevLett.103.211101 | ⚠ **arXiv preprint, not the PRL version** |
| `Marsch2018_AnnGeophys36_1607.pdf` | Marsch (2018), *Ann. Geophys.* **36**, 1607–1630, doi:10.5194/angeo-36-1607-2018 | Open-access version of record |

**Recommended minimal set: Hellinger 2006 + Matteini 2007 + Pierrard 2016**, with Marsch 1982
for the foundational observation at abstract level only and Bale 2009 optional.

**Hellinger et al. (2006)** — protons, Wind SWE + MFI, 1995–2001, ~1 AU. Supplies the
four-instability threshold formula and coefficients that everyone else reuses. Two scope facts
we must respect: (i) it is the **boundary of the occurrence distribution** that follows the
thresholds, not the distribution — their ¶10 says *"a majority of observations lies outside the
regions unstable"*; (ii) in the slow wind the operative bounds are the **oblique** modes
(mirror, oblique firehose), explicitly *"in contradiction with the results of the linear
theory"* which favours the parallel modes. If we name modes, name the oblique ones.

**Pierrard et al. (2016)** — **this is the anisotropy ⊗ Kappa tie-in.** The measured electron
VDF is fitted as a bi-Maxwellian core plus an **anisotropic bi-Kappa halo** with independent
`T_h,∥`, `T_h,⊥` and index κ. ~124 000 events, Helios 1 / Cluster II / Ulysses, 0.3–3.95 AU.
κ falls from 7.57 at 0.35 AU to 3.16 at 3.0 AU, and *"Deviations from isotropy decrease with
increasing κ"*. This is the single best citation for why an anisotropic Kappa with independent
θ⊥, θ∥ is the right object to sample.

### Two caveats that must not be papered over

1. **Pierrard 2016 is ELECTRONS.** A search specifically for a primary in-situ observation
   fitting an anisotropic bi-Kappa to solar-wind **proton** VDFs found none.
   `yoon2023bikappa`, already in our bib, is **theory**, not observation. The
   anisotropy⊗κ claim **for ions** is **BLOCKED**.
2. **Marsch et al. (1982) supports "anisotropic AND non-Maxwellian", but the non-Maxwellian
   character is beams, high-energy shoulders and heat flux — NOT a κ power-law tail.** Nothing
   in it licenses a Kappa claim for ions. Also: the publisher-deposited abstract carries
   typos, so **quote no numeric detail from it**; and the author order is Marsch, Mühlhäuser,
   Schwenn, Rosenbauer, Pilipp, Neubauer (Crossref + Hellinger's own reference list), not the
   order given by several web sources.

### Verdicts

| # | Planned statement | Verdict |
|---|---|---|
| S1 | "Temperature anisotropy is a persistent feature of the expanding solar wind" | **MODIFY** — Marsch's "persistent" is scoped to *high-speed streams*; low/intermediate-speed distributions show the opposite sign. Matteini supplies the radial ("expanding") half. |
| S2 | "the observed distribution of T⊥/T∥ is bounded by kinetic instability thresholds" | **MODIFY** — it is the *boundary of the occurrence distribution*; and name the oblique modes or none. |
| S3 | "ion distributions are commonly anisotropic **and** non-Maxwellian/suprathermal" | **MODIFY** — "non-Maxwellian (beams, shoulders, heat flux)" is supported for ions; a κ power-law tail is **not**. |
| S4 | "Simulations intended to represent such anisotropic suprathermal populations…" | **PASS**, once S3 is fixed upstream |
| S5 | any ion observation fitted with an anisotropic bi-Kappa | **BLOCKED** — no primary source found |
| S6 | any page/figure-level detail from Marsch et al. (1982) | **BLOCKED** — abstract only |

## 7.3 Net effect of round 2

1. **One planned claim inverted.** We were preparing to say prior loaders had not been shown in
   multidimensional PIC. The opposite is true and published. Removing this *before* a returning
   referee finds it is worth more than any claim it would have supported.
2. **The concession gets blunter again.** A&M 2015 Eq. (22) is the same scale mixture as
   Z&N 2022 and ZUM2026 Eq. (6). Three independent groups, 2015 / 2022 / 2026.
3. **The surviving differentiator sharpened** — see §7.1, "The differentiator this leaves".
4. **R1.7 is answerable** with Hellinger + Matteini + Pierrard, subject to the electron/ion
   caveat, which must appear in the text rather than being finessed.
5. **Two new blockers recorded** (ion anisotropy⊗κ; Marsch 1982 full text). Both are
   containable by narrowing the wording rather than by dropping the paragraph.
</content>
</invoke>
