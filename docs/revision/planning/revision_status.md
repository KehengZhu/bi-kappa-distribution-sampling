# JTJ1001 — revision status after the manuscript rewrite

**Date:** 2026-08-18. **State:** manuscript rewritten, cap-geometry section corrected and
compressed, building clean; rebuttal letter not yet written, per plan §14.

**Publication-hygiene pass (2026-08-18, later same day).** The pre-writing gate was added to
`CLAUDE.md` and to plan §*Mandatory pre-writing publication-hygiene gate*, then run against the
manuscript. Prose changes, no evidence or figures touched:

- **Abstract rewritten a third time** under the gate. Full account in plan §1, status note dated
  2026-08-18. Both operating-range numbers retained; λ now given units ("the cap set at 50
  thermal speeds") because the λ^−(2κ−1) decay law moved out of the Abstract.
- **§I contribution item 5** no longer says the performance measurement was "used to retire the
  performance claims of the previous version of this work"; it states what the measurement shows.
- **§VIII opening** no longer opens on "the previous version of this work claimed…"; it states
  the two properties commonly assumed of the construction and that neither holds. A first-time
  reader had no referent for the old wording.
- **§IX** now says "earlier releases of the package" rather than "the previous version of this
  work" for the cap-default flip, which is a software changelog fact and belongs in §IX.
- **Terminology:** "rejected mass" → **"rejected fraction"** throughout `main.tex` (7 sites,
  including the §IV.C subsection title, which had contradicted its own first sentence);
  "second-moment anisotropy" → "second-moment ratio" (1 site).
- **Left alone, flagged for a decision:** §X *Prior art and independent development* still states
  that the prior publications "were identified after this work was first submitted". Gate rule 4
  argues against it; publication ethics argues for it. Not an editorial call to make silently.

**Introduction rewrite (2026-08-18, same pass).** §I rewritten against
`manuscript_revision_plan.md` §2 and the matrix's R2.A1/A2/A3 and §9.3 ledger. Decisions
recorded in plan §2's dated status block; the prose changes are:

- The difficulty paragraph no longer opens on "the earlier literature, **including our own
  previous framing**", and no longer says the inverse-CDF point is what "earlier framings
  treated as decisive". Both were unreadable without the revision history.
- The prior-art paragraph no longer opens "and the chronology matters", and the standalone
  independent-development disclaimer is deleted from §I. §X carries that disclosure once.
- Four difficulties became three; the heavy-tail-diagnostic failure is now subordinate to
  moment non-existence, per plan §0b rule B.
- The inverse-CDF sentence moved into the prior-art paragraph and is stated positively, so it
  can no longer be read as evidence of a sampling gap.
- The contribution list lost the "in the following order of weight" framing (the §9.3 tier
  ledger showing through), the `1.35×10⁷` draw count, and item 3's RNG/default/workflow
  bookkeeping. Order is unchanged and still matches plan §0a.
- The closing line no longer defends the package ("not offered as a contribution independent
  of them"); it states the three conditions a reliable loader must satisfy, matching the
  Abstract's closing.
- **New paragraph answering R2.A3**, which the matrix still lists HIGH/UNADDRESSED: anisotropy
  adds no new stochastic problem, only implementation and verification burden. Forwards to
  §II.E. R2.A3 can now be marked addressed in the matrix.
- **Terminology:** Summary's "the scale-implied anisotropy" → "the scale-implied second-moment
  ratio", closing the last seam left by the earlier normalization.

Build after the pass: 18 pages, 0 undefined references, 0 undefined citations.

**Build:** `paper/overleaf/latex-build/main.pdf`, 18 pages, **0 undefined references,
0 undefined citations, 0 overfull boxes**. Artifacts confined to `latex-build/`.
(16 underfull hboxes remain, all cosmetic looseness in two-column text; three are the
unbreakable GitHub/DOI URLs in §IX.)

**Provenance:** all four experiments re-verified by `make verify` against committed
`checksums.sha256` — every file matches. The released header
`cpp/bi_kappa_distribution.H` hashes to `6b138af5…`, which is the `sampler_header_sha256`
recorded in all four `results/exp*_results.json`. **No experiment was rerun**; none needed to be.

**Tests:** `cpp/` suite builds warning-free under `-Wall -Wextra -std=c++11` and all 11 groups
pass (`test_transform`, `test_define_*` ×5, `test_max_normalized_velocity_*` ×2, `test_seed_api`,
`test_radius_formation`, `test_no_cap_semantics` incl. B5–B8), exit 0.

**Analytic provenance:** `paper/figures/verify_cap_geometry.py` — 34 deterministic quadrature
assertions covering every closed form in §IV.F, no sampling and no experiment input, exit 0. It
runs from `make_manuscript_assets.py` before any asset is written, so asset generation fails if a
closed form and independent quadrature ever disagree.

---

## 0. What changed on 2026-08-18

An audit of the cap-geometry work found two boundary errors in
`../literature/cap_geometry_novelty_audit.md` and a strategic overreach in the §IV.F that had been
written from it. Both are fixed; details in that file's new §0.

**Mathematical corrections.**

1. **The speed-cap plateau lives on the open interval `1/2 < κ < 3/2`, not `κ ≤ 3/2`.** At
   `κ = 3/2` the truncated radial second moment diverges only logarithmically, and the leading
   `log C` coefficient is direction-independent, so the limit is exactly **1** and the bias decays
   as `O(1/log v_max)` — slowly, but to zero. The manuscript said "for κ ≤ 3/2 it does not
   [vanish]" in two places (speed cap, cylinder). Both fixed.
2. **`0.645` is a wide-cap asymptote, not the value at every cap width.** The finite-cap ratio at
   `κ = 0.75`, `θ∥/θ⊥ = 2` runs 0.498 → 0.607 → 0.643 → 0.6446 as `v_max/θ⊥` runs 3 → 10 → 10² →
   10⁴, approaching from **below**: a narrow cap is worse, not better. "No matter how wide the
   bound is made" is gone.
3. **Terminology.** On `1/2 < κ < 3/2` there is no untruncated `T∥/T⊥` for a cap to be biased
   against, so nothing is described as biasing a "true" or "imposed" temperature anisotropy. The
   manuscript compares the **second-moment ratio of the capped law** against the **scale-implied**
   `(θ∥/θ⊥)²`. Both capped moments diverge; only their ratio has a finite limit, and that limit is
   a property of the chosen regularization. Two geometries can select two different limits while
   both removing vanishing probability mass — which is the point, not a paradox.

**Strategic corrections to §IV.F.** It had grown into a small cap-geometry study, which conflicts
with the minimal-closure decision and opens reviewer surface (other anisotropy ratios? cylinder
closed form? matched ε? simulation consequence?). Compressed back to an analytic caution plus a
property map:

- **Deleted** the 5 × 3 numerical sweep over (κ, θ∥/θ⊥). It read as "we systematically studied the
  physical-speed cap parameter space."
- **Added** Table III, a four-row conceptual property table (box / normalized ellipsoid / physical
  sphere / gyrotropic cylinder). Its job is the conceptual distinction that is actually worth
  having: **gyrotropy of `f` ≠ angular-law preservation ≠ second-moment-ratio preservation**.
- **Kept** one equation (the wide-cap limit) and one illustrative number (`\SpeedCapLimit` = 0.645).
- **Replaced** the global "a code … should prefer it" with a conditional recommendation keyed to
  what a simulation requires — consistent with this repo's own "no universal best cap" finding.
- **Not** promoted to the Abstract headline or the numbered contribution list. It supports the
  existing thesis; it is not a new pillar.

**Thesis unchanged:** *reliable bi-Kappa loading requires explicit control of the target law, any
finite-velocity bounding policy, and the numerical implementation.* The cap results are three
supporting examples of "bounding policy matters", not a study of optimal bounding geometry.

**Deliberately still not done:** no Exp 5, no sampler-geometry change, no PIC demonstration, no
Monte-Carlo run for the plateau (unnecessary once the manuscript quotes an analytic limit rather
than a measured value).

---

## 1. Manuscript structure, old → new

| Old | New | Change |
|---|---|---|
| Abstract | Abstract | Rewritten from scratch |
| I Introduction | §I Introduction | Substantially rewritten; prior-art chronology + 5-item contribution list |
| — | §II Background and related work | **New** — answers R2.A1, R2.A3, R2.B |
| II Notations and probabilistic background | §II.B, §III.A–B | Split: physical content → Background, mathematical machinery → Construction |
| III Sampling method (A–D) | §III Probabilistic construction | Rewritten; full Jacobian, geometry, pseudocode |
| III.E Truncation | §IV The optional component-wise cap | **Rewritten as a characterization**; 6 subsections |
| IV Sampling and verification | §VI Validation | Restructured around the test hierarchy |
| — | §VII Numerical considerations and operating range | **New** — R1.3 low-κ, Exp 4 |
| — | §VIII Performance | **New** — R1.4, and it reports against us |
| V Magnetic field-aligned frame | §V Magnetic field-aligned frame transform | Retained, validation + edge case added |
| VI Code availability | §IX Software availability | Rewritten and corrected |
| — | §X Discussion and limitations | **New** |
| VII Summary | §XI Summary | Rewritten |

---

## 2. Reviewer coverage — where each comment is answered

### First referee

| # | Comment | Answered at | Evidence |
|---|---|---|---|
| **R1.1** | Separate exact untruncated from cutoff sampling | §III.F (scoping of "distributionally exact", no outer rejection loop); §II.C (three target laws); §IV.A (Eq. 20 capped law); §IX (API corrected, default flipped) | Exp 2 subsequence check 240/240 |
| **R1.2** | Disclose truncation settings and validation provenance | §VI.A (mode, κ, θ, N, seeds, RNG, precision for every dataset); every figure/table caption states mode; §IV.F (forward reference deleted, replaced) | Exp 2 manifest |
| **R1.3** | Expand validation to radial law, low κ, one truncated case | §VI.B (radial, Table IV); §VI.C (direction, independence, MAD); §VI.D (frame); §VI.G (bounded scope statement); §IV (capped example w/ λ and rejected fraction, Table II); §VII (operating range) | Exp 1, 2, 4 |
| **R1.4** | Support or moderate speed and constant-cost claims | §VIII entire, Table VII. All speed claims **deleted**; measured result is unfavorable and reported as such | Exp 3 |
| **R1.5** | Discuss the larger κ=2 variance errors | §VI.F — the two largest stated explicitly (−3.6% $v_x$, −6.7% $v_\perp$), replicate spread reported, finite-sample downward bias explained, **and** the fourth-moment caveat (spread is not a standard error at κ=2) | Exp 1, regenerated Table V |
| **R1.6** | Make Figs. 2–4 self-contained and quantitatively labeled | Figs. 4–5 captions carry θ⊥, θ∥, N, κ, mode, seeds, precision, curve identities; ordinate "Probability density"; axes normalized and stated dimensionless; histogram normalization stated; bi-Maxwellian-limit purpose stated; **clipping disclosed and made consistent** | regenerated by `paper/figures/make_manuscript_assets.py` |
| **R1.7** | Add foundational, sampling, and PDF references/definitions | §II.A (PDF definition: non-negativity, unit integral, Lebesgue measure, variables named); §I ¶2 (Hellinger 2006, Matteini 2007, Chandran 2010, Pierrard 2016); §I ¶3 (Birdsall, Hockney, Cartwright, Verboncoeur); §II.D (Devroye, and acceptance quantified rather than asserted) | 8 new primary-source-verified bib entries |

### Second referee

| # | Comment | Answered at |
|---|---|---|
| **R2.A1** | State of Kappa-sampling literature | §I ¶5 (chronology, primary-source accurate); §II.D (survey by method class); §X ¶1 |
| **R2.A2** | Difficulty, solved problems, novelty, performance evidence | §I ¶4 (four real difficulties, inverse-CDF demoted and stated precisely); §I ¶6 (explicit non-novelty statement); §I ¶7 (5-item contribution list); §VIII (measured: we do **not** outperform) |
| **R2.A3** | Isotropic vs bi-Kappa challenges | §II.E — affine rescaling ⇒ same radial law; added burden is parameterization, frame, cutoff geometry, joint validation |
| **R2.B** | Background, methods, parameter meanings | §II.A–D; Table I gives meaning, domain **and units**; thermal-speed vs temperature convention stated explicitly; θ never called a temperature |
| **R2.C1** | Purpose of the projected 1-D sample | §VI.E ¶3 — computed *from* the final 3-D vectors, never an input; states what it can and cannot establish |
| **R2.C2** | Role of the transformation equation | §III.A — explicitly a derivational density rule, **not** an executable step; applied three times (Eqs. 12, 14, 15) |
| **R2.C3** | Isotropization derivation and ellipsoid/sphere language | §III.C — `V = A U`, quadratic-form collapse (Eq. 11), `|det A|` (Eq. 12), cancellation shown (Eq. 13); geometry paragraph clarifies **ellipsoids in V-space, spheres in U-space** by naming the coordinates rather than contradicting the referee; §V states $Q$ orthogonal ⇒ preserves the radial law |

---

## 3. Figures and tables

| Item | Content | Source | Answers |
|---|---|---|---|
| Fig. 1 | Sphere-sampling bias (retained) | existing | — |
| Fig. 2 | Pseudocode for one uncapped draw | — | R2.C2 |
| Fig. 3 | Cap: TV vs λ; p99.9 ratio; azimuthal $|z_4|$ | Exp 2 | R1.1, R1.2 |
| Fig. 4 | Cartesian marginals vs exact marginal and bi-Maxwellian limit | Exp 1 raw | R1.6 |
| Fig. 5 | Q–Q vs exact marginal (scaled Student-$t$) | Exp 1 raw | R1.3, R1.6 |
| Table I | Symbols, meanings, domains, **units** | — | R2.B |
| Table II | TV / p99.9 ratio, κ × λ | Exp 2 | R1.1, R1.2 |
| Table III | **New** — what four bounding regions preserve (conceptual, no sweep) | analytic | R1.1 |
| Table IV | Validation summary, κ = 0.51…10 | Exp 1 | R1.3 |
| Table V | Second moments, κ > 3/2 only, with replicate spread | Exp 1 raw | R1.5 |
| Table VI | Non-finite fraction by precision | Exp 4 | R1.3 |
| Table VII | Per-sample cost, three implementations | Exp 3 | R1.4 |

**Tables II and IV–VII are `\input` fragments generated by
`paper/figures/make_manuscript_assets.py`.** No number is typed into the manuscript by hand.
Table III carries no measured numbers; the one computed number quoted in §IV.F prose comes from
`tables/capgeom.tex`, which defines `\SpeedCapLimit` and is `\input` in the preamble — so even
that number is generated, not typed. Rerun with
`uv run --project ../../python python make_manuscript_assets.py` from `paper/figures/`; the
cap-geometry self-test runs first and aborts the run on any disagreement.

**Table numbering changed** on 2026-08-18: the new §IV.F table takes III, so the old III–VI are now
IV–VII. Any section/table cross-reference in the response letter must be read off the *current*
build.

**Figures dropped:** the three separate `marginal-normal-comparison-k=*.pdf` (merged into one
3×3 panel) and the two `marginal-qqplot-k=*.png` (merged into one 2×3 panel, and now vector).
The KDE curve was removed: standard bandwidth rules are calibrated to a normal reference
through the sample standard deviation, which does not exist for κ ≤ 3/2 and is not stable at
κ = 2. Histogram + exact marginal + bi-Maxwellian reference is the defensible triple.

---

## 4. Mechanical audits — all passed

*(Re-run 2026-08-18 after the §IV.F rewrite; all six still pass. Counts below re-verified.)*

1. **Terminology.** `rejection-free` 0; `constant time` 0; `prohibit*` 0; `outperform` 0.
   `exact` survives only as mathematical exactness or as "distributionally exact for the
   uncapped target in exact arithmetic". `novel`/`new` survive only in negations. `fast`
   survives only where a *competitor* is faster, plus "fast wind" (physical). `negligible`
   survives once, describing the cap's *timing* overhead. `confirm` survives once, meaning
   "confirms the implementation, not the fact".
2. **Mode.** Every figure caption, every table caption, and every numerical statement declares
   uncapped or capped.
3. **Moments.** No second-moment or variance interpretation of the *untruncated* law at κ ≤ 3/2
   anywhere; Table V starts at κ = 2; §IV.D explicitly warns against the capped variance at
   κ ≤ 3/2. **§IV.F is the one place that computes second moments in that range**, and it is
   allowed because the *capped* law has bounded support and therefore always has them — the
   section says so in as many words, compares against the scale-implied `(θ∥/θ⊥)²` rather than any
   temperature, and states that no untruncated `T∥/T⊥` exists to be biased against.
3a. **Gyrotropy, two-level rule.** No unqualified "does not preserve gyrotropy" (0 occurrences);
   no "agyrotropic" (0); every one of the five sites states both halves — continuous gyrotropy of
   `f` broken, pressure-tensor gyrotropy preserved exactly. See §4b.
3b. **Plateau boundary.** No `κ ≤ 3/2` in any plateau or speed-cap context; the two that existed
   (speed cap, cylinder) are now `1/2 < κ < 3/2`, with the `κ = 3/2` logarithmic case stated
   separately. No "no matter how wide" (0). No global "should prefer" (0).
4. **Citations.** A&M 2014 described as one-dimensional (Bailey polar transform); "multivariate"
   appears once, correctly, for the Student-$t$ class; Zenitani & **Umeda** authorship correct in
   `refs.bib`; ZUM 2026 Algorithms 3.1/3.2 described as writing the anisotropic loader out
   explicitly; §I states multidimensional PIC use of prior loaders **is** established.
5. **Empirical range.** All three occurrences use "Within the implementations and test protocol
   examined here, reliable operation was validated down to…".
6. **Compilation.** 0 undefined refs, 0 undefined citations, 0 overfull boxes.

---

## 4a. Novelty audit of the replacement claims — 2026-08-17

Run **after** the rewrite, because retiring a false claim is only half the job: whatever replaces
it is a new claim with the same burden of proof. Full record: `../literature/step1_claim_audit.md`
§8; standing procedure: `manuscript_revision_plan.md` §13 pass 6.

| Claim (plan §0a) | Verdict |
|---|---|
| 1. Direct 3-D validation incl. `1/2 < κ ≤ 3/2` | clear |
| 2a. `TV = rejected mass`, `λ^−(2κ−1)` | clear — framed as characterization, not discovery |
| 2b. Cap does not preserve gyrotropy | **superseded 2026-08-18** — see below |
| 3. Arbitrary-**B**-frame loading | clear **as scoped** — survives only because §V concedes it is elementary linear algebra |
| 4. Finite-precision characterization | ⚠ **one defect, fixed** |
| 5. Performance characterization | clear |

**The defect.** §X claimed the log-domain small-shape Gamma mitigation with **no citation**,
reading as though the idea were ours. It is an established remedy (Liu, Martin & Syring 2017;
Devroye 1986). Fixed in `main.tex` §X, which now cites both and frames our step as an adaptation.
This is the retired Gamma-ratio failure mode one level down, and it is why pass 6 is now standing.

**Two accuracy corrections**, neither a novelty defect: regularized-Kappa samplers *do* exist
(ZUM 2026 gives two), so §II C and §X no longer imply a gap; and Scherer et al. 2020 (the
κ-cookbook) is now cited at §II C and §X as the systematic survey of κ-variants.

**Bibliography:** `liu2017gamma` and `scherer2020cookbook` added — authors/title/journal/year/DOI
verified, **volume and pages deliberately left blank** rather than guessed. See §5.

## 4b. Claim 2b superseded — 2026-08-18

The 2026-08-17 audit called "the cap does not preserve gyrotropy" the **strongest genuinely new
result**. That verdict rested on an incomplete statement of the result, and the correction cuts
both ways.

**Weaker than believed, in two respects.** The damage is `O(ε)` in the rejected mass — any bounded
directional statistic obeys `|Δ| ≤ 2Gε/(1−ε)`, and `a₄ ∝ ε` empirically — so unlike the tail
quantiles it is *already bounded by a number the paper reports*. And it is invisible to the
diagnostic a referee would reach for: the cube's permutation and sign symmetries make the capped
pressure tensor **exactly** gyrotropic, with `E[v∥²]/E[vx²] = (θ∥/θ⊥)²` exactly, for every λ and
every κ > 1/2. Left unqualified, the old wording invited a referee to compute the pressure tensor,
find it perfectly gyrotropic, and conclude we overstated.

**Stronger than believed, as stated correctly.** The result is not "the cap breaks gyrotropy" but
*the cap breaks continuous gyrotropy of `f` while preserving pressure-tensor gyrotropy and the
scale-implied second-moment ratio exactly* — i.e. **second-order diagnostics cannot detect this
artifact at all**. That is a more useful thing to tell a code author than the original claim, and
it is harder to attack because both halves are proved rather than measured.

All four stale sites are now fixed (Abstract, §I contribution list, §X Discussion, §XI Summary),
and `manuscript_revision_plan.md` rule **D** now forbids the unqualified form so it cannot
regenerate. The `a₄` statistic is **not** to be called a standard agyrotropy measure; it is a
higher-order angular diagnostic, and every standard pressure-tensor agyrotropy measure returns
exactly zero here.

**New claim introduced by §IV.F, and its scope.** The physical-speed-cap plateau is stated as a
derivation, with no novelty claim attached, and deliberately kept out of the Abstract and the
contribution list. Its defensible narrow form is: *when the corresponding untruncated moment does
not exist, a vanishing total-variation distortion does not determine the limiting value of a
cutoff-dependent moment ratio.* Do **not** generalize this to "small rejected mass ≠ faithful
moment-derived structure" — for bounded statistics and convergent moments the rejected mass *does*
control the error, and the broad version would be indefensible.

## 5. Remaining manuscript-only cleanup

- **Point-by-point response letter and list of changes.** Deliberately not started; plan §14
  requires the manuscript to be stable first. It now is, so this is the next task.
- **Marked-up PDF** for resubmission (editor requires one).
- **Zenodo DOI/version.** `bikappa_code` cites the concept DOI. Before submission, mint a
  release matching the tested tree and decide whether to cite the version DOI instead — §IX
  currently says "concept DOI resolving to the latest archived version", which is accurate but
  weaker than citing the exact tested release.
- **`omelchenko2023agu`** is an AGU abstract; check the journal's tolerance for it.
- **Volume/page numbers for `liu2017gamma` and `scherer2020cookbook`** — added by the novelty
  audit with DOI verified but volume/pages left blank on purpose. Complete before submission;
  this is the only known incomplete bibliographic data in the file.
- **Section/equation cross-references in the response letter** must be taken from the *final*
  numbering, which changed again on 2026-08-18: **Eqs. 1–29, Figs. 1–5, Tables I–VII**. The new
  §IV.F table is III, so the old Tables III–VI are now IV–VII. Do not copy numbers from the
  2026-08-17 draft of this file.
- **The letter must state the gyrotropy result in its two-level form** (§4b). Reusing the
  2026-08-17 shorthand from `reviewer_response_matrix.md` would put an attackable claim in front of
  the same referees; those three sites in that file are now corrected, but check before quoting.
- **APS context.** Per the Chief Editor (2026-08-04), the revised materials return to the *original*
  referees, with the 2026-07-30 prior-art-overlap explanation forwarded alongside. Per Y. Omelchenko
  the same day, the simulation application is a separate paper and this one is method + software.
  Both facts argue against letting the cap-geometry material expand any further.

## 6. Deliberately not done

- Landing the log-domain Gamma mitigation (authors' call; breaks the RNG stream). Disclosed in §X.
- Making `bi_maxwellian_distribution` symmetric with the new bi-Kappa default.
- Any new simulation application.
- Any claim resting on the retired items in `reviewer_response_matrix.md` §9.3.
- **Exp 5 / any cap-geometry experiment.** Not needed: §IV.F quotes an analytic limit, so there is
  no measured value requiring Monte-Carlo support. `verify_cap_geometry.py` supplies the provenance
  deterministically.
- **Implementing any alternative bounding geometry.** §IV.F says so explicitly and leaves the
  matched comparison to future work.
