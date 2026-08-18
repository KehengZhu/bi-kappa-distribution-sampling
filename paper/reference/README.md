# Primary sources for the JTJ1001 revision

Naming: `AuthorsYear_Journal_Volume_Page.pdf`. Every file below was identified by reading its
first page, not by trusting the download filename — two of them were mislabelled at the source.

> **The PDFs themselves are deliberately NOT tracked** — they are publisher copies, and the
> repository-wide `*.pdf` ignore rule covers them. This README is the tracked record of what
> was read and what it supports. Anyone reproducing the audit re-downloads from the DOIs below.

> **Round 2, 2026-08-17.** Five further sources obtained; see the two new tables. Verdicts and
> consequences are in `docs/revision/literature/step1_claim_audit.md` §7.

## Prior-art / related-work sources (sequencing step 1)

| File | Full citation | What it is needed for |
|---|---|---|
| `ZenitaniUsamiMatsukiyo2026_JGR131_e2025JA034669.pdf` | Zenitani, Usami & Matsukiyo (2026), *JGR Space Phys.* **131**, e2025JA034669, doi:10.1029/2025JA034669 | **The conflicting paper.** `(r,q)` → bi-Kappa at `r=0, q=κ+1`; Algorithm 3.1 = our core construction. Read in full. |
| `AbdulMace2014_CPC185_2383.pdf` | Abdul & Mace (2014), *Comput. Phys. Commun.* **185**(10), 2383–2386, doi:10.1016/j.cpc.2014.05.006 | Student-*t* route; first practical Kappa loader. Verify P4's characterization and whether 2014 or the 2015 companion covers multi-dimensional loading. |
| `ZenitaniNakano2022_PoP29_113904.pdf` | Zenitani & Nakano (2022), *Phys. Plasmas* **29**, 113904, doi:10.1063/5.0117628 | Beta-prime/Gamma-ratio route, relativistic Kappa. Confirm P4 attributes the machinery correctly. |
| `ZenitaniNakano2023_JGR128_arXiv.pdf` | Zenitani & Nakano (2023), *JGR Space Phys.* **128**(10), doi:10.1029/2023JA031983 | Loss-cone / KLC. ⚠ **This copy is arXiv:2309.06879v2, not the published version.** Fine for reading; cite the JGR version, and check any quoted detail against it. |
| `Zenitani2025_RNAAS.pdf` | Zenitani (2025), *Res. Notes AAS*, doi:10.3847/2515-5172/ae1c41 | Pareto-envelope rejection, acceptance ≈0.73–0.8 — the quantitative counter to "prohibitively low". ✅ **κ range verified from primary source:** the quoted efficiencies are *"0.806 at κ = 1.5, 0.785 at κ = 2, 0.750 at κ = 5"*, asymptotic to `√(πe)/4 ≈ 0.731`, under the recommended envelope index `n = κ/2`. Since the envelope requires `0 < n < κ − 1/2`, `n = κ/2` is **undefined for κ ≤ 1** — the range is not merely unquoted there, the recommended method does not apply. ✅ **§2 procedure transcribed and independently implemented** in `experiments/exp3_benchmark/`; our measured acceptance reproduces 0.806 / 0.786 / 0.751 at κ = 1.5 / 2 / 5. |
| `ZenitaniUmeda2026_EPS78_119.pdf` | **Zenitani & Umeda** (2026), *Earth, Planets and Space* **78**:119, doi:10.1186/s40623-026-02465-0 | Approximate Kappa generator, q-exponential inverse transform, GPU-oriented. ⚠ **Consensus listed this as "Zenitani et al."; the paper is Zenitani & Umeda.** Fix before it reaches `refs.bib`. |
| `An2022_JGR127_e2021JA030031.pdf` | An, Artemyev, Angelopoulos, Lu & Pritchett (2022), *JGR Space Phys.* **127**, e2021JA030031, doi:10.1029/2021JA030031 | Chebsampling. Kills the "inversion is intractable" premise. **Verify how broad a class it actually covers and in how many dimensions** — P3 currently says "broad classes of one- and low-dimensional distributions". |

## Regularized-Kappa sources (§1 fact 6)

The original Zenitani et al. cite is **Scherer, Fichtner & Lazar (2017), *Europhysics Letters*
**120**(5), 50002, doi:10.1209/0295-5075/120/50002** — EPL, not ApJ, which is why an ApJ search
for "Scherer 2017" came up empty. It is **still missing**; see the note below on whether we need it.

| File | Full citation | Relevance |
|---|---|---|
| `SchererLazarHusidicFichtner2019_ApJ880_118.pdf` | Scherer, Lazar, Husidic & Fichtner (2019), *ApJ* **880**, 118, doi:10.3847/1538-4357/ab1ea1 | *Moments of the Anisotropic Regularized κ-distributions.* **Arguably better for our purposes than the 2017 original**: this is the *anisotropic* regularized Kappa, i.e. the principled smooth-cutoff counterpart to our component-wise box cap on a bi-Kappa, and it gives the moments. |
| `SchererFichtnerFahrLazar2019_ApJ881_93.pdf` | Scherer, Fichtner, Fahr & Lazar (2019), *ApJ* **881**, 93, doi:10.3847/1538-4357/ab2df9 | *On the Applicability of κ-distributions.* Bears on the κ ≤ 3/2 moment-divergence discussion and on how to describe the normalization domain. |

## R1.7 foundational sources

| File | Full citation | Relevance |
|---|---|---|
| `Chandran2010_ApJ720_503.pdf` | Chandran, Li, Rogers, Quataert & Germaschewski (2010), *ApJ* **720**, 503–515, doi:10.1088/0004-637X/720/1/503 | Turbulence-driven perpendicular ion heating; the foundational reference R1.7 asks for in Introduction ¶2. |

## Abdul group — the prior-art chain (obtained 2026-08-17)

| File | Full citation | What it establishes |
|---|---|---|
| `AbdulMace2015_PoP22_102107.pdf` | Abdul & Mace (2015), *Phys. Plasmas* **22**(10), 102107, doi:10.1063/1.4933005 | **The multivariate-loader blocker, resolved.** Samples the **isotropic 3-D** Kappa via a genuinely trivariate normal-over-χ² scale mixture, Eq. (22), with a **single χ²_ν shared across all three components**. "One-dimensional" in the title is 1D3V *simulation geometry* (§VI: "Results of 1D3V particle-in-cell simulations"), not velocity-space dimensionality. ⚠ It is a Bernstein-wave PIC study whose §III B (≈2 pp.) contains the method — not a dedicated methods paper. ⚠ It never says how the non-integer-ν χ² deviate is generated. |
| `AbdulMatthewsMace2021_PoP28_062104.pdf` | Abdul, Matthews & Mace (2021), *Phys. Plasmas* **28**(6), 062104, doi:10.1063/5.0047638 | **Kills any "not demonstrated in multidimensional PIC" claim.** 2D3V GPU PIC, standard bi-Kappa (Summers & Thorne form), T⊥/T∥ = 3.0, 6.7×10⁷ particles/species, explicitly using the 2015 loader. |
| `Abdul2013_MScThesis_UKZN10413-12288.pdf` | Abdul (2013), MSc, UKZN, hdl:10413/12288 | Open access. Contains the 2014 CPC method (§2.1.2) and the Student-*t* ↔ Kappa equivalence derivation. |
| `Abdul2018_PhDThesis_UKZN10413-22458.pdf` | Abdul (2018), PhD, UKZN, hdl:10413/22458 | Open access. **§2.3 (pp. 23–28) is the fullest written statement of the loader anywhere** — and it is *isotropic only*, which is the point. |

**The load-bearing negative result:** the anisotropic bi-Kappa scale mixture is *used* by this
group (2021) and by Zenitani (ZUM2026 Algorithm 3.1, *"We can easily extend it for θ∥ ≠ θ⊥"*),
but is **written out by nobody**. That, not the construction, is what we can document.

## R1.7 foundational sources — round 2 (obtained 2026-08-17)

| File | Full citation | Exact proposition it supports |
|---|---|---|
| `HellingerTravnicekKasperLazarus2006_GRL33_L09101.pdf` | Hellinger, Trávníček, Kasper & Lazarus (2006), *GRL* **33**, L09101, doi:10.1029/2006GL025925 | Protons, Wind SWE+MFI, 1995–2001, ~1 AU. Supplies the four-instability threshold formula + coefficients (Table 1). ⚠ It is the **boundary of the occurrence distribution** that follows the thresholds — ¶10: *"a majority of observations lies outside the regions unstable"*. ⚠ In the slow wind the operative modes are the **oblique** ones, explicitly contradicting linear theory's preference for the parallel ones. |
| `Matteini2007_GRL34_L20105.pdf` | Matteini et al. (2007), *GRL* **34**, L20105, doi:10.1029/2007GL030920 | Helios + Ulysses, 0.3–2.5 AU. **The source for the word "expanding":** fast wind moves from low-β, T⊥ > T∥ at 0.3 AU to β∥ ≈ 1, T∥ > T⊥ at 1 AU. |
| `Pierrard2016_SolPhys291_2165_arXiv.pdf` | Pierrard, Lazar, Poedts, Štverák, Maksimovic & Trávníček (2016), *Solar Phys.* **291**(7), 2165–2179, doi:10.1007/s11207-016-0961-7 | **The anisotropy ⊗ Kappa tie-in.** Electron VDF fitted as bi-Maxwellian core + **anisotropic bi-Kappa halo** with independent T∥, T⊥ and index κ. ~124 000 events, 0.3–3.95 AU; κ 7.57 → 3.16; *"Deviations from isotropy decrease with increasing κ"*. ⚠ **ELECTRONS, not ions.** ⚠ On-disk copy is the arXiv preprint (titled "… I. Comparison"); the published title drops the "I." |
| `Bale2009_PRL103_211101_arXiv.pdf` | Bale, Kasper, Howes, Quataert, Salem & Sundkvist (2009), *PRL* **103**, 211101, doi:10.1103/PhysRevLett.103.211101 | Optional. Gyroscale fluctuation power *enhanced along* the mirror/IC and oblique-firehose thresholds — evidence the instabilities are **actively regulating**, not merely coinciding. Says nothing about non-Maxwellian tails. ⚠ arXiv preprint. |
| `Marsch2018_AnnGeophys36_1607.pdf` | Marsch (2018), *Ann. Geophys.* **36**, 1607–1630, doi:10.5194/angeo-36-1607-2018 | Open-access fallback by the author of the 1982 paper, reproducing the 1982 Helios VDFs. Fig. 12b: measured VDF at 0.35 AU with *"anisotropic core and extended tail along the magnetic field"*. It is a review, so second choice as a foundational citation, but verifiable end to end. |

## Still outstanding

- **Marsch et al. (1982), *JGR* 87(A1), 52–72, doi:10.1029/JA087iA01p00052** — full text **not
  obtained** (not open access; Wiley Cloudflare-gated). The publisher-deposited abstract was
  verified via Crossref and supports abstract-level propositions only: *"A marked anisotropy in
  the core of proton distributions … is a persistent feature of high speed streams and becomes
  most pronounced in the perihelion (≈0.3 AU)"*. ⚠ **Quote no numeric detail** — the deposited
  abstract carries typos. ⚠ Author order is Marsch, **Mühlhäuser, Schwenn, Rosenbauer**, Pilipp,
  Neubauer (Crossref + Hellinger's reference list), not the order several web sources give.
- **A primary ion observation fitted with an anisotropic bi-Kappa** — **none found.** Pierrard
  2016 covers electrons. `yoon2023bikappa` in our bib is *theory*, not observation. Do not
  generalize the electron result to ions.
- **Published versions of Bale 2009 and Pierrard 2016** — on-disk copies are arXiv preprints.
  Verify before any quoted detail reaches the manuscript.
- **Scherer, Fichtner & Lazar (2017), EPL 120, 50002** — optional. Needed only if we want to cite
  the regularized Kappa at its origin; the 2019 anisotropic paper may serve better on its own.
- **Small-shape Gamma prior art**, now relevant because Experiment 4 promoted C5: Ahrens & Dieter
  (1974) *Computing* **12**(3) 223–246; Best (1983) *Computing* **30**(2) 185–188; Devroye (1986)
  ch. VII §2.6. Needed only if the log-domain construction is landed.
