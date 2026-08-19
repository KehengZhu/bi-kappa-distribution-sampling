# Cap-geometry novelty, prior-art and feasibility audit

**Date:** 2026-08-17. **Status:** audit only. No cap implemented, no manuscript claim added,
no API changed, no cap-comparison experiment run.

**Method.** Repository state recovered first (`revision_status.md`, `step1_claim_audit.md`,
`reviewer_response_matrix.md`, Experiment 2, `main.tex` §IV, and the shipped predicate
`bi_kappa_distribution.H:withinNormalizedVelocityCap`). All geometry mathematics below was
derived independently against the repository's own convention and then **verified numerically**
by short audit probes (§3.6). Prior art was checked against the 19 primary PDFs already in
`paper/reference/` plus targeted literature search.

> **Verification caveat, stated up front.** Six parallel literature agents were launched and all
> six died on a shared upstream API rate limit before reporting. The prior-art sections below
> rest on (a) the local primary PDFs, which I read directly, and (b) a reduced set of targeted
> searches I ran myself. Several statistics citations in §5 are named from the standard
> literature but **were not opened as primary sources**. They are tagged accordingly and listed
> in §15. The mathematics in §3–§4 does not depend on any of them.

---

# 0. Corrections applied 2026-08-18

Two of this audit's statements about the speed-cap plateau were wrong at the boundary and were
corrected before anything went into the manuscript. Both corrections are now enforced by
`paper/figures/verify_cap_geometry.py`, a deterministic quadrature self-test (34 assertions, no
sampling) that runs from `make_manuscript_assets.py`, so asset generation fails if either drifts.

**C1. The plateau exists on the open interval $1/2<\kappa<3/2$, not on $\kappa\le3/2$.**
At $\kappa=3/2$ exactly the truncated radial second moment diverges *logarithmically*,
$\int^{R_{\max}}r^4(1+r^2)^{-5/2}dr\sim\log R_{\max}$. Writing $R_{\max}(\mu)=C/h(\mu)$ gives
$\log R_{\max}=\log C-\log h(\mu)$: the direction-dependent piece is an $O(1)$ correction to a
direction-*independent* leading term, so the leading coefficients cancel from the ratio and the
limit is exactly $1$. The bias does not plateau there — it decays as $O(1/\log v_{\max})$.
Verified: at $\theta_\parallel/\theta_\perp=2$ the ratio runs $0.9644$, $0.9772$, $0.9833$,
$0.9891$ for $v_{\max}/\theta_\perp=10^{4},10^{6},10^{8},10^{12}$, with
$(1-r)\log v_{\max}$ constant to $\pm4\%$ over that range. Consistently, the closed form of §3.5
returns exactly $1$ at $\kappa=3/2$ (the exponent $(3-2\kappa)/2$ vanishes, so $W\equiv1$), and
the plateau closes continuously from below: $0.9434$, $0.9713$, $0.9942$, $0.99942$ at
$\kappa=1.4,1.45,1.49,1.4999$. There is no discontinuity at the endpoint — which is why the
correct statement is a strict inequality and not a change of regime.

**C2. $0.6446$ is a wide-cap asymptote, not the value at every cap width.**
This audit's §1 said "at *every* cap width" and §3.5 said "independent of $v_{\max}$"; both
overstate. The *limit* is independent of $v_{\max}$; the finite-cap ratio is not. Deterministic
quadrature at $\kappa=0.75$, $\theta_\parallel/\theta_\perp=2$ gives $0.4981$, $0.6070$,
$0.6429$, $0.64462$ at $v_{\max}/\theta_\perp=3,10,10^{2},10^{4}$ — approaching $0.644623$ from
*below*, with the residual falling as $v_{\max}^{-(3-2\kappa)}$ (the exponent is confirmed to
within 25% across three decades for every $\kappa$ tested). A narrow cap is worse than the
asymptote, not better. The correct manuscript phrasing is "approaches $0.645$, rather than
unity, as $v_{\max}\to\infty$."

**C3. Terminology.** On $1/2<\kappa<3/2$ the untruncated law has no second moment, so there is no
"true" or "imposed" $T_\parallel/T_\perp$ for a cap to be biased against. This audit's use of
"realized temperature anisotropy" is therefore unsafe. Read every such phrase below as *the
second-moment ratio of the capped law*, compared against the *scale-implied*
$(\theta_\parallel/\theta_\perp)^2$. Both capped moments diverge as $v_{\max}\to\infty$; only the
ratio has a finite limit, and that limit is a property of the chosen regularization. This is not
a paradox — it is the cleanest statement of the result, and the manuscript now says it that way.

The $\le$ vs $<$ occurrences below have been corrected in place. Occurrences of $\kappa\le3/2$
that refer to *non-existence of the untruncated second moment* (§2, §4.3, §4.5) are correct as
written and were left alone.

---

# 1. Executive verdict

**Is cap-geometry comparison scientifically worthwhile?** Yes, but not for the reason the
handoff assumed, and not at the size it assumed.

**Is anything plausibly new?** One thing, and it is *not* the box/gyrotropy result that the
prior audit flagged as strongest. The strongest surviving candidate is:

> A hard cap on **physical speed** — the geometrically obvious and simulation-motivated choice —
> corrupts the second-moment ratio of the capped law relative to the scale-implied
> $(\theta_\parallel/\theta_\perp)^2$, and for $1/2<\kappa<3/2$ it does so by a **finite amount that
> does not vanish as the cap is widened**. The bias tends to a computable constant (§3.5, verified
> by deterministic quadrature). At $\kappa=0.75$, $\theta_\parallel/\theta_\perp=2$, the capped
> $E[v_\parallel^2]/E[v_x^2]$ approaches $0.6446$ times the scale-implied value as
> $v_{\max}\to\infty$ — see correction C2 in §0: this is an asymptote, not the value at every
> cap width. At $\kappa=3/2$ exactly the limit is instead $1$, reached as $O(1/\log v_{\max})$ (C1).

This is a new *hazard class*. The paper's existing thesis is "a small rejected mass does not
imply a faithful tail." This extends it to "…and for a speed cap it does not imply a faithful
*second-moment ratio* either, however wide the cap is made, over most of the $\kappa$ range Kappa
distributions are used for." That is the same argument one level deeper, not a second paper.

The precise form, which is narrow enough to be hard to attack: **when the corresponding
untruncated moment does not exist, a vanishing total-variation distortion does not determine the
limiting value of a cutoff-dependent moment ratio.** Note what this does *not* say — it makes no
claim about bounded statistics or convergent moments, where a small rejected mass does control the
error (§4.5). Stating it generally ("small rejected mass ≠ faithful moment-derived structure")
would be too broad to defend.

**Three findings that change the current manuscript regardless of what is added.**

1. **§IV.E's gyrotropy claim is imprecise in a way a referee can attack.** The box cap preserves
   the **pressure tensor's** gyrotropy *exactly* (§3.4, §4.2): $P_{xx}=P_{yy}$, all off-diagonals
   vanish, and $T_\parallel/T_\perp$ is exactly right. Every standard agyrotropy measure
   (Scudder–Daughton, Aunai, Swisdak) returns **zero**. The breaking is real but lives at fourth
   order in the distribution function. As written, "the cap does not preserve gyrotropy" invites
   a referee to check the pressure tensor, find it perfectly gyrotropic, and conclude we
   overstated. **The correct — and stronger — statement is that the distortion is invisible to
   the standard diagnostic.**
2. **§IV.E overstates detectability against the paper's own table.** The text says the modulation
   "at small $\kappa$ remains detectable across the whole tested ladder." Experiment 2 §4 gives
   $|z_4| = 2.3$ at $\kappa=0.75,\lambda=50$ and $1.6$ at $\lambda=100$ — both **below** the
   $3\sigma$ line the Fig. 3 caption itself draws. Detection holds for $\lambda\le20$ at
   $\kappa=0.75$, and $\lambda\le5$ at $\kappa=1$. This needs one sentence rewritten.
3. **The box's agyrotropy is bounded by the rejected mass** (§4.5): any bounded directional
   statistic obeys $|\Delta| \le 2G\varepsilon/(1-\varepsilon)$, and $a_4 \propto \varepsilon$ is
   confirmed empirically. So the angular damage, unlike the *quantile* damage, is already
   controlled by the number the paper reports. This tempers claim 2b's standing as an
   independent hazard — and is the honest framing.

**Should it enter the current revision?** **Recommendation B** — a small, self-contained
addition now; the full comparison as a follow-up. Reasoning in §13.

---

# 2. Existing measured motivation

Only the Experiment-2 facts that bear on the question, from
`experiments/exp2_cap_characterization/results/exp2_table.md`:

| Fact | Value | Why it matters here |
|---|---|---|
| Capped law is a distinct conditional target | verified bitwise, 240/240 | Cap geometry is part of the target definition, not an implementation detail |
| $\mathrm{TV} = 1-p$ exactly | closed form | **This is geometry-blind** (§9). It cannot rank caps |
| Rejected mass decays as $\lambda^{-(2\kappa-1)}$ | slopes to 0.5% | Exponent is set by the tail, not the geometry — no geometry fixes heavy tails |
| Small TV ≠ faithful tail | $\kappa=1.5,\lambda=50$: TV $=6.3\times10^{-4}$, p99.9 speed $-24\%$ | The template for the new anisotropy result |
| Box is a cube in $\mathbf{U}$, $\theta$ cancels | exact | Confirmed against code, §3.1 |
| Four-fold azimuthal modulation | $a_4=-0.052$, $z_4=-11.6$ at $\kappa=0.75,\lambda=3$ | Real, but only where $\varepsilon$ is already large (§4.5) |
| Capped variance at $\kappa\le3/2$ is a property of $\lambda$ | $1.28\to178$ | Same logic recurs for the sphere's anisotropy bias |

---

# 3. Mathematical formulation

## 3.1 Repository convention (checked against code, not assumed)

From `bi_kappa_distribution.H:282–305` and `main.tex` Eqs. (16)–(19):

$$T=R^2\sim\beta'(3/2,\kappa-1/2),\qquad \mathbf{U}=R\,\mathbf{n},\ \ \mathbf{n}\sim\mathrm{Unif}(S^2)\perp R,$$
$$\mathbf{V}=A\mathbf{U},\qquad A=\sqrt{\kappa}\,\mathrm{diag}(\theta_\perp,\theta_\perp,\theta_\parallel),\qquad
p_{\mathbf U}(\mathbf u)\propto(1+|\mathbf u|^2)^{-(\kappa+1)}.$$

$p_{\mathbf U}$ is isotropic; that isotropy is the entire structure a cap can preserve or destroy.

## 3.2 Geometry A — component-wise box (what ships today)

Shipped predicate: $|v_x|/\theta_\perp\le\lambda$, $|v_y|/\theta_\perp\le\lambda$,
$|v_\parallel|/\theta_\parallel\le\lambda$. Substituting $v_x=\sqrt\kappa\,\theta_\perp u_x$ etc.:

$$\boxed{\ \max_i|u_i|\le c,\qquad c=\lambda/\sqrt{\kappa}\ }$$

A **cube** of half-side $c$. The $\theta$'s cancel identically — confirming Exp-2 finding 4 from
the source. In $\mathbf V$-space: a rectangular box with half-widths
$(\lambda\theta_\perp,\lambda\theta_\perp,\lambda\theta_\parallel)$.
Max speed $\lambda\sqrt{2\theta_\perp^2+\theta_\parallel^2}$.

**Symmetry group:** $O_h$ (octahedral) — finite, *not* containing rotations about $\hat{\mathbf b}$.
This single fact drives everything: it retains the 4-fold subgroup $D_{4h}$ about $\hat{\mathbf b}$,
which is why the leading surviving azimuthal harmonic is $m=4$ and not $m=2$.

## 3.3 Geometry B — normalized-radial / ellipsoidal cap

$$\frac{v_x^2+v_y^2}{\theta_\perp^2}+\frac{v_\parallel^2}{\theta_\parallel^2}\le\lambda^2
\iff \kappa|\mathbf u|^2\le\lambda^2 \iff \boxed{\ R\le R_{\max}=\lambda/\sqrt\kappa\ },\quad T_{\max}=\lambda^2/\kappa.$$

Note the elegant relation: **at equal $\lambda$, geometry B is exactly the ball inscribed in
geometry A's cube.** Ellipsoid in $\mathbf V$ with semi-axes
$(\lambda\theta_\perp,\lambda\theta_\perp,\lambda\theta_\parallel)$; sphere in $\mathbf U$.
Max speed $\lambda\max(\theta_\perp,\theta_\parallel)$; per-component bounds also hold.

**Acceptance in closed form** — no angular integral needed, unlike the box:
$$p=P(T\le T_{\max})=I_z\!\left(\tfrac32,\kappa-\tfrac12\right),\qquad z=\frac{T_{\max}}{1+T_{\max}}=\frac{\lambda^2}{\kappa+\lambda^2}.$$
Verified: $\kappa=0.75,\lambda=3$ gives $0.402136$ analytic vs $0.402209$ empirical ($N=4\times10^6$).

The box's acceptance requires $E_M[I_z(3/2,\kappa-1/2)]$ over the density of $M=\max_i|n_i|$
(Exp 2). **The ellipsoid's is a single incomplete beta.** That is a real, if minor, simplification.

As $\lambda\to\infty$, $1-p\propto\lambda^{-(2\kappa-1)}$ — the *same exponent as the box*. Any
region squeezed between two balls has the same exponent; **only the prefactor is
geometry-dependent.** No geometry rescues heavy tails.

## 3.4 Geometry C — physical-speed sphere $|\mathbf V|\le v_{\max}$

$|\mathbf V|^2=\kappa R^2(\theta_\perp^2\sin^2\Theta+\theta_\parallel^2\cos^2\Theta)$, $\Theta$ the
polar angle from $\hat{\mathbf b}$ in $\mathbf U$-space. So

$$R\le R_{\max}(\Theta)=\frac{v_{\max}}{\sqrt\kappa\,\sqrt{\theta_\perp^2\sin^2\Theta+\theta_\parallel^2\cos^2\Theta}}.$$

The acceptance radius depends on $\Theta$ but **not on $\varphi$**. Hence gyrotropy is exact and
the angular law is distorted in the polar direction only. Radius and direction become dependent.
Reduces to geometry B iff $\theta_\perp=\theta_\parallel$.

## 3.5 The non-vanishing anisotropy bias (the key derivation)

For a cap with $\Theta$-dependent radius $R_{\max}(\Theta)$ and $\kappa<3/2$ (where
$E[R^2]=\infty$), the truncated second moment is tail-dominated,
$\int^{R_{\max}}r^2\cdot r^2(1+r^2)^{-(\kappa+1)}dr\sim R_{\max}^{3-2\kappa}$, giving

$$\lim_{\varepsilon\to0}\frac{E[u_z^2]}{E[u_x^2]}
=\frac{\displaystyle\int_0^1 \mu^2\,W(\mu)\,d\mu}{\displaystyle\int_0^1 \tfrac12(1-\mu^2)\,W(\mu)\,d\mu},
\qquad W(\mu)=\left(\theta_\perp^2(1-\mu^2)+\theta_\parallel^2\mu^2\right)^{-(3-2\kappa)/2}.$$

The overall scale cancels from the *limit*: **the limiting bias is independent of $v_{\max}$**, though the finite-cap ratio is not (§0, C2). The derivation above requires $3-2\kappa>0$; at $\kappa=3/2$ the divergence is logarithmic and the limit is $1$ instead (§0, C1). Predicted vs measured
($\theta_\parallel/\theta_\perp=2$; target value is exactly 1):

| $\kappa$ | 0.60 | 0.75 | 1.00 | 1.25 | 1.40 | 1.49 |
|---|---|---|---|---|---|---|
| predicted plateau | 0.5907 | **0.6446** | 0.7462 | 0.8641 | 0.9434 | 0.9942 |

Measured at $\kappa=0.75$: $0.6443\ (\varepsilon=10^{-2})$, $0.6426\ (10^{-3})$,
$0.6622\ (10^{-4})$, $0.6380\ (10^{-5})$ — a plateau, not a decay. **Confirmed.**

For $\kappa>3/2$ the moment is core-dominated and the bias *does* vanish as $\varepsilon\to0$
(measured $\kappa=5$: $0.694\to0.918\to0.983$ for $\varepsilon=10^{-1},10^{-2},10^{-3}$).
The hazard is specific to $1/2<\kappa<3/2$ — precisely the regime this paper opened up. At the
endpoint $\kappa=3/2$ it degenerates to an $O(1/\log v_{\max})$ decay (§0, C1).

## 3.6 Geometry D — gyrotropic cylinder

$v_\perp\le v_{\perp,\max}$, $|v_\parallel|\le v_{\parallel,\max}$ gives
$R\le\min\!\big(a/\sin\Theta,\ b/|\cos\Theta|\big)$ with $a=v_{\perp\max}/(\sqrt\kappa\theta_\perp)$,
$b=v_{\parallel\max}/(\sqrt\kappa\theta_\parallel)$ — again $\Theta$-dependent, $\varphi$-free.
Same structural class as C: gyrotropy exact, angular law and anisotropy distorted, and the same
non-vanishing plateau on $1/2<\kappa<3/2$ (measured $1.0918/1.0991/1.0997$ at
$\varepsilon=10^{-1},10^{-2},10^{-3}$, $\kappa=0.75$). Its unique feature is **independent**
$\perp$/$\parallel$ bounds; B bounds both but locks their ratio to $\theta_\perp:\theta_\parallel$.

## 3.7 Geometry E — regularized Kappa (comparison class only, not a cap)

Scherer, Lazar, Husidic & Fichtner 2019 (ApJ 880, 118), **Eq. (3)**, read from the PDF:

$$f_{\rm RBK}=n_0N_{\rm RBK}\left(1+\frac{w_\parallel^2}{\kappa}+\frac{w_\perp^2}{\kappa}\right)^{-\kappa-1}
e^{-\alpha_\parallel^2 w_\parallel^2-\alpha_\perp^2 w_\perp^2},\qquad
w_\parallel=\frac{v_\parallel}{\Theta_\parallel},\ w_\perp=\frac{v_\perp}{\Theta_\perp},$$

with "two distinct positive cutoff parameters $\alpha_\parallel\ne\alpha_\perp$."

**This is directly load-bearing for the audit.** The regularizer is a function of the
**normalized** components, so its level sets are **ellipsoids in physical velocity space** — not
speed spheres, not component boxes. At $\alpha_\parallel=\alpha_\perp$ they coincide exactly with
the bi-Kappa's own density level sets: the smooth analogue of geometry B.

Differences from every hard cap: support stays $\mathbb{R}^3$; the law is a *different physical
model* with its own moments, not a conditioned version of the bi-Kappa; there is no rejected
mass; and samplers are published (ZUM 2026 gives two rejection methods).

## 3.8 Numerical verification of §3.2–3.6

Audit probe, $N=4$–$20\times10^6$ per case, isotropic-$\mathbf U$ draws with each region applied.
Target values: $a_4=0$, $E[u_z^2]/E[u_x^2]=1$.

At $\kappa=2$, $\lambda=3$, $(\theta_\perp,\theta_\parallel)=(1,2)$:

| region | accept | $E[u_z^2]/E[u_x^2]$ | $z_4$ | $E[u_x^4-6u_x^2u_y^2+u_y^4]$ |
|---|---|---|---|---|
| A cube | 0.919 | **1.0007** | **−25.5** | −3.18e−01 |
| B ball | 0.876 | **1.0005** | +0.1 | −4.7e−04 |
| C sphere | 0.876 | **1.4486** | +0.1 | −2.3e−04 |
| D cylinder | 0.904 | **0.9361** | −0.2 | −1.5e−03 |

The pattern is exactly complementary and holds at every $\kappa$ tested (0.75, 1.5, 2, 5).

---

# 4. Preservation-property analysis

## 4.1 Gyrotropy of the distribution function

A cap region $C$ preserves the continuous gyrotropy of $f$ **iff $C$ is invariant under rotation
about $\hat{\mathbf b}$.** ($f_C=f\mathbf 1_C/p$, $f$ already gyrotropic, so $f_C$ is gyrotropic iff
$\mathbf 1_C$ is.) One line, and it is the entire design rule. B, C, D pass; **A is the only one
of the four that fails**, retaining $D_{4h}$ and hence a leading $m=4$ harmonic.

## 4.2 Gyrotropy of the pressure tensor — where the current manuscript is imprecise

The cube is invariant under **permutation of the three normalized axes** (equal half-side $c$ in
all three, because $\theta$ cancels) and under sign flips. With $p_{\mathbf U}$ isotropic this
forces $E[u_x^2]=E[u_y^2]=E[u_z^2]$ exactly and all off-diagonals zero. Therefore under the box:

$$P_{xx}=P_{yy},\quad P_{ij}=0\ (i\ne j),\quad \frac{E[v_\parallel^2]}{E[v_x^2]}=\left(\frac{\theta_\parallel}{\theta_\perp}\right)^2\ \textbf{exactly}.$$

**The box-capped distribution is exactly gyrotropic at the pressure-tensor level and reproduces
the imposed anisotropy exactly.** Verified: $E[u_z^2]/E[u_x^2]=1.0007$ (κ=2), $1.0016$ (κ=0.75).
Every standard agyrotropy measure returns zero. The breaking appears first at fourth order.

## 4.3 The ellipsoidal cap's moment-ratio theorem (proved, as §4 of the handoff asked)

Under B the conditional law of $\mathbf U$ is isotropic (ball ∩ isotropic density), so
$E[u_x^2]=E[u_y^2]=E[u_z^2]=E[R^2\mid R\le R_{\max}]/3$ and

$$\frac{E[v_\parallel^2]}{E[v_{x}^2]}=\left(\frac{\theta_\parallel}{\theta_\perp}\right)^2\quad\textbf{exactly, for every }R_{\max}<\infty\text{ and every }\kappa>1/2,$$

**including $\kappa\le3/2$, where the untruncated moments do not exist.** *Proved.*

**Per the handoff's explicit request, the distinction matters:** for the perpendicular
*magnitude*, $E[v_\perp^2]=2E[v_x^2]$, so
$E[v_\parallel^2]/E[v_\perp^2]=\tfrac12(\theta_\parallel/\theta_\perp)^2$. Always name which.

**Disproof of the handoff's Claim 4 as stated:** the box preserves this ratio *equally exactly*
(§4.2). The ellipsoid is **not** more faithful in anisotropy than the box. It is more faithful
than the **sphere and cylinder** — which is the real result, and the opposite of the assumption.

## 4.4 Summary of what each geometry costs

| | preserves $m=4$-free angular law | preserves $T_\parallel/T_\perp$ |
|---|---|---|
| A box | ✗ (at $O(\varepsilon)$) | ✓ exactly |
| B ellipsoid | ✓ exactly | ✓ exactly |
| C sphere | ✓ gyrotropy; ✗ polar law | ✗ (plateau on $1/2<\kappa<3/2$; $O(1/\log v_{\max})$ at $\kappa=3/2$) |
| D cylinder | ✓ gyrotropy; ✗ polar law | ✗ (plateau on $1/2<\kappa<3/2$) |

**B is the only geometry that pays nothing structurally.** Its whole distortion is radial.

## 4.5 The $O(\varepsilon)$ bound on directional damage

For any bounded $g$ with $|g|\le G$ and any region $C$ with $P(C)=1-\varepsilon$:
$$\left|E_C[g]-E[g]\right|\le \frac{2G\varepsilon}{1-\varepsilon}.$$
So **every** bounded directional statistic — $a_4$, any gyrotropy index — is $O(\varepsilon)$ for
*every* geometry. Confirmed: at $\kappa=5$, $a_4=-0.036,-0.0044,-0.0004$ for
$\varepsilon=10^{-1},10^{-2},10^{-3}$: linear in $\varepsilon$.

**Consequence, and it cuts against claim 2b:** the box's agyrotropy is *bounded by the number the
paper already reports*. It is not an independent hazard the way the tail-quantile error is — that
error is unbounded because quantiles are unbounded functionals. The anisotropy plateau of §3.5
is a genuine third category: **bounded functional, yet $O(1)$ error**, because at $\kappa\le3/2$
the reference moment is infinite and the bound above compares against a divergent quantity.

---

# 5. Closest mathematical / statistical prior art

**The stochastic representation is the whole story and it is textbook.** An elliptically
contoured $\mathbf X=\boldsymbol\mu+A R\mathbf D$ with $\mathbf D\sim\mathrm{Unif}(S^{d-1})$
independent of $R$ — Cambanis, Huang & Simons (1981); Fang, Kotz & Ng (1990). Conditioning on
$R\le R_{\max}$ touches only the radial factor. Independence and the uniform angular law survive
by construction. **Claims 2a and 3 are immediate corollaries.** Nothing there is new.

Searching the truncated-elliptical literature is instructive mainly for what it shows is *hard*:
essentially all of it addresses **rectangular / polytope** truncation — Ho, Lin, Chen & Wang
(2011) [6]; Li & Ghosh (2015) [5]; Botev & L'Ecuyer (2015) [12,17]; Galarza et al. (2020, 2021)
[1,8]; Morán-Vásquez et al. (2019) [4]; Valeriano et al. (2021) [2], with the R package
`relliptical`; Wu et al. (2024) [3]. Arismendi & Broda (2017) [9] handle moments of **quadratic
forms**, the closest to Mahalanobis-radius truncation.

**The reading is unambiguous: ellipsoidal truncation gets no dedicated literature because it is
the trivial case.** Effort goes to rectangular/polytope truncation precisely because it destroys
the radial–angular factorization and therefore needs MCMC, tilting or slice sampling. Chiron et
al. (2023) [11] use the radial × von Mises–Fisher decomposition explicitly as a routine tool.

**On direct radial sampling (Claim 5).** $T\sim\beta'(a,b)\Rightarrow Y=T/(1+T)\sim\mathrm{Beta}(a,b)$
is standard (Johnson, Kotz & Balakrishnan, Pearson Type VI). Inverse-CDF for a truncated law,
$X=F^{-1}(UF(b))$, is Devroye (1986). Composing them gives the truncated radial sampler in two
textbook steps via `betaincinv`. **This is not a new algorithm and must not be presented as one.**

One genuine numerical caveat, already known in this repo: $Y=T/(1+T)$ rounds to exactly $1$ in
double precision at small $\kappa$ (Experiment 1), so a naive implementation would need the
complement branch (`betainccinv` / $I_{1-z}(b,a)$) — the same rule Exp 2 already follows.

---

# 6. Closest plasma / PIC prior art

**Read from the local primary PDFs (`paper/reference/`), grepped for cutoff / truncation /
maximum velocity / gyrotropy / azimuthal.**

| Source | Hard velocity cap? | Geometry | Finding |
|---|---|---|---|
| Abdul & Mace 2014 (CPC 185, 2383) | none found | — | 1-D Bailey polar transform; no cap discussion |
| Abdul & Mace 2015 (PoP 22, 102107) | none found | — | no cutoff-geometry discussion |
| Abdul 2013 MSc / 2018 PhD theses | none found | — | only $v_{\min},v_{\max}$ as *integration limits*; "cut-off" is always a **wave cut-off frequency** |
| Zenitani & Nakano 2022, 2023; Zenitani 2025 | none found | — | rejection method + beta-prime; no hard cap geometry |
| Zenitani, Usami & Matsukiyo 2026 | **smooth** cutoff | regularized-Kappa | "Kappa distribution with a high-energy cutoff"; acceptance efficiency vs $v_c$; enables $0<\kappa<3/2$ |
| Scherer, Fichtner, Fahr & Lazar 2019 (ApJ 881, 93) | **smooth** | $e^{-\alpha^2w^2}$ on speed | isotropic RKD; truncates superluminal contribution |
| **Scherer, Lazar, Husidic & Fichtner 2019 (ApJ 880, 118)** | **smooth** | **normalized/ellipsoidal, Eq. (3)** | anisotropic RKD: $e^{-\alpha_\parallel^2w_\parallel^2-\alpha_\perp^2w_\perp^2}$, $w=v/\Theta$ per component |

**Two conclusions, and the second is the important one.**

1. **No hard velocity cap of any geometry appears in the Kappa-loading literature I hold.** The
   component-wise box in this repository appears to be a local implementation choice with no
   published precedent — and correspondingly no published analysis. Consistent with
   `step1_claim_audit.md` §8.3.
2. **The anisotropic-Kappa literature already uses the normalized/ellipsoidal geometry** for its
   principled smooth cutoff (Scherer et al. 2019b Eq. 3). So "the right bounding geometry for an
   anisotropic Kappa is the normalized ellipsoid, not the physical speed sphere" is **already
   implicitly established practice** in the smooth setting. This substantially deflates any
   novelty claim for geometry B and simultaneously *strengthens* the practical recommendation:
   the hard-cap analogue of accepted practice is B, and the shipped box is not it.

Scherer et al. 2019b also state gyrotropy as a structural requirement ("gyrotropic distributions
typical in magnetized plasmas... $P_{22}=P_{33}$"), which is the right citation for why breaking
it matters.

**Agyrotropy as a numerical artifact:** searched; **not found**. The agyrotropy literature
(Scudder & Daughton, Aunai et al., Swisdak) treats it as a *physical* reconnection diagnostic. No
paper reporting spurious agyrotropy from particle initialization surfaced. Given §4.2 — the box
produces *zero* pressure-tensor agyrotropy — this negative result is also less consequential than
the prior audit assumed.

---

# 7. Claim-by-claim novelty table

| # | Claim | Verdict | Basis |
|---|---|---|---|
| 1 | Box truncation induces $m=4$ modulation, breaks gyrotropy | **PLAUSIBLY NEW NUMERICAL CHARACTERIZATION** — mathematically obvious ($O_h\not\supset$ rotations); no published plasma statement found; **but** damage is $O(\varepsilon)$ (§4.5) and invisible to the pressure tensor (§4.2), so it must be scoped far more carefully than §IV.E currently does | §3.2, §4.1–4.2, §4.5, §6 |
| 2a | Ellipsoidal cap preserves gyrotropy — mathematics | **STANDARD CONSEQUENCE / NOT NOVEL** | Cambanis et al. 1981; Fang–Kotz–Ng 1990 |
| 2b | …in Kappa/bi-Kappa literature | **ESTABLISHED PRIOR ART** | Scherer et al. 2019b Eq. (3) already regularizes in normalized components |
| 2c | …as a hard-cap PIC loading policy | **PLAUSIBLY NEW IN PLASMA LOADING** — but trivial given 2a+2b; not worth a claim | §6 |
| 3 | Radial conditioning preserves angular law + $R$–direction independence | **STANDARD CONSEQUENCE / NOT NOVEL** — immediate from radial–angular factorization | §5 |
| 4 | Radial cap preserves anisotropy better than box/sphere/cylinder | **Partly FALSE as stated.** vs **box**: false — both exact (§4.2, §4.3). vs **sphere/cylinder**: true, and the non-vanishing $1/2<\kappa<3/2$ plateau is a **STRONG NOVELTY CANDIDATE** | §3.5, §4.3 |
| 5 | Direct conditional radial sampler, no outer redraw | **STANDARD CONSEQUENCE / NOT NOVEL** — Beta-prime→Beta + inverse-CDF, both textbook | §5 |
| 6 | Systematic matched comparison of cap geometries | **PLAUSIBLY NEW IN PLASMA LOADING** — none found; modest weight on its own | §6 |
| 7 | Radial cap is "optimal" | **NOT VIABLE AS STATED** — all $f$-divergences are geometry-blind (§11); a precise multi-objective statement is available instead | §11 |

**Headline correction to the prior audit.** `step1_claim_audit.md` §8.3 called the box/gyrotropy
result "the strongest genuinely new empirical result in the paper." This audit finds it weaker
than believed (bounded by $\varepsilon$; invisible to the standard diagnostic) and finds the
**speed-cap anisotropy plateau** stronger.

---

# 8. Geometry property table

| Property | A box | B ellipsoid | C phys. sphere | D cylinder | E reg-Kappa |
|---|---|---|---|---|---|
| Gyrotropy of $f$ | **✗** ($D_{4h}$, $m=4$) | ✓ exact | ✓ exact | ✓ exact | ✓ exact |
| Pressure tensor gyrotropic | ✓ exact | ✓ exact | ✓ exact | ✓ exact | ✓ exact |
| Full normalized angular law | ✗ | **✓ exact** | ✗ (polar) | ✗ (polar) | ✗ (polar unless $\alpha_\parallel=\alpha_\perp$) |
| $R$ ⟂ direction independence | ✗ | **✓ exact** | ✗ | ✗ | ✗ |
| Scale-implied second-moment ratio | ✓ exact | ✓ exact | **✗ plateau on $1/2<\kappa<3/2$** | **✗ plateau** | ✗ by design |
| Hard bound on $\vert\mathbf V\vert$ | ✓ $\lambda\sqrt{2\theta_\perp^2+\theta_\parallel^2}$ | ✓ $\lambda\max\theta$ | **✓ exactly $v_{\max}$** | ✓ $\sqrt{v_{\perp\max}^2+v_{\parallel\max}^2}$ | **✗ unbounded support** |
| Independent ⊥/∥ bounds | ✓ (ratio locked) | ✓ (ratio locked) | ✗ | **✓ free** | ✓ free |
| Closed-form acceptance | angular integral | **single $I_z$** | angular integral | angular integral | n/a |
| Direct conditional sampler | ✗ | **✓ trivial** | ✗ | ✗ | published (ZUM 2026) |
| Known prior art | **none found** | smooth analogue: Scherer 2019b | none found | none found | Scherer; ZUM 2026 |

---

# 9. Fair comparison framework

**The crucial negative result first.** For $f_C=f\mathbf 1_C/p$ and *any* $f$-divergence
$D_\phi$:
$$D_\phi(f_C\Vert f)=\int f\,\phi\!\left(\tfrac{f_C}{f}\right)=p\,\phi(1/p)+(1-p)\,\phi(0),$$
a function of $p$ **alone**. In particular $\mathrm{TV}=1-p$ and $\mathrm{KL}(f_C\Vert f)=-\log p$
for every geometry; $\mathrm{KL}(f\Vert f_C)=\infty$ always.

> **No $f$-divergence can rank cap geometries at matched rejected mass.** The paper's own
> Eq. (20) result is a special case. Any comparison must use **structural** statistics —
> directional harmonics, moment ratios, quantiles — or a transport metric.

| Match | Question it answers | Verdict |
|---|---|---|
| **A — equal rejected mass $\varepsilon$** | "Same cost, which geometry preserves most structure?" | **Primary.** Only sound at fixed $\varepsilon$; must use structural metrics per above |
| **B — equal hard $\vert\mathbf V\vert\le v_{\max}$** | "Simulation demands a true speed limit; which policy distorts least?" | **Primary and the practically decisive one.** All four can satisfy it; see below |
| C — equal cost | timing | **Secondary.** Mean attempts $=1/p$, so equal cost *is* Match A. Redundant |
| D — equal $R_{\max}$ | theory only | **Weak.** Corresponds to no physical constraint; use only for exposition |

**Match B has a clean trade-off theorem.** Among regions $C\subseteq\{|\mathbf V|\le v_{\max}\}$,
the sphere is the *maximal* one, so it uniquely minimizes rejected mass. An ellipsoid meeting the
same guarantee must be inscribed, $\lambda=v_{\max}/\max(\theta_\perp,\theta_\parallel)$, and pays
more rejection. Measured price ($\theta_\parallel/\theta_\perp=2$):

| $\kappa$ | $\varepsilon_{\rm ellip}/\varepsilon_{\rm sphere}$ | sphere $T$-ratio error | ellipsoid |
|---|---|---|---|
| 0.75 | 1.21 | **0.646** | 1.000 |
| 1.5 | 2.00 | 0.898 | 0.998 |
| 2 | 2.5–2.8 | 0.981 | 1.002 |
| 5 | 6.1 | 0.995 | 0.999 |

**So: you cannot simultaneously minimize rejected mass at fixed $v_{\max}$ and preserve the
angular law and anisotropy, unless $\theta_\perp=\theta_\parallel$. The price of correctness is a
factor 1.2–6 in rejected mass — and it buys removal of a 35% anisotropy error at $\kappa=0.75$.**
That is a genuine, quantified, practically actionable design result.

---

# 10. Minimal experiment design (justified only at Recommendation-B scale)

Not run. Designed to be the smallest thing that settles the question.

- **$\kappa\in\{0.75,\,1.5,\,2,\,5\}$** — spans divergent / boundary / finite second moment.
  Drop 0.51 and 10: the physics is monotone and already covered by Exp 1/2.
- **$\theta_\parallel/\theta_\perp=2$** only, plus a single $(1,1)$ isotropic control (where
  B ≡ C, a free correctness check).
- **Four geometries** at **Match A** ($\varepsilon\in\{10^{-1},10^{-2},10^{-3}\}$) and
  **Match B** ($v_{\max}\in\{10,50,200\}$).
- $N=10^6$ accepted × 5 seeds. $\hat b=\hat z$ only — the cap is applied pre-rotation
  (`bi_kappa_distribution.H:307` before `:321`), so multiple $\hat{\mathbf b}$ adds nothing.
  **Exp 1 already validated the rotation.**
- **Outputs (5 only):** $a_4$ with $z_4$; $E[v_\parallel^2]/E[v_x^2]$; p99/p99.9 speed ratio;
  rejected mass; max realized $|\mathbf V|$.
- **Cost:** ≈ $2.4\times10^8$ draws — smaller than Experiment 2. Roughly 1 day including analysis.

**Deliverable: one figure (3 panels — $a_4$, anisotropy ratio, p99.9 — vs $\varepsilon$, four
geometries) and one table.** Not more.

Drop from the handoff's list: runtime (it is $1/p$, already known), finite-precision failures
(Exp 4 covers it, and it is geometry-independent), radial-law distortion (analytic for B,
uninteresting for others), multiple $\hat{\mathbf b}$, and $\kappa=0.75$–$10$ full sweeps.

---

# 11. The optimality question

**What is known and applies:**
- Density level sets of an elliptical law are ellipsoids; the HDR at any probability level is
  therefore an ellipsoid, and it is the **minimum-volume** region carrying that mass (standard;
  Anderson's lemma / Hyndman 1996 — *not verified from primary source*, §15).
- Minimum volume is **not** minimum distributional distortion. Do not conflate them.

**What is not true:** "the ellipsoid minimizes the distortion at fixed rejected mass" is
**vacuous** — §9 shows every geometry ties in every $f$-divergence.

**What *can* be stated precisely**, and should be, if anything is:

> *(Symmetry, not optimality.)* Among cap regions with a given rejected mass, the normalized ball
> is the unique one (up to null sets) invariant under the full symmetry group $O(3)$ of
> $p_{\mathbf U}$; it is therefore the only choice leaving the angular law and $R$–direction
> independence exactly intact, and the only one whose distortion is purely radial.

> *(Constrained trade-off — §9.)* Among regions guaranteeing $|\mathbf V|\le v_{\max}$, the sphere
> uniquely minimizes rejected mass but incurs an anisotropy bias that does not vanish with
> $v_{\max}$ for $1/2<\kappa<3/2$; the inscribed ellipsoid removes that bias at a bounded,
> computable cost in rejected mass.

**Verdict: do not claim optimality.** Claim a symmetry characterization and a quantified
trade-off. Both are defensible; neither needs a new theorem.

---

# 12. Effect on the current paper's thesis

The current thesis — *reliable bi-Kappa loading requires validating the sampler, the bounding
policy, and the implementation* — **already contains the bounding policy.** §IV is a
characterization of one policy. A short comparison makes §IV land better because it converts a
purely negative section ("our cap is bad") into a constructive one ("…and here is the geometry
that isn't, and why").

- **Small addition (§13-B): stronger and more coherent.** §IV currently leaves the reader with a
  problem and no alternative, which is a weak place for a methods paper to end a section.
- **Full comparison inside this revision: broader, defensible, but unfocused and risky.** It adds
  a second thesis (cap-geometry design) to a paper whose contribution list is already five items,
  and it invites a fresh round of refereeing on material no reviewer asked for.

**Relevant fact for scoping:** `reviewer_response_matrix.md` line 60 records the referee point
that "Sec. III.E implements a **component-wise box cutoff, not a radial speed cutoff**," and line
633 lists "anisotropic cutoff geometry" among the genuine added implementation concerns. **Cap
geometry is inside the reviewers' expressed concern surface.** A short, decisive treatment is
*responsive*; a full comparative study is not what they asked for.

---

# 13. Recommendation

## **B. SMALL ADDITION ONLY; SAVE THE FULL STUDY FOR A FOLLOW-UP.**

**Why not A (full study now).** The revision is evidence-complete, builds clean, and is one
rebuttal letter from resubmission. A full comparison needs new code, a new experiment, a new
figure and table, ~3 manuscript pages, and a second thesis. Against that: two of the seven claims
are standard mathematics, one is already established practice (Scherer 2019b), and one is false
as stated. The genuinely new piece — the $1/2<\kappa<3/2$ speed-cap plateau — does **not** require
the full four-geometry apparatus to state.

**Why not C/D (nothing now).** Three items are cheap, responsive, and fix real defects:
two of them are **corrections to claims already in the manuscript**, which must happen regardless.

**What to add now — ~2/3 page, no new experiment, no new code, no API change:**

1. **Fix §IV.E (required, independent of everything else).** State that the box preserves the
   pressure tensor's gyrotropy and the anisotropy ratio *exactly*, and that the modulation is a
   fourth-order property invisible to standard agyrotropy measures. This is more defensible
   *and* a better result. Repair the "detectable across the whole tested ladder" sentence against
   Exp 2 §4 ($|z_4|=2.3,\,1.6$ at $\lambda=50,100$).
2. **Add a short §IV.F, "The geometry is a choice"** (~2 paragraphs + 4-row table): the box is a
   cube in $\mathbf U$; the normalized ellipsoid $\kappa|\mathbf u|^2\le\lambda^2$ preserves
   gyrotropy, the full angular law and $R$–direction independence, with acceptance
   $I_z(3/2,\kappa-1/2)$, $z=\lambda^2/(\kappa+\lambda^2)$; a physical-speed sphere preserves
   gyrotropy but biases the capped second-moment ratio away from the scale-implied value, and for
   $1/2<\kappa<3/2$ that bias does
   not vanish as $v_{\max}\to\infty$ (quote 0.645 at $\kappa=0.75$). One sentence noting this is
   the hard-cap analogue of the geometry Scherer et al. (2019b) already use smoothly.
   Cite `scherer2019anisotropic` (already in `refs.bib`).
3. **One sentence in §X** flagging the systematic matched comparison as future work.

All three are derivable analytically and verified (§3.8); **no new experiment is required to
state any of them.** The numbers quoted come from closed forms plus the audit probes, which
should be promoted into a small reproducible script if item 2 is adopted.

**The follow-up.** "Symmetry-preserving finite-velocity bounding for heavy-tailed kinetic
particle loading" is viable *only* if it carries more than the comparison. Minimum for
substance: (i) the four-geometry matched study of §10; (ii) the $1/2<\kappa<3/2$ plateau derived
in closed form for all three $\Theta$-dependent geometries; (iii) a PIC/hybrid demonstration that
one of these biases measurably changes a physical outcome — an instability threshold or growth
rate sensitive to $T_\parallel/T_\perp$. **Without (iii) it is a note, not a paper**, and the
$O(\varepsilon)$ bound of §4.5 means the box-agyrotropy angle alone will not carry it.

---

# 14. If added: manuscript structure

Only §IV changes. No section moves; no renumbering beyond §IV.

```
IV. The optional component-wise cap
    A. The capped law                                  unchanged
    B. Cost and distortion are the same number         + one line: all f-divergences are
                                                         geometry-blind, so the comparison
                                                         must be structural
    C. A small rejected mass ≠ a faithful tail         unchanged
    D. A finite variance that belongs to the box       unchanged
    E. What the cap does to gyrotropy                  RETITLED + CORRECTED (§13 item 1)
    F. The geometry is a choice                        NEW, ~2 paragraphs + 4-row table
    G. Consequences for the released implementation    was F; unchanged
```

The handoff's proposed ten-section restructure (separate "Finite-velocity bounding strategies"
part) is **rejected for this revision**: it elevates the cap to a co-equal thesis, contradicts
§IX's framing of the cap as optional and discouraged, and reopens sections the referees accepted.

---

# 15. References and unresolved verification gaps

**Verified from primary source (PDFs in `paper/reference/`, read directly this session):**
- Scherer, Lazar, Husidic & Fichtner 2019, ApJ 880, 118 — **Eq. (3) read from the PDF page**;
  anisotropic RKD regularizes in normalized components. Load-bearing for §6 and §7 claim 2b.
- Scherer, Fichtner, Fahr & Lazar 2019, ApJ 881, 93 — isotropic RKD, $e^{-\alpha^2w^2}$ on speed.
- Zenitani, Usami & Matsukiyo 2026, JGR 131 — high-energy-cutoff Kappa; acceptance vs $v_c$.
- Abdul & Mace 2014, 2015; Abdul 2013/2018 theses; Zenitani & Nakano 2022/2023; Zenitani 2025 —
  grepped for cutoff/truncation/gyrotropy: **no hard velocity cap of any geometry**.

**Named but NOT opened as primary sources — must be verified before any of this is cited:**
1. Cambanis, Huang & Simons (1981), *J. Multivar. Anal.* — elliptical stochastic representation.
2. Fang, Kotz & Ng (1990), *Symmetric Multivariate and Related Distributions*.
3. Tallis (1963), *Ann. Math. Statist.* — elliptical and radial truncation in normal populations.
   **Potentially the closest single prior-art item to claims 2a/3/4; unverified.**
4. Anderson (1955) lemma; Hyndman (1996) HDR — for §11's minimum-volume statement.
5. Devroye (1986) — inverse-CDF for truncated laws; chapter/page not pinned.
6. Johnson, Kotz & Balakrishnan Vol. 2 — Beta-prime ↔ Beta (Pearson Type VI).
7. Scudder & Daughton (2008); Aunai et al. (2013); Swisdak (2016) — agyrotropy measures. Needed
   only if §IV.E's corrected claim names them.

**Searched, nothing found (negative results, recorded deliberately):**
- Spurious agyrotropy from particle initialization or velocity-space truncation — no hit.
- Any comparison of hard velocity-cap geometries in plasma particle loading — no hit.
- Any report of $m=4$ azimuthal structure from Cartesian velocity-space truncation — no hit.

**Blocked / not completed:** PIC-code loader survey (WarpX, Smilei, OSIRIS, EPOCH, VPIC —
does any expose a velocity cap, and with what geometry?), and the Vlasov Cartesian
velocity-grid-boundary literature (Vlasiator). Both agents died on the API limit. **Neither
affects the recommendation**, but both should be closed before a follow-up paper is drafted.

**Audit probe scripts** (`/tmp/audit_check.py`, `/tmp/audit_probe.py`, `/tmp/audit_plateau.py`,
`/tmp/audit_matchb.py`) are scratch, not committed. If §13 item 2 is adopted, the plateau formula
and the acceptance closed form should be promoted into a small tracked script so the quoted
numbers have provenance to the same standard as Exps. 1–4.
