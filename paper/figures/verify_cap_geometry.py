#!/usr/bin/env python3
"""Deterministic verification of the closed-form cap-geometry results in Sec. IV F.

Every claim Sec. IV F makes about a bounding region is analytic, so it is checked
here by quadrature rather than by sampling.  No random numbers are used and no
experiment output is read, which makes this file a self-contained provenance
record: run it and every printed number in Sec. IV F is either reproduced or the
script exits non-zero.

The point of the design is *independence*.  The manuscript derives the wide-cap
limit of a physical-speed bound by working in the normalized coordinates
:math:`\\mathbf{U}` of Eq. (10), taking the radial integral to its tail-dominated
asymptote, and cancelling :math:`v_{\\max}`.  Re-integrating that same asymptotic
expression would verify nothing but the arithmetic.  So the checks below come at
the quantities from a different direction:

* ``speed_cap_ratio_vspace`` integrates the bi-Kappa density directly in
  *velocity* space, in spherical coordinates aligned with :math:`\\hat{b}`, over
  the literal region :math:`|\\mathbf{v}|\\le v_{\\max}` at *finite*
  :math:`v_{\\max}`.  It never forms the normalized coordinates and never takes a
  limit.  Check 1 then asks whether that finite-cap sequence converges to the
  analytic limit as :math:`v_{\\max}` grows.
* ``box_ratio_vspace`` does the same for the component-wise box, in Cartesian
  velocity coordinates over the unequal-sided rectangular region the released
  code actually tests, so the exactness claim of Eq. (24) is a numerical result
  about an anisotropic integral rather than a restatement of a symmetry.
* ``ellipsoid_acceptance_quad`` integrates the Beta-prime radial density to
  check the incomplete-beta closed form of Eq. (25).

Checks
------
1. The finite-:math:`v_{\\max}` velocity-space ratio converges to
   ``speed_cap_anisotropy_limit`` for ``1/2 < kappa < 3/2``, and that limit is
   *not* unity: the bias survives an arbitrarily wide cap.
2. The limit is a wide-cap asymptote, not a value attained at every finite cap.
3. At ``kappa = 3/2`` exactly the limit is unity, and the finite-cap ratio
   approaches unity like ``1 / log v_max`` -- so the boundary case is a
   logarithmically slow bias, not a plateau.  This is the case the first draft of
   Sec. IV F got wrong by writing ``kappa <= 3/2``.
4. For ``kappa > 3/2`` the ratio converges to unity at the ordinary rate.
5. The component-wise box reproduces the scale-implied ratio exactly, at every
   finite cap and every ``kappa > 1/2``.
6. The normalized-radial ellipsoid's acceptance probability equals
   ``I_z(3/2, kappa - 1/2)`` with ``z = lambda^2 / (kappa + lambda^2)``.

``speed_cap_anisotropy_limit`` lives here rather than in
``make_manuscript_assets.py`` so that the function the manuscript's number is
printed from is the same one these checks exercise; the asset script imports it.
"""

import math
import sys

from scipy.integrate import quad
from scipy.special import betainc, beta as beta_fn


# --------------------------------------------------------------------------
# the analytic wide-cap limit quoted in Sec. IV F
# --------------------------------------------------------------------------

def speed_cap_anisotropy_limit(kappa, ratio):
    """Wide-cap limit of E[v_par^2]/E[v_x^2] under a physical-speed cap
    |v| <= v_max, expressed in units of the scale-implied (theta_par/theta_perp)^2
    so that 1 means the scale anisotropy is realized.

    A cap on |v| accepts radii up to R_max(mu) = v_max / (sqrt(kappa) *
    sqrt(theta_perp^2 (1 - mu^2) + theta_par^2 mu^2)) in the normalized
    coordinates, which depends on the polar cosine mu.  For 1/2 < kappa < 3/2 the
    second moment of the radial law diverges, the truncated moment is
    tail-dominated,

        int^{R_max} r^2 . r^2 (1 + r^2)^{-(kappa+1)} dr ~ R_max^{3-2 kappa},

    and the overall scale v_max cancels from the ratio.  What survives is a pure
    angular average, and it does not tend to 1 as the cap is relaxed.

    At kappa = 3/2 the exponent (3 - 2 kappa)/2 vanishes, the angular weight
    becomes constant, and the expression returns exactly 1 -- correctly, because
    there the divergence is logarithmic and its leading coefficient is
    direction-independent.  For kappa > 3/2 the truncated moment is core-dominated
    and the limit is likewise 1, but this expression is not the right asymptotics
    for that case and is not used there.
    """
    expo = (3.0 - 2.0 * kappa) / 2.0

    def w(mu):
        return (1.0 + (ratio ** 2 - 1.0) * mu ** 2) ** (-expo)

    num = quad(lambda mu: mu ** 2 * w(mu), 0.0, 1.0, limit=200)[0]
    den = quad(lambda mu: 0.5 * (1.0 - mu ** 2) * w(mu), 0.0, 1.0, limit=200)[0]
    return num / den


# --------------------------------------------------------------------------
# independent finite-cap quadratures, carried out in velocity space
# --------------------------------------------------------------------------

def _decade_quad(f, upper, limit=200):
    """Integrate f over [0, upper], splitting at decade boundaries.

    The radial integrands here are concentrated within a few thermal speeds of
    the origin while ``upper`` may be $10^{12}$.  A single adaptive call on such
    an interval silently misses the peak and returns a wrong answer that is not
    flagged, which is exactly the failure mode a provenance check must not have.
    Splitting into decades bounds the dynamic range each call sees.
    """
    total, lo, hi = 0.0, 0.0, 1.0
    while lo < upper:
        h = min(hi, upper)
        total += quad(f, lo, h, limit=limit)[0]
        lo, hi = h, hi * 10.0
    return total


def _bikappa_shape(vperp_sq, vpar_sq, kappa, tperp, tpar):
    """The bi-Kappa density of Eq. (2) up to its constant normalization, which
    cancels from every ratio computed here."""
    q = vperp_sq / (kappa * tperp ** 2) + vpar_sq / (kappa * tpar ** 2)
    return (1.0 + q) ** (-(kappa + 1.0))


def speed_cap_ratio_vspace(kappa, ratio, vmax, tperp=1.0):
    """E[v_par^2]/E[v_x^2] over |v| <= vmax, in units of (tpar/tperp)^2.

    Spherical coordinates in *velocity* space: v = s (sin T cos p, sin T sin p,
    cos T), volume element s^2 ds sin T dT dp.  Azimuthally the density is
    already symmetric, so <v_x^2> = <v_perp^2>/2 = s^2 sin^2(T)/2 and the phi
    integral contributes a common factor that cancels.  Finite vmax throughout;
    no asymptotics, no normalized coordinates.
    """
    tpar = ratio * tperp

    def radial(mu, weight_par):
        # weight_par: True -> mu^2 (parallel), False -> (1 - mu^2)/2 (one
        # perpendicular Cartesian component)
        ang = mu ** 2 if weight_par else 0.5 * (1.0 - mu ** 2)
        if ang == 0.0:
            return 0.0
        return ang * _decade_quad(
            lambda s: s ** 4 * _bikappa_shape(
                s ** 2 * (1.0 - mu ** 2), s ** 2 * mu ** 2, kappa, tperp, tpar),
            vmax)

    num = quad(lambda mu: radial(mu, True), 0.0, 1.0, limit=200)[0]
    den = quad(lambda mu: radial(mu, False), 0.0, 1.0, limit=200)[0]
    return (num / den) / ratio ** 2


def box_ratio_vspace(kappa, ratio, lam, tperp=1.0):
    """E[v_par^2]/E[v_x^2] over the component-wise box, in units of
    (tpar/tperp)^2.

    Integrated in Cartesian *velocity* coordinates over the literal region the
    released predicate tests, |v_x| <= lam tperp, |v_y| <= lam tperp,
    |v_z| <= lam tpar -- three unequal half-widths against an anisotropic
    density.  That the answer is exactly 1 is therefore a numerical statement
    about that integral, not a property of the quadrature grid.
    """
    tpar = ratio * tperp
    ax, az = lam * tperp, lam * tpar

    def moment(component):
        def inner(vz):
            def mid(vy):
                return quad(
                    lambda vx: (vz ** 2 if component == "par" else vx ** 2)
                    * _bikappa_shape(vx ** 2 + vy ** 2, vz ** 2,
                                     kappa, tperp, tpar),
                    0.0, ax, limit=120)[0]
            return quad(mid, 0.0, ax, limit=120)[0]
        return quad(inner, 0.0, az, limit=120)[0]

    return (moment("par") / moment("perp")) / ratio ** 2


# --------------------------------------------------------------------------
# ellipsoid acceptance
# --------------------------------------------------------------------------

def ellipsoid_acceptance_closed(kappa, lam):
    """Eq. (25): p = I_z(3/2, kappa - 1/2), z = lam^2 / (kappa + lam^2)."""
    z = lam ** 2 / (kappa + lam ** 2)
    return betainc(1.5, kappa - 0.5, z)


def ellipsoid_acceptance_quad(kappa, lam):
    """P(T <= lam^2 / kappa) by direct quadrature of the Beta-prime density
    p_T(t) = t^{1/2} (1 + t)^{-(kappa+1)} / B(3/2, kappa - 1/2)."""
    tmax = lam ** 2 / kappa
    integral = _decade_quad(
        lambda t: math.sqrt(t) * (1.0 + t) ** (-(kappa + 1.0)), tmax)
    return integral / beta_fn(1.5, kappa - 0.5)


# --------------------------------------------------------------------------
# checks
# --------------------------------------------------------------------------

FAILURES = []


def check(label, condition, detail=""):
    tag = "ok  " if condition else "FAIL"
    print(f"  [{tag}] {label}" + (f"   {detail}" if detail else ""))
    if not condition:
        FAILURES.append(label)


def close(a, b, tol):
    return abs(a - b) <= tol


def main():
    print("1. wide-cap plateau for 1/2 < kappa < 3/2 "
          "(velocity-space quadrature vs analytic limit)")
    print("   the finite-cap deficit should fall as v_max^-(3-2 kappa), so "
          "deficit * v_max^(3-2 kappa)")
    print("   is checked for constancy: that pins both the limit and the rate "
          "at which it is reached.")
    for kappa, ratio, expect in [(0.75, 2.0, 0.644623),
                                 (0.6, 2.0, None),
                                 (1.0, 2.0, None),
                                 (1.25, 2.0, None),
                                 (1.4, 2.0, None),
                                 (0.75, 0.5, None)]:
        lim = speed_cap_anisotropy_limit(kappa, ratio)
        q = 3.0 - 2.0 * kappa
        scaled = []
        for v in (1.0e4, 1.0e5, 1.0e6):
            f = speed_cap_ratio_vspace(kappa, ratio, v)
            scaled.append((lim - f) * v ** q)
        tag = f"kappa={kappa:g}, tpar/tperp={ratio:g}"
        # For theta_par < theta_perp the limit is above unity and is approached
        # from above, so the residual changes sign with the anisotropy; only its
        # magnitude and consistency of sign are asserted.
        mags = [abs(s) for s in scaled]
        check(f"{tag}: residual decays as v_max^-{q:.2f} toward the limit",
              len({s > 0.0 for s in scaled}) == 1
              and max(mags) / min(mags) < 1.25,
              f"limit={lim:.6f}, residual*v_max^q = "
              + ", ".join(f"{s:.4f}" for s in scaled))
        check(f"{tag}: the limit is not unity",
              abs(lim - 1.0) > 0.02, f"limit={lim:.6f}")
        if expect is not None:
            check(f"{tag}: limit equals the value quoted in Sec. IV F",
                  close(lim, expect, 5.0e-6), f"{lim:.6f} vs {expect:.6f}")

    print("   the plateau closes continuously as kappa -> 3/2 from below, "
          "so 1/2 < kappa < 3/2 is")
    print("   an open interval and not a discontinuity at the endpoint:")
    approach = [speed_cap_anisotropy_limit(k, 2.0)
                for k in (1.4, 1.45, 1.49, 1.499, 1.4999)]
    print("     kappa = 1.4, 1.45, 1.49, 1.499, 1.4999 -> "
          + ", ".join(f"{a:.6f}" for a in approach))
    check("the limit increases monotonically to 1 as kappa -> 3/2^-",
          all(approach[i] < approach[i + 1] for i in range(len(approach) - 1))
          and approach[-1] < 1.0 and approach[-1] > 0.999)

    print("\n2. the plateau is an asymptote, not a value attained at every cap")
    lim = speed_cap_anisotropy_limit(0.75, 2.0)
    seq = [speed_cap_ratio_vspace(0.75, 2.0, v) for v in (3.0, 10.0, 1.0e2, 1.0e4)]
    print("     v_max = 3, 10, 1e2, 1e4 -> "
          + ", ".join(f"{s:.6f}" for s in seq))
    check("narrow caps sit well below the asymptote",
          abs(seq[0] - lim) > 0.1, f"v_max=3 gives {seq[0]:.6f}, limit {lim:.6f}")
    check("the sequence increases monotonically toward the asymptote",
          all(seq[i] < seq[i + 1] for i in range(len(seq) - 1))
          and seq[-1] < lim + 1.0e-5)

    print("\n3. kappa = 3/2 exactly: unity, approached like 1/log(v_max)")
    for ratio in (0.5, 2.0, 4.0):
        check(f"analytic limit is exactly 1 at kappa=3/2, tpar/tperp={ratio:g}",
              close(speed_cap_anisotropy_limit(1.5, ratio), 1.0, 1.0e-9))
    vs = [1.0e4, 1.0e6, 1.0e8, 1.0e12]
    rs = [speed_cap_ratio_vspace(1.5, 2.0, v) for v in vs]
    for v, r in zip(vs, rs):
        print(f"     v_max={v:>8.0e} -> {r:.6f}   "
              f"(1-r)*log(v_max) = {(1.0 - r) * math.log(v):.4f}")
    check("the finite-cap ratio increases toward 1",
          all(rs[i] < rs[i + 1] for i in range(len(rs) - 1)) and rs[-1] > 0.98)
    prods = [(1.0 - r) * math.log(v) for v, r in zip(vs, rs)]
    check("the deficit scales as 1/log(v_max), i.e. no plateau",
          max(prods) / min(prods) < 1.15,
          f"(1-r)log(v_max) spans {min(prods):.4f}-{max(prods):.4f}")

    print("\n4. kappa > 3/2: unity at the ordinary rate")
    for kappa in (2.0, 5.0):
        r = speed_cap_ratio_vspace(kappa, 2.0, 1.0e4)
        check(f"kappa={kappa:g}: ratio at v_max=1e4 is 1", close(r, 1.0, 1.0e-3),
              f"{r:.6f}")

    print("\n5. the component-wise box realizes the scale-implied ratio exactly")
    for kappa, ratio, lam in [(0.75, 2.0, 3.0), (0.75, 2.0, 20.0),
                              (1.5, 4.0, 5.0), (3.0, 0.5, 10.0), (0.55, 2.0, 8.0)]:
        r = box_ratio_vspace(kappa, ratio, lam)
        check(f"kappa={kappa:g}, tpar/tperp={ratio:g}, lambda={lam:g}",
              close(r, 1.0, 1.0e-6), f"ratio/(tpar/tperp)^2 = {r:.9f}")

    print("\n6. ellipsoid acceptance: closed form vs direct quadrature")
    for kappa, lam in [(0.75, 3.0), (0.75, 50.0), (1.5, 10.0),
                       (2.0, 5.0), (10.0, 3.0), (0.51, 20.0)]:
        a = ellipsoid_acceptance_closed(kappa, lam)
        b = ellipsoid_acceptance_quad(kappa, lam)
        check(f"kappa={kappa:g}, lambda={lam:g}", close(a, b, 1.0e-9),
              f"closed={a:.9f} quad={b:.9f}")

    print()
    if FAILURES:
        print(f"{len(FAILURES)} check(s) FAILED:")
        for f in FAILURES:
            print(f"  - {f}")
        return 1
    print("all cap-geometry checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
