"""The representability boundary of a loader cell, and the two F5 upper-tail members.

Why this module exists
----------------------
Amendment 1.3.0 gave the F5 loader battery an upper-tail statistic, because a battery
certifying a heavy-tailed law had nothing in it that looked at the tail.  Both of its
members were stated against the *untruncated* law:

  * the COUNT above ``z0 = -log q0`` against ``Binomial(n_attempted, q0)``, and
  * a Kolmogorov-Smirnov test of the excesses ``Z - z0`` against ``Exp(1)``.

Neither null is the one the data obey, and the first holdout failed a correct candidate on
both.  The reason is the same in each case, and it is not a defect of the statistics: near
``kappa = 1/2`` the bi-Kappa law puts non-zero probability outside every finite
floating-point range, so the loader cannot return the far tail.  It is *right* not to.

The excess member was the worse of the two.  Honest overflow removes exactly the largest
draws, so the surviving excesses are right-censored, and a Kolmogorov-Smirnov test against an
uncensored ``Exp(1)`` rejects with probability approaching one.  Measured over 2000
replicates of a perfectly correct loader at the frozen production size, the frozen member's
rejection rate was **1.000** on cases C3 and C4 against a nominal 0.010.  The count member
had the same problem in a milder form, and amendment 1.3.0 had already patched around it
with a side condition -- ``count_member_requires_q0_above_honest_floor`` -- that declared the
member "not applicable" wherever it would have misfired.

The censoring is also **direction-dependent**, which is what makes a single cutoff
insufficient.  A draw is returned iff every component is representable,

    R * M(n) <= max(),     M(n) = max_j |g_j|,   g = Q(ub) diag(s_perp, s_perp, s_par) n,

with ``s = sqrt(kappa) theta``.  ``M`` varies over the sphere by up to
``sqrt(3) theta_max / theta_min``, so the tail probability at which a draw stops being
returnable is a property of its own direction, not of the cell.  A cutoff inferred from the
cell's total overflow rate -- ``z_f = -log f`` -- is the average of that boundary, not the
boundary.

What replaces them
------------------
**Member 1, the count.**  The probe records ``log_r_ref``, the radius the attempt carried,
for *every* attempt, the overflowed ones included (``src/exp7_probe.cpp``, pushed
unconditionally; an uncapped cell runs its core mapping once per attempt, so its record count
equals its attempt count).  The intended ``Z`` of every attempt is therefore already on disk,
uncensored.  Against that sample the count above ``z0`` is ``Binomial(n_attempted, q0)``
**exactly**, at every threshold and in every uncapped cell.  The frozen null is recovered as
it was written, and the side condition it needed is no longer needed: the observability
problem is solved by data that was always there rather than by declining to test.

**Member 2, the excess.**  The count member reads the law the loader *intended*.  The
fidelity claim is about the population it *returned*, so the excess member stays on the
returned draws -- and is given the null those draws actually obey.  ``Z`` is independent of
the direction, so conditional on its own direction a returned draw above ``z0`` is
``Exp(1)`` truncated to ``(z0, C_i]``, where

    C_i = Z(log max() - log M(n_i))

is that draw's own representability threshold, computed from its own recovered direction.
Its probability integral transform

    U_i = (1 - exp(-(Z_i - z0))) / (1 - exp(-(C_i - z0)))

is i.i.d. ``Uniform(0,1)`` under the null, and the member is a Kolmogorov-Smirnov test of it.

This is not a different test.  Where no attempt can overflow, ``C_i`` is effectively infinite,
``U_i`` reduces to ``1 - exp(-(Z_i - z0))``, and because a Kolmogorov-Smirnov statistic is
invariant under a common monotone transform of the data and the null CDF, the number it
returns is *identical* to the frozen member's.  Measured on case C0, the two agree to
1.1e-16.  The replacement equals the frozen test wherever the frozen test was valid, and is
defined where it was not.  It is also the construction the protocol already uses for the
capped cells, whose accepted radius is truncated at a direction-dependent bound in exactly
the same way (``loader_tests``' cap-law member).

Neither member needs quadrature, a grid, or any numerical integration: both are per-draw
arithmetic in the frozen ``Z`` transform, whose evaluator gate G0 validates against an
arbitrary-precision incomplete beta.
"""

from __future__ import annotations

import math

import numpy as np

import exp7_io as IO
import exp7_stats as S

#: The largest finite value of each precision the experiment runs in.
MAX_FINITE = {"double": float(np.finfo(np.float64).max),
              "float": float(np.finfo(np.float32).max)}


def g_from_direction(n_hat, kappa: float, theta_perp: float, theta_par: float, ub):
    """``g`` with ``V = R g``, for unit directions in the field-aligned frame.

    The scaling and the rotation are the released header's, in its order:
    ``g = Q(ub) diag(sqrt(kappa) theta_perp, sqrt(kappa) theta_perp, sqrt(kappa) theta_par) n``.
    ``IO.field_basis`` re-derives ``Q`` independently of the C++ rather than importing it, so
    that a frame error would show up here rather than cancel.
    """
    n_hat = np.atleast_2d(np.asarray(n_hat, dtype=float))
    sk = math.sqrt(float(kappa))
    scale = np.array([sk * theta_perp, sk * theta_perp, sk * theta_par], dtype=float)
    return (n_hat * scale) @ IO.field_basis(ub).T


def max_component(n_hat, kappa: float, theta_perp: float, theta_par: float, ub):
    """``M(n) = max_j |g_j|`` -- the component that decides representability."""
    return np.max(np.abs(g_from_direction(n_hat, kappa, theta_perp, theta_par, ub)), axis=1)


def censoring_threshold_z(n_hat, kappa: float, theta_perp: float, theta_par: float, ub,
                          precision: str):
    """``C(n)``: the ``Z`` at which this direction's largest component reaches ``max()``.

    Returned in the same ``Z`` units as the data, through the same frozen evaluator, so that
    any bias in the transform cancels between the observation and its boundary instead of
    being compared across two different routes.
    """
    a = float(kappa) - 0.5
    m = max_component(n_hat, kappa, theta_perp, theta_par, ub)
    with np.errstate(divide="ignore"):
        log_rc = math.log(MAX_FINITE[precision]) - np.log(m)
    return S.z_from_log_w(S.log_w_from_log_r(log_rc), a)


def recovery_tolerance(log_r, precision: str) -> np.ndarray:
    """Bound on ``|recovered log R - log_r_ref|`` for a returned draw.

    The log representation carries an absolute error of order ``eps |log R|`` (the header
    documents this), the recovery adds an exponential, a rotation and a logarithm, and eight
    times the dominant term covers all of it.  It is derived from that error model, not fitted
    to what came out: the observed worst cases are 1.7e-6 in ``float`` and 5.7e-14 in
    ``double``, against bounds of 8.6e-5 and 1.3e-12.
    """
    eps = float(np.finfo(np.float32 if precision == "float" else np.float64).eps)
    return 8.0 * eps * (np.abs(np.asarray(log_r, dtype=float)) + 2.0)
