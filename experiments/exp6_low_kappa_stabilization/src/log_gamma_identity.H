// LOG-ID: small-shape Gamma variate generated directly on the logarithmic scale via the
// shape-boosting identity.
//
// Attribution.  For 0 < a < 1,
//
//     X ~ Gamma(a, 1)   <=>   X = Y U^{1/a},   Y ~ Gamma(a+1, 1),  U ~ Uniform(0,1),
//
// is the classical shape-boosting relation of Ahrens and Dieter (1974); it is the route
// by which Marsaglia-Tsang style generators - and hence both libc++ and libstdc++ - reach
// shapes below one.  Writing it on the logarithmic scale,
//
//     log X = log Y + (log U) / a,
//
// is likewise established rather than new: Xi, Tan and Liu (2013, J. Stat. Softw. 55(4),
// doi:10.18637/jss.v055.i04) return log-scale output from their Gamma generators at small
// shape, and Liu, Martin and Syring (2017, Comput. Stat. 32, 1767-1775,
// doi:10.1007/s00180-016-0692-0) generate small-shape Gamma variates on the log scale for
// exactly this reason.  Nothing in this file is claimed as a new sampling method.  What is
// under test here is its integration into a complete anisotropic bi-Kappa loader.
//
// What the identity buys.  The product form materializes U^{1/a}, which for small a
// underflows to exactly zero long before the mathematical variate is negligible: at
// a = 0.01, U^{1/a} = U^{100}, so U < 0.000724 already drives a double to zero.  The sum
// form has no such loss: log Y is order unity and (log U)/a is a moderate negative number.
//
// What it does not buy.  log X is a stable representation of X, not an exact one: both
// log Y and (log U)/a are rounded, so log X carries roughly one ulp of the larger term.
// At a = 0.01 that is an absolute error near 1e-14 in the log, i.e. a relative error near
// 1e-14 in X.  The experiment measures this rather than asserting it away.

#ifndef EXP6_LOG_GAMMA_IDENTITY_H
#define EXP6_LOG_GAMMA_IDENTITY_H

#include <cmath>
#include <cstdint>
#include <limits>
#include <random>

namespace exp6
{

/// One small-shape Gamma draw, reported on every scale the experiment needs to compare.
template <typename T> struct LogGammaDraw
{
    T log_x;      ///< log of the variate; finite whenever the draw exists
    T x_product;  ///< what the product form Y U^{1/a} would have materialized (may be 0)
    T y;          ///< the Gamma(a+1,1) variate actually drawn
    T u;          ///< the Uniform(0,1) variate actually drawn
    uint32_t retries;   ///< open-interval endpoint redraws consumed
    uint32_t attempts;  ///< acceptance-rejection attempts (1 for a rejection-free method)
    uint32_t uniforms;  ///< uniform variates consumed, including redraws
    uint32_t gammas;    ///< Gamma variates consumed
    bool out_of_domain; ///< the primitive's published domain did not cover this shape
    bool ok;            ///< false only if log_x came out non-finite
};

/// LOG-ID.  `gBoost` must already be parameterised Gamma(a+1, 1) and `unif` Uniform(0,1).
///
/// The uniform must land in the open interval (0,1) because its logarithm is taken.
/// std::uniform_real_distribution yields [0,1), so the single reachable bad endpoint is
/// exactly 0; it is redrawn, and the redraw is counted, because the redraw rule is part of
/// the implemented discrete sampler and not an implementation detail.  The endpoint is
/// never clamped to the smallest positive value: clamping would silently replace the
/// tail of the sampled law with an atom.
template <typename T, typename Gen>
LogGammaDraw<T> drawLogGammaIdentity(T a, Gen &gen, std::gamma_distribution<T> &gBoost,
                                     std::uniform_real_distribution<T> &unif)
{
    LogGammaDraw<T> d;
    d.retries= 0;
    d.attempts= 1;
    d.uniforms= 1;
    d.gammas= 1;
    d.out_of_domain= false;

    d.y= gBoost(gen);
    T u= unif(gen);
    while (!(u > T(0)))
    {
        ++d.retries;
        ++d.uniforms;
        u= unif(gen);
    }
    d.u= u;

    d.log_x= std::log(d.y) + std::log(u) / a;
    d.x_product= d.y * std::pow(u, T(1) / a);
    d.ok= (d.log_x == d.log_x) && (d.log_x != std::numeric_limits<T>::infinity()) &&
          (d.log_x != -std::numeric_limits<T>::infinity());
    return d;
}

/// The same identity evaluated in double from the recorded primitives.  This is the
/// "intended" log X2 that the mechanism decomposition compares a working-precision result
/// against; the independent arbitrary-precision oracle re-derives it again from the same
/// primitives at 100 decimal digits.
inline double logGammaIdentityReference(double a, double y, double u)
{
    return std::log(y) + std::log(u) / a;
}

} // namespace exp6

#endif // EXP6_LOG_GAMMA_IDENTITY_H
