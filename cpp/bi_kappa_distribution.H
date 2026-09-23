#ifndef _BI_KAPPA_DISTRIBUTION_H_
#define _BI_KAPPA_DISTRIBUTION_H_

/// Release version of this header.  A consumer that vendors the file can test these to tell
/// which sampler it has; 2.x is the first line whose radius is built in the log domain, and
/// its random stream differs from 1.x for every seed.
///
/// 2.1.0 supersedes 2.0.0, which was never released.  The two differ in which vectors they
/// return at the representability boundary and in what `n_nonfinite()` counts under a cap,
/// so they are two samplers and are versioned as two; see CHANGELOG.md.
///
/// 2.2.0 carries the log-domain radius of a `float` instantiation in `double` and rounds
/// once, when the component is materialized.  It returns the same numbers as 2.1.0 in
/// `double`, bit for bit, and more accurate ones in `float`; see
/// @ref bikappa_detail::log_accumulator.
#define BI_KAPPA_VERSION_MAJOR 2
#define BI_KAPPA_VERSION_MINOR 2
#define BI_KAPPA_VERSION_PATCH 1
#define BI_KAPPA_VERSION_STRING "2.2.1"

#include <array>
#include <algorithm>
#include <cmath>
#include <limits>
#include <random>
#include <stdexcept>

// This sampler decides what it may return by materializing a value and asking whether it is
// finite.  `-ffast-math` (and `-ffinite-math-only` alone) licenses a compiler to assume no
// infinity or NaN ever occurs, and so to fold `mag <= max()` to `true` -- which would return
// an overflowed component as if it were a number AND leave n_nonfinite() reading zero: the
// exact failure this file is built to prevent, with no symptom.  Refusing to compile is the
// only safe response, because there is nothing to check at run time once the check is gone.
// Define BI_KAPPA_ALLOW_FAST_MATH to override, and do not trust the counters if you do.
#if defined(__FAST_MATH__) && !defined(BI_KAPPA_ALLOW_FAST_MATH)
#error "bi_kappa_distribution.H requires IEEE-754 infinity semantics: -ffast-math / \
-ffinite-math-only would silently disable its representability test. Define \
BI_KAPPA_ALLOW_FAST_MATH to compile anyway (the overflow counters then mean nothing)."
#endif

/**
 * @brief Sampling primitives used to build the bi-Kappa radius.
 *
 * These live here, rather than being delegated to `<random>`, because the radius is a Gamma
 * ratio whose denominator shape goes to zero as kappa -> 1/2.  In that limit the denominator
 * variate underflows to exactly zero and the draw is lost, so the loader carries `log X2`
 * and never materializes `X2`.  Taking a logarithm forces three requirements that `<random>`
 * does not meet:
 *
 *  - the uniform must lie in the *open* interval, while `std::uniform_real_distribution`
 *    yields `[0,1)` -- and, in one of the two standard libraries this loader is tested
 *    against, can yield exactly `1` for `float`, whose logarithm is zero;
 *  - `std::gamma_distribution` uses a different algorithm in each standard library, which
 *    would make the stream, the variate count and the finite-precision failure envelope
 *    properties of the library rather than of the method; and
 *  - `std::normal_distribution` caches the second variate of each generated pair, and the
 *    two libraries return that pair in opposite orders.
 *
 * Fixing all three here makes the sampled sequence a function of the engine alone.
 */
namespace bikappa_detail
{

/**
 * @brief The type in which the log-domain radius is accumulated.
 *
 * `log R` is a stable representation of the radius, not an exact one: its absolute error
 * is about `eps |log R|`, so the *relative* error of the radius it stands for is also about
 * `eps |log R|`.  Near the top of a type's range `|log R|` is the logarithm of the largest
 * finite number -- 88.7 in `float`, 709.8 in `double` -- and the shape-boosting identity
 * divides one of its two terms by `a = kappa - 1/2`, which amplifies the error of `log U` by
 * `1/a`: 200-fold at `kappa = 0.505`.  In `float` the result is an error of a few parts in
 * a million on a quantity whose scale is 88.7, which is larger than the gap between
 * neighbouring `float`s at the overflow threshold.  The radius is then not resolved finely
 * enough to say which side of the threshold a component near it falls on, and a draw that
 * is representable can be exponentiated one ulp high and lost.  One such draw appeared in
 * 14.2 million audited attempts at `float kappa = 0.505`: its largest component lay
 * 3.1e-08 natural-log units -- about half an ulp -- below `FLT_MAX`, and the `float` log
 * radius carried 3.5e-06 of error, a hundred times that margin.
 *
 * Widening only the accumulator fixes it.  The variates are still drawn in `RealType`, from
 * the same engine bits, in the same order, so the sampled sequence does not change; what
 * changes is that the sum `log Y + (log U)/a` and the radius built from it are formed in a
 * type whose resolution is far finer than the margin being decided, and the value is rounded
 * to `RealType` exactly once, where the component is materialized.  For `float` that moves
 * the resolution of `log R` from about 5.3e-06 to about 1.6e-14 and leaves a single rounding
 * between the mathematics and the returned number.
 *
 * For `double` (and for `long double`, where the two coincide on the platforms this header
 * targets) the accumulator *is* the working type, so the arithmetic is unchanged and a
 * `double` instantiation of 2.2.0 returns the same bits as 2.1.0.
 */
template <typename RealType> struct log_accumulator
{
    typedef RealType type;
};

/// `float` accumulates in `double`: 53 bits against 24, over a range that easily holds
/// `log(FLT_MAX)` and the largest `float` alike.
template <> struct log_accumulator<float>
{
    typedef double type;
};

/// Whole bits of entropy the engine delivers per call (32 for `std::mt19937`).
template <typename GeneratorType> inline int engineBits()
{
    unsigned long long range= static_cast<unsigned long long>(GeneratorType::max()) -
                              static_cast<unsigned long long>(GeneratorType::min());
    int bits= 0;
    while (range >= 1ULL)
    {
        ++bits;
        range>>= 1;
    }
    return bits;
}

/**
 * @brief Uniform on the **open** interval (0,1), by an explicit mapping from engine bits.
 *
 * `b = digits - 1` random bits give an integer `k` uniform on `[0, 2^b)`, and the value
 * returned is `(k + 1/2) / 2^b`.  Spending one bit fewer than the significand is what makes
 * the arithmetic exact: `k + 1/2` is then exactly representable and the division by a power
 * of two is exact, so the result lies strictly inside (0,1) for every `k` -- the smallest is
 * `2^-(b+1)`, the largest `1 - 2^-(b+1)`.  No endpoint is ever produced, so none is redrawn
 * and none is clamped.  Clamping is not a neutral choice: it replaces the end of the sampled
 * law with an atom, which is what one standard library does for the value 1.
 *
 * The bit costs one part in `2^53` (`double`) or `2^24` (`float`) of resolution, far below
 * what `log(U)/a` can carry at any shape.
 */
template <typename RealType, typename GeneratorType> RealType openCanonical(GeneratorType &gen)
{
    const int bits= std::numeric_limits<RealType>::digits - 1;
    const int per= engineBits<GeneratorType>();
    unsigned long long k= 0;
    int have= 0;
    while (have < bits)
    {
        const unsigned long long word= static_cast<unsigned long long>(gen()) -
                                       static_cast<unsigned long long>(GeneratorType::min());
        const int take= (bits - have < per) ? (bits - have) : per;
        k= (k << take) | (word & ((1ULL << take) - 1ULL));
        have+= take;
    }
    return (static_cast<RealType>(k) + RealType(0.5)) * std::ldexp(RealType(1.0), -bits);
}

/// Standard normal by the Marsaglia polar method.  The second variate of each generated pair
/// is cached, so `cached` must be cleared whenever the engine is reseeded -- otherwise a
/// reseeded stream opens with a variate belonging to the previous one.
template <typename RealType, typename GeneratorType>
RealType polarNormal(GeneratorType &gen, RealType &cache, bool &cached)
{
    if (cached)
    {
        cached= false;
        return cache;
    }
    RealType x, y, s;
    do
    {
        x= RealType(2.0) * openCanonical<RealType, GeneratorType>(gen) - RealType(1.0);
        y= RealType(2.0) * openCanonical<RealType, GeneratorType>(gen) - RealType(1.0);
        s= x * x + y * y;
    } while (s >= RealType(1.0));
    // s == 0 is unreachable: openCanonical never returns exactly 1/2, so neither x nor y is
    // ever exactly zero.
    const RealType f= std::sqrt(RealType(-2.0) * std::log(s) / s);
    cache= y * f;
    cached= true;
    return x * f;
}

/// Gamma(shape, 1) for shape >= 1, by Marsaglia and Tsang (2000), ACM TOMS 26(3), 363-372.
/// Expected attempts stay below 1.01 for every shape >= 1, so the cost does not depend on
/// the shape -- which is why the boosted shape below is always taken above one.
template <typename RealType, typename GeneratorType>
RealType marsagliaTsangGamma(RealType shape, GeneratorType &gen, RealType &ncache, bool &ncached)
{
    const RealType d= shape - RealType(1.0) / RealType(3.0);
    const RealType c= RealType(1.0) / std::sqrt(RealType(9.0) * d);
    for (;;)
    {
        RealType x, v;
        do
        {
            x= polarNormal<RealType, GeneratorType>(gen, ncache, ncached);
            v= RealType(1.0) + c * x;
        } while (v <= RealType(0.0));
        v= v * v * v;
        const RealType u= openCanonical<RealType, GeneratorType>(gen);
        const RealType xx= x * x;
        if (u < RealType(1.0) - RealType(0.0331) * xx * xx)
        {
            return d * v;
        }
        if (std::log(u) < RealType(0.5) * xx + d * (RealType(1.0) - v + std::log(v)))
        {
            return d * v;
        }
    }
}

/**
 * @brief `log X` for `X ~ Gamma(a, 1)`, for every `a > 0`, without ever forming `X`.
 *
 * Uses the shape-boosting relation of Ahrens and Dieter (1974),
 * `X = Y U^{1/a}` with `Y ~ Gamma(a+1, 1)` and `U ~ Uniform(0,1)`, written additively:
 *
 *     log X = log Y + (log U) / a.
 *
 * The product form is what underflows: at `a = 0.01`, `U^{1/a} = U^{100}` is already zero in
 * `double` for `U < 7.2e-4`.  The sum has no such loss, because `log Y` is order unity and
 * `(log U)/a` is a moderate negative number.  The relation holds for every positive shape, so
 * this is one code path with no domain split and no fallback, and it is rejection-free at the
 * boosting step: exactly one Gamma and one uniform per call.
 *
 * `log X` is a *stable* representation of `X`, not an exact one.  `log Y` carries about one
 * unit in the last place, and `(log U)/a` about one unit in the last place of `|log U|/a`, so
 * the absolute error of `log X` grows like `eps * |log X|` -- equivalently the relative error
 * of the radius grows like `eps * |log R|`, losing roughly one decimal digit per decade of
 * dynamic range.  The uniform also has finite resolution, so `|log U|` cannot exceed about
 * 36.7 (`double`) or 16.6 (`float`), and the sampled law is truncated below `exp(-36.7/a)`.
 * That region carries about `1.1e-16` of the probability and lies far outside the range in
 * which any velocity component is representable.
 */
template <typename RealType, typename GeneratorType>
typename log_accumulator<RealType>::type logGamma(RealType a, GeneratorType &gen,
                                                  RealType &ncache, bool &ncached)
{
    typedef typename log_accumulator<RealType>::type AccumType;
    const RealType y=
        marsagliaTsangGamma<RealType, GeneratorType>(a + RealType(1.0), gen, ncache, ncached);
    const RealType u= openCanonical<RealType, GeneratorType>(gen);
    // Y, U and the shape are the working-precision values the sampler actually drew and the
    // caller actually asked for; only the logarithms and the division are widened.  Widening
    // the shape instead of the variates would change which distribution is sampled.
    return std::log(static_cast<AccumType>(y)) +
           std::log(static_cast<AccumType>(u)) / static_cast<AccumType>(a);
}

/**
 * @brief Materialize `V_j = R g_j` from `log R`, and report whether any component overflowed.
 *
 * This is the one place the loader decides what is representable, and it decides it from the
 * value it is about to return rather than from that value's logarithm.
 *
 * The obvious test -- compare `log|V_j| = log R + log|g_j|` against `log(max())` -- is wrong
 * by up to a rounding step, in both directions, and neither error is academic:
 *
 *  - `log(max())` is itself a rounded value and can land *above* the true logarithm, so
 *    `exp` of it need not be finite.  In `float`, `exp(log(FLT_MAX))` is exactly `+inf`.  A
 *    draw whose `log R` equals `log(FLT_MAX)` bit for bit therefore passes a
 *    `log R <= logMax` test, and `exp(log R)` then overflows -- while every component it
 *    should have produced is perfectly representable, because `|g_j| < 1` pulls each of them
 *    back under the limit.  A confirmatory run found exactly such a draw in 14.2 million
 *    audited attempts, at `float kappa = 0.51`: `log R` was `0x1.62e43p+6`, identical to
 *    `log(FLT_MAX)`, the largest component's logarithm was 88.151 against a limit of 88.723,
 *    and the loader returned `(-inf, +inf, +inf)` -- and, the overflow never having been
 *    detected, counted nothing.
 *  - The same rounding in the other direction makes a logarithmic test discard a draw that
 *    *is* representable, which is an avoidable loss rather than an honest one.
 *
 * Comparing the materialized magnitude against `max()` cannot have either defect: it is the
 * number that will be returned, so the test and the result agree by construction.  It is
 * also cheaper, because the common path takes no logarithm at all -- one per component is
 * needed only where `R` alone overflows, which is where the log domain is doing the work
 * anyway and where there is no accuracy left to lose.
 *
 * @param logR   logarithm of the radius, in the accumulator type of
 *               @ref bikappa_detail::log_accumulator.
 * @param g      the order-unity vector, already rotated into the frame the caller returns,
 *               in the same accumulator type.  `R` and `g_j` are both order-unity-to-moderate
 *               quantities whose product spans the whole range of the returned type, so the
 *               product is where the one rounding belongs: the working-precision `g_j` of a
 *               `float` run carries about 6e-8 of relative error, which is larger than the
 *               margin by which a component near `FLT_MAX` is representable, and rounding it
 *               before the multiply throws away a draw the type can hold.
 * @param point  receives the materialized components, `+/-inf` where one has no
 *               representation, with the sign the mathematics gives it.
 * @return       true iff at least one component had no representation in @c RealType.
 */
template <typename RealType>
bool materializeComponents(typename log_accumulator<RealType>::type logR,
                           const std::array<typename log_accumulator<RealType>::type, 3> &g,
                           std::array<RealType, 3> &point)
{
    typedef typename log_accumulator<RealType>::type AccumType;
    const AccumType accumMax= std::numeric_limits<AccumType>::max();
    const RealType maxFinite= std::numeric_limits<RealType>::max();
    const AccumType rTry= std::exp(logR);
    // `false` for +infinity and for NaN, both of which fail every comparison.  The test is
    // in the accumulator, because that is where the product below is formed; a `float` run
    // reaches this branch only for a radius beyond `exp(709.78)`, far above anything a
    // `float` component could hold.
    const bool radiusRepresentable= (rTry <= accumMax);
    const AccumType R= radiusRepresentable ? rTry : AccumType(0.0);

    bool anyOverflow= false;
    for (int j= 0; j < 3; ++j)
    {
        if (g[j] == AccumType(0.0))
        {
            point[j]= RealType(0.0);
            continue;
        }
        const AccumType absG= std::abs(g[j]);
        const AccumType wide= radiusRepresentable ? R * absG
                                                  : std::exp(logR + std::log(absG));
        // The single rounding.  What decides representability is the number the caller will
        // receive, which is `wide` rounded to `RealType` -- so a magnitude that lies within
        // half an ulp of `max()` rounds to `max()` and is returned, and one that rounds past
        // it becomes an infinity here and is counted.  Deciding from `wide` instead would
        // discard a draw the type can hold; deciding from a logarithm would do it in both
        // directions at once.
        const RealType mag= static_cast<RealType>(wide);
        if (!(mag <= maxFinite))
        {
            anyOverflow= true;
            point[j]= (g[j] < AccumType(0.0)) ? -std::numeric_limits<RealType>::infinity()
                                              : std::numeric_limits<RealType>::infinity();
            continue;
        }
        point[j]= (g[j] < AccumType(0.0)) ? -mag : mag;
    }
    return anyOverflow;
}

} // namespace bikappa_detail

/**
 * @brief Bi-Kappa distribution class for sampling 3D velocity vectors.
 *
 * This class follows the interface pattern of standard library distributions
 * like std::uniform_real_distribution<T>.
 *
 * @code{.cpp}
 * bi_kappa_distribution<double> dist;
 * dist.define(2.0, 1.0, 2.0, {0, 0, 1});   // kappa, theta_perp, theta_par, ub
 * auto v = dist();                          // {vx, vy, vz}, global frame
 * @endcode
 *
 * @note **The cap is off by default**, so the default target law is the bi-Kappa
 * distribution itself.  Passing a finite `max_normalized_velocity` opts in to a
 * *different*, bounded target: the bi-Kappa law conditioned on the component-wise
 * box.  See no_cap() and param_type::capped().
 *
 * @note **The radius is built in the logarithmic domain.**  `R = sqrt(X1)/sqrt(X2)` with
 * `X2 ~ Gamma(kappa - 1/2, 1)`, so as `kappa -> 1/2` the denominator's shape goes to zero
 * and `X2` underflows to exactly zero.  `X2` is therefore never formed: the sampler carries
 * `log X2` and `log R`, and reaches the linear domain only at the last step, once per
 * component.  Whether a component is representable is then decided from the component
 * itself -- the number about to be returned -- and never from its logarithm, because
 * `log(max())` is a rounded value that can sit on either side of the true limit; see
 * @ref bikappa_detail::materializeComponents.  Two consequences a caller should know about:
 *
 *  - **The random stream differs from releases before 2.0.0.**  The Gamma, normal and
 *    uniform primitives are implemented in @ref bikappa_detail rather than taken from
 *    `<random>`, and the direction is drawn by a rejection method that uses no trigonometry,
 *    so a given seed produces a different sequence than it used to.  In exchange the sequence
 *    no longer depends on which standard library is in use: the same seed gives bit-identical
 *    output under libc++ and libstdc++, which it did not before -- `std::gamma_distribution`
 *    uses a different algorithm in each, and a `sin`/`cos` pair on one argument is fused by
 *    one compiler and not the other.  What the sequence does still depend on is whether the
 *    compiler contracts a multiply and an add; pass `-ffp-contract=off` to pin it to what the
 *    source text says.  Across architectures nothing is bitwise: `log` and `exp` are not
 *    correctly rounded and their implementations differ.
 *  - **Non-finite output is still possible, and is now countable.**  Sufficiently near
 *    `kappa = 1/2` the bi-Kappa law puts non-zero probability outside every finite
 *    floating-point range; those draws exist as mathematics and cannot be returned as
 *    numbers.  Untruncated mode returns them with non-finite components, as before, and
 *    n_nonfinite() counts them.  What has changed is how many there are: the loss is now the
 *    representability floor of @c RealType rather than the underflow limit of an
 *    intermediate.  Measured over 2x10^6 draws at `kappa = 0.505` in `double`, the rate falls
 *    from 2.4e-2 to 8.4e-4, which matches the analytic floor.
 *
 * @see \ref usage for more recipes (seeding, external RNG, reconfiguration).
 *
 * @tparam RealType The floating-point type (e.g., double, float)
 */
template <typename RealType= double> class bi_kappa_distribution
{
  public:
    using result_type= std::array<RealType, 3>;
    using point_type= result_type;

    /// The type in which the sampler forms `R g` before rounding once; see
    /// @ref bikappa_detail::log_accumulator.  Identical to @c point_type unless
    /// @c RealType is @c float.  Internal: no public entry point takes or returns it.
    using accum_point=
        std::array<typename bikappa_detail::log_accumulator<RealType>::type, 3>;

    /**
     * @brief Parameter class to store distribution parameters
     */
    class param_type
    {
      public:
        using distribution_type= bi_kappa_distribution<RealType>;

        // The cap defaults to +infinity, i.e. off -- the same value no_cap() returns.
        // It is spelled out here rather than as no_cap() because the enclosing class is
        // still incomplete while this nested default argument is parsed.
        param_type(RealType kappa= RealType(2.0), RealType theta_perp= RealType(1.0),
                   RealType theta_par= RealType(1.0),
                   point_type ub= point_type({RealType(0.0), RealType(0.0), RealType(1.0)}),
                   RealType max_normalized_velocity=
                       std::numeric_limits<RealType>::infinity())
            : m_kappa(kappa), m_theta_perp(theta_perp), m_theta_par(theta_par), m_ub(ub),
              m_maxNormalizedVelocity(max_normalized_velocity)
        {
            // Each test is written so that NaN fails it.  `x <= 0.5` is false for a NaN,
            // so the natural spelling admits one, and a NaN parameter then propagates into
            // log R, where it can no longer be told from an unrepresentable velocity: the
            // sampler would report a non-finite draw and count it as honest overflow, which
            // it is not.  Rejecting the parameter is the only point at which the two are
            // still distinguishable.
            if (!(m_kappa > RealType(0.5)))
            {
                throw std::invalid_argument("kappa must be > 0.5 for normalization");
            }
            if (!(m_theta_perp > RealType(0.0)) || !(m_theta_par > RealType(0.0)) ||
                !(m_theta_perp <= std::numeric_limits<RealType>::max()) ||
                !(m_theta_par <= std::numeric_limits<RealType>::max()))
            {
                throw std::invalid_argument(
                    "Temperature parameters must be positive and finite");
            }
            {
                // The field direction is validated here rather than inside the rotation, so
                // that an axis-aligned-looking NaN cannot slip past the fast path.
                RealType n2= RealType(0.0);
                for (int j= 0; j < 3; ++j)
                {
                    if (!(std::abs(m_ub[j]) <= std::numeric_limits<RealType>::max()))
                    {
                        throw std::invalid_argument("ub components must be finite");
                    }
                    n2+= m_ub[j] * m_ub[j];
                }
                if (!(n2 > RealType(0.0)))
                {
                    throw std::invalid_argument("ub must be a non-zero vector");
                }
            }
            // A positive finite value selects the capped mode.  no_cap() (= +inf)
            // selects the untruncated mode, in which no acceptance-rejection loop
            // runs at all and samples follow the full bi-Kappa law.
            if (!(m_maxNormalizedVelocity > RealType(0.0)))
            {
                throw std::invalid_argument(
                    "max_normalized_velocity must be positive, or no_cap() for the "
                    "untruncated distribution");
            }
        }

        RealType kappa() const { return m_kappa; }
        RealType theta_perp() const { return m_theta_perp; }
        RealType theta_par() const { return m_theta_par; }
        point_type ub() const { return m_ub; }
        RealType max_normalized_velocity() const { return m_maxNormalizedVelocity; }

        /// True when a finite component-wise cap is in force, i.e. when the
        /// sampled law is the bi-Kappa distribution conditioned on the cap event
        /// rather than the bi-Kappa distribution itself.
        bool capped() const { return std::isfinite(m_maxNormalizedVelocity); }

        bool operator==(const param_type &other) const
        {
            return m_kappa == other.m_kappa && m_theta_perp == other.m_theta_perp &&
                   m_theta_par == other.m_theta_par && m_ub == other.m_ub &&
                   m_maxNormalizedVelocity == other.m_maxNormalizedVelocity;
        }

        bool operator!=(const param_type &other) const { return !(*this == other); }

      private:
        RealType m_kappa;
        RealType m_theta_perp;
        RealType m_theta_par;
        point_type m_ub;
        RealType m_maxNormalizedVelocity;
    };

    /// Sentinel selecting the untruncated sampler, and **the default**.  In this
    /// mode the draw uses a fixed sequence of high-level variates with no outer
    /// acceptance-rejection loop, and the sampled law is the full bi-Kappa
    /// distribution.
    ///
    /// Naming it explicitly is equivalent to omitting the argument, and is worth
    /// doing where the target law matters to a reader:
    ///
    ///     bi_kappa_distribution<double> d(2.0, 1.0, 2.0, {0,0,1},
    ///                                     bi_kappa_distribution<double>::no_cap());
    ///
    /// With a finite cap the sampled law is instead that distribution conditioned
    /// on all normalized components lying inside the component-wise box -- a
    /// different, bounded law with its own normalization.  That is opt-in: it is
    /// never selected by default, because at kappa <= 3/2 no practical cap value
    /// is a good approximation to the bi-Kappa law and the conditioned law is not
    /// even axisymmetric about B.
    ///
    /// The cap is tested in the field-aligned frame, before the rotation into the
    /// global frame, so the returned components are *not* individually bounded by
    /// cap*theta unless ub is {0,0,1}; the box rotates with ub.  The speed is
    /// bounded in every frame, by cap*sqrt(2*theta_perp^2 + theta_par^2).
    ///
    /// \see \ref usage_bikappa_cap for what the capped law is, what it costs, and
    /// why a small rejected fraction does not imply faithful tails.
    static constexpr RealType no_cap()
    {
        return std::numeric_limits<RealType>::infinity();
    }

    // ====================================================================
    // Construction
    // ====================================================================
    bi_kappa_distribution()
        : m_params(), m_gen(std::random_device{}()), m_normalCache(RealType(0.0)),
          m_normalCached(false), m_attempts(0), m_nonfinite(0)
    {
    }

    explicit bi_kappa_distribution(RealType kappa, RealType theta_perp= RealType(1.0),
                                   RealType theta_par= RealType(1.0),
                                   point_type ub= point_type({RealType(0.0), RealType(0.0),
                                                              RealType(1.0)}),
                                   RealType max_normalized_velocity= no_cap())
        : m_params(kappa, theta_perp, theta_par, ub, max_normalized_velocity),
          m_gen(std::random_device{}()), m_normalCache(RealType(0.0)), m_normalCached(false),
          m_attempts(0), m_nonfinite(0)
    {
        define(kappa, theta_perp, theta_par, ub, max_normalized_velocity);
    }

    explicit bi_kappa_distribution(const param_type &params)
        : m_params(params), m_gen(std::random_device{}()), m_normalCache(RealType(0.0)),
          m_normalCached(false), m_attempts(0), m_nonfinite(0)
    {
        define(params);
    }

    // ====================================================================
    // Setters
    // ====================================================================

    // Bulk configuration
    void define(const param_type &params)
    {
        m_params= params;
        reset(); // just to clear cache
    }

    /**
     * @brief Configure the sampler.  Validates everything; throws rather than
     * deferring a bad parameter to sampling time.
     *
     * @param kappa   Spectral index.  Smaller means heavier suprathermal tails;
     *                kappa -> infinity recovers the bi-Maxwellian.  Must be > 0.5.
     *                Note that moments of order l exist only for kappa > (l+1)/2,
     *                so for kappa <= 1.5 the loaded population has no finite
     *                variance and no convergent temperature.
     * @param theta_perp  Thermal speed transverse to @p ub, > 0.  A speed, not a
     *                temperature: theta = sqrt(2*k_B*T/m).
     * @param theta_par   Thermal speed along @p ub, > 0.
     * @param ub      Magnetic-field direction {bx,by,bz} in the global frame.
     *                Need not be normalized; must be non-zero.  Samples are drawn
     *                in the field-aligned frame and rotated into the global frame
     *                before return.
     * @param max_normalized_velocity  Velocity cap in thermal speeds, or no_cap()
     *                (the default) for the untruncated bi-Kappa law.  A finite
     *                value selects a *different*, bounded target law -- see
     *                \ref usage_bikappa_cap.  Must be > 0; passing 0 to mean
     *                "no cap" throws.
     *
     * @note To seed reproducibly, either call seed() afterwards or use the
     * overload taking a trailing int: >= 0 is reproducible, < 0 reseeds from
     * std::random_device.  Without either, the stream stays as seeded at
     * construction, which is *not* reproducible.
     *
     * @throws std::invalid_argument if kappa <= 0.5, either theta is non-positive or
     *                non-finite, `ub` is zero or has a non-finite component, any of these is
     *                NaN,
     *         @p max_normalized_velocity is non-positive, or @p ub is the zero vector.
     *
     * @see \ref parameters for the full table, \ref choosing to compare samplers.
     */
    void define(RealType kappa, RealType theta_perp= RealType(1.0),
                RealType theta_par= RealType(1.0),
                point_type ub= point_type({RealType(0.0), RealType(0.0), RealType(1.0)}),
                RealType max_normalized_velocity= no_cap())
    {
        define(param_type(kappa, theta_perp, theta_par, ub, max_normalized_velocity));
    }

    void define(const param_type &params, int seed) { applySeed(m_gen, seed); define(params); }

    void define(RealType kappa, RealType theta_perp, RealType theta_par,
                point_type ub, RealType max_normalized_velocity, int seed)
    {
        applySeed(m_gen, seed);
        define(kappa, theta_perp, theta_par, ub, max_normalized_velocity);
    }

    // Individual parameter setters
    void kappa(RealType value)
    {
        define(value, m_params.theta_perp(), m_params.theta_par(), m_params.ub(),
               m_params.max_normalized_velocity());
    }

    void theta_perp(RealType value)
    {
        define(m_params.kappa(), value, m_params.theta_par(), m_params.ub(),
               m_params.max_normalized_velocity());
    }

    void theta_par(RealType value)
    {
        define(m_params.kappa(), m_params.theta_perp(), value, m_params.ub(),
               m_params.max_normalized_velocity());
    }

    void ub(const point_type &value)
    {
        define(m_params.kappa(), m_params.theta_perp(), m_params.theta_par(), value,
               m_params.max_normalized_velocity());
    }

    /// Unlike the other setters, this one changes *which distribution* is sampled,
    /// not merely its parameters.  Pass no_cap() to return to the bi-Kappa law
    /// itself; passing 0 to mean "no cap" throws, as every non-positive value does.
    ///
    ///     dist.max_normalized_velocity(50.0);                  // capped law
    ///     dist.max_normalized_velocity(decltype(dist)::no_cap());  // bi-Kappa law
    void max_normalized_velocity(RealType value)
    {
        define(m_params.kappa(), m_params.theta_perp(), m_params.theta_par(), m_params.ub(),
               value);
    }

    // ====================================================================
    // Getters
    // ====================================================================

    // Bulk getter
    param_type param() const { return m_params; }

    // Individual parameter getter
    RealType kappa() const { return m_params.kappa(); }
    RealType theta_perp() const { return m_params.theta_perp(); }
    RealType theta_par() const { return m_params.theta_par(); }
    point_type ub() const { return m_params.ub(); }
    RealType max_normalized_velocity() const { return m_params.max_normalized_velocity(); }

    // ====================================================================
    // Random engine control
    // ====================================================================
    /// Reseed the internal engine.  This also calls reset(), so a reseeded stream restarts
    /// cleanly: without it the cached second normal variate of the previous stream would be
    /// consumed by the first draw of the new one, and `seed(s); draw()` would not reproduce.
    void seed(int s)
    {
        applySeed(m_gen, s);
        reset();
    }

    /// Discard any variate cached inside the sampler and zero the observability counters.
    /// Called by seed() and by every define() overload.
    void reset()
    {
        m_normalCached= false;
        m_normalCache= RealType(0.0);
        m_attempts= 0;
        m_nonfinite= 0;
    }

    // ====================================================================
    // Observability
    // ====================================================================

    /// Core-mapping attempts made since the last reset(), define() or seed().  In
    /// untruncated mode this equals the number of returned samples; in capped mode it also
    /// counts the draws the box rejected.
    unsigned long long n_attempts() const { return m_attempts; }

    /// Attempts since the last reset(), define() or seed() whose intended velocity vector
    /// had no representation in @c RealType.
    ///
    /// Sufficiently near kappa = 1/2 the bi-Kappa law puts non-zero probability outside
    /// every finite floating-point range, so such draws exist as mathematics and cannot be
    /// returned as numbers.  Untruncated mode returns them as non-finite components -- the
    /// honest answer, and unchanged from earlier releases -- and counts them here.  Capped
    /// mode rejects them, as it rejects anything outside the box, and still counts them
    /// here, so that a caller can tell "the box rejected a representable sample" from "the
    /// box rejected a sample that could not have been returned at all".  In both modes the
    /// vector whose representability is counted is the *velocity* the caller would have
    /// received, rotation and thermal scaling included -- not the normalized coordinate the
    /// cap predicate is written in, which understates the largest component whenever
    /// `theta_par` and `theta_perp` differ.  Without this
    /// counter the second case is invisible: the returned samples are all finite and the
    /// acceptance rate alone does not distinguish the two.
    unsigned long long n_nonfinite() const { return m_nonfinite; }

    // ====================================================================
    // Sampling
    // ====================================================================
    template <typename GeneratorType> result_type operator()(GeneratorType &gen)
    {
        return (*this)(gen, m_params);
    }

    result_type operator()()
    {
        return (*this)(m_gen);
    }

    template <typename GeneratorType>
    result_type operator()(GeneratorType &gen, const param_type &params)
    {
        typedef typename bikappa_detail::log_accumulator<RealType>::type AccumType;
        const RealType kHalf= RealType(0.5);
        const RealType kOne= RealType(1.0);
        const AccumType kHalfA= AccumType(0.5);
        const AccumType kOneA= AccumType(1.0);
        const int kMaxCapRejectTries= 1000000;

        RealType kappa= params.kappa();
        RealType theta_perp= params.theta_perp();
        RealType theta_par= params.theta_par();
        point_type ub= params.ub();
        RealType maxNormalizedVelocity= params.max_normalized_velocity();

        const RealType a= kappa - kHalf;
        // sqrt(kappa) is not a variate; it is part of the deterministic map from the drawn
        // variates to the returned velocity, so it is formed in the accumulator with
        // everything else on that path.
        const AccumType sqrtKappa= std::sqrt(static_cast<AccumType>(kappa));

        // Untruncated mode runs the core mapping exactly once; capped mode redraws
        // until the sample falls inside the box.
        const bool capped= params.capped();
        const int maxTries= capped ? kMaxCapRejectTries : 1;

        const bool axisAligned= approximatelyEqual(ub[0], RealType(0.0)) &&
                                approximatelyEqual(ub[1], RealType(0.0)) &&
                                approximatelyEqual(ub[2], RealType(1.0));

        // Four logarithms that only the capped branch reads, and so are only taken when
        // there is a cap.  Untruncated mode decides representability from the components
        // themselves and needs no logarithm at all outside step 1; computing these anyway
        // cost two `std::log` calls on every untruncated draw, which `-fmath-errno` (the
        // default) forbids a compiler from removing.
        //
        // `logMax` is a screening quantity only: the cap branch uses it to decide whether a
        // rejected draw is near enough to the type's limit to be worth materializing.
        // Representability itself is never decided from a logarithm -- see
        // bikappa_detail::materializeComponents.
        //
        // `logGBound` bounds the *velocity* from its radius alone:
        // |g_j| <= |g| <= sqrt(kappa) max(theta_perp, theta_par) for every direction and
        // every rotation, because a rotation preserves the norm.  That is what the capped
        // branch needs and what the normalized coordinate the cap test is written in cannot
        // give it.  It is finite because define() rejects a non-positive or non-finite
        // theta, so the screen below is never comparing against a NaN.
        // All four are formed in the accumulator, so that the screen and the cap predicate
        // resolve the same margins the radius does.
        const AccumType logMax=
            capped ? std::log(static_cast<AccumType>(std::numeric_limits<RealType>::max()))
                   : AccumType(0.0);
        const AccumType logCap=
            capped ? std::log(static_cast<AccumType>(maxNormalizedVelocity)) : AccumType(0.0);
        const AccumType logSqrtKappa= capped ? std::log(sqrtKappa) : AccumType(0.0);
        const AccumType logGBound=
            capped ? logSqrtKappa +
                         std::log(static_cast<AccumType>(std::max(theta_perp, theta_par)))
                   : AccumType(0.0);

        int tries= 0;
        for (tries= 0; tries < maxTries; ++tries)
        {
            ++m_attempts;

            // Step 1: the radius, carried as its logarithm.
            //
            // R = sqrt(X1)/sqrt(X2) with X1 ~ Gamma(3/2,1) and X2 ~ Gamma(a,1), a = kappa-1/2.
            // As kappa -> 1/2 the shape a -> 0 and X2 underflows to exactly zero, which no
            // rearrangement of sqrt(X1)/sqrt(X2) can survive.  X2 is therefore never formed:
            // log X2 comes from the shape-boosting identity and
            //
            //     log R = (log X1 - log X2) / 2
            //
            // stays finite and accurate wherever the mathematical radius exists.  X1 has
            // shape 3/2 > 1, so it is drawn directly and its logarithm is well conditioned.
            const RealType x1=
                bikappa_detail::marsagliaTsangGamma<RealType, GeneratorType>(
                    RealType(1.5), gen, m_normalCache, m_normalCached);
            const AccumType logX2=
                bikappa_detail::logGamma<RealType, GeneratorType>(a, gen, m_normalCache,
                                                                  m_normalCached);
            const AccumType logR=
                kHalfA * (std::log(static_cast<AccumType>(x1)) - logX2);

            // Step 2: a direction uniform on S^2, by the rejection method of Marsaglia
            // (1972), Ann. Math. Statist. 43(2), 645-646.  Draw a point uniformly in the
            // unit disc and lift it; the lift is area-preserving, so the result is exactly
            // uniform on the sphere.
            //
            // The textbook alternative -- cos(theta) uniform on (-1,1) and phi uniform on
            // (0, 2pi) -- is equally correct as mathematics, and is what this loader used
            // before 2.0.0.  It needs sin(phi) and cos(phi) of the same argument, and a
            // compiler is free to fuse that pair into one combined routine.  Clang does, on
            // Darwin, into __sincos_stret, whose sine differs from the standalone sin by one
            // unit in the last place for roughly one argument in a thousand; GCC calls the
            // two separately.  Two builds of this header then return different numbers for
            // the same seed, and the difference belongs to the compiler rather than to the
            // sampler.  The rejection method uses only sqrt, which IEEE-754 requires to be
            // correctly rounded, so the direction is identical under every compiler.  It is
            // also faster here, despite averaging 4/pi = 1.273 attempts per direction.
            //
            // The two coordinates are the variates and are drawn in RealType; the disc test
            // and the lift are the deterministic part and are evaluated in the accumulator.
            // For a `float` run that also makes the test exact -- `d1*d1 + d2*d2` of two
            // `float`s needs 49 bits, which `double` holds and `float` does not, so the
            // `float` form could accept a point outside the disc or reject one inside it.
            RealType d1, d2;
            AccumType dsq;
            do
            {
                d1= RealType(2.0) *
                    bikappa_detail::openCanonical<RealType, GeneratorType>(gen) - kOne;
                d2= RealType(2.0) *
                    bikappa_detail::openCanonical<RealType, GeneratorType>(gen) - kOne;
                dsq= static_cast<AccumType>(d1) * static_cast<AccumType>(d1) +
                     static_cast<AccumType>(d2) * static_cast<AccumType>(d2);
            } while (dsq >= kOneA);
            // dsq == 0 is unreachable: openCanonical never returns exactly 1/2, so neither
            // d1 nor d2 is ever exactly zero.
            const AccumType dscale= AccumType(2.0) * std::sqrt(kOneA - dsq);

            const AccumType nx= static_cast<AccumType>(d1) * dscale;
            const AccumType ny= static_cast<AccumType>(d2) * dscale;
            const AccumType nz= kOneA - AccumType(2.0) * dsq;

            // Step 3: the cap, tested in the field-aligned frame before anything is
            // exponentiated.  This is the same predicate the linear form applies --
            // |point_j| / theta_j = sqrt(kappa) R |n_j| <= lambda -- moved into the log
            // domain so that an intended sample far outside the box is rejected rather
            // than overflowed first and rejected afterwards.  Only the largest component
            // has to be examined, so this costs one logarithm rather than three.
            if (capped)
            {
                const AccumType maxAbsN=
                    std::max(std::abs(nx), std::max(std::abs(ny), std::abs(nz)));
                const AccumType logLargest= logR + logSqrtKappa + std::log(maxAbsN);
                if (logLargest > logCap)
                {
                    // An intended sample that is outside the box and also unrepresentable is
                    // still counted, so a caller can tell the two apart.  Two things have to
                    // be right here, and the released 2.0.0 form had neither.
                    //
                    // *What* is tested.  Representability is a property of the velocity
                    // `R g_j`, not of the normalized coordinate `sqrt(kappa) R |n_j|` the cap
                    // test just used: the two differ by `theta` and by the rotation.  The
                    // ratio is bounded by `sqrt(3) max(theta_perp, theta_par)/min(...)`,
                    // since `|w|/sqrt(3) <= max_j |w_j| <= |w|` -- at
                    // `theta_par/theta_perp = 2` that is 3.46, and 2.33 is the largest value
                    // the C5/C6 geometry actually reaches.  The rotation concentrates the
                    // largest component as readily as it spreads it, which is why the bound
                    // is not simply the theta ratio.  The counter must answer the question
                    // n_nonfinite() documents, which is about the velocity.
                    //
                    // *How* it is tested.  By materializing, for the same reason as step 5.
                    // The exponential is worth avoiding here -- in capped mode this is the
                    // hot branch -- so `logGBound` screens first.  It is an upper bound on
                    // log|V_j| that holds for every direction and every rotation, so a draw a
                    // whole natural-log unit below the limit by that bound cannot overflow.
                    //
                    // The screen is sound, not merely generous.  Its error is dominated by
                    // half an ulp of `logR + logGBound` at the limit, which is 88.72 * 2^-24
                    // = 5.3e-6 in `float` and 709.78 * 2^-53 = 7.9e-14 in `double`; with the
                    // libm and normalization terms the worst case is about 1.1e-5 and
                    // 1.0e-13.  One natural-log unit exceeds those by 9.4e4 and 1.0e13.
                    //
                    // How often it fires is a property of the configuration, not a constant:
                    // about 7e-7 of attempts at double kappa = 0.51, but 17 % at float
                    // kappa = 0.51 and 87 % at double kappa = 0.5001, where most of the law
                    // genuinely lies outside the type.  Where it fires often, it is doing
                    // the work it exists for.
                    if (logR + logGBound > logMax - kOneA)
                    {
                        accum_point gRej= fieldAlignedG(sqrtKappa, theta_perp, theta_par,
                                                        nx, ny, nz);
                        if (!axisAligned)
                        {
                            gRej= rotateFromFieldAlignedFrameIn<AccumType>(
                                gRej, widen<AccumType>(ub));
                        }
                        point_type scratch;
                        if (bikappa_detail::materializeComponents<RealType>(logR, gRej,
                                                                            scratch))
                        {
                            ++m_nonfinite;
                        }
                    }
                    continue;
                }
            }

            // Step 4: the order-unity vector g, with V = R g.  Forming g first is what makes
            // the representability test meaningful: an orthogonal rotation can spread one
            // over-large local component across three representable global ones, so the test
            // has to be applied to the components that are actually returned.
            accum_point g= fieldAlignedG(sqrtKappa, theta_perp, theta_par, nx, ny, nz);
            if (!axisAligned)
            {
                g= rotateFromFieldAlignedFrameIn<AccumType>(g, widen<AccumType>(ub));
            }

            // Step 5: materialize the components and decide representability from the
            // components themselves, never from their logarithms.  The whole of that
            // decision lives in bikappa_detail::materializeComponents, which documents why
            // a logarithmic test is wrong in both directions; the capped branch above calls
            // the same routine, so the counter and the return can never disagree.
            point_type point;
            const bool anyOverflow=
                bikappa_detail::materializeComponents<RealType>(logR, g, point);

            if (anyOverflow)
            {
                // Honest overflow: the intended velocity has no representation in RealType.
                // No formation can rescue it, so it is returned as a non-finite vector and
                // counted, never silently redrawn.  In capped mode this is unreachable for
                // every cap a caller would write: an accepted draw has
                // sqrt(kappa) R max_j|n_j| <= lambda, hence
                // max_j |g_j| R <= sqrt(3) lambda max(theta_perp, theta_par), so reaching it
                // would take lambda max(theta) above 1.97e38 in `float` or 1.04e308 in
                // `double` -- a box larger than the type.
                ++m_nonfinite;
            }

            return point;
        }

        throw std::runtime_error("failed to sample velocity under max_normalized_velocity cap: "
                                 "rejection sampling exceeded internal attempt limit");
    }

    // ====================================================================
    // Geometry helper
    // ====================================================================
    // Transform from velocity components (v_perp, v_par) to vector velocity (vx, vy, vz)
    // given the parallel direction ub.
    point_type rotate_from_fieldAligned_frame(const point_type &local, const point_type &ub) const
    {
        return rotateFromFieldAlignedFrameIn<RealType>(local, widen<RealType>(ub));
    }

    /// The same rotation, carried out in @p ValueType.
    ///
    /// The basis is a function of `ub` alone, so evaluating it in the accumulator changes
    /// no decision and no stream -- it only stops the rotation from spending an ulp of the
    /// returned type before the component is formed.  With `ValueType = RealType` this is
    /// the released 2.1.0 routine, character for character, which is why a `double`
    /// instantiation is unaffected.
    template <typename ValueType>
    static std::array<ValueType, 3> rotateFromFieldAlignedFrameIn(
        const std::array<ValueType, 3> &local, const std::array<ValueType, 3> &ub)
    {
        typedef std::array<ValueType, 3> vec_type;
        const ValueType kZero= ValueType(0.0);
        const ValueType kOne= ValueType(1.0);

        // Normalize ub to get e3
        ValueType ubNorm2= ub[0] * ub[0] + ub[1] * ub[1] + ub[2] * ub[2];
        if (approximatelyEqualIn<ValueType>(ubNorm2, kZero))
        {
            throw std::invalid_argument("ub must be a non-zero vector");
        }

        ValueType ubNorm= std::sqrt(ubNorm2);

        vec_type e3;
        e3[0]= ub[0] / ubNorm;
        e3[1]= ub[1] / ubNorm;
        e3[2]= ub[2] / ubNorm;

        // Build e2 perpendicular to e3
        vec_type e2;
        e2[0]= kOne;
        e2[1]= kOne;
        e2[2]= kOne;

        int maxcomp= 0;
        if (std::abs(e3[1]) > std::abs(e3[maxcomp]))
        {
            maxcomp= 1;
        }
        if (std::abs(e3[2]) > std::abs(e3[maxcomp]))
        {
            maxcomp= 2;
        }

        ValueType sum= e3[0] + e3[1] + e3[2];
        e2[maxcomp]= kOne - sum / e3[maxcomp];

        ValueType e2Norm2= e2[0] * e2[0] + e2[1] * e2[1] + e2[2] * e2[2];
        ValueType e2Norm= std::sqrt(e2Norm2);
        e2[0]/= e2Norm;
        e2[1]/= e2Norm;
        e2[2]/= e2Norm;

        // e1 = e2 x e3
        vec_type e1;
        e1[0]= e2[1] * e3[2] - e2[2] * e3[1];
        e1[1]= e2[2] * e3[0] - e2[0] * e3[2];
        e1[2]= e2[0] * e3[1] - e2[1] * e3[0];

        // Transform local coordinates to Cartesian coordinates
        vec_type point;
        point[0]= local[0] * e1[0] + local[1] * e2[0] + local[2] * e3[0];
        point[1]= local[0] * e1[1] + local[1] * e2[1] + local[2] * e3[1];
        point[2]= local[0] * e1[2] + local[1] * e2[2] + local[2] * e3[2];

        return point;
    }

    // ====================================================================
    // Comparison
    // ====================================================================
    bool operator==(const bi_kappa_distribution &other) const { return m_params == other.m_params; }

    bool operator!=(const bi_kappa_distribution &other) const { return !(*this == other); }

  private:
    // ====================================================================
    // Internal helpers
    // ====================================================================
    static void applySeed(std::mt19937 &gen, int s)
    {
        if (s < 0)
            gen.seed(std::random_device{}());
        else
            gen.seed(static_cast<unsigned int>(s));
    }

    /// The order-unity vector in the field-aligned frame, `V = R g` before any rotation,
    /// in the accumulator of @ref bikappa_detail::log_accumulator.
    ///
    /// Shared by step 4 and by the capped branch of step 3, so the vector whose
    /// representability is counted is the vector that would have been returned.  The
    /// direction and the parameters are the working-precision values the sampler drew and
    /// the caller supplied; only the three products are widened, which is what lets the
    /// component carry a single rounding.  With `RealType = double` the accumulator is
    /// `double` and this is the released 2.1.0 expression unchanged.
    static accum_point fieldAlignedG(
        typename bikappa_detail::log_accumulator<RealType>::type sqrtKappa,
        RealType theta_perp, RealType theta_par,
        typename bikappa_detail::log_accumulator<RealType>::type nx,
        typename bikappa_detail::log_accumulator<RealType>::type ny,
        typename bikappa_detail::log_accumulator<RealType>::type nz)
    {
        typedef typename bikappa_detail::log_accumulator<RealType>::type A;
        accum_point g;
        g[0]= sqrtKappa * static_cast<A>(theta_perp) * nx;
        g[1]= sqrtKappa * static_cast<A>(theta_perp) * ny;
        g[2]= sqrtKappa * static_cast<A>(theta_par) * nz;
        return g;
    }

    /// `ub`, widened componentwise.  Widening is exact, so the rotated basis is the basis
    /// of the `ub` the caller gave, evaluated more finely -- not a different frame.
    template <typename ValueType>
    static std::array<ValueType, 3> widen(const point_type &v)
    {
        std::array<ValueType, 3> w;
        w[0]= static_cast<ValueType>(v[0]);
        w[1]= static_cast<ValueType>(v[1]);
        w[2]= static_cast<ValueType>(v[2]);
        return w;
    }

    template <typename ValueType> static bool approximatelyEqualIn(ValueType a, ValueType b)
    {
        ValueType scale= std::max(std::max(ValueType(1.0), std::abs(a)), std::abs(b));
        return std::abs(a - b) <= std::numeric_limits<ValueType>::epsilon() * scale;
    }

    static bool approximatelyEqual(RealType a, RealType b)
    {
        return approximatelyEqualIn<RealType>(a, b);
    }

    /// The cap predicate in linear form: every normalized component inside the box.
    ///
    /// The sampler itself no longer calls this -- it applies the algebraically identical
    /// test `sqrt(kappa) R |n_j| <= lambda` in the log domain, before exponentiating, so that
    /// a sample far outside the box is rejected rather than overflowed first and rejected
    /// afterwards.  The linear form is kept because it is the published definition of the
    /// capped law, it is what `experiments/exp2_cap_characterization` transcribes, and a
    /// regression test checks that the two predicates accept exactly the same set wherever
    /// both are representable.
    static bool withinNormalizedVelocityCap(const point_type &point, RealType theta_perp,
                                            RealType theta_par, RealType max_normalized_velocity)
    {
        RealType norm_x= std::abs(point[0]) / theta_perp;
        RealType norm_y= std::abs(point[1]) / theta_perp;
        RealType norm_z= std::abs(point[2]) / theta_par;
        return norm_x <= max_normalized_velocity && norm_y <= max_normalized_velocity &&
               norm_z <= max_normalized_velocity;
    }

    param_type m_params;
    std::mt19937 m_gen;
    RealType m_normalCache;      ///< second variate of the current polar-normal pair
    bool m_normalCached;         ///< whether m_normalCache holds an unconsumed variate
    unsigned long long m_attempts;
    unsigned long long m_nonfinite;
};

#endif // _BI_KAPPA_DISTRIBUTION_H_
