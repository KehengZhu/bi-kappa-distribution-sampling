#ifndef _BI_KAPPA_DISTRIBUTION_H_
#define _BI_KAPPA_DISTRIBUTION_H_

/// Release version of this header.  A consumer that vendors the file can test these to tell
/// which sampler it has; 2.0.0 is the first release whose radius is built in the log domain,
/// and its random stream differs from 1.x for every seed.
#define BI_KAPPA_VERSION_MAJOR 2
#define BI_KAPPA_VERSION_MINOR 0
#define BI_KAPPA_VERSION_PATCH 0
#define BI_KAPPA_VERSION_STRING "2.0.0"

#include <array>
#include <algorithm>
#include <cmath>
#include <limits>
#include <random>
#include <stdexcept>

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
RealType logGamma(RealType a, GeneratorType &gen, RealType &ncache, bool &ncached)
{
    const RealType y=
        marsagliaTsangGamma<RealType, GeneratorType>(a + RealType(1.0), gen, ncache, ncached);
    const RealType u= openCanonical<RealType, GeneratorType>(gen);
    return std::log(y) + std::log(u) / a;
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
 * `log X2` and `log R`, and decides whether each returned component is representable from
 * `log|V_j| = log R + log|g_j|` before exponentiating anything.  Two consequences a caller
 * should know about:
 *
 *  - **The random stream differs from releases before 2.0.0.**  The Gamma, normal and
 *    uniform primitives are implemented in @ref bikappa_detail rather than taken from
 *    `<random>`, so a given seed produces a different sequence than it used to.  In exchange
 *    the sequence no longer depends on which standard library is in use: the same seed gives
 *    bit-identical output under libc++ and libstdc++.
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
            if (m_kappa <= RealType(0.5))
            {
                throw std::invalid_argument("kappa must be > 0.5 for normalization");
            }
            if (m_theta_perp <= RealType(0.0) || m_theta_par <= RealType(0.0))
            {
                throw std::invalid_argument("Temperature parameters must be positive");
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
     * @throws std::invalid_argument if kappa <= 0.5, either theta is non-positive,
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
    /// box rejected a sample that could not have been returned at all".  Without this
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
        const RealType kHalf= RealType(0.5);
        const RealType kOne= RealType(1.0);
        const int kMaxCapRejectTries= 1000000;

        RealType kappa= params.kappa();
        RealType theta_perp= params.theta_perp();
        RealType theta_par= params.theta_par();
        point_type ub= params.ub();
        RealType maxNormalizedVelocity= params.max_normalized_velocity();

        const RealType a= kappa - kHalf;
        RealType sqrtKappa= std::sqrt(kappa);

        // Untruncated mode runs the core mapping exactly once; capped mode redraws
        // until the sample falls inside the box.
        const bool capped= params.capped();
        const int maxTries= capped ? kMaxCapRejectTries : 1;

        const bool axisAligned= approximatelyEqual(ub[0], RealType(0.0)) &&
                                approximatelyEqual(ub[1], RealType(0.0)) &&
                                approximatelyEqual(ub[2], RealType(1.0));

        // The largest finite magnitude, and its logarithm.  A component is representable
        // exactly when log|V_j| is below this; the comparison is made in the log domain,
        // before anything is exponentiated.
        const RealType logMax= std::log(std::numeric_limits<RealType>::max());
        const RealType logCap= capped ? std::log(maxNormalizedVelocity)
                                      : std::numeric_limits<RealType>::infinity();
        const RealType logSqrtKappa= std::log(sqrtKappa);
        // |g_j| never exceeds sqrt(kappa) max(theta_perp, theta_par), because a rotation
        // preserves the norm; see step 5.
        const RealType logGBound= logSqrtKappa + std::log(std::max(theta_perp, theta_par));

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
            const RealType logX2=
                bikappa_detail::logGamma<RealType, GeneratorType>(a, gen, m_normalCache,
                                                                  m_normalCached);
            const RealType logR= kHalf * (std::log(x1) - logX2);

            // Step 2: a direction uniform on S^2, from two open-interval uniforms.
            const RealType cosTheta=
                RealType(2.0) * bikappa_detail::openCanonical<RealType, GeneratorType>(gen) - kOne;
            const RealType sinTheta=
                std::sqrt(std::max(RealType(0.0), kOne - cosTheta * cosTheta));
            const RealType phi= RealType(2.0) * piValue() *
                                bikappa_detail::openCanonical<RealType, GeneratorType>(gen);

            const RealType nx= sinTheta * std::cos(phi);
            const RealType ny= sinTheta * std::sin(phi);
            const RealType nz= cosTheta;

            // Step 3: the cap, tested in the field-aligned frame before anything is
            // exponentiated.  This is the same predicate the linear form applies --
            // |point_j| / theta_j = sqrt(kappa) R |n_j| <= lambda -- moved into the log
            // domain so that an intended sample far outside the box is rejected rather
            // than overflowed first and rejected afterwards.  Only the largest component
            // has to be examined, so this costs one logarithm rather than three.
            if (capped)
            {
                const RealType maxAbsN=
                    std::max(std::abs(nx), std::max(std::abs(ny), std::abs(nz)));
                const RealType logLargest= logR + logSqrtKappa + std::log(maxAbsN);
                if (logLargest > logCap)
                {
                    // An intended sample that is outside the box and also unrepresentable
                    // is still counted, so a caller can tell the two apart.
                    if (logLargest >= logMax)
                    {
                        ++m_nonfinite;
                    }
                    continue;
                }
            }

            // Step 4: the order-unity vector g, with V = R g.  Forming g first is what makes
            // the representability test meaningful: an orthogonal rotation can spread one
            // over-large local component across three representable global ones, so the test
            // has to be applied to the components that are actually returned.
            point_type g;
            g[0]= sqrtKappa * theta_perp * nx;
            g[1]= sqrtKappa * theta_perp * ny;
            g[2]= sqrtKappa * theta_par * nz;
            if (!axisAligned)
            {
                g= rotate_from_fieldAligned_frame(g, ub);
            }

            // Step 5: decide representability from log|V_j| = log R + log|g_j|, and only
            // then exponentiate.  Where the radius itself is representable the product
            // R * g_j keeps a full working-precision result; the per-component exponential
            // is used only where it is not.
            //
            // A rotation is orthogonal, so |g_j| <= |g| <= sqrt(kappa) max(theta_perp,
            // theta_par) =: gBound whatever ub is.  When log R + log gBound is below the
            // type's limit every component is representable, and the three per-component
            // logarithms can be skipped entirely.  That is the overwhelmingly common case;
            // the exact per-component test below runs only near the limit.
            const bool radiusRepresentable= (logR <= logMax);
            const RealType R= radiusRepresentable ? std::exp(logR) : RealType(0.0);

            point_type point;
            bool anyOverflow= false;
            if (radiusRepresentable && logR + logGBound < logMax)
            {
                point[0]= R * g[0];
                point[1]= R * g[1];
                point[2]= R * g[2];
            }
            else
            {
                for (int j= 0; j < 3; ++j)
                {
                    if (g[j] == RealType(0.0))
                    {
                        point[j]= RealType(0.0);
                        continue;
                    }
                    const RealType logComponent= logR + std::log(std::abs(g[j]));
                    if (logComponent >= logMax)
                    {
                        anyOverflow= true;
                        point[j]= (g[j] < RealType(0.0))
                                      ? -std::numeric_limits<RealType>::infinity()
                                      : std::numeric_limits<RealType>::infinity();
                        continue;
                    }
                    if (radiusRepresentable)
                    {
                        point[j]= R * g[j];
                    }
                    else
                    {
                        const RealType mag= std::exp(logComponent);
                        point[j]= (g[j] < RealType(0.0)) ? -mag : mag;
                    }
                }
            }

            if (anyOverflow)
            {
                // Honest overflow: the intended velocity has no representation in RealType.
                // No formation can rescue it, so it is returned as a non-finite vector and
                // counted, never silently redrawn.  In capped mode it is unreachable, having
                // already been rejected by the box test above.
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
        const RealType kZero= RealType(0.0);
        const RealType kOne= RealType(1.0);

        // Normalize ub to get e3
        RealType ubNorm2= ub[0] * ub[0] + ub[1] * ub[1] + ub[2] * ub[2];
        if (approximatelyEqual(ubNorm2, kZero))
        {
            throw std::invalid_argument("ub must be a non-zero vector");
        }

        RealType ubNorm= std::sqrt(ubNorm2);

        point_type e3;
        e3[0]= ub[0] / ubNorm;
        e3[1]= ub[1] / ubNorm;
        e3[2]= ub[2] / ubNorm;

        // Build e2 perpendicular to e3
        point_type e2({kOne, kOne, kOne});

        int maxcomp= 0;
        if (std::abs(e3[1]) > std::abs(e3[maxcomp]))
        {
            maxcomp= 1;
        }
        if (std::abs(e3[2]) > std::abs(e3[maxcomp]))
        {
            maxcomp= 2;
        }

        RealType sum= e3[0] + e3[1] + e3[2];
        e2[maxcomp]= kOne - sum / e3[maxcomp];

        RealType e2Norm2= e2[0] * e2[0] + e2[1] * e2[1] + e2[2] * e2[2];
        RealType e2Norm= std::sqrt(e2Norm2);
        e2[0]/= e2Norm;
        e2[1]/= e2Norm;
        e2[2]/= e2Norm;

        // e1 = e2 x e3
        point_type e1;
        e1[0]= e2[1] * e3[2] - e2[2] * e3[1];
        e1[1]= e2[2] * e3[0] - e2[0] * e3[2];
        e1[2]= e2[0] * e3[1] - e2[1] * e3[0];

        // Transform local coordinates to Cartesian coordinates
        point_type point;
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

    static RealType piValue() { return RealType(3.14159265358979323846264338327950288L); }

    static bool approximatelyEqual(RealType a, RealType b)
    {
        RealType scale= std::max({RealType(1.0), std::abs(a), std::abs(b)});
        return std::abs(a - b) <= std::numeric_limits<RealType>::epsilon() * scale;
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
