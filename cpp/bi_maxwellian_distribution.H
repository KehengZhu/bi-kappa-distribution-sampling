#ifndef _BI_MAXWELLIAN_DISTRIBUTION_H_
#define _BI_MAXWELLIAN_DISTRIBUTION_H_

#include <array>
#include <algorithm>
#include <cmath>
#include <limits>
#include <random>
#include <stdexcept>

/**
 * @brief Bi-Maxwellian distribution for 3D velocity sampling.
 *
 * Components are independent Gaussians with anisotropic thermal widths:
 * Vx ~ N(0, theta_perp^2/2), Vy ~ N(0, theta_perp^2/2), Vz ~ N(0, theta_par^2/2).
 *
 * @code{.cpp}
 * bi_maxwellian_distribution<double> dist;
 * dist.define(1.0, 2.0, {0, 0, 1});   // theta_perp, theta_par, ub
 * auto v = dist();                    // {vx, vy, vz}, global frame
 * @endcode
 *
 * @note **The cap defaults to 20 thermal speeds**, as in bi_kappa_distribution.  For a
 * bi-Maxwellian the probability beyond it is negligible.  Pass no_cap() to sample the bi-Maxwellian distribution itself; see
 * no_cap() and param_type::capped().
 *
 * @see \ref usage for more recipes (seeding, external RNG, reconfiguration).
 */
template <typename RealType= double> class bi_maxwellian_distribution
{
  public:
    using result_type= std::array<RealType, 3>;   ///< A velocity {vx, vy, vz}.
    using point_type= result_type;                  ///< Same as result_type.

    /// @brief The parameters of a bi_maxwellian_distribution: theta_perp, theta_par, ub and
    /// the velocity cap.
    class param_type
    {
      public:
        using distribution_type= bi_maxwellian_distribution<RealType>;   ///< The owning sampler.

        /// Store the parameters.  Throws std::invalid_argument if either theta is zero or
        /// negative, or if the cap is not positive; @p ub is not checked here.
        param_type(RealType theta_perp= RealType(1.0), RealType theta_par= RealType(1.0),
                   point_type ub= point_type({RealType(0.0), RealType(0.0), RealType(1.0)}),
                   RealType max_normalized_velocity= RealType(20.0))
            : m_theta_perp(theta_perp), m_theta_par(theta_par), m_ub(ub),
              m_maxNormalizedVelocity(max_normalized_velocity)
        {
            if (m_theta_perp <= RealType(0.0) || m_theta_par <= RealType(0.0))
            {
                throw std::invalid_argument("Bi-Maxwellian temperatures must be positive");
            }
            // A positive finite value selects the capped mode.  no_cap() (= +inf)
            // selects the untruncated mode, in which no acceptance-rejection loop
            // runs at all and samples follow the bi-Maxwellian distribution itself.
            if (!(m_maxNormalizedVelocity > RealType(0.0)))
            {
                throw std::invalid_argument(
                    "max_normalized_velocity must be positive, or no_cap() for the "
                    "untruncated distribution");
            }
        }

        RealType theta_perp() const { return m_theta_perp; }   ///< Perpendicular thermal speed.
        RealType theta_par() const { return m_theta_par; }     ///< Parallel thermal speed.
        point_type ub() const { return m_ub; }                 ///< Field direction, as given.
        /// Velocity cap in thermal speeds; +infinity for no_cap().
        RealType max_normalized_velocity() const { return m_maxNormalizedVelocity; }

        /// True when a finite component-wise cap is in force, that is, when the sampled
        /// distribution is the bi-Maxwellian distribution conditioned on the cap box rather
        /// than the bi-Maxwellian distribution itself.
        bool capped() const { return std::isfinite(m_maxNormalizedVelocity); }

        /// True if every parameter compares equal.
        bool operator==(const param_type &other) const
        {
            return m_theta_perp == other.m_theta_perp && m_theta_par == other.m_theta_par &&
                   m_ub == other.m_ub && m_maxNormalizedVelocity == other.m_maxNormalizedVelocity;
        }

        /// Negation of operator==.
        bool operator!=(const param_type &other) const { return !(*this == other); }

      private:
        RealType m_theta_perp;
        RealType m_theta_par;
        point_type m_ub;
        RealType m_maxNormalizedVelocity;
    };

    /// Value of the cap that selects the uncapped sampler: `+infinity`.  In this mode each
    /// draw is three Gaussians with no outer acceptance-rejection loop, and the samples
    /// follow the bi-Maxwellian distribution itself.
    ///
    /// It must be passed explicitly, because the default cap is 20:
    ///
    ///     bi_maxwellian_distribution<double> d(1.0, 2.0, {0,0,1},
    ///                                          bi_maxwellian_distribution<double>::no_cap());
    ///
    /// With a finite cap the sampled distribution is instead that distribution conditioned
    /// on all normalized components lying inside the component-wise box -- a
    /// different, bounded distribution with its own normalization.  The cap is tested in the
    /// field-aligned frame, before the rotation into the global frame, so the
    /// returned components are *not* individually bounded by cap*theta unless ub is
    /// {0,0,1}; the box rotates with ub.
    ///
    /// \see \ref usage_bikappa_cap, which describes the same construction for
    ///      bi_kappa_distribution.
    static constexpr RealType no_cap()
    {
        return std::numeric_limits<RealType>::infinity();
    }

    // ====================================================================
    // Construction
    // ====================================================================
    /// Construct with theta_perp = theta_par = 1, ub = {0, 0, 1} and a cap of 20.  The
    /// internal engine is seeded from `std::random_device`.
    bi_maxwellian_distribution() : m_params(), m_gen(std::random_device{}()), m_standardNormal(RealType(0.0), RealType(1.0)) {}

    /// Construct and call define() with these arguments.  The internal engine is seeded from
    /// `std::random_device`; call seed() for a reproducible stream.
    explicit bi_maxwellian_distribution(RealType theta_perp, RealType theta_par= RealType(1.0),
                                        point_type ub= point_type({RealType(0.0), RealType(0.0),
                                                                   RealType(1.0)}),
                                        RealType max_normalized_velocity= RealType(20.0))
        : m_params(theta_perp, theta_par, ub, max_normalized_velocity),
          m_gen(std::random_device{}()), m_standardNormal(RealType(0.0), RealType(1.0))
    {
        define(theta_perp, theta_par, ub, max_normalized_velocity);
    }

    /// Construct from a parameter set.  The internal engine is seeded from
    /// `std::random_device`.
    explicit bi_maxwellian_distribution(const param_type &params)
        : m_params(params), m_gen(std::random_device{}()), m_standardNormal(RealType(0.0), RealType(1.0))
    {
        define(params);
    }

    // ====================================================================
    // Setters
    // ====================================================================

    /// Replace all parameters with @p params, then call reset().
    void define(const param_type &params)
    {
        m_params= params;
        reset(); // just to clear cache
    }

    /**
     * @brief Configure the sampler.
     *
     * @param theta_perp  Thermal speed of the two components transverse to @p ub,
     *                > 0.  Each is N(0, theta_perp^2/2) -- a speed, not a
     *                temperature: theta = sqrt(2*k_B*T/m).
     * @param theta_par   Thermal speed along @p ub, > 0.
     * @param ub      Magnetic-field direction {bx,by,bz} in the global frame.
     *                Need not be normalized; must be non-zero, with a length above
     *                about 1.5e-8 in `double` (3.5e-4 in `float`); see \ref params_ub.
     * @param max_normalized_velocity  Velocity cap in thermal speeds (default 20), or
     *                no_cap() for the bi-Maxwellian distribution itself.  A finite value
     *                selects a *different*, bounded target distribution; see
     *                \ref usage_bikappa_cap.  Must be > 0; passing 0 to mean "no cap"
     *                throws.
     *
     * @note To seed reproducibly, either call seed() afterwards or use the
     * overload taking a trailing int: >= 0 is reproducible, < 0 reseeds from
     * std::random_device.
     *
     * @throws std::invalid_argument if either theta is zero or negative, or
     *         @p max_normalized_velocity is not positive (NaN included).
     *
     * @note define() does not check @p ub, and it does not reject a NaN or infinite theta.
     * A @p ub that is zero, or shorter than the limit above, passes define() and makes
     * every draw throw std::invalid_argument("ub must be a non-zero vector").
     *
     * @see \ref parameters for the full table, \ref choosing to compare samplers.
     */
    void define(RealType theta_perp, RealType theta_par= RealType(1.0),
                point_type ub= point_type({RealType(0.0), RealType(0.0), RealType(1.0)}),
                RealType max_normalized_velocity= RealType(20.0))
    {
        define(param_type(theta_perp, theta_par, ub, max_normalized_velocity));
    }

    /// define(params), after seeding the internal engine with @p seed as seed() does.
    void define(const param_type &params, int seed) { applySeed(m_gen, seed); define(params); }

    /// define() with all four parameters, after seeding the internal engine with @p seed as
    /// seed() does (`seed >= 0` reproducible, `seed < 0` from `std::random_device`).
    void define(RealType theta_perp, RealType theta_par,
                point_type ub, RealType max_normalized_velocity, int seed)
    {
        applySeed(m_gen, seed);
        define(theta_perp, theta_par, ub, max_normalized_velocity);
    }

    // Individual parameter setters.  Each calls define() with the other parameters
    // unchanged.

    /// Set theta_perp, keeping the other parameters; throws as define() does.
    void theta_perp(RealType value)
    {
        define(value, m_params.theta_par(), m_params.ub(), m_params.max_normalized_velocity());
    }

    /// Set theta_par, keeping the other parameters; throws as define() does.
    void theta_par(RealType value)
    {
        define(m_params.theta_perp(), value, m_params.ub(), m_params.max_normalized_velocity());
    }

    /// Set the field direction, keeping the other parameters.
    void ub(const point_type &value)
    {
        define(m_params.theta_perp(), m_params.theta_par(), value,
               m_params.max_normalized_velocity());
    }

    /// Set the velocity cap, keeping the other parameters.  Unlike the other setters, this
    /// one changes *which distribution* is sampled, not merely its parameters.  Pass
    /// no_cap() to return to the bi-Maxwellian distribution itself; passing 0 to mean
    /// "no cap" throws, as every non-positive value does.
    void max_normalized_velocity(RealType value)
    {
        define(m_params.theta_perp(), m_params.theta_par(), m_params.ub(), value);
    }

    // ====================================================================
    // Getters
    // ====================================================================

    /// The current parameters, for example to test param().capped().
    param_type param() const { return m_params; }

    RealType theta_perp() const { return m_params.theta_perp(); }   ///< Perpendicular thermal speed.
    RealType theta_par() const { return m_params.theta_par(); }     ///< Parallel thermal speed.
    point_type ub() const { return m_params.ub(); }                 ///< Field direction, as given.
    /// Velocity cap in thermal speeds; +infinity for no_cap().
    RealType max_normalized_velocity() const { return m_params.max_normalized_velocity(); }

    // ====================================================================
    // Random engine control
    // ====================================================================
    /// Reseed the internal engine and discard any variate cached from before, so that
    /// seed(s) followed by n draws gives the same values however many draws came first.
    void seed(int s)
    {
        applySeed(m_gen, s);
        reset();
    }

    /// Discard the normal variate cached by `std::normal_distribution`.  Called by seed()
    /// and by every define() overload.
    void reset() { m_standardNormal.reset(); }

    // ====================================================================
    // Sampling
    // ====================================================================
    /// Draw one velocity from the engine @p gen instead of the internal one.
    /// @return {vx, vy, vz} in the global frame.
    template <typename GeneratorType> result_type operator()(GeneratorType &gen)
    {
        return (*this)(gen, m_params);
    }

    /// Draw one velocity from the internal engine.
    /// @return {vx, vy, vz} in the global frame.
    result_type operator()()
    {
        return (*this)(m_gen);
    }

    /// Draw one velocity from @p gen with the parameters @p params instead of the stored
    /// ones.
    /// @throws std::invalid_argument if the field direction of @p params is numerically
    ///         zero; std::runtime_error after 10^6 consecutive rejections by a finite cap.
    template <typename GeneratorType>
    result_type operator()(GeneratorType &gen, const param_type &params)
    {
        const RealType sigmaScale= RealType(1.0) / std::sqrt(RealType(2.0));
        const int kMaxCapRejectTries= 1000000;
        RealType maxNormalizedVelocity= params.max_normalized_velocity();

        // Untruncated mode runs the core draw exactly once; capped mode redraws
        // until the sample falls inside the box.
        const bool capped= params.capped();
        const int maxTries= capped ? kMaxCapRejectTries : 1;

        int tries= 0;
        for (tries= 0; tries < maxTries; ++tries)
        {
            RealType vx= params.theta_perp() * sigmaScale * m_standardNormal(gen);
            RealType vy= params.theta_perp() * sigmaScale * m_standardNormal(gen);
            RealType vz= params.theta_par() * sigmaScale * m_standardNormal(gen);

            point_type point;
            point[0]= vx;
            point[1]= vy;
            point[2]= vz;

            if (capped && !withinNormalizedVelocityCap(point, params.theta_perp(),
                                                       params.theta_par(),
                                                       maxNormalizedVelocity))
            {
                continue;
            }

            // {0,0,1} is the identity rotation, so skip the transform for it.
            point_type ub= params.ub();
            if (approximatelyEqual(ub[0], RealType(0.0)) &&
                approximatelyEqual(ub[1], RealType(0.0)) &&
                approximatelyEqual(ub[2], RealType(1.0)))
            {
                return point;
            }
            return rotate_from_fieldAligned_frame(point, ub);
        }

        throw std::runtime_error("failed to sample velocity under max_normalized_velocity cap: "
                                 "rejection sampling exceeded internal attempt limit");
    }

    // ====================================================================
    // Geometry helper
    // ====================================================================
    /// Rotate @p local = {v_perp1, v_perp2, v_par} from field-aligned components to the
    /// global frame, for the field direction @p ub.  The perpendicular axes are the ones the
    /// sampler uses.  Throws std::invalid_argument if @p ub is numerically zero.
    point_type rotate_from_fieldAligned_frame(const point_type &local, const point_type &ub) const
    {
        const RealType kZero= RealType(0.0);
        const RealType kOne= RealType(1.0);

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

        point_type e1;
        e1[0]= e2[1] * e3[2] - e2[2] * e3[1];
        e1[1]= e2[2] * e3[0] - e2[0] * e3[2];
        e1[2]= e2[0] * e3[1] - e2[1] * e3[0];

        point_type point;
        point[0]= local[0] * e1[0] + local[1] * e2[0] + local[2] * e3[0];
        point[1]= local[0] * e1[1] + local[1] * e2[1] + local[2] * e3[1];
        point[2]= local[0] * e1[2] + local[1] * e2[2] + local[2] * e3[2];
        return point;
    }

    // ====================================================================
    // Comparison
    // ====================================================================
    /// True if the two samplers have equal parameters.  Engine state is not compared.
    bool operator==(const bi_maxwellian_distribution &other) const
    {
        return m_params == other.m_params;
    }

    /// Negation of operator==.
    bool operator!=(const bi_maxwellian_distribution &other) const { return !(*this == other); }

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

    static bool approximatelyEqual(RealType a, RealType b)
    {
        RealType scale= std::max({RealType(1.0), std::abs(a), std::abs(b)});
        return std::abs(a - b) <= std::numeric_limits<RealType>::epsilon() * scale;
    }

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
    std::normal_distribution<RealType> m_standardNormal;
};

#endif // _BI_MAXWELLIAN_DISTRIBUTION_H_
