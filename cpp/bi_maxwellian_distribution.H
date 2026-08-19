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
 * @note **The cap is off by default**, so the default target law is the
 * bi-Maxwellian distribution itself.  Passing a finite `max_normalized_velocity`
 * opts in to a *different*, bounded target: the bi-Maxwellian law conditioned on
 * the component-wise box.  See no_cap() and param_type::capped().  This matches
 * bi_kappa_distribution; both samplers default to the law they are named after.
 *
 * @see \ref usage for more recipes (seeding, external RNG, reconfiguration).
 */
template <typename RealType= double> class bi_maxwellian_distribution
{
  public:
    using result_type= std::array<RealType, 3>;
    using point_type= result_type;

    class param_type
    {
      public:
        using distribution_type= bi_maxwellian_distribution<RealType>;

        // The cap defaults to +infinity, i.e. off -- the same value no_cap() returns.
        // It is spelled out here rather than as no_cap() because the enclosing class is
        // still incomplete while this nested default argument is parsed.
        param_type(RealType theta_perp= RealType(1.0), RealType theta_par= RealType(1.0),
                   point_type ub= point_type({RealType(0.0), RealType(0.0), RealType(1.0)}),
                   RealType max_normalized_velocity=
                       std::numeric_limits<RealType>::infinity())
            : m_theta_perp(theta_perp), m_theta_par(theta_par), m_ub(ub),
              m_maxNormalizedVelocity(max_normalized_velocity)
        {
            if (m_theta_perp <= RealType(0.0) || m_theta_par <= RealType(0.0))
            {
                throw std::invalid_argument("Bi-Maxwellian temperatures must be positive");
            }
            // A positive finite value selects the capped mode.  no_cap() (= +inf)
            // selects the untruncated mode, in which no acceptance-rejection loop
            // runs at all and samples follow the full bi-Maxwellian law.
            if (!(m_maxNormalizedVelocity > RealType(0.0)))
            {
                throw std::invalid_argument(
                    "max_normalized_velocity must be positive, or no_cap() for the "
                    "untruncated distribution");
            }
        }

        RealType theta_perp() const { return m_theta_perp; }
        RealType theta_par() const { return m_theta_par; }
        point_type ub() const { return m_ub; }
        RealType max_normalized_velocity() const { return m_maxNormalizedVelocity; }

        /// True when a finite component-wise cap is in force, i.e. when the
        /// sampled law is the bi-Maxwellian distribution conditioned on the cap
        /// event rather than the bi-Maxwellian distribution itself.
        bool capped() const { return std::isfinite(m_maxNormalizedVelocity); }

        bool operator==(const param_type &other) const
        {
            return m_theta_perp == other.m_theta_perp && m_theta_par == other.m_theta_par &&
                   m_ub == other.m_ub && m_maxNormalizedVelocity == other.m_maxNormalizedVelocity;
        }

        bool operator!=(const param_type &other) const { return !(*this == other); }

      private:
        RealType m_theta_perp;
        RealType m_theta_par;
        point_type m_ub;
        RealType m_maxNormalizedVelocity;
    };

    /// Sentinel selecting the untruncated sampler, and **the default**.  In this
    /// mode each draw is three Gaussians with no outer acceptance-rejection loop,
    /// and the sampled law is the full bi-Maxwellian distribution.
    ///
    /// Naming it explicitly is equivalent to omitting the argument, and is worth
    /// doing where the target law matters to a reader:
    ///
    ///     bi_maxwellian_distribution<double> d(1.0, 2.0, {0,0,1},
    ///                                          bi_maxwellian_distribution<double>::no_cap());
    ///
    /// With a finite cap the sampled law is instead that distribution conditioned
    /// on all normalized components lying inside the component-wise box -- a
    /// different, bounded law with its own normalization.  The cap is tested in the
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
    bi_maxwellian_distribution() : m_params(), m_gen(std::random_device{}()), m_standardNormal(RealType(0.0), RealType(1.0)) {}

    explicit bi_maxwellian_distribution(RealType theta_perp, RealType theta_par= RealType(1.0),
                                        point_type ub= point_type({RealType(0.0), RealType(0.0),
                                                                   RealType(1.0)}),
                                        RealType max_normalized_velocity= no_cap())
        : m_params(theta_perp, theta_par, ub, max_normalized_velocity),
          m_gen(std::random_device{}()), m_standardNormal(RealType(0.0), RealType(1.0))
    {
        define(theta_perp, theta_par, ub, max_normalized_velocity);
    }

    explicit bi_maxwellian_distribution(const param_type &params)
        : m_params(params), m_gen(std::random_device{}()), m_standardNormal(RealType(0.0), RealType(1.0))
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
     * @param theta_perp  Thermal speed of the two components transverse to @p ub,
     *                > 0.  Each is N(0, theta_perp^2/2) -- a speed, not a
     *                temperature: theta = sqrt(2*k_B*T/m).
     * @param theta_par   Thermal speed along @p ub, > 0.
     * @param ub      Magnetic-field direction {bx,by,bz} in the global frame.
     *                Need not be normalized; must be non-zero.
     * @param max_normalized_velocity  Velocity cap in thermal speeds, or no_cap()
     *                (the default) for the untruncated bi-Maxwellian law.  A finite
     *                value selects a *different*, bounded target law -- see
     *                \ref usage_bikappa_cap.  Must be > 0; passing 0 to mean
     *                "no cap" throws.
     *
     * @note To seed reproducibly, either call seed() afterwards or use the
     * overload taking a trailing int: >= 0 is reproducible, < 0 reseeds from
     * std::random_device.
     *
     * @throws std::invalid_argument if either theta is non-positive,
     *         @p max_normalized_velocity is non-positive, or @p ub is the zero vector.
     *
     * @see \ref parameters for the full table, \ref choosing to compare samplers.
     */
    void define(RealType theta_perp, RealType theta_par= RealType(1.0),
                point_type ub= point_type({RealType(0.0), RealType(0.0), RealType(1.0)}),
                RealType max_normalized_velocity= no_cap())
    {
        define(param_type(theta_perp, theta_par, ub, max_normalized_velocity));
    }

    void define(const param_type &params, int seed) { applySeed(m_gen, seed); define(params); }

    void define(RealType theta_perp, RealType theta_par,
                point_type ub, RealType max_normalized_velocity, int seed)
    {
        applySeed(m_gen, seed);
        define(theta_perp, theta_par, ub, max_normalized_velocity);
    }

    // Individual parameter setters
    void theta_perp(RealType value)
    {
        define(value, m_params.theta_par(), m_params.ub(), m_params.max_normalized_velocity());
    }

    void theta_par(RealType value)
    {
        define(m_params.theta_perp(), value, m_params.ub(), m_params.max_normalized_velocity());
    }

    void ub(const point_type &value)
    {
        define(m_params.theta_perp(), m_params.theta_par(), value,
               m_params.max_normalized_velocity());
    }

    /// Unlike the other setters, this one changes *which distribution* is sampled,
    /// not merely its parameters.  Pass no_cap() to return to the bi-Maxwellian law
    /// itself; passing 0 to mean "no cap" throws, as every non-positive value does.
    void max_normalized_velocity(RealType value)
    {
        define(m_params.theta_perp(), m_params.theta_par(), m_params.ub(), value);
    }

    // ====================================================================
    // Getters
    // ====================================================================

    // Bulk getter
    param_type param() const { return m_params; }

    // Individual parameter getter
    RealType theta_perp() const { return m_params.theta_perp(); }
    RealType theta_par() const { return m_params.theta_par(); }
    point_type ub() const { return m_params.ub(); }
    RealType max_normalized_velocity() const { return m_params.max_normalized_velocity(); }

    // ====================================================================
    // Random engine control
    // ====================================================================
    void seed(int s) { applySeed(m_gen, s); }

    void reset() { m_standardNormal.reset(); }

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
    // Transform local velocity components to global Cartesian coordinates given the parallel
    // direction ub.
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
    bool operator==(const bi_maxwellian_distribution &other) const
    {
        return m_params == other.m_params;
    }

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
