#ifndef _BI_KAPPA_DISTRIBUTION_H_
#define _BI_KAPPA_DISTRIBUTION_H_

#include <array>
#include <algorithm>
#include <cmath>
#include <limits>
#include <random>
#include <stdexcept>

/**
 * @brief Bi-Kappa distribution class for sampling 3D velocity vectors.
 *
 * This class follows the interface pattern of standard library distributions
 * like std::uniform_real_distribution<T>.
 *
 * @code{.cpp}
 * bi_kappa_distribution<double> dist;
 * dist.define(2.0, 1.0, 2.0, {0, 0, 1}, 20.0, 12345); // kappa, theta_perp, theta_par, ub, cap, seed
 * auto v = dist();                                     // {vx, vy, vz}, global frame
 * @endcode
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

        param_type(RealType kappa= RealType(2.0), RealType theta_perp= RealType(1.0),
                   RealType theta_par= RealType(1.0),
                   point_type ub= point_type({RealType(0.0), RealType(0.0), RealType(1.0)}),
                   RealType max_normalized_velocity= RealType(20.0))
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

    /// Sentinel selecting the untruncated sampler.  Pass as
    /// max_normalized_velocity to disable the component-wise cap entirely:
    ///
    ///     bi_kappa_distribution<double> d(2.0, 1.0, 2.0, {0,0,1},
    ///                                     bi_kappa_distribution<double>::no_cap());
    ///
    /// In this mode the draw uses a fixed sequence of high-level variates with no
    /// outer acceptance-rejection loop, and the sampled law is the full bi-Kappa
    /// distribution.  With a finite cap the sampled law is instead that
    /// distribution conditioned on all normalized components lying inside the box.
    static constexpr RealType no_cap()
    {
        return std::numeric_limits<RealType>::infinity();
    }

    // ====================================================================
    // Construction
    // ====================================================================
    bi_kappa_distribution()
        : m_params(), m_gen(std::random_device{}()), m_distX1(RealType(1.5), RealType(1.0)),
          m_distX2(RealType(1.0), RealType(1.0)), m_distCosTheta(RealType(-1.0), RealType(1.0)),
          m_distPhi(RealType(0.0), RealType(2.0) * piValue())
    {
    }

    explicit bi_kappa_distribution(RealType kappa, RealType theta_perp= RealType(1.0),
                                   RealType theta_par= RealType(1.0),
                                   point_type ub= point_type({RealType(0.0), RealType(0.0),
                                                              RealType(1.0)}),
                                   RealType max_normalized_velocity= RealType(20.0))
        : m_params(kappa, theta_perp, theta_par, ub, max_normalized_velocity),
          m_gen(std::random_device{}()),
          m_distX1(RealType(1.5), RealType(1.0)), m_distX2(RealType(1.0), RealType(1.0)),
          m_distCosTheta(RealType(-1.0), RealType(1.0)),
          m_distPhi(RealType(0.0), RealType(2.0) * piValue())
    {
        define(kappa, theta_perp, theta_par, ub, max_normalized_velocity);
    }

    explicit bi_kappa_distribution(const param_type &params)
        : m_params(params), m_gen(std::random_device{}()),
          m_distX1(RealType(1.5), RealType(1.0)),
          m_distX2(RealType(1.0), RealType(1.0)), m_distCosTheta(RealType(-1.0), RealType(1.0)),
          m_distPhi(RealType(0.0), RealType(2.0) * piValue())
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

    void define(RealType kappa, RealType theta_perp= RealType(1.0),
                RealType theta_par= RealType(1.0),
                point_type ub= point_type({RealType(0.0), RealType(0.0), RealType(1.0)}),
                RealType max_normalized_velocity= RealType(20.0))
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
    void seed(int s) { applySeed(m_gen, s); }

    void reset()
    {
        m_distX1.reset();
        m_distX2.reset();
        m_distCosTheta.reset();
        m_distPhi.reset();
    }

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

        typename std::gamma_distribution<RealType>::param_type distX2Param(kappa - kHalf, kOne);
        RealType sqrtKappa= std::sqrt(kappa);

        // Untruncated mode runs the core mapping exactly once; capped mode redraws
        // until the sample falls inside the box.
        const bool capped= params.capped();
        const int maxTries= capped ? kMaxCapRejectTries : 1;

        int tries= 0;
        for (tries= 0; tries < maxTries; ++tries)
        {
            // Step 1: Sample radii.
            // Take the two square roots separately rather than forming the ratio
            // first: T = x1/x2 overflows for x2 near the underflow limit even where
            // R = sqrt(T) is perfectly representable, which discards valid draws at
            // small kappa (shape(x2) = kappa - 1/2 << 1).
            RealType x1= m_distX1(gen);
            RealType x2= m_distX2(gen, distX2Param);
            RealType r= std::sqrt(x1) / std::sqrt(x2);

            // Step 2: Sample directions uniformly on S^2
            RealType cosTheta= m_distCosTheta(gen);
            RealType sinTheta= std::sqrt(std::max(RealType(0.0), kOne - cosTheta * cosTheta));
            RealType phi= m_distPhi(gen);

            // Step 3: Compute unit vectors (u)
            RealType ux= r * sinTheta * std::cos(phi);
            RealType uy= r * sinTheta * std::sin(phi);
            RealType uz= r * cosTheta;

            // Step 4: Map back to v-space and apply velocity cap in local frame
            point_type point;
            point[0]= sqrtKappa * theta_perp * ux;
            point[1]= sqrtKappa * theta_perp * uy;
            point[2]= sqrtKappa * theta_par * uz;

            if (capped &&
                !withinNormalizedVelocityCap(point, theta_perp, theta_par, maxNormalizedVelocity))
            {
                continue;
            }

            if (approximatelyEqual(ub[0], RealType(0.0)) &&
                approximatelyEqual(ub[1], RealType(0.0)) &&
                approximatelyEqual(ub[2], RealType(1.0)))
            {
                return point;
            }

            // Step 5: rotate from local field-aligned frame to global frame
            return rotate_from_fieldAligned_frame(point, ub);
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
    std::gamma_distribution<RealType> m_distX1;
    std::gamma_distribution<RealType> m_distX2;
    std::uniform_real_distribution<RealType> m_distCosTheta;
    std::uniform_real_distribution<RealType> m_distPhi;
};

#endif // _BI_KAPPA_DISTRIBUTION_H_
