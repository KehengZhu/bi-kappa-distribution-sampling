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
            if (m_maxNormalizedVelocity <= RealType(0.0))
            {
                throw std::invalid_argument("max_normalized_velocity must be positive");
            }
        }

        RealType theta_perp() const { return m_theta_perp; }
        RealType theta_par() const { return m_theta_par; }
        point_type ub() const { return m_ub; }
        RealType max_normalized_velocity() const { return m_maxNormalizedVelocity; }

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

    // ====================================================================
    // Construction
    // ====================================================================
    bi_maxwellian_distribution() : m_params(), m_gen(std::random_device{}()), m_standardNormal(RealType(0.0), RealType(1.0)) {}

    explicit bi_maxwellian_distribution(RealType theta_perp, RealType theta_par= RealType(1.0),
                                        point_type ub= point_type({RealType(0.0), RealType(0.0),
                                                                   RealType(1.0)}),
                                        RealType max_normalized_velocity= RealType(20.0))
        : m_params(theta_perp, theta_par, ub), m_gen(std::random_device{}()), m_standardNormal(RealType(0.0), RealType(1.0))
    {
        define(theta_perp, theta_par, ub, max_normalized_velocity);
    }

    explicit bi_maxwellian_distribution(const param_type &params)
        : m_params(params), m_gen(std::random_device{}()), m_standardNormal(RealType(0.0), RealType(1.0))
    {
        define(params);
    }

    // ====================================================================
    // Bulk configuration
    // ====================================================================
    void define(const param_type &params)
    {
        m_params= params;
        reset(); // just to clear cache
    }

    void define(RealType theta_perp, RealType theta_par= RealType(1.0),
                point_type ub= point_type({RealType(0.0), RealType(0.0), RealType(1.0)}),
                RealType max_normalized_velocity= RealType(20.0))
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

    // ====================================================================
    // Individual parameter setters
    // ====================================================================
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

    void max_normalized_velocity(RealType value)
    {
        define(m_params.theta_perp(), m_params.theta_par(), m_params.ub(), value);
    }

    // ====================================================================
    // Parameter access
    // ====================================================================
    param_type param() const { return m_params; }
    void param(const param_type &params) { define(params); }

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

        int tries= 0;
        for (tries= 0; tries < kMaxCapRejectTries; ++tries)
        {
            RealType vx= params.theta_perp() * sigmaScale * m_standardNormal(gen);
            RealType vy= params.theta_perp() * sigmaScale * m_standardNormal(gen);
            RealType vz= params.theta_par() * sigmaScale * m_standardNormal(gen);

            point_type point;
            point[0]= vx;
            point[1]= vy;
            point[2]= vz;

            if (!withinNormalizedVelocityCap(point, params.theta_perp(), params.theta_par(),
                                             maxNormalizedVelocity))
            {
                continue;
            }

            // Rotate only when ub is explicitly set (non-zero).
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
