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
                   point_type ub= point_type({RealType(0.0), RealType(0.0), RealType(1.0)}))
            : m_theta_perp(theta_perp), m_theta_par(theta_par), m_ub(ub)
        {
            if (m_theta_perp <= RealType(0.0) || m_theta_par <= RealType(0.0))
            {
                throw std::invalid_argument("Bi-Maxwellian temperatures must be positive");
            }
        }

        RealType theta_perp() const { return m_theta_perp; }
        RealType theta_par() const { return m_theta_par; }
        point_type ub() const { return m_ub; }

        bool operator==(const param_type &other) const
        {
            return m_theta_perp == other.m_theta_perp && m_theta_par == other.m_theta_par &&
                   m_ub == other.m_ub;
        }

        bool operator!=(const param_type &other) const { return !(*this == other); }

      private:
        RealType m_theta_perp;
        RealType m_theta_par;
        point_type m_ub;
    };

    bi_maxwellian_distribution() : m_params(), m_standardNormal(RealType(0.0), RealType(1.0)) {}

    explicit bi_maxwellian_distribution(RealType theta_perp, RealType theta_par= RealType(1.0),
                                        point_type ub= point_type({RealType(0.0), RealType(0.0),
                                                                   RealType(1.0)}))
        : m_params(theta_perp, theta_par, ub), m_standardNormal(RealType(0.0), RealType(1.0))
    {
    }

    explicit bi_maxwellian_distribution(const param_type &params)
        : m_params(params), m_standardNormal(RealType(0.0), RealType(1.0))
    {
    }

    void reset() { m_standardNormal.reset(); }

    template <typename GeneratorType> result_type operator()(GeneratorType &gen)
    {
        return (*this)(gen, m_params);
    }

    template <typename GeneratorType>
    result_type operator()(GeneratorType &gen, const param_type &params)
    {
        const RealType sigmaScale= RealType(1.0) / std::sqrt(RealType(2.0));

        RealType vx= params.theta_perp() * sigmaScale * m_standardNormal(gen);
        RealType vy= params.theta_perp() * sigmaScale * m_standardNormal(gen);
        RealType vz= params.theta_par() * sigmaScale * m_standardNormal(gen);

        point_type point;
        point[0]= vx;
        point[1]= vy;
        point[2]= vz;

        // Rotate only when ub is explicitly set (non-zero).
        point_type ub= params.ub();
        if (approximatelyEqual(ub[0], RealType(0.0)) && approximatelyEqual(ub[1], RealType(0.0)) &&
            approximatelyEqual(ub[2], RealType(1.0)))
        {
            return point;
        }
        return rotate_from_fieldAligned_frame(point, ub);
    }

    param_type param() const { return m_params; }
    void param(const param_type &params) { m_params= params; }

    RealType theta_perp() const { return m_params.theta_perp(); }
    RealType theta_par() const { return m_params.theta_par(); }

    // Helper function: transform local velocity components to global Cartesian coordinates
    // given parallel direction ub.
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

    bool operator==(const bi_maxwellian_distribution &other) const
    {
        return m_params == other.m_params;
    }

    bool operator!=(const bi_maxwellian_distribution &other) const { return !(*this == other); }

  private:
    static bool approximatelyEqual(RealType a, RealType b)
    {
        RealType scale= std::max({RealType(1.0), std::abs(a), std::abs(b)});
        return std::abs(a - b) <= std::numeric_limits<RealType>::epsilon() * scale;
    }

    param_type m_params;
    std::normal_distribution<RealType> m_standardNormal;
};

#endif // _BI_MAXWELLIAN_DISTRIBUTION_H_
