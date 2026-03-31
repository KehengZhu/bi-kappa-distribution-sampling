#ifndef _BI_MAXWELLIAN_DISTRIBUTION_H_
#define _BI_MAXWELLIAN_DISTRIBUTION_H_

#include <array>
#include <cmath>
#include <random>
#include <stdexcept>

/**
 * @brief Bi-Maxwellian distribution for 3D velocity sampling.
 *
 * Components are independent Gaussians with anisotropic thermal widths:
 * Vx ~ N(0, theta_perp^2/2), Vy ~ N(0, theta_perp^2/2), Vz ~ N(0, theta_para^2/2).
 */
template <typename RealType= double> class BiMaxwellianDistribution
{
  public:
    typedef std::array<RealType, 3> Point3D;
    typedef Point3D ResultType;
    typedef ResultType result_type;

    class ParamType
    {
      public:
        typedef BiMaxwellianDistribution<RealType> DistributionType;

        ParamType(RealType theta_perp= RealType(1.0), RealType theta_para= RealType(1.0),
                  Point3D ub= Point3D({RealType(0.0), RealType(0.0), RealType(0.0)}))
            : m_theta_perp(theta_perp), m_theta_para(theta_para), m_ub(ub)
        {
            if (m_theta_perp <= RealType(0.0) || m_theta_para <= RealType(0.0))
            {
                throw std::invalid_argument("Bi-Maxwellian temperatures must be positive");
            }
        }

        RealType theta_perp() const { return m_theta_perp; }
        RealType theta_para() const { return m_theta_para; }
        Point3D ub() const { return m_ub; }

        bool operator==(const ParamType &other) const
        {
            return m_theta_perp == other.m_theta_perp && m_theta_para == other.m_theta_para &&
                   m_ub == other.m_ub;
        }

        bool operator!=(const ParamType &other) const { return !(*this == other); }

      private:
        RealType m_theta_perp;
        RealType m_theta_para;
        Point3D m_ub;
    };
    typedef ParamType param_type;

    BiMaxwellianDistribution()
        : m_params(), m_standardNormal(RealType(0.0), RealType(1.0))
    {
    }

    explicit BiMaxwellianDistribution(RealType theta_perp, RealType theta_para= RealType(1.0),
                                      Point3D ub= Point3D({RealType(0.0), RealType(0.0),
                                                           RealType(0.0)}))
        : m_params(theta_perp, theta_para, ub), m_standardNormal(RealType(0.0), RealType(1.0))
    {
    }

    explicit BiMaxwellianDistribution(const ParamType &params)
        : m_params(params), m_standardNormal(RealType(0.0), RealType(1.0))
    {
    }

    void reset() { m_standardNormal.reset(); }

    template <typename GeneratorType> ResultType operator()(GeneratorType &gen)
    {
        return (*this)(gen, m_params);
    }

    template <typename GeneratorType>
    ResultType operator()(GeneratorType &gen, const ParamType &params)
    {
        const RealType kOne= RealType(1.0);
        const RealType sigmaScale= kOne / std::sqrt(RealType(2.0));

        RealType vx= params.theta_perp() * sigmaScale * m_standardNormal(gen);
        RealType vy= params.theta_perp() * sigmaScale * m_standardNormal(gen);
        RealType vz= params.theta_para() * sigmaScale * m_standardNormal(gen);

        Point3D point;
        point[0]= vx;
        point[1]= vy;
        point[2]= vz;

        // Rotate only when ub is explicitly set (non-zero).
        Point3D ub= params.ub();
        RealType ubNorm2= ub[0] * ub[0] + ub[1] * ub[1] + ub[2] * ub[2];
        if (ubNorm2 > RealType(0.0))
        {
            return rotate_from_fieldAligned_frame(point, ub);
        }
        return point;
    }

    ParamType param() const { return m_params; }
    void param(const ParamType &params) { m_params= params; }

    RealType theta_perp() const { return m_params.theta_perp(); }
    RealType theta_para() const { return m_params.theta_para(); }

    // Helper function: transform local velocity components to global Cartesian coordinates
    // given parallel direction ub.
    Point3D rotate_from_fieldAligned_frame(const Point3D &local, const Point3D &ub) const
    {
        const RealType kZero= RealType(0.0);
        const RealType kOne= RealType(1.0);

        RealType ubNorm2= ub[0] * ub[0] + ub[1] * ub[1] + ub[2] * ub[2];
        if (ubNorm2 <= kZero)
        {
            throw std::invalid_argument("ub must be a non-zero vector");
        }

        RealType ubNorm= std::sqrt(ubNorm2);

        Point3D e3;
        e3[0]= ub[0] / ubNorm;
        e3[1]= ub[1] / ubNorm;
        e3[2]= ub[2] / ubNorm;

        Point3D e2({kOne, kOne, kOne});

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

        Point3D e1;
        e1[0]= e2[1] * e3[2] - e2[2] * e3[1];
        e1[1]= e2[2] * e3[0] - e2[0] * e3[2];
        e1[2]= e2[0] * e3[1] - e2[1] * e3[0];

        Point3D point;
        point[0]= local[0] * e1[0] + local[1] * e2[0] + local[2] * e3[0];
        point[1]= local[0] * e1[1] + local[1] * e2[1] + local[2] * e3[1];
        point[2]= local[0] * e1[2] + local[1] * e2[2] + local[2] * e3[2];
        return point;
    }

    bool operator==(const BiMaxwellianDistribution &other) const
    {
        return m_params == other.m_params;
    }

    bool operator!=(const BiMaxwellianDistribution &other) const { return !(*this == other); }

  private:
    ParamType m_params;
        std::normal_distribution<RealType> m_standardNormal;
};

#endif // _BI_MAXWELLIAN_DISTRIBUTION_H_
