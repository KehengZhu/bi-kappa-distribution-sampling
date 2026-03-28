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
 * Vx ~ N(0, theta_perp^2), Vy ~ N(0, theta_perp^2), Vz ~ N(0, theta_para^2).
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
        typedef DistributionType distribution_type;

        ParamType(RealType theta_perp= static_cast<RealType>(1.0),
                  RealType theta_para= static_cast<RealType>(1.0))
            : m_theta_perp(theta_perp), m_theta_para(theta_para)
        {
            if (m_theta_perp <= static_cast<RealType>(0.0) ||
                m_theta_para <= static_cast<RealType>(0.0))
            {
                throw std::invalid_argument("Bi-Maxwellian temperatures must be positive");
            }
        }

        RealType theta_perp() const { return m_theta_perp; }
        RealType theta_para() const { return m_theta_para; }

        bool operator==(const ParamType &other) const
        {
            return m_theta_perp == other.m_theta_perp && m_theta_para == other.m_theta_para;
        }

        bool operator!=(const ParamType &other) const { return !(*this == other); }

      private:
        RealType m_theta_perp;
        RealType m_theta_para;
    };

    BiMaxwellianDistribution() : m_params() {}

    explicit BiMaxwellianDistribution(RealType theta_perp,
                                      RealType theta_para= static_cast<RealType>(1.0))
        : m_params(theta_perp, theta_para)
    {
    }

    explicit BiMaxwellianDistribution(const ParamType &params) : m_params(params) {}

    void reset() {}

    template <typename GeneratorType> ResultType operator()(GeneratorType &gen)
    {
        return (*this)(gen, m_params);
    }

    template <typename GeneratorType>
    ResultType operator()(GeneratorType &gen, const ParamType &params)
    {
        std::normal_distribution<RealType> standardNormal(static_cast<RealType>(0.0),
                                                          static_cast<RealType>(1.0));

        RealType vx= params.theta_perp() * standardNormal(gen);
        RealType vy= params.theta_perp() * standardNormal(gen);
        RealType vz= params.theta_para() * standardNormal(gen);

        Point3D point;
        point[0]= vx;
        point[1]= vy;
        point[2]= vz;
        return point;
    }

    ParamType param() const { return m_params; }
    void param(const ParamType &params) { m_params= params; }

    RealType theta_perp() const { return m_params.theta_perp(); }
    RealType theta_para() const { return m_params.theta_para(); }

    bool operator==(const BiMaxwellianDistribution &other) const
    {
        return m_params == other.m_params;
    }

    bool operator!=(const BiMaxwellianDistribution &other) const { return !(*this == other); }

  private:
    ParamType m_params;
};

#endif // _BI_MAXWELLIAN_DISTRIBUTION_H_
