#ifndef _BI_KAPPA_DISTRIBUTION_H_
#define _BI_KAPPA_DISTRIBUTION_H_

#include <array>
#include <algorithm>
#include <cmath>
#include <random>
#include <stdexcept>

/**
 * @brief Bi-Kappa distribution class for sampling 3D velocity vectors.
 *
 * This class follows the interface pattern of standard library distributions
 * like std::uniform_real_distribution<T>.
 *
 * @tparam RealType The floating-point type (e.g., double, float)
 */
template <typename RealType= double> class BiKappaDistribution
{
  public:
    typedef std::array<RealType, 3> Point3D;
    typedef Point3D ResultType;

    /**
     * @brief Parameter class to store distribution parameters
     */
    class ParamType
    {
      public:
        typedef BiKappaDistribution<RealType> DistributionType;
        typedef DistributionType distribution_type;

        ParamType(RealType kappa= static_cast<RealType>(2.0),
                  RealType theta_perp= static_cast<RealType>(1.0),
                  RealType theta_para= static_cast<RealType>(1.0))
            : m_kappa(kappa), m_theta_perp(theta_perp), m_theta_para(theta_para)
        {
            if (m_kappa <= static_cast<RealType>(0.5))
            {
                throw std::invalid_argument("kappa must be > 0.5 for normalization");
            }
            if (m_theta_perp <= static_cast<RealType>(0.0) ||
                m_theta_para <= static_cast<RealType>(0.0))
            {
                throw std::invalid_argument("Temperature parameters must be positive");
            }
        }

        RealType kappa() const { return m_kappa; }
        RealType theta_perp() const { return m_theta_perp; }
        RealType theta_para() const { return m_theta_para; }

        bool operator==(const ParamType &other) const
        {
            return m_kappa == other.m_kappa && m_theta_perp == other.m_theta_perp &&
                   m_theta_para == other.m_theta_para;
        }

        bool operator!=(const ParamType &other) const { return !(*this == other); }

      private:
        RealType m_kappa;
        RealType m_theta_perp;
        RealType m_theta_para;
    };

    // Constructors
    BiKappaDistribution() : m_params() {}

    explicit BiKappaDistribution(RealType kappa, RealType theta_perp= static_cast<RealType>(1.0),
                                 RealType theta_para= static_cast<RealType>(1.0))
        : m_params(kappa, theta_perp, theta_para)
    {
    }

    explicit BiKappaDistribution(const ParamType &params) : m_params(params) {}

    // Reset (no-op for this distribution, but required by standard interface)
    void reset() {}

    // Generate a sample
    template <typename GeneratorType> ResultType operator()(GeneratorType &gen)
    {
        return (*this)(gen, m_params);
    }

    template <typename GeneratorType>
    ResultType operator()(GeneratorType &gen, const ParamType &params)
    {
        RealType kappa= params.kappa();
        RealType theta_perp= params.theta_perp();
        RealType theta_para= params.theta_para();
        RealType pi= static_cast<RealType>(3.14159265358979323846264338327950288L);

        // Gamma distributions for the radius (R = sqrt(X1/X2))
        std::gamma_distribution<RealType> distX1(static_cast<RealType>(1.5),
                                                 static_cast<RealType>(1.0));
        std::gamma_distribution<RealType> distX2(kappa - static_cast<RealType>(0.5),
                                                 static_cast<RealType>(1.0));

        // Uniform distributions for spherical directions
        std::uniform_real_distribution<RealType> distCosTheta(static_cast<RealType>(-1.0),
                                                              static_cast<RealType>(1.0));
        std::uniform_real_distribution<RealType> distPhi(static_cast<RealType>(0.0),
                                                         static_cast<RealType>(2.0) * pi);

        // Step 1: Sample radii
        RealType x1= distX1(gen);
        RealType x2= distX2(gen);
        RealType r= std::sqrt(x1 / x2);

        // Step 2: Sample directions uniformly on S^2
        RealType cosTheta= distCosTheta(gen);
        RealType sinTheta= std::sqrt(
            std::max(static_cast<RealType>(0.0), static_cast<RealType>(1.0) - cosTheta * cosTheta));
        RealType phi= distPhi(gen);

        // Step 3: Compute unit vectors (u)
        RealType ux= r * sinTheta * std::cos(phi);
        RealType uy= r * sinTheta * std::sin(phi);
        RealType uz= r * cosTheta;

        // Step 4: Map back to v-space and return as array
        RealType sqrtKappa= std::sqrt(kappa);

        Point3D point;
        point[0]= sqrtKappa * theta_perp * ux;
        point[1]= sqrtKappa * theta_perp * uy;
        point[2]= sqrtKappa * theta_para * uz;
        return point;
    }

    // Parameter access
    ParamType param() const { return m_params; }
    void param(const ParamType &params) { m_params= params; }

    RealType kappa() const { return m_params.kappa(); }
    RealType theta_perp() const { return m_params.theta_perp(); }
    RealType theta_para() const { return m_params.theta_para(); }

    // Comparison operators
    bool operator==(const BiKappaDistribution &other) const { return m_params == other.m_params; }

    bool operator!=(const BiKappaDistribution &other) const { return !(*this == other); }

  private:
    ParamType m_params;
};

#endif // _BI_KAPPA_DISTRIBUTION_H_
