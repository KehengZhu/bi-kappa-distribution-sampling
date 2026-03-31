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
    typedef ResultType result_type;

    /**
     * @brief Parameter class to store distribution parameters
     */
    class ParamType
    {
      public:
        typedef BiKappaDistribution<RealType> DistributionType;
        typedef DistributionType distribution_type;

        ParamType(RealType kappa= RealType(2.0), RealType theta_perp= RealType(1.0),
              RealType theta_para= RealType(1.0),
              Point3D ub= Point3D({RealType(0.0), RealType(0.0), RealType(0.0)}))
            : m_kappa(kappa), m_theta_perp(theta_perp), m_theta_para(theta_para), m_ub(ub)
        {
            if (m_kappa <= RealType(0.5))
            {
                throw std::invalid_argument("kappa must be > 0.5 for normalization");
            }
            if (m_theta_perp <= RealType(0.0) || m_theta_para <= RealType(0.0))
            {
                throw std::invalid_argument("Temperature parameters must be positive");
            }
        }

        RealType kappa() const { return m_kappa; }
        RealType theta_perp() const { return m_theta_perp; }
        RealType theta_para() const { return m_theta_para; }
        Point3D ub() const { return m_ub; }

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
        Point3D m_ub;
    };

    // Constructors
    BiKappaDistribution()
        : m_params(), m_distX1(RealType(1.5), RealType(1.0)),
          m_distX2(RealType(1.0), RealType(1.0)), m_distCosTheta(RealType(-1.0), RealType(1.0)),
          m_distPhi(RealType(0.0), RealType(2.0) * piValue())
    {
    }

    explicit BiKappaDistribution(RealType kappa, RealType theta_perp= RealType(1.0),
                                 RealType theta_para= RealType(1.0),
                                 Point3D ub= Point3D({RealType(0.0), RealType(0.0), RealType(0.0)}))
        : m_params(kappa, theta_perp, theta_para, ub), m_distX1(RealType(1.5), RealType(1.0)),
          m_distX2(RealType(1.0), RealType(1.0)), m_distCosTheta(RealType(-1.0), RealType(1.0)),
          m_distPhi(RealType(0.0), RealType(2.0) * piValue())
    {
    }

    explicit BiKappaDistribution(const ParamType &params)
        : m_params(params), m_distX1(RealType(1.5), RealType(1.0)),
          m_distX2(RealType(1.0), RealType(1.0)), m_distCosTheta(RealType(-1.0), RealType(1.0)),
          m_distPhi(RealType(0.0), RealType(2.0) * piValue())
    {
    }

    typedef ParamType param_type;

    void reset()
    {
        m_distX1.reset();
        m_distX2.reset();
        m_distCosTheta.reset();
        m_distPhi.reset();
    }

    // Generate a sample
    template <typename GeneratorType> ResultType operator()(GeneratorType &gen)
    {
        return (*this)(gen, m_params);
    }

    template <typename GeneratorType>
    ResultType operator()(GeneratorType &gen, const ParamType &params)
    {
        const RealType kHalf= RealType(0.5);
        const RealType kOne= RealType(1.0);

        RealType kappa= params.kappa();
        RealType theta_perp= params.theta_perp();
        RealType theta_para= params.theta_para();
        Point3D ub= params.ub();

        typename std::gamma_distribution<RealType>::param_type distX2Param(kappa - kHalf, kOne);

        // Step 1: Sample radii
        RealType x1= m_distX1(gen);
        RealType x2= m_distX2(gen, distX2Param);
        RealType r= std::sqrt(x1 / x2);

        // Step 2: Sample directions uniformly on S^2
        RealType cosTheta= m_distCosTheta(gen);
        RealType sinTheta= std::sqrt(std::max(RealType(0.0), kOne - cosTheta * cosTheta));
        RealType phi= m_distPhi(gen);

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

        // Step 5: rotate from local field-aligned frame to global frame
        // only when ub is explicitly set (non-zero).
        RealType ubNorm2= ub[0] * ub[0] + ub[1] * ub[1] + ub[2] * ub[2];
        if (ubNorm2 > RealType(0.0))
        {
            point = rotate_from_fieldAligned_frame(point, ub);
        }
        
        return point;
    }

    // Parameter access
    ParamType param() const { return m_params; }
    void param(const ParamType &params) { m_params= params; }

    RealType kappa() const { return m_params.kappa(); }
    RealType theta_perp() const { return m_params.theta_perp(); }
    RealType theta_para() const { return m_params.theta_para(); }

    // Helper function: transform from velocity components (v_perp, v_par) to
    // vector velocity (vx,vy,vz) given parallel direction ub
    Point3D rotate_from_fieldAligned_frame(const Point3D &local, const Point3D &ub) const
    {
        const RealType kZero= RealType(0.0);
        const RealType kOne= RealType(1.0);

        // Normalize ub to get e3
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

        // Build e2 perpendicular to e3
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

        // e1 = e2 x e3
        Point3D e1;
        e1[0]= e2[1] * e3[2] - e2[2] * e3[1];
        e1[1]= e2[2] * e3[0] - e2[0] * e3[2];
        e1[2]= e2[0] * e3[1] - e2[1] * e3[0];

        // Transform local coordinates to Cartesian coordinates
        Point3D point;
        point[0]= local[0] * e1[0] + local[1] * e2[0] + local[2] * e3[0];
        point[1]= local[0] * e1[1] + local[1] * e2[1] + local[2] * e3[1];
        point[2]= local[0] * e1[2] + local[1] * e2[2] + local[2] * e3[2];

        return point;
    }

    // Comparison operators
    bool operator==(const BiKappaDistribution &other) const { return m_params == other.m_params; }

    bool operator!=(const BiKappaDistribution &other) const { return !(*this == other); }

  private:
    static RealType piValue() { return RealType(3.14159265358979323846264338327950288L); }

    ParamType m_params;
    std::gamma_distribution<RealType> m_distX1;
    std::gamma_distribution<RealType> m_distX2;
    std::uniform_real_distribution<RealType> m_distCosTheta;
    std::uniform_real_distribution<RealType> m_distPhi;
};

#endif // _BI_KAPPA_DISTRIBUTION_H_
