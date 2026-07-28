#ifndef _GENERAL_VELOCITY_GENERATOR_H_
#define _GENERAL_VELOCITY_GENERATOR_H_

#include <array>
#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <random>
#include <stdexcept>

/**
 * @brief Samples velocities from a user-defined speed-squared distribution g(w).
 *
 * The sampled value is w = |V|^2, so the speed is |V| = sqrt(w) and the
 * direction is isotropic in 3D. (For a distribution originally posed in energy,
 * w = |V|^2 = 2E/m, so the mass drops out of this generator entirely.)
 *
 * @code{.cpp}
 * auto gW = [](double w){ return std::exp(-0.5 * w); };  // any positive g(w), w = |v|^2
 * general_velocity_generator<double> gen;
 * gen.define(gW, 0.0, 40.0);                              // g, v2_min, v2_max
 * gen.seed(12345);
 * auto v = gen();
 * @endcode
 *
 * @see \ref usage for more recipes.
 */
template <typename RealType= double> class general_velocity_generator
{
  public:
    using result_type= std::array<RealType, 3>;
    using point_type= result_type;
    using speed_sq_pdf_function= std::function<RealType(RealType)>;

    general_velocity_generator()
        : m_speedSqPdf(), m_speedSqMin(RealType(0.0)), m_speedSqMax(RealType(1.0)),
          m_maxRejectTries(100000),
          m_pdfUpperBound(RealType(0.0)), m_isDefined(false), m_gen(std::random_device{}()),
          m_distSpeedSq(RealType(0.0), RealType(1.0)), m_distAccept(RealType(0.0), RealType(1.0)),
          m_distCosTheta(RealType(-1.0), RealType(1.0)),
          m_distPhi(RealType(0.0), RealType(2.0) * piValue())
    {
    }

    general_velocity_generator(speed_sq_pdf_function speed_sq_pdf, RealType speed_sq_min,
                               RealType speed_sq_max, int probe_points= 1024,
                               int max_reject_tries= 100000)
        : general_velocity_generator()
    {
        define(speed_sq_pdf, speed_sq_min, speed_sq_max, probe_points, max_reject_tries);
    }

    void define(speed_sq_pdf_function speed_sq_pdf, RealType speed_sq_min, RealType speed_sq_max,
                int probe_points= 1024, int max_reject_tries= 100000)
    {
        m_isDefined= false;

        if (!speed_sq_pdf)
        {
            throw std::invalid_argument("speed_sq_pdf callable must be provided");
        }
        if (speed_sq_min < RealType(0.0))
        {
            throw std::invalid_argument("speed_sq_min must be non-negative (w = |v|^2 >= 0)");
        }
        if (!(speed_sq_min < speed_sq_max))
        {
            throw std::invalid_argument("speed_sq_min must be strictly less than speed_sq_max");
        }
        if (probe_points <= 1)
        {
            throw std::invalid_argument("probe_points must be positive");
        }
        if (max_reject_tries <= 0)
        {
            throw std::invalid_argument("max_reject_tries must be positive");
        }

        m_speedSqPdf= speed_sq_pdf;
        m_speedSqMin= speed_sq_min;
        m_speedSqMax= speed_sq_max;
        m_maxRejectTries= max_reject_tries;
        m_pdfUpperBound= RealType(0.0);

        m_distSpeedSq.param(typename std::uniform_real_distribution<RealType>::param_type(
            m_speedSqMin, m_speedSqMax));
        m_distAccept.param(
            typename std::uniform_real_distribution<RealType>::param_type(RealType(0.0),
                                                                           RealType(1.0)));
        m_distCosTheta.param(typename std::uniform_real_distribution<RealType>::param_type(
            RealType(-1.0), RealType(1.0)));
        m_distPhi.param(typename std::uniform_real_distribution<RealType>::param_type(
            RealType(0.0), RealType(2.0) * piValue()));

        m_pdfUpperBound= estimatePdfUpperBound(probe_points);
        if (m_pdfUpperBound <= RealType(0.0))
        {
            throw std::invalid_argument(
                "speed_sq_pdf must be positive somewhere within [speed_sq_min, speed_sq_max]");
        }

        m_distAccept.param(typename std::uniform_real_distribution<RealType>::param_type(
            RealType(0.0), m_pdfUpperBound));

        reset();
        m_isDefined= true;
    }

    void seed(int s) { applySeed(m_gen, s); }

    void reset()
    {
        m_distSpeedSq.reset();
        m_distAccept.reset();
        m_distCosTheta.reset();
        m_distPhi.reset();
    }

    template <typename GeneratorType> point_type operator()(GeneratorType &gen)
    {
        if (!m_isDefined)
        {
            throw std::runtime_error(
                "general_velocity_generator must be define(...)d before sampling");
        }

        const RealType kOne= RealType(1.0);

        RealType speed= std::sqrt(sampleSpeedSq(gen));   // |V| = sqrt(w), w = |V|^2

        RealType cosTheta= m_distCosTheta(gen);
        RealType sinTheta= std::sqrt(std::max(RealType(0.0), kOne - cosTheta * cosTheta));
        RealType phi= m_distPhi(gen);

        RealType vx= speed * sinTheta * std::cos(phi);
        RealType vy= speed * sinTheta * std::sin(phi);
        RealType vz= speed * cosTheta;

        point_type point;
        point[0]= vx;
        point[1]= vy;
        point[2]= vz;
        return point;
    }

    point_type operator()()
    {
        return (*this)(m_gen);
    }

    RealType speedSqMin() const { return m_speedSqMin; }
    RealType speedSqMax() const { return m_speedSqMax; }
    RealType pdfUpperBound() const { return m_pdfUpperBound; }

  private:
    RealType evaluateSpeedSqPdf(RealType speed_sq) const
    {
        RealType value= m_speedSqPdf(speed_sq);
        if (!std::isfinite(value))
        {
            throw std::runtime_error("speed_sq_pdf returned a non-finite value");
        }
        if (value < RealType(0.0))
        {
            throw std::runtime_error(
                "speed_sq_pdf returned a negative value; expected a non-negative PDF");
        }
        return value;
    }

    RealType estimatePdfUpperBound(int probe_points) const
    {
        RealType max_value= RealType(0.0);
        int i= 0;
        for (i= 0; i < probe_points; ++i)
        {
            RealType ratio= RealType(i) / RealType(probe_points - 1);
            RealType speed_sq= m_speedSqMin + ratio * (m_speedSqMax - m_speedSqMin);
            RealType value= evaluateSpeedSqPdf(speed_sq);
            if (value > max_value)
            {
                max_value= value;
            }
        }

        // Add safety margin in case peak lies between probes.
        return max_value * RealType(1.05);
    }

    template <typename GeneratorType> RealType sampleSpeedSq(GeneratorType &gen)
    {
        int tries= 0;
        for (tries= 0; tries < m_maxRejectTries; ++tries)
        {
            RealType speed_sq= m_distSpeedSq(gen);
            RealType y= m_distAccept(gen);
            RealType pdf_value= evaluateSpeedSqPdf(speed_sq);
            if (y <= pdf_value)
            {
                return speed_sq;
            }
        }

        throw std::runtime_error("failed to sample speed_sq from speed_sq_pdf: rejection sampling "
                                 "exceeded max_reject_tries");
    }

    static void applySeed(std::mt19937 &gen, int s)
    {
        if (s < 0)
            gen.seed(std::random_device{}());
        else
            gen.seed(static_cast<unsigned int>(s));
    }

    static RealType piValue() { return RealType(3.14159265358979323846264338327950288L); }

    speed_sq_pdf_function m_speedSqPdf;
    RealType m_speedSqMin;
    RealType m_speedSqMax;
    int m_maxRejectTries;
    RealType m_pdfUpperBound;
    bool m_isDefined;
    std::mt19937 m_gen;
    std::uniform_real_distribution<RealType> m_distSpeedSq;
    std::uniform_real_distribution<RealType> m_distAccept;
    std::uniform_real_distribution<RealType> m_distCosTheta;
    std::uniform_real_distribution<RealType> m_distPhi;
};

#endif // _GENERAL_VELOCITY_GENERATOR_H_
