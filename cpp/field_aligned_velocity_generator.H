#ifndef _FIELD_ALIGNED_VELOCITY_GENERATOR_H_
#define _FIELD_ALIGNED_VELOCITY_GENERATOR_H_

#include <array>
#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <random>
#include <stdexcept>

/**
 * @brief Samples field-aligned velocities: an arbitrary parallel speed-squared
 *        distribution along B, with a Maxwellian perpendicular spread.
 *
 * The parallel component is w_par = v_par^2, drawn from a user-defined
 * distribution g(w_par); the parallel speed is |v_par| = sqrt(w_par) and its
 * direction is set by a fixed sign (+1 along B, -1 opposite B). (For a
 * distribution originally posed in energy, w_par = v_par^2 = 2E_par/m, so the
 * mass drops out.) The two perpendicular components are independent Gaussians
 * with N(0, theta_perp^2/2), i.e. a Maxwellian with thermal speed theta_perp
 * (same convention as bi_maxwellian_distribution). The result is returned in
 * the global (lab) frame.
 *
 * @code{.cpp}
 * auto gW = [](double w){ return std::exp(-0.5 * w); }; // any positive g(w), w = v_par^2
 * field_aligned_velocity_generator<double> gen;
 * gen.define(gW, 0.0, 40.0, 1.0, {0, 0, 1}, +1);        // g, v2_min, v2_max,
 *                                                        // theta_perp, B, sign
 * gen.seed(12345);
 * auto v = gen();                                        // {vx, vy, vz}
 * @endcode
 *
 * @see \ref usage for more recipes.
 */
template <typename RealType= double> class field_aligned_velocity_generator
{
  public:
    using result_type= std::array<RealType, 3>;
    using point_type= result_type;
    using speed_sq_pdf_function= std::function<RealType(RealType)>;

    field_aligned_velocity_generator()
        : m_speedSqPdf(), m_speedSqMin(RealType(0.0)), m_speedSqMax(RealType(1.0)),
          m_thetaPerp(RealType(1.0)),
          m_ub(point_type({RealType(0.0), RealType(0.0), RealType(1.0)})),
          m_parallelSign(RealType(1.0)), m_maxRejectTries(100000),
          m_pdfUpperBound(RealType(0.0)), m_isDefined(false), m_gen(std::random_device{}()),
          m_distSpeedSq(RealType(0.0), RealType(1.0)), m_distAccept(RealType(0.0), RealType(1.0)),
          m_standardNormal(RealType(0.0), RealType(1.0))
    {
    }

    field_aligned_velocity_generator(speed_sq_pdf_function speed_sq_pdf, RealType speed_sq_min,
                                     RealType speed_sq_max, RealType theta_perp, point_type ub,
                                     int parallel_sign, int probe_points= 1024,
                                     int max_reject_tries= 100000)
        : field_aligned_velocity_generator()
    {
        define(speed_sq_pdf, speed_sq_min, speed_sq_max, theta_perp, ub, parallel_sign,
               probe_points, max_reject_tries);
    }

    void define(speed_sq_pdf_function speed_sq_pdf, RealType speed_sq_min, RealType speed_sq_max,
                RealType theta_perp, point_type ub, int parallel_sign,
                int probe_points= 1024, int max_reject_tries= 100000)
    {
        m_isDefined= false;

        if (!speed_sq_pdf)
        {
            throw std::invalid_argument("speed_sq_pdf callable must be provided");
        }
        if (speed_sq_min < RealType(0.0))
        {
            throw std::invalid_argument("speed_sq_min must be non-negative (w = v_par^2 >= 0)");
        }
        if (!(speed_sq_min < speed_sq_max))
        {
            throw std::invalid_argument("speed_sq_min must be strictly less than speed_sq_max");
        }
        if (theta_perp <= RealType(0.0))
        {
            throw std::invalid_argument("theta_perp must be positive");
        }
        if (parallel_sign != 1 && parallel_sign != -1)
        {
            throw std::invalid_argument("parallel_sign must be exactly +1 or -1");
        }
        if (probe_points <= 1)
        {
            throw std::invalid_argument("probe_points must be positive");
        }
        if (max_reject_tries <= 0)
        {
            throw std::invalid_argument("max_reject_tries must be positive");
        }

        RealType ubNorm2= ub[0] * ub[0] + ub[1] * ub[1] + ub[2] * ub[2];
        if (approximatelyEqual(ubNorm2, RealType(0.0)))
        {
            throw std::invalid_argument("ub (magnetic field) must be a non-zero vector");
        }

        m_speedSqPdf= speed_sq_pdf;
        m_speedSqMin= speed_sq_min;
        m_speedSqMax= speed_sq_max;
        m_thetaPerp= theta_perp;
        m_ub= ub;
        m_parallelSign= static_cast<RealType>(parallel_sign);
        m_maxRejectTries= max_reject_tries;
        m_pdfUpperBound= RealType(0.0);

        m_distSpeedSq.param(typename std::uniform_real_distribution<RealType>::param_type(
            m_speedSqMin, m_speedSqMax));
        m_distAccept.param(
            typename std::uniform_real_distribution<RealType>::param_type(RealType(0.0),
                                                                           RealType(1.0)));

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
        m_standardNormal.reset();
    }

    template <typename GeneratorType> point_type operator()(GeneratorType &gen)
    {
        if (!m_isDefined)
        {
            throw std::runtime_error(
                "field_aligned_velocity_generator must be define(...)d before sampling");
        }

        const RealType sigmaScale= RealType(1.0) / std::sqrt(RealType(2.0));

        // Parallel: v_par^2 -> signed speed along B.
        RealType vPar= m_parallelSign * std::sqrt(sampleSpeedSq(gen));

        // Perpendicular: 2D Maxwellian, N(0, theta_perp^2/2) per component.
        RealType vPerp1= m_thetaPerp * sigmaScale * m_standardNormal(gen);
        RealType vPerp2= m_thetaPerp * sigmaScale * m_standardNormal(gen);

        // Field-aligned frame: (perp1 -> e1, perp2 -> e2, par -> e3 = B/|B|).
        point_type local;
        local[0]= vPerp1;
        local[1]= vPerp2;
        local[2]= vPar;
        return rotate_from_fieldAligned_frame(local, m_ub);
    }

    point_type operator()()
    {
        return (*this)(m_gen);
    }

    RealType parallelSpeedSqMin() const { return m_speedSqMin; }
    RealType parallelSpeedSqMax() const { return m_speedSqMax; }
    RealType thetaPerp() const { return m_thetaPerp; }
    point_type ub() const { return m_ub; }
    RealType parallelSign() const { return m_parallelSign; }
    RealType pdfUpperBound() const { return m_pdfUpperBound; }

    // ====================================================================
    // Geometry helper
    // ====================================================================
    // Transform local velocity components (perp1, perp2, par) to global
    // Cartesian coordinates given the parallel direction ub.
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

    static bool approximatelyEqual(RealType a, RealType b)
    {
        RealType scale= std::max({RealType(1.0), std::abs(a), std::abs(b)});
        return std::abs(a - b) <= std::numeric_limits<RealType>::epsilon() * scale;
    }

    speed_sq_pdf_function m_speedSqPdf;
    RealType m_speedSqMin;
    RealType m_speedSqMax;
    RealType m_thetaPerp;
    point_type m_ub;
    RealType m_parallelSign;
    int m_maxRejectTries;
    RealType m_pdfUpperBound;
    bool m_isDefined;
    std::mt19937 m_gen;
    std::uniform_real_distribution<RealType> m_distSpeedSq;
    std::uniform_real_distribution<RealType> m_distAccept;
    std::normal_distribution<RealType> m_standardNormal;
};

#endif // _FIELD_ALIGNED_VELOCITY_GENERATOR_H_
