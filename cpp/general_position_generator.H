#ifndef _GENERAL_POSITION_GENERATOR_H_
#define _GENERAL_POSITION_GENERATOR_H_

#include <cmath>
#include <cstddef>
#include <functional>
#include <random>
#include <stdexcept>
#include <vector>

/**
 * @brief Samples positions in 1D/2D/3D from a user density rho(x), rho(x,y), or rho(x,y,z).
 *
 * The generator uses rejection sampling in an axis-aligned box domain.
 */
template <typename RealType= double> class general_position_generator
{
  public:
    using result_type= std::vector<RealType>;
    using point_type= result_type;
    using density_function= std::function<RealType(const point_type &)>;

    general_position_generator(int dimension, const point_type &lower_bounds,
                               const point_type &upper_bounds, density_function density_function,
                               int probe_points= 4096, int max_reject_tries= 200000)
        : m_dimension(dimension), m_lowerBounds(lower_bounds), m_upperBounds(upper_bounds),
          m_rangeWidths(static_cast<std::size_t>(dimension), RealType(0.0)),
          m_densityFunction(density_function), m_maxRejectTries(max_reject_tries),
          m_densityUpperBound(RealType(0.0)), m_distAccept(RealType(0.0), RealType(1.0)),
          m_distUnit(RealType(0.0), RealType(1.0))
    {
        if (!(m_dimension >= 1 && m_dimension <= 3))
        {
            throw std::invalid_argument("dimension must be 1, 2, or 3");
        }
        if (!m_densityFunction)
        {
            throw std::invalid_argument("density_function callable must be provided");
        }
        if (m_lowerBounds.size() != static_cast<std::size_t>(m_dimension) ||
            m_upperBounds.size() != static_cast<std::size_t>(m_dimension))
        {
            throw std::invalid_argument("lower_bounds and upper_bounds sizes must match dimension");
        }
        if (probe_points <= 0)
        {
            throw std::invalid_argument("probe_points must be positive");
        }
        if (m_maxRejectTries <= 0)
        {
            throw std::invalid_argument("max_reject_tries must be positive");
        }

        int d= 0;
        for (d= 0; d < m_dimension; ++d)
        {
            std::size_t idx= static_cast<std::size_t>(d);
            if (!(m_lowerBounds[idx] < m_upperBounds[idx]))
            {
                throw std::invalid_argument("each lower bound must be strictly less than its upper "
                                            "bound");
            }
            m_rangeWidths[idx]= m_upperBounds[idx] - m_lowerBounds[idx];
        }

        m_densityUpperBound= estimateDensityUpperBound(probe_points);
        if (m_densityUpperBound <= RealType(0.0))
        {
            throw std::invalid_argument("density_function must be positive somewhere inside the "
                                        "domain");
        }

        m_distAccept.param(typename std::uniform_real_distribution<RealType>::param_type(
            RealType(0.0), m_densityUpperBound));
    }

    void reset()
    {
        m_distAccept.reset();
        m_distUnit.reset();
    }

    int dimension() const { return m_dimension; }
    const point_type &lowerBounds() const { return m_lowerBounds; }
    const point_type &upperBounds() const { return m_upperBounds; }
    RealType densityUpperBound() const { return m_densityUpperBound; }

    template <typename GeneratorType> point_type operator()(GeneratorType &gen) const
    {
        int tries= 0;
        for (tries= 0; tries < m_maxRejectTries; ++tries)
        {
            point_type point= sampleUniformPoint(gen);
            RealType rho= evaluateDensity(point);
            RealType y= m_distAccept(gen);
            if (y <= rho)
            {
                return point;
            }
        }

        throw std::runtime_error("failed to sample position from density_function: rejection "
                                 "sampling exceeded max_reject_tries");
    }

  private:
    RealType evaluateDensity(const point_type &point) const
    {
        RealType value= m_densityFunction(point);
        if (!std::isfinite(value))
        {
            throw std::runtime_error("density_function returned a non-finite value");
        }
        if (value < RealType(0.0))
        {
            throw std::runtime_error("density_function returned a negative value; expected a "
                                     "non-negative density");
        }
        return value;
    }

    template <typename GeneratorType> point_type sampleUniformPoint(GeneratorType &gen) const
    {
        point_type point(static_cast<std::size_t>(m_dimension), RealType(0.0));
        int d= 0;
        for (d= 0; d < m_dimension; ++d)
        {
            std::size_t idx= static_cast<std::size_t>(d);
            point[idx]= m_lowerBounds[idx] + m_rangeWidths[idx] * m_distUnit(gen);
        }
        return point;
    }

    RealType estimateDensityUpperBound(int probe_points) const
    {
        std::mt19937 probeGen(1337u);
        RealType max_value= RealType(0.0);

        int i= 0;
        for (i= 0; i < probe_points; ++i)
        {
            point_type point= sampleUniformPoint(probeGen);
            RealType value= evaluateDensity(point);
            if (value > max_value)
            {
                max_value= value;
            }
        }

        // Add safety margin in case peak lies between probes.
        return max_value * RealType(1.05);
    }

    int m_dimension;
    point_type m_lowerBounds;
    point_type m_upperBounds;
    point_type m_rangeWidths;
    density_function m_densityFunction;
    int m_maxRejectTries;
    RealType m_densityUpperBound;
    mutable std::uniform_real_distribution<RealType> m_distAccept;
    mutable std::uniform_real_distribution<RealType> m_distUnit;
};

#endif // _GENERAL_POSITION_GENERATOR_H_
