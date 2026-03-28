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
template <typename RealType= double> class GeneralPositionGenerator
{
  public:
    typedef std::vector<RealType> PositionVector; // modern form: using PositionVector =
                                                  // std::vector<RealType>;
    typedef std::function<RealType(const PositionVector &)> DensityFunction; // modern form: using
                                                                             // DensityFunction =
                                                                             // std::function<RealType(const
                                                                             // PositionVector&)>;

    GeneralPositionGenerator(int dimension, const PositionVector &lower_bounds,
                             const PositionVector &upper_bounds, DensityFunction density_function,
                             int probe_points= 4096, int max_reject_tries= 200000)
        : m_dimension(dimension), m_lowerBounds(lower_bounds), m_upperBounds(upper_bounds),
          m_densityFunction(density_function), m_maxRejectTries(max_reject_tries)
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
            if (!(m_lowerBounds[static_cast<std::size_t>(d)] <
                  m_upperBounds[static_cast<std::size_t>(d)]))
            {
                throw std::invalid_argument("each lower bound must be strictly less than its upper "
                                            "bound");
            }
        }

        m_densityUpperBound= estimateDensityUpperBound(probe_points);
        if (m_densityUpperBound <= static_cast<RealType>(0.0))
        {
            throw std::invalid_argument("density_function must be positive somewhere inside the "
                                        "domain");
        }
    }

    int dimension() const
    {
        return m_dimension;
    }
    const PositionVector &lowerBounds() const
    {
        return m_lowerBounds;
    }
    const PositionVector &upperBounds() const
    {
        return m_upperBounds;
    }
    RealType densityUpperBound() const
    {
        return m_densityUpperBound;
    }

    template <typename GeneratorType> PositionVector operator()(GeneratorType &gen) const
    {
        std::uniform_real_distribution<RealType> distAccept(static_cast<RealType>(0.0),
                                                            m_densityUpperBound);

        int tries= 0;
        for (tries= 0; tries < m_maxRejectTries; ++tries)
        {
            PositionVector point= sampleUniformPoint(gen);
            RealType rho= evaluateDensity(point);
            RealType y= distAccept(gen);
            if (y <= rho)
            {
                return point;
            }
        }

        throw std::runtime_error("failed to sample position from density_function: rejection "
                                 "sampling exceeded max_reject_tries");
    }

  private:
    RealType evaluateDensity(const PositionVector &point) const
    {
        RealType value= m_densityFunction(point);
        if (!std::isfinite(value))
        {
            throw std::runtime_error("density_function returned a non-finite value");
        }
        if (value < static_cast<RealType>(0.0))
        {
            throw std::runtime_error("density_function returned a negative value; expected a "
                                     "non-negative density");
        }
        return value;
    }

    template <typename GeneratorType> PositionVector sampleUniformPoint(GeneratorType &gen) const
    {
        PositionVector point(static_cast<std::size_t>(m_dimension), static_cast<RealType>(0.0));
        int d= 0;
        for (d= 0; d < m_dimension; ++d)
        {
            std::uniform_real_distribution<RealType> distAxis(
                m_lowerBounds[static_cast<std::size_t>(d)],
                m_upperBounds[static_cast<std::size_t>(d)]);
            point[static_cast<std::size_t>(d)]= distAxis(gen);
        }
        return point;
    }

    RealType estimateDensityUpperBound(int probe_points) const
    {
        std::mt19937 probeGen(1337u);
        RealType max_value= static_cast<RealType>(0.0);

        int i= 0;
        for (i= 0; i < probe_points; ++i)
        {
            PositionVector point= sampleUniformPoint(probeGen);
            RealType value= evaluateDensity(point);
            if (value > max_value)
            {
                max_value= value;
            }
        }

        // Add safety margin in case peak lies between probes.
        return max_value * static_cast<RealType>(1.05);
    }

    int m_dimension;
    PositionVector m_lowerBounds;
    PositionVector m_upperBounds;
    DensityFunction m_densityFunction;
    int m_maxRejectTries;
    RealType m_densityUpperBound;
};

#endif // _GENERAL_POSITION_GENERATOR_H_
