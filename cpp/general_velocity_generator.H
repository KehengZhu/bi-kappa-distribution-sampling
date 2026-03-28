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
 * @brief Samples velocities from a user-defined energy distribution f(E).
 *
 * The sampled velocity satisfies E = m * |V|^2 / 2 and direction is isotropic in 3D.
 */
template<typename RealType = double>
class GeneralVelocityGenerator
{
public:
    typedef std::array<RealType, 3> Point3D; // modern form: using Point3D = std::array<RealType, 3>;
    typedef std::function<RealType(RealType)> EnergyPdfFunction; // modern form: using EnergyPdfFunction = std::function<RealType(RealType)>;

    GeneralVelocityGenerator(EnergyPdfFunction energy_pdf,
                             RealType energy_min,
                             RealType energy_max,
                             RealType particle_mass,
                             int probe_points = 1024,
                             int max_reject_tries = 100000)
        : m_energyPdf(energy_pdf),
          m_energyMin(energy_min),
          m_energyMax(energy_max),
          m_particleMass(particle_mass),
          m_maxRejectTries(max_reject_tries)
    {
        if (!m_energyPdf) {
            throw std::invalid_argument("energy_pdf callable must be provided");
        }
        if (!(m_energyMin < m_energyMax)) {
            throw std::invalid_argument("energy_min must be strictly less than energy_max");
        }
        if (m_particleMass <= static_cast<RealType>(0.0)) {
            throw std::invalid_argument("particle_mass must be positive");
        }
        if (probe_points <= 1) {
            throw std::invalid_argument("probe_points must be positive");
        }
        if (m_maxRejectTries <= 0) {
            throw std::invalid_argument("max_reject_tries must be positive");
        }

        m_pdfUpperBound= estimatePdfUpperBound(probe_points);
        if (m_pdfUpperBound <= static_cast<RealType>(0.0)) {
            throw std::invalid_argument("energy_pdf must be positive somewhere within [energy_min, energy_max]");
        }
    }

    template<typename GeneratorType>
    Point3D operator()(GeneratorType& gen) const
    {
        RealType energy= sampleEnergy(gen);
        RealType speed= std::sqrt(static_cast<RealType>(2.0) * energy / m_particleMass);

        RealType pi= static_cast<RealType>(3.14159265358979323846264338327950288L);
        std::uniform_real_distribution<RealType> distCosTheta(static_cast<RealType>(-1.0), static_cast<RealType>(1.0));
        std::uniform_real_distribution<RealType> distPhi(static_cast<RealType>(0.0), static_cast<RealType>(2.0) * pi);

        RealType cosTheta= distCosTheta(gen);
        RealType sinTheta= std::sqrt(std::max(static_cast<RealType>(0.0), static_cast<RealType>(1.0) - cosTheta * cosTheta));
        RealType phi= distPhi(gen);

        RealType vx= speed * sinTheta * std::cos(phi);
        RealType vy= speed * sinTheta * std::sin(phi);
        RealType vz= speed * cosTheta;

        Point3D point;
        point[0]= vx;
        point[1]= vy;
        point[2]= vz;
        // modern aggregate alternative: return {vx, vy, vz};
        return point;
    }

    RealType energyMin() const { return m_energyMin; }
    RealType energyMax() const { return m_energyMax; }
    RealType particleMass() const { return m_particleMass; }
    RealType pdfUpperBound() const { return m_pdfUpperBound; }

private:
    RealType evaluateEnergyPdf(RealType energy) const
    {
        RealType value= m_energyPdf(energy);
        if (!std::isfinite(value)) {
            throw std::runtime_error("energy_pdf returned a non-finite value");
        }
        if (value < static_cast<RealType>(0.0)) {
            throw std::runtime_error("energy_pdf returned a negative value; expected a non-negative PDF");
        }
        return value;
    }

    RealType estimatePdfUpperBound(int probe_points) const
    {
        RealType max_value= static_cast<RealType>(0.0);
        int i= 0;
        for (i = 0; i < probe_points; ++i) {
            RealType ratio= static_cast<RealType>(i) / static_cast<RealType>(probe_points - 1);
            RealType energy= m_energyMin + ratio * (m_energyMax - m_energyMin);
            RealType value= evaluateEnergyPdf(energy);
            if (value > max_value) {
                max_value= value;
            }
        }

        // Add safety margin in case peak lies between probes.
        return max_value * static_cast<RealType>(1.05);
    }

    template<typename GeneratorType>
    RealType sampleEnergy(GeneratorType& gen) const
    {
        std::uniform_real_distribution<RealType> distEnergy(m_energyMin, m_energyMax);
        std::uniform_real_distribution<RealType> distAccept(static_cast<RealType>(0.0), m_pdfUpperBound);

        int tries= 0;
        for (tries = 0; tries < m_maxRejectTries; ++tries) {
            RealType energy= distEnergy(gen);
            RealType y= distAccept(gen);
            RealType pdf_value= evaluateEnergyPdf(energy);
            if (y <= pdf_value) {
                return energy;
            }
        }

        throw std::runtime_error("failed to sample energy from energy_pdf: rejection sampling exceeded max_reject_tries");
    }

    EnergyPdfFunction m_energyPdf;
    RealType m_energyMin;
    RealType m_energyMax;
    RealType m_particleMass;
    int m_maxRejectTries;
    RealType m_pdfUpperBound;
};

#endif // _GENERAL_VELOCITY_GENERATOR_H_
