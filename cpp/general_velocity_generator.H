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
class general_velocity_generator
{
public:
        using result_type = std::array<RealType, 3>;
        using point_type = result_type;
        using energy_pdf_function = std::function<RealType(RealType)>;

        general_velocity_generator(energy_pdf_function energy_pdf,
                                                             RealType energy_min,
                                                             RealType energy_max,
                                                             RealType particle_mass,
                                                             int probe_points = 1024,
                                                             int max_reject_tries = 100000)
        : m_energyPdf(energy_pdf),
          m_energyMin(energy_min),
          m_energyMax(energy_max),
          m_particleMass(particle_mass),
            m_maxRejectTries(max_reject_tries),
            m_pdfUpperBound(RealType(0.0)),
            m_distEnergy(energy_min, energy_max),
            m_distAccept(RealType(0.0), RealType(1.0)),
            m_distCosTheta(RealType(-1.0), RealType(1.0)),
            m_distPhi(RealType(0.0), RealType(2.0) * piValue())
    {
        if (!m_energyPdf) {
            throw std::invalid_argument("energy_pdf callable must be provided");
        }
        if (!(m_energyMin < m_energyMax)) {
            throw std::invalid_argument("energy_min must be strictly less than energy_max");
        }
        if (m_particleMass <= RealType(0.0)) {
            throw std::invalid_argument("particle_mass must be positive");
        }
        if (probe_points <= 1) {
            throw std::invalid_argument("probe_points must be positive");
        }
        if (m_maxRejectTries <= 0) {
            throw std::invalid_argument("max_reject_tries must be positive");
        }

        m_pdfUpperBound= estimatePdfUpperBound(probe_points);
        if (m_pdfUpperBound <= RealType(0.0)) {
            throw std::invalid_argument("energy_pdf must be positive somewhere within [energy_min, energy_max]");
        }

        m_distAccept.param(
            typename std::uniform_real_distribution<RealType>::param_type(RealType(0.0), m_pdfUpperBound));
    }

    void reset()
    {
        m_distEnergy.reset();
        m_distAccept.reset();
        m_distCosTheta.reset();
        m_distPhi.reset();
    }

    template<typename GeneratorType>
    point_type operator()(GeneratorType& gen) const
    {
        const RealType kOne= RealType(1.0);
        const RealType kTwo= RealType(2.0);

        RealType energy= sampleEnergy(gen);
        RealType speed= std::sqrt(kTwo * energy / m_particleMass);

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
        if (value < RealType(0.0)) {
            throw std::runtime_error("energy_pdf returned a negative value; expected a non-negative PDF");
        }
        return value;
    }

    RealType estimatePdfUpperBound(int probe_points) const
    {
        RealType max_value= RealType(0.0);
        int i= 0;
        for (i = 0; i < probe_points; ++i) {
            RealType ratio= RealType(i) / RealType(probe_points - 1);
            RealType energy= m_energyMin + ratio * (m_energyMax - m_energyMin);
            RealType value= evaluateEnergyPdf(energy);
            if (value > max_value) {
                max_value= value;
            }
        }

        // Add safety margin in case peak lies between probes.
        return max_value * RealType(1.05);
    }

    template<typename GeneratorType>
    RealType sampleEnergy(GeneratorType& gen) const
    {
        int tries= 0;
        for (tries = 0; tries < m_maxRejectTries; ++tries) {
            RealType energy= m_distEnergy(gen);
            RealType y= m_distAccept(gen);
            RealType pdf_value= evaluateEnergyPdf(energy);
            if (y <= pdf_value) {
                return energy;
            }
        }

        throw std::runtime_error("failed to sample energy from energy_pdf: rejection sampling exceeded max_reject_tries");
    }

    static RealType piValue() { return RealType(3.14159265358979323846264338327950288L); }

    energy_pdf_function m_energyPdf;
    RealType m_energyMin;
    RealType m_energyMax;
    RealType m_particleMass;
    int m_maxRejectTries;
    RealType m_pdfUpperBound;
    mutable std::uniform_real_distribution<RealType> m_distEnergy;
    mutable std::uniform_real_distribution<RealType> m_distAccept;
    mutable std::uniform_real_distribution<RealType> m_distCosTheta;
    mutable std::uniform_real_distribution<RealType> m_distPhi;
};

#endif // _GENERAL_VELOCITY_GENERATOR_H_
