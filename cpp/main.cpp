#include <random>
#include <cmath>
#include <iostream>
#include <fstream>

typedef double Real;

// Local headers
#include "bi_kappa_distribution.H"
#include "bi_maxwellian_distribution.H"
#include "general_position_generator.H"
#include "general_velocity_generator.H"
#include "test_suite.H"

int main()
{
    const Real kappa = 2.0;
    const Real theta_perp = 1.0;
    const Real theta_par = 2.0;
    const int n_particle = 100000;

    const unsigned int seed = 20030410u;
    std::mt19937 gen(seed);

    if (run_all_tests() != 0)
    {
        return 1;
    }
    std::cout << std::endl;

    std::cout << "=== Example 1: BiKappaDistribution ===" << std::endl;
    bi_kappa_distribution<Real> biKappa;
    biKappa.define(kappa, theta_perp, theta_par, {0,0,1}, 100.0);
    {
        std::ofstream out("samples_bikappa.txt");
        for (int i = 0; i < n_particle; ++i) {
            const bi_kappa_distribution<Real>::point_type v = biKappa(gen);
            out << v[0] << " " << v[1] << " " << v[2] << "\n";
        }
        std::cout << "wrote " << n_particle << " samples to samples_bikappa.txt" << std::endl;
    }

    std::cout << "\n=== Example 2: BiMaxwellianDistribution ===" << std::endl;
    bi_maxwellian_distribution<Real> biMaxwell;
    biMaxwell.define(theta_perp, theta_par);
    {
        std::ofstream out("samples_bimaxwellian.txt");
        for (int i = 0; i < n_particle; ++i) {
            const bi_maxwellian_distribution<Real>::point_type v = biMaxwell(gen);
            out << v[0] << " " << v[1] << " " << v[2] << "\n";
        }
        std::cout << "wrote " << n_particle << " samples to samples_bimaxwellian.txt" << std::endl;
    }

    std::cout << "\n=== Example 3: general_velocity_generator (f(E)=exp(-E)) ===" << std::endl;
    auto energyPdf = [](Real E) -> Real { return std::exp(-E); };
    general_velocity_generator<Real> velocityGen;
    velocityGen.define(energyPdf,
                       static_cast<Real>(0.0),
                       static_cast<Real>(20.0),
                       static_cast<Real>(1.0));
    {
        std::ofstream out("samples_general_velocity.txt");
        for (int i = 0; i < n_particle; ++i) {
            const general_velocity_generator<Real>::point_type v = velocityGen(gen);
            out << v[0] << " " << v[1] << " " << v[2] << "\n";
        }
        std::cout << "wrote " << n_particle << " samples to samples_general_velocity.txt"
                  << std::endl;
    }

    std::cout << "\n=== Example 4: general_position_generator (rho=1+sin(x)sin(y)) ===" << std::endl;
    typedef general_position_generator<Real>::point_type position_point_type;
    auto rho = [](const position_point_type &x) -> Real {
        return 1 + std::sin(x[0]) * std::sin(x[1]);
    };
    static double pi = 3.14159;
    general_position_generator<Real> positionGen;
    positionGen.define(2, {-pi, -pi}, {pi, pi}, rho);
    {
        std::ofstream out("samples_general_position.txt");
        for (int i = 0; i < n_particle; ++i) {
            const position_point_type x = positionGen(gen);
            out << x[0] << " " << x[1] << "\n";
        }
        std::cout << "wrote " << n_particle << " samples to samples_general_position.txt"
                  << std::endl;
    }

    return 0;
}
