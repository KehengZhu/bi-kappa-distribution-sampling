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
#include "field_aligned_velocity_generator.H"
#include "test_suite.H"

int main()
{
    const Real kappa = 2.0;
    const Real theta_perp = 1.0;
    const Real theta_par = 2.0;
    const int n_particle = 200000;

    if (run_all_tests() != 0)
    {
        return 1;
    }
    std::cout << std::endl;

    std::cout << "=== Example 1: BiKappaDistribution ===" << std::endl;
    bi_kappa_distribution<Real> biKappa;
    biKappa.define(kappa, theta_perp, theta_par, {0,0,1}, 100.0, 20030410);
    {
        std::ofstream out("samples_bikappa.txt");
        for (int i = 0; i < n_particle; ++i) {
            const bi_kappa_distribution<Real>::point_type v = biKappa();
            out << v[0] << " " << v[1] << " " << v[2] << "\n";
        }
        std::cout << "wrote " << n_particle << " samples to samples_bikappa.txt" << std::endl;
    }

    std::cout << "\n=== Example 2: BiMaxwellianDistribution ===" << std::endl;
    bi_maxwellian_distribution<Real> biMaxwell;
    biMaxwell.define(theta_perp, theta_par, {0,0,1}, 20.0, 20030410);
    {
        std::ofstream out("samples_bimaxwellian.txt");
        for (int i = 0; i < n_particle; ++i) {
            const bi_maxwellian_distribution<Real>::point_type v = biMaxwell();
            out << v[0] << " " << v[1] << " " << v[2] << "\n";
        }
        std::cout << "wrote " << n_particle << " samples to samples_bimaxwellian.txt" << std::endl;
    }

    std::cout << "\n=== Example 3: general_velocity_generator (isotropic Maxwellian via w=|v|^2) ==="
              << std::endl;
    // g(w) = sqrt(w) * exp(-w/theta^2) with w = |v|^2 is the speed-squared density
    // of an isotropic Maxwellian of thermal speed theta: with dw = 2|v| d|v| it maps
    // to the Maxwell speed pdf ~ |v|^2 exp(-|v|^2/theta^2), i.e. each Cartesian
    // component is N(0, theta^2/2). theta = theta_perp = 1 here, and w_max = 25
    // truncates at |v| = 5 theta.
    const Real theta_iso = theta_perp;
    auto speedSqPdf = [theta_iso](Real w) -> Real {
        return std::sqrt(w) * std::exp(-w / (theta_iso * theta_iso));
    };
    general_velocity_generator<Real> velocityGen;
    velocityGen.define(speedSqPdf,
                       static_cast<Real>(0.0),
                       static_cast<Real>(25.0));
    velocityGen.seed(20030410);
    {
        std::ofstream out("samples_general_velocity.txt");
        for (int i = 0; i < n_particle; ++i) {
            const general_velocity_generator<Real>::point_type v = velocityGen();
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
    positionGen.seed(20030410);
    {
        std::ofstream out("samples_general_position.txt");
        for (int i = 0; i < n_particle; ++i) {
            const position_point_type x = positionGen();
            out << x[0] << " " << x[1] << "\n";
        }
        std::cout << "wrote " << n_particle << " samples to samples_general_position.txt"
                  << std::endl;
    }

    std::cout << "\n=== Example 5: field_aligned_velocity_generator (beam along B=(1,1,1)) ==="
              << std::endl;
    // g(w_par) = exp(-w_par/theta_par^2) with w_par = v_par^2 gives the Rayleigh
    // parallel beam p(v_par) = (2 v_par/theta_par^2) exp(-v_par^2/theta_par^2),
    // v_par >= 0, peaking at v_par = theta_par/sqrt(2). The perpendicular plane is
    // Maxwellian with thermal speed theta_perp. B is tilted off the axes so the
    // field-aligned rotation is exercised; w_max = 100 truncates at v_par = 5 theta_par.
    auto parallelSpeedSqPdf = [theta_par](Real w) -> Real {
        return std::exp(-w / (theta_par * theta_par));
    };
    field_aligned_velocity_generator<Real> fieldAlignedGen;
    fieldAlignedGen.define(parallelSpeedSqPdf,
                           static_cast<Real>(0.0),
                           static_cast<Real>(100.0),
                           theta_perp,
                           {1, 1, 1},
                           +1);
    fieldAlignedGen.seed(20030410);
    {
        std::ofstream out("samples_field_aligned.txt");
        for (int i = 0; i < n_particle; ++i) {
            const field_aligned_velocity_generator<Real>::point_type v = fieldAlignedGen();
            out << v[0] << " " << v[1] << " " << v[2] << "\n";
        }
        std::cout << "wrote " << n_particle << " samples to samples_field_aligned.txt"
                  << std::endl;
    }

    return 0;
}
