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

int test_transform()
{
    auto abs_real = [](Real x) -> Real { return (x < Real(0.0)) ? -x : x; };
    auto approx_equal = [&](Real a, Real b, Real tol) -> bool {
        return abs_real(a - b) <= tol;
    };
    auto dot3 = [](const std::array<Real, 3> &a, const std::array<Real, 3> &b) -> Real {
        return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    };
    auto norm3 = [](const std::array<Real, 3> &v) -> Real {
        return std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    };
    auto normalize3 = [&](const std::array<Real, 3> &v) -> std::array<Real, 3> {
        Real n = norm3(v);
        return {v[0] / n, v[1] / n, v[2] / n};
    };
    auto add3 = [](const std::array<Real, 3> &a,
                   const std::array<Real, 3> &b) -> std::array<Real, 3> {
        return {a[0] + b[0], a[1] + b[1], a[2] + b[2]};
    };

    {
        std::cout << "=== Transform self-test (BiKappa/BiMaxwellian) ===" << std::endl;

        const Real tol = Real(1e-10);
        bool all_ok = true;

        bi_kappa_distribution<Real> kappaTest(Real(2.0), Real(1.0), Real(2.0));
        bi_maxwellian_distribution<Real> maxwellTest(Real(1.0), Real(2.0));

        // Test 1: For any ub, parallel projection equals local[2],
        // and perpendicular magnitude equals sqrt(local[0]^2 + local[1]^2).
        {
            bi_kappa_distribution<Real>::point_type local = {Real(0.7), Real(-1.2), Real(2.5)};
            bi_kappa_distribution<Real>::point_type ub = {Real(1.0), Real(2.0), Real(2.0)};
            bi_kappa_distribution<Real>::point_type out =
                kappaTest.rotate_from_fieldAligned_frame(local, ub);
            bi_kappa_distribution<Real>::point_type ubHat = normalize3(ub);

            Real parallel = dot3(out, ubHat);
            Real outNorm = norm3(out);
            Real perpNormSq = outNorm * outNorm - parallel * parallel;
            Real localPerpNormSq = local[0] * local[0] + local[1] * local[1];

            bool ok = approx_equal(parallel, local[2], Real(1e-9)) &&
                      approx_equal(perpNormSq, localPerpNormSq, Real(1e-8));
            all_ok = all_ok && ok;
            std::cout << "BiKappa component decomposition: " << (ok ? "PASS" : "FAIL")
                   << std::endl;
            if (!ok)
            {
                std::cout << "  local: (" << local[0] << ", " << local[1] << ", " << local[2] << ")\n";
                std::cout << "  ub: (" << ub[0] << ", " << ub[1] << ", " << ub[2] << ")\n";
                std::cout << "  out: (" << out[0] << ", " << out[1] << ", " << out[2] << ")\n";
                std::cout << "  parallel(out, ubHat): " << parallel << " expected " << local[2]
                       << "\n";
                std::cout << "  perp^2(out): " << perpNormSq << " expected " << localPerpNormSq
                       << "\n";
            }
        }

        // Test 2: parallel component maps along ub direction (non-normalized ub)
        {
            bi_maxwellian_distribution<Real>::point_type local = {Real(0.0), Real(0.0), Real(3.0)};
            bi_maxwellian_distribution<Real>::point_type ub = {Real(0.0), Real(5.0), Real(0.0)};
            bi_maxwellian_distribution<Real>::point_type out =
                maxwellTest.rotate_from_fieldAligned_frame(local, ub);

            bool ok = approx_equal(out[0], Real(0.0), tol) &&
                      approx_equal(out[1], Real(3.0), tol) &&
                      approx_equal(out[2], Real(0.0), tol);
            all_ok = all_ok && ok;
            std::cout << "BiMaxwellian parallel-axis mapping: " << (ok ? "PASS" : "FAIL")
                   << std::endl;
        }

        // Test 3: Dot products are preserved under a rotation.
        {
            bi_kappa_distribution<Real>::point_type a = {Real(1.2), Real(-0.4), Real(2.7)};
            bi_kappa_distribution<Real>::point_type b = {Real(-0.8), Real(1.5), Real(0.3)};
            bi_kappa_distribution<Real>::point_type ub = {Real(2.0), Real(-1.0), Real(3.0)};
            bi_kappa_distribution<Real>::point_type ra =
                kappaTest.rotate_from_fieldAligned_frame(a, ub);
            bi_kappa_distribution<Real>::point_type rb =
                kappaTest.rotate_from_fieldAligned_frame(b, ub);

            bool ok = approx_equal(dot3(a, b), dot3(ra, rb), Real(1e-9));
            all_ok = all_ok && ok;
            std::cout << "BiKappa dot-product preservation: " << (ok ? "PASS" : "FAIL")
                   << std::endl;
        }

        // Test 4: Linearity T(a+b)=T(a)+T(b) for the same ub.
        {
            bi_maxwellian_distribution<Real>::point_type a = {Real(0.4), Real(1.1), Real(-0.2)};
            bi_maxwellian_distribution<Real>::point_type b = {Real(-1.3), Real(0.7), Real(2.2)};
            bi_maxwellian_distribution<Real>::point_type ab = add3(a, b);
            bi_maxwellian_distribution<Real>::point_type ub = {Real(1.0), Real(3.0), Real(2.0)};

            bi_maxwellian_distribution<Real>::point_type ta =
                maxwellTest.rotate_from_fieldAligned_frame(a, ub);
            bi_maxwellian_distribution<Real>::point_type tb =
                maxwellTest.rotate_from_fieldAligned_frame(b, ub);
            bi_maxwellian_distribution<Real>::point_type tab =
                maxwellTest.rotate_from_fieldAligned_frame(ab, ub);
            bi_maxwellian_distribution<Real>::point_type tSum = add3(ta, tb);

            bool ok = approx_equal(tab[0], tSum[0], Real(1e-9)) &&
                      approx_equal(tab[1], tSum[1], Real(1e-9)) &&
                      approx_equal(tab[2], tSum[2], Real(1e-9));
            all_ok = all_ok && ok;
            std::cout << "BiMaxwellian linearity: " << (ok ? "PASS" : "FAIL") << std::endl;
        }

        if (!all_ok)
        {
            std::cout << "Transform self-test failed. Exiting." << std::endl;
            return 1;
        }
        std::cout << "Transform self-test: ALL PASS" << std::endl;
    }
    return 0;
}

int main()
{
    const Real kappa = 2.0;
    const Real theta_perp = 1.0;
    const Real theta_par = 2.0;
    const int n_particle = 100000;

    std::random_device rd;
    std::mt19937 gen(rd());

    std::cout << "=== Test: Rotation transform ===" << std::endl;
    if (test_transform() != 0)
    {
        return 1;
    }
    std::cout << std::endl;

    std::cout << "=== Example 1: BiKappaDistribution ===" << std::endl;
    bi_kappa_distribution<Real> biKappa(kappa, theta_perp, theta_par);
    {
        std::ofstream out("samples_bikappa.txt");
        for (int i = 0; i < n_particle; ++i) {
            const bi_kappa_distribution<Real>::point_type v = biKappa(gen);
            out << v[0] << " " << v[1] << " " << v[2] << "\n";
        }
        std::cout << "wrote " << n_particle << " samples to samples_bikappa.txt" << std::endl;
    }

    std::cout << "\n=== Example 2: BiMaxwellianDistribution ===" << std::endl;
    bi_maxwellian_distribution<Real> biMaxwell(theta_perp, theta_par);
    {
        std::ofstream out("samples_bimaxwellian.txt");
        for (int i = 0; i < n_particle; ++i) {
            const bi_maxwellian_distribution<Real>::point_type v = biMaxwell(gen);
            out << v[0] << " " << v[1] << " " << v[2] << "\n";
        }
        std::cout << "wrote " << n_particle << " samples to samples_bimaxwellian.txt"
               << std::endl;
    }

    std::cout << "\n=== Example 3: general_velocity_generator (f(E)=exp(-E)) ===" << std::endl;
    auto energyPdf = [](Real E) -> Real { return std::exp(-E); };
    general_velocity_generator<Real> velocityGen(energyPdf,
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
        return 1 + std::sin(x[0])*std::sin(x[1]);
    };
    static double pi = 3.14159;
    general_position_generator<Real> positionGen(2,
                                                 {-pi, -pi},
                                                 {pi, pi},
                                                 rho);
    {
        std::ofstream out("samples_general_position.txt");
        for (int i = 0; i < n_particle; ++i) {
            const position_point_type x = positionGen(gen);
            out << x[0] << " " << x[1] << " " << x[2] << "\n";
        }
        std::cout << "wrote " << n_particle << " samples to samples_general_position.txt"
               << std::endl;
    }

    return 0;
}
