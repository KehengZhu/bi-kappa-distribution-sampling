#ifndef _TEST_SUITE_H_
#define _TEST_SUITE_H_

#include <array>
#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <vector>

#include "bi_kappa_distribution.H"
#include "bi_maxwellian_distribution.H"
#include "general_position_generator.H"
#include "general_velocity_generator.H"

// All test functions are static to avoid ODR issues if this header is ever
// included from more than one translation unit.

typedef double TestReal;

// ================================================================
// test_transform
// Rotation-frame consistency tests (previously in main.cpp).
// ================================================================
static int test_transform()
{
    auto abs_real= [](TestReal x) -> TestReal { return (x < TestReal(0.0)) ? -x : x; };
    auto approx_equal= [&](TestReal a, TestReal b, TestReal tol) -> bool {
        return abs_real(a - b) <= tol;
    };
    auto dot3= [](const std::array<TestReal, 3> &a, const std::array<TestReal, 3> &b) -> TestReal {
        return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    };
    auto norm3= [](const std::array<TestReal, 3> &v) -> TestReal {
        return std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    };
    auto normalize3= [&](const std::array<TestReal, 3> &v) -> std::array<TestReal, 3> {
        TestReal n= norm3(v);
        std::array<TestReal, 3> r= {v[0] / n, v[1] / n, v[2] / n};
        return r;
    };
    auto add3= [](const std::array<TestReal, 3> &a,
                  const std::array<TestReal, 3> &b) -> std::array<TestReal, 3> {
        std::array<TestReal, 3> r= {a[0] + b[0], a[1] + b[1], a[2] + b[2]};
        return r;
    };

    std::cout << "=== test_transform ===" << std::endl;
    const TestReal tol= TestReal(1e-10);
    bool all_ok= true;

    bi_kappa_distribution<TestReal> kappaTest(TestReal(2.0), TestReal(1.0), TestReal(2.0));
    bi_maxwellian_distribution<TestReal> maxwellTest(TestReal(1.0), TestReal(2.0));

    // T1: parallel projection equals local[2]; perp magnitude preserved
    {
        bi_kappa_distribution<TestReal>::point_type local= {TestReal(0.7), TestReal(-1.2),
                                                            TestReal(2.5)};
        bi_kappa_distribution<TestReal>::point_type ub= {TestReal(1.0), TestReal(2.0),
                                                         TestReal(2.0)};
        bi_kappa_distribution<TestReal>::point_type out=
            kappaTest.rotate_from_fieldAligned_frame(local, ub);
        bi_kappa_distribution<TestReal>::point_type ubHat= normalize3(ub);

        TestReal parallel= dot3(out, ubHat);
        TestReal outNorm= norm3(out);
        TestReal perpNormSq= outNorm * outNorm - parallel * parallel;
        TestReal localPerpNormSq= local[0] * local[0] + local[1] * local[1];

        bool ok= approx_equal(parallel, local[2], TestReal(1e-9)) &&
                 approx_equal(perpNormSq, localPerpNormSq, TestReal(1e-8));
        all_ok= all_ok && ok;
        std::cout << "  T1 BiKappa component decomposition: " << (ok ? "PASS" : "FAIL")
                  << std::endl;
        if (!ok)
        {
            std::cout << "    parallel=" << parallel << " expected=" << local[2] << "\n";
            std::cout << "    perpNormSq=" << perpNormSq << " expected=" << localPerpNormSq << "\n";
        }
    }

    // T2: pure parallel local vector maps exactly along ub
    {
        bi_maxwellian_distribution<TestReal>::point_type local= {TestReal(0.0), TestReal(0.0),
                                                                  TestReal(3.0)};
        bi_maxwellian_distribution<TestReal>::point_type ub= {TestReal(0.0), TestReal(5.0),
                                                               TestReal(0.0)};
        bi_maxwellian_distribution<TestReal>::point_type out=
            maxwellTest.rotate_from_fieldAligned_frame(local, ub);

        bool ok= approx_equal(out[0], TestReal(0.0), tol) &&
                 approx_equal(out[1], TestReal(3.0), tol) &&
                 approx_equal(out[2], TestReal(0.0), tol);
        all_ok= all_ok && ok;
        std::cout << "  T2 BiMaxwellian parallel-axis mapping: " << (ok ? "PASS" : "FAIL")
                  << std::endl;
    }

    // T3: inner product is preserved under the rotation
    {
        bi_kappa_distribution<TestReal>::point_type a= {TestReal(1.2), TestReal(-0.4),
                                                        TestReal(2.7)};
        bi_kappa_distribution<TestReal>::point_type b= {TestReal(-0.8), TestReal(1.5),
                                                        TestReal(0.3)};
        bi_kappa_distribution<TestReal>::point_type ub= {TestReal(2.0), TestReal(-1.0),
                                                         TestReal(3.0)};
        bi_kappa_distribution<TestReal>::point_type ra=
            kappaTest.rotate_from_fieldAligned_frame(a, ub);
        bi_kappa_distribution<TestReal>::point_type rb=
            kappaTest.rotate_from_fieldAligned_frame(b, ub);

        bool ok= approx_equal(dot3(a, b), dot3(ra, rb), TestReal(1e-9));
        all_ok= all_ok && ok;
        std::cout << "  T3 BiKappa dot-product preservation: " << (ok ? "PASS" : "FAIL")
                  << std::endl;
    }

    // T4: linearity T(a+b) = T(a) + T(b)
    {
        bi_maxwellian_distribution<TestReal>::point_type a= {TestReal(0.4), TestReal(1.1),
                                                             TestReal(-0.2)};
        bi_maxwellian_distribution<TestReal>::point_type b= {TestReal(-1.3), TestReal(0.7),
                                                             TestReal(2.2)};
        bi_maxwellian_distribution<TestReal>::point_type ab= add3(a, b);
        bi_maxwellian_distribution<TestReal>::point_type ub= {TestReal(1.0), TestReal(3.0),
                                                              TestReal(2.0)};

        bi_maxwellian_distribution<TestReal>::point_type ta=
            maxwellTest.rotate_from_fieldAligned_frame(a, ub);
        bi_maxwellian_distribution<TestReal>::point_type tb=
            maxwellTest.rotate_from_fieldAligned_frame(b, ub);
        bi_maxwellian_distribution<TestReal>::point_type tab=
            maxwellTest.rotate_from_fieldAligned_frame(ab, ub);
        bi_maxwellian_distribution<TestReal>::point_type tSum= add3(ta, tb);

        bool ok= approx_equal(tab[0], tSum[0], TestReal(1e-9)) &&
                 approx_equal(tab[1], tSum[1], TestReal(1e-9)) &&
                 approx_equal(tab[2], tSum[2], TestReal(1e-9));
        all_ok= all_ok && ok;
        std::cout << "  T4 BiMaxwellian linearity: " << (ok ? "PASS" : "FAIL") << std::endl;
    }

    std::cout << "test_transform: " << (all_ok ? "PASS" : "FAIL") << std::endl;
    return all_ok ? 0 : 1;
}

// ================================================================
// test_define_bikappa
// Tests for the define() two-phase initialization API.
// ================================================================
static int test_define_bikappa()
{
    std::cout << "=== test_define_bikappa ===" << std::endl;
    bool all_ok= true;

    // D1: default construct, then define(individual args) – param getters agree
    {
        bi_kappa_distribution<TestReal> dist;
        dist.define(TestReal(3.0), TestReal(1.5), TestReal(2.5));
        bool ok= (dist.kappa() == TestReal(3.0)) && (dist.theta_perp() == TestReal(1.5)) &&
                 (dist.theta_par() == TestReal(2.5));
        all_ok= all_ok && ok;
        std::cout << "  D1 default construct + define(args) param check: "
                  << (ok ? "PASS" : "FAIL") << std::endl;
    }

    // D2: define(param_type) overload
    {
        bi_kappa_distribution<TestReal>::param_type p(TestReal(4.0), TestReal(0.5),
                                                       TestReal(3.0));
        bi_kappa_distribution<TestReal> dist;
        dist.define(p);
        bool ok= (dist.kappa() == TestReal(4.0)) && (dist.theta_perp() == TestReal(0.5)) &&
                 (dist.theta_par() == TestReal(3.0));
        all_ok= all_ok && ok;
        std::cout << "  D2 define(param_type): " << (ok ? "PASS" : "FAIL") << std::endl;
    }

    // D3: redefine with larger theta_par – mean |vz| should scale proportionally
    {
        const unsigned int seed= 42u;
        const int n= 2000;
        bi_kappa_distribution<TestReal> dist;

        std::mt19937 gen(seed);
        dist.define(TestReal(2.0), TestReal(1.0), TestReal(1.0));
        TestReal vz_sum_1= TestReal(0.0);
        for (int i= 0; i < n; ++i)
        {
            vz_sum_1+= std::abs(dist(gen)[2]);
        }

        gen.seed(seed);
        dist.define(TestReal(2.0), TestReal(1.0), TestReal(5.0));
        TestReal vz_sum_2= TestReal(0.0);
        for (int i= 0; i < n; ++i)
        {
            vz_sum_2+= std::abs(dist(gen)[2]);
        }

        // Expect ratio ≈ 5; accept [3, 8] to accommodate sampling variance
        TestReal ratio= (vz_sum_1 > TestReal(0.0)) ? vz_sum_2 / vz_sum_1 : TestReal(0.0);
        bool ok= (ratio > TestReal(3.0)) && (ratio < TestReal(8.0));
        all_ok= all_ok && ok;
        std::cout << "  D3 redefine scales vz (ratio=" << ratio << "): "
                  << (ok ? "PASS" : "FAIL") << std::endl;
    }

    std::cout << "test_define_bikappa: " << (all_ok ? "PASS" : "FAIL") << std::endl;
    return all_ok ? 0 : 1;
}

// ================================================================
// test_define_bimaxwellian
// ================================================================
static int test_define_bimaxwellian()
{
    std::cout << "=== test_define_bimaxwellian ===" << std::endl;
    bool all_ok= true;

    // D1: default construct + define(individual args)
    {
        bi_maxwellian_distribution<TestReal> dist;
        dist.define(TestReal(2.0), TestReal(4.0));
        bool ok= (dist.theta_perp() == TestReal(2.0)) && (dist.theta_par() == TestReal(4.0));
        all_ok= all_ok && ok;
        std::cout << "  D1 default construct + define(args) param check: "
                  << (ok ? "PASS" : "FAIL") << std::endl;
    }

    // D2: define(param_type) overload
    {
        bi_maxwellian_distribution<TestReal>::param_type p(TestReal(1.5), TestReal(3.0));
        bi_maxwellian_distribution<TestReal> dist;
        dist.define(p);
        bool ok= (dist.theta_perp() == TestReal(1.5)) && (dist.theta_par() == TestReal(3.0));
        all_ok= all_ok && ok;
        std::cout << "  D2 define(param_type): " << (ok ? "PASS" : "FAIL") << std::endl;
    }

    // D3: redefine with larger theta_perp – mean |vx| should scale proportionally
    {
        const unsigned int seed= 42u;
        const int n= 2000;
        bi_maxwellian_distribution<TestReal> dist;

        std::mt19937 gen(seed);
        dist.define(TestReal(1.0), TestReal(1.0));
        TestReal vx_sum_1= TestReal(0.0);
        for (int i= 0; i < n; ++i)
        {
            vx_sum_1+= std::abs(dist(gen)[0]);
        }

        gen.seed(seed);
        dist.define(TestReal(5.0), TestReal(1.0));
        TestReal vx_sum_2= TestReal(0.0);
        for (int i= 0; i < n; ++i)
        {
            vx_sum_2+= std::abs(dist(gen)[0]);
        }

        TestReal ratio= (vx_sum_1 > TestReal(0.0)) ? vx_sum_2 / vx_sum_1 : TestReal(0.0);
        bool ok= (ratio > TestReal(3.0)) && (ratio < TestReal(8.0));
        all_ok= all_ok && ok;
        std::cout << "  D3 redefine scales vx (ratio=" << ratio << "): "
                  << (ok ? "PASS" : "FAIL") << std::endl;
    }

    std::cout << "test_define_bimaxwellian: " << (all_ok ? "PASS" : "FAIL") << std::endl;
    return all_ok ? 0 : 1;
}

// ================================================================
// test_define_general_velocity
// ================================================================
static int test_define_general_velocity()
{
    std::cout << "=== test_define_general_velocity ===" << std::endl;
    bool all_ok= true;

    auto energyPdf= [](TestReal E) -> TestReal { return std::exp(-E); };
    std::mt19937 gen(99u);

    // D1: default construct + define + sample succeeds with finite output
    {
        bool ok= false;
        try
        {
            general_velocity_generator<TestReal> vgen;
            vgen.define(energyPdf, TestReal(0.0), TestReal(10.0), TestReal(1.0));
            auto v= vgen(gen);
            ok= std::isfinite(v[0]) && std::isfinite(v[1]) && std::isfinite(v[2]);
        }
        catch (...)
        {
        }
        all_ok= all_ok && ok;
        std::cout << "  D1 default construct + define + sample: " << (ok ? "PASS" : "FAIL")
                  << std::endl;
    }

    // D2: sampling before define() throws std::runtime_error
    {
        bool threw= false;
        try
        {
            general_velocity_generator<TestReal> vgen;
            vgen(gen);
        }
        catch (const std::runtime_error &)
        {
            threw= true;
        }
        catch (...)
        {
        }
        all_ok= all_ok && threw;
        std::cout << "  D2 sample before define throws: " << (threw ? "PASS" : "FAIL")
                  << std::endl;
    }

    // D3: redefine with different parameters, sampling still succeeds
    {
        bool ok= false;
        try
        {
            general_velocity_generator<TestReal> vgen;
            vgen.define(energyPdf, TestReal(0.0), TestReal(10.0), TestReal(1.0));
            vgen.define(energyPdf, TestReal(0.0), TestReal(5.0), TestReal(2.0));
            auto v= vgen(gen);
            ok= std::isfinite(v[0]) && std::isfinite(v[1]) && std::isfinite(v[2]);
        }
        catch (...)
        {
        }
        all_ok= all_ok && ok;
        std::cout << "  D3 redefine then sample: " << (ok ? "PASS" : "FAIL") << std::endl;
    }

    std::cout << "test_define_general_velocity: " << (all_ok ? "PASS" : "FAIL") << std::endl;
    return all_ok ? 0 : 1;
}

// ================================================================
// test_define_general_position
// ================================================================
static int test_define_general_position()
{
    std::cout << "=== test_define_general_position ===" << std::endl;
    bool all_ok= true;

    typedef general_position_generator<TestReal>::point_type pvec;
    std::mt19937 gen(77u);

    // D1: default construct + define(1D) + sample succeeds
    {
        bool ok= false;
        try
        {
            auto rho1d= [](const pvec &x) -> TestReal {
                return TestReal(1.0) + std::sin(x[0]);
            };
            general_position_generator<TestReal> pgen;
            pgen.define(1, {TestReal(0.0)}, {TestReal(6.28)}, rho1d);
            auto p= pgen(gen);
            ok= (p.size() == 1) && std::isfinite(p[0]);
        }
        catch (...)
        {
        }
        all_ok= all_ok && ok;
        std::cout << "  D1 default construct + define(1D) + sample: " << (ok ? "PASS" : "FAIL")
                  << std::endl;
    }

    // D2: sampling before define() throws std::runtime_error
    {
        bool threw= false;
        try
        {
            general_position_generator<TestReal> pgen;
            pgen(gen);
        }
        catch (const std::runtime_error &)
        {
            threw= true;
        }
        catch (...)
        {
        }
        all_ok= all_ok && threw;
        std::cout << "  D2 sample before define throws: " << (threw ? "PASS" : "FAIL")
                  << std::endl;
    }

    // D3: redefine with a different dimension, sample returns correct dimension
    {
        bool ok= false;
        try
        {
            auto rho1d= [](const pvec &x) -> TestReal {
                return TestReal(1.0) + std::sin(x[0]);
            };
            auto rho2d= [](const pvec &x) -> TestReal {
                return TestReal(1.0) + std::sin(x[0]) * std::cos(x[1]);
            };
            general_position_generator<TestReal> pgen;
            pgen.define(1, {TestReal(0.0)}, {TestReal(6.28)}, rho1d);
            pgen.define(2, {TestReal(0.0), TestReal(0.0)}, {TestReal(6.28), TestReal(6.28)},
                        rho2d);
            auto p= pgen(gen);
            ok= (p.size() == 2) && std::isfinite(p[0]) && std::isfinite(p[1]);
        }
        catch (...)
        {
        }
        all_ok= all_ok && ok;
        std::cout << "  D3 redefine to 2D then sample: " << (ok ? "PASS" : "FAIL") << std::endl;
    }

    std::cout << "test_define_general_position: " << (all_ok ? "PASS" : "FAIL") << std::endl;
    return all_ok ? 0 : 1;
}

// ================================================================
// test_max_normalized_velocity_bikappa
// ================================================================
static int test_max_normalized_velocity_bikappa()
{
    std::cout << "=== test_max_normalized_velocity_bikappa ===" << std::endl;
    bool all_ok= true;

    typedef bi_kappa_distribution<TestReal>::point_type K3;

    // V1: getter returns the value set via define()
    {
        K3 ub_default= {TestReal(0.0), TestReal(0.0), TestReal(1.0)};
        bi_kappa_distribution<TestReal> dist;
        dist.define(TestReal(2.0), TestReal(1.0), TestReal(1.0), ub_default, TestReal(3.5));
        bool ok= (dist.max_normalized_velocity() == TestReal(3.5));
        all_ok= all_ok && ok;
        std::cout << "  V1 cap getter: " << (ok ? "PASS" : "FAIL") << std::endl;
    }

    // V2: all samples satisfy the normalized cap (no rotation, so local == global frame)
    //     Cap is applied as |vi| / theta_i <= cap for each component.
    {
        const TestReal cap= TestReal(6.0);
        const TestReal theta_perp= TestReal(1.0);
        const TestReal theta_par= TestReal(2.0);
        K3 ub_default= {TestReal(0.0), TestReal(0.0), TestReal(1.0)};

        bi_kappa_distribution<TestReal> dist;
        dist.define(TestReal(2.0), theta_perp, theta_par, ub_default, cap);

        std::mt19937 gen(12345u);
        const int n= 500;
        bool all_capped= true;
        for (int i= 0; i < n; ++i)
        {
            auto v= dist(gen);
            if (std::abs(v[0]) / theta_perp > cap + TestReal(1e-10) ||
                std::abs(v[1]) / theta_perp > cap + TestReal(1e-10) ||
                std::abs(v[2]) / theta_par  > cap + TestReal(1e-10))
            {
                all_capped= false;
                std::cout << "    failed at sample " << i << ": |v|/theta = ("
                          << std::abs(v[0]) / theta_perp << ", "
                          << std::abs(v[1]) / theta_perp << ", "
                          << std::abs(v[2]) / theta_par  << ")\n";
                break;
            }
        }
        all_ok= all_ok && all_capped;
        std::cout << "  V2 " << n << " samples satisfy cap=" << cap << ": "
                  << (all_capped ? "PASS" : "FAIL") << std::endl;
    }

    // V3: cap still enforced after redefine with a tighter bound
    {
        const TestReal cap= TestReal(4.0);
        const TestReal theta_perp= TestReal(1.0);
        const TestReal theta_par= TestReal(1.0);
        K3 ub_default= {TestReal(0.0), TestReal(0.0), TestReal(1.0)};

        bi_kappa_distribution<TestReal> dist;
        dist.define(TestReal(2.0), theta_perp, theta_par);            // no cap first
        dist.define(TestReal(2.0), theta_perp, theta_par, ub_default, cap); // redefine with cap

        std::mt19937 gen(99999u);
        const int n= 500;
        bool all_capped= true;
        for (int i= 0; i < n; ++i)
        {
            auto v= dist(gen);
            if (std::abs(v[0]) / theta_perp > cap + TestReal(1e-10) ||
                std::abs(v[1]) / theta_perp > cap + TestReal(1e-10) ||
                std::abs(v[2]) / theta_par  > cap + TestReal(1e-10))
            {
                all_capped= false;
                break;
            }
        }
        all_ok= all_ok && all_capped;
        std::cout << "  V3 cap enforced after redefine (" << n << " samples, cap=" << cap
                  << "): " << (all_capped ? "PASS" : "FAIL") << std::endl;
    }

    // V4: invalid cap (<= 0) throws std::invalid_argument
    {
        bool threw= false;
        try
        {
            K3 ub_default= {TestReal(0.0), TestReal(0.0), TestReal(1.0)};
            bi_kappa_distribution<TestReal> dist;
            dist.define(TestReal(2.0), TestReal(1.0), TestReal(1.0), ub_default,
                        TestReal(-1.0));
        }
        catch (const std::invalid_argument &)
        {
            threw= true;
        }
        catch (...)
        {
        }
        all_ok= all_ok && threw;
        std::cout << "  V4 negative cap throws: " << (threw ? "PASS" : "FAIL") << std::endl;
    }

    std::cout << "test_max_normalized_velocity_bikappa: " << (all_ok ? "PASS" : "FAIL")
              << std::endl;
    return all_ok ? 0 : 1;
}

// ================================================================
// test_max_normalized_velocity_bimaxwellian
// ================================================================
static int test_max_normalized_velocity_bimaxwellian()
{
    std::cout << "=== test_max_normalized_velocity_bimaxwellian ===" << std::endl;
    bool all_ok= true;

    typedef bi_maxwellian_distribution<TestReal>::point_type M3;

    // V1: getter returns the value set via define()
    {
        M3 ub_default= {TestReal(0.0), TestReal(0.0), TestReal(1.0)};
        bi_maxwellian_distribution<TestReal> dist;
        dist.define(TestReal(1.0), TestReal(1.0), ub_default, TestReal(4.5));
        bool ok= (dist.max_normalized_velocity() == TestReal(4.5));
        all_ok= all_ok && ok;
        std::cout << "  V1 cap getter: " << (ok ? "PASS" : "FAIL") << std::endl;
    }

    // V2: all samples satisfy the normalized cap (no rotation)
    {
        const TestReal cap= TestReal(6.0);
        const TestReal theta_perp= TestReal(2.0);
        const TestReal theta_par= TestReal(1.0);
        M3 ub_default= {TestReal(0.0), TestReal(0.0), TestReal(1.0)};

        bi_maxwellian_distribution<TestReal> dist;
        dist.define(theta_perp, theta_par, ub_default, cap);

        std::mt19937 gen(54321u);
        const int n= 500;
        bool all_capped= true;
        for (int i= 0; i < n; ++i)
        {
            auto v= dist(gen);
            if (std::abs(v[0]) / theta_perp > cap + TestReal(1e-10) ||
                std::abs(v[1]) / theta_perp > cap + TestReal(1e-10) ||
                std::abs(v[2]) / theta_par  > cap + TestReal(1e-10))
            {
                all_capped= false;
                std::cout << "    failed at sample " << i << ": |v|/theta = ("
                          << std::abs(v[0]) / theta_perp << ", "
                          << std::abs(v[1]) / theta_perp << ", "
                          << std::abs(v[2]) / theta_par  << ")\n";
                break;
            }
        }
        all_ok= all_ok && all_capped;
        std::cout << "  V2 " << n << " samples satisfy cap=" << cap << ": "
                  << (all_capped ? "PASS" : "FAIL") << std::endl;
    }

    // V3: cap still enforced after redefine with a tighter bound
    {
        const TestReal cap= TestReal(4.0);
        const TestReal theta_perp= TestReal(1.0);
        const TestReal theta_par= TestReal(1.0);
        M3 ub_default= {TestReal(0.0), TestReal(0.0), TestReal(1.0)};

        bi_maxwellian_distribution<TestReal> dist;
        dist.define(theta_perp, theta_par);                           // no cap first
        dist.define(theta_perp, theta_par, ub_default, cap);          // redefine with cap

        std::mt19937 gen(11111u);
        const int n= 500;
        bool all_capped= true;
        for (int i= 0; i < n; ++i)
        {
            auto v= dist(gen);
            if (std::abs(v[0]) / theta_perp > cap + TestReal(1e-10) ||
                std::abs(v[1]) / theta_perp > cap + TestReal(1e-10) ||
                std::abs(v[2]) / theta_par  > cap + TestReal(1e-10))
            {
                all_capped= false;
                break;
            }
        }
        all_ok= all_ok && all_capped;
        std::cout << "  V3 cap enforced after redefine (" << n << " samples, cap=" << cap
                  << "): " << (all_capped ? "PASS" : "FAIL") << std::endl;
    }

    // V4: zero cap throws std::invalid_argument
    {
        bool threw= false;
        try
        {
            M3 ub_default= {TestReal(0.0), TestReal(0.0), TestReal(1.0)};
            bi_maxwellian_distribution<TestReal> dist;
            dist.define(TestReal(1.0), TestReal(1.0), ub_default, TestReal(0.0));
        }
        catch (const std::invalid_argument &)
        {
            threw= true;
        }
        catch (...)
        {
        }
        all_ok= all_ok && threw;
        std::cout << "  V4 zero cap throws: " << (threw ? "PASS" : "FAIL") << std::endl;
    }

    std::cout << "test_max_normalized_velocity_bimaxwellian: " << (all_ok ? "PASS" : "FAIL")
              << std::endl;
    return all_ok ? 0 : 1;
}

// ================================================================
// run_all_tests – call this from main()
// Returns 0 if all tests pass, 1 if any test fails.
// ================================================================
static int run_all_tests()
{
    int result= 0;
    result|= test_transform();
    result|= test_define_bikappa();
    result|= test_define_bimaxwellian();
    result|= test_define_general_velocity();
    result|= test_define_general_position();
    result|= test_max_normalized_velocity_bikappa();
    result|= test_max_normalized_velocity_bimaxwellian();
    std::cout << "\n=== Test suite: " << (result == 0 ? "ALL PASS" : "SOME FAILED") << " ==="
              << std::endl;
    return result;
}

#endif // _TEST_SUITE_H_
