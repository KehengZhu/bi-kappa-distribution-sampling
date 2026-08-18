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
#include "field_aligned_velocity_generator.H"

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

    // D3: individual setters update the canonical parameter state
    {
        bi_kappa_distribution<TestReal> dist;
        bi_kappa_distribution<TestReal>::point_type ub= {TestReal(1.0), TestReal(2.0),
                                                         TestReal(3.0)};
        dist.define(TestReal(2.0), TestReal(1.0), TestReal(1.0));
        dist.kappa(TestReal(5.0));
        dist.theta_perp(TestReal(1.5));
        dist.theta_par(TestReal(2.5));
        dist.ub(ub);
        dist.max_normalized_velocity(TestReal(7.0));

        bi_kappa_distribution<TestReal>::param_type params= dist.param();
        bool ok= (dist.kappa() == TestReal(5.0)) && (dist.theta_perp() == TestReal(1.5)) &&
                 (dist.theta_par() == TestReal(2.5)) && (dist.ub() == ub) &&
                 (dist.max_normalized_velocity() == TestReal(7.0)) &&
                 (params.kappa() == TestReal(5.0)) &&
                 (params.theta_perp() == TestReal(1.5)) &&
                 (params.theta_par() == TestReal(2.5)) && (params.ub() == ub) &&
                 (params.max_normalized_velocity() == TestReal(7.0));
        all_ok= all_ok && ok;
        std::cout << "  D3 individual setters update params: " << (ok ? "PASS" : "FAIL")
                  << std::endl;
    }

    // D4: redefine with larger theta_par – mean |vz| should scale proportionally
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
        std::cout << "  D4 redefine scales vz (ratio=" << ratio << "): "
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

    // D3: individual setters update the canonical parameter state
    {
        bi_maxwellian_distribution<TestReal> dist;
        bi_maxwellian_distribution<TestReal>::point_type ub= {TestReal(3.0), TestReal(1.0),
                                                              TestReal(4.0)};
        dist.define(TestReal(1.0), TestReal(1.0));
        dist.theta_perp(TestReal(2.5));
        dist.theta_par(TestReal(4.5));
        dist.ub(ub);
        dist.max_normalized_velocity(TestReal(8.0));

        bi_maxwellian_distribution<TestReal>::param_type params= dist.param();
        bool ok= (dist.theta_perp() == TestReal(2.5)) &&
                 (dist.theta_par() == TestReal(4.5)) && (dist.ub() == ub) &&
                 (dist.max_normalized_velocity() == TestReal(8.0)) &&
                 (params.theta_perp() == TestReal(2.5)) &&
                 (params.theta_par() == TestReal(4.5)) && (params.ub() == ub) &&
                 (params.max_normalized_velocity() == TestReal(8.0));
        all_ok= all_ok && ok;
        std::cout << "  D3 individual setters update params: " << (ok ? "PASS" : "FAIL")
                  << std::endl;
    }

    // D4: redefine with larger theta_perp – mean |vx| should scale proportionally
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
    std::cout << "  D4 redefine scales vx (ratio=" << ratio << "): "
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

    auto speedSqPdf= [](TestReal w) -> TestReal { return std::exp(-TestReal(0.5) * w); };
    std::mt19937 gen(99u);

    // D1: default construct + define + sample succeeds with finite output
    {
        bool ok= false;
        try
        {
            general_velocity_generator<TestReal> vgen;
            vgen.define(speedSqPdf, TestReal(0.0), TestReal(20.0));
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
            vgen.define(speedSqPdf, TestReal(0.0), TestReal(20.0));
            vgen.define(speedSqPdf, TestReal(0.0), TestReal(10.0));
            auto v= vgen(gen);
            ok= std::isfinite(v[0]) && std::isfinite(v[1]) && std::isfinite(v[2]);
        }
        catch (...)
        {
        }
        all_ok= all_ok && ok;
        std::cout << "  D3 redefine then sample: " << (ok ? "PASS" : "FAIL") << std::endl;
    }

    // D4: negative speed_sq_min throws std::invalid_argument
    {
        bool threw= false;
        try
        {
            general_velocity_generator<TestReal> vgen;
            vgen.define(speedSqPdf, TestReal(-1.0), TestReal(10.0));
        }
        catch (const std::invalid_argument &)
        {
            threw= true;
        }
        catch (...)
        {
        }
        all_ok= all_ok && threw;
        std::cout << "  D4 negative speed_sq_min throws: " << (threw ? "PASS" : "FAIL")
                  << std::endl;
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
// test_seed_api
// Tests the internal seeded generator: define(..., seed), seed(s),
// and the no-arg operator()().  The external operator()(gen) path
// is already exercised in the D3 and V2/V3 tests above.
// ================================================================
static int test_seed_api()
{
    std::cout << "=== test_seed_api ===" << std::endl;
    bool all_ok= true;

    // S1: bi_kappa define(..., seed) + operator()() is reproducible
    {
        const int seed= 55555;
        const int n= 100;
        bi_kappa_distribution<TestReal>::point_type ub= {TestReal(0.0), TestReal(0.0),
                                                         TestReal(1.0)};
        bi_kappa_distribution<TestReal> dist;
        dist.define(TestReal(2.0), TestReal(1.0), TestReal(2.0), ub, TestReal(20.0), seed);
        TestReal sum1= TestReal(0.0);
        for (int i= 0; i < n; ++i)
        {
            auto v= dist();
            sum1+= v[0] + v[1] + v[2];
        }

        dist.define(TestReal(2.0), TestReal(1.0), TestReal(2.0), ub, TestReal(20.0), seed);
        TestReal sum2= TestReal(0.0);
        for (int i= 0; i < n; ++i)
        {
            auto v= dist();
            sum2+= v[0] + v[1] + v[2];
        }

        bool ok= (sum1 == sum2);
        all_ok= all_ok && ok;
        std::cout << "  S1 bi_kappa define(..., seed): " << (ok ? "PASS" : "FAIL") << std::endl;
    }

    // S2: bi_kappa seed(s) reseeds the internal generator reproducibly
    {
        const int seed= 33333;
        const int n= 100;
        bi_kappa_distribution<TestReal> dist;
        dist.define(TestReal(2.0), TestReal(1.0), TestReal(2.0));

        dist.seed(seed);
        TestReal sum1= TestReal(0.0);
        for (int i= 0; i < n; ++i)
        {
            auto v= dist();
            sum1+= v[0] + v[1] + v[2];
        }

        dist.seed(seed);
        TestReal sum2= TestReal(0.0);
        for (int i= 0; i < n; ++i)
        {
            auto v= dist();
            sum2+= v[0] + v[1] + v[2];
        }

        bool ok= (sum1 == sum2);
        all_ok= all_ok && ok;
        std::cout << "  S2 bi_kappa seed(s): " << (ok ? "PASS" : "FAIL") << std::endl;
    }

    // S3: bi_maxwellian define(..., seed) + operator()() is reproducible
    {
        const int seed= 77777;
        const int n= 100;
        bi_maxwellian_distribution<TestReal>::point_type ub= {TestReal(0.0), TestReal(0.0),
                                                               TestReal(1.0)};
        bi_maxwellian_distribution<TestReal> dist;
        dist.define(TestReal(1.0), TestReal(2.0), ub, TestReal(20.0), seed);
        TestReal sum1= TestReal(0.0);
        for (int i= 0; i < n; ++i)
        {
            auto v= dist();
            sum1+= v[0] + v[1] + v[2];
        }

        dist.define(TestReal(1.0), TestReal(2.0), ub, TestReal(20.0), seed);
        TestReal sum2= TestReal(0.0);
        for (int i= 0; i < n; ++i)
        {
            auto v= dist();
            sum2+= v[0] + v[1] + v[2];
        }

        bool ok= (sum1 == sum2);
        all_ok= all_ok && ok;
        std::cout << "  S3 bi_maxwellian define(..., seed): " << (ok ? "PASS" : "FAIL")
                  << std::endl;
    }

    // S4: bi_maxwellian seed(s) reseeds reproducibly
    {
        const int seed= 88888;
        const int n= 100;
        bi_maxwellian_distribution<TestReal> dist;
        dist.define(TestReal(1.0), TestReal(2.0));

        dist.seed(seed);
        TestReal sum1= TestReal(0.0);
        for (int i= 0; i < n; ++i)
        {
            auto v= dist();
            sum1+= v[0] + v[1] + v[2];
        }

        dist.seed(seed);
        TestReal sum2= TestReal(0.0);
        for (int i= 0; i < n; ++i)
        {
            auto v= dist();
            sum2+= v[0] + v[1] + v[2];
        }

        bool ok= (sum1 == sum2);
        all_ok= all_ok && ok;
        std::cout << "  S4 bi_maxwellian seed(s): " << (ok ? "PASS" : "FAIL") << std::endl;
    }

    // S5: general_velocity_generator seed(s) + operator()() is reproducible
    {
        const int seed= 44444;
        const int n= 50;
        auto speedSqPdf= [](TestReal w) -> TestReal { return std::exp(-TestReal(0.5) * w); };
        general_velocity_generator<TestReal> vgen;
        vgen.define(speedSqPdf, TestReal(0.0), TestReal(20.0));

        vgen.seed(seed);
        TestReal sum1= TestReal(0.0);
        for (int i= 0; i < n; ++i)
        {
            auto v= vgen();
            sum1+= v[0] + v[1] + v[2];
        }

        vgen.seed(seed);
        TestReal sum2= TestReal(0.0);
        for (int i= 0; i < n; ++i)
        {
            auto v= vgen();
            sum2+= v[0] + v[1] + v[2];
        }

        bool ok= (sum1 == sum2);
        all_ok= all_ok && ok;
        std::cout << "  S5 general_velocity_generator seed(s): " << (ok ? "PASS" : "FAIL")
                  << std::endl;
    }

    // S6: general_position_generator seed(s) + operator()() is reproducible
    {
        const int seed= 22222;
        const int n= 50;
        typedef general_position_generator<TestReal>::point_type pvec;
        auto rho= [](const pvec &x) -> TestReal { return TestReal(1.0) + std::sin(x[0]); };
        general_position_generator<TestReal> pgen;
        pgen.define(1, {TestReal(0.0)}, {TestReal(6.28)}, rho);

        pgen.seed(seed);
        TestReal sum1= TestReal(0.0);
        for (int i= 0; i < n; ++i)
            sum1+= pgen()[0];

        pgen.seed(seed);
        TestReal sum2= TestReal(0.0);
        for (int i= 0; i < n; ++i)
            sum2+= pgen()[0];

        bool ok= (sum1 == sum2);
        all_ok= all_ok && ok;
        std::cout << "  S6 general_position_generator seed(s): " << (ok ? "PASS" : "FAIL")
                  << std::endl;
    }

    // S7: negative seed on bi_kappa routes through random_device; output is finite
    {
        bool ok= false;
        try
        {
            bi_kappa_distribution<TestReal> dist;
            dist.define(TestReal(2.0), TestReal(1.0), TestReal(2.0));
            dist.seed(-1);
            auto v= dist();
            ok= std::isfinite(v[0]) && std::isfinite(v[1]) && std::isfinite(v[2]);
        }
        catch (...)
        {
        }
        all_ok= all_ok && ok;
        std::cout << "  S7 negative seed (random_device) gives finite output: "
                  << (ok ? "PASS" : "FAIL") << std::endl;
    }

    // S8: two different fixed seeds produce different sequences for bi_maxwellian
    {
        const int n= 100;
        bi_maxwellian_distribution<TestReal> dist;

        dist.define(TestReal(1.0), TestReal(1.0));
        dist.seed(1111);
        TestReal sumA= TestReal(0.0);
        for (int i= 0; i < n; ++i)
        {
            auto v= dist();
            sumA+= v[0] + v[1] + v[2];
        }

        dist.seed(9999);
        TestReal sumB= TestReal(0.0);
        for (int i= 0; i < n; ++i)
        {
            auto v= dist();
            sumB+= v[0] + v[1] + v[2];
        }

        bool ok= (sumA != sumB);
        all_ok= all_ok && ok;
        std::cout << "  S8 different seeds produce different sequences: "
                  << (ok ? "PASS" : "FAIL") << std::endl;
    }

    std::cout << "test_seed_api: " << (all_ok ? "PASS" : "FAIL") << std::endl;
    return all_ok ? 0 : 1;
}

// ================================================================
// test_define_field_aligned
// Statistical checks for the field-aligned velocity generator:
// signed parallel speed along B + Maxwellian perpendicular spread.
// ================================================================
static int test_define_field_aligned()
{
    std::cout << "=== test_define_field_aligned ===" << std::endl;
    bool all_ok= true;

    typedef field_aligned_velocity_generator<TestReal>::point_type F3;
    auto speedSqPdf= [](TestReal w) -> TestReal { return std::exp(-TestReal(0.5) * w); };

    // D1: default construct + define + sample succeeds with finite output
    {
        std::mt19937 gen(99u);
        bool ok= false;
        try
        {
            field_aligned_velocity_generator<TestReal> fgen;
            fgen.define(speedSqPdf, TestReal(0.0), TestReal(40.0), TestReal(1.0),
                        F3({TestReal(0.0), TestReal(0.0), TestReal(1.0)}), +1);
            auto v= fgen(gen);
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
        std::mt19937 gen(99u);
        bool threw= false;
        try
        {
            field_aligned_velocity_generator<TestReal> fgen;
            fgen(gen);
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

    // S1: with B=z, vz is one-signed (matches +1 sign) and vx,vy Maxwellian
    //     with mean~0 and variance ~ theta_perp^2/2.
    {
        const TestReal theta_perp= TestReal(2.0);
        field_aligned_velocity_generator<TestReal> fgen;
        fgen.define(speedSqPdf, TestReal(0.0), TestReal(40.0), theta_perp,
                    F3({TestReal(0.0), TestReal(0.0), TestReal(1.0)}), +1);
        std::mt19937 gen(12345u);

        const int n= 200000;
        TestReal sum_x= 0, sum_y= 0, sum_xx= 0, sum_yy= 0, sum_xy= 0;
        bool vz_all_positive= true;
        for (int i= 0; i < n; ++i)
        {
            auto v= fgen(gen);
            sum_x+= v[0];
            sum_y+= v[1];
            sum_xx+= v[0] * v[0];
            sum_yy+= v[1] * v[1];
            sum_xy+= v[0] * v[1];
            if (v[2] <= TestReal(0.0))
            {
                vz_all_positive= false;
            }
        }
        TestReal mean_x= sum_x / n;
        TestReal mean_y= sum_y / n;
        TestReal var_x= sum_xx / n - mean_x * mean_x;
        TestReal var_y= sum_yy / n - mean_y * mean_y;
        TestReal cov_xy= sum_xy / n - mean_x * mean_y;
        TestReal expected_var= theta_perp * theta_perp / TestReal(2.0);

        bool sign_ok= vz_all_positive;
        bool mean_ok= std::abs(mean_x) < TestReal(0.05) && std::abs(mean_y) < TestReal(0.05);
        bool var_ok= std::abs(var_x - expected_var) < TestReal(0.1) &&
                     std::abs(var_y - expected_var) < TestReal(0.1);
        bool iso_ok= std::abs(var_x - var_y) < TestReal(0.1) &&
                     std::abs(cov_xy) < TestReal(0.1);

        all_ok= all_ok && sign_ok && mean_ok && var_ok && iso_ok;
        std::cout << "  S1 vz>0 for sign=+1: " << (sign_ok ? "PASS" : "FAIL") << std::endl;
        std::cout << "  S1 perp mean~0: " << (mean_ok ? "PASS" : "FAIL") << " (" << mean_x << ", "
                  << mean_y << ")" << std::endl;
        std::cout << "  S1 perp var~theta_perp^2/2=" << expected_var << ": "
                  << (var_ok ? "PASS" : "FAIL") << " (" << var_x << ", " << var_y << ")"
                  << std::endl;
        std::cout << "  S1 perp isotropy (var_x~var_y, cov~0): " << (iso_ok ? "PASS" : "FAIL")
                  << std::endl;
    }

    // S2: sign=-1 flips the parallel direction (vz all negative for B=z).
    {
        field_aligned_velocity_generator<TestReal> fgen;
        fgen.define(speedSqPdf, TestReal(0.0), TestReal(40.0), TestReal(1.0),
                    F3({TestReal(0.0), TestReal(0.0), TestReal(1.0)}), -1);
        std::mt19937 gen(2024u);
        bool vz_all_negative= true;
        for (int i= 0; i < 20000; ++i)
        {
            auto v= fgen(gen);
            if (v[2] >= TestReal(0.0))
            {
                vz_all_negative= false;
                break;
            }
        }
        all_ok= all_ok && vz_all_negative;
        std::cout << "  S2 vz<0 for sign=-1: " << (vz_all_negative ? "PASS" : "FAIL")
                  << std::endl;
    }

    // S3: with a tilted B, the perpendicular part is orthogonal to B, i.e. the
    //     parallel projection v.b_hat is one-signed and removing it leaves
    //     v_perp with v_perp.B ~ 0 for every sample.
    {
        F3 B= {TestReal(1.0), TestReal(1.0), TestReal(1.0)};
        TestReal Bnorm= std::sqrt(B[0] * B[0] + B[1] * B[1] + B[2] * B[2]);
        F3 bhat= {B[0] / Bnorm, B[1] / Bnorm, B[2] / Bnorm};

        field_aligned_velocity_generator<TestReal> fgen;
        fgen.define(speedSqPdf, TestReal(0.0), TestReal(40.0), TestReal(1.5), B, +1);
        std::mt19937 gen(7u);

        bool ortho_ok= true;
        bool proj_positive= true;
        for (int i= 0; i < 20000; ++i)
        {
            auto v= fgen(gen);
            TestReal vpar= v[0] * bhat[0] + v[1] * bhat[1] + v[2] * bhat[2];
            F3 vperp= {v[0] - vpar * bhat[0], v[1] - vpar * bhat[1], v[2] - vpar * bhat[2]};
            TestReal vperp_dot_B= vperp[0] * B[0] + vperp[1] * B[1] + vperp[2] * B[2];
            if (std::abs(vperp_dot_B) > TestReal(1e-9) * (TestReal(1.0) + std::abs(vpar) * Bnorm))
            {
                ortho_ok= false;
            }
            if (vpar <= TestReal(0.0))
            {
                proj_positive= false;
            }
        }
        all_ok= all_ok && ortho_ok && proj_positive;
        std::cout << "  S3 v_perp orthogonal to tilted B: " << (ortho_ok ? "PASS" : "FAIL")
                  << std::endl;
        std::cout << "  S3 v.b_hat>0 for sign=+1: " << (proj_positive ? "PASS" : "FAIL")
                  << std::endl;
    }

    // E1: invalid arguments throw std::invalid_argument.
    {
        int threw= 0;
        int expected= 5;
        // zero B
        try { field_aligned_velocity_generator<TestReal> f;
            f.define(speedSqPdf, TestReal(0.0), TestReal(1.0), TestReal(1.0),
                     F3({TestReal(0.0), TestReal(0.0), TestReal(0.0)}), +1); }
        catch (const std::invalid_argument &) { ++threw; } catch (...) {}
        // theta_perp <= 0
        try { field_aligned_velocity_generator<TestReal> f;
            f.define(speedSqPdf, TestReal(0.0), TestReal(1.0), TestReal(0.0),
                     F3({TestReal(0.0), TestReal(0.0), TestReal(1.0)}), +1); }
        catch (const std::invalid_argument &) { ++threw; } catch (...) {}
        // speed_sq_min >= speed_sq_max
        try { field_aligned_velocity_generator<TestReal> f;
            f.define(speedSqPdf, TestReal(1.0), TestReal(1.0), TestReal(1.0),
                     F3({TestReal(0.0), TestReal(0.0), TestReal(1.0)}), +1); }
        catch (const std::invalid_argument &) { ++threw; } catch (...) {}
        // speed_sq_min < 0
        try { field_aligned_velocity_generator<TestReal> f;
            f.define(speedSqPdf, TestReal(-1.0), TestReal(1.0), TestReal(1.0),
                     F3({TestReal(0.0), TestReal(0.0), TestReal(1.0)}), +1); }
        catch (const std::invalid_argument &) { ++threw; } catch (...) {}
        // parallel_sign not +/-1
        try { field_aligned_velocity_generator<TestReal> f;
            f.define(speedSqPdf, TestReal(0.0), TestReal(1.0), TestReal(1.0),
                     F3({TestReal(0.0), TestReal(0.0), TestReal(1.0)}), 2); }
        catch (const std::invalid_argument &) { ++threw; } catch (...) {}

        bool ok= (threw == expected);
        all_ok= all_ok && ok;
        std::cout << "  E1 invalid args throw (" << threw << "/" << expected << "): "
                  << (ok ? "PASS" : "FAIL") << std::endl;
    }

    std::cout << "test_define_field_aligned: " << (all_ok ? "PASS" : "FAIL") << std::endl;
    return all_ok ? 0 : 1;
}

// ================================================================
// test_radius_formation
//
// Regression tests for the radial formation in bi_kappa_distribution.
//
// The radius is R = sqrt(X1/X2) with X1 ~ Ga(3/2,1) and X2 ~ Ga(kappa-1/2,1).
// The two square roots must be taken SEPARATELY:
//
//     sqrt(x1) / sqrt(x2)        not        sqrt(x1 / x2)
//
// The two are algebraically identical but not numerically identical.  For small
// kappa the shape of X2 is far below 1, X2 is frequently driven into the
// subnormal range, and the intermediate quotient x1/x2 overflows to +inf even
// where the final radius R is comfortably representable.  Taking the roots first
// halves the exponent before the division and keeps those draws.
//
// Measured effect (experiments/exp4_precision, 5e6 draws per configuration,
// identical variates through both formations, libc++ and libstdc++ agreeing to
// within seed noise):
//
//     float,  kappa = 0.55   loss 1.22e-2 -> 5.7e-3   (fix recovers 6.5e-3)
//     float,  kappa = 0.60   loss 1.5e-4  -> 2.9e-5   (about 5x)
//     double, kappa = 0.51   loss 8.2e-4  -> 5.7e-4
//
// A1 is a direct, deterministic check of the overflow mechanism itself; A2 and A3
// guard the released sampler at the two ends of the range.
// ================================================================
static int test_radius_formation()
{
    bool all_ok= true;
    std::cout << "\n--- test_radius_formation ---" << std::endl;

    // ---- A1: the mechanism, deterministically ----------------------------
    // Pick x2 subnormal enough that x1/x2 overflows in float while
    // sqrt(x1)/sqrt(x2) does not.  No RNG involved.
    {
        const float x1= 1.0f;
        const float x2= std::numeric_limits<float>::denorm_min(); // ~1.4e-45
        const float r_ratio= std::sqrt(x1 / x2);
        const float r_split= std::sqrt(x1) / std::sqrt(x2);

        const bool ratio_overflows= !std::isfinite(r_ratio);
        const bool split_survives= std::isfinite(r_split) && r_split > 0.0f;
        const bool ok= ratio_overflows && split_survives;
        all_ok= all_ok && ok;
        std::cout << "  A1 sqrt(x1)/sqrt(x2) survives where sqrt(x1/x2) overflows"
                  << " (ratio=" << r_ratio << ", split=" << r_split << "): "
                  << (ok ? "PASS" : "FAIL") << std::endl;
    }

    // ---- A2: the released sampler is clean in double at kappa >= 0.55 ----
    // Experiment 4 measures exactly zero non-finite draws over 5e6 draws per
    // configuration for double at every kappa >= 0.55, on both standard
    // libraries.  Any regression here means the formation was changed back.
    {
        const int n= 200000;
        int nonfinite= 0;
        bi_kappa_distribution<double> dist(0.55, 1.0, 2.0, {0.0, 0.0, 1.0},
                                           bi_kappa_distribution<double>::no_cap());
        dist.seed(4001);
        for (int i= 0; i < n; ++i)
        {
            const auto v= dist();
            if (!std::isfinite(v[0]) || !std::isfinite(v[1]) || !std::isfinite(v[2]))
                ++nonfinite;
        }
        const bool ok= (nonfinite == 0);
        all_ok= all_ok && ok;
        std::cout << "  A2 double, kappa=0.55, " << n << " draws, non-finite=" << nonfinite
                  << " (expect 0): " << (ok ? "PASS" : "FAIL") << std::endl;
    }

    // ---- A3: float at kappa = 0.6 stays under the measured bound ---------
    // Measured non-finite fraction with the split formation is 2.9e-5 (libc++)
    // and 3.1e-5 (libstdc++); with the ratio formation it is 1.5e-4.  The bound
    // below sits between the two, so it fails if the formation regresses and
    // passes with ordinary seed-to-seed variation.
    {
        const int n= 500000;
        int nonfinite= 0;
        bi_kappa_distribution<float> dist(0.6f, 1.0f, 1.0f, {0.0f, 0.0f, 1.0f},
                                          bi_kappa_distribution<float>::no_cap());
        dist.seed(4001);
        for (int i= 0; i < n; ++i)
        {
            const auto v= dist();
            if (!std::isfinite(v[0]) || !std::isfinite(v[1]) || !std::isfinite(v[2]))
                ++nonfinite;
        }
        const double frac= static_cast<double>(nonfinite) / n;
        const bool ok= (frac < 8.0e-5);
        all_ok= all_ok && ok;
        std::cout << "  A3 float, kappa=0.6, non-finite fraction=" << frac
                  << " (bound 8e-5, ratio-formation value 1.5e-4): "
                  << (ok ? "PASS" : "FAIL") << std::endl;
    }

    // ---- A4: the fix does not perturb the law at ordinary kappa ----------
    // At kappa = 2 neither formation can overflow, so the two must describe the
    // same real number.  The right quantity to bound is the MAGNITUDE of their
    // disagreement, not its frequency: `sqrt(x1/x2)` rounds twice and
    // `sqrt(x1)/sqrt(x2)` rounds three times, at different points, so the two
    // land on different ulps for a large fraction of inputs (measured: 37%)
    // while never differing by more than a couple of ulp.  The frequency
    // therefore carries no information about correctness and the magnitude does.
    //
    // Measured over 2e5 draws at kappa = 2: max relative difference 4.34e-16,
    // i.e. 1.95 * DBL_EPSILON, max 2 ulp.  The bound below is 4 ulp-equivalent,
    // which is tight enough to catch any real change of law and loose enough to
    // be insensitive to the order of operations a compiler chooses.
    {
        const int n= 200000;
        std::mt19937 gen(4001);
        std::gamma_distribution<double> gX1(1.5, 1.0);
        std::gamma_distribution<double> gX2(2.0 - 0.5, 1.0);
        double max_rel= 0.0;
        for (int i= 0; i < n; ++i)
        {
            const double x1= gX1(gen);
            const double x2= gX2(gen);
            const double a= std::sqrt(x1 / x2);
            const double b= std::sqrt(x1) / std::sqrt(x2);
            const double rel= std::abs(a - b) / std::abs(a);
            if (rel > max_rel)
                max_rel= rel;
        }
        const double bound= 4.0 * std::numeric_limits<double>::epsilon();
        const bool ok= (max_rel <= bound);
        all_ok= all_ok && ok;
        std::cout << "  A4 kappa=2, max relative difference between formations=" << max_rel
                  << " (bound " << bound << " = 4 eps): " << (ok ? "PASS" : "FAIL")
                  << std::endl;
    }

    std::cout << "test_radius_formation: " << (all_ok ? "PASS" : "FAIL") << std::endl;
    return all_ok ? 0 : 1;
}

// ================================================================
// test_no_cap_semantics
//
// no_cap() must select the untruncated law and must not run the
// acceptance-rejection loop.  A finite cap selects a different (conditional)
// law.  These are different distributions and the API must keep them distinct.
// ================================================================
static int test_no_cap_semantics()
{
    bool all_ok= true;
    std::cout << "\n--- test_no_cap_semantics ---" << std::endl;

    // ---- B1: no_cap() is accepted and reports itself as uncapped ---------
    {
        bool ok= false;
        try
        {
            bi_kappa_distribution<double> dist(2.0, 1.0, 2.0, {0.0, 0.0, 1.0},
                                               bi_kappa_distribution<double>::no_cap());
            ok= !dist.param().capped() &&
                !std::isfinite(dist.max_normalized_velocity());
        }
        catch (...)
        {
            ok= false;
        }
        all_ok= all_ok && ok;
        std::cout << "  B1 no_cap() constructs and reports capped()==false: "
                  << (ok ? "PASS" : "FAIL") << std::endl;
    }

    // ---- B2: a finite cap reports itself as capped ------------------------
    {
        bi_kappa_distribution<double> dist(2.0, 1.0, 2.0, {0.0, 0.0, 1.0}, 5.0);
        const bool ok= dist.param().capped();
        all_ok= all_ok && ok;
        std::cout << "  B2 finite cap reports capped()==true: " << (ok ? "PASS" : "FAIL")
                  << std::endl;
    }

    // ---- B3: uncapped mode really is untruncated --------------------------
    // With cap = 5 no normalized component may exceed 5.  Uncapped, some must.
    {
        const int n= 200000;
        bi_kappa_distribution<double> dist(1.0, 1.0, 1.0, {0.0, 0.0, 1.0},
                                           bi_kappa_distribution<double>::no_cap());
        dist.seed(4002);
        int beyond= 0;
        for (int i= 0; i < n; ++i)
        {
            const auto v= dist();
            if (std::abs(v[0]) > 5.0 || std::abs(v[1]) > 5.0 || std::abs(v[2]) > 5.0)
                ++beyond;
        }
        const bool ok= (beyond > 0);
        all_ok= all_ok && ok;
        std::cout << "  B3 uncapped draws exceed the box (" << beyond << "/" << n
                  << " beyond 5): " << (ok ? "PASS" : "FAIL") << std::endl;
    }

    // ---- B4: capped mode respects the box exactly -------------------------
    {
        const int n= 100000;
        const double lambda= 5.0;
        const double theta_perp= 1.0, theta_par= 2.0;
        bi_kappa_distribution<double> dist(1.0, theta_perp, theta_par, {0.0, 0.0, 1.0},
                                           lambda);
        dist.seed(4003);
        int violations= 0;
        for (int i= 0; i < n; ++i)
        {
            const auto v= dist();
            if (std::abs(v[0]) / theta_perp > lambda ||
                std::abs(v[1]) / theta_perp > lambda ||
                std::abs(v[2]) / theta_par > lambda)
                ++violations;
        }
        const bool ok= (violations == 0);
        all_ok= all_ok && ok;
        std::cout << "  B4 capped draws never leave the box (" << violations
                  << " violations): " << (ok ? "PASS" : "FAIL") << std::endl;
    }

    // ---- B5: the DEFAULT is uncapped, on every construction path ----------
    // This is the regression guard for the default flip.  A silently capped
    // default means a caller who omits the argument samples the conditional
    // box law while believing they sampled bi-Kappa -- which at kappa <= 3/2
    // is neither a good approximation nor even axisymmetric about B
    // (Experiment 2), and which hides non-finite draws (Experiment 4).
    {
        bi_kappa_distribution<double> dflt;                          // default ctor
        bi_kappa_distribution<double> ctor(2.0, 1.0, 2.0);           // value ctor
        bi_kappa_distribution<double> defined;
        defined.define(2.0, 1.0, 2.0);                               // define(...)
        bi_kappa_distribution<double> seeded;
        seeded.define(bi_kappa_distribution<double>::param_type(2.0, 1.0, 2.0), 7);
        bi_kappa_distribution<double>::param_type p;                 // param_type default

        const bool ok= !dflt.param().capped() && !ctor.param().capped() &&
                       !defined.param().capped() && !seeded.param().capped() &&
                       !p.capped();
        all_ok= all_ok && ok;
        std::cout << "  B5 default is uncapped on every construction path: "
                  << (ok ? "PASS" : "FAIL") << std::endl;
    }

    // ---- B6: a defaulted sampler really produces untruncated draws --------
    // B5 checks the reported flag; this checks the behaviour behind it.
    {
        const int n= 200000;
        bi_kappa_distribution<double> dist(1.0, 1.0, 1.0);   // cap omitted
        dist.seed(4002);
        int beyond= 0;
        for (int i= 0; i < n; ++i)
        {
            const auto v= dist();
            if (std::abs(v[0]) > 5.0 || std::abs(v[1]) > 5.0 || std::abs(v[2]) > 5.0)
                ++beyond;
        }
        const bool ok= (beyond > 0);
        all_ok= all_ok && ok;
        std::cout << "  B6 defaulted sampler draws exceed the box (" << beyond << "/" << n
                  << " beyond 5): " << (ok ? "PASS" : "FAIL") << std::endl;
    }

    // ---- B7: omitting the cap is identical to naming no_cap() -------------
    // Same seed, so the two must agree draw for draw, not just in law.
    {
        const int n= 20000;
        bi_kappa_distribution<double> omitted(1.5, 1.0, 2.0);
        bi_kappa_distribution<double> named(1.5, 1.0, 2.0, {0.0, 0.0, 1.0},
                                            bi_kappa_distribution<double>::no_cap());
        omitted.seed(4004);
        named.seed(4004);
        int mismatches= 0;
        for (int i= 0; i < n; ++i)
        {
            const auto a= omitted();
            const auto b= named();
            if (a[0] != b[0] || a[1] != b[1] || a[2] != b[2])
                ++mismatches;
        }
        const bool ok= (mismatches == 0);
        all_ok= all_ok && ok;
        std::cout << "  B7 omitted cap == explicit no_cap() bitwise (" << mismatches
                  << " mismatches): " << (ok ? "PASS" : "FAIL") << std::endl;
    }

    // ---- B8: the two target laws are measurably different ------------------
    // Guards against a future "simplification" that quietly makes the cap a
    // no-op.  At kappa = 0.75, lambda = 3 the cap rejects ~55% of attempts
    // (Experiment 2), so the capped sample's spread must be visibly smaller.
    {
        const int n= 50000;
        const double kappa= 0.75, lambda= 3.0;
        bi_kappa_distribution<double> uncapped(kappa, 1.0, 1.0);
        bi_kappa_distribution<double> capped(kappa, 1.0, 1.0, {0.0, 0.0, 1.0}, lambda);
        uncapped.seed(4005);
        capped.seed(4005);
        double max_uncapped= 0.0, max_capped= 0.0;
        for (int i= 0; i < n; ++i)
        {
            const auto a= uncapped();
            const auto b= capped();
            const double ra= std::max(std::abs(a[0]), std::max(std::abs(a[1]),
                                                               std::abs(a[2])));
            const double rb= std::max(std::abs(b[0]), std::max(std::abs(b[1]),
                                                               std::abs(b[2])));
            if (std::isfinite(ra) && ra > max_uncapped) max_uncapped= ra;
            if (rb > max_capped) max_capped= rb;
        }
        // The capped run is bounded by lambda by construction; the uncapped run
        // at kappa = 0.75 must reach far past it.
        const bool ok= (max_capped <= lambda) && (max_uncapped > 10.0 * lambda);
        all_ok= all_ok && ok;
        std::cout << "  B8 capped and uncapped are different laws (max |v_i|/theta: "
                  << max_capped << " vs " << max_uncapped << "): "
                  << (ok ? "PASS" : "FAIL") << std::endl;
    }

    std::cout << "test_no_cap_semantics: " << (all_ok ? "PASS" : "FAIL") << std::endl;
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
    result|= test_define_field_aligned();
    result|= test_max_normalized_velocity_bikappa();
    result|= test_max_normalized_velocity_bimaxwellian();
    result|= test_seed_api();
    result|= test_radius_formation();
    result|= test_no_cap_semantics();
    std::cout << "\n=== Test suite: " << (result == 0 ? "ALL PASS" : "SOME FAILED") << " ==="
              << std::endl;
    return result;
}

#endif // _TEST_SUITE_H_
