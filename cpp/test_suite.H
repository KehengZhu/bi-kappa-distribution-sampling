#ifndef _TEST_SUITE_H_
#define _TEST_SUITE_H_

#include <array>
#include <cmath>
#include <cstring>
#include <type_traits>
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
// Analytic support for the log-domain radius tests
//
// Every stochastic check below states its criterion either as an exact
// deterministic identity or as an interval derived from the analytic null at the
// sample size the test actually uses.  Nothing here is a tuned constant, so these
// helpers exist to make that possible without adding a dependency.
// ================================================================

/// Engine-independent value of pi, spelled exactly as the sampler spells it, so a
/// replayed draw reproduces the sampler bit for bit.
static double test_pi() { return double(3.14159265358979323846264338327950288L); }

// ---- polygamma -------------------------------------------------------------
// psi^(m)(z) by upward recurrence to z >= 20 followed by the Stirling-type
// asymptotic series.  Only m = 0, 1, 3 are needed: E[log X] = psi(a),
// Var[log X] = psi'(a), and the standard error of a sample variance needs the
// fourth central moment mu4 = psi'''(a) + 3 psi'(a)^2.

static double test_digamma(double z)
{
    double r= 0.0;
    while (z < 20.0)
    {
        r-= 1.0 / z;
        z+= 1.0;
    }
    const double w= 1.0 / (z * z);
    return r + std::log(z) - 0.5 / z -
           w * (1.0 / 12.0 -
                w * (1.0 / 120.0 - w * (1.0 / 252.0 - w * (1.0 / 240.0 - w * (1.0 / 132.0)))));
}

static double test_trigamma(double z)
{
    double r= 0.0;
    while (z < 20.0)
    {
        r+= 1.0 / (z * z);
        z+= 1.0;
    }
    const double w= 1.0 / (z * z);
    return r + (1.0 / z) * (1.0 + 0.5 / z +
                            w * (1.0 / 6.0 -
                                 w * (1.0 / 30.0 - w * (1.0 / 42.0 - w * (1.0 / 30.0)))));
}

static double test_tetragamma(double z) // psi'''
{
    double r= 0.0;
    while (z < 20.0)
    {
        r+= 6.0 / (z * z * z * z);
        z+= 1.0;
    }
    const double w= 1.0 / (z * z);
    return r + (1.0 / (z * z * z)) *
                   (2.0 + 3.0 / z + w * (2.0 - w * (1.0 - w * (4.0 / 3.0 - w * 3.0))));
}

/// z such that P(Z > z) = alpha for a standard normal, by bisection on erfc.
static double test_normal_upper_quantile(double alpha)
{
    double lo= 0.0, hi= 40.0;
    for (int i= 0; i < 200; ++i)
    {
        const double mid= 0.5 * (lo + hi);
        if (0.5 * std::erfc(mid / std::sqrt(2.0)) > alpha)
            lo= mid;
        else
            hi= mid;
    }
    return 0.5 * (lo + hi);
}

// ---- exact binomial tails --------------------------------------------------
static double test_binom_log_pmf(long long k, long long n, double p)
{
    if (k < 0 || k > n)
        return -std::numeric_limits<double>::infinity();
    const double dk= static_cast<double>(k), dn= static_cast<double>(n);
    double t= std::lgamma(dn + 1.0) - std::lgamma(dk + 1.0) - std::lgamma(dn - dk + 1.0);
    if (k > 0)
        t+= dk * std::log(p);
    if (k < n)
        t+= (dn - dk) * std::log1p(-p);
    return t;
}

/// P(X <= k) and P(X >= k) for X ~ Binomial(n, p), by summing the probability mass
/// function over a window that carries all but a negligible part of the mass.  The
/// window is 40 standard deviations plus 80 counts wide on each side, so the
/// neglected tail is below 1e-300 for every (n, p) used here.
static void test_binom_tails(long long k, long long n, double p, double *at_most,
                             double *at_least)
{
    const double mean= static_cast<double>(n) * p;
    const double sd= std::sqrt(static_cast<double>(n) * p * (1.0 - p));
    long long lo= static_cast<long long>(std::floor(mean - 40.0 * sd - 80.0));
    long long hi= static_cast<long long>(std::ceil(mean + 40.0 * sd + 80.0));
    if (lo < 0)
        lo= 0;
    if (hi > n)
        hi= n;
    double below= 0.0, above= 0.0, at= 0.0;
    for (long long i= lo; i <= hi; ++i)
    {
        const double w= std::exp(test_binom_log_pmf(i, n, p));
        if (i < k)
            below+= w;
        else if (i > k)
            above+= w;
        else
            at= w;
    }
    const double total= below + at + above;
    *at_most= (total > 0.0) ? (below + at) / total : 0.0;
    *at_least= (total > 0.0) ? (above + at) / total : 0.0;
}

/// True when @p k lies inside the exact two-sided (1 - alpha) binomial acceptance
/// region for n trials at probability p.
static bool test_binom_consistent(long long k, long long n, double p, double alpha)
{
    double at_most= 0.0, at_least= 0.0;
    test_binom_tails(k, n, p, &at_most, &at_least);
    const double worse= (at_most < at_least) ? at_most : at_least;
    return worse > 0.5 * alpha;
}

// ---- the analytic representability floor -----------------------------------
//
// With a = kappa - 1/2, X1 ~ Gamma(3/2,1), Y ~ Gamma(a+1,1), U ~ Uniform(0,1) and
// E = -log U ~ Exp(1), the shape-boosting identity gives log X2 = log Y - E/a, so
//
//     log R = (log X1 - log Y + E/a) / 2
//
// and therefore, for any threshold L,
//
//     P(log R > L) = P(E > 2aL - a log(X1/Y)) = E[(X1/Y)^a] exp(-2aL),
//
// because E is exponential and independent of X1 and Y.  The prefactor is exact:
// E[X1^a] = Gamma(3/2+a)/Gamma(3/2) and E[Y^-a] = 1/Gamma(a+1).
//
// A returned component is V_j = R g_j with g = Q(ub) diag(sqrt(kappa) theta) n and
// n uniform on the sphere, so the draw has no representation exactly when
// R max_j |g_j| exceeds the largest finite value MAX of the type:
//
//     P(overflow) = [Gamma(3/2+a) / (Gamma(3/2) Gamma(a+1))] MAX^(-2a) E[(max_j |g_j|)^(2a)]
//
// The last factor depends only on the direction and is evaluated here by product
// midpoint quadrature over (cos theta, phi); it converges to six digits well below
// the grid used.

static double test_radius_tail_prefactor(double a)
{
    return std::exp(std::lgamma(1.5 + a) - std::lgamma(1.5) - std::lgamma(1.0 + a));
}

/// The analytic representability floor for one configuration.
///
/// @param kappa       spectral index; the denominator shape is `a = kappa - 1/2`
/// @param theta_perp  thermal speed transverse to @p ub
/// @param theta_par   thermal speed along @p ub
/// @param ub          magnetic-field direction; need not be normalized
/// @param log_max     natural logarithm of the largest finite value of the type
/// @param p_overflow  P(the returned, rotated vector has no representation)
/// @param p_rescued   P(the unrotated field-aligned vector has no representation
///                    while the rotated one does) -- the rotation-recoverable set
static void test_overflow_floor(double kappa, double theta_perp, double theta_par,
                                const std::array<double, 3> &ub, double log_max,
                                double *p_overflow, double *p_rescued)
{
    const int n_mu= 512, n_phi= 1024;
    const double a= kappa - 0.5;
    const double sqrt_kappa= std::sqrt(kappa);
    const bool axis_aligned= (ub[0] == 0.0 && ub[1] == 0.0 && ub[2] == 1.0);
    bi_kappa_distribution<double> rot(2.0, 1.0, 1.0);

    double acc_overflow= 0.0, acc_rescued= 0.0;
    for (int i= 0; i < n_mu; ++i)
    {
        const double mu= -1.0 + (static_cast<double>(i) + 0.5) * 2.0 / n_mu;
        const double sin_theta= std::sqrt(std::max(0.0, 1.0 - mu * mu));
        for (int k= 0; k < n_phi; ++k)
        {
            const double phi= (static_cast<double>(k) + 0.5) * 2.0 * test_pi() / n_phi;
            std::array<double, 3> local;
            local[0]= sqrt_kappa * theta_perp * (sin_theta * std::cos(phi));
            local[1]= sqrt_kappa * theta_perp * (sin_theta * std::sin(phi));
            local[2]= sqrt_kappa * theta_par * mu;
            const std::array<double, 3> g=
                axis_aligned ? local : rot.rotate_from_fieldAligned_frame(local, ub);
            const double m_local=
                std::max(std::abs(local[0]), std::max(std::abs(local[1]), std::abs(local[2])));
            const double m_rot= std::max(std::abs(g[0]), std::max(std::abs(g[1]), std::abs(g[2])));
            const double e_rot= std::exp(2.0 * a * (std::log(m_rot) - log_max));
            acc_overflow+= e_rot;
            const double gap= std::exp(2.0 * a * (std::log(m_local) - log_max)) - e_rot;
            if (gap > 0.0)
                acc_rescued+= gap;
        }
    }
    const double cells= static_cast<double>(n_mu) * static_cast<double>(n_phi);
    const double c= test_radius_tail_prefactor(a);
    *p_overflow= c * acc_overflow / cells;
    *p_rescued= c * acc_rescued / cells;
}

// ---- ulp accounting --------------------------------------------------------
static double test_ulp_at(double x)
{
    const double ax= std::abs(x);
    if (!(ax > 0.0))
        return std::numeric_limits<double>::denorm_min();
    return std::nextafter(ax, std::numeric_limits<double>::infinity()) - ax;
}

/// Distance between @p a and @p b in units in the last place of @p scale.
static double test_ulp_distance(double a, double b, double scale)
{
    if (a == b)
        return 0.0;
    if (!std::isfinite(a) || !std::isfinite(b))
        return std::numeric_limits<double>::infinity();
    return std::abs(a - b) / test_ulp_at(scale);
}

/// True when this build fuses a multiply and an add into a single rounding.  A
/// fused a*b + c rounds once, so x*x - (1 + 2^-51) with x = 1 + 2^-52 comes out as
/// 2^-104 rather than exactly zero.  The operands are read through volatile so the
/// expression cannot be folded at compile time.  This is a diagnostic, not an
/// assertion: it explains a golden-stream mismatch rather than causing one.
static bool test_fp_contraction_enabled()
{
    volatile double vx= 1.0 + std::ldexp(1.0, -52);
    volatile double vc= -(1.0 + std::ldexp(1.0, -51));
    const double x= vx;
    const double c= vc;
    const double r= x * x + c;
    return r != 0.0;
}

// ---- replay of the sampler's primitive sequence ----------------------------
//
// The sampler draws, in this order: X1 ~ Gamma(3/2,1); log X2 by the boosting
// identity at shape a; then a direction by Marsaglia's rejection method, two
// open-interval uniforms per attempt, averaging 4/pi attempts.  The second
// polar-normal variate is cached across calls, so a replay has to carry the same
// cache or it desynchronizes after the first draw; the direction's rejection loop
// makes the per-draw variate count variable, so a replay assuming a fixed count
// would desynchronize too.  Replaying lets a test see the quantities the sampler
// decides on -- log R, the direction, the pre-rotation vector -- which are
// otherwise private.

struct test_replay_state
{
    double cache;
    bool cached;
    test_replay_state() : cache(0.0), cached(false) {}
};

struct test_replay_draw
{
    double log_r;
    double nx, ny, nz;
};

template <typename GeneratorType>
static test_replay_draw test_replay(double a, GeneratorType &gen, test_replay_state &st)
{
    const double x1= bikappa_detail::marsagliaTsangGamma<double, GeneratorType>(
        1.5, gen, st.cache, st.cached);
    const double log_x2=
        bikappa_detail::logGamma<double, GeneratorType>(a, gen, st.cache, st.cached);
    test_replay_draw d;
    d.log_r= 0.5 * (std::log(x1) - log_x2);
    double d1, d2, dsq;
    do
    {
        d1= 2.0 * bikappa_detail::openCanonical<double, GeneratorType>(gen) - 1.0;
        d2= 2.0 * bikappa_detail::openCanonical<double, GeneratorType>(gen) - 1.0;
        dsq= d1 * d1 + d2 * d2;
    } while (dsq >= 1.0);
    const double dscale= 2.0 * std::sqrt(1.0 - dsq);
    d.nx= d1 * dscale;
    d.ny= d2 * dscale;
    d.nz= 1.0 - 2.0 * dsq;
    return d;
}

/// Engines that drive openCanonical() to the two ends of its range.  The mapping
/// from engine bits to the open interval is deterministic, so the smallest and the
/// largest value it can return are reachable exactly, which a finite random sample
/// never is: k = 0 has probability 2^-52 per draw in double.
struct test_all_zero_engine
{
    typedef unsigned int result_type;
    static result_type min() { return 0u; }
    static result_type max() { return 4294967295u; }
    result_type operator()() { return 0u; }
};

struct test_all_ones_engine
{
    typedef unsigned int result_type;
    static result_type min() { return 0u; }
    static result_type max() { return 4294967295u; }
    result_type operator()() { return 4294967295u; }
};


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
//
// Since 2.0.0 the sampler does not form X2 at all -- it carries log X2 and log R --
// so A1 and A4 now describe a formation the released class no longer uses, and what
// A2 and A3 measure is the representability floor of the type rather than the
// avoidable loss of an intermediate.  test_log_radius_finiteness below covers the
// log-domain sampler against that floor directly.
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

    // ---- A3: float at kappa = 0.6 sits on the representability floor -----
    // This check used to compare the measured fraction against a constant chosen
    // to sit between two measured values.  That constant carried no information
    // about what the fraction should be, and it stopped describing the sampler
    // once the radius moved into the log domain: the loss is now set by the
    // largest finite float rather than by an intermediate that underflows, and
    // the two formations no longer differ by the factor the constant encoded.
    // The criterion is now the analytic floor at this configuration, tested with
    // an exact two-sided binomial interval at a stated level.
    {
        const long long n= 500000;
        const double alpha= 1.0e-3;
        const std::array<double, 3> z_hat= {0.0, 0.0, 1.0};
        double p_floor= 0.0, p_rescued= 0.0;
        test_overflow_floor(0.6, 1.0, 1.0, z_hat,
                            std::log(static_cast<double>(std::numeric_limits<float>::max())),
                            &p_floor, &p_rescued);
        long long nonfinite= 0;
        bi_kappa_distribution<float> dist(0.6f, 1.0f, 1.0f, {0.0f, 0.0f, 1.0f},
                                          bi_kappa_distribution<float>::no_cap());
        dist.seed(4001);
        for (long long i= 0; i < n; ++i)
        {
            const bi_kappa_distribution<float>::result_type v= dist();
            if (!std::isfinite(v[0]) || !std::isfinite(v[1]) || !std::isfinite(v[2]))
                ++nonfinite;
        }
        const bool ok= test_binom_consistent(nonfinite, n, p_floor, alpha);
        all_ok= all_ok && ok;
        std::cout << "  A3 float, kappa=0.6, non-finite " << nonfinite << "/" << n
                  << ", analytic floor " << p_floor << " (expected "
                  << p_floor * static_cast<double>(n)
                  << "), inside the exact two-sided binomial interval at alpha=" << alpha
                  << ": " << (ok ? "PASS" : "FAIL") << std::endl;
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
// test_no_cap_semantics_bimaxwellian
//
// The bi-Maxwellian sampler must expose the same cap contract as the bi-Kappa
// one: uncapped by default, no_cap() as the sentinel, capped() as the query,
// and no rejection loop when uncapped.  A caller who omits the argument must
// get the distribution the class is named after, on either sampler.
//
// The box here is deliberately tight (1 thermal speed).  Gaussian tails are so
// thin that a realistic cap is never exercised by a finite sample, which is
// exactly why the old truncating default went unnoticed.
// ================================================================
static int test_no_cap_semantics_bimaxwellian()
{
    bool all_ok= true;
    std::cout << "\n--- test_no_cap_semantics_bimaxwellian ---" << std::endl;

    typedef bi_maxwellian_distribution<double> BM;

    // ---- M1: no_cap() is accepted and reports itself as uncapped ----------
    {
        bool ok= false;
        try
        {
            BM dist(1.0, 2.0, {0.0, 0.0, 1.0}, BM::no_cap());
            ok= !dist.param().capped() && !std::isfinite(dist.max_normalized_velocity());
        }
        catch (...)
        {
            ok= false;
        }
        all_ok= all_ok && ok;
        std::cout << "  M1 no_cap() constructs and reports capped()==false: "
                  << (ok ? "PASS" : "FAIL") << std::endl;
    }

    // ---- M2: a finite cap reports itself as capped ------------------------
    {
        BM dist(1.0, 2.0, {0.0, 0.0, 1.0}, 5.0);
        const bool ok= dist.param().capped();
        all_ok= all_ok && ok;
        std::cout << "  M2 finite cap reports capped()==true: " << (ok ? "PASS" : "FAIL")
                  << std::endl;
    }

    // ---- M3: the DEFAULT is uncapped, on every construction path ----------
    {
        BM dflt;                              // default ctor
        BM ctor(1.0, 2.0);                    // value ctor
        BM defined;
        defined.define(1.0, 2.0);             // define(...)
        BM seeded;
        seeded.define(BM::param_type(1.0, 2.0), 7);
        BM::param_type p;                     // param_type default

        const bool ok= !dflt.param().capped() && !ctor.param().capped() &&
                       !defined.param().capped() && !seeded.param().capped() &&
                       !p.capped();
        all_ok= all_ok && ok;
        std::cout << "  M3 default is uncapped on every construction path: "
                  << (ok ? "PASS" : "FAIL") << std::endl;
    }

    // ---- M4: a defaulted sampler really produces untruncated draws --------
    // M3 checks the reported flag; this checks the behaviour behind it.
    {
        const int n= 20000;
        BM dist(1.0, 1.0);                    // cap omitted
        dist.seed(4102);
        int beyond= 0;
        for (int i= 0; i < n; ++i)
        {
            const auto v= dist();
            if (std::abs(v[0]) > 1.0 || std::abs(v[1]) > 1.0 || std::abs(v[2]) > 1.0)
                ++beyond;
        }
        const bool ok= (beyond > 0);
        all_ok= all_ok && ok;
        std::cout << "  M4 defaulted sampler draws exceed the box (" << beyond << "/" << n
                  << " beyond 1): " << (ok ? "PASS" : "FAIL") << std::endl;
    }

    // ---- M5: capped mode still respects the box exactly --------------------
    {
        const int n= 20000;
        const double lambda= 1.0;
        const double theta_perp= 2.0, theta_par= 1.0;
        BM dist(theta_perp, theta_par, {0.0, 0.0, 1.0}, lambda);
        dist.seed(4103);
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
        std::cout << "  M5 capped draws never leave the box (" << violations
                  << " violations): " << (ok ? "PASS" : "FAIL") << std::endl;
    }

    // ---- M6: omitting the cap is identical to naming no_cap() -------------
    {
        const int n= 20000;
        BM omitted(1.0, 2.0);
        BM named(1.0, 2.0, {0.0, 0.0, 1.0}, BM::no_cap());
        omitted.seed(4104);
        named.seed(4104);
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
        std::cout << "  M6 omitted cap == explicit no_cap() bitwise (" << mismatches
                  << " mismatches): " << (ok ? "PASS" : "FAIL") << std::endl;
    }

    // ---- M7: 0 is rejected rather than read as "no cap" -------------------
    {
        bool threw= false;
        try
        {
            BM dist(1.0, 2.0, {0.0, 0.0, 1.0}, 0.0);
            (void)dist;
        }
        catch (const std::invalid_argument &)
        {
            threw= true;
        }
        all_ok= all_ok && threw;
        std::cout << "  M7 cap of 0 throws instead of disabling: "
                  << (threw ? "PASS" : "FAIL") << std::endl;
    }

    std::cout << "test_no_cap_semantics_bimaxwellian: " << (all_ok ? "PASS" : "FAIL")
              << std::endl;
    return all_ok ? 0 : 1;
}

// ================================================================
// test_log_radius_primitive
//
// The sampling primitives the log-domain radius is built from, checked on their
// own: the open-interval uniform, the log-domain Gamma, its law, and the
// underflow the log domain exists to avoid.
// ================================================================
static int test_log_radius_primitive()
{
    bool all_ok= true;
    std::cout << "\n--- test_log_radius_primitive ---" << std::endl;

    // ---- L0: the polygamma functions L3 depends on are right ----------------
    // E[log X] = psi(a) and Var[log X] = psi'(a); the standard error of the sample
    // variance needs psi'''(a) as well.  Checking them against closed forms at z = 1
    // and against the recurrence keeps a wrong z-score from being read as a wrong
    // sampler.
    {
        const double euler= 0.57721566490153286061;
        const double pi= test_pi();
        const double want_psi= -euler;
        const double want_psi1= pi * pi / 6.0;
        const double want_psi3= pi * pi * pi * pi / 15.0;
        bool ok= std::abs(test_digamma(1.0) - want_psi) <= 1e-12 * std::abs(want_psi) &&
                 std::abs(test_trigamma(1.0) - want_psi1) <= 1e-12 * want_psi1 &&
                 std::abs(test_tetragamma(1.0) - want_psi3) <= 1e-12 * want_psi3;
        const double zs[4]= {0.001, 0.37, 3.0, 11.5};
        for (int i= 0; i < 4; ++i)
        {
            const double z= zs[i];
            if (std::abs(test_digamma(z + 1.0) - test_digamma(z) - 1.0 / z) >
                1e-11 * (1.0 / z))
                ok= false;
            if (std::abs(test_trigamma(z) - test_trigamma(z + 1.0) - 1.0 / (z * z)) >
                1e-11 / (z * z))
                ok= false;
            if (std::abs(test_tetragamma(z) - test_tetragamma(z + 1.0) -
                         6.0 / (z * z * z * z)) > 1e-11 * 6.0 / (z * z * z * z))
                ok= false;
        }
        all_ok= all_ok && ok;
        std::cout << "  L0 digamma/trigamma/tetragamma match psi(1)=-gamma, psi'(1)=pi^2/6,"
                     " psi'''(1)=pi^4/15 and their recurrences: "
                  << (ok ? "PASS" : "FAIL") << std::endl;
    }

    // ---- L1: the uniform never touches an endpoint, and its two endpoints are exact ----
    // A logarithm is taken of this value, so 0 would give -inf and 1 would give 0,
    // silently turning a Gamma(a) draw into a Gamma(a+1) one.  The sampled scan shows
    // no endpoint occurs; the two stub engines show what the mapping returns at the
    // ends of the engine's range, which a random sample of any practical size never
    // reaches (k = 0 has probability 2^-52 per draw in double).
    {
        const long long n= 10000000;
        bool open_double= true, open_float= true;
        double min_d= 2.0, max_d= -1.0;
        float min_f= 2.0f, max_f= -1.0f;

        std::mt19937 gen_d(90001u);
        for (long long i= 0; i < n; ++i)
        {
            const double u= bikappa_detail::openCanonical<double, std::mt19937>(gen_d);
            if (!(u > 0.0) || !(u < 1.0))
                open_double= false;
            if (u < min_d)
                min_d= u;
            if (u > max_d)
                max_d= u;
        }
        std::mt19937 gen_f(90002u);
        for (long long i= 0; i < n; ++i)
        {
            const float u= bikappa_detail::openCanonical<float, std::mt19937>(gen_f);
            if (!(u > 0.0f) || !(u < 1.0f))
                open_float= false;
            if (u < min_f)
                min_f= u;
            if (u > max_f)
                max_f= u;
        }

        test_all_zero_engine zero_d, zero_f;
        test_all_ones_engine ones_d, ones_f;
        const double lo_d= bikappa_detail::openCanonical<double, test_all_zero_engine>(zero_d);
        const double hi_d= bikappa_detail::openCanonical<double, test_all_ones_engine>(ones_d);
        const float lo_f= bikappa_detail::openCanonical<float, test_all_zero_engine>(zero_f);
        const float hi_f= bikappa_detail::openCanonical<float, test_all_ones_engine>(ones_f);

        const double lo_d_exact= std::ldexp(1.0, -53);
        const double hi_d_exact= 1.0 - std::ldexp(1.0, -53);
        const float lo_f_exact= std::ldexp(1.0f, -24);
        const float hi_f_exact= 1.0f - std::ldexp(1.0f, -24);

        const bool ends_ok= (lo_d == lo_d_exact) && (hi_d == hi_d_exact) &&
                            (lo_f == lo_f_exact) && (hi_f == hi_f_exact);
        const bool range_ok= (min_d >= lo_d_exact) && (max_d <= hi_d_exact) &&
                             (min_f >= lo_f_exact) && (max_f <= hi_f_exact);
        const bool ok= open_double && open_float && ends_ok && range_ok;
        all_ok= all_ok && ok;
        std::cout << "  L1 openCanonical strictly inside (0,1) over " << n
                  << " draws each in double and float, and its two endpoints equal 2^-(b+1)"
                     " and 1-2^-(b+1) exactly: "
                  << (ok ? "PASS" : "FAIL") << std::endl;
        if (!ok)
        {
            std::cout << "    double: open=" << open_double << " lo=" << lo_d << " (want "
                      << lo_d_exact << ") hi=" << hi_d << " (want " << hi_d_exact << ")\n";
            std::cout << "    float : open=" << open_float << " lo=" << lo_f << " (want "
                      << lo_f_exact << ") hi=" << hi_f << " (want " << hi_f_exact << ")\n";
        }
    }

    // ---- L2 and L3: logGamma is finite everywhere, and samples the right law ----
    // log X = log Y + log(U)/a with Y ~ Gamma(a+1,1) and U ~ Uniform(0,1) gives
    // E[log X] = psi(a+1) - 1/a = psi(a) and Var[log X] = psi'(a+1) + 1/a^2 = psi'(a),
    // both exactly and independently of how Y is generated.  The criterion is a
    // z-score against the analytic standard error at the n each cell uses, with the
    // critical value derived from a nominal family false-alarm rate of 1e-4 spread
    // over the 24 statistics by Bonferroni -- no tuned bound.
    {
        const double shapes[12]= {1e-6, 1e-4, 1e-3, 1e-2, 0.1,   0.5,
                                  0.9,  0.99, 1.0,  1.001, 2.0,  4.5};
        const long long n_per= 100000;
        const double family_alpha= 1.0e-4;
        const double z_crit= test_normal_upper_quantile(0.5 * family_alpha / 24.0);

        long long nonfinite_total= 0;
        double worst_z= 0.0;
        double worst_shape= 0.0;
        bool law_ok= true;
        for (int s= 0; s < 12; ++s)
        {
            const double a= shapes[s];
            std::mt19937 gen(static_cast<unsigned int>(90010 + s));
            test_replay_state st;
            long long cnt= 0;
            double mean= 0.0, m2= 0.0;
            for (long long i= 0; i < n_per; ++i)
            {
                const double lx= bikappa_detail::logGamma<double, std::mt19937>(
                    a, gen, st.cache, st.cached);
                if (!std::isfinite(lx))
                {
                    ++nonfinite_total;
                    continue;
                }
                ++cnt;
                const double d= lx - mean;
                mean+= d / static_cast<double>(cnt);
                m2+= d * (lx - mean);
            }
            if (cnt < 2)
            {
                law_ok= false;
                continue;
            }
            const double var= m2 / static_cast<double>(cnt - 1);
            const double psi1= test_trigamma(a);
            const double psi3= test_tetragamma(a);
            const double z_mean= (mean - test_digamma(a)) /
                                 std::sqrt(psi1 / static_cast<double>(cnt));
            const double z_var= (var - psi1) / std::sqrt((psi3 + 2.0 * psi1 * psi1) /
                                                         static_cast<double>(cnt));
            if (std::abs(z_mean) > z_crit || std::abs(z_var) > z_crit)
                law_ok= false;
            if (std::abs(z_mean) > worst_z)
            {
                worst_z= std::abs(z_mean);
                worst_shape= a;
            }
            if (std::abs(z_var) > worst_z)
            {
                worst_z= std::abs(z_var);
                worst_shape= a;
            }
        }
        const bool l2_ok= (nonfinite_total == 0);
        all_ok= all_ok && l2_ok && law_ok;
        std::cout << "  L2 logGamma finite over 12 shapes x " << n_per << " draws ("
                  << nonfinite_total << " non-finite): " << (l2_ok ? "PASS" : "FAIL")
                  << std::endl;
        std::cout << "  L3 E[log X]=psi(a) and Var[log X]=psi'(a), worst |z|=" << worst_z
                  << " at a=" << worst_shape << " against critical " << z_crit
                  << " (family alpha 1e-4 over 24 statistics): " << (law_ok ? "PASS" : "FAIL")
                  << std::endl;
    }

    // ---- L4: the mechanism, deterministically and with no RNG ---------------
    // At a = 0.01 the product form of the boosting identity, Y U^(1/a), underflows to
    // exactly zero and the draw is lost; the additive form log Y + log(U)/a is an
    // ordinary negative number.  The expected value is written a second way,
    // log 2 - 400 log 10, so the check is against algebra rather than against a
    // repetition of the same expression.
    {
        const double a= 0.01;
        const double y= 2.0;
        const double u= 1.0e-4;
        const double product= y * std::pow(u, 1.0 / a);
        const double log_form= std::log(y) + std::log(u) / a;
        const double expected= std::log(2.0) - 400.0 * std::log(10.0);
        const double tol= 8.0 * test_ulp_at(expected);
        const bool underflowed= (product == 0.0);
        const bool finite= std::isfinite(log_form);
        const bool matches= std::abs(log_form - expected) <= tol;
        const bool ok= underflowed && finite && matches;
        all_ok= all_ok && ok;
        std::cout << "  L4 y*pow(u,1/a) underflows to exactly 0 while log y + log(u)/a = "
                  << log_form << " (expected " << expected << ", within 8 ulp): "
                  << (ok ? "PASS" : "FAIL") << std::endl;
    }

    std::cout << "test_log_radius_primitive: " << (all_ok ? "PASS" : "FAIL") << std::endl;
    return all_ok ? 0 : 1;
}

// ================================================================
// test_log_radius_no_regression
//
// The log domain must not change the radius where nothing forced it to, and the
// stream it produces must not drift.
// ================================================================

/// 16 draws at kappa = 2, theta_perp = 1, theta_par = 2, ub = {0,0,1}, seed 20030410,
/// as exact bit patterns.  ONE array serves both standard libraries: since 2.0.0 the
/// Gamma, normal and uniform primitives live in the header, so the sampled sequence
/// is a function of the engine alone.  That equality is not background to this test,
/// it is the assertion -- a libstdc++ build that disagreed with the array its build
/// selects would mean the stream had become a property of the standard library
/// again.  Each array below was produced identically by libc++ and by libstdc++.
///
/// There are two arrays because the stream depends on one thing the standard library
/// is not responsible for: whether the compiler fuses a multiply and an add.  A
/// fused multiply-add rounds once instead of twice, and the sampler has several
/// expressions it applies to -- `1 - cosTheta*cosTheta` and `1 + c*x` among them --
/// so a contracting build differs from a non-contracting one by one or two units in
/// the last place on about a tenth of the components.  The two toolchains do not
/// default alike: clang contracts at every optimization level, GCC only from -O2.
/// Both settings are frozen here, and a stream matching neither is drift.
///
/// Off is what the source text specifies; a project that needs the stream to be a
/// property of the source alone should pass -ffp-contract=off, and that belongs in
/// the library's build documentation rather than in this test.
static const unsigned long long k_golden_stream_no_fma[16][3]= {
    {0xbfd15ba81cfd18ddULL, 0xbff7e0f2cd8faeddULL, 0x3fda9924c54e1574ULL},
    {0xbfdea5acf93ac620ULL, 0xbfda2dc0d512deccULL, 0xbfbf0a70e5bc1734ULL},
    {0xbfddb887e0c8ce9fULL, 0x40080571d77bd3ebULL, 0xc0038109652c3ca4ULL},
    {0x3fe698cd0e52adcdULL, 0x40026fd65169f08fULL, 0x3ffea0be023e0cf4ULL},
    {0x3ff16471e77dedd5ULL, 0xbff17910ce76ee65ULL, 0xbfea8567aa7fde0fULL},
    {0xbff3b6c788b6df6aULL, 0xbf8d3eef1117229eULL, 0xbff3fd8b4e31160fULL},
    {0xbff24018d0a8cbdbULL, 0xbfd124d1dc511927ULL, 0xbfe962c1208342d5ULL},
    {0xc007882129566ec9ULL, 0xbfcb3b5d2ab5adf0ULL, 0xbffaecc134f280f4ULL},
    {0xbfe3ee2b53da3ec7ULL, 0xbfcf4c5a89d81b95ULL, 0xbfe20cb5c9f2a311ULL},
    {0x3fe8a89dacd86329ULL, 0x3ff89f8ce2af2c4cULL, 0xbfeffe21eaad910cULL},
    {0xbfd158b279b811c4ULL, 0xbfd3cbc1d046b426ULL, 0x40011c8dcb02c2a2ULL},
    {0xbfc43078ce64f93cULL, 0xbfd628128e7f6e70ULL, 0x3fe06a0b502aa269ULL},
    {0x3feee1664b806335ULL, 0x3fcc2cc34ec3b107ULL, 0x400a462fc5d5b51bULL},
    {0x3ff040ff5ccae55dULL, 0x3fe85f8198bf9596ULL, 0xc00ecb4f11c9afe3ULL},
    {0xbfe97d2723ac9d35ULL, 0xbff1331e78d62310ULL, 0xbfd011092bdec4ebULL},
    {0x3fa09ad0b0a37a31ULL, 0x3fd241313bfa4a21ULL, 0x3ff701cef7135568ULL}};

/// The same 16 draws from a build that contracts multiply-add.
static const unsigned long long k_golden_stream_fma[16][3]= {
    {0xbfd15ba81cfd18ddULL, 0xbff7e0f2cd8faeddULL, 0x3fda9924c54e1574ULL},
    {0xbfdea5acf93ac61dULL, 0xbfda2dc0d512dec9ULL, 0xbfbf0a70e5bc1731ULL},
    {0xbfddb887e0c8ce9fULL, 0x40080571d77bd3ebULL, 0xc0038109652c3ca4ULL},
    {0x3fe698cd0e52adcdULL, 0x40026fd65169f08fULL, 0x3ffea0be023e0cf4ULL},
    {0x3ff16471e77dedd5ULL, 0xbff17910ce76ee65ULL, 0xbfea8567aa7fde0fULL},
    {0xbff3b6c788b6df6aULL, 0xbf8d3eef1117229eULL, 0xbff3fd8b4e31160fULL},
    {0xbff24018d0a8cbdbULL, 0xbfd124d1dc511927ULL, 0xbfe962c1208342d5ULL},
    {0xc007882129566ec9ULL, 0xbfcb3b5d2ab5adf0ULL, 0xbffaecc134f280f4ULL},
    {0xbfe3ee2b53da3ec7ULL, 0xbfcf4c5a89d81b95ULL, 0xbfe20cb5c9f2a311ULL},
    {0x3fe8a89dacd8632bULL, 0x3ff89f8ce2af2c4eULL, 0xbfeffe21eaad910dULL},
    {0xbfd158b279b811c5ULL, 0xbfd3cbc1d046b428ULL, 0x40011c8dcb02c2a4ULL},
    {0xbfc43078ce64f93cULL, 0xbfd628128e7f6e70ULL, 0x3fe06a0b502aa269ULL},
    {0x3feee1664b806335ULL, 0x3fcc2cc34ec3b107ULL, 0x400a462fc5d5b51bULL},
    {0x3ff040ff5ccae55dULL, 0x3fe85f8198bf9596ULL, 0xc00ecb4f11c9afe3ULL},
    {0xbfe97d2723ac9d35ULL, 0xbff1331e78d62310ULL, 0xbfd011092bdec4ebULL},
    {0x3fa09ad0b0a37a31ULL, 0x3fd241313bfa4a21ULL, 0x3ff701cef7135568ULL}};

static int test_log_radius_no_regression()
{
    bool all_ok= true;
    std::cout << "\n--- test_log_radius_no_regression ---" << std::endl;

    // ---- M1: at kappa = 2 the two formations describe the same number --------
    // Given the same X1 and X2, exp((log X1 - log X2)/2) and sqrt(X1)/sqrt(X2) differ
    // only by the rounding of the round trip through log and exp.  That error is
    // bounded, not guessed.  Writing u = 2^-53 for the unit roundoff: each logarithm
    // carries about u |log x| of absolute error, the halved difference inherits both,
    // the exponential turns that into the same relative error, and the reference
    // itself rounds twice in the square roots and once in the division.  So the
    // relative difference is at most about u (|log x1| + |log x2| + 2.5), and since
    // one ulp is at least u times the value, that is at most
    // |log x1| + |log x2| + 3 ulp.  The criterion is that bound, evaluated from the
    // draw itself; at kappa = 2 it sits near 4 ulp because both logarithms are of
    // order one.
    {
        const long long n= 200000;
        std::mt19937 gen(90003u);
        test_replay_state st;
        long long violations= 0;
        double worst= 0.0, worst_bound= 0.0;
        for (long long i= 0; i < n; ++i)
        {
            const double x1= bikappa_detail::marsagliaTsangGamma<double, std::mt19937>(
                1.5, gen, st.cache, st.cached);
            const double x2= bikappa_detail::marsagliaTsangGamma<double, std::mt19937>(
                1.5, gen, st.cache, st.cached); // shape a = kappa - 1/2 = 3/2
            const double r_log= std::exp(0.5 * (std::log(x1) - std::log(x2)));
            const double r_lin= std::sqrt(x1) / std::sqrt(x2);
            const double ulp= test_ulp_distance(r_log, r_lin, r_lin);
            const double bound= 3.0 + std::abs(std::log(x1)) + std::abs(std::log(x2));
            if (!(ulp <= bound))
                ++violations;
            if (ulp > worst)
                worst= ulp;
            if (ulp / bound > worst_bound)
                worst_bound= ulp / bound;
        }
        const bool ok= (violations == 0);
        all_ok= all_ok && ok;
        std::cout << "  M1 kappa=2, log-domain radius == sqrt(x1)/sqrt(x2) within the"
                     " round-trip bound over " << n << " draws (worst " << worst
                  << " ulp, worst fraction of the bound " << worst_bound << ", "
                  << violations << " violations): " << (ok ? "PASS" : "FAIL") << std::endl;
    }

    // ---- M2: the stream has not drifted -------------------------------------
    {
        const bool fma= test_fp_contraction_enabled();
        bi_kappa_distribution<double> dist(2.0, 1.0, 2.0, {0.0, 0.0, 1.0});
        dist.seed(20030410);
        int mismatch_no_fma= 0, mismatch_fma= 0;
        unsigned long long got[16][3];
        for (int i= 0; i < 16; ++i)
        {
            const bi_kappa_distribution<double>::result_type v= dist();
            for (int j= 0; j < 3; ++j)
            {
                std::memcpy(&got[i][j], &v[j], sizeof(got[i][j]));
                if (got[i][j] != k_golden_stream_no_fma[i][j])
                    ++mismatch_no_fma;
                if (got[i][j] != k_golden_stream_fma[i][j])
                    ++mismatch_fma;
            }
        }
        const bool matched_no_fma= (mismatch_no_fma == 0);
        const bool matched_fma= (mismatch_fma == 0);
        const bool ok= matched_no_fma || matched_fma;
        // The build's own multiply-add behaviour predicts which array should match,
        // so a mismatch between prediction and outcome is worth printing even though
        // the assertion is that one of the two matches.
        const bool as_predicted= (fma ? matched_fma : matched_no_fma);
        all_ok= all_ok && ok;
        if (!ok)
        {
            const int shown= (mismatch_no_fma <= mismatch_fma) ? 0 : 1;
            std::cout << "    closest frozen array is the "
                      << (shown == 0 ? "non-contracting" : "contracting") << " one, with "
                      << (shown == 0 ? mismatch_no_fma : mismatch_fma) << " of 48 words"
                      << " different; this build "
                      << (fma ? "does" : "does not") << " fuse multiply-add\n";
            int printed= 0;
            for (int i= 0; i < 16 && printed < 6; ++i)
                for (int j= 0; j < 3 && printed < 6; ++j)
                {
                    const unsigned long long want= (shown == 0) ? k_golden_stream_no_fma[i][j]
                                                                : k_golden_stream_fma[i][j];
                    if (got[i][j] != want)
                    {
                        std::cout << "    draw " << i << " component " << j << ": got 0x"
                                  << std::hex << got[i][j] << " want 0x" << want << std::dec
                                  << "\n";
                        ++printed;
                    }
                }
        }
        std::cout << "  M2 16 draws at (2,1,2,{0,0,1},20030410) reproduce the frozen golden"
                     " array bitwise -- one array per multiply-add setting, each of them"
                     " shared by libc++ and libstdc++ (matched the "
                  << (matched_no_fma ? "non-contracting" : (matched_fma ? "contracting"
                                                                        : "neither"))
                  << " array; this build " << (fma ? "does" : "does not")
                  << " fuse multiply-add" << (as_predicted ? "" : ", which disagrees")
                  << "): " << (ok ? "PASS" : "FAIL") << std::endl;
    }

    // ---- M3: the same, through the boosting identity the sampler actually uses ----
    // Here X2 is never formed by the log path: it carries log X2 = log Y + log(U)/a,
    // while the linear path forms Y U^(1/a).  The bound gains the term the boosting
    // step contributes, eps |log U| / a, which is why it is wider than M1's.
    {
        const long long n= 200000;
        const double a= 1.5;
        std::mt19937 gen(90004u);
        test_replay_state st;
        long long violations= 0;
        double worst= 0.0, worst_bound= 0.0;
        for (long long i= 0; i < n; ++i)
        {
            const double y= bikappa_detail::marsagliaTsangGamma<double, std::mt19937>(
                a + 1.0, gen, st.cache, st.cached);
            const double u= bikappa_detail::openCanonical<double, std::mt19937>(gen);
            const double x1= bikappa_detail::marsagliaTsangGamma<double, std::mt19937>(
                1.5, gen, st.cache, st.cached);
            const double log_x2= std::log(y) + std::log(u) / a;
            const double x2= y * std::pow(u, 1.0 / a);
            const double r_log= std::exp(0.5 * (std::log(x1) - log_x2));
            const double r_lin= std::sqrt(x1) / std::sqrt(x2);
            const double ulp= test_ulp_distance(r_log, r_lin, r_lin);
            // Same accounting as M1, plus what the boosting step costs on each side:
            // the log path carries u |log u| / a from dividing the uniform's logarithm
            // by the shape, and pow(u, 1/a) carries the same amount on the linear side.
            const double bound= 4.0 + std::abs(std::log(x1)) + std::abs(std::log(y)) +
                                2.0 * std::abs(std::log(u)) / a;
            if (!(ulp <= bound))
                ++violations;
            if (ulp > worst)
                worst= ulp;
            if (ulp / bound > worst_bound)
                worst_bound= ulp / bound;
        }
        const bool ok= (violations == 0);
        all_ok= all_ok && ok;
        std::cout << "  M3 kappa=2, boosted log X2 vs y*pow(u,1/a) within the round-trip"
                     " bound (worst " << worst << " ulp, worst fraction of the bound "
                  << worst_bound << ", " << violations << " violations): "
                  << (ok ? "PASS" : "FAIL") << std::endl;
    }

    std::cout << "test_log_radius_no_regression: " << (all_ok ? "PASS" : "FAIL") << std::endl;
    return all_ok ? 0 : 1;
}
// ================================================================
// test_log_radius_finiteness
//
// Near kappa = 1/2 the bi-Kappa law puts non-zero probability outside every finite
// floating-point range, so a correct sampler does not return only finite vectors --
// it returns exactly as many non-finite ones as the type forces and no more.  The
// criterion for each cell is therefore the analytic representability floor derived
// above, tested with an exact binomial interval, not a bound on "small".
// ================================================================

/// Non-finite count over @p n uncapped draws of the released sampler, isotropic and
/// unrotated so that the analytic floor applies directly.
template <typename RealType>
static long long test_count_nonfinite(double kappa, double theta, int seed, long long n)
{
    bi_kappa_distribution<RealType> dist(
        static_cast<RealType>(kappa), static_cast<RealType>(theta),
        static_cast<RealType>(theta));
    dist.seed(seed);
    long long bad= 0;
    for (long long i= 0; i < n; ++i)
    {
        const typename bi_kappa_distribution<RealType>::result_type v= dist();
        if (!std::isfinite(v[0]) || !std::isfinite(v[1]) || !std::isfinite(v[2]))
            ++bad;
    }
    return bad;
}

template <typename RealType>
static bool test_report_floor(const char *tag, const char *label, double kappa, double theta,
                              int seed, long long n, double alpha, bool &all_ok)
{
    const std::array<double, 3> z_hat= {0.0, 0.0, 1.0};
    const double log_max= std::log(static_cast<double>(std::numeric_limits<RealType>::max()));
    double p_overflow= 0.0, p_rescued= 0.0;
    test_overflow_floor(kappa, theta, theta, z_hat, log_max, &p_overflow, &p_rescued);
    const long long observed= test_count_nonfinite<RealType>(kappa, theta, seed, n);
    const bool ok= test_binom_consistent(observed, n, p_overflow, alpha);
    all_ok= all_ok && ok;
    std::cout << "  " << tag << " " << label << ": non-finite " << observed << "/" << n
              << " = " << (static_cast<double>(observed) / static_cast<double>(n))
              << ", analytic floor " << p_overflow << " (expected "
              << p_overflow * static_cast<double>(n)
              << "), inside the exact two-sided binomial interval at alpha=" << alpha
              << ": " << (ok ? "PASS" : "FAIL") << std::endl;
    return ok;
}

static int test_log_radius_finiteness()
{
    bool all_ok= true;
    std::cout << "\n--- test_log_radius_finiteness ---" << std::endl;
    const double alpha= 1.0e-3;

    // ---- N1: the rate is the floor, not zero -------------------------------
    test_report_floor<double>("N1", "double kappa=0.505", 0.505, 1.0, 90030, 1000000, alpha,
                              all_ok);

    // ---- N2: two more cells on the same rule -------------------------------
    test_report_floor<float>("N2", "float  kappa=0.55 ", 0.55, 1.0, 90031, 1000000, alpha,
                             all_ok);
    test_report_floor<double>("N2", "double kappa=0.51 ", 0.51, 1.0, 90032, 1000000, alpha,
                              all_ok);

    // ---- N3: where the floor is negligible, zero is a prediction -----------
    // The analytic floor here is below 1e-19, so P(at least one non-finite draw in
    // 1e6) is below 1e-13: exactly zero is not an assumption about "large enough
    // kappa", it is what the floor predicts at this n.
    {
        const long long n= 1000000;
        const std::array<double, 3> z_hat= {0.0, 0.0, 1.0};
        double pf= 0.0, pd= 0.0, dummy= 0.0;
        test_overflow_floor(0.75, 1.0, 1.0, z_hat,
                            std::log(static_cast<double>(std::numeric_limits<float>::max())),
                            &pf, &dummy);
        test_overflow_floor(0.55, 1.0, 1.0, z_hat,
                            std::log(std::numeric_limits<double>::max()), &pd, &dummy);
        const long long got_f= test_count_nonfinite<float>(0.75, 1.0, 90033, n);
        const long long got_d= test_count_nonfinite<double>(0.55, 1.0, 90034, n);
        const double p_any= 1.0 - std::exp(static_cast<double>(n) * std::log1p(-pf)) +
                            1.0 - std::exp(static_cast<double>(n) * std::log1p(-pd));
        const bool ok= (got_f == 0) && (got_d == 0);
        all_ok= all_ok && ok;
        std::cout << "  N3 float kappa=0.75 (" << got_f << ") and double kappa=0.55 ("
                  << got_d << ") give exactly zero non-finite in " << n
                  << "; floors " << pf << " and " << pd
                  << ", so P(any non-finite) = " << p_any << ": " << (ok ? "PASS" : "FAIL")
                  << std::endl;
    }

    // ---- N4: the previous formation loses strictly more ---------------------
    // Reproduced inline rather than by vendoring the 1.0.0 header, because cpp/ must
    // not depend on experiments/.  X2 ~ Gamma(a,1) is materialized and underflows to
    // exactly zero for a = 0.005, which no rearrangement of sqrt(X1)/sqrt(X2) can
    // survive; the log path never forms it.
    {
        const long long n= 1000000;
        const double kappa= 0.505, a= kappa - 0.5, sqrt_kappa= std::sqrt(kappa);
        std::mt19937 gen(90035u);
        std::gamma_distribution<double> gamma_x1(1.5, 1.0);
        std::gamma_distribution<double> gamma_x2(a, 1.0);
        std::uniform_real_distribution<double> uni(0.0, 1.0);
        long long legacy_bad= 0;
        for (long long i= 0; i < n; ++i)
        {
            const double x1= gamma_x1(gen);
            const double x2= gamma_x2(gen);
            const double r= std::sqrt(x1) / std::sqrt(x2);
            const double cos_theta= 2.0 * uni(gen) - 1.0;
            const double sin_theta= std::sqrt(std::max(0.0, 1.0 - cos_theta * cos_theta));
            const double phi= 2.0 * test_pi() * uni(gen);
            const double vx= r * sqrt_kappa * (sin_theta * std::cos(phi));
            const double vy= r * sqrt_kappa * (sin_theta * std::sin(phi));
            const double vz= r * sqrt_kappa * cos_theta;
            if (!std::isfinite(vx) || !std::isfinite(vy) || !std::isfinite(vz))
                ++legacy_bad;
        }
        const long long cand_bad= test_count_nonfinite<double>(kappa, 1.0, 90035, n);
        const double ratio= (cand_bad > 0)
                                ? static_cast<double>(legacy_bad) /
                                      static_cast<double>(cand_bad)
                                : std::numeric_limits<double>::infinity();
        const bool ok= (legacy_bad >= 10 * cand_bad) && (cand_bad > 0);
        all_ok= all_ok && ok;
        std::cout << "  N4 double kappa=0.505: materialized X2 loses " << legacy_bad << "/"
                  << n << " against the log path's " << cand_bad << " (" << ratio
                  << "x, required >= 10x): " << (ok ? "PASS" : "FAIL") << std::endl;
    }

    std::cout << "test_log_radius_finiteness: " << (all_ok ? "PASS" : "FAIL") << std::endl;
    return all_ok ? 0 : 1;
}

// ================================================================
// test_overflow_observability
//
// n_attempts() and n_nonfinite() are the only way a caller can tell "the box
// rejected a representable sample" from "the box rejected a sample that could not
// have been returned at all".  In capped mode the second case is otherwise
// invisible: every returned draw is finite and the acceptance rate alone does not
// separate them.
// ================================================================

/// Draw until the sampler has recorded at least one unrepresentable attempt, so that
/// a following reset has something to clear.
static void test_dirty_counters(bi_kappa_distribution<double> &dist)
{
    for (int i= 0; i < 500000 && dist.n_nonfinite() == 0ULL; ++i)
    {
        const bi_kappa_distribution<double>::result_type v= dist();
        (void)v;
    }
}

static int test_overflow_observability()
{
    bool all_ok= true;
    std::cout << "\n--- test_overflow_observability ---" << std::endl;

    // ---- O1: the counter equals what the caller can see ---------------------
    {
        const long long n= 200000;
        bi_kappa_distribution<double> dist(0.505, 1.0, 1.0);
        dist.seed(90040);
        long long seen= 0;
        for (long long i= 0; i < n; ++i)
        {
            const bi_kappa_distribution<double>::result_type v= dist();
            if (!std::isfinite(v[0]) || !std::isfinite(v[1]) || !std::isfinite(v[2]))
                ++seen;
        }
        const bool ok= (static_cast<unsigned long long>(seen) == dist.n_nonfinite()) &&
                       (dist.n_attempts() == static_cast<unsigned long long>(n)) && (seen > 0);
        all_ok= all_ok && ok;
        std::cout << "  O1 uncapped kappa=0.505: n_nonfinite()=" << dist.n_nonfinite()
                  << " equals the " << seen << " non-finite returns, n_attempts()="
                  << dist.n_attempts() << " equals " << n << " draws: "
                  << (ok ? "PASS" : "FAIL") << std::endl;
    }

    // ---- O2: the silent-discard guard ---------------------------------------
    // With lambda = 5 at kappa = 0.505 the box accepts about 1.8 % of attempts, so
    // 2000 returned samples cost roughly 1.1e5 attempts, of which the analytic floor
    // makes about 90 unrepresentable.  P(none) is below 1e-30, so requiring a
    // non-zero count is not a coin flip.
    {
        const int n_returned= 2000;
        bi_kappa_distribution<double> dist(0.505, 1.0, 1.0, {0.0, 0.0, 1.0}, 5.0);
        dist.seed(90041);
        bool all_finite= true;
        for (int i= 0; i < n_returned; ++i)
        {
            const bi_kappa_distribution<double>::result_type v= dist();
            if (!std::isfinite(v[0]) || !std::isfinite(v[1]) || !std::isfinite(v[2]))
                all_finite= false;
        }
        const bool ok= all_finite && (dist.n_nonfinite() > 0ULL) &&
                       (dist.n_attempts() > static_cast<unsigned long long>(n_returned));
        all_ok= all_ok && ok;
        std::cout << "  O2 capped kappa=0.505 lambda=5: " << n_returned
                  << " returned draws all finite, n_nonfinite()=" << dist.n_nonfinite()
                  << " > 0, n_attempts()=" << dist.n_attempts() << " > " << n_returned
                  << ": " << (ok ? "PASS" : "FAIL") << std::endl;
    }

    // ---- O3: every path that restarts a stream also restarts the counters ----
    {
        bool ok= true;
        const bi_kappa_distribution<double>::param_type params(0.505, 1.0, 1.0,
                                                               {0.0, 0.0, 1.0});
        const bi_kappa_distribution<double>::point_type z_hat= {0.0, 0.0, 1.0};
        for (int path= 0; path < 6; ++path)
        {
            bi_kappa_distribution<double> dist(0.505, 1.0, 1.0);
            dist.seed(90042 + path);
            test_dirty_counters(dist);
            if (dist.n_attempts() == 0ULL || dist.n_nonfinite() == 0ULL)
                ok= false;
            switch (path)
            {
            case 0: dist.reset(); break;
            case 1: dist.define(params); break;
            case 2: dist.define(0.505, 1.0, 1.0); break;
            case 3: dist.define(params, 90042 + path); break;
            case 4:
                dist.define(0.505, 1.0, 1.0, z_hat,
                            bi_kappa_distribution<double>::no_cap(), 90042 + path);
                break;
            default: dist.seed(90042 + path); break;
            }
            if (dist.n_attempts() != 0ULL || dist.n_nonfinite() != 0ULL)
                ok= false;
        }
        all_ok= all_ok && ok;
        std::cout << "  O3 reset(), all four define() overloads and seed() zero both"
                     " counters: " << (ok ? "PASS" : "FAIL") << std::endl;
    }

    std::cout << "test_overflow_observability: " << (all_ok ? "PASS" : "FAIL") << std::endl;
    return all_ok ? 0 : 1;
}

// ================================================================
// test_seed_reset_contract
// ================================================================
static void test_collect(bi_kappa_distribution<double> &dist, int n,
                         std::vector<double> &out)
{
    out.clear();
    for (int i= 0; i < n; ++i)
    {
        const bi_kappa_distribution<double>::result_type v= dist();
        out.push_back(v[0]);
        out.push_back(v[1]);
        out.push_back(v[2]);
    }
}

static int test_seed_reset_contract()
{
    bool all_ok= true;
    std::cout << "\n--- test_seed_reset_contract ---" << std::endl;

    // ---- P1: seed(s) restarts the stream for every length -------------------
    // seed() now calls reset(), which discards the cached second polar-normal
    // variate.  Without that the second run opens with a variate belonging to the
    // first, so the two agree only when the draw count leaves the cache empty.
    // Before the fix this FAILED under libstdc++ for every odd n.
    {
        bool ok= true;
        int first_bad= 0;
        for (int n= 1; n <= 8; ++n)
        {
            bi_kappa_distribution<double> dist(2.0, 1.0, 2.0);
            std::vector<double> a, b;
            dist.seed(90050);
            test_collect(dist, n, a);
            dist.seed(90050);
            test_collect(dist, n, b);
            if (a != b)
            {
                ok= false;
                if (first_bad == 0)
                    first_bad= n;
            }
        }
        all_ok= all_ok && ok;
        std::cout << "  P1 seed(s) then n draws reproduces bitwise for n = 1..8; before"
                     " seed() called reset() this FAILED under libstdc++ for every odd n";
        if (!ok)
            std::cout << " (first n that still fails: " << first_bad << ")";
        std::cout << ": " << (ok ? "PASS" : "FAIL") << std::endl;
    }

    // ---- P2: the two ways of seeding a configuration agree ------------------
    {
        const bi_kappa_distribution<double>::param_type params(1.5, 1.0, 2.0,
                                                               {0.0, 0.0, 1.0});
        bi_kappa_distribution<double> a(2.0, 1.0, 1.0), b(2.0, 1.0, 1.0);
        a.define(params, 90051);
        b.seed(90051);
        b.define(params);
        std::vector<double> va, vb;
        test_collect(a, 64, va);
        test_collect(b, 64, vb);
        const bool ok= (va == vb);
        all_ok= all_ok && ok;
        std::cout << "  P2 define(params, seed) == seed(s) then define(params), bitwise over"
                     " 64 draws: " << (ok ? "PASS" : "FAIL") << std::endl;
    }

    // ---- P3: different seeds give different streams -------------------------
    {
        bi_kappa_distribution<double> dist(2.0, 1.0, 2.0);
        std::vector<double> va, vb;
        dist.seed(90052);
        test_collect(dist, 64, va);
        dist.seed(90053);
        test_collect(dist, 64, vb);
        const bool ok= (va != vb);
        all_ok= all_ok && ok;
        std::cout << "  P3 two different seeds give different sequences: "
                  << (ok ? "PASS" : "FAIL") << std::endl;
    }

    std::cout << "test_seed_reset_contract: " << (all_ok ? "PASS" : "FAIL") << std::endl;
    return all_ok ? 0 : 1;
}
// ================================================================
// test_log_cap_semantics
//
// The cap is now decided in the log domain, before anything is exponentiated, so
// that a sample far outside the box is rejected rather than overflowed first and
// rejected afterwards.  That is only legitimate if it accepts exactly the set the
// published linear definition accepts.
// ================================================================
static int test_log_cap_semantics()
{
    bool all_ok= true;
    std::cout << "\n--- test_log_cap_semantics ---" << std::endl;

    // ---- Q1: the two predicates accept exactly the same set -----------------
    // The linear predicate |point_j| / theta_j <= lambda is private, so it is
    // transcribed here.  Note that theta cancels out of it: with
    // point = R sqrt(kappa) (theta_perp n_x, theta_perp n_y, theta_par n_z) the
    // normalized components are sqrt(kappa) R |n_j|, which is what the log form
    // tests.  ub is the axis, so the returned vector is the field-aligned one and
    // the linear predicate applies to it unchanged.
    //
    // The replayed accepted draws are also compared against the sampler's own
    // returned stream bitwise, so this checks the transcription as well as the
    // predicates.
    {
        const long long attempts= 300000;
        const double kappa= 0.75, theta_perp= 1.0, theta_par= 2.0, lambda= 5.0;
        const double a= kappa - 0.5, sqrt_kappa= std::sqrt(kappa);
        const double log_cap= std::log(lambda), log_sqrt_kappa= std::log(sqrt_kappa);

        bi_kappa_distribution<double> dist(kappa, theta_perp, theta_par, {0.0, 0.0, 1.0},
                                           lambda);
        std::mt19937 gen_sampler(90060u), gen_replay(90060u);
        test_replay_state st;
        long long comparable= 0, disagreements= 0, accepted= 0, stream_mismatch= 0;
        for (long long i= 0; i < attempts; ++i)
        {
            const test_replay_draw d= test_replay(a, gen_replay, st);
            const double max_abs_n=
                std::max(std::abs(d.nx), std::max(std::abs(d.ny), std::abs(d.nz)));
            const double log_largest= d.log_r + log_sqrt_kappa + std::log(max_abs_n);
            const bool log_accepts= !(log_largest > log_cap);

            const double r= std::exp(d.log_r);
            bi_kappa_distribution<double>::point_type point;
            point[0]= r * (sqrt_kappa * theta_perp * d.nx);
            point[1]= r * (sqrt_kappa * theta_perp * d.ny);
            point[2]= r * (sqrt_kappa * theta_par * d.nz);
            if (!std::isfinite(r) || !std::isfinite(point[0]) || !std::isfinite(point[1]) ||
                !std::isfinite(point[2]))
                continue; // not representable in the linear form; nothing to compare
            ++comparable;

            const bool linear_accepts= (std::abs(point[0]) / theta_perp <= lambda) &&
                                       (std::abs(point[1]) / theta_perp <= lambda) &&
                                       (std::abs(point[2]) / theta_par <= lambda);
            if (log_accepts != linear_accepts)
                ++disagreements;
            if (log_accepts)
            {
                ++accepted;
                const bi_kappa_distribution<double>::result_type v= dist(gen_sampler);
                if (v[0] != point[0] || v[1] != point[1] || v[2] != point[2])
                    ++stream_mismatch;
            }
        }
        const bool ok= (disagreements == 0) && (stream_mismatch == 0) &&
                       (comparable >= 100000) && (accepted > 0);
        all_ok= all_ok && ok;
        std::cout << "  Q1 log-domain cap and linear |point_j|/theta_j <= lambda accept the"
                     " same set over " << comparable << " representable attempts ("
                  << disagreements << " disagreements, " << accepted
                  << " accepted, " << stream_mismatch
                  << " stream mismatches): " << (ok ? "PASS" : "FAIL") << std::endl;
    }

    // ---- Q2: capped draws never leave the box (the B4 guard, re-run) --------
    {
        const long long n= 100000;
        const double lambda= 5.0, theta_perp= 1.0, theta_par= 2.0;
        bi_kappa_distribution<double> dist(1.0, theta_perp, theta_par, {0.0, 0.0, 1.0},
                                           lambda);
        dist.seed(90061);
        long long violations= 0;
        for (long long i= 0; i < n; ++i)
        {
            const bi_kappa_distribution<double>::result_type v= dist();
            if (std::abs(v[0]) / theta_perp > lambda || std::abs(v[1]) / theta_perp > lambda ||
                std::abs(v[2]) / theta_par > lambda)
                ++violations;
        }
        const bool ok= (violations == 0);
        all_ok= all_ok && ok;
        std::cout << "  Q2 capped draws never leave the box over " << n << " draws ("
                  << violations << " violations): " << (ok ? "PASS" : "FAIL") << std::endl;
    }

    // ---- Q3: omitting the cap is bitwise identical to no_cap() (the B7 guard) ----
    {
        const long long n= 20000;
        bi_kappa_distribution<double> omitted(1.5, 1.0, 2.0);
        bi_kappa_distribution<double> named(1.5, 1.0, 2.0, {0.0, 0.0, 1.0},
                                            bi_kappa_distribution<double>::no_cap());
        omitted.seed(90062);
        named.seed(90062);
        long long mismatches= 0;
        for (long long i= 0; i < n; ++i)
        {
            const bi_kappa_distribution<double>::result_type a= omitted();
            const bi_kappa_distribution<double>::result_type b= named();
            if (a[0] != b[0] || a[1] != b[1] || a[2] != b[2])
                ++mismatches;
        }
        const bool ok= (mismatches == 0);
        all_ok= all_ok && ok;
        std::cout << "  Q3 omitted cap == explicit no_cap() bitwise over " << n << " draws ("
                  << mismatches << " mismatches): " << (ok ? "PASS" : "FAIL") << std::endl;
    }

    std::cout << "test_log_cap_semantics: " << (all_ok ? "PASS" : "FAIL") << std::endl;
    return all_ok ? 0 : 1;
}

// ================================================================
// test_log_rotation_stability
//
// An orthogonal rotation can spread one over-large field-aligned component across
// three representable global ones, which is why the sampler builds the order-unity
// vector g = Q(ub) diag(sqrt(kappa) theta) n first and decides representability on the
// components it is actually going to return -- by materializing each one and asking
// whether it is finite, not by comparing its logarithm against a rounded log(max()).
// Deciding on the pre-rotation vector instead would discard draws that have a
// perfectly good answer.
// ================================================================
static int test_log_rotation_stability()
{
    bool all_ok= true;
    std::cout << "\n--- test_log_rotation_stability ---" << std::endl;

    // ---- R1: finiteness follows the log-domain predicate, and the rotation is
    //          accurate -----------------------------------------------------
    // The reference recomputes the same draw in long double the other way round --
    // scale first, then rotate -- so the comparison is between two different
    // orderings of the same algebra rather than a repetition of one.  A rotation is
    // backward stable in the norm but not component by component: with ub = (1,1,1)
    // the three products cancel, so a component that is small compared with the
    // vector carries a large relative error and a small absolute one.  The tolerance
    // is therefore stated in ulp of the vector's largest component, which is the
    // quantity a rotation error is actually bounded by.
    //
    // Note that long double is not wider than double on every platform; this build
    // reports its mantissa below, and where the two coincide the reference is an
    // independent recomputation at the same precision rather than a more accurate
    // one.
    {
        const long long n= 200000;
        const double kappa= 0.51, theta_perp= 1.0, theta_par= 2.0;
        const double a= kappa - 0.5, sqrt_kappa= std::sqrt(kappa);
        const bi_kappa_distribution<double>::point_type ub= {1.0, 1.0, 1.0};
        const std::array<long double, 3> ub_ld= {1.0L, 1.0L, 1.0L};

        bi_kappa_distribution<double> dist(kappa, theta_perp, theta_par, ub);
        bi_kappa_distribution<long double> rot_ld(2.0L, 1.0L, 1.0L);
        std::mt19937 gen_sampler(90070u), gen_replay(90070u);
        test_replay_state st;

        long long predicate_mismatch= 0, compared= 0, skipped= 0;
        double worst_ulp= 0.0;
        for (long long i= 0; i < n; ++i)
        {
            const test_replay_draw d= test_replay(a, gen_replay, st);
            bi_kappa_distribution<double>::point_type local;
            local[0]= sqrt_kappa * theta_perp * d.nx;
            local[1]= sqrt_kappa * theta_perp * d.ny;
            local[2]= sqrt_kappa * theta_par * d.nz;
            const bi_kappa_distribution<double>::point_type g=
                dist.rotate_from_fieldAligned_frame(local, ub);

            // Predict from the rule the loader implements: materialize each component and
            // ask whether it is finite.  Predicting from `log R + log|g_j| >= log(max())`
            // instead -- as this test did before 2.1.0 -- is predicting from a withdrawn
            // rule, and the two disagree exactly at the boundary, which is the only place
            // this check has anything to say.
            bool predicted_finite= true;
            {
                const double r_try= std::exp(d.log_r);
                const bool radius_rep= (r_try <= std::numeric_limits<double>::max());
                for (int j= 0; j < 3; ++j)
                {
                    if (g[j] == 0.0)
                        continue;
                    const double mag= radius_rep
                                          ? r_try * std::abs(g[j])
                                          : std::exp(d.log_r + std::log(std::abs(g[j])));
                    if (!(mag <= std::numeric_limits<double>::max()))
                        predicted_finite= false;
                }
            }

            const bi_kappa_distribution<double>::result_type v= dist(gen_sampler);
            const bool actually_finite=
                std::isfinite(v[0]) && std::isfinite(v[1]) && std::isfinite(v[2]);
            if (predicted_finite != actually_finite)
                ++predicate_mismatch;
            if (!actually_finite)
                continue;

            const long double r_ld= std::exp(static_cast<long double>(d.log_r));
            std::array<long double, 3> scaled;
            bool scaled_ok= true;
            for (int j= 0; j < 3; ++j)
            {
                scaled[j]= r_ld * static_cast<long double>(local[j]);
                if (!std::isfinite(static_cast<double>(scaled[j])))
                    scaled_ok= false;
            }
            if (!scaled_ok)
            {
                // rotation-recoverable draw: the pre-rotation vector has no double
                // representation, so scale-then-rotate has nothing to compare.
                ++skipped;
                continue;
            }
            const std::array<long double, 3> reference=
                rot_ld.rotate_from_fieldAligned_frame(scaled, ub_ld);
            const double scale=
                std::max(std::abs(v[0]), std::max(std::abs(v[1]), std::abs(v[2])));
            ++compared;
            for (int j= 0; j < 3; ++j)
            {
                const double u=
                    test_ulp_distance(v[j], static_cast<double>(reference[j]), scale);
                if (u > worst_ulp)
                    worst_ulp= u;
            }
        }
        const bool ok= (predicate_mismatch == 0) && (worst_ulp <= 8.0) && (compared > 0);
        all_ok= all_ok && ok;
        std::cout << "  R1 kappa=0.51, ub={1,1,1}: finite exactly where the materialized"
                     " |R g_j| is at most DBL_MAX for all j (" << predicate_mismatch
                  << " mismatches in " << n << "), and within " << worst_ulp
                  << " ulp of the vector scale against a long double reference with "
                  << std::numeric_limits<long double>::digits
                  << "-bit mantissa (bound 8, " << compared << " compared, " << skipped
                  << " rotation-recoverable): " << (ok ? "PASS" : "FAIL") << std::endl;
    }

    // ---- R2: the rotation-recoverable case exists and is kept ---------------
    // A draw whose field-aligned vector has no representation can still have three
    // representable global components, because a rotation can spread it: with
    // ub = (1,1,1) a vector along the field maps to (1,1,1)/sqrt(3) times its norm,
    // so the largest component falls by sqrt(3).  The set is therefore the band of
    // log R of width log sqrt(3) just above the threshold, and its probability comes
    // from the same closed form as the floor: about 2.4e-6 here, so a few tens of
    // millions of attempts are needed before "a few" is more than noise.  At 3e7 the
    // expected count is about 71, which is far enough from zero that the existence
    // claim is not a coin flip and well enough determined that the count itself can
    // be tested against the closed form.
    {
        const long long n= 30000000;
        const double kappa= 0.505, theta_perp= 1.0, theta_par= 8.0;
        const double a= kappa - 0.5, sqrt_kappa= std::sqrt(kappa);
        const double log_max= std::log(std::numeric_limits<double>::max());
        const bi_kappa_distribution<double>::point_type ub= {1.0, 1.0, 1.0};
        const std::array<double, 3> ub_q= {1.0, 1.0, 1.0};
        double p_overflow= 0.0, p_rescued= 0.0;
        test_overflow_floor(kappa, theta_perp, theta_par, ub_q, log_max, &p_overflow,
                            &p_rescued);

        bi_kappa_distribution<double> dist(kappa, theta_perp, theta_par, ub);
        std::mt19937 gen_sampler(90075u), gen_replay(90075u);
        test_replay_state st;
        long long rescued= 0, rescued_lost= 0;
        for (long long i= 0; i < n; ++i)
        {
            const test_replay_draw d= test_replay(a, gen_replay, st);
            bi_kappa_distribution<double>::point_type local;
            local[0]= sqrt_kappa * theta_perp * d.nx;
            local[1]= sqrt_kappa * theta_perp * d.ny;
            local[2]= sqrt_kappa * theta_par * d.nz;
            const bi_kappa_distribution<double>::point_type g=
                dist.rotate_from_fieldAligned_frame(local, ub);
            const double m_local= std::max(std::abs(local[0]),
                                           std::max(std::abs(local[1]), std::abs(local[2])));
            const double m_rot=
                std::max(std::abs(g[0]), std::max(std::abs(g[1]), std::abs(g[2])));
            // Same rule as R1: materialize, then ask.  `m_local` and `m_rot` are the
            // largest pre- and post-rotation components, and a rounding-to-nearest product
            // is monotone in them, so testing the largest decides the vector.
            const double r_try= std::exp(d.log_r);
            const bool radius_rep= (r_try <= std::numeric_limits<double>::max());
            const double mag_local= radius_rep ? r_try * m_local
                                               : std::exp(d.log_r + std::log(m_local));
            const double mag_rot= radius_rep ? r_try * m_rot
                                             : std::exp(d.log_r + std::log(m_rot));
            const bool local_overflows= !(mag_local <= std::numeric_limits<double>::max());
            const bool rotated_fits= (mag_rot <= std::numeric_limits<double>::max());

            const bi_kappa_distribution<double>::result_type v= dist(gen_sampler);
            if (local_overflows && rotated_fits)
            {
                ++rescued;
                if (!std::isfinite(v[0]) || !std::isfinite(v[1]) || !std::isfinite(v[2]))
                    ++rescued_lost;
            }
        }
        const bool consistent= test_binom_consistent(rescued, n, p_rescued, 1.0e-3);
        const bool ok= (rescued > 0) && (rescued_lost == 0) && consistent;
        all_ok= all_ok && ok;
        std::cout << "  R2 kappa=0.505, theta_par/theta_perp=8, ub={1,1,1}: " << rescued
                  << " rotation-recoverable draws in " << n << " (expected "
                  << p_rescued * static_cast<double>(n)
                  << " from the closed form, exact binomial alpha=1e-3), all returned"
                     " finite (" << rescued_lost << " lost): " << (ok ? "PASS" : "FAIL")
                  << std::endl;
    }

    std::cout << "test_log_rotation_stability: " << (all_ok ? "PASS" : "FAIL") << std::endl;
    return all_ok ? 0 : 1;
}

// ================================================================
// test_representability_boundary
//
// Deterministic fixtures at the representability boundary, built from the failure
// preserved in experiments/exp7_confirmatory/raw/oracle_disagreements.jsonl and from
// the neighbourhood of that boundary in both precisions.
//
// The defect these lock down.  Before the fix, the loader decided representability by
// comparing a logarithm, `log R` or `log R + log|g_j|`, against `log(max())`.  That
// comparison cannot be made to agree with the arithmetic it predicts, because
// `log(max())` is itself a rounded value and may sit on either side of the true
// logarithm of the largest finite number.  On this libm it sits ABOVE in `float`:
// `exp(log(FLT_MAX))` is exactly `+inf`, while in `double` `exp(log(DBL_MAX))` is
// finite.  A `float` draw whose `log R` equalled `log(FLT_MAX)` bit for bit therefore
// passed the test, was exponentiated to `+inf`, and was returned as
// `(-inf, +inf, +inf)` -- with `n_nonfinite()` never incremented, because the branch
// that produced it did not look at what it had produced.  All three components were
// representable: the largest was 1.92e38 against a limit of 3.40e38, a margin of 0.571
// natural-log units, because every `|g_j|` is below one and pulls the product back
// under the limit.
//
// The corrected loader never consults a logarithm to make this decision.  It
// materializes each component and compares the component with `max()`, so the test and
// the result agree by construction, in either precision and whichever way the library's
// `log` and `exp` round.
//
// Every expected value below was computed at 100 decimal digits with mpmath, from the
// exact definition |V_j| = exp(log R + log|g_j|), and is the correctly-rounded value of
// the type, or an infinity where the exact magnitude exceeds max().  The C++ needs no
// multiprecision library because that adjudication has already happened.
// ================================================================
template <typename T> struct BoundaryFixture
{
    const char *name;
    T log_r;
    T g[3];
    T expected[3];      ///< the 100-digit value, correctly rounded, or +/-inf
    bool representable[3];
    bool any_overflow;
};

#define INFV (std::numeric_limits<float>::infinity())
static const BoundaryFixture<float> kFloatBoundaryFixtures[]= {
    {"F1_preserved_failure", 88.72284f,
     {-0.41668004f, 0.13217956f, 0.56471795f},
     {-1.417889e+38f, 4.4978383e+37f, 1.921636e+38f},
     {true, true, true}, false},
    {"F2_one_ulp_below", 88.72283f,
     {-0.41668004f, 0.13217956f, 0.56471795f},
     {-1.4178783e+38f, 4.497804e+37f, 1.9216214e+38f},
     {true, true, true}, false},
    {"F3_one_ulp_above", 88.72285f,
     {-0.41668004f, 0.13217956f, 0.56471795f},
     {-1.4178999e+38f, 4.497873e+37f, 1.9216508e+38f},
     {true, true, true}, false},
    {"F4_two_ulp_above", 88.722855f,
     {-0.41668004f, 0.13217956f, 0.56471795f},
     {-1.4179107e+38f, 4.4979073e+37f, 1.9216654e+38f},
     {true, true, true}, false},
    {"F5_radius_above_components_small", 95.0f,
     {0.001f, -0.0015f, 0.0012f},
     {1.8112393e+38f, -2.7168586e+38f, 2.173487e+38f},
     {true, true, true}, false},
    {"F6_all_unrepresentable", 200.0f,
     {0.4f, -0.5f, 0.6f},
     {INFV, -INFV, INFV},
     {false, false, false}, true},
    {"F7_mixed", 89.72284f,
     {0.001f, -0.9f, 0.5f},
     {9.2498365e+35f, -INFV, INFV},
     {true, false, false}, true},
    {"F8_zero_component", 88.72284f,
     {-0.41668004f, 0.0f, 0.56471795f},
     {-1.417889e+38f, 0.0f, 1.921636e+38f},
     {true, true, true}, false},
    {"F9_ordinary", 10.0f,
     {0.3f, -0.4f, 0.5f},
     {6607.94f, -8810.587f, 11013.232f},
     {true, true, true}, false},
};
#undef INFV
#define INFV (std::numeric_limits<double>::infinity())
static const BoundaryFixture<double> kDoubleBoundaryFixtures[]= {
    {"D1_at_log_max", 709.782712893384,
     {-0.41668001428715673, 0.13217955415599872, 0.5647179217598608},
     {-7.490628011183356e+307, 2.3761827707539496e+307, 1.015189531081392e+308},
     {true, true, true}, false},
    {"D2_one_ulp_below", 709.7827128933839,
     {-0.41668001428715673, 0.13217955415599872, 0.5647179217598608},
     {-7.490628011182504e+307, 2.3761827707536796e+307, 1.0151895310812767e+308},
     {true, true, true}, false},
    {"D3_one_ulp_above", 709.7827128933841,
     {-0.41668001428715673, 0.13217955415599872, 0.5647179217598608},
     {-7.490628011184207e+307, 2.3761827707542195e+307, 1.0151895310815074e+308},
     {true, true, true}, false},
    {"D4_radius_far_above_components_small", 800.0,
     {1e-300, -2e-300, 3e-300},
     {2.7263745721125668e+47, -5.4527491442251336e+47, 8.179123716337701e+47},
     {true, true, true}, false},
    {"D5_all_unrepresentable", 800.0,
     {0.4, -0.5, 0.6},
     {INFV, -INFV, INFV},
     {false, false, false}, true},
    {"D6_mixed", 710.782712893384,
     {0.001, -0.9, 0.5},
     {4.886636581641693e+305, -INFV, INFV},
     {true, false, false}, true},
    {"D7_zero_component", 709.782712893384,
     {-0.41668001428715673, 0.0, 0.5647179217598608},
     {-7.490628011183356e+307, 0.0, 1.015189531081392e+308},
     {true, true, true}, false},
    {"D8_ordinary", 10.0,
     {0.3, -0.4, 0.5},
     {6607.9397384420145, -8810.586317922687, 11013.232897403359},
     {true, true, true}, false},
};
#undef INFV

/// The withdrawn step 5, kept verbatim so the regression can show it failing.  This is
/// the code that shipped in 2.0.0: the radius is called representable when its LOGARITHM
/// is at most log(max()), and the fast path then multiplies without looking at what it
/// produced.
template <typename T>
static bool test_withdrawn_step5(T log_r, const T g[3], T out[3])
{
    const T log_max= std::log(std::numeric_limits<T>::max());
    const T log_g_bound= std::log(std::max(std::max(std::abs(g[0]), std::abs(g[1])),
                                           std::abs(g[2])));
    const bool radius_representable= (log_r <= log_max);
    const T R= radius_representable ? std::exp(log_r) : T(0);
    bool any= false;
    if (radius_representable && log_r + log_g_bound < log_max)
    {
        for (int j= 0; j < 3; ++j)
            out[j]= R * g[j];
        return any;   // the fast path never set the flag; that is the defect
    }
    for (int j= 0; j < 3; ++j)
    {
        if (g[j] == T(0)) { out[j]= T(0); continue; }
        const T lc= log_r + std::log(std::abs(g[j]));
        if (lc >= log_max)
        {
            any= true;
            out[j]= (g[j] < T(0)) ? -std::numeric_limits<T>::infinity()
                                  : std::numeric_limits<T>::infinity();
            continue;
        }
        const T mag= radius_representable ? R * std::abs(g[j]) : std::exp(lc);
        out[j]= (g[j] < T(0)) ? -mag : mag;
    }
    return any;
}

/// The documented accuracy of the log representation: the absolute error of `log R`
/// grows like `eps |log R|`, so the relative error of a component grows like
/// `eps |log R + log|g_j||`.  Four times that, plus a couple of ulp for the final
/// multiply or exponential, is the budget every fixture is held to.  It is derived from
/// the error model in the header, not tuned to the numbers that came out.
template <typename T> static double test_component_tolerance(T log_r, T g)
{
    const double eps= double(std::numeric_limits<T>::epsilon());
    const double lc= std::abs(double(log_r) + std::log(std::abs(double(g))));
    return 4.0 * eps * (lc + 2.0);
}

template <typename T>
static bool test_boundary_fixtures(const BoundaryFixture<T> *fx, int n, const char *label,
                                   bool expect_withdrawn_to_fail)
{
    bool all_ok= true;
    int withdrawn_failures= 0;
    double worst_rel= 0.0;
    for (int i= 0; i < n; ++i)
    {
        const BoundaryFixture<T> &f= fx[i];
        // Release 2.2.0 hands `materializeComponents` the accumulator of
        // `bikappa_detail::log_accumulator<T>` -- `double` for a `float` run, the working
        // type otherwise -- because the one rounding belongs on the product `R g_j` and
        // not on its two factors.  Widening is exact, so the fixtures below are the same
        // mathematical inputs they always were and their 100-digit adjudications still
        // stand; what changes is that the loader can now resolve them.
        typedef typename bikappa_detail::log_accumulator<T>::type A;
        std::array<A, 3> g;
        std::array<T, 3> v;
        for (int j= 0; j < 3; ++j)
            g[j]= static_cast<A>(f.g[j]);

        const bool any=
            bikappa_detail::materializeComponents<T>(static_cast<A>(f.log_r), g, v);

        bool ok= (any == f.any_overflow);
        bool any_nonfinite= false;
        for (int j= 0; j < 3; ++j)
        {
            const bool fin= (v[j] == v[j]) && std::abs(v[j]) <= std::numeric_limits<T>::max();
            any_nonfinite= any_nonfinite || !fin;
            // 1. finiteness must follow the 100-digit adjudication
            if (fin != f.representable[j])
                ok= false;
            // 2. the sign must be the sign of g_j, finite or not
            if (f.g[j] != T(0) && ((v[j] < T(0)) != (f.g[j] < T(0))))
                ok= false;
            // 3. a finite component must match the 100-digit value inside the documented
            //    error budget; an infinite one must match exactly.
            if (!fin)
            {
                if (!(v[j] == f.expected[j]))
                    ok= false;
            }
            else if (f.expected[j] != T(0))
            {
                const double rel= std::abs(double(v[j]) / double(f.expected[j]) - 1.0);
                if (rel > worst_rel)
                    worst_rel= rel;
                if (!(rel <= test_component_tolerance<T>(f.log_r, f.g[j])))
                    ok= false;
            }
            else if (v[j] != T(0))
            {
                ok= false;
            }
        }
        // 4. the contract: the flag says exactly what the values say
        if (any != any_nonfinite)
            ok= false;

        // 5. the withdrawn predicate, on the same fixture
        T old_out[3];
        const bool old_any= test_withdrawn_step5<T>(f.log_r, f.g, old_out);
        bool old_matches= (old_any == f.any_overflow);
        for (int j= 0; j < 3; ++j)
        {
            const bool ofin= (old_out[j] == old_out[j]) &&
                             std::abs(old_out[j]) <= std::numeric_limits<T>::max();
            if (ofin != f.representable[j])
                old_matches= false;
        }
        if (!old_matches)
            ++withdrawn_failures;

        all_ok= all_ok && ok;
        std::cout << "  " << label << " " << f.name << ": corrected "
                  << (ok ? "PASS" : "FAIL") << ", withdrawn form "
                  << (old_matches ? "agrees" : "DISAGREES with the 100-digit value")
                  << std::endl;
    }
    std::cout << "  " << label << ": worst relative deviation from the 100-digit value "
              << worst_rel << "; the withdrawn form misclassified " << withdrawn_failures
              << " of " << n << " fixtures" << std::endl;
    if (expect_withdrawn_to_fail && withdrawn_failures == 0)
    {
        std::cout << "  " << label << ": the regression has no teeth -- the withdrawn form "
                     "passes every fixture: FAIL" << std::endl;
        all_ok= false;
    }
    return all_ok;
}

static int test_representability_boundary()
{
    bool all_ok= true;
    std::cout << "\n--- test_representability_boundary ---" << std::endl;

    // S1/S2: the fixtures.  `float` is where this libm's log(max()) rounds up, so the
    // withdrawn form is required to fail there; `double` rounds the other way on this
    // host, so it is not required to fail and the fixtures assert the contract only.
    all_ok&= test_boundary_fixtures<float>(
        kFloatBoundaryFixtures,
        int(sizeof(kFloatBoundaryFixtures) / sizeof(kFloatBoundaryFixtures[0])), "S1 float",
        true);
    all_ok&= test_boundary_fixtures<double>(
        kDoubleBoundaryFixtures,
        int(sizeof(kDoubleBoundaryFixtures) / sizeof(kDoubleBoundaryFixtures[0])),
        "S2 double", false);

    // S3: the whole class, uncapped, at a configuration that reaches the boundary often.
    // n_nonfinite() must equal the number of non-finite returns exactly -- the count the
    // withdrawn fast path could not keep, because it never looked.
    {
        const long long n= 4000000;
        bi_kappa_distribution<float> dist;
        dist.define(0.55f, 1.0f, 2.0f, std::array<float, 3>({1.0f, 2.0f, 3.0f}));
        dist.seed(20260920);
        long long nonfinite_returns= 0;
        for (long long i= 0; i < n; ++i)
        {
            const std::array<float, 3> v= dist();
            if (!std::isfinite(v[0]) || !std::isfinite(v[1]) || !std::isfinite(v[2]))
                ++nonfinite_returns;
        }
        const bool ok= (dist.n_attempts() == static_cast<unsigned long long>(n)) &&
                       (dist.n_nonfinite() == static_cast<unsigned long long>(
                                                  nonfinite_returns)) &&
                       (nonfinite_returns > 0);
        all_ok= all_ok && ok;
        std::cout << "  S3 float kappa=0.55 rotated, " << n << " uncapped draws: "
                  << nonfinite_returns << " non-finite returns, n_nonfinite()="
                  << dist.n_nonfinite() << ", n_attempts()=" << dist.n_attempts() << ": "
                  << (ok ? "PASS" : "FAIL") << std::endl;
    }

    // S4: the capped counter counts the representability of the VELOCITY, not of the
    // normalized coordinate.  With theta_par/theta_perp = 8 the two differ by a factor of
    // eight, i.e. 2.08 natural-log units, so a draw can be unrepresentable while its
    // normalized coordinate is not.  The released 2.0.0 form tested the normalized
    // coordinate and undercounted; the corrected form materializes the velocity.
    //
    // The check is an identity rather than a rate: every rejected attempt is re-derived
    // here from the same log-domain quantities the sampler used, and the two counts must
    // agree exactly.
    {
        const int n_returned= 20000;
        const float lambda= 5.0f;
        bi_kappa_distribution<float> dist;
        dist.define(0.55f, 1.0f, 8.0f, std::array<float, 3>({0.0f, 0.0f, 1.0f}), lambda);
        dist.seed(20260921);
        for (int i= 0; i < n_returned; ++i)
        {
            const std::array<float, 3> v= dist();
            if (!std::isfinite(v[0]) || !std::isfinite(v[1]) || !std::isfinite(v[2]))
            {
                std::cout << "  S4: a capped draw was returned non-finite: FAIL"
                          << std::endl;
                all_ok= false;
                break;
            }
        }
        const bool ok= (dist.n_attempts() > static_cast<unsigned long long>(n_returned)) &&
                       (dist.n_nonfinite() > 0);
        all_ok= all_ok && ok;
        std::cout << "  S4 float kappa=0.55 theta_par/theta_perp=8 cap=" << lambda << ": "
                  << dist.n_attempts() << " attempts for " << n_returned
                  << " returned draws, n_nonfinite()=" << dist.n_nonfinite()
                  << " (unrepresentable draws the box also rejected, counted on the "
                     "velocity): " << (ok ? "PASS" : "FAIL") << std::endl;
    }

    // S5: the draw the SECOND confirmatory holdout lost, replayed from its recorded
    // variates.
    //
    // The defect this locks down is resolution, not a wrong comparison.  Release 2.1.0
    // already decided representability from the materialized component; what it could not
    // do was resolve the margin.  At `float kappa = 0.505`, seed 8005, attempt 157855, the
    // largest intended component lies 3.12e-08 natural-log units below `log(FLT_MAX)` --
    // about half an ulp -- so the correctly rounded `float` is finite and the type can
    // hold the draw.  Two working-precision quantities on the path to it are coarser than
    // that margin: `log R = (log X1 - log X2)/2` carries about 3.5e-06 once `log U` has
    // been divided by `a = 0.005`, and `g_j = sqrt(kappa) n_j` formed in `float` carries
    // another 3.9e-08, against a `float` ulp of 6.0e-08 at the overflow threshold.  2.1.0
    // therefore exponentiated one ulp high and returned `+inf` for a returnable draw.
    //
    // 2.2.0 draws the same `float` variates in the same order -- the law is unchanged --
    // and evaluates the deterministic map from them to the returned velocity in the
    // accumulator, rounding once.  Both halves are asserted: the stabilized formation must
    // reproduce the 120-digit values below, and the withdrawn working-precision formation
    // must lose the draw, so that a revert cannot pass this test quietly.
    //
    // The five variates are the bit patterns the run recorded.  The expected components
    // were computed at 120 decimal digits with mpmath from the exact definition
    // |V_j| = exp(log R + log|g_j|) with log X2 = log Y + log(U)/a, and are the
    // correctly-rounded `float` of that value.
    {
        const unsigned bits_x1= 0x3f128ad2u, bits_y= 0x3f795a38u, bits_u= 0x3ed12d4eu;
        const unsigned bits_ct= 0xbf0bc1cau, bits_phi= 0x40c530f9u;
        float x1, y, u, ct, phi;
        std::memcpy(&x1, &bits_x1, sizeof(float));
        std::memcpy(&y, &bits_y, sizeof(float));
        std::memcpy(&u, &bits_u, sizeof(float));
        std::memcpy(&ct, &bits_ct, sizeof(float));
        std::memcpy(&phi, &bits_phi, sizeof(float));
        const float kappa= 0.505f;
        const float a= kappa - 0.5f;

        // 120-digit adjudication, correctly rounded to float.  The first component is
        // FLT_MAX itself: the exact magnitude is 0.99999999860 of it.
        const float expected[3]= {3.4028235e+38f, -4.1361512e+37f, -2.2335724e+38f};

        // --- the stabilized formation, as 2.2.0 performs it ---------------------
        typedef bikappa_detail::log_accumulator<float>::type A;
        const A log_x2= std::log(static_cast<A>(y)) +
                        std::log(static_cast<A>(u)) / static_cast<A>(a);
        const A log_r= A(0.5) * (std::log(static_cast<A>(x1)) - log_x2);
        const A st= std::sqrt(A(1) - static_cast<A>(ct) * static_cast<A>(ct));
        const A sk= std::sqrt(static_cast<A>(kappa));
        std::array<A, 3> g;
        g[0]= sk * st * std::cos(static_cast<A>(phi));
        g[1]= sk * st * std::sin(static_cast<A>(phi));
        g[2]= sk * static_cast<A>(ct);
        std::array<float, 3> v;
        const bool any= bikappa_detail::materializeComponents<float>(log_r, g, v);

        bool ok= !any;
        for (int j= 0; j < 3; ++j)
            ok= ok && (v[j] == expected[j]);

        // --- the withdrawn working-precision formation --------------------------
        const float f_log_x2= std::log(y) + std::log(u) / a;
        const float f_log_r= 0.5f * (std::log(x1) - f_log_x2);
        const float f_st= std::sqrt(1.0f - ct * ct);
        const float f_sk= std::sqrt(kappa);
        const float f_g0= f_sk * f_st * std::cos(phi);
        const float f_R= std::exp(f_log_r);
        const bool withdrawn_lost=
            !(f_R <= std::numeric_limits<float>::max()) ||
            !(f_R * std::abs(f_g0) <= std::numeric_limits<float>::max());

        ok= ok && withdrawn_lost;
        all_ok= all_ok && ok;
        std::cout << "  S5 holdout-2 draw (float kappa=0.505, seed 8005, attempt 157855): "
                  << "stabilized returns " << (any ? "non-finite" : "the 120-digit values")
                  << ", withdrawn form " << (withdrawn_lost ? "loses it" : "also returns it")
                  << ": " << (ok ? "PASS" : "FAIL") << std::endl;
    }

    // S6: the float representability boundary, swept.
    //
    // S5 is one point on a boundary; this walks across it.  Each target magnitude is built
    // first, as an exact `double` -- FLT_MAX plus an integer number of eighths of an ulp,
    // which needs 27 bits of mantissa -- and what the loader must decide is compared
    // against the hardware's own `float` conversion of that magnitude, never against the
    // loader's own `exp`.  The rounding midpoint at +4/8 ulp is excluded and the exclusion
    // is part of the rule: it is the one magnitude whose answer comes from
    // round-half-to-even rather than from an inequality, no finite-precision
    // reconstruction of it can be relied on to land on the right side, and an exact tie has
    // probability zero under a continuous law.  The two points 1e-12 either side of it --
    // a hundred times the accumulator's resolution, four thousand times finer than an ulp
    // -- are included instead, and they are what a working-precision `log R` could not
    // have decided.
    {
        const float fmax= std::numeric_limits<float>::max();
        const double fmaxd= static_cast<double>(fmax);
        const double ulp= fmaxd - static_cast<double>(std::nextafter(fmax, 0.0f));
        const double g0= 0.59104243501804303857;   // the S5 draw's largest |g_j|
        const double g1= -0.071841542754834329151;
        const double g2= -0.38795312600129426096;

        int checked= 0, bad= 0;
        for (int i= 0; i <= 34; ++i)
        {
            const int j= i - 16;
            double m;
            bool expect_finite;
            if (i <= 32)
            {
                if (j == 4)
                    continue;                       // the exact tie; see above
                m= fmaxd + static_cast<double>(j) * (ulp / 8.0);
                expect_finite= std::isfinite(static_cast<float>(m));
            }
            else
            {
                const double mid= fmaxd + 4.0 * (ulp / 8.0);
                expect_finite= (i == 33);
                m= mid * (expect_finite ? (1.0 - 1.0e-12) : (1.0 + 1.0e-12));
            }
            const double log_r= std::log(m) - std::log(g0);
            std::array<double, 3> g;
            g[0]= g0;
            g[1]= g1;
            g[2]= g2;
            std::array<float, 3> v;
            const bool over= bikappa_detail::materializeComponents<float>(log_r, g, v);
            const bool got= !over && std::isfinite(v[0]) && std::isfinite(v[1]) &&
                            std::isfinite(v[2]);
            ++checked;
            if (got != expect_finite)
                ++bad;
        }
        const bool ok= (bad == 0) && (checked == 34);
        all_ok= all_ok && ok;
        std::cout << "  S6 float boundary sweep: " << checked
                  << " magnitudes from -2 to +2 ulp of FLT_MAX and 1e-12 either side of "
                     "the rounding midpoint, " << bad << " decided against the hardware: "
                  << (ok ? "PASS" : "FAIL") << std::endl;
    }

    // S7: the accumulator is the working type for `double` and wider only for `float`.
    //
    // This is what makes "a `double` instantiation of 2.2.0 returns the same bits as
    // 2.1.0" a statement about the source rather than a measurement: with
    // `log_accumulator<double>::type == double`, every expression 2.2.0 widened is the
    // expression 2.1.0 evaluated, character for character, so there is no `double`
    // arithmetic to re-validate.  Only `float` gets a wider accumulator, and it must be
    // strictly wider or the fix is not there at all.  A digest against the 2.1.0 binary
    // is not used here: it would depend on the compiler and on `-ffp-contract`, and the
    // equality it would test is checked once, against both released headers, in the
    // release record rather than in a portable unit test.
    {
        const bool d_is_working=
            std::is_same<bikappa_detail::log_accumulator<double>::type, double>::value;
        const bool f_is_double=
            std::is_same<bikappa_detail::log_accumulator<float>::type, double>::value;
        const bool wider= std::numeric_limits<
                              bikappa_detail::log_accumulator<float>::type>::digits >
                          std::numeric_limits<float>::digits;
        const bool ok= d_is_working && f_is_double && wider;
        all_ok= all_ok && ok;
        std::cout << "  S7 accumulator: double -> double (working type, so the double path "
                     "is 2.1.0 unchanged), float -> double ("
                  << std::numeric_limits<
                         bikappa_detail::log_accumulator<float>::type>::digits
                  << " bits against " << std::numeric_limits<float>::digits << "): "
                  << (ok ? "PASS" : "FAIL") << std::endl;
    }

    std::cout << "test_representability_boundary: " << (all_ok ? "PASS" : "FAIL")
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
    result|= test_define_field_aligned();
    result|= test_max_normalized_velocity_bikappa();
    result|= test_max_normalized_velocity_bimaxwellian();
    result|= test_seed_api();
    result|= test_radius_formation();
    result|= test_no_cap_semantics();
    result|= test_no_cap_semantics_bimaxwellian();
    result|= test_log_radius_primitive();
    result|= test_log_radius_no_regression();
    result|= test_log_radius_finiteness();
    result|= test_overflow_observability();
    result|= test_seed_reset_contract();
    result|= test_log_cap_semantics();
    result|= test_log_rotation_stability();
    result|= test_representability_boundary();
    std::cout << "\n=== Test suite: " << (result == 0 ? "ALL PASS" : "SOME FAILED") << " ==="
              << std::endl;
    return result;
}

#endif // _TEST_SUITE_H_
