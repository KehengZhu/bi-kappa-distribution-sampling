// Experiment 4 — finite-precision / low-kappa audit of the released bi-Kappa sampler.
//
// The question this program answers is NOT "can we invent a new low-kappa algorithm".
// It is: what is the actual supported numerical range of the released implementation as a
// function of kappa, precision and standard-library implementation, and how much of the
// observed failure is avoidable (an artifact of how the radius is formed) versus honest
// (the mathematical sample itself is not representable at that precision).
//
// Four modes, selected by argv[1]:
//
//   released   Draw from the released class bi_kappa_distribution<T> in no_cap mode and
//              count non-finite outputs.  This is the ground truth for "what the shipped
//              code does".  No instrumentation is added to the released header.
//
//   variates   The mechanism decomposition.  A reference small-shape Gamma generator is
//              driven so that the SAME underlying variates feed four different radius
//              formations, which makes the comparison exact draw-by-draw rather than
//              distributional:
//                  r_ratio = sqrt(x1 / x2)      the pre-fix formation
//                  r_split = sqrt(x1)/sqrt(x2)  the current formation
//                  log_r   = (log x1 - log x2)/2, carried in the log domain throughout
//              log_r never underflows or overflows, so it is the arbiter that decides
//              whether a lost draw was recoverable.
//
//   logdump    Dump log R from the log-domain path for distributional validation.
//
//   capped     Q5.  Does the component-wise box cap prevent the invalid computation, or
//              merely resample it away after the fact?
//
// Build: see GNUmakefile.  Every mode prints one JSON object per configuration to stdout.

#include "../../cpp/bi_kappa_distribution.H"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <random>
#include <string>
#include <vector>

// ---------------------------------------------------------------------------
// Environment identification
// ---------------------------------------------------------------------------
static const char *compilerId()
{
#if defined(__clang__)
    return "clang " __clang_version__;
#elif defined(__GNUC__)
    return "gcc " __VERSION__;
#else
    return "unknown";
#endif
}

static const char *stdlibId()
{
#if defined(_LIBCPP_VERSION)
    return "libc++";
#elif defined(__GLIBCXX__)
    return "libstdc++";
#else
    return "unknown";
#endif
}

static const char *archId()
{
#if defined(__aarch64__) || defined(_M_ARM64)
    return "arm64";
#elif defined(__x86_64__)
    return "x86_64";
#else
    return "unknown";
#endif
}

template <typename T> const char *precisionId();
template <> const char *precisionId<float>() { return "float"; }
template <> const char *precisionId<double>() { return "double"; }

// ---------------------------------------------------------------------------
// Reference small-shape Gamma, in both plain and log form, from shared variates.
//
// For shape a in (0,1) the standard boost identity is
//     X ~ Ga(a,1)   <=>   X = Y * U^(1/a),   Y ~ Ga(a+1,1),  U ~ Uniform(0,1),
// which is what Marsaglia-Tsang style generators (and hence both libc++ and libstdc++)
// use to reach a < 1.  Shape a+1 lies in (1,2) and is numerically benign, so ALL of the
// small-shape difficulty is concentrated in the single factor U^(1/a):
//
//     X     = Y * pow(U, 1/a)          underflows to exactly 0 for small a
//     log X = log(Y) + log(U)/a        exact, and cannot underflow
//
// Computing both from one (Y, U) pair is what makes the spurious/honest split exact.
// ---------------------------------------------------------------------------
template <typename T> struct GammaDraw
{
    T x;     // the value a standard implementation would produce (may be 0 or inf)
    T log_x; // the same draw carried in the log domain (always finite)
};

template <typename T, typename Gen>
GammaDraw<T> drawGammaSmallShape(T a, Gen &gen, std::gamma_distribution<T> &gBoost,
                                 std::uniform_real_distribution<T> &unif)
{
    GammaDraw<T> out;
    if (a >= T(1))
    {
        // No boost needed; shape >= 1 is well conditioned.
        std::gamma_distribution<T> g(a, T(1));
        out.x= g(gen);
        out.log_x= std::log(out.x);
        return out;
    }
    const T y= gBoost(gen); // Ga(a+1, 1)
    T u= unif(gen);
    while (!(u > T(0)))
    {
        u= unif(gen); // guard the open-interval endpoint
    }
    out.x= y * std::pow(u, T(1) / a);
    out.log_x= std::log(y) + std::log(u) / a;
    return out;
}

// ---------------------------------------------------------------------------
// Per-configuration accumulators
// ---------------------------------------------------------------------------
struct Counters
{
    long long n= 0;
    long long x2_zero= 0;         // denominator Gamma underflowed to exactly 0
    long long x2_subnormal= 0;    // denominator in the subnormal range (still non-zero)
    long long ratio_nonfinite= 0; // sqrt(x1/x2) not finite
    long long split_nonfinite= 0; // sqrt(x1)/sqrt(x2) not finite
    long long recovered= 0;       // ratio lost it, split kept it
    long long honest= 0;          // true |log r| exceeds what the precision can represent
    long long spurious_ratio= 0;  // representable, but the ratio formation lost it
    long long spurious_split= 0;  // representable, but the split formation still lost it
    double max_finite_log_r= -std::numeric_limits<double>::infinity();
    double max_log_r= -std::numeric_limits<double>::infinity(); // over ALL draws, log domain
};

// ---------------------------------------------------------------------------
// Mode: variates
// ---------------------------------------------------------------------------
template <typename T>
Counters runVariates(double kappa, unsigned seed, long long nDraws)
{
    Counters c;
    const T a= T(kappa) - T(0.5);
    const T logMax= std::log(std::numeric_limits<T>::max());

    std::mt19937 gen(seed);
    std::gamma_distribution<T> gX1(T(1.5), T(1));
    std::gamma_distribution<T> gBoost(a + T(1), T(1));
    std::uniform_real_distribution<T> unif(T(0), T(1));

    for (long long i= 0; i < nDraws; ++i)
    {
        const T x1= gX1(gen);
        const GammaDraw<T> d2= drawGammaSmallShape<T>(a, gen, gBoost, unif);
        const T x2= d2.x;

        const T rRatio= std::sqrt(x1 / x2);
        const T rSplit= std::sqrt(x1) / std::sqrt(x2);
        const T logR= T(0.5) * (std::log(x1) - d2.log_x);

        ++c.n;
        if (x2 == T(0))
            ++c.x2_zero;
        else if (std::fabs(x2) < std::numeric_limits<T>::min())
            ++c.x2_subnormal;

        const bool ratioLost= !std::isfinite(rRatio);
        const bool splitLost= !std::isfinite(rSplit);
        // "Representable" means the mathematical radius itself fits in this precision.
        // Judged in the log domain, which neither underflows nor overflows.
        const bool representable= (logR <= logMax);

        if (ratioLost)
            ++c.ratio_nonfinite;
        if (splitLost)
            ++c.split_nonfinite;
        if (ratioLost && !splitLost)
            ++c.recovered;
        if (!representable)
            ++c.honest;
        if (ratioLost && representable)
            ++c.spurious_ratio;
        if (splitLost && representable)
            ++c.spurious_split;

        const double logRd= static_cast<double>(logR);
        if (logRd > c.max_log_r)
            c.max_log_r= logRd;
        if (std::isfinite(rSplit))
        {
            const double lr= std::log(static_cast<double>(rSplit));
            if (lr > c.max_finite_log_r)
                c.max_finite_log_r= lr;
        }
    }
    return c;
}

// ---------------------------------------------------------------------------
// Mode: released — the shipped class, uncapped, counting non-finite components
// ---------------------------------------------------------------------------
template <typename T>
void runReleased(double kappa, unsigned seed, long long nDraws, long long *nonFinite,
                 long long *thrown, double *maxFiniteLogR)
{
    *nonFinite= 0;
    *thrown= 0;
    *maxFiniteLogR= -std::numeric_limits<double>::infinity();

    bi_kappa_distribution<T> dist(T(kappa), T(1), T(1),
                                  {T(0), T(0), T(1)},
                                  bi_kappa_distribution<T>::no_cap());
    dist.seed(static_cast<int>(seed));

    for (long long i= 0; i < nDraws; ++i)
    {
        try
        {
            const auto v= dist();
            const bool bad= !std::isfinite(v[0]) || !std::isfinite(v[1]) || !std::isfinite(v[2]);
            if (bad)
            {
                ++(*nonFinite);
            }
            else
            {
                // overflow-safe norm: hypot chains, never forms the sum of squares
                const double r= std::hypot(std::hypot(static_cast<double>(v[0]),
                                                      static_cast<double>(v[1])),
                                           static_cast<double>(v[2]));
                if (r > 0 && std::isfinite(r))
                {
                    const double lr= std::log(r);
                    if (lr > *maxFiniteLogR)
                        *maxFiniteLogR= lr;
                }
            }
        }
        catch (const std::exception &)
        {
            ++(*thrown);
        }
    }
}

// ---------------------------------------------------------------------------
// Mode: capped — Q5.  Does the cap prevent the bad computation or hide it?
// ---------------------------------------------------------------------------
template <typename T>
void runCapped(double kappa, double lambda, unsigned seed, long long nDraws,
               long long *nonFiniteOut, long long *thrown, long long *nonFiniteRedrawn,
               long long *totalAttempts)
{
    *nonFiniteOut= 0;
    *thrown= 0;
    *nonFiniteRedrawn= 0;
    *totalAttempts= 0;

    // Reproduce the sampler's inner loop faithfully so that attempts can be counted
    // without instrumenting the released header.  Same variate order, same predicate.
    const T a= T(kappa) - T(0.5);
    const T sqrtKappa= std::sqrt(T(kappa));
    const T thetaPerp= T(1), thetaPar= T(1);

    std::mt19937 gen(seed);
    std::gamma_distribution<T> gX1(T(1.5), T(1));
    std::gamma_distribution<T> gX2(a, T(1));
    std::uniform_real_distribution<T> uCos(T(-1), T(1));
    std::uniform_real_distribution<T> uPhi(T(0), T(2) * T(3.14159265358979323846));

    const int kMaxTries= 1000000;
    for (long long i= 0; i < nDraws; ++i)
    {
        bool accepted= false;
        for (int t= 0; t < kMaxTries; ++t)
        {
            ++(*totalAttempts);
            const T x1= gX1(gen);
            const T x2= gX2(gen);
            const T r= std::sqrt(x1) / std::sqrt(x2);
            const T cosTheta= uCos(gen);
            const T sinTheta= std::sqrt(std::max(T(0), T(1) - cosTheta * cosTheta));
            const T phi= uPhi(gen);

            const T px= sqrtKappa * thetaPerp * r * sinTheta * std::cos(phi);
            const T py= sqrtKappa * thetaPerp * r * sinTheta * std::sin(phi);
            const T pz= sqrtKappa * thetaPar * r * cosTheta;

            const bool bad= !std::isfinite(px) || !std::isfinite(py) || !std::isfinite(pz);
            const bool inBox= std::fabs(px) / thetaPerp <= T(lambda) &&
                              std::fabs(py) / thetaPerp <= T(lambda) &&
                              std::fabs(pz) / thetaPar <= T(lambda);
            if (bad)
            {
                // A non-finite draw can never satisfy the box predicate, so the cap
                // silently redraws it.  Count how often that happens.
                ++(*nonFiniteRedrawn);
            }
            if (inBox)
            {
                if (bad)
                    ++(*nonFiniteOut);
                accepted= true;
                break;
            }
        }
        if (!accepted)
            ++(*thrown);
    }
}

// ---------------------------------------------------------------------------
// Mode: logdump — log R samples for distributional validation
// ---------------------------------------------------------------------------
template <typename T>
void runLogDump(double kappa, unsigned seed, long long nDraws, const char *path)
{
    const T a= T(kappa) - T(0.5);
    std::mt19937 gen(seed);
    std::gamma_distribution<T> gX1(T(1.5), T(1));
    std::gamma_distribution<T> gBoost(a + T(1), T(1));
    std::uniform_real_distribution<T> unif(T(0), T(1));

    std::vector<double> buf;
    buf.reserve(static_cast<size_t>(nDraws));
    for (long long i= 0; i < nDraws; ++i)
    {
        const T x1= gX1(gen);
        const GammaDraw<T> d2= drawGammaSmallShape<T>(a, gen, gBoost, unif);
        buf.push_back(0.5 * (static_cast<double>(std::log(x1)) - static_cast<double>(d2.log_x)));
    }
    FILE *f= std::fopen(path, "wb");
    if (!f)
    {
        std::fprintf(stderr, "cannot open %s\n", path);
        std::exit(1);
    }
    std::fwrite(buf.data(), sizeof(double), buf.size(), f);
    std::fclose(f);
}

// ---------------------------------------------------------------------------

static void emitEnv(const char *precision)
{
    std::printf("\"compiler\":\"%s\",\"stdlib\":\"%s\",\"arch\":\"%s\",\"precision\":\"%s\","
                "\"rng\":\"mt19937\"",
                compilerId(), stdlibId(), archId(), precision);
}

template <typename T> void variatesSweep(const std::vector<double> &kappas,
                                         const std::vector<unsigned> &seeds, long long n)
{
    for (double k : kappas)
    {
        for (unsigned s : seeds)
        {
            const Counters c= runVariates<T>(k, s, n);
            std::printf("{\"mode\":\"variates\",");
            emitEnv(precisionId<T>());
            std::printf(",\"kappa\":%.6g,\"shape_x2\":%.6g,\"seed\":%u,\"n\":%lld,"
                        "\"x2_zero\":%lld,\"x2_subnormal\":%lld,"
                        "\"ratio_nonfinite\":%lld,\"split_nonfinite\":%lld,"
                        "\"recovered_by_split\":%lld,\"honest_overflow\":%lld,"
                        "\"spurious_ratio\":%lld,\"spurious_split\":%lld,"
                        "\"max_finite_log_r\":%.10g,\"max_log_r\":%.10g}\n",
                        k, k - 0.5, s, c.n, c.x2_zero, c.x2_subnormal, c.ratio_nonfinite,
                        c.split_nonfinite, c.recovered, c.honest, c.spurious_ratio,
                        c.spurious_split, c.max_finite_log_r, c.max_log_r);
        }
    }
}

template <typename T> void releasedSweep(const std::vector<double> &kappas,
                                         const std::vector<unsigned> &seeds, long long n)
{
    for (double k : kappas)
    {
        for (unsigned s : seeds)
        {
            long long nf= 0, th= 0;
            double mx= 0;
            runReleased<T>(k, s, n, &nf, &th, &mx);
            std::printf("{\"mode\":\"released\",");
            emitEnv(precisionId<T>());
            std::printf(",\"kappa\":%.6g,\"seed\":%u,\"n\":%lld,\"nonfinite\":%lld,"
                        "\"thrown\":%lld,\"max_finite_log_r\":%.10g}\n",
                        k, s, n, nf, th, mx);
        }
    }
}

template <typename T> void cappedSweep(const std::vector<double> &kappas,
                                       const std::vector<double> &lambdas,
                                       const std::vector<unsigned> &seeds, long long n)
{
    for (double k : kappas)
    {
        for (double lam : lambdas)
        {
            for (unsigned s : seeds)
            {
                long long nfo= 0, th= 0, nfr= 0, att= 0;
                runCapped<T>(k, lam, s, n, &nfo, &th, &nfr, &att);
                std::printf("{\"mode\":\"capped\",");
                emitEnv(precisionId<T>());
                std::printf(",\"kappa\":%.6g,\"lambda\":%.6g,\"seed\":%u,\"n\":%lld,"
                            "\"attempts\":%lld,\"nonfinite_redrawn\":%lld,"
                            "\"nonfinite_returned\":%lld,\"exhausted\":%lld}\n",
                            k, lam, s, n, att, nfr, nfo, th);
            }
        }
    }
}

int main(int argc, char **argv)
{
    if (argc < 2)
    {
        std::fprintf(stderr,
                     "usage: %s (released|variates|capped) [n_draws]\n"
                     "       %s logdump <kappa> <seed> <n> <out.bin> <precision:float|double>\n",
                     argv[0], argv[0]);
        return 2;
    }
    const std::string mode= argv[1];

    // The kappa ladder is fixed here, in the committed source, so that a run is fully
    // determined by the command line recorded in the README.
    const std::vector<double> kappas= {0.5001, 0.501, 0.505, 0.51, 0.55,
                                       0.60,   0.75,  1.0,   1.5};
    const std::vector<unsigned> seeds= {4001, 4002, 4003, 4004, 4005};

    if (mode == "logdump")
    {
        if (argc < 7)
        {
            std::fprintf(stderr, "logdump needs <kappa> <seed> <n> <out.bin> <precision>\n");
            return 2;
        }
        const double k= std::atof(argv[2]);
        const unsigned s= static_cast<unsigned>(std::atoi(argv[3]));
        const long long n= std::atoll(argv[4]);
        const char *path= argv[5];
        if (std::strcmp(argv[6], "float") == 0)
            runLogDump<float>(k, s, n, path);
        else
            runLogDump<double>(k, s, n, path);
        return 0;
    }

    const long long n= (argc >= 3) ? std::atoll(argv[2]) : 1000000LL;

    if (mode == "variates")
    {
        variatesSweep<float>(kappas, seeds, n);
        variatesSweep<double>(kappas, seeds, n);
    }
    else if (mode == "released")
    {
        releasedSweep<float>(kappas, seeds, n);
        releasedSweep<double>(kappas, seeds, n);
    }
    else if (mode == "capped")
    {
        // The capped sweep is deliberately much smaller than the uncapped ones.  Near
        // kappa = 1/2 the box acceptance probability is ~4e-4, so each accepted sample
        // costs thousands of attempts; 2000 accepted samples per configuration already
        // costs ~5e6 attempts and is ample for the question being asked, which is
        // whether non-finite draws reach the caller or are silently redrawn.
        const std::vector<double> lambdas= {5.0, 20.0};
        const std::vector<double> capKappas= {0.5001, 0.51, 0.55, 0.75, 1.5};
        const long long nCap= 2000;
        cappedSweep<float>(capKappas, lambdas, seeds, nCap);
        cappedSweep<double>(capKappas, lambdas, seeds, nCap);
    }
    else
    {
        std::fprintf(stderr, "unknown mode: %s\n", mode.c_str());
        return 2;
    }
    return 0;
}
