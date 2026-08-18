// Experiment 3 -- reproducible absolute and comparative performance benchmark.
//
// Answers R1.4 and the "how does it outperform" part of R2.A2.
//
// THE RULE THIS FILE EXISTS TO OBEY: no timing number is meaningful until every
// compared method is verified to sample the intended target law.  So this program
// has two modes, and `validate` must pass before `time` is believed.
//
//   validate   dump normalized radii + directions per method -> distributional check
//   time       repeated timing batches per (method, kappa) -> ns/sample
//
// ---------------------------------------------------------------------------
// The three methods, each transcribed from its primary source.  Two of them are
// NOT independent algorithms and are labelled as such -- inflating an
// implementation variant into a rival algorithm is exactly what R2.A2 warns about.
//
// M1  gamma_ratio_spherical   THE RELEASED IMPLEMENTATION.
//     cpp/bi_kappa_distribution.H in no_cap() mode.  X1 ~ Ga(3/2,1),
//     X2 ~ Ga(kappa-1/2,1), R = sqrt(X1)/sqrt(X2), then an explicit uniform
//     direction on S^2.  Equivalent to Zenitani & Nakano (2022) Algorithm 1-1 and
//     to ZUM (2026) Algorithm 3.1.
//
// M2  scale_mixture_normals   Abdul & Mace (2015), Phys. Plasmas 22, 102107,
//     Eq. (22) with Eqs. (19)-(20):
//         X = mu + sigma * sqrt(nu / chisq_nu) * Z,
//         nu = 2*kappa - 1,   sigma^2 = kappa*theta^2/(2*kappa - 1),
//     with Z three independent standard normals and ONE chi-squared deviate of
//     nu degrees of freedom shared across all three components (that shared
//     denominator is what makes it a genuine trivariate scale mixture).
//     Substituting sigma and nu collapses to
//         v_i = theta*sqrt(kappa) * Z_i / sqrt(chisq_nu).
//     SAME TARGET LAW AS M1, and provably the same construction: |Z|^2 ~ chisq_3
//     = 2*Ga(3/2,1) and chisq_nu = 2*Ga(kappa-1/2,1), so |v| has exactly M1's
//     Gamma-ratio radius.  It is therefore an IMPLEMENTATION VARIANT of M1 -- it
//     buys the direction out of three normals instead of two uniforms -- and the
//     benchmark treats it as one.
//     ** Disclosure: the paper never says how the non-integer-nu chi-squared
//     deviate is generated.  chisq_nu = 2*Ga(nu/2,1) via std::gamma_distribution
//     is OUR choice, not theirs, and any cost difference from a different
//     chi-squared route is not attributable to Abdul & Mace.
//
// M3  pareto_rejection        Zenitani (2025), Res. Notes AAS 9, 299.
//     A GENUINELY DIFFERENT ALGORITHM: rejection sampling of the beta-prime
//     radial law under a Pareto envelope, using uniform variates ONLY (the point
//     of the paper is portability -- no Gamma or normal generator needed).
//     Transcribed verbatim from its Section 2 procedure:
//         Step 0: D <- (2k-2n-1)^(k-n-1/2) * (2k-2n)^(n-k)
//         Step 1: U1, U2 ~ U(0,1)
//         Step 2: W <- sqrt((1-U1)^(-1/n) - 1)
//         Step 3: if W*(1-U1)^((k-n)/n) < D*U2  -> back to Step 1
//         Step 4: U3, U4 ~ U(0,1)
//         Step 5: V <- sqrt(k*theta^2)*W
//                 v_x <- V*(2*U3-1)
//                 v_y <- 2*V*sqrt(U3*(1-U3))*cos(2*pi*U4)
//                 v_z <- 2*V*sqrt(U3*(1-U3))*sin(2*pi*U4)
//     Envelope index n must satisfy 0 < n < kappa - 1/2.  We use the author's
//     recommendation n = kappa/2, which he reports gives efficiency ~0.73-0.8 and
//     simplifies the Step 2/3 exponents to -2/kappa and 1.  NOTE n = kappa/2
//     requires kappa > 1: at kappa = 1 it hits the n = kappa - 1/2 bound and D
//     degenerates.  The author's quoted efficiencies start at kappa = 1.5.
//     ** So M3 is INAPPLICABLE for kappa <= 1 under its recommended setting, and
//     the harness records that as a result rather than silently skipping it.
// ---------------------------------------------------------------------------
//
// Fairness rules, all deliberate:
//   * one RNG type for every method (std::mt19937), seeded identically per run;
//   * the isotropic core (theta_perp = theta_par = 1) is what is compared, since
//     M2 and M3 as published are isotropic;  the released anisotropic + field
//     rotation path is timed separately as M1's own variants, not against them;
//   * timing loops accumulate a checksum so the optimizer cannot delete the work;
//   * every configuration is timed in several independent batches and the
//     DISTRIBUTION of batch times is reported, never a single number.

#include "../../cpp/bi_kappa_distribution.H"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <random>
#include <string>
#include <vector>

// --------------------------------------------------------------------------- //
// environment capture
// --------------------------------------------------------------------------- //
static const char *compiler_id()
{
#if defined(__clang__)
    return "clang " __clang_version__;
#elif defined(__GNUC__)
    return "gcc " __VERSION__;
#else
    return "unknown";
#endif
}

static const char *stdlib_id()
{
#if defined(_LIBCPP_VERSION)
    return "libc++";
#elif defined(__GLIBCXX__)
    return "libstdc++";
#else
    return "unknown";
#endif
}

static const char *arch_id()
{
#if defined(__aarch64__) || defined(_M_ARM64)
    return "arm64";
#elif defined(__x86_64__)
    return "x86_64";
#else
    return "unknown";
#endif
}

static double pi_value() { return 3.14159265358979323846264338327950288; }

// --------------------------------------------------------------------------- //
// M2 -- Abdul & Mace (2015) Eq. (22)
// --------------------------------------------------------------------------- //
struct ScaleMixtureNormals
{
    double kappa, theta;
    std::normal_distribution<double> z;
    std::gamma_distribution<double> g;   // chisq_nu = 2 * Ga(nu/2, 1)

    ScaleMixtureNormals(double k, double th)
        : kappa(k), theta(th), z(0.0, 1.0), g(k - 0.5, 1.0)
    {
    }

    // Returns the three velocity components.
    void draw(std::mt19937 &gen, double out[3])
    {
        const double chisq= 2.0 * g(gen);
        const double scale= theta * std::sqrt(kappa) / std::sqrt(chisq);
        out[0]= scale * z(gen);
        out[1]= scale * z(gen);
        out[2]= scale * z(gen);
    }
};

// --------------------------------------------------------------------------- //
// M3 -- Zenitani (2025) Pareto-envelope rejection
// --------------------------------------------------------------------------- //
struct ParetoRejection
{
    double kappa, theta, n, D, inv_n, exp3;
    bool applicable;
    long attempts;   // counts EVERY envelope draw, so acceptance is measurable
    long accepted;   // counted here rather than inferred, so warm-up batches and
                     // any future change to the batch schedule cannot corrupt it

    ParetoRejection(double k, double th)
        : kappa(k), theta(th), attempts(0), accepted(0)
    {
        n= 0.5 * k;                       // author's recommendation
        applicable= (n > 0.0) && (n < k - 0.5);
        // Step 0.  D = (2k-2n-1)^(k-n-1/2) * (2k-2n)^(n-k)
        const double a= 2.0 * k - 2.0 * n - 1.0;
        const double b= 2.0 * k - 2.0 * n;
        D= std::pow(a, k - n - 0.5) * std::pow(b, n - k);
        inv_n= -1.0 / n;                  // Step 2 exponent
        exp3= (k - n) / n;                // Step 3 exponent
    }

    void draw(std::mt19937 &gen, std::uniform_real_distribution<double> &u, double out[3])
    {
        double W;
        for (;;)
        {
            ++attempts;
            const double u1= u(gen);
            const double u2= u(gen);
            const double om= 1.0 - u1;
            W= std::sqrt(std::pow(om, inv_n) - 1.0);           // Step 2
            if (!(W * std::pow(om, exp3) < D * u2))            // Step 3
                break;
        }
        ++accepted;
        const double u3= u(gen);                                // Step 4
        const double u4= u(gen);
        const double V= std::sqrt(kappa * theta * theta) * W;   // Step 5
        const double s= 2.0 * V * std::sqrt(u3 * (1.0 - u3));
        out[0]= V * (2.0 * u3 - 1.0);
        out[1]= s * std::cos(2.0 * pi_value() * u4);
        out[2]= s * std::sin(2.0 * pi_value() * u4);
    }
};

// --------------------------------------------------------------------------- //
// validate -- dump log|v| and the direction cosine for a distributional check
// --------------------------------------------------------------------------- //
// log|v| is dumped, never |v|^2: at low kappa the square overflows where |v| is
// still representable, which would discard exactly the heavy-tail draws the check
// exists to test.  The norm uses hypot, never a naive sum of squares.
static double log_norm(const double v[3])
{
    const double r= std::hypot(std::hypot(v[0], v[1]), v[2]);
    if (r > 0.0 && std::isfinite(r))
        return std::log(r);
    // Fall back to a scaled evaluation only if hypot itself overflowed.
    const double m= std::fmax(std::fabs(v[0]), std::fmax(std::fabs(v[1]), std::fabs(v[2])));
    if (!(m > 0.0) || !std::isfinite(m))
        return std::numeric_limits<double>::quiet_NaN();
    const double a= v[0] / m, b= v[1] / m, c= v[2] / m;
    return std::log(m) + 0.5 * std::log(a * a + b * b + c * c);
}

static int run_validate(const char *outpath, double kappa, double theta, int seed, long n)
{
    FILE *fh= std::fopen(outpath, "wb");
    if (!fh)
    {
        std::fprintf(stderr, "cannot open %s\n", outpath);
        return 1;
    }

    // Three records per draw: method id, log|v|, cos(theta_dir).  Written as
    // doubles so the analysis never re-derives anything the sampler already knew.
    std::vector<double> buf;
    buf.reserve(3 * 3 * 1024);

    std::mt19937 gen1(static_cast<unsigned>(seed));
    std::mt19937 gen2(static_cast<unsigned>(seed));
    std::mt19937 gen3(static_cast<unsigned>(seed));

    bi_kappa_distribution<double> m1(kappa, theta, theta, {0.0, 0.0, 1.0},
                                     bi_kappa_distribution<double>::no_cap());
    ScaleMixtureNormals m2(kappa, theta);
    ParetoRejection m3(kappa, theta);
    std::uniform_real_distribution<double> uni(0.0, 1.0);

    for (long i= 0; i < n; ++i)
    {
        double v[3];

        const auto a= m1(gen1);
        v[0]= a[0]; v[1]= a[1]; v[2]= a[2];
        double lr= log_norm(v);
        buf.push_back(1.0); buf.push_back(lr);
        buf.push_back(std::isfinite(lr) ? v[2] / std::exp(lr) : NAN);

        m2.draw(gen2, v);
        lr= log_norm(v);
        buf.push_back(2.0); buf.push_back(lr);
        buf.push_back(std::isfinite(lr) ? v[2] / std::exp(lr) : NAN);

        if (m3.applicable)
        {
            m3.draw(gen3, uni, v);
            lr= log_norm(v);
            buf.push_back(3.0); buf.push_back(lr);
            buf.push_back(std::isfinite(lr) ? v[2] / std::exp(lr) : NAN);
        }

        if (buf.size() > 3 * 3 * 1000)
        {
            std::fwrite(buf.data(), sizeof(double), buf.size(), fh);
            buf.clear();
        }
    }
    if (!buf.empty())
        std::fwrite(buf.data(), sizeof(double), buf.size(), fh);
    std::fclose(fh);

    // Acceptance for the one rejection-based method, to stdout as JSON.
    std::printf("{\"kappa\": %.6g, \"theta\": %.6g, \"seed\": %d, \"n\": %ld, "
                "\"m3_applicable\": %s, \"m3_n_envelope\": %.6g, \"m3_attempts\": %ld, "
                "\"m3_accepted\": %ld, \"m3_acceptance\": %.9g}\n",
                kappa, theta, seed, n, m3.applicable ? "true" : "false", m3.n,
                m3.attempts, m3.accepted,
                m3.attempts > 0 ? static_cast<double>(m3.accepted) /
                                      static_cast<double>(m3.attempts)
                                : 0.0);
    return 0;
}

// --------------------------------------------------------------------------- //
// time -- repeated batches, checksum-guarded
// --------------------------------------------------------------------------- //
#include <chrono>

static double now_seconds()
{
    using clock= std::chrono::steady_clock;
    return std::chrono::duration<double>(clock::now().time_since_epoch()).count();
}

static int run_time(double kappa, double theta, int seed, long n, int repeats,
                    const char *method, const char *label)
{
    std::mt19937 gen(static_cast<unsigned>(seed));
    std::uniform_real_distribution<double> uni(0.0, 1.0);

    // Distinguish the released implementation's own variants so that the cost of
    // anisotropy and of the field rotation is attributable, rather than folded
    // into a single "our method" number.
    const bool aniso= (std::strcmp(label, "aniso") == 0);
    const bool rotate= (std::strcmp(label, "rotated") == 0);
    const bool capped= (std::strcmp(label, "capped20") == 0);

    const double tperp= theta;
    const double tpar= aniso || rotate ? 2.0 * theta : theta;
    bi_kappa_distribution<double>::point_type ub= {0.0, 0.0, 1.0};
    if (rotate)
        ub= {0.3, -0.5, 0.8};
    const double cap= capped ? 20.0 : bi_kappa_distribution<double>::no_cap();

    bi_kappa_distribution<double> m1(kappa, tperp, tpar, ub, cap);
    ScaleMixtureNormals m2(kappa, theta);
    ParetoRejection m3(kappa, theta);

    if (std::strcmp(method, "pareto_rejection") == 0 && !m3.applicable)
    {
        std::printf("{\"method\": \"%s\", \"variant\": \"%s\", \"kappa\": %.6g, "
                    "\"applicable\": false, \"reason\": "
                    "\"n = kappa/2 violates 0 < n < kappa - 1/2\"}\n",
                    method, label, kappa);
        return 0;
    }

    // "cap" is emitted as JSON null when the cap is off: `inf` is not valid JSON
    // and would make the whole record unparseable.
    char capbuf[32];
    if (std::isfinite(cap))
        std::snprintf(capbuf, sizeof(capbuf), "%.6g", cap);
    else
        std::snprintf(capbuf, sizeof(capbuf), "null");

    std::printf("{\"method\": \"%s\", \"variant\": \"%s\", \"kappa\": %.6g, "
                "\"theta_perp\": %.6g, \"theta_par\": %.6g, \"cap\": %s, "
                "\"seed\": %d, \"n_per_batch\": %ld, \"repeats\": %d, "
                "\"applicable\": true, \"compiler\": \"%s\", \"stdlib\": \"%s\", "
                "\"arch\": \"%s\", \"batches_ns_per_sample\": [",
                method, label, kappa, tperp, tpar, capbuf, seed, n, repeats,
                compiler_id(), stdlib_id(), arch_id());

    // One untimed warm-up batch: first-call allocation and page faults are not
    // steady-state per-sample cost.
    double sink= 0.0;
    for (long i= 0; i < n / 10 + 1; ++i)
    {
        double v[3];
        if (std::strcmp(method, "gamma_ratio_spherical") == 0)
        { const auto a= m1(gen); v[0]=a[0]; v[1]=a[1]; v[2]=a[2]; }
        else if (std::strcmp(method, "scale_mixture_normals") == 0)
            m2.draw(gen, v);
        else
            m3.draw(gen, uni, v);
        sink+= v[0] + v[1] + v[2];
    }

    for (int r= 0; r < repeats; ++r)
    {
        const double t0= now_seconds();
        for (long i= 0; i < n; ++i)
        {
            double v[3];
            if (std::strcmp(method, "gamma_ratio_spherical") == 0)
            { const auto a= m1(gen); v[0]=a[0]; v[1]=a[1]; v[2]=a[2]; }
            else if (std::strcmp(method, "scale_mixture_normals") == 0)
                m2.draw(gen, v);
            else
                m3.draw(gen, uni, v);
            // Checksum: keeps the loop observable.  Adding the components is
            // cheap next to a Gamma draw and is applied identically to all
            // methods, so it cannot bias the comparison.
            sink+= v[0] + v[1] + v[2];
        }
        const double t1= now_seconds();
        std::printf("%s%.6f", r ? ", " : "", (t1 - t0) * 1e9 / static_cast<double>(n));
    }

    std::printf("], \"acceptance\": %.9g, \"n_attempts\": %ld, \"n_accepted\": %ld, "
                "\"checksum_finite\": %s}\n",
                (std::strcmp(method, "pareto_rejection") == 0 && m3.attempts > 0)
                    ? static_cast<double>(m3.accepted) / static_cast<double>(m3.attempts)
                    : 1.0,
                m3.attempts, m3.accepted,
                std::isfinite(sink) ? "true" : "false");
    return 0;
}

// --------------------------------------------------------------------------- //
int main(int argc, char **argv)
{
    if (argc < 2)
    {
        std::fprintf(stderr,
                     "usage:\n"
                     "  %s validate <out.bin> <kappa> <theta> <seed> <n>\n"
                     "  %s time <method> <variant> <kappa> <theta> <seed> <n> <repeats>\n"
                     "    method  = gamma_ratio_spherical | scale_mixture_normals | "
                     "pareto_rejection\n"
                     "    variant = iso | aniso | rotated | capped20\n",
                     argv[0], argv[0]);
        return 2;
    }

    if (std::strcmp(argv[1], "validate") == 0 && argc == 7)
        return run_validate(argv[2], std::atof(argv[3]), std::atof(argv[4]),
                            std::atoi(argv[5]), std::atol(argv[6]));

    if (std::strcmp(argv[1], "time") == 0 && argc == 9)
        return run_time(std::atof(argv[4]), std::atof(argv[5]), std::atoi(argv[6]),
                        std::atol(argv[7]), std::atoi(argv[8]), argv[2], argv[3]);

    std::fprintf(stderr, "bad arguments\n");
    return 2;
}
