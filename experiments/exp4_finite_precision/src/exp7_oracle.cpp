// Experiment 7 - independent arbitrary-precision oracle.
//
// This program exists to disagree with the probe.  It reads the probe's audit stream, which
// contains only the *primitives* of each audited attempt - x1, y, u, cos(theta), phi, kappa
// and the precision - and re-derives every classification from scratch at 100 decimal
// digits, using its own code.  It shares no arithmetic with exp7_loaders.H; the only things
// it borrows from exp7_common.H are the binary record layout, because it has to read the
// file, and the frozen accuracy thresholds, because those are the protocol's.
//
// (cos(theta), phi) is the *paired diagnostic layer's* declared direction primitive, not a
// record of what either released sampler drew: release 1.0.0 draws its angles through
// std::uniform_real_distribution, and release 2.0.0 draws no angle at all - its direction
// comes from a rejection method on the unit disc.  What this program adjudicates is the
// radius formations on the layer's declared draw, which is what the audit stream is for.
// One consequence, bounded rather than ignored: the working-precision recomputation below
// evaluates sin(phi) and cos(phi) of one argument, and a compiler may fuse that pair into a
// combined routine whose sine differs by an ulp, so this program's working-precision
// direction can differ in the last place from the direction the probe formed when the two
// were built by different compilers.  A classification can only turn on that where a
// component sits within an ulp of the type's limit, which is inside the fully audited
// within-margin stratum and is adjudicated at 100 digits from the primitives themselves.
//
// C++14, because Boost.Multiprecision requires it.  That is also why the oracle is a
// separate program: the probe and the sampling headers stay strictly C++11, and the
// dependency cannot leak into them.
//
// Four questions per audited attempt:
//
//   1. Is the intended final velocity representable in the run's floating-point type?
//      Decided against the exact overflow threshold (2 - 2^-p) 2^emax, which is half an ulp
//      above max() and differs from it only in the 16th digit - the reason the band near it
//      is audited at all.
//   2. What does each formation actually produce in working precision?  Recomputed here
//      from the primitives.
//   3. Is a returned radius *accurate*?  Amendment 1.1.0 makes a finite return a success
//      only when its relative error against this program's 100-digit value is at or below
//      2^-(digits/2) for its type.  Experiment 6 had no such question, and so scored a
//      radius that had lost 48 % of its value as a success.
//   4. Does the terminal category the probe recorded match the one that follows?
//
// Two defects of Experiment 6's oracle are fixed here.  It recorded the audit file's path
// but not its content hash, so nothing tied an adjudication to the bytes it adjudicated;
// every summary line now carries the SHA-256 of the file it read.  And `make oracle`
// swallowed a non-zero exit, so a disagreement could not fail the target; this program
// still exits non-zero on any disagreement, and the GNUmakefile no longer hides it.
//
// Build:  see GNUmakefile target exp7_oracle.exe.

#include "exp7_common.H"

#include <boost/multiprecision/cpp_dec_float.hpp>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

using boost::multiprecision::cpp_dec_float_100;
using namespace exp7;

namespace
{

/// Exact conversion check.  A double is a dyadic rational; for every magnitude this
/// experiment produces, 100 decimal digits hold it exactly.  Round-tripping back and
/// demanding bitwise equality turns that assumption into a test.
bool exactFromDouble(double v, cpp_dec_float_100 *out)
{
    *out= cpp_dec_float_100(v);
    return static_cast<double>(*out) == v;
}

cpp_dec_float_100 logExact(const cpp_dec_float_100 &x) { return log(x); }

/// log of the exact real threshold at which a value rounds to infinity in a type with
/// `digits` bits of precision and maximum exponent `emax`.
cpp_dec_float_100 logOverflowExact(int digits, int emax)
{
    cpp_dec_float_100 two(2);
    cpp_dec_float_100 half_ulp= cpp_dec_float_100(1);
    for (int i= 0; i < digits; ++i)
        half_ulp/= two; // 2^-p
    cpp_dec_float_100 mant= two - half_ulp;
    return log(mant) + cpp_dec_float_100(emax) * log(two);
}

/// Working-precision recomputation, written independently of the probe.  The audit stream
/// covers the isotropic, unrotated configuration only, which is what P2 and P3 sample.
template <typename T> struct WorkingOutcome
{
    bool qf_finite;
    bool legacy_finite;
    bool candidate_finite;
    bool x2_zero;
    bool x2_subnormal;
    T log_x2;
    T r_legacy;
    T r_qf;
    T log_r;
};

template <typename T>
WorkingOutcome<T> recomputeWorking(double kappa, double x1d, double yd, double ud, double ctd,
                                   double phd)
{
    const T a= static_cast<T>(kappa) - T(0.5);
    const T x1= static_cast<T>(x1d);
    const T y= static_cast<T>(yd);
    const T u= static_cast<T>(ud);
    const T ct= static_cast<T>(ctd);
    const T ph= static_cast<T>(phd);

    const T x2= y * std::pow(u, T(1) / a);
    const T log_x2= std::log(y) + std::log(u) / a;
    const T log_x1= std::log(x1);

    const T r_qf= std::sqrt(x1 / x2);
    const T r_legacy= std::sqrt(x1) / std::sqrt(x2);
    const T log_r= T(0.5) * (log_x1 - log_x2);

    const T sk= std::sqrt(static_cast<T>(kappa));
    const T st= std::sqrt(std::max(T(0), T(1) - ct * ct));
    const T n0= st * std::cos(ph), n1= st * std::sin(ph), n2= ct;

    WorkingOutcome<T> o;
    o.x2_zero= (x2 == T(0));
    o.x2_subnormal= (x2 != T(0)) && (std::fabs(x2) < std::numeric_limits<T>::min());
    o.log_x2= log_x2;
    o.r_legacy= r_legacy;
    o.r_qf= r_qf;
    o.log_r= log_r;

    // LEGACY operation order: radius times unit direction, then the thermal scaling.
    const T qf0= sk * (r_qf * n0), qf1= sk * (r_qf * n1), qf2= sk * (r_qf * n2);
    const T sp0= sk * (r_legacy * n0), sp1= sk * (r_legacy * n1), sp2= sk * (r_legacy * n2);
    o.qf_finite= std::isfinite(qf0) && std::isfinite(qf1) && std::isfinite(qf2);
    o.legacy_finite= std::isfinite(sp0) && std::isfinite(sp1) && std::isfinite(sp2);

    // CANDIDATE: build g first, then materialize each component and decide from the
    // component, which is the predicate the released header applies.  Written out here
    // rather than included from exp7_loaders.H: the oracle must not share code with the
    // thing it audits, and an independent transcription of the same rule is the point.
    //
    // The released header used to compare `log R` and `log R + log|g_j|` against
    // `log(max())`.  That rule is not the arithmetic it predicts -- `log(max())` is a
    // rounded value and `exp` of it need not be finite -- so an oracle carrying it would
    // model code that no longer exists, and every draw within a rounding step of the
    // boundary would be adjudicated against the wrong predicate.  Those draws are exactly
    // the ones the margin rule of PROTOCOL.md 2.3 audits at rate 1.
    const T max_finite= std::numeric_limits<T>::max();
    bool ok= (log_r == log_r);
    if (ok)
    {
        const T gs[3]= {sk * n0, sk * n1, sk * n2};
        const T r_try= std::exp(log_r);
        const bool radius_rep= (r_try <= max_finite);
        for (int j= 0; j < 3; ++j)
        {
            if (gs[j] == T(0))
                continue;
            const T mag= radius_rep ? r_try * std::fabs(gs[j])
                                    : std::exp(log_r + std::log(std::fabs(gs[j])));
            if (!(mag <= max_finite))
                ok= false;
        }
    }
    o.candidate_finite= ok;
    return o;
}

/// Method indices, restated here rather than included from exp7_loaders.H.  The oracle must
/// not share code with the thing it audits; what it does share is the record layout, and a
/// three-value enum is part of that layout.  A mismatch would show up immediately as a
/// wholesale disagreement rather than as a subtle one.
enum
{
    kOracleMethodLegacy= 0,
    kOracleMethodCandidate= 1,
    kOracleMethodQF= 2
};

/// The terminal classifier, re-derived.
int categoryFor(int method, bool finite, bool representable, double rel_err,
                double max_rel_err, bool x2_zero)
{
    if (!representable)
        return finite ? kCatOverflowReturnedFinite : kCatHonestOverflow;
    if (finite)
    {
        if (rel_err == rel_err && rel_err > max_rel_err)
            return kCatFiniteButWrong;
        return kCatFinite;
    }
    if (method == kOracleMethodCandidate)
        return kCatLogPrimFail;
    if (x2_zero)
        return kCatDenomZero;
    return (method == kOracleMethodQF) ? kCatQuotientLoss : kCatLegacyLoss;
}

/// |r/R - 1| where r is a working-precision radius and R the exact one, both given as logs.
/// Evaluated in the log domain so that nothing overflows.
double relErrFromLogs(double log_r_working, const cpp_dec_float_100 &log_r_exact)
{
    if (!(log_r_working == log_r_working) || std::isinf(log_r_working))
        return std::numeric_limits<double>::quiet_NaN();
    const double d= log_r_working - static_cast<double>(log_r_exact);
    if (d > 1.0 || d < -1.0)
        return 1.0;
    return std::fabs(std::expm1(d));
}

/// Accuracy, tracked per precision and never pooled across them.  Amendment 1.1.0 exists
/// because Experiment 6's headline contrast put a float worst case next to a float figure
/// without labelling either, while its double worst cases were six orders of magnitude
/// apart from both.  Subnormal denominators are tracked separately inside each precision,
/// because that is where the legacy radius is degraded on exactly the draws it survives.
struct Acc
{
    long long n= 0, n_subnormal= 0;
    double max_log_x2_err= 0.0, sum_log_x2_err= 0.0;
    long long n_cmp= 0, n_cmp_sub= 0;
    double max_rel_legacy= 0.0, sum_rel_legacy= 0.0;
    double max_rel_candidate= 0.0, sum_rel_candidate= 0.0;
    double max_rel_legacy_sub= 0.0, max_rel_candidate_sub= 0.0;
    long long fbw_legacy= 0, fbw_candidate= 0;
    long long returned_finite_legacy= 0, returned_finite_candidate= 0;
};

} // namespace

int main(int argc, char **argv)
{
    if (argc < 2)
    {
        std::fprintf(stderr,
                     "usage: %s <audit.bin> [--disagreements out.jsonl] [--progress]\n"
                     "       [--shard I --shards N]\n",
                     argv[0]);
        return 2;
    }
    const char *path= argv[1];
    const char *disagree_path= 0;
    bool progress= false;
    // Sharding exists only to use more than one core: at about 1.1 ms a record, adjudicating
    // the frozen matrix takes hours on one.  A shard reads the whole file and adjudicates the
    // records congruent to `shard` modulo `shards`, so the union of all shards is exactly the
    // unsharded run, every record is adjudicated once, and the record index is a property of
    // the file rather than of how the work was divided.
    long long shard= 0, shards= 1;
    for (int i= 2; i < argc; ++i)
    {
        if (std::strcmp(argv[i], "--disagreements") == 0 && i + 1 < argc)
            disagree_path= argv[++i];
        else if (std::strcmp(argv[i], "--progress") == 0)
            progress= true;
        else if (std::strcmp(argv[i], "--shard") == 0 && i + 1 < argc)
            shard= std::atoll(argv[++i]);
        else if (std::strcmp(argv[i], "--shards") == 0 && i + 1 < argc)
            shards= std::atoll(argv[++i]);
    }
    if (shards < 1 || shard < 0 || shard >= shards)
    {
        std::fprintf(stderr, "exp7_oracle: --shard must satisfy 0 <= shard < shards\n");
        return 2;
    }

    // The hash of the bytes being adjudicated.  Experiment 6 recorded the path and not the
    // content, so nothing tied an adjudication to a file that might since have changed.
    const std::string file_sha= Sha256::file(path);

    FILE *f= std::fopen(path, "rb");
    if (!f)
    {
        std::fprintf(stderr, "exp7_oracle: cannot open %s\n", path);
        return 1;
    }
    FileHeader hdr;
    if (std::fread(&hdr, sizeof(hdr), 1, f) != 1 || std::memcmp(hdr.magic, "EXP7REC", 7) != 0)
    {
        std::fprintf(stderr, "exp7_oracle: %s is not an exp7 record file\n", path);
        std::fclose(f);
        return 1;
    }
    if (hdr.record_kind != kKindAudit || hdr.record_size != sizeof(AuditRecord) ||
        hdr.schema != kAuditSchemaVersion)
    {
        std::fprintf(stderr,
                     "exp7_oracle: %s has kind %u / size %u / schema %u, expected %u / %zu / %u\n",
                     path, hdr.record_kind, hdr.record_size, hdr.schema, kKindAudit,
                     sizeof(AuditRecord), kAuditSchemaVersion);
        std::fclose(f);
        return 1;
    }

    FILE *dj= disagree_path ? std::fopen(disagree_path, "wb") : 0;
    if (dj)
    {
        // A header record, written whether or not anything disagrees.  An empty file is
        // otherwise indistinguishable from an oracle that never ran, which is what
        // `require_header_record_even_when_zero_disagreements` in the protocol is about.
        std::fprintf(dj,
                     "{\"kind\":\"header\",\"tool\":\"exp7_oracle\",\"file\":\"%s\","
                     "\"file_sha256\":\"%s\",\"n_records_declared\":%lld,"
                     "\"protocol_sha256\":\"%s\",\"max_rel_error_double\":%.17g,"
                     "\"max_rel_error_float\":%.17g}\n",
                     path, file_sha.c_str(), static_cast<long long>(hdr.n_records),
                     protocol::kProtocolSha256, protocol::kMaxRelErrorDouble,
                     protocol::kMaxRelErrorFloat);
    }

    long long n= 0, record_index= 0, disagreements= 0, conversion_failures= 0, near_limit= 0;
    long long cat_checked= 0;
    double worst_margin= 0.0;
    Acc acc[2]; // [0] double, [1] float

    const cpp_dec_float_100 log_ovf_double= logOverflowExact(53, 1023);
    const cpp_dec_float_100 log_ovf_float= logOverflowExact(24, 127);

    std::vector<AuditRecord> buf(65536);
    for (;;)
    {
        const size_t got= std::fread(&buf[0], sizeof(AuditRecord), buf.size(), f);
        if (got == 0)
            break;
        for (size_t i= 0; i < got; ++i)
        {
            const AuditRecord &a= buf[i];
            ++record_index;
            if (shards > 1 && (record_index - 1) % shards != shard)
                continue;
            ++n;
            if (progress && (n % 100000) == 0)
                std::fprintf(stderr, "  exp7_oracle: %lld records, %lld disagreements\n", n,
                             disagreements);
            const bool is_float= a.precision_is_float != 0;
            const double max_rel= is_float ? protocol::kMaxRelErrorFloat
                                           : protocol::kMaxRelErrorDouble;

            cpp_dec_float_100 x1, y, u, ct, ph, kap;
            bool ok= exactFromDouble(a.x1, &x1) && exactFromDouble(a.y, &y) &&
                     exactFromDouble(a.u, &u) && exactFromDouble(a.cos_theta, &ct) &&
                     exactFromDouble(a.phi, &ph) && exactFromDouble(a.kappa, &kap);
            if (!ok)
            {
                ++conversion_failures;
                continue;
            }

            // The shape the run actually used, formed in the working precision.
            const double a_work= is_float
                                     ? static_cast<double>(static_cast<float>(a.kappa) - 0.5f)
                                     : (a.kappa - 0.5);
            cpp_dec_float_100 alpha;
            if (!exactFromDouble(a_work, &alpha))
            {
                ++conversion_failures;
                continue;
            }

            const cpp_dec_float_100 log_x1= logExact(x1);
            const cpp_dec_float_100 log_x2= logExact(y) + logExact(u) / alpha;
            const cpp_dec_float_100 log_r= (log_x1 - log_x2) / cpp_dec_float_100(2);

            const cpp_dec_float_100 one(1);
            const cpp_dec_float_100 st= sqrt(one - ct * ct);
            const cpp_dec_float_100 sk= sqrt(kap);
            const cpp_dec_float_100 g0= sk * st * cos(ph);
            const cpp_dec_float_100 g1= sk * st * sin(ph);
            const cpp_dec_float_100 g2= sk * ct;

            const cpp_dec_float_100 &log_ovf= is_float ? log_ovf_float : log_ovf_double;
            cpp_dec_float_100 max_comp= log_r + logExact(abs(g0));
            {
                const cpp_dec_float_100 c1= log_r + logExact(abs(g1));
                const cpp_dec_float_100 c2= log_r + logExact(abs(g2));
                if (c1 > max_comp) max_comp= c1;
                if (c2 > max_comp) max_comp= c2;
            }
            const bool representable= (max_comp < log_ovf);
            const double margin= static_cast<double>(max_comp - log_ovf);
            if (std::fabs(margin) <= oracleBandLogUnits())
                ++near_limit;
            if (std::fabs(margin) > std::fabs(worst_margin))
                worst_margin= margin;

            // Independent working-precision recomputation.
            int cl_legacy, cl_candidate, cl_qf;
            Acc &A= acc[is_float ? 1 : 0];
            {
                double wlog_x2, wlog_r, wlog_r_legacy, wlog_r_qf;
                bool x2_zero, x2_sub, qf_fin, lg_fin, cd_fin;
                if (is_float)
                {
                    const WorkingOutcome<float> wo=
                        recomputeWorking<float>(a.kappa, a.x1, a.y, a.u, a.cos_theta, a.phi);
                    wlog_x2= wo.log_x2;
                    wlog_r= wo.log_r;
                    wlog_r_legacy= (std::isfinite(wo.r_legacy) && wo.r_legacy > 0.0f)
                                       ? std::log(static_cast<double>(wo.r_legacy))
                                       : std::numeric_limits<double>::quiet_NaN();
                    wlog_r_qf= (std::isfinite(wo.r_qf) && wo.r_qf > 0.0f)
                                   ? std::log(static_cast<double>(wo.r_qf))
                                   : std::numeric_limits<double>::quiet_NaN();
                    x2_zero= wo.x2_zero;
                    x2_sub= wo.x2_subnormal;
                    qf_fin= wo.qf_finite;
                    lg_fin= wo.legacy_finite;
                    cd_fin= wo.candidate_finite;
                }
                else
                {
                    const WorkingOutcome<double> wo=
                        recomputeWorking<double>(a.kappa, a.x1, a.y, a.u, a.cos_theta, a.phi);
                    wlog_x2= wo.log_x2;
                    wlog_r= wo.log_r;
                    wlog_r_legacy= (std::isfinite(wo.r_legacy) && wo.r_legacy > 0.0)
                                       ? std::log(wo.r_legacy)
                                       : std::numeric_limits<double>::quiet_NaN();
                    wlog_r_qf= (std::isfinite(wo.r_qf) && wo.r_qf > 0.0)
                                   ? std::log(wo.r_qf)
                                   : std::numeric_limits<double>::quiet_NaN();
                    x2_zero= wo.x2_zero;
                    x2_sub= wo.x2_subnormal;
                    qf_fin= wo.qf_finite;
                    lg_fin= wo.legacy_finite;
                    cd_fin= wo.candidate_finite;
                }

                const double rel_legacy= relErrFromLogs(wlog_r_legacy, log_r);
                const double rel_candidate= relErrFromLogs(wlog_r, log_r);
                const double rel_qf= relErrFromLogs(wlog_r_qf, log_r);

                cl_legacy= categoryFor(kOracleMethodLegacy, lg_fin, representable, rel_legacy,
                                       max_rel, x2_zero);
                cl_candidate= categoryFor(kOracleMethodCandidate, cd_fin, representable,
                                          rel_candidate, max_rel, x2_zero);
                cl_qf= categoryFor(kOracleMethodQF, qf_fin, representable, rel_qf, max_rel, x2_zero);

                ++A.n;
                if (x2_sub)
                    ++A.n_subnormal;
                const double lerr= std::fabs(wlog_x2 - static_cast<double>(log_x2));
                A.sum_log_x2_err+= lerr;
                if (lerr > A.max_log_x2_err)
                    A.max_log_x2_err= lerr;

                if (lg_fin)
                    ++A.returned_finite_legacy;
                if (cd_fin)
                    ++A.returned_finite_candidate;
                if (cl_legacy == kCatFiniteButWrong)
                    ++A.fbw_legacy;
                if (cl_candidate == kCatFiniteButWrong)
                    ++A.fbw_candidate;

                if (rel_legacy == rel_legacy && rel_candidate == rel_candidate &&
                    rel_legacy < 1.0 && rel_candidate < 1.0)
                {
                    ++A.n_cmp;
                    A.sum_rel_legacy+= rel_legacy;
                    A.sum_rel_candidate+= rel_candidate;
                    if (rel_legacy > A.max_rel_legacy) A.max_rel_legacy= rel_legacy;
                    if (rel_candidate > A.max_rel_candidate) A.max_rel_candidate= rel_candidate;
                    if (x2_sub)
                    {
                        ++A.n_cmp_sub;
                        if (rel_legacy > A.max_rel_legacy_sub)
                            A.max_rel_legacy_sub= rel_legacy;
                        if (rel_candidate > A.max_rel_candidate_sub)
                            A.max_rel_candidate_sub= rel_candidate;
                    }
                }
            }

            ++cat_checked;
            const bool same= (cl_legacy == static_cast<int>(a.cat_legacy)) &&
                             (cl_candidate == static_cast<int>(a.cat_candidate)) &&
                             (cl_qf == static_cast<int>(a.cat_qf));
            if (!same)
            {
                ++disagreements;
                if (dj)
                    std::fprintf(dj,
                                 "{\"kind\":\"disagreement\",\"file\":\"%s\","
                                 "\"file_sha256\":\"%s\",\"kappa\":%.17g,\"precision\":\"%s\","
                                 "\"x1\":%.17g,\"y\":%.17g,\"u\":%.17g,\"cos_theta\":%.17g,"
                                 "\"phi\":%.17g,\"probe\":[%u,%u,%u],\"oracle\":[%d,%d,%d],"
                                 "\"margin_log_units\":%.17g}\n",
                                 path, file_sha.c_str(), a.kappa,
                                 is_float ? "float" : "double", a.x1, a.y, a.u, a.cos_theta,
                                 a.phi, a.cat_legacy, a.cat_candidate, a.cat_qf, cl_legacy,
                                 cl_candidate, cl_qf, margin);
            }
        }
        if (got < buf.size())
            break;
    }
    std::fclose(f);

    if (dj)
    {
        std::fprintf(dj,
                     "{\"kind\":\"footer\",\"file\":\"%s\",\"file_sha256\":\"%s\","
                     "\"n_records\":%lld,\"disagreements\":%lld,\"conversion_failures\":%lld}\n",
                     path, file_sha.c_str(), n, disagreements, conversion_failures);
        std::fclose(dj);
    }

    for (int pi= 0; pi < 2; ++pi)
    {
        const Acc &A= acc[pi];
        if (A.n == 0)
            continue;
        std::printf(
            "{\"tool\":\"exp7_oracle\","
            "\"oracle\":\"boost::multiprecision::cpp_dec_float_100\","
            "\"file\":\"%s\",\"file_sha256\":\"%s\",\"protocol_sha256\":\"%s\","
            "\"precision\":\"%s\",\"n_records\":%lld,\"n_records_total\":%lld,"
            "\"n_classified\":%lld,\"disagreements\":%lld,\"conversion_failures\":%lld,"
            "\"near_limit_records\":%lld,\"worst_overflow_margin_log_units\":%.6g,"
            "\"n_subnormal_denominator\":%lld,\"max_abs_log_x2_error\":%.6g,"
            "\"mean_abs_log_x2_error\":%.6g,\"n_radius_compared\":%lld,"
            "\"n_radius_compared_subnormal\":%lld,\"max_rel_error_threshold\":%.17g,"
            "\"returned_finite_legacy\":%lld,\"returned_finite_candidate\":%lld,"
            "\"finite_but_wrong_legacy\":%lld,\"finite_but_wrong_candidate\":%lld,"
            "\"max_rel_err_legacy\":%.6g,\"mean_rel_err_legacy\":%.6g,"
            "\"max_rel_err_candidate\":%.6g,\"mean_rel_err_candidate\":%.6g,"
            "\"max_rel_err_legacy_subnormal_denominator\":%.6g,"
            "\"max_rel_err_candidate_subnormal_denominator\":%.6g}\n",
            path, file_sha.c_str(), protocol::kProtocolSha256, pi ? "float" : "double", A.n,
            n, cat_checked, disagreements, conversion_failures, near_limit, worst_margin,
            A.n_subnormal, A.max_log_x2_err, A.n ? A.sum_log_x2_err / A.n : 0.0, A.n_cmp,
            A.n_cmp_sub, pi ? protocol::kMaxRelErrorFloat : protocol::kMaxRelErrorDouble,
            A.returned_finite_legacy, A.returned_finite_candidate, A.fbw_legacy,
            A.fbw_candidate, A.max_rel_legacy,
            A.n_cmp ? A.sum_rel_legacy / A.n_cmp : 0.0, A.max_rel_candidate,
            A.n_cmp ? A.sum_rel_candidate / A.n_cmp : 0.0, A.max_rel_legacy_sub,
            A.max_rel_candidate_sub);
    }

    if (n == 0)
    {
        // An audit file with no records adjudicates nothing.  Saying so is the honest
        // answer; a silent success here is how an absent audit passes for a clean one.
        std::printf("{\"tool\":\"exp7_oracle\",\"file\":\"%s\",\"file_sha256\":\"%s\","
                    "\"n_records\":0,\"note\":\"no audit records in this file\"}\n",
                    path, file_sha.c_str());
    }

    return (disagreements == 0 && conversion_failures == 0) ? 0 : 1;
}
