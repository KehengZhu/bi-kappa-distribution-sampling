// Experiment 6 - independent arbitrary-precision oracle.
//
// This program exists to disagree with the probe.  It reads the probe's audit stream,
// which contains only the *primitives* of each audited attempt - x1, y, u, cos(theta),
// phi, kappa and the precision - and re-derives every classification from scratch at 100
// decimal digits, using its own code.  It shares no arithmetic with exp6_loaders.H; the
// only thing it borrows from exp6_common.H is the binary record layout, because it has to
// read the file.
//
// C++14, because Boost.Multiprecision requires it.  That is also why the oracle is a
// separate program: the probe and the sampling headers stay strictly C++11, and the
// dependency cannot leak into them.
//
// Three questions per audited attempt:
//
//   1. Is the intended final velocity representable in the run's floating-point type?
//      Decided against the exact overflow threshold (2 - 2^-p) 2^emax, which is half an
//      ulp above max() and differs from it only in the 16th digit - the reason the band
//      near it is audited at all.
//   2. What does each formation actually produce in working precision?  Recomputed here
//      from the primitives.
//   3. Does the terminal category the probe recorded match the one that follows?
//
// Build:  see GNUmakefile target exp6_oracle.exe.

#include "exp6_common.H"

#include <boost/multiprecision/cpp_dec_float.hpp>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

using boost::multiprecision::cpp_dec_float_100;
using namespace exp6;

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
        half_ulp/= two;                      // 2^-p
    cpp_dec_float_100 mant= two - half_ulp;  // 2 - 2^-p
    return log(mant) + cpp_dec_float_100(emax) * log(two);
}

/// Working-precision recomputation, written independently of the probe.
template <typename T> struct WorkingOutcome
{
    bool qf_finite;
    bool split_finite;
    bool log_finite;
    bool x2_zero;
    bool x2_subnormal;
    T log_x2;
    T r_split;   ///< the radius the split form materialized (may be non-finite)
    T log_r;     ///< the radius the log path carried
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
    const T r_split= std::sqrt(x1) / std::sqrt(x2);
    const T log_r= T(0.5) * (log_x1 - log_x2);

    const T sk= std::sqrt(static_cast<T>(kappa));
    const T st= std::sqrt(std::max(T(0), T(1) - ct * ct));
    const T n0= st * std::cos(ph), n1= st * std::sin(ph), n2= ct;

    WorkingOutcome<T> o;
    o.x2_zero= (x2 == T(0));
    o.x2_subnormal= (x2 != T(0)) && (std::fabs(x2) < std::numeric_limits<T>::min());
    o.log_x2= log_x2;
    o.r_split= r_split;
    o.log_r= log_r;

    // Released operation order: radius times unit direction, then the thermal scaling.
    // Isotropic, unrotated: the audit stream covers exactly that configuration.
    const T qf0= sk * (r_qf * n0), qf1= sk * (r_qf * n1), qf2= sk * (r_qf * n2);
    const T sp0= sk * (r_split * n0), sp1= sk * (r_split * n1), sp2= sk * (r_split * n2);
    o.qf_finite= std::isfinite(qf0) && std::isfinite(qf1) && std::isfinite(qf2);
    o.split_finite= std::isfinite(sp0) && std::isfinite(sp1) && std::isfinite(sp2);

    // Log path: decide before exponentiating.
    const int digits= std::numeric_limits<T>::digits;
    const int emax= std::numeric_limits<T>::max_exponent - 1;
    const double kLn2= 0.69314718055994530942;
    const T log_ovf= static_cast<T>(emax * kLn2 + std::log(2.0 - std::ldexp(1.0, -digits)));
    bool log_ok= (log_r == log_r);
    if (log_ok)
    {
        const T gs[3]= {sk * n0, sk * n1, sk * n2};
        for (int j= 0; j < 3; ++j)
            if (gs[j] != T(0) && log_r + std::log(std::fabs(gs[j])) >= log_ovf)
                log_ok= false;
    }
    o.log_finite= log_ok;
    return o;
}

int categoryFor(int method, bool finite, bool representable, bool x2_zero)
{
    if (!representable)
        return finite ? kCatFinite : kCatHonestOverflow;
    if (finite)
        return kCatFinite;
    if (method == 2)
        return kCatLogPrimFail;
    if (x2_zero)
        return kCatDenomZero;
    return (method == 0) ? kCatQuotientLoss : kCatSplitLoss;
}

} // namespace

int main(int argc, char **argv)
{
    if (argc < 2)
    {
        std::fprintf(stderr, "usage: %s <audit.bin> [--disagreements out.jsonl]\n", argv[0]);
        return 2;
    }
    const char *path= argv[1];
    const char *disagree_path= 0;
    for (int i= 2; i + 1 < argc; ++i)
        if (std::strcmp(argv[i], "--disagreements") == 0)
            disagree_path= argv[i + 1];

    FILE *f= std::fopen(path, "rb");
    if (!f)
    {
        std::fprintf(stderr, "exp6_oracle: cannot open %s\n", path);
        return 1;
    }
    FileHeader hdr;
    if (std::fread(&hdr, sizeof(hdr), 1, f) != 1 || std::memcmp(hdr.magic, "EXP6REC", 7) != 0)
    {
        std::fprintf(stderr, "exp6_oracle: %s is not an exp6 record file\n", path);
        return 1;
    }
    if (hdr.record_kind != 4 || hdr.record_size != sizeof(AuditRecord) ||
        hdr.schema != kAuditSchemaVersion)
    {
        std::fprintf(stderr,
                     "exp6_oracle: %s has kind %u / size %u / schema %u, expected 4 / %zu / %u\n",
                     path, hdr.record_kind, hdr.record_size, hdr.schema, sizeof(AuditRecord),
                     kAuditSchemaVersion);
        return 1;
    }

    FILE *dj= disagree_path ? std::fopen(disagree_path, "wb") : 0;

    long long n= 0, disagreements= 0, conversion_failures= 0, near_limit= 0;
    long long cat_checked= 0;
    double worst_margin= 0.0;

    // Accuracy is tracked separately for the two precisions, and separately for draws whose
    // denominator went subnormal.  Pooling them hides the point: a subnormal denominator has
    // already lost most of its significand, so the split form's radius is degraded on exactly
    // the draws it barely survives, while the log path never materializes the denominator.
    struct Acc
    {
        long long n= 0, n_subnormal= 0;
        double max_log_x2_err= 0.0, sum_log_x2_err= 0.0;
        long long n_cmp= 0, n_cmp_sub= 0;
        double max_rel_split= 0.0, sum_rel_split= 0.0;
        double max_rel_log= 0.0, sum_rel_log= 0.0;
        double max_rel_split_sub= 0.0, max_rel_log_sub= 0.0;
    };
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
            ++n;
            const bool is_float= a.precision_is_float != 0;

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
            if (std::fabs(margin) <= kOracleBandLogUnits)
                ++near_limit;
            if (std::fabs(margin) > std::fabs(worst_margin))
                worst_margin= margin;

            // Independent working-precision recomputation.
            int cq, cs, cl;
            Acc &A= acc[is_float ? 1 : 0];
            {
                double wlog_x2, wr_split, wlog_r;
                bool x2_zero, x2_sub, qf_fin, sp_fin, lg_fin;
                if (is_float)
                {
                    const WorkingOutcome<float> wo=
                        recomputeWorking<float>(a.kappa, a.x1, a.y, a.u, a.cos_theta, a.phi);
                    wlog_x2= wo.log_x2; wr_split= wo.r_split; wlog_r= wo.log_r;
                    x2_zero= wo.x2_zero; x2_sub= wo.x2_subnormal;
                    qf_fin= wo.qf_finite; sp_fin= wo.split_finite; lg_fin= wo.log_finite;
                }
                else
                {
                    const WorkingOutcome<double> wo=
                        recomputeWorking<double>(a.kappa, a.x1, a.y, a.u, a.cos_theta, a.phi);
                    wlog_x2= wo.log_x2; wr_split= wo.r_split; wlog_r= wo.log_r;
                    x2_zero= wo.x2_zero; x2_sub= wo.x2_subnormal;
                    qf_fin= wo.qf_finite; sp_fin= wo.split_finite; lg_fin= wo.log_finite;
                }
                cq= categoryFor(0, qf_fin, representable, x2_zero);
                cs= categoryFor(1, sp_fin, representable, x2_zero);
                cl= categoryFor(2, lg_fin, representable, x2_zero);

                ++A.n;
                if (x2_sub)
                    ++A.n_subnormal;
                const double lerr= std::fabs(wlog_x2 - static_cast<double>(log_x2));
                A.sum_log_x2_err+= lerr;
                if (lerr > A.max_log_x2_err)
                    A.max_log_x2_err= lerr;

                // Relative error of each formation's radius against the exact radius, on
                // draws where both formations resolve it.  Evaluated in the log domain so
                // that no comparison overflows: |r/R - 1| = |expm1(log r - log R)|.
                if (std::isfinite(wr_split) && wr_split > 0.0 && std::isfinite(wlog_r))
                {
                    const double rel_split=
                        std::fabs(std::expm1(std::log(wr_split) - static_cast<double>(log_r)));
                    const double rel_log=
                        std::fabs(std::expm1(wlog_r - static_cast<double>(log_r)));
                    if (rel_split < 1.0 && rel_log < 1.0)
                    {
                        ++A.n_cmp;
                        A.sum_rel_split+= rel_split;
                        A.sum_rel_log+= rel_log;
                        if (rel_split > A.max_rel_split) A.max_rel_split= rel_split;
                        if (rel_log > A.max_rel_log) A.max_rel_log= rel_log;
                        if (x2_sub)
                        {
                            ++A.n_cmp_sub;
                            if (rel_split > A.max_rel_split_sub) A.max_rel_split_sub= rel_split;
                            if (rel_log > A.max_rel_log_sub) A.max_rel_log_sub= rel_log;
                        }
                    }
                }
            }

            ++cat_checked;
            const bool same= (cq == static_cast<int>(a.cat_qf)) &&
                             (cs == static_cast<int>(a.cat_split)) &&
                             (cl == static_cast<int>(a.cat_log));
            if (!same)
            {
                ++disagreements;
                if (dj)
                    std::fprintf(dj,
                                 "{\"kappa\":%.17g,\"precision\":\"%s\",\"x1\":%.17g,"
                                 "\"y\":%.17g,\"u\":%.17g,\"cos_theta\":%.17g,\"phi\":%.17g,"
                                 "\"probe\":[%u,%u,%u],\"oracle\":[%d,%d,%d],"
                                 "\"margin_log_units\":%.17g}\n",
                                 a.kappa, is_float ? "float" : "double", a.x1, a.y, a.u,
                                 a.cos_theta, a.phi, a.cat_qf, a.cat_split, a.cat_log, cq, cs,
                                 cl, margin);
            }
        }
        if (got < buf.size())
            break;
    }
    std::fclose(f);
    if (dj)
        std::fclose(dj);

    for (int pi= 0; pi < 2; ++pi)
    {
        const Acc &A= acc[pi];
        if (A.n == 0)
            continue;
        std::printf(
            "{\"tool\":\"exp6_oracle\","
            "\"oracle\":\"boost::multiprecision::cpp_dec_float_100\","
            "\"file\":\"%s\",\"precision\":\"%s\",\"n_records\":%lld,"
            "\"n_records_total\":%lld,\"n_classified\":%lld,\"disagreements\":%lld,"
            "\"conversion_failures\":%lld,\"near_limit_records\":%lld,"
            "\"worst_overflow_margin_log_units\":%.6g,\"n_subnormal_denominator\":%lld,"
            "\"max_abs_log_x2_error\":%.6g,\"mean_abs_log_x2_error\":%.6g,"
            "\"n_radius_compared\":%lld,\"n_radius_compared_subnormal\":%lld,"
            "\"max_rel_err_split\":%.6g,\"mean_rel_err_split\":%.6g,"
            "\"max_rel_err_log\":%.6g,\"mean_rel_err_log\":%.6g,"
            "\"max_rel_err_split_subnormal_denominator\":%.6g,"
            "\"max_rel_err_log_subnormal_denominator\":%.6g}\n",
            path, pi ? "float" : "double", A.n, n, cat_checked, disagreements,
            conversion_failures, near_limit, worst_margin, A.n_subnormal, A.max_log_x2_err,
            A.n ? A.sum_log_x2_err / A.n : 0.0, A.n_cmp, A.n_cmp_sub, A.max_rel_split,
            A.n_cmp ? A.sum_rel_split / A.n_cmp : 0.0, A.max_rel_log,
            A.n_cmp ? A.sum_rel_log / A.n_cmp : 0.0, A.max_rel_split_sub, A.max_rel_log_sub);
    }
    return (disagreements == 0 && conversion_failures == 0) ? 0 : 1;
}
