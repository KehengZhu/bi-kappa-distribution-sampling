// Experiment 7 - confirmatory probe.
//
// Strict C++11, header-only sampling, no dependency outside the standard library.  The
// arbitrary-precision oracle is a separate program; this one records the primitives it will
// need and never shares code with it.
//
// Phases (argv[1]) are the frozen matrix of PROTOCOL.md Sec. 4:
//
//   p1   scalar        LEGACY and CANDIDATE, full ladder, float and double, 10^6 x 5 seeds
//   p2   mechanism     paired diagnostic layer + native layer, full ladder, uncapped
//   p3   conditioning  per-attempt records: intended state vs success
//   p4   loader        complete anisotropic, rotated, capped 3-D loader, cases C0-C6
//   p5   portability   CANDIDATE only, full ladder, one tagged environment per build
//   p6   performance   timed blocks, cases B0-B4, seeds 7006-7010
//   selftest           accounting identity, replica-vs-class equality, error paths, digests
//   env                the floating-point environment the run actually executes in
//
// Every phase accepts --n, --out, --seeds, --smoke, --tag and --execution.  Smoke mode
// writes only under <out>/smoke/, including its manifest, so that it can never be swept
// into a production manifest, a checksum bundle or an oracle run.

#include "exp7_common.H"
#include "exp7_loaders.H"

#include <algorithm>
#include <cfenv>
#include <cfloat>
#include <chrono>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

using namespace exp7;

// ---------------------------------------------------------------------------
// Options
// ---------------------------------------------------------------------------
struct Options
{
    long long n;
    std::string out;
    std::string tag;
    std::string execution;
    std::string only_precision;
    std::string case_id;
    bool smoke;
    int n_seeds;

    Options()
        : n(0), out("raw"), tag("primary"), execution("native"), only_precision(""),
          case_id(""), smoke(false), n_seeds(5)
    {
    }
};

static Options parseOptions(int argc, char **argv, int start)
{
    Options o;
    for (int i= start; i < argc; ++i)
    {
        const std::string a= argv[i];
        const bool has_next= (i + 1 < argc);
        if (a == "--n" && has_next)
            o.n= std::atoll(argv[++i]);
        else if (a == "--out" && has_next)
            o.out= argv[++i];
        else if (a == "--tag" && has_next)
            o.tag= argv[++i];
        else if (a == "--execution" && has_next)
            o.execution= argv[++i];
        else if (a == "--precision" && has_next)
            o.only_precision= argv[++i];
        else if (a == "--case" && has_next)
            o.case_id= argv[++i];
        else if (a == "--seeds" && has_next)
            o.n_seeds= std::atoi(argv[++i]);
        else if (a == "--smoke")
            o.smoke= true;
        else
        {
            std::fprintf(stderr, "exp7: unknown option %s\n", a.c_str());
            std::exit(2);
        }
    }
    if (o.smoke)
    {
        o.n= (o.n > 0 && o.n < 2000) ? o.n : 2000;
        o.n_seeds= 1;
    }
    return o;
}

/// The declared production seeds, truncated only by --seeds for a development run.  No seed
/// is ever computed; the vector is the one declared in exp7_common.H.
static std::vector<unsigned> seedsFor(const Options &o)
{
    std::vector<unsigned> all= productionSeeds();
    if (o.n_seeds < static_cast<int>(all.size()))
        all.resize(static_cast<size_t>(o.n_seeds));
    return all;
}

static std::vector<unsigned> benchSeedsFor(const Options &o)
{
    std::vector<unsigned> all= performanceSeeds();
    if (o.n_seeds < static_cast<int>(all.size()))
        all.resize(static_cast<size_t>(o.n_seeds));
    return all;
}

static std::string joinPath(const std::string &dir, const std::string &name)
{
    if (dir.empty())
        return name;
    if (dir[dir.size() - 1] == '/')
        return dir + name;
    return dir + "/" + name;
}

/// Smoke output is diverted wholesale, manifest included.  Experiment 6 diverted the data
/// files but left the manifest, the checksum sweep and `make oracle` looking at the same
/// tree, and all three picked smoke files up.
static std::string phaseDir(const Options &o, const char *phase)
{
    if (o.smoke)
        return joinPath(joinPath(o.out, "smoke"), phase);
    return joinPath(o.out, phase);
}

static std::string manifestPath(const Options &o)
{
    if (o.smoke)
        return joinPath(joinPath(o.out, "smoke"), "manifest.csv");
    return joinPath(o.out, "manifest.csv");
}

// ---------------------------------------------------------------------------
// Atomic line-oriented output, registered in the manifest on finish.
// ---------------------------------------------------------------------------
class JsonlWriter
{
  public:
    explicit JsonlWriter(const std::string &path)
        : m_path(path), m_part(path + ".part"), m_f(0), m_lines(0)
    {
        m_f= std::fopen(m_part.c_str(), "wb");
        if (!m_f)
        {
            std::fprintf(stderr, "exp7: cannot open %s\n", m_part.c_str());
            std::exit(1);
        }
    }
    FILE *f() { return m_f; }
    long long lines() const { return m_lines; }
    void countLine() { ++m_lines; }

    std::string finish(const RunContext &rc, const char *phase)
    {
        std::fclose(m_f);
        m_f= 0;
        if (std::rename(m_part.c_str(), m_path.c_str()) != 0)
        {
            std::fprintf(stderr, "exp7: cannot rename %s\n", m_part.c_str());
            std::exit(1);
        }
        const long long bytes= verifyJsonlContract(m_path, m_lines);
        manifestAppend(rc, phase, "summary", "all", "all", "all",
                       std::numeric_limits<double>::quiet_NaN(), -1, m_path, "jsonl", m_lines,
                       -1, -1, -1, bytes, Sha256::file(m_path));
        return m_path;
    }

  private:
    std::string m_path, m_part;
    FILE *m_f;
    long long m_lines;
};

// ---------------------------------------------------------------------------
// Variate accounting.  A generator that counts is the only honest way to report "random
// variates consumed", which is a cost the protocol requires and which differs between the
// two methods.  Its min() and max() mirror std::mt19937 exactly, so wrapping the engine
// does not change a single variate.
// ---------------------------------------------------------------------------
class CountingEngine
{
  public:
    typedef std::mt19937::result_type result_type;
    explicit CountingEngine(unsigned s) : m_g(s), m_calls(0) {}
    static constexpr result_type min() { return std::mt19937::min(); }
    static constexpr result_type max() { return std::mt19937::max(); }
    result_type operator()()
    {
        ++m_calls;
        return m_g();
    }
    long long calls() const { return m_calls; }

  private:
    std::mt19937 m_g;
    long long m_calls;
};

// ---------------------------------------------------------------------------
// Shared per-configuration emission
// ---------------------------------------------------------------------------
static void emitConfigFields(FILE *f, const RunContext &rc, const char *phase,
                             const char *layer, const char *method, const char *precision,
                             double kappa, double a, unsigned seed, long long n)
{
    std::fprintf(f,
                 "{\"phase\":\"%s\",\"layer\":\"%s\",\"method\":\"%s\",\"tag\":\"%s\","
                 "\"execution\":\"%s\",\"protocol_sha256\":\"%s\",",
                 phase, layer, method, rc.tag.c_str(), rc.execution.c_str(),
                 protocol::kProtocolSha256);
    emitEnvFields(f, precision);
    std::fprintf(f, ",\"kappa\":%.17g,\"shape_a\":%.17g,\"seed\":%u,\"n_attempted\":%lld",
                 kappa, a, seed, n);
}

static void emitCategories(FILE *f, const MethodCounters &mc)
{
    for (int i= 0; i < kNumCategories; ++i)
        std::fprintf(f, ",\"cat_%s\":%lld", categoryName(i), mc.cat[i]);
    std::fprintf(f, ",\"n_finite\":%lld,\"n_avoidable\":%lld,\"n_honest\":%lld",
                 mc.cat[kCatFinite], mc.avoidable(), mc.honest());
    std::fprintf(f, ",\"nonfinite_output\":%lld", mc.nonfinite_output);
    std::fprintf(f, ",\"n_rel_err\":%lld,", mc.n_rel_err);
    emitDouble(f, "mean_rel_err_radius", mc.n_rel_err ? mc.sum_rel_err / mc.n_rel_err : 0.0);
    std::fprintf(f, ",");
    emitDouble(f, "max_rel_err_radius", mc.max_rel_err);
    std::fprintf(f, ",");
    emitDouble(f, "max_finite_log_component", mc.max_finite_log_component);
}

static void emitDigest(FILE *f, const char *prefix, const StreamDigest &d)
{
    std::fprintf(f, ",\"%s_sha256\":\"%s\",\"%s_fnv1a64\":\"%s\",\"%s_bytes\":%llu", prefix,
                 d.sha256().c_str(), prefix, d.fnv1a().c_str(), prefix,
                 static_cast<unsigned long long>(d.byteCount()));
}

static void emitAuditStrata(FILE *f, const AttemptCounters &ac)
{
    for (int i= 1; i < kNumAuditStrata; ++i)
        std::fprintf(f, ",\"audit_%s_total\":%lld,\"audit_%s_audited\":%lld,"
                        "\"audit_%s_rate\":%.17g",
                     auditStratumName(i), ac.stratum_total[i], auditStratumName(i),
                     ac.stratum_audited[i], auditStratumName(i), auditStratumRate(i));
    std::fprintf(f,
                 ",\"audit_margin_log_units\":%.17g,\"audit_margin_assertions\":%lld,"
                 "\"audit_margin_assertion_failures\":%lld,",
                 oracleBandLogUnits(), ac.margin_assertions, ac.margin_assertion_failures);
    emitDouble(f, "audit_min_unambiguous_margin", ac.min_unambiguous_margin);
}

/// Fold every counter of a configuration into the digest, so that "bitwise equality of the
/// returned vectors and of every counter" is a single comparison rather than a list of them.
static void digestCounters(StreamDigest *d, const MethodCounters &mc, const VariateCounts &vc,
                           long long n)
{
    d->counter("n_attempted", n);
    for (int i= 0; i < kNumCategories; ++i)
        d->counter(categoryName(i), mc.cat[i]);
    d->counter("nonfinite_output", mc.nonfinite_output);
    d->counter("gamma_variates", vc.gamma);
    d->counter("uniform_variates", vc.uniform);
    d->counter("engine_calls", vc.engine_calls);
}

// ===========================================================================
// The two released classes, run end to end.  This is what a user receives, and it is what
// the portability claim is about: the digest below is over the bytes the class returned.
// ===========================================================================
template <typename T> struct ClassRun
{
    long long n_finite;
    long long n_nonfinite;
    long long n_thrown;
    unsigned long long n_attempts_reported;  ///< CANDIDATE only; 0 for LEGACY
    unsigned long long n_nonfinite_reported; ///< CANDIDATE only; 0 for LEGACY
    double seconds;
    std::string sha256;
    std::string fnv1a;
    unsigned long long digest_bytes;
};

template <typename T>
static ClassRun<T> runLegacyClass(double kappa, double theta_perp, double theta_par,
                                  const double ub[3], double cap, unsigned seed, long long n,
                                  RecordWriter<LoaderRecord> *dump)
{
    typedef typename bikappa_v1::bi_kappa_distribution<T>::point_type P;
    P ubt;
    ubt[0]= static_cast<T>(ub[0]);
    ubt[1]= static_cast<T>(ub[1]);
    ubt[2]= static_cast<T>(ub[2]);
    const T capT= (cap == std::numeric_limits<double>::infinity())
                      ? bikappa_v1::bi_kappa_distribution<T>::no_cap()
                      : static_cast<T>(cap);
    bikappa_v1::bi_kappa_distribution<T> dist(static_cast<T>(kappa),
                                              static_cast<T>(theta_perp),
                                              static_cast<T>(theta_par), ubt, capT);
    dist.seed(static_cast<int>(seed));

    ClassRun<T> r;
    r.n_finite= r.n_nonfinite= r.n_thrown= 0;
    r.n_attempts_reported= 0;
    r.n_nonfinite_reported= 0;
    StreamDigest d;
    const std::chrono::steady_clock::time_point t0= std::chrono::steady_clock::now();
    for (long long i= 0; i < n; ++i)
    {
        try
        {
            const P v= dist();
            bool finite= true;
            for (int j= 0; j < 3; ++j)
            {
                d.real<T>(v[j]);
                if (!std::isfinite(v[j]))
                    finite= false;
            }
            if (finite)
                ++r.n_finite;
            else
                ++r.n_nonfinite;
            if (dump)
            {
                LoaderRecord lr;
                for (int j= 0; j < 3; ++j)
                    lr.v[j]= static_cast<double>(v[j]);
                lr.log_r_ref= std::numeric_limits<double>::quiet_NaN();
                lr.status= finite ? kCatFinite : kCatLegacyLoss;
                lr.attempts= 1;
                dump->push(lr);
            }
        }
        catch (const std::exception &)
        {
            ++r.n_thrown;
            d.label("throw");
        }
    }
    r.seconds= std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    d.counter("n_finite", r.n_finite);
    d.counter("n_nonfinite", r.n_nonfinite);
    d.counter("n_thrown", r.n_thrown);
    r.sha256= d.sha256();
    r.fnv1a= d.fnv1a();
    r.digest_bytes= d.byteCount();
    return r;
}

template <typename T>
static ClassRun<T> runCandidateClass(double kappa, double theta_perp, double theta_par,
                                     const double ub[3], double cap, unsigned seed,
                                     long long n, RecordWriter<LoaderRecord> *dump)
{
    typedef typename bi_kappa_distribution<T>::point_type P;
    P ubt;
    ubt[0]= static_cast<T>(ub[0]);
    ubt[1]= static_cast<T>(ub[1]);
    ubt[2]= static_cast<T>(ub[2]);
    const T capT= (cap == std::numeric_limits<double>::infinity())
                      ? bi_kappa_distribution<T>::no_cap()
                      : static_cast<T>(cap);
    bi_kappa_distribution<T> dist(static_cast<T>(kappa), static_cast<T>(theta_perp),
                                  static_cast<T>(theta_par), ubt, capT);
    dist.seed(static_cast<int>(seed));

    ClassRun<T> r;
    r.n_finite= r.n_nonfinite= r.n_thrown= 0;
    StreamDigest d;
    const std::chrono::steady_clock::time_point t0= std::chrono::steady_clock::now();
    for (long long i= 0; i < n; ++i)
    {
        try
        {
            const P v= dist();
            bool finite= true;
            for (int j= 0; j < 3; ++j)
            {
                d.real<T>(v[j]);
                if (!std::isfinite(v[j]))
                    finite= false;
            }
            if (finite)
                ++r.n_finite;
            else
                ++r.n_nonfinite;
            if (dump)
            {
                LoaderRecord lr;
                for (int j= 0; j < 3; ++j)
                    lr.v[j]= static_cast<double>(v[j]);
                lr.log_r_ref= std::numeric_limits<double>::quiet_NaN();
                lr.status= finite ? kCatFinite : kCatHonestOverflow;
                lr.attempts= 1;
                dump->push(lr);
            }
        }
        catch (const std::exception &)
        {
            ++r.n_thrown;
            d.label("throw");
        }
    }
    r.seconds= std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    r.n_attempts_reported= dist.n_attempts();
    r.n_nonfinite_reported= dist.n_nonfinite();
    d.counter("n_finite", r.n_finite);
    d.counter("n_nonfinite", r.n_nonfinite);
    d.counter("n_thrown", r.n_thrown);
    d.counter("n_attempts_reported", static_cast<long long>(r.n_attempts_reported));
    d.counter("n_nonfinite_reported", static_cast<long long>(r.n_nonfinite_reported));
    r.sha256= d.sha256();
    r.fnv1a= d.fnv1a();
    r.digest_bytes= d.byteCount();
    return r;
}

template <typename T>
static void emitClassRow(FILE *f, const RunContext &rc, const char *phase, int method,
                         double kappa, unsigned seed, long long n, const ClassRun<T> &r)
{
    emitConfigFields(f, rc, phase, "class", methodName(method), precisionId<T>(), kappa,
                     static_cast<double>(shapeFromKappa<T>(kappa)), seed, n);
    std::fprintf(f,
                 ",\"n_finite\":%lld,\"nonfinite_output\":%lld,\"thrown\":%lld,"
                 "\"n_attempts_reported\":%llu,\"n_nonfinite_reported\":%llu,"
                 "\"seconds\":%.9g,\"digest_sha256\":\"%s\",\"digest_fnv1a64\":\"%s\","
                 "\"digest_bytes\":%llu}\n",
                 r.n_finite, r.n_nonfinite, r.n_thrown, r.n_attempts_reported,
                 r.n_nonfinite_reported, r.seconds, r.sha256.c_str(), r.fnv1a.c_str(),
                 r.digest_bytes);
}

// ===========================================================================
// P1 - scalar radial law
//
// The native replica of each method is run for 10^6 attempts and its log radius kept, then
// reduced to exactly what families F1-F4 need.  The full per-draw sample is NOT written:
// Experiment 6's raw/ reached 4.1 GB, almost all of it per-draw dumps, and none of F1-F4
// reads a draw twice.  What is written instead is
//
//   - the order statistics at the frozen F2 index brackets, plus the point estimate, which
//     is everything F2 and F3 can use;
//   - the exact Anderson-Darling, Kolmogorov-Smirnov and Cramer-von Mises statistics of
//     Z = -log I_W(a, 3/2) against Exp(1), computed here from the sorted sample;
//   - the empirical CDF of Z on a fixed grid, for the figures and as an independent check
//     on those three statistics;
//   - the exceedance count above each protocol threshold, and *the excesses themselves* in
//     a small companion file, because amendment 1.1.0 Sec. 2.4 measures F4 as the only
//     family with power against survivor conditioning below a 10^-3 loss fraction and F4
//     needs the excesses, not a grid.
//
// A draw whose radius the method could not resolve is sorted to +infinity rather than
// dropped.  That is not a convenience: failure here is overflow, so an unresolved draw *is*
// the largest draw, and dropping it would silently turn the sample into its own survivors.
// The resolved count is reported with every cell so a conditional cell is labelled as one.
// ===========================================================================
struct ScalarSummary
{
    long long n;
    long long n_resolved;   ///< draws whose log R is finite
    long long n_nonresolved;///< draws lost to overflow, sorted to +inf
    double os_lo[16], os_hi[16], os_point[16];
    double ad, ks, cvm;
    long long z_grid[kZGridPoints];
    long long exceed[8];
    long long exceed_nonresolved;
    long long n_tail_records;
    std::string tail_file;
};

/// Reduce the log-radius sample.  `log_r` is consumed (sorted in place).
static ScalarSummary reduceScalar(std::vector<double> &log_r, double a)
{
    ScalarSummary s;
    s.n= static_cast<long long>(log_r.size());
    s.n_tail_records= 0;
    const double kInf= std::numeric_limits<double>::infinity();

    for (size_t i= 0; i < log_r.size(); ++i)
        if (!(log_r[i] == log_r[i]))
            log_r[i]= kInf; // a NaN radius is not a small draw; it is a lost one
    std::sort(log_r.begin(), log_r.end());

    s.n_nonresolved= 0;
    for (size_t i= log_r.size(); i > 0; --i)
    {
        if (log_r[i - 1] == kInf)
            ++s.n_nonresolved;
        else
            break;
    }
    s.n_resolved= s.n - s.n_nonresolved;

    for (int q= 0; q < protocol::kNumQuantileLevels; ++q)
    {
        const long long lo= protocol::kQuantileLoIndex[q];
        const long long hi= protocol::kQuantileHiIndex[q];
        const double p= protocol::kQuantileLevels[q];
        long long pt= static_cast<long long>(std::ceil(p * static_cast<double>(s.n))) - 1;
        if (pt < 0)
            pt= 0;
        if (pt >= s.n)
            pt= s.n - 1;
        s.os_lo[q]= (lo >= 0 && lo < s.n) ? log_r[static_cast<size_t>(lo)]
                                          : std::numeric_limits<double>::quiet_NaN();
        s.os_hi[q]= (hi >= 0 && hi < s.n) ? log_r[static_cast<size_t>(hi)]
                                          : std::numeric_limits<double>::quiet_NaN();
        s.os_point[q]= (s.n > 0) ? log_r[static_cast<size_t>(pt)]
                                 : std::numeric_limits<double>::quiet_NaN();
    }

    // Z is a strictly increasing function of log R, so the resolved prefix of the sorted
    // log-radius sample maps straight onto the sorted Z sample; the re-sort below only
    // removes the last-bit non-monotonicity of the evaluator itself.
    std::vector<double> z;
    z.reserve(static_cast<size_t>(s.n_resolved));
    for (long long i= 0; i < s.n_resolved; ++i)
        z.push_back(zFromLogR(log_r[static_cast<size_t>(i)], a));
    std::sort(z.begin(), z.end());

    for (int i= 0; i < kZGridPoints; ++i)
        s.z_grid[i]= 0;
    {
        // Cumulative counts of Z <= grid point, by a single merge.
        size_t j= 0;
        long long running= 0;
        for (int i= 0; i < kZGridPoints; ++i)
        {
            const double zi= zGridPoint(i);
            while (j < z.size() && z[j] <= zi)
            {
                ++running;
                ++j;
            }
            s.z_grid[i]= running;
        }
    }

    for (int t= 0; t < protocol::kNumTailQ0; ++t)
        s.exceed[t]= 0;
    s.exceed_nonresolved= s.n_nonresolved;
    for (int t= 0; t < protocol::kNumTailQ0; ++t)
    {
        const double z0= -std::log(protocol::kTailQ0[t]);
        long long c= 0;
        for (size_t i= z.size(); i > 0; --i)
        {
            if (z[i - 1] > z0)
                ++c;
            else
                break;
        }
        s.exceed[t]= c;
    }

    // Exact statistics against Exp(1), on the resolved sample.  The formulas are the ones
    // the analysis module uses, including its clipping of the CDF, so the two agree to the
    // last bit rather than approximately.
    const long long nz= static_cast<long long>(z.size());
    if (nz >= 8)
    {
        const double dn= static_cast<double>(nz);
        double sum_ad= 0.0, sum_cvm= 0.0, dmax= 0.0;
        for (long long i= 0; i < nz; ++i)
        {
            double u= -std::expm1(-z[static_cast<size_t>(i)]);
            if (u < 1e-300)
                u= 1e-300;
            if (u > 1.0 - 1e-16)
                u= 1.0 - 1e-16;
            double ur= -std::expm1(-z[static_cast<size_t>(nz - 1 - i)]);
            if (ur < 1e-300)
                ur= 1e-300;
            if (ur > 1.0 - 1e-16)
                ur= 1.0 - 1e-16;
            const double w= 2.0 * static_cast<double>(i + 1) - 1.0;
            sum_ad+= w * (std::log(u) + std::log1p(-ur));
            const double d= u - w / (2.0 * dn);
            sum_cvm+= d * d;
            const double dplus= static_cast<double>(i + 1) / dn - u;
            const double dminus= u - static_cast<double>(i) / dn;
            if (dplus > dmax)
                dmax= dplus;
            if (dminus > dmax)
                dmax= dminus;
        }
        s.ad= -dn - sum_ad / dn;
        s.cvm= 1.0 / (12.0 * dn) + sum_cvm;
        s.ks= dmax;
    }
    else
    {
        s.ad= s.cvm= s.ks= std::numeric_limits<double>::quiet_NaN();
    }
    return s;
}

static void emitScalarSummary(FILE *f, const ScalarSummary &s)
{
    std::fprintf(f, ",\"n_resolved\":%lld,\"n_nonresolved\":%lld,", s.n_resolved,
                 s.n_nonresolved);
    emitDouble(f, "loss_fraction",
               s.n ? static_cast<double>(s.n_nonresolved) / static_cast<double>(s.n) : 0.0);
    std::fprintf(f, ",\"conditional\":%s", s.n_nonresolved ? "true" : "false");

    std::fprintf(f, ",\"quantiles\":[");
    for (int q= 0; q < protocol::kNumQuantileLevels; ++q)
    {
        if (q)
            std::fprintf(f, ",");
        std::fprintf(f, "{\"p\":%.17g,\"lo_index\":%lld,\"hi_index\":%lld,"
                        "\"achieved_coverage\":%.17g,\"lo_value\":",
                     protocol::kQuantileLevels[q], protocol::kQuantileLoIndex[q],
                     protocol::kQuantileHiIndex[q], protocol::kQuantileAchievedCoverage[q]);
        emitDoubleValue(f, s.os_lo[q]);
        std::fprintf(f, ",\"hi_value\":");
        emitDoubleValue(f, s.os_hi[q]);
        std::fprintf(f, ",\"point_value\":");
        emitDoubleValue(f, s.os_point[q]);
        std::fprintf(f, "}");
    }
    std::fprintf(f, "],");

    emitDouble(f, "ad_statistic", s.ad);
    std::fprintf(f, ",");
    emitDouble(f, "ks_statistic", s.ks);
    std::fprintf(f, ",");
    emitDouble(f, "cvm_statistic", s.cvm);

    std::fprintf(f, ",\"tail\":[");
    for (int t= 0; t < protocol::kNumTailQ0; ++t)
    {
        if (t)
            std::fprintf(f, ",");
        std::fprintf(f,
                     "{\"q0\":%.17g,\"z0\":%.17g,\"observed_resolved\":%lld,"
                     "\"observed_nonresolved\":%lld,\"expected\":%.17g}",
                     protocol::kTailQ0[t], -std::log(protocol::kTailQ0[t]), s.exceed[t],
                     s.exceed_nonresolved,
                     static_cast<double>(s.n) * protocol::kTailQ0[t]);
    }
    std::fprintf(f, "]");

    std::fprintf(f, ",\"z_grid_points\":%d,\"z_grid_step\":%.17g,\"z_ecdf_counts\":[",
                 kZGridPoints, kZGridStep);
    for (int i= 0; i < kZGridPoints; ++i)
        std::fprintf(f, "%s%lld", i ? "," : "", s.z_grid[i]);
    std::fprintf(f, "]");

    std::fprintf(f, ",\"tail_records\":%lld,\"tail_file\":\"%s\"", s.n_tail_records,
                 s.tail_file.c_str());
}

template <typename T>
static void p1Native(const Options &o, const RunContext &rc, JsonlWriter &w, int method,
                     double kappa, unsigned seed, long long n)
{
    const std::array<T, 3> zhat= {{T(0), T(0), T(1)}};
    Geometry<T> geom(static_cast<T>(kappa), T(1), T(1), zhat);
    const double log_ovf= logOverflowThreshold<T>();
    const double a_ref= shapeReference<T>(kappa);

    CountingEngine gen(seed);
    LegacyNative<T> legacy(kappa);
    CandidateNative<T> candidate(kappa);

    MethodCounters mc;
    AttemptCounters ac;
    VariateCounts vc;
    StreamDigest d;
    std::vector<double> log_r;
    log_r.reserve(static_cast<size_t>(n));

    const std::chrono::steady_clock::time_point t0= std::chrono::steady_clock::now();
    for (long long i= 0; i < n; ++i)
    {
        NativeAttempt<T> at= (method == kMethodLegacy) ? legacy.attempt(gen, geom)
                                                       : candidate.attempt(gen, geom);
        for (int j= 0; j < 3; ++j)
            d.real<T>(at.v[j]);
        log_r.push_back(at.log_r);
        ++ac.n;

        // Two Gamma variates either way.  The uniform count is taken from the attempt
        // rather than assumed: the candidate's direction is a rejection loop, so it draws
        // one uniform for the boosting identity plus two per try of that loop.
        vc.gamma+= 2;
        vc.uniform+= at.uniform_variates;
        if (method == kMethodLegacy)
        {
            if (at.x2 == T(0))
                ++ac.x2_zero;
            else if (std::fabs(at.x2) < std::numeric_limits<T>::min())
                ++ac.x2_subnormal;
        }

        // Native classification.  The native layer cannot separate avoidable from honest
        // loss on a draw whose denominator materialized as zero, because the intended value
        // is then unknowable from the output; what it *can* do is use the radius the method
        // itself carried.  The authoritative split is P2's paired layer.
        int cat;
        if (at.finite)
            cat= kCatFinite;
        else if (method == kMethodLegacy && at.x2 == T(0))
            cat= kCatDenomZero;
        else if (at.log_r == at.log_r && at.log_r < log_ovf)
            cat= (method == kMethodLegacy) ? kCatLegacyLoss : kCatLogPrimFail;
        else
            cat= kCatHonestOverflow;
        mc.observe(cat, at.finite, std::numeric_limits<double>::quiet_NaN(),
                   at.finite ? at.log_r : -std::numeric_limits<double>::infinity());
        if (at.log_r > ac.max_log_r_ref && at.log_r != std::numeric_limits<double>::infinity())
            ac.max_log_r_ref= at.log_r;
    }
    const double secs=
        std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    vc.engine_calls= gen.calls();
    digestCounters(&d, mc, vc, n);

    ScalarSummary s= reduceScalar(log_r, a_ref);

    // The excesses above the lowest protocol threshold, which nests the other two.
    {
        char buf[320];
        std::snprintf(buf, sizeof(buf), "p1_tail_%s_%s_k%g_s%u_%s.bin", rc.tag.c_str(),
                      methodName(method), kappa, seed, precisionId<T>());
        const double ub[3]= {0.0, 0.0, 1.0};
        RecordWriter<TailRecord> tw(joinPath(phaseDir(o, "p1"), buf), kTailSchemaVersion,
                                    kKindTail, kappa, 1.0, 1.0, ub,
                                    std::numeric_limits<double>::infinity(), seed,
                                    sizeof(T) == sizeof(float));
        double z0min= -std::log(protocol::kTailQ0[0]);
        for (int t= 1; t < protocol::kNumTailQ0; ++t)
        {
            const double z0= -std::log(protocol::kTailQ0[t]);
            if (z0 < z0min)
                z0min= z0;
        }
        // log_r has been sorted in place by reduceScalar, so the resolved prefix is sorted
        // and Z is increasing along it: walk back from the largest until it drops below z0.
        for (long long i= s.n_resolved; i > 0; --i)
        {
            const double z= zFromLogR(log_r[static_cast<size_t>(i - 1)], a_ref);
            if (!(z > z0min))
                break;
            TailRecord tr;
            tr.z= z;
            tw.push(tr);
        }
        s.n_tail_records= tw.count();
        s.tail_file= tw.finish(rc, "p1", "tail", methodName(method), precisionId<T>(), "-",
                               kappa, static_cast<long long>(seed));
    }

    emitConfigFields(w.f(), rc, "p1", "native", methodName(method), precisionId<T>(), kappa,
                     a_ref, seed, n);
    emitCategories(w.f(), mc);
    std::fprintf(w.f(),
                 ",\"x2_zero\":%lld,\"x2_subnormal\":%lld,\"gamma_variates\":%lld,"
                 "\"uniform_variates\":%lld,\"engine_calls\":%lld,\"seconds\":%.9g,",
                 ac.x2_zero, ac.x2_subnormal, vc.gamma, vc.uniform, vc.engine_calls, secs);
    emitDouble(w.f(), "max_finite_log_r", ac.max_log_r_ref);
    emitScalarSummary(w.f(), s);
    emitDigest(w.f(), "digest", d);
    std::fprintf(w.f(), "}\n");
    w.countLine();
}

template <typename T>
static void p1ForPrecision(const Options &o, const RunContext &rc, JsonlWriter &w,
                           const std::vector<double> &kappas,
                           const std::vector<unsigned> &seeds, long long n)
{
    const double zhat[3]= {0.0, 0.0, 1.0};
    for (size_t ik= 0; ik < kappas.size(); ++ik)
    {
        const double k= kappas[ik];
        for (size_t is= 0; is < seeds.size(); ++is)
        {
            const unsigned s= seeds[is];
            const ClassRun<T> lg=
                runLegacyClass<T>(k, 1.0, 1.0, zhat, std::numeric_limits<double>::infinity(),
                                  s, n, 0);
            emitClassRow<T>(w.f(), rc, "p1", kMethodLegacy, k, s, n, lg);
            w.countLine();
            const ClassRun<T> cd= runCandidateClass<T>(
                k, 1.0, 1.0, zhat, std::numeric_limits<double>::infinity(), s, n, 0);
            emitClassRow<T>(w.f(), rc, "p1", kMethodCandidate, k, s, n, cd);
            w.countLine();

            p1Native<T>(o, rc, w, kMethodLegacy, k, s, n);
            p1Native<T>(o, rc, w, kMethodCandidate, k, s, n);
        }
    }
}

static int phaseP1(const Options &o, const RunContext &rc)
{
    const std::vector<double> kappas= kappaLadder();
    const std::vector<unsigned> seeds= seedsFor(o);
    const long long n= o.n > 0 ? o.n : protocol::kNScalar;
    const std::string path= joinPath(phaseDir(o, "p1"), "p1_" + o.tag + ".jsonl");
    JsonlWriter w(path);
    if (o.only_precision != "float")
        p1ForPrecision<double>(o, rc, w, kappas, seeds, n);
    if (o.only_precision != "double")
        p1ForPrecision<float>(o, rc, w, kappas, seeds, n);
    const std::string out= w.finish(rc, "p1");
    std::fprintf(stderr, "exp7 P1: wrote %lld rows to %s\n", w.lines(), out.c_str());
    return 0;
}

// ===========================================================================
// P2 - mechanism
// ===========================================================================
/// Which stratum an attempt belongs to.  Strata are ordered by decision value and an attempt
/// is assigned to the first one it matches, so the counts partition the sample and the
/// reported coverage of each stratum is a coverage of a disjoint class.
inline int auditStratumOf(const PairedResult &pr)
{
    if (pr.ref.near_limit)
        return kAuditWithinMargin;
    const bool all_same= pr.cat[kMethodLegacy] == pr.cat[kMethodCandidate] &&
                         pr.cat[kMethodCandidate] == pr.cat[kMethodQF];
    if (!all_same)
        return kAuditDisagree;
    if (pr.cat[kMethodCandidate] == kCatFiniteButWrong)
        return kAuditFiniteButWrong;
    if (isAvoidable(pr.cat[kMethodCandidate]))
        return kAuditAvoidableCandidate;
    if ((pr.flags & (kFlagX2Zero | kFlagX2Subnormal)) && pr.ref.representable)
        return kAuditSubnormalRepresentable;
    if (!pr.finite[kMethodLegacy] || !pr.finite[kMethodCandidate])
        return kAuditUnambiguous;
    return kAuditSample;
}

template <typename T>
static void maybeAudit(RecordWriter<AuditRecord> *aw, const PairedResult &pr,
                       const SharedPrimitives<T> &s, AttemptCounters *ac, double kappa)
{
    const int stratum= auditStratumOf(pr);
    ++ac->stratum_total[stratum];
    if (stratum == kAuditUnambiguous)
    {
        // The margin, asserted rather than assumed.  `near_limit` already decided this, but
        // the protocol asks for the inequality to be evaluated and reported for every
        // attempt the run declines to adjudicate, so it is evaluated and reported.
        const double m= std::fabs(pr.ref.overflow_margin);
        ++ac->margin_assertions;
        if (!(m > oracleBandLogUnits()))
            ++ac->margin_assertion_failures;
        if (m < ac->min_unambiguous_margin)
            ac->min_unambiguous_margin= m;
    }
    if (!aw)
        return;
    const long long stride= auditStride(stratum);
    if (stride <= 0)
        return;
    if ((ac->stratum_total[stratum] - 1) % stride != 0)
        return;
    ++ac->stratum_audited[stratum];
    AuditRecord a;
    a.x1= static_cast<double>(s.x1);
    a.y= static_cast<double>(s.y);
    a.u= static_cast<double>(s.u);
    a.cos_theta= static_cast<double>(s.cos_theta);
    a.phi= static_cast<double>(s.phi);
    a.log_x2_working= static_cast<double>(s.log_x2);
    a.kappa= kappa;
    a.flags= pr.flags | kFlagAudited;
    a.precision_is_float= (sizeof(T) == sizeof(float)) ? 1u : 0u;
    a.pad= 0;
    a.cat_legacy= static_cast<uint8_t>(pr.cat[kMethodLegacy]);
    a.cat_candidate= static_cast<uint8_t>(pr.cat[kMethodCandidate]);
    a.cat_qf= static_cast<uint8_t>(pr.cat[kMethodQF]);
    a.reserved= 0;
    aw->push(a);
    ++ac->audited;
}

template <typename T>
static void p2Paired(const Options &o, const RunContext &rc, JsonlWriter &w, double kappa,
                     unsigned seed, long long n, RecordWriter<AuditRecord> *aw)
{
    (void)o;
    const std::array<T, 3> zhat= {{T(0), T(0), T(1)}};
    Geometry<T> geom(static_cast<T>(kappa), T(1), T(1), zhat);
    const T a= shapeFromKappa<T>(kappa);
    const double log_ovf= logOverflowThreshold<T>();

    CountingEngine gen(seed);
    T ncache= T(0);
    bool ncached= false;

    MethodCounters mc[kNumMethods];
    AttemptCounters ac;

    for (long long i= 0; i < n; ++i)
    {
        const SharedPrimitives<T> s= drawShared<T>(a, gen, ncache, ncached);
        const PairedResult pr= runPairedAttempt<T>(kappa, s, geom, log_ovf);

        ++ac.n;
        if (pr.flags & kFlagX2Zero)
            ++ac.x2_zero;
        if (pr.flags & kFlagX2Subnormal)
            ++ac.x2_subnormal;
        if (pr.flags & kFlagHonestOverflowRef)
            ++ac.ref_nonrepresentable;
        if (pr.flags & kFlagRotationRecoverable)
            ++ac.rotation_recoverable;
        if (pr.flags & kFlagNearLimit)
            ++ac.near_limit;
        if (pr.ref.log_r > ac.max_log_r_ref)
            ac.max_log_r_ref= pr.ref.log_r;
        for (int m= 0; m < kNumMethods; ++m)
            mc[m].observe(pr.cat[m], pr.finite[m], pr.rel_err[m], pr.ref.max_log_component);
        maybeAudit<T>(aw, pr, s, &ac, kappa);
    }

    for (int m= 0; m < kNumMethods; ++m)
    {
        emitConfigFields(w.f(), rc, "p2", "paired", methodName(m), precisionId<T>(), kappa,
                         static_cast<double>(a), seed, n);
        emitCategories(w.f(), mc[m]);
        std::fprintf(w.f(),
                     ",\"diagnostic_only\":%s,\"x2_zero\":%lld,\"x2_subnormal\":%lld,"
                     "\"ref_nonrepresentable\":%lld,\"rotation_recoverable\":%lld,"
                     "\"near_limit\":%lld,\"audited\":%lld,",
                     isTestedMethod(m) ? "false" : "true", ac.x2_zero, ac.x2_subnormal,
                     ac.ref_nonrepresentable, ac.rotation_recoverable, ac.near_limit,
                     ac.audited);
        emitDouble(w.f(), "max_rel_error_threshold", maxRelErrorFor<T>());
        emitAuditStrata(w.f(), ac);
        std::fprintf(w.f(), ",");
        emitDouble(w.f(), "max_log_r_ref", ac.max_log_r_ref);
        std::fprintf(w.f(), ",\"accounting_ok\":%s}\n",
                     mc[m].accountingHolds(n) ? "true" : "false");
        w.countLine();
    }
}

template <typename T>
static void p2Native(const Options &o, const RunContext &rc, JsonlWriter &w, int method,
                     double kappa, unsigned seed, long long n)
{
    (void)o;
    const std::array<T, 3> zhat= {{T(0), T(0), T(1)}};
    Geometry<T> geom(static_cast<T>(kappa), T(1), T(1), zhat);
    const double log_ovf= logOverflowThreshold<T>();

    CountingEngine gen(seed);
    LegacyNative<T> legacy(kappa);
    CandidateNative<T> candidate(kappa);
    MethodCounters mc;
    AttemptCounters ac;
    VariateCounts vc;
    StreamDigest d;

    const std::chrono::steady_clock::time_point t0= std::chrono::steady_clock::now();
    for (long long i= 0; i < n; ++i)
    {
        NativeAttempt<T> at= (method == kMethodLegacy) ? legacy.attempt(gen, geom)
                                                       : candidate.attempt(gen, geom);
        for (int j= 0; j < 3; ++j)
            d.real<T>(at.v[j]);
        ++ac.n;
        // Two Gamma variates either way.  The uniform count is taken from the attempt
        // rather than assumed: the candidate's direction is a rejection loop, so it draws
        // one uniform for the boosting identity plus two per try of that loop.
        vc.gamma+= 2;
        vc.uniform+= at.uniform_variates;
        if (method == kMethodLegacy)
        {
            if (at.x2 == T(0))
                ++ac.x2_zero;
            else if (std::fabs(at.x2) < std::numeric_limits<T>::min())
                ++ac.x2_subnormal;
        }
        int cat;
        if (at.finite)
            cat= kCatFinite;
        else if (method == kMethodLegacy && at.x2 == T(0))
            cat= kCatDenomZero;
        else if (at.log_r == at.log_r && at.log_r < log_ovf)
            cat= (method == kMethodLegacy) ? kCatLegacyLoss : kCatLogPrimFail;
        else
            cat= kCatHonestOverflow;
        mc.observe(cat, at.finite, std::numeric_limits<double>::quiet_NaN(),
                   at.finite ? at.log_r : -std::numeric_limits<double>::infinity());
    }
    const double secs=
        std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    vc.engine_calls= gen.calls();
    digestCounters(&d, mc, vc, n);

    emitConfigFields(w.f(), rc, "p2", "native", methodName(method), precisionId<T>(), kappa,
                     static_cast<double>(shapeFromKappa<T>(kappa)), seed, n);
    emitCategories(w.f(), mc);
    std::fprintf(w.f(),
                 ",\"x2_zero\":%lld,\"x2_subnormal\":%lld,\"gamma_variates\":%lld,"
                 "\"uniform_variates\":%lld,\"engine_calls\":%lld,\"seconds\":%.9g,"
                 "\"accounting_ok\":%s",
                 ac.x2_zero, ac.x2_subnormal, vc.gamma, vc.uniform, vc.engine_calls, secs,
                 mc.accountingHolds(n) ? "true" : "false");
    emitDigest(w.f(), "digest", d);
    std::fprintf(w.f(), "}\n");
    w.countLine();
}

static int phaseP2(const Options &o, const RunContext &rc)
{
    const std::vector<double> kappas= kappaLadder();
    const std::vector<unsigned> seeds= seedsFor(o);
    const long long n= o.n > 0 ? o.n : protocol::kNMechanism;

    const std::string dir= phaseDir(o, "p2");
    const std::string path= joinPath(dir, "p2_" + o.tag + ".jsonl");
    JsonlWriter w(path);

    const double ubz[3]= {0.0, 0.0, 1.0};
    const std::string apath= joinPath(dir, "audit_p2_" + o.tag + ".bin");
    RecordWriter<AuditRecord> aw(apath, kAuditSchemaVersion, kKindAudit, 0.0, 1.0, 1.0, ubz,
                                 std::numeric_limits<double>::infinity(), 0, 0);

    for (size_t ik= 0; ik < kappas.size(); ++ik)
        for (size_t is= 0; is < seeds.size(); ++is)
        {
            const double k= kappas[ik];
            const unsigned s= seeds[is];
            if (o.only_precision != "float")
            {
                p2Paired<double>(o, rc, w, k, s, n, &aw);
                p2Native<double>(o, rc, w, kMethodLegacy, k, s, n);
                p2Native<double>(o, rc, w, kMethodCandidate, k, s, n);
            }
            if (o.only_precision != "double")
            {
                p2Paired<float>(o, rc, w, k, s, n, &aw);
                p2Native<float>(o, rc, w, kMethodLegacy, k, s, n);
                p2Native<float>(o, rc, w, kMethodCandidate, k, s, n);
            }
        }

    const std::string aout= aw.finish(rc, "p2", "audit", "paired", "mixed", "-",
                                      std::numeric_limits<double>::quiet_NaN(), -1);
    const std::string out= w.finish(rc, "p2");
    std::fprintf(stderr, "exp7 P2: wrote %lld rows to %s and %lld audit records to %s\n",
                 w.lines(), out.c_str(), aw.count(), aout.c_str());
    return 0;
}

// ===========================================================================
// P3 - conditioning
// ===========================================================================
template <typename T>
static void p3One(const Options &o, const RunContext &rc, JsonlWriter &w, double kappa,
                  unsigned seed, long long n, RecordWriter<AuditRecord> *aw)
{
    const std::array<T, 3> zhat= {{T(0), T(0), T(1)}};
    Geometry<T> geom(static_cast<T>(kappa), T(1), T(1), zhat);
    const T a= shapeFromKappa<T>(kappa);
    const double log_ovf= logOverflowThreshold<T>();

    CountingEngine gen(seed);
    T ncache= T(0);
    bool ncached= false;

    char buf[320];
    std::snprintf(buf, sizeof(buf), "p3_cond_%s_k%g_s%u_%s.bin", rc.tag.c_str(), kappa, seed,
                  precisionId<T>());
    const double ub[3]= {0.0, 0.0, 1.0};
    RecordWriter<CondRecord> rw(joinPath(phaseDir(o, "p3"), buf), kCondSchemaVersion,
                                kKindCond, kappa, 1.0, 1.0, ub,
                                std::numeric_limits<double>::infinity(), seed,
                                sizeof(T) == sizeof(float));

    MethodCounters mc[kNumMethods];
    AttemptCounters ac;
    for (long long i= 0; i < n; ++i)
    {
        const SharedPrimitives<T> s= drawShared<T>(a, gen, ncache, ncached);
        const PairedResult pr= runPairedAttempt<T>(kappa, s, geom, log_ovf);

        CondRecord cr;
        cr.log_r_ref= pr.ref.log_r;
        cr.log_speed_ref= pr.ref.log_speed;
        cr.max_log_component_ref= pr.ref.max_log_component;
        cr.flags= pr.flags;
        cr.cat_legacy= static_cast<uint8_t>(pr.cat[kMethodLegacy]);
        cr.cat_candidate= static_cast<uint8_t>(pr.cat[kMethodCandidate]);
        cr.cat_qf= static_cast<uint8_t>(pr.cat[kMethodQF]);
        cr.reserved= 0;
        rw.push(cr);

        ++ac.n;
        if (pr.flags & kFlagX2Zero)
            ++ac.x2_zero;
        if (pr.flags & kFlagX2Subnormal)
            ++ac.x2_subnormal;
        if (pr.flags & kFlagHonestOverflowRef)
            ++ac.ref_nonrepresentable;
        if (pr.flags & kFlagNearLimit)
            ++ac.near_limit;
        if (pr.ref.log_r > ac.max_log_r_ref)
            ac.max_log_r_ref= pr.ref.log_r;
        for (int m= 0; m < kNumMethods; ++m)
            mc[m].observe(pr.cat[m], pr.finite[m], pr.rel_err[m], pr.ref.max_log_component);
        maybeAudit<T>(aw, pr, s, &ac, kappa);
    }
    const std::string out= rw.finish(rc, "p3", "paired", "both", precisionId<T>(), "-", kappa,
                                     static_cast<long long>(seed));

    for (int m= 0; m < kNumMethods; ++m)
    {
        emitConfigFields(w.f(), rc, "p3", "paired", methodName(m), precisionId<T>(), kappa,
                         static_cast<double>(a), seed, n);
        emitCategories(w.f(), mc[m]);
        std::fprintf(w.f(),
                     ",\"diagnostic_only\":%s,\"x2_zero\":%lld,\"x2_subnormal\":%lld,"
                     "\"ref_nonrepresentable\":%lld,\"near_limit\":%lld,\"audited\":%lld,",
                     isTestedMethod(m) ? "false" : "true", ac.x2_zero, ac.x2_subnormal,
                     ac.ref_nonrepresentable, ac.near_limit, ac.audited);
        emitDouble(w.f(), "max_rel_error_threshold", maxRelErrorFor<T>());
        emitAuditStrata(w.f(), ac);
        std::fprintf(w.f(), ",\"accounting_ok\":%s,\"raw_file\":\"%s\"}\n",
                     mc[m].accountingHolds(n) ? "true" : "false", out.c_str());
        w.countLine();
    }
}

static int phaseP3(const Options &o, const RunContext &rc)
{
    // PROTOCOL.md Sec. 4: double at the rescue-relevant shapes, float at its own boundary,
    // plus a benign control in each.
    const double dk[]= {0.501, 0.505, 0.51, 0.75};
    const double fk[]= {0.55, 0.60, 0.75};
    const std::vector<unsigned> seeds= seedsFor(o);
    const long long n= o.n > 0 ? o.n : protocol::kNConditioning;

    const std::string dir= phaseDir(o, "p3");
    const std::string path= joinPath(dir, "p3_" + o.tag + ".jsonl");
    JsonlWriter w(path);
    const double ubz[3]= {0.0, 0.0, 1.0};
    RecordWriter<AuditRecord> aw(joinPath(dir, "audit_p3_" + o.tag + ".bin"),
                                 kAuditSchemaVersion, kKindAudit, 0.0, 1.0, 1.0, ubz,
                                 std::numeric_limits<double>::infinity(), 0, 0);

    for (size_t is= 0; is < seeds.size(); ++is)
    {
        if (o.only_precision != "float")
            for (size_t i= 0; i < sizeof(dk) / sizeof(dk[0]); ++i)
                p3One<double>(o, rc, w, dk[i], seeds[is], n, &aw);
        if (o.only_precision != "double")
            for (size_t i= 0; i < sizeof(fk) / sizeof(fk[0]); ++i)
                p3One<float>(o, rc, w, fk[i], seeds[is], n, &aw);
    }
    const std::string aout= aw.finish(rc, "p3", "audit", "paired", "mixed", "-",
                                      std::numeric_limits<double>::quiet_NaN(), -1);
    const std::string out= w.finish(rc, "p3");
    std::fprintf(stderr, "exp7 P3: wrote %lld rows to %s and %lld audit records to %s\n",
                 w.lines(), out.c_str(), aw.count(), aout.c_str());
    return 0;
}

// ===========================================================================
// P4 - complete anisotropic, rotated, capped 3-D loader
// ===========================================================================
struct CaseSpec
{
    const char *id;
    bool is_float;
    double kappa;
    double theta_ratio; ///< theta_par / theta_perp, with theta_perp = 1
    double ub[3];
    double cap; ///< +inf for uncapped
    const char *role;
};

static std::vector<CaseSpec> caseTable()
{
    const double s14= 1.0 / std::sqrt(14.0);
    const double inf= std::numeric_limits<double>::infinity();
    std::vector<CaseSpec> v;
    CaseSpec c;
    c.id= "C0"; c.is_float= false; c.kappa= 2.0;   c.theta_ratio= 2.0;
    c.ub[0]= 0; c.ub[1]= 0; c.ub[2]= 1; c.cap= inf; c.role= "benign reference"; v.push_back(c);
    c.id= "C1"; c.is_float= false; c.kappa= 0.51;  c.theta_ratio= 2.0;
    c.ub[0]= s14; c.ub[1]= 2 * s14; c.ub[2]= 3 * s14; c.cap= inf;
    c.role= "primary rescued case"; v.push_back(c);
    c.id= "C2"; c.is_float= false; c.kappa= 0.51;  c.theta_ratio= 0.5;
    c.ub[0]= -2.0 / 3.0; c.ub[1]= 1.0 / 3.0; c.ub[2]= 2.0 / 3.0; c.cap= inf;
    c.role= "inverse anisotropy + rotation"; v.push_back(c);
    c.id= "C3"; c.is_float= true;  c.kappa= 0.55;  c.theta_ratio= 2.0;
    c.ub[0]= s14; c.ub[1]= 2 * s14; c.ub[2]= 3 * s14; c.cap= inf;
    c.role= "single-precision boundary"; v.push_back(c);
    c.id= "C4"; c.is_float= false; c.kappa= 0.505; c.theta_ratio= 1.0;
    c.ub[0]= 0; c.ub[1]= 0; c.ub[2]= 1; c.cap= inf;
    c.role= "honest-overflow semantics"; v.push_back(c);
    c.id= "C5"; c.is_float= false; c.kappa= 0.51;  c.theta_ratio= 2.0;
    c.ub[0]= s14; c.ub[1]= 2 * s14; c.ub[2]= 3 * s14; c.cap= 5.0;
    c.role= "strong bounded target"; v.push_back(c);
    c.id= "C6"; c.is_float= false; c.kappa= 0.51;  c.theta_ratio= 2.0;
    c.ub[0]= s14; c.ub[1]= 2 * s14; c.ub[2]= 3 * s14; c.cap= 20.0;
    c.role= "weak bounded target"; v.push_back(c);
    return v;
}

/// The released headers' own internal cap-rejection limit.  Reproduced, not chosen.
static const int kMaxCapTries= 1000000;

/// One method's complete loader, run through its verified native replica so that attempts,
/// rejections and the intended radius are observable.
///
/// Uncapped cases stop after exactly `n` *intended attempts*.  Capped cases stop after `n`
/// *returned samples* and keep every attempt and rejection counter, which is how the
/// protocol words it: an uncapped run has no rejection, so its attempt count is the honest
/// unit, while a capped run's returned count is what a caller asked for.
template <typename T>
static void p4One(const Options &o, const RunContext &rc, JsonlWriter &w, const CaseSpec &cs,
                  int method, unsigned seed, long long n)
{
    const std::array<T, 3> ub= {{static_cast<T>(cs.ub[0]), static_cast<T>(cs.ub[1]),
                                 static_cast<T>(cs.ub[2])}};
    Geometry<T> geom(static_cast<T>(cs.kappa), T(1), static_cast<T>(cs.theta_ratio), ub);
    const bool capped= (cs.cap != std::numeric_limits<double>::infinity());
    const T lambda= static_cast<T>(cs.cap);
    const T log_cap= capped ? static_cast<T>(std::log(cs.cap)) : T(0);

    CountingEngine gen(seed);
    LegacyNative<T> legacy(cs.kappa);
    CandidateNative<T> candidate(cs.kappa);

    char buf[320];
    std::snprintf(buf, sizeof(buf), "p4_loader_%s_%s_%s_s%u.bin", rc.tag.c_str(), cs.id,
                  methodName(method), seed);
    RecordWriter<LoaderRecord> rw(joinPath(phaseDir(o, "p4"), buf), kLoaderSchemaVersion,
                                  kKindLoader, cs.kappa, 1.0, cs.theta_ratio, cs.ub, cs.cap,
                                  seed, cs.is_float);

    long long attempts= 0, returned= 0, nonfinite_attempt= 0, cap_reject= 0, exhausted= 0,
              nonfinite_returned= 0, x2_zero= 0, cap_reject_unrepresentable= 0;
    long long max_run= 0;
    StreamDigest d;
    const std::chrono::steady_clock::time_point t0= std::chrono::steady_clock::now();

    while ((capped && returned < n) || (!capped && attempts < n))
    {
        long long this_attempts= 0;
        bool emitted= false;
        for (int t= 0; t < (capped ? kMaxCapTries : 1); ++t)
        {
            ++attempts;
            ++this_attempts;

            std::array<T, 3> v;
            double log_r_ref= std::numeric_limits<double>::quiet_NaN();
            bool finite= false;
            uint32_t status= kCatFinite;

            if (method == kMethodLegacy)
            {
                // The legacy loader materializes the radius, forms the local vector, tests
                // the box on it, and rotates last.  Reproduced in that order.
                NativeAttempt<T> at= legacy.attempt(gen, geom);
                const T r= legacy.lastRadius();
                if (at.x2 == T(0))
                    ++x2_zero;
                log_r_ref= at.log_r;
                if (capped && !legacyCapAccept<T>(r, legacy.lastCosTheta(), legacy.lastPhi(),
                                                  geom, lambda))
                {
                    ++cap_reject;
                    if (!at.finite)
                        ++cap_reject_unrepresentable;
                    continue;
                }
                v= at.v;
                finite= at.finite;
                if (!finite)
                    ++nonfinite_attempt;
                status= finite ? kCatFinite
                               : (at.x2 == T(0) ? kCatDenomZero : kCatLegacyLoss);
            }
            else
            {
                // The candidate tests the box in the log domain, before exponentiating.
                // Its attempt() does not, so the test is applied here exactly as the header
                // applies it and the attempt is only completed on acceptance.
                NativeAttempt<T> at= candidate.attempt(gen, geom);
                log_r_ref= at.log_r;
                if (capped)
                {
                    bool unrep= false;
                    if (!candidateCapAcceptLog<T>(static_cast<T>(at.log_r),
                                                  candidate.lastDirection(), geom, log_cap,
                                                  &unrep))
                    {
                        ++cap_reject;
                        if (unrep)
                            ++cap_reject_unrepresentable;
                        continue;
                    }
                }
                v= at.v;
                finite= at.finite;
                if (!finite)
                    ++nonfinite_attempt;
                status= finite ? kCatFinite : kCatHonestOverflow;
            }

            LoaderRecord lr;
            for (int j= 0; j < 3; ++j)
            {
                lr.v[j]= static_cast<double>(v[j]);
                d.real<T>(v[j]);
            }
            lr.log_r_ref= log_r_ref;
            lr.status= status;
            lr.attempts= static_cast<uint32_t>(this_attempts);
            rw.push(lr);
            ++returned;
            if (!finite)
                ++nonfinite_returned;
            if (this_attempts > max_run)
                max_run= this_attempts;
            emitted= true;
            break;
        }
        if (!emitted)
        {
            ++exhausted;
            ++returned; // an exhausted request still consumes one of the n slots
        }
    }
    const double secs=
        std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    d.counter("attempts", attempts);
    d.counter("returned", returned);
    d.counter("cap_reject", cap_reject);
    d.counter("nonfinite_returned", nonfinite_returned);
    const std::string out= rw.finish(rc, "p4", capped ? "capped" : "uncapped",
                                     methodName(method), precisionId<T>(), cs.id, cs.kappa,
                                     static_cast<long long>(seed));

    emitConfigFields(w.f(), rc, "p4", capped ? "capped" : "uncapped", methodName(method),
                     precisionId<T>(),
                     cs.kappa, static_cast<double>(shapeFromKappa<T>(cs.kappa)), seed,
                     attempts);
    std::fprintf(w.f(), ",\"case\":\"%s\",\"role\":\"%s\",\"theta_ratio\":%.17g,"
                        "\"ub\":[%.17g,%.17g,%.17g],",
                 cs.id, cs.role, cs.theta_ratio, cs.ub[0], cs.ub[1], cs.ub[2]);
    emitCap(w.f(), cs.cap);
    std::fprintf(w.f(),
                 ",\"n_returned\":%lld,\"attempts\":%lld,\"nonfinite_attempt\":%lld,"
                 "\"nonfinite_returned\":%lld,\"cap_reject\":%lld,"
                 "\"cap_reject_unrepresentable\":%lld,\"cap_exhausted\":%lld,"
                 "\"max_attempt_run\":%lld,\"x2_zero\":%lld,\"engine_calls\":%lld,"
                 "\"seconds\":%.9g,\"raw_file\":\"%s\"",
                 returned, attempts, nonfinite_attempt, nonfinite_returned, cap_reject,
                 cap_reject_unrepresentable, exhausted, max_run, x2_zero, gen.calls(), secs,
                 out.c_str());
    emitDigest(w.f(), "digest", d);
    std::fprintf(w.f(), "}\n");
    w.countLine();
}

/// The released classes themselves on the uncapped cases, dumped so the analysis can check
/// the replicas against them on the very configurations the loader claim is about.
template <typename T>
static void p4ClassDump(const Options &o, const RunContext &rc, JsonlWriter &w,
                        const CaseSpec &cs, int method, unsigned seed, long long n)
{
    char buf[320];
    std::snprintf(buf, sizeof(buf), "p4_class_%s_%s_%s_s%u.bin", rc.tag.c_str(), cs.id,
                  methodName(method), seed);
    RecordWriter<LoaderRecord> rw(joinPath(phaseDir(o, "p4"), buf), kLoaderSchemaVersion,
                                  kKindLoader, cs.kappa, 1.0, cs.theta_ratio, cs.ub, cs.cap,
                                  seed, cs.is_float);
    ClassRun<T> r;
    if (method == kMethodLegacy)
        r= runLegacyClass<T>(cs.kappa, 1.0, cs.theta_ratio, cs.ub, cs.cap, seed, n, &rw);
    else
        r= runCandidateClass<T>(cs.kappa, 1.0, cs.theta_ratio, cs.ub, cs.cap, seed, n, &rw);
    const std::string out= rw.finish(rc, "p4", "class", methodName(method), precisionId<T>(),
                                     cs.id, cs.kappa, static_cast<long long>(seed));

    emitConfigFields(w.f(), rc, "p4", "class", methodName(method), precisionId<T>(), cs.kappa,
                     static_cast<double>(shapeFromKappa<T>(cs.kappa)), seed, n);
    std::fprintf(w.f(), ",\"case\":\"%s\",\"role\":\"%s\",\"theta_ratio\":%.17g,", cs.id,
                 cs.role, cs.theta_ratio);
    emitCap(w.f(), cs.cap);
    std::fprintf(w.f(),
                 ",\"n_finite\":%lld,\"nonfinite_output\":%lld,\"thrown\":%lld,"
                 "\"n_attempts_reported\":%llu,\"n_nonfinite_reported\":%llu,"
                 "\"seconds\":%.9g,\"digest_sha256\":\"%s\",\"digest_fnv1a64\":\"%s\","
                 "\"digest_bytes\":%llu,\"raw_file\":\"%s\"}\n",
                 r.n_finite, r.n_nonfinite, r.n_thrown, r.n_attempts_reported,
                 r.n_nonfinite_reported, r.seconds, r.sha256.c_str(), r.fnv1a.c_str(),
                 r.digest_bytes, out.c_str());
    w.countLine();
}

static int phaseP4(const Options &o, const RunContext &rc)
{
    const std::vector<CaseSpec> cases= caseTable();
    const std::vector<unsigned> seeds= seedsFor(o);
    const long long n= o.n > 0 ? o.n : protocol::kNLoader;

    const std::string path= joinPath(phaseDir(o, "p4"), "p4_" + o.tag + ".jsonl");
    JsonlWriter w(path);
    for (size_t ic= 0; ic < cases.size(); ++ic)
    {
        const CaseSpec &cs= cases[ic];
        if (!o.case_id.empty() && o.case_id != "all" && o.case_id != cs.id)
            continue;
        const bool uncapped= (cs.cap == std::numeric_limits<double>::infinity());
        for (size_t is= 0; is < seeds.size(); ++is)
        {
            const unsigned s= seeds[is];
            if (cs.is_float)
            {
                p4One<float>(o, rc, w, cs, kMethodLegacy, s, n);
                p4One<float>(o, rc, w, cs, kMethodCandidate, s, n);
                if (uncapped)
                {
                    p4ClassDump<float>(o, rc, w, cs, kMethodLegacy, s, n);
                    p4ClassDump<float>(o, rc, w, cs, kMethodCandidate, s, n);
                }
            }
            else
            {
                p4One<double>(o, rc, w, cs, kMethodLegacy, s, n);
                p4One<double>(o, rc, w, cs, kMethodCandidate, s, n);
                if (uncapped)
                {
                    p4ClassDump<double>(o, rc, w, cs, kMethodLegacy, s, n);
                    p4ClassDump<double>(o, rc, w, cs, kMethodCandidate, s, n);
                }
            }
        }
    }
    const std::string out= w.finish(rc, "p4");
    std::fprintf(stderr, "exp7 P4: wrote %lld rows to %s\n", w.lines(), out.c_str());
    return 0;
}

// ===========================================================================
// P5 - portability
//
// CANDIDATE only, as PROTOCOL.md Sec. 4 freezes it: the released 2.0.0 class run end to end
// on the full P1 ladder, with a digest over the returned bit patterns and over every
// counter.  Two builds of this phase against different standard libraries on one
// architecture must produce identical digests; that is an equality, not an overlap.
//
// The negative half of the contrast -- the same comparison for LEGACY, where it does not
// hold -- comes from P1's class rows, which cover the same ladder, seeds and geometry.  It
// is not duplicated here, because P5's method list is frozen.
// ===========================================================================
template <typename T>
static void p5ForPrecision(const Options &o, const RunContext &rc, JsonlWriter &w,
                           const std::vector<double> &kappas,
                           const std::vector<unsigned> &seeds, long long n)
{
    (void)o;
    const double zhat[3]= {0.0, 0.0, 1.0};
    for (size_t ik= 0; ik < kappas.size(); ++ik)
        for (size_t is= 0; is < seeds.size(); ++is)
        {
            const ClassRun<T> cd= runCandidateClass<T>(
                kappas[ik], 1.0, 1.0, zhat, std::numeric_limits<double>::infinity(),
                seeds[is], n, 0);
            emitClassRow<T>(w.f(), rc, "p5", kMethodCandidate, kappas[ik], seeds[is], n, cd);
            w.countLine();
        }
}

static int phaseP5(const Options &o, const RunContext &rc)
{
    const std::vector<double> kappas= kappaLadder();
    const std::vector<unsigned> seeds= seedsFor(o);
    const long long n= o.n > 0 ? o.n : protocol::kNScalar;
    const std::string path= joinPath(phaseDir(o, "p5"), "p5_" + o.tag + ".jsonl");
    JsonlWriter w(path);
    if (o.only_precision != "float")
        p5ForPrecision<double>(o, rc, w, kappas, seeds, n);
    if (o.only_precision != "double")
        p5ForPrecision<float>(o, rc, w, kappas, seeds, n);
    const std::string out= w.finish(rc, "p5");
    std::fprintf(stderr, "exp7 P5: wrote %lld rows to %s\n", w.lines(), out.c_str());
    return 0;
}

// ===========================================================================
// P6 - performance
// ===========================================================================
struct BenchCase
{
    const char *id;
    bool is_float;
    double kappa;
    double cap;
};

static std::vector<BenchCase> benchTable()
{
    const double inf= std::numeric_limits<double>::infinity();
    std::vector<BenchCase> v;
    BenchCase b;
    b.id= "B0"; b.is_float= false; b.kappa= 0.51; b.cap= inf;  v.push_back(b);
    b.id= "B1"; b.is_float= true;  b.kappa= 0.55; b.cap= inf;  v.push_back(b);
    b.id= "B2"; b.is_float= false; b.kappa= 0.75; b.cap= inf;  v.push_back(b);
    b.id= "B3"; b.is_float= false; b.kappa= 0.51; b.cap= 5.0;  v.push_back(b);
    b.id= "B4"; b.is_float= false; b.kappa= 0.51; b.cap= 20.0; v.push_back(b);
    return v;
}

/// Seeds the randomization of the *method order within a block*.  It is not a sampling
/// seed: it never reaches a variate, and the order it produces is recorded in every row, so
/// the randomization is auditable rather than merely asserted.  It is written out here,
/// like every other constant in this experiment, rather than derived from anything.
static const unsigned kBenchOrderRandomization= 917u;

struct BlockResult
{
    double seconds;
    long long returned;
    long long finite;
    long long nonfinite;
    unsigned long long attempts_reported;
};

template <typename T>
static BlockResult benchBlock(int method, double kappa, double cap, unsigned seed,
                              long long n_return)
{
    const double zhat[3]= {0.0, 0.0, 1.0};
    BlockResult br;
    br.returned= br.finite= br.nonfinite= 0;
    br.attempts_reported= 0;

    if (method == kMethodLegacy)
    {
        typedef typename bikappa_v1::bi_kappa_distribution<T>::point_type P;
        P ubt;
        ubt[0]= T(0); ubt[1]= T(0); ubt[2]= T(1);
        const T capT= (cap == std::numeric_limits<double>::infinity())
                          ? bikappa_v1::bi_kappa_distribution<T>::no_cap()
                          : static_cast<T>(cap);
        bikappa_v1::bi_kappa_distribution<T> dist(static_cast<T>(kappa), T(1), T(1), ubt,
                                                  capT);
        dist.seed(static_cast<int>(seed));
        const std::chrono::steady_clock::time_point t0= std::chrono::steady_clock::now();
        for (long long i= 0; i < n_return; ++i)
        {
            const P v= dist();
            if (std::isfinite(v[0]) && std::isfinite(v[1]) && std::isfinite(v[2]))
                ++br.finite;
            else
                ++br.nonfinite;
            ++br.returned;
        }
        br.seconds=
            std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    }
    else
    {
        typedef typename bi_kappa_distribution<T>::point_type P;
        P ubt;
        ubt[0]= T(0); ubt[1]= T(0); ubt[2]= T(1);
        const T capT= (cap == std::numeric_limits<double>::infinity())
                          ? bi_kappa_distribution<T>::no_cap()
                          : static_cast<T>(cap);
        bi_kappa_distribution<T> dist(static_cast<T>(kappa), T(1), T(1), ubt, capT);
        dist.seed(static_cast<int>(seed));
        const std::chrono::steady_clock::time_point t0= std::chrono::steady_clock::now();
        for (long long i= 0; i < n_return; ++i)
        {
            const P v= dist();
            if (std::isfinite(v[0]) && std::isfinite(v[1]) && std::isfinite(v[2]))
                ++br.finite;
            else
                ++br.nonfinite;
            ++br.returned;
        }
        br.seconds=
            std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        br.attempts_reported= dist.n_attempts();
    }
    (void)zhat;
    return br;
}

template <typename T>
static void benchOne(const Options &o, const RunContext &rc, JsonlWriter &w,
                     const BenchCase &bc, std::mt19937 &order_rng)
{
    const std::vector<unsigned> bseeds= benchSeedsFor(o);
    const int kBlocks= o.smoke ? 2 : 10;
    const double kMinSeconds= o.smoke ? 0.0 : 2.0;
    const long long kMinAttempts= o.smoke ? 100 : 1000000;

    // One untimed warm-up per method, which also calibrates the block size so that every
    // timed block clears both floors.  The warm-up is never reported.
    long long n_return[2];
    for (int mi= 0; mi < 2; ++mi)
    {
        const int method= (mi == 0) ? kMethodLegacy : kMethodCandidate;
        long long guess= o.smoke ? (o.n > 0 ? o.n : 200) : 200000;
        const BlockResult warm= benchBlock<T>(method, bc.kappa, bc.cap, bseeds[0], guess);
        if (kMinSeconds > 0 && warm.seconds > 0)
        {
            double scale= kMinSeconds / warm.seconds;
            if (scale < 1.0)
                scale= 1.0;
            guess= static_cast<long long>(static_cast<double>(guess) * scale * 1.15) + 1;
        }
        if (guess < kMinAttempts)
            guess= kMinAttempts;
        n_return[mi]= guess;
    }

    for (int b= 0; b < kBlocks; ++b)
    {
        // Randomized method order, drawn from a declared generator and recorded per row, so
        // that neither method systematically occupies a warmed or a cold cache slot and the
        // realized order can be checked rather than assumed.
        const bool candidate_first= ((order_rng() & 1u) != 0u);
        const unsigned bseed= bseeds[static_cast<size_t>(b) % bseeds.size()];
        for (int k= 0; k < 2; ++k)
        {
            const bool doing_candidate= candidate_first ? (k == 0) : (k == 1);
            const int method= doing_candidate ? kMethodCandidate : kMethodLegacy;
            const int mi= doing_candidate ? 1 : 0;
            const BlockResult r=
                benchBlock<T>(method, bc.kappa, bc.cap, bseed, n_return[mi]);
            emitConfigFields(w.f(), rc, "p6",
                             bc.cap == std::numeric_limits<double>::infinity() ? "uncapped"
                                                                               : "capped",
                             methodName(method), precisionId<T>(), bc.kappa,
                             static_cast<double>(shapeFromKappa<T>(bc.kappa)), bseed,
                             r.returned);
            std::fprintf(w.f(), ",\"case\":\"%s\",", bc.id);
            emitCap(w.f(), bc.cap);
            std::fprintf(w.f(),
                         ",\"block\":%d,\"slot\":%d,\"candidate_first\":%s,"
                         "\"n_returned\":%lld,\"n_finite\":%lld,\"n_nonfinite\":%lld,"
                         "\"attempts_reported\":%llu,\"seconds\":%.9g,",
                         b, k, candidate_first ? "true" : "false", r.returned, r.finite,
                         r.nonfinite, r.attempts_reported, r.seconds);
            emitDouble(w.f(), "seconds_per_returned",
                       r.returned ? r.seconds / static_cast<double>(r.returned) : 0.0);
            std::fprintf(w.f(), "}\n");
            w.countLine();
        }
    }

    // Same-build, same-method, same-seed reproducibility.
    for (int mi= 0; mi < 2; ++mi)
    {
        const int method= (mi == 0) ? kMethodLegacy : kMethodCandidate;
        const long long nrep= o.smoke ? 200 : 10000;
        const BlockResult a1= benchBlock<T>(method, bc.kappa, bc.cap, bseeds[0], nrep);
        const BlockResult a2= benchBlock<T>(method, bc.kappa, bc.cap, bseeds[0], nrep);
        emitConfigFields(w.f(), rc, "p6", "reproducibility", methodName(method),
                         precisionId<T>(), bc.kappa,
                         static_cast<double>(shapeFromKappa<T>(bc.kappa)), bseeds[0],
                         a1.returned);
        std::fprintf(w.f(), ",\"case\":\"%s\",", bc.id);
        emitCap(w.f(), bc.cap);
        std::fprintf(w.f(),
                     ",\"repeat_identical\":%s,\"finite_1\":%lld,\"finite_2\":%lld,"
                     "\"nonfinite_1\":%lld,\"nonfinite_2\":%lld,"
                     "\"attempts_1\":%llu,\"attempts_2\":%llu}\n",
                     (a1.finite == a2.finite && a1.nonfinite == a2.nonfinite &&
                      a1.attempts_reported == a2.attempts_reported)
                         ? "true"
                         : "false",
                     a1.finite, a2.finite, a1.nonfinite, a2.nonfinite, a1.attempts_reported,
                     a2.attempts_reported);
        w.countLine();
    }
}

static int phaseP6(const Options &o, const RunContext &rc)
{
    const std::vector<BenchCase> cases= benchTable();
    const std::string path= joinPath(phaseDir(o, "p6"), "p6_" + o.tag + ".jsonl");
    JsonlWriter w(path);
    std::mt19937 order_rng(kBenchOrderRandomization);
    for (size_t i= 0; i < cases.size(); ++i)
    {
        if (!o.case_id.empty() && o.case_id != "all" && o.case_id != cases[i].id)
            continue;
        if (cases[i].is_float)
            benchOne<float>(o, rc, w, cases[i], order_rng);
        else
            benchOne<double>(o, rc, w, cases[i], order_rng);
    }
    const std::string out= w.finish(rc, "p6");
    std::fprintf(stderr, "exp7 P6: wrote %lld rows to %s\n", w.lines(), out.c_str());
    return 0;
}

// ===========================================================================
// selftest - run before every production phase
// ===========================================================================
static int g_fail= 0;
static void check(bool ok, const char *what)
{
    std::printf("%-6s %s\n", ok ? "PASS" : "FAIL", what);
    if (!ok)
        ++g_fail;
}

/// The LEGACY native replica against bikappa_v1::bi_kappa_distribution, component by
/// component.  Both are seeded identically and neither has drawn before the loop starts.
template <typename T>
static bool legacyReplicaMatchesClass(double kappa, double ratio, const double ub[3],
                                      unsigned seed, long long n)
{
    typedef typename bikappa_v1::bi_kappa_distribution<T>::point_type P;
    P ubt;
    ubt[0]= static_cast<T>(ub[0]);
    ubt[1]= static_cast<T>(ub[1]);
    ubt[2]= static_cast<T>(ub[2]);
    const std::array<T, 3> ubarr= {{ubt[0], ubt[1], ubt[2]}};
    Geometry<T> geom(static_cast<T>(kappa), T(1), static_cast<T>(ratio), ubarr);

    bikappa_v1::bi_kappa_distribution<T> dist(static_cast<T>(kappa), T(1),
                                              static_cast<T>(ratio), ubt,
                                              bikappa_v1::bi_kappa_distribution<T>::no_cap());
    dist.seed(static_cast<int>(seed));

    std::mt19937 gen(seed);
    LegacyNative<T> rep(kappa);
    for (long long i= 0; i < n; ++i)
    {
        const P vc= dist();
        const NativeAttempt<T> vr= rep.attempt(gen, geom);
        for (int j= 0; j < 3; ++j)
        {
            const bool both_nan= !(vc[j] == vc[j]) && !(vr.v[j] == vr.v[j]);
            if (!both_nan && !(vc[j] == vr.v[j]))
                return false;
        }
    }
    return true;
}

/// The CANDIDATE native replica against the released 2.0.0 class.
template <typename T>
static bool candidateReplicaMatchesClass(double kappa, double ratio, const double ub[3],
                                         unsigned seed, long long n)
{
    typedef typename bi_kappa_distribution<T>::point_type P;
    P ubt;
    ubt[0]= static_cast<T>(ub[0]);
    ubt[1]= static_cast<T>(ub[1]);
    ubt[2]= static_cast<T>(ub[2]);
    const std::array<T, 3> ubarr= {{ubt[0], ubt[1], ubt[2]}};
    Geometry<T> geom(static_cast<T>(kappa), T(1), static_cast<T>(ratio), ubarr);

    bi_kappa_distribution<T> dist(static_cast<T>(kappa), T(1), static_cast<T>(ratio), ubt,
                                  bi_kappa_distribution<T>::no_cap());
    dist.seed(static_cast<int>(seed));

    std::mt19937 gen(seed);
    CandidateNative<T> rep(kappa);
    long long nonfinite= 0;
    for (long long i= 0; i < n; ++i)
    {
        const P vc= dist();
        const NativeAttempt<T> vr= rep.attempt(gen, geom);
        if (!vr.finite)
            ++nonfinite;
        for (int j= 0; j < 3; ++j)
        {
            const bool both_nan= !(vc[j] == vc[j]) && !(vr.v[j] == vr.v[j]);
            if (!both_nan && !(vc[j] == vr.v[j]))
                return false;
        }
    }
    // The class's own counter must agree with the replica's observation, which is what makes
    // n_nonfinite() usable as evidence rather than as a comment.
    return dist.n_attempts() == static_cast<unsigned long long>(n) &&
           dist.n_nonfinite() == static_cast<unsigned long long>(nonfinite);
}

static int phaseSelftest(const Options &o)
{
    (void)o;
    // Fixture seeds, declared in exp7_common.H and disjoint from both production blocks.
    // Several checks below are gate predicates; running them on a holdout seed before the
    // holdout is read would be reading the holdout.
    const std::vector<unsigned> fx= selftestSeeds();
    std::printf("exp7 selftest  [%s / %s / %s]  protocol %s %s\n", compilerId(), stdlibId(),
                archId(), protocol::kProtocolVersion, protocol::kProtocolSha256);

    // 0. The digest the portability claim rests on must be SHA-256 and not something that
    //    merely looks like it.  FIPS 180-4 test vectors.
    {
        Sha256 h;
        const bool empty_ok=
            h.hex() == "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855";
        Sha256 h2;
        h2.update("abc", 3);
        const bool abc_ok=
            h2.hex() == "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad";
        Sha256 h3;
        const char *msg= "abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq";
        h3.update(msg, std::strlen(msg));
        const bool long_ok=
            h3.hex() == "248d6a61d20638b8e5c026930c3e6039a33ce45964ff2167f6ecedd419db06c1";
        check(empty_ok && abc_ok && long_ok, "SHA-256 matches the FIPS 180-4 test vectors");
    }

    // 1. Each native replica reproduces the class it stands for, exactly.  That is what
    //    licenses using a replica wherever an intermediate has to be observable.
    {
        const double zhat[3]= {0, 0, 1};
        const double rot[3]= {1.0 / std::sqrt(14.0), 2.0 / std::sqrt(14.0),
                              3.0 / std::sqrt(14.0)};
        check(legacyReplicaMatchesClass<double>(2.0, 2.0, zhat, fx[0], 20000),
              "LEGACY replica == v1 class, double, kappa=2, unrotated");
        check(legacyReplicaMatchesClass<double>(0.51, 2.0, rot, fx[1], 20000),
              "LEGACY replica == v1 class, double, kappa=0.51, rotated");
        check(legacyReplicaMatchesClass<float>(0.55, 2.0, rot, fx[2], 20000),
              "LEGACY replica == v1 class, float, kappa=0.55, rotated");
        check(candidateReplicaMatchesClass<double>(2.0, 2.0, zhat, fx[0], 20000),
              "CANDIDATE replica == 2.0.0 class, double, kappa=2, unrotated");
        check(candidateReplicaMatchesClass<double>(0.501, 2.0, rot, fx[1], 20000),
              "CANDIDATE replica == 2.0.0 class, double, kappa=0.501, rotated");
        check(candidateReplicaMatchesClass<float>(0.55, 2.0, rot, fx[2], 20000),
              "CANDIDATE replica == 2.0.0 class, float, kappa=0.55, rotated");
    }

    // 2. Accounting identity on a run guaranteed to exercise every branch, and the
    //    candidate's avoidable loss at the primary rescued double shape.
    {
        const std::array<double, 3> zhat= {{0, 0, 1}};
        Geometry<double> geom(0.501, 1.0, 1.0, zhat);
        const double a= shapeReference<double>(0.501);
        std::mt19937 gen(fx[0]);
        double ncache= 0.0;
        bool ncached= false;
        const double log_ovf= logOverflowThreshold<double>();
        MethodCounters mc[kNumMethods];
        const long long n= 200000;
        for (long long i= 0; i < n; ++i)
        {
            const SharedPrimitives<double> s= drawShared<double>(a, gen, ncache, ncached);
            const PairedResult pr= runPairedAttempt<double>(0.501, s, geom, log_ovf);
            for (int m= 0; m < kNumMethods; ++m)
                mc[m].observe(pr.cat[m], pr.finite[m], pr.rel_err[m],
                              pr.ref.max_log_component);
        }
        bool ok= true, saw_avoidable= false, saw_honest= false, saw_finite= false;
        for (int m= 0; m < kNumMethods; ++m)
        {
            ok= ok && mc[m].accountingHolds(n);
            saw_avoidable= saw_avoidable || mc[m].avoidable() > 0;
            saw_honest= saw_honest || mc[m].honest() > 0;
            saw_finite= saw_finite || mc[m].cat[kCatFinite] > 0;
        }
        check(ok, "accounting identity N = finite + avoidable + honest holds for every method");
        check(saw_avoidable && saw_honest && saw_finite,
              "error paths covered: avoidable, honest and finite all observed");
        char msg[220];
        std::snprintf(msg, sizeof(msg),
                      "CANDIDATE has zero avoidable loss at kappa=0.501 (double): "
                      "%lld avoidable, of which %lld finite-but-wrong, against LEGACY's %lld",
                      mc[kMethodCandidate].avoidable(),
                      mc[kMethodCandidate].cat[kCatFiniteButWrong],
                      mc[kMethodLegacy].avoidable());
        check(mc[kMethodCandidate].avoidable() == 0, msg);
        check(mc[kMethodLegacy].avoidable() > 0 &&
                  mc[kMethodQF].avoidable() >= mc[kMethodLegacy].avoidable(),
              "QF loses at least as much as LEGACY");
    }

    // 3. Deterministic rerun: the same seed gives byte-identical output and counters.
    {
        const std::array<double, 3> zhat= {{0, 0, 1}};
        Geometry<double> geom(0.51, 1.0, 1.0, zhat);
        StreamDigest d1, d2;
        for (int pass= 0; pass < 2; ++pass)
        {
            std::mt19937 gen(fx[0]);
            CandidateNative<double> rep(0.51);
            StreamDigest &d= pass ? d2 : d1;
            for (long long i= 0; i < 50000; ++i)
            {
                const NativeAttempt<double> at= rep.attempt(gen, geom);
                for (int j= 0; j < 3; ++j)
                    d.real<double>(at.v[j]);
            }
        }
        check(d1.sha256() == d2.sha256() && d1.fnv1a() == d2.fnv1a(),
              "deterministic rerun: identical digest from the same seed");
    }

    // 4. seed(int) clears the sampler's cached normal variate.  Before 2.0.0 it did not, so
    //    seed(s) followed by an odd number of draws did not reproduce; the fix is a
    //    behaviour change and is tested as one.
    {
        bi_kappa_distribution<double> dist(2.0, 1.0, 2.0);
        dist.seed(static_cast<int>(fx[3]));
        const std::array<double, 3> a1= dist();
        const std::array<double, 3> a2= dist();
        const std::array<double, 3> a3= dist();
        (void)a2;
        dist.seed(static_cast<int>(fx[3]));
        const std::array<double, 3> b1= dist();
        const std::array<double, 3> b2= dist();
        const std::array<double, 3> b3= dist();
        (void)b2;
        bool same= true;
        for (int j= 0; j < 3; ++j)
            same= same && a1[j] == b1[j] && a3[j] == b3[j];
        // And after an odd number of draws, which is the case 1.x got wrong.
        dist.seed(static_cast<int>(fx[4]));
        const std::array<double, 3> c1= dist();
        dist.seed(static_cast<int>(fx[4]));
        const std::array<double, 3> d1= dist();
        for (int j= 0; j < 3; ++j)
            same= same && c1[j] == d1[j];
        check(same, "seed(int) calls reset(), so a reseeded stream reproduces");
    }

    // 5. The two formations agree where neither can fail.
    {
        const double a= shapeReference<double>(2.0);
        std::mt19937 gen(fx[0]);
        double ncache= 0.0;
        bool ncached= false;
        double worst= 0.0;
        for (int i= 0; i < 200000; ++i)
        {
            const SharedPrimitives<double> s= drawShared<double>(a, gen, ncache, ncached);
            const double r_legacy= std::sqrt(s.x1) / std::sqrt(s.x2);
            const double log_r= 0.5 * (s.log_x1 - s.log_x2);
            const double rel= std::fabs(std::expm1(log_r - std::log(r_legacy)));
            if (rel > worst)
                worst= rel;
        }
        char msg[180];
        std::snprintf(msg, sizeof(msg),
                      "LEGACY and log-domain radii agree at kappa=2 (max rel. diff %.3g, "
                      "threshold %.3g)",
                      worst, maxRelErrorFor<double>());
        check(worst < maxRelErrorFor<double>(), msg);
    }

    // 6. The overflow threshold classifies correctly on either side of the boundary.
    {
        const double lo= logOverflowThreshold<double>();
        const double lof= logOverflowThreshold<float>();
        check(std::exp(lo - 1e-9) <= std::numeric_limits<double>::max() * 1.0000001 &&
                  !std::isfinite(static_cast<double>(std::exp(lo + 1.0))),
              "double overflow threshold brackets max()");
        check(static_cast<float>(std::exp(lof - 1.0)) <= std::numeric_limits<float>::max() &&
                  !std::isfinite(static_cast<float>(std::exp(lof + 1.0))),
              "float overflow threshold brackets max()");
    }

    // 7. The log-domain cap predicate accepts exactly what the linear predicate accepts, on
    //    draws where the linear one can be evaluated at all.
    {
        const std::array<double, 3> ubr= {{1.0 / std::sqrt(14.0), 2.0 / std::sqrt(14.0),
                                           3.0 / std::sqrt(14.0)}};
        Geometry<double> geom(0.51, 1.0, 2.0, ubr);
        const double a= shapeReference<double>(0.51);
        std::mt19937 gen(fx[3]);
        double ncache= 0.0;
        bool ncached= false;
        long long disagree= 0, comparable= 0;
        for (int i= 0; i < 200000; ++i)
        {
            const SharedPrimitives<double> s= drawShared<double>(a, gen, ncache, ncached);
            if (!(s.x2 > 0))
                continue;
            const double r= std::sqrt(s.x1) / std::sqrt(s.x2);
            if (!std::isfinite(r))
                continue;
            const std::array<double, 3> n= Geometry<double>::direction(s.cos_theta, s.phi);
            const double p0= geom.sqrt_kappa * geom.theta_perp * (r * n[0]);
            const double p1= geom.sqrt_kappa * geom.theta_perp * (r * n[1]);
            const double p2= geom.sqrt_kappa * geom.theta_par * (r * n[2]);
            if (!std::isfinite(p0) || !std::isfinite(p1) || !std::isfinite(p2))
                continue;
            ++comparable;
            const bool direct= legacyCapAccept<double>(r, s.cos_theta, s.phi, geom, 5.0);
            bool unrep= false;
            const bool logged= candidateCapAcceptLog<double>(std::log(r), s.cos_theta, s.phi,
                                                             geom, std::log(5.0), &unrep);
            if (direct != logged)
                ++disagree;
        }
        char msg[200];
        std::snprintf(msg, sizeof(msg),
                      "log-domain cap predicate matches the linear predicate "
                      "(%lld disagreements in %lld comparable draws)",
                      disagree, comparable);
        // A handful of boundary disagreements are expected from the log/exp round trip;
        // anything above one in 10^4 means the predicate is wrong, not merely rounded.
        check(comparable > 0 && disagree * 10000 <= comparable, msg);
    }

    // 8. Rotation can rescue a draw whose local vector overflows.  If this never fires, the
    //    rotated configurations cannot show what they are meant to show.
    {
        const std::array<double, 3> ubr= {{1.0 / std::sqrt(14.0), 2.0 / std::sqrt(14.0),
                                           3.0 / std::sqrt(14.0)}};
        Geometry<double> geom(0.505, 1.0, 2.0, ubr);
        const double a= shapeReference<double>(0.505);
        std::mt19937 gen(fx[4]);
        double ncache= 0.0;
        bool ncached= false;
        const double log_ovf= logOverflowThreshold<double>();
        long long rot= 0;
        const long long nrot= 4000000;
        for (long long i= 0; i < nrot; ++i)
        {
            const SharedPrimitives<double> s= drawShared<double>(a, gen, ncache, ncached);
            const PairedResult pr= runPairedAttempt<double>(0.505, s, geom, log_ovf);
            if (pr.flags & kFlagRotationRecoverable)
                ++rot;
        }
        char msg[180];
        std::snprintf(msg, sizeof(msg),
                      "rotation-recoverable draws observed (%lld in 4e6 at kappa=0.505)", rot);
        check(rot > 0, msg);
    }

    // 9. The open-interval uniform never reaches an endpoint, in either precision.  This is
    //    the property that lets log(U) be taken with neither a redraw nor a clamp.
    {
        std::mt19937 gen(fx[0]);
        bool ok_d= true, ok_f= true;
        for (int i= 0; i < 500000; ++i)
        {
            const double u= bikappa_detail::openCanonical<double, std::mt19937>(gen);
            if (!(u > 0.0 && u < 1.0))
                ok_d= false;
        }
        for (int i= 0; i < 500000; ++i)
        {
            const float u= bikappa_detail::openCanonical<float, std::mt19937>(gen);
            if (!(u > 0.0f && u < 1.0f))
                ok_f= false;
        }
        check(ok_d && ok_f, "openCanonical stays strictly inside (0,1) for float and double");
    }

    // 10. Error paths: every argument the loaders reject must throw, in both versions.
    {
        int thrown= 0;
        try { bi_kappa_distribution<double> d(0.5); (void)d; } catch (const std::invalid_argument &) { ++thrown; }
        try { bi_kappa_distribution<double> d(2.0, 0.0); (void)d; } catch (const std::invalid_argument &) { ++thrown; }
        try { bi_kappa_distribution<double> d(2.0, 1.0, 1.0, {{0, 0, 1}}, 0.0); (void)d; } catch (const std::invalid_argument &) { ++thrown; }
        try {
            bi_kappa_distribution<double> d(2.0, 1.0, 1.0, {{0, 0, 0}});
            d();
        } catch (const std::invalid_argument &) { ++thrown; }
        try { bikappa_v1::bi_kappa_distribution<double> d(0.5); (void)d; } catch (const std::invalid_argument &) { ++thrown; }
        check(thrown == 5, "error paths: kappa<=1/2, theta<=0, cap<=0 and a zero ub all throw");
    }

    // 11. The two-piece log I_W evaluator is continuous across its frozen switch and gives a
    //     uniform I_W, which is what makes Z exponential.  A coarse but independent check;
    //     the arbitrary-precision comparison under gate G0 is the analysis's.
    {
        double worst_jump= 0.0;
        for (int i= 0; i < 6; ++i)
        {
            const double a= std::ldexp(1.0, -4 * i - 1);
            // Straddle the switch closely enough that the true variation of log I across
            // the interval (its derivative is `a`, so at most 3e-12 here) cannot be mistaken
            // for a step between the two formulas.
            const double lo= logRegularizedIncompleteBeta(protocol::kLogQSwitch - 1e-12, a, 1.5);
            const double hi= logRegularizedIncompleteBeta(protocol::kLogQSwitch + 1e-12, a, 1.5);
            const double jump= std::fabs(hi - lo);
            if (jump > worst_jump)
                worst_jump= jump;
        }
        // Z of the median draw must be near log 2 for any shape, since I_W is uniform.
        std::mt19937 gen(fx[1]);
        double ncache= 0.0;
        bool ncached= false;
        bool mean_ok= true;
        double worst_mean= 0.0;
        const double shapes[]= {1e-4, 1e-2, 0.25, 1.5};
        for (size_t si= 0; si < sizeof(shapes) / sizeof(shapes[0]); ++si)
        {
            const double a= shapes[si];
            const long long n= 40000;
            double sum= 0.0;
            for (long long i= 0; i < n; ++i)
            {
                const SharedPrimitives<double> s= drawShared<double>(a, gen, ncache, ncached);
                const double log_r= 0.5 * (s.log_x1 - s.log_x2);
                const double z= zFromLogR(log_r, a);
                if (std::isfinite(z))
                    sum+= z;
            }
            const double mean= sum / static_cast<double>(n);
            const double zscore= (mean - 1.0) / (1.0 / std::sqrt(static_cast<double>(n)));
            if (std::fabs(zscore) > std::fabs(worst_mean))
                worst_mean= zscore;
            if (std::fabs(zscore) > 5.0)
                mean_ok= false;
        }
        char msg[220];
        std::snprintf(msg, sizeof(msg),
                      "Z = -log I_W(a,3/2) is Exp(1): E[Z] = 1 (worst z = %.2f), evaluator "
                      "continuous at the switch (worst jump %.3g)",
                      worst_mean, worst_jump);
        check(mean_ok && worst_jump < 1e-9, msg);
    }

    // 12. The declared matrix agrees with the frozen protocol.  checkAgainstProtocol() exits
    //     on failure, so reaching this line is the result.
    {
        checkAgainstProtocol();
        check(true, "declared ladder and seed blocks agree with config/protocol.json");
    }

    std::printf("%s: %d failure(s)\n", g_fail ? "SELFTEST FAILED" : "selftest ok", g_fail);
    return g_fail ? 1 : 0;
}

// ===========================================================================
// env - the floating-point environment the run actually executes in.
//
// Recorded rather than assumed: subnormal flushing (FTZ/DAZ on x86, the FPCR FZ bit on
// AArch64) silently turns the small-shape denominator into an exact zero earlier than IEEE
// arithmetic would, which is the mechanism under study.
// ===========================================================================
static int phaseEnv(const Options &o)
{
    (void)o;
    volatile double dsub= std::numeric_limits<double>::denorm_min();
    volatile float fsub= std::numeric_limits<float>::denorm_min();
    volatile double dprod= dsub * 1.0;
    volatile float fprod= fsub * 1.0f;
    const bool d_flushes= !(dprod == dsub) || dprod == 0.0;
    const bool f_flushes= !(fprod == fsub) || fprod == 0.0f;

    const char *round= "unknown";
    switch (std::fegetround())
    {
    case FE_TONEAREST: round= "FE_TONEAREST"; break;
    case FE_DOWNWARD: round= "FE_DOWNWARD"; break;
    case FE_UPWARD: round= "FE_UPWARD"; break;
    case FE_TOWARDZERO: round= "FE_TOWARDZERO"; break;
    default: break;
    }

    std::printf("{");
    emitEnvFields(stdout, "n/a");
    std::printf(",\"bi_kappa_version\":\"%s\",\"protocol_version\":\"%s\","
                "\"protocol_sha256\":\"%s\",\"cxxflags\":\"%s\",",
                BI_KAPPA_VERSION_STRING, protocol::kProtocolVersion,
                protocol::kProtocolSha256, cxxFlagsId());
    std::printf("\"rounding_mode\":\"%s\",\"flt_eval_method\":%d,"
                "\"double_subnormals_flushed\":%s,\"float_subnormals_flushed\":%s,"
                "\"sizeof_long_double\":%zu,\"fp_contract\":\"%s\",",
                round, static_cast<int>(FLT_EVAL_METHOD), d_flushes ? "true" : "false",
                f_flushes ? "true" : "false", sizeof(long double),
#if defined(__FP_FAST_FMA) || defined(FP_FAST_FMA)
                "fma available; built with -ffp-contract=off"
#else
                "built with -ffp-contract=off"
#endif
    );
    emitDouble(stdout, "log_overflow_threshold_double", logOverflowThreshold<double>());
    std::printf(",");
    emitDouble(stdout, "log_overflow_threshold_float", logOverflowThreshold<float>());
    std::printf(",");
    emitDouble(stdout, "log_max_double", logMaxFinite<double>());
    std::printf(",");
    emitDouble(stdout, "log_max_float", logMaxFinite<float>());
    std::printf(",");
    emitDouble(stdout, "max_rel_error_double", protocol::kMaxRelErrorDouble);
    std::printf(",");
    emitDouble(stdout, "max_rel_error_float", protocol::kMaxRelErrorFloat);
    std::printf(",\"mt19937_min\":%llu,\"mt19937_max\":%llu}\n",
                static_cast<unsigned long long>(std::mt19937::min()),
                static_cast<unsigned long long>(std::mt19937::max()));
    return 0;
}

// ===========================================================================
static void makeDirs(const Options &o, const char *phase)
{
    const std::string d= phaseDir(o, phase);
    std::string cmd= "mkdir -p '" + d + "'";
    if (std::system(cmd.c_str()) != 0)
    {
        std::fprintf(stderr, "exp7: cannot create %s\n", d.c_str());
        std::exit(1);
    }
}

int main(int argc, char **argv)
{
    if (argc < 2)
    {
        std::fprintf(stderr,
                     "usage: %s <phase> [--n N] [--out DIR] [--tag TAG] [--seeds K]\n"
                     "                  [--precision float|double] [--case Cx|Bx|all]\n"
                     "                  [--execution native|translated] [--smoke]\n"
                     "phases: p1 p2 p3 p4 p5 p6 selftest env\n",
                     argv[0]);
        return 2;
    }
    const std::string phase= argv[1];
    const Options o= parseOptions(argc, argv, 2);

    checkAgainstProtocol();

    if (phase == "selftest")
        return phaseSelftest(o);
    if (phase == "env")
        return phaseEnv(o);

    RunContext rc;
    rc.manifest_path= manifestPath(o);
    rc.run_utc= nowUtcIso();
    rc.probe_path= argv[0];
    rc.probe_sha256= Sha256::file(argv[0]);
    rc.tag= o.tag;
    rc.execution= o.execution;
    rc.smoke= o.smoke;
    if (rc.probe_sha256.empty())
    {
        std::fprintf(stderr, "exp7: cannot hash the running probe at %s; refusing to write "
                             "raw data whose provenance cannot be recorded\n",
                     argv[0]);
        return 1;
    }

    if (phase == "p1")
    {
        makeDirs(o, "p1");
        return phaseP1(o, rc);
    }
    if (phase == "p2")
    {
        makeDirs(o, "p2");
        return phaseP2(o, rc);
    }
    if (phase == "p3")
    {
        makeDirs(o, "p3");
        return phaseP3(o, rc);
    }
    if (phase == "p4")
    {
        makeDirs(o, "p4");
        return phaseP4(o, rc);
    }
    if (phase == "p5")
    {
        makeDirs(o, "p5");
        return phaseP5(o, rc);
    }
    if (phase == "p6")
    {
        makeDirs(o, "p6");
        return phaseP6(o, rc);
    }

    std::fprintf(stderr, "exp7: unknown phase %s\n", phase.c_str());
    return 2;
}
