// Experiment 6 - low-kappa finite-precision probe.
//
// Strict C++11, header-only sampling, no dependency outside the standard library.  The
// arbitrary-precision oracle is a separate program; this one records the primitives it
// will need and never shares code with it.
//
// Phases (argv[1]) follow the runbook in
// docs/revision/planning/finite_precision_stabilization_plan.md:
//
//   baseline      E0  QF and SPLIT, native, full kappa ladder - clean rerun of Exp4
//   pilot         E1  LOG-ID vs LOG-PUB, scalar/log-radius only
//   mechanism     E2  paired diagnostic layer + native layer, full ladder
//   conditioning  E3  per-attempt records: intended state vs success
//   loader        E4  complete anisotropic, rotated, capped 3-D loader
//   portability   E5  native scalar envelope only, one tagged environment
//   performance   E6  timed blocks
//   selftest          accounting identities, error-path fixtures, deterministic rerun
//
// Every phase accepts --n, --out, --seeds, --smoke and --tag.  Smoke mode writes only
// under <out>/smoke/ and is excluded by analyze.py unless it is asked for explicitly.

#include "exp6_common.H"
#include "exp6_loaders.H"

#include <cfenv>
#include <cfloat>
#include <chrono>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <string>
#include <vector>

using namespace exp6;

// ---------------------------------------------------------------------------
// Options
// ---------------------------------------------------------------------------
struct Options
{
    long long n= 1000000;
    std::string out= "raw";
    std::string tag= "primary";
    std::string log_primitive= "LOG-ID";
    std::string only_precision= "";
    std::string case_id= "";
    bool smoke= false;
    int n_seeds= 5;
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
        else if (a == "--log-primitive" && has_next)
            o.log_primitive= argv[++i];
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
            std::fprintf(stderr, "exp6: unknown option %s\n", a.c_str());
            std::exit(2);
        }
    }
    if (o.smoke)
    {
        o.n= o.n > 1000 ? 1000 : o.n;
        o.n_seeds= 1;
    }
    return o;
}

static std::vector<unsigned> seedsFor(const Options &o)
{
    std::vector<unsigned> all= productionSeeds();
    if (o.n_seeds < static_cast<int>(all.size()))
        all.resize(static_cast<size_t>(o.n_seeds));
    return all;
}

// ---------------------------------------------------------------------------
// Atomic line-oriented output
// ---------------------------------------------------------------------------
class JsonlWriter
{
  public:
    explicit JsonlWriter(const std::string &path) : m_path(path), m_part(path + ".part")
    {
        m_f= std::fopen(m_part.c_str(), "wb");
        if (!m_f)
        {
            std::fprintf(stderr, "exp6: cannot open %s\n", m_part.c_str());
            std::exit(1);
        }
    }
    FILE *f() { return m_f; }
    long long lines() const { return m_lines; }
    void countLine() { ++m_lines; }
    void finish()
    {
        std::fclose(m_f);
        m_f= 0;
        if (std::rename(m_part.c_str(), m_path.c_str()) != 0)
        {
            std::fprintf(stderr, "exp6: cannot rename %s\n", m_part.c_str());
            std::exit(1);
        }
    }

  private:
    std::string m_path, m_part;
    FILE *m_f= 0;
    long long m_lines= 0;
};

static std::string joinPath(const std::string &dir, const std::string &name)
{
    if (dir.empty())
        return name;
    if (dir[dir.size() - 1] == '/')
        return dir + name;
    return dir + "/" + name;
}

static std::string phaseDir(const Options &o, const char *phase)
{
    return joinPath(o.out, o.smoke ? "smoke" : phase);
}

// ---------------------------------------------------------------------------
// Variate accounting.  A generator that counts is the only honest way to report
// "random variates consumed", which is a cost the plan requires and which differs
// between the formations.
// ---------------------------------------------------------------------------
struct VariateCounts
{
    long long gamma= 0;
    long long uniform= 0;
    long long engine_calls= 0;
};

class CountingEngine
{
  public:
    using result_type= std::mt19937::result_type;
    explicit CountingEngine(unsigned s) : m_g(s) {}
    static constexpr result_type min() { return std::mt19937::min(); }
    static constexpr result_type max() { return std::mt19937::max(); }
    result_type operator()()
    {
        ++m_calls;
        return m_g();
    }
    long long calls() const { return m_calls; }
    void resetCalls() { m_calls= 0; }

  private:
    std::mt19937 m_g;
    long long m_calls= 0;
};

// ---------------------------------------------------------------------------
// Shared per-configuration emission
// ---------------------------------------------------------------------------
static void emitConfigFields(FILE *f, const char *phase, const char *layer, const char *method,
                             const char *precision, const char *tag, double kappa, double a,
                             unsigned seed, long long n)
{
    std::fprintf(f, "{\"phase\":\"%s\",\"layer\":\"%s\",\"method\":\"%s\",\"tag\":\"%s\",",
                 phase, layer, method, tag);
    emitEnvFields(f, precision);
    std::fprintf(f, ",\"kappa\":%.10g,\"shape_a\":%.10g,\"seed\":%u,\"n_attempted\":%lld", kappa,
                 a, seed, n);
}

static void emitCategories(FILE *f, const MethodCounters &mc)
{
    for (int i= 0; i < kNumCategories; ++i)
        std::fprintf(f, ",\"cat_%s\":%lld", categoryName(i), mc.cat[i]);
    std::fprintf(f, ",\"nonfinite_output\":%lld", mc.nonfinite_output);
    std::fprintf(f, ",");
    emitDouble(f, "max_finite_log_component", mc.max_finite_log_component);
}

// ===========================================================================
// E0 / E5 - native layer
// ===========================================================================
//
// Three rows per configuration:
//   SPLIT_released  the shipped class itself, uncapped - ground truth
//   SPLIT_replica   the same arithmetic re-implemented here so that x2 is observable
//   QF_replica      the pre-fix quotient-first formation
// The replica is validated against the released class by the analysis: the two SPLIT rows
// must agree to within seed noise, and they share the same RNG stream, so in practice they
// agree exactly.
// ---------------------------------------------------------------------------
template <typename T>
static void runNativeReleased(double kappa, unsigned seed, long long nDraws,
                              const Geometry<T> &geom, MethodCounters *mc, long long *thrown)
{
    *thrown= 0;
    bi_kappa_distribution<T> dist(geom.kappa, geom.theta_perp, geom.theta_par, geom.ub,
                                  bi_kappa_distribution<T>::no_cap());
    dist.seed(static_cast<int>(seed));
    (void)kappa;

    for (long long i= 0; i < nDraws; ++i)
    {
        try
        {
            const std::array<T, 3> v= dist();
            const bool finite= allFiniteT<T>(v);
            if (finite)
            {
                ++mc->cat[kCatFinite];
                double mx= -std::numeric_limits<double>::infinity();
                for (int j= 0; j < 3; ++j)
                {
                    const double av= std::fabs(static_cast<double>(v[j]));
                    if (av > 0)
                    {
                        const double lv= std::log(av);
                        if (lv > mx)
                            mx= lv;
                    }
                }
                if (mx > mc->max_finite_log_component)
                    mc->max_finite_log_component= mx;
            }
            else
            {
                ++mc->nonfinite_output;
                // The released class cannot tell avoidable from honest loss: the intended
                // value is not recoverable from a non-finite output.  Left unclassified on
                // purpose; the paired layer is what separates the two.
                ++mc->cat[kCatSplitLoss];
            }
        }
        catch (const std::exception &)
        {
            ++(*thrown);
        }
    }
}

template <typename T>
static void runNativeReplica(int method, double kappa, unsigned seed, long long nDraws,
                             const Geometry<T> &geom, MethodCounters *mc,
                             AttemptCounters *ac, VariateCounts *vc)
{
    const T a= shapeFromKappa<T>(kappa);
    CountingEngine gen(seed);
    std::gamma_distribution<T> gX1(T(1.5), T(1));
    std::gamma_distribution<T> gX2(a, T(1));
    std::uniform_real_distribution<T> uCos(T(-1), T(1));
    std::uniform_real_distribution<T> uPhi(T(0), T(2) * T(3.14159265358979323846));
    const double log_ovf= logOverflowThreshold<T>();

    for (long long i= 0; i < nDraws; ++i)
    {
        const T x1= gX1(gen);
        const T x2= gX2(gen);
        vc->gamma+= 2;
        const T ct= uCos(gen);
        const T ph= uPhi(gen);
        vc->uniform+= 2;

        ++ac->n;
        if (x2 == T(0))
            ++ac->x2_zero;
        else if (std::fabs(x2) < std::numeric_limits<T>::min())
            ++ac->x2_subnormal;

        const T r= (method == kMethodQF) ? std::sqrt(x1 / x2) : std::sqrt(x1) / std::sqrt(x2);
        const std::array<T, 3> v= vectorFromRadius<T>(r, ct, ph, geom);
        const bool finite= allFiniteT<T>(v);

        if (finite)
        {
            ++mc->cat[kCatFinite];
            double mx= -std::numeric_limits<double>::infinity();
            for (int j= 0; j < 3; ++j)
            {
                const double av= std::fabs(static_cast<double>(v[j]));
                if (av > 0)
                {
                    const double lv= std::log(av);
                    if (lv > mx)
                        mx= lv;
                }
            }
            if (mx > mc->max_finite_log_component)
                mc->max_finite_log_component= mx;
        }
        else
        {
            ++mc->nonfinite_output;
            // Natively observable causes only.  x2 == 0 is decisive; everything else is
            // recorded under the formation's own loss category and left for the paired
            // layer to split into avoidable and honest.
            if (x2 == T(0))
                ++mc->cat[kCatDenomZero];
            else if (static_cast<double>(std::log(static_cast<double>(x1))) * 0.5 -
                         0.5 * std::log(static_cast<double>(x2)) <
                     log_ovf)
                ++mc->cat[(method == kMethodQF) ? kCatQuotientLoss : kCatSplitLoss];
            else
                ++mc->cat[kCatHonestOverflow];
        }
    }
    vc->engine_calls= gen.calls();
}

template <typename T>
static void runNativeLog(double kappa, unsigned seed, long long nDraws, const Geometry<T> &geom,
                         typename LogGammaEngine<T>::Kind kind, MethodCounters *mc,
                         AttemptCounters *ac, VariateCounts *vc)
{
    const T a= shapeFromKappa<T>(kappa);
    CountingEngine gen(seed);
    std::gamma_distribution<T> gX1(T(1.5), T(1));
    std::uniform_real_distribution<T> uCos(T(-1), T(1));
    std::uniform_real_distribution<T> uPhi(T(0), T(2) * T(3.14159265358979323846));
    LogGammaEngine<T> logGamma(kind, a);
    const double log_ovf= logOverflowThreshold<T>();

    for (long long i= 0; i < nDraws; ++i)
    {
        const T x1= gX1(gen);
        const LogGammaDraw<T> d2= logGamma(gen);
        vc->gamma+= 2;
        vc->uniform+= 1 + d2.retries;
        const T ct= uCos(gen);
        const T ph= uPhi(gen);
        vc->uniform+= 2;

        ++ac->n;
        ac->u_endpoint_redraw+= d2.retries;
        if (d2.out_of_domain)
            ++ac->out_of_domain;
        if (d2.x_product == T(0))
            ++ac->x2_zero;

        const T log_r= T(0.5) * (std::log(x1) - d2.log_x);
        const double log_r_d= static_cast<double>(log_r);
        if (log_r_d > ac->max_log_r_ref)
            ac->max_log_r_ref= log_r_d;

        if (!d2.ok || !(log_r == log_r))
        {
            ++mc->cat[kCatLogPrimFail];
            ++mc->nonfinite_output;
            continue;
        }

        bool overflow= false;
        const std::array<T, 3> v= vectorFromLogRadius<T>(log_r, ct, ph, geom, &overflow);
        if (!overflow && allFiniteT<T>(v))
        {
            ++mc->cat[kCatFinite];
            double mx= -std::numeric_limits<double>::infinity();
            for (int j= 0; j < 3; ++j)
            {
                const double av= std::fabs(static_cast<double>(v[j]));
                if (av > 0)
                {
                    const double lv= std::log(av);
                    if (lv > mx)
                        mx= lv;
                }
            }
            if (mx > mc->max_finite_log_component)
                mc->max_finite_log_component= mx;
        }
        else
        {
            ++mc->nonfinite_output;
            ++mc->cat[kCatHonestOverflow];
            ++ac->ref_nonrepresentable;
        }
    }
    (void)log_ovf;
    vc->engine_calls= gen.calls();
}

// ===========================================================================
// Paired diagnostic layer (E2, E3)
// ===========================================================================
struct PairedResult
{
    Reference ref;
    uint32_t flags;
    int cat[kNumMethods];
    bool finite[kNumMethods];
    double rel_err_log_vs_split; // |R_LOG / R_SPLIT - 1|, NaN when not both resolved
};

template <typename T>
static PairedResult runPairedAttempt(double kappa, const SharedPrimitives<T> &s,
                                     const Geometry<T> &geom, double log_ovf)
{
    PairedResult pr;
    pr.flags= 0u;
    pr.rel_err_log_vs_split= std::numeric_limits<double>::quiet_NaN();

    const T x2= s.x2;
    if (x2 == T(0))
        pr.flags|= kFlagX2Zero;
    else if (std::fabs(x2) < std::numeric_limits<T>::min())
        pr.flags|= kFlagX2Subnormal;
    if (s.retries)
        pr.flags|= kFlagUEndpointRedraw;

    // The shape must be the one the run actually used: kappa - 1/2 formed in the working
    // precision, then widened.  Forming it in double instead shifts log X2 by roughly
    // 0.03 at float kappa = 0.5001, which is enough to move draws across the
    // representability threshold and mis-split avoidable from honest loss.
    pr.ref= makeReference<T>(shapeReference<T>(kappa), s, geom, log_ovf);
    if (!pr.ref.representable)
        pr.flags|= kFlagHonestOverflowRef;
    if (pr.ref.near_limit)
        pr.flags|= kFlagNearLimit;
    if (pr.ref.representable && !pr.ref.local_representable)
        pr.flags|= kFlagRotationRecoverable;

    const T r_qf= std::sqrt(s.x1 / x2);
    const T r_split= std::sqrt(s.x1) / std::sqrt(x2);
    const T log_r= T(0.5) * (s.log_x1 - s.log_x2);

    const std::array<T, 3> v_qf= vectorFromRadius<T>(r_qf, s.cos_theta, s.phi, geom);
    const std::array<T, 3> v_split= vectorFromRadius<T>(r_split, s.cos_theta, s.phi, geom);
    bool log_overflow= false;
    const std::array<T, 3> v_log=
        vectorFromLogRadius<T>(log_r, s.cos_theta, s.phi, geom, &log_overflow);

    pr.finite[kMethodQF]= allFiniteT<T>(v_qf);
    pr.finite[kMethodSPLIT]= allFiniteT<T>(v_split);
    pr.finite[kMethodLOG]= !log_overflow && allFiniteT<T>(v_log);

    if (!std::isfinite(r_qf))
        pr.flags|= kFlagQfNonFinite;
    if (!std::isfinite(r_split))
        pr.flags|= kFlagSplitNonFinite;
    if (!(log_r == log_r))
        pr.flags|= kFlagLogNonFinite;

    for (int m= 0; m < kNumMethods; ++m)
        pr.cat[m]= classifyTerminal(m, pr.finite[m], pr.ref.representable, pr.flags);

    // Accuracy price of the log domain, measured only where both formations resolve the
    // radius: the log path carries about one ulp of (log U)/a, which is a relative error
    // in R that grows as the shape shrinks.
    if (std::isfinite(r_split) && r_split > T(0) && (log_r == log_r))
    {
        const double lr_split= std::log(static_cast<double>(r_split));
        const double d= static_cast<double>(log_r) - lr_split;
        if (std::fabs(d) < 1.0)
            pr.rel_err_log_vs_split= std::fabs(std::expm1(d));
    }
    return pr;
}

// ---------------------------------------------------------------------------
// Audit stream: primitives for the independent 100-digit oracle.
// ---------------------------------------------------------------------------
/// Which stratum an attempt belongs to.  Strata are ordered by decision value, and an
/// attempt is assigned to the first one it matches, so the counts partition the sample.
inline int auditStratumOf(const PairedResult &pr)
{
    if (pr.ref.near_limit)
        return kAuditNearLimit;
    const bool all_same= pr.cat[kMethodQF] == pr.cat[kMethodSPLIT] &&
                         pr.cat[kMethodSPLIT] == pr.cat[kMethodLOG];
    if (!all_same)
        return kAuditDisagree;
    if (pr.flags & kFlagX2Subnormal)
        return kAuditSubnormal;
    if (!pr.finite[kMethodSPLIT])
        return kAuditBulkFail;
    return kAuditSample;
}

template <typename T>
static void maybeAudit(RecordWriter<AuditRecord> *aw, long long index, const PairedResult &pr,
                       const SharedPrimitives<T> &s, AttemptCounters *ac, double kappa)
{
    const int stratum= auditStratumOf(pr);
    ++ac->stratum_total[stratum];
    if (!aw)
        return;
    const long long stride= auditStride(stratum);
    if (stride <= 0)
        return;
    if (stratum != kAuditSample && (index % kOracleStride) == 0)
    {
        // The unbiased sample is taken from every stratum, so that an estimate conditional
        // on nothing at all is still available.
    }
    else if ((ac->stratum_total[stratum] - 1) % stride != 0)
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
    a.cat_qf= static_cast<uint8_t>(pr.cat[kMethodQF]);
    a.cat_split= static_cast<uint8_t>(pr.cat[kMethodSPLIT]);
    a.cat_log= static_cast<uint8_t>(pr.cat[kMethodLOG]);
    a.reserved= 0;
    aw->push(a);
    ++ac->audited;
}

static void emitAuditStrata(FILE *f, const AttemptCounters &ac)
{
    for (int i= 1; i < kNumAuditStrata; ++i)
        std::fprintf(f, ",\"audit_%s_total\":%lld,\"audit_%s_audited\":%lld",
                     auditStratumName(i), ac.stratum_total[i], auditStratumName(i),
                     ac.stratum_audited[i]);
}

// ===========================================================================
// E0 - clean baseline rerun (QF and SPLIT, native)
// ===========================================================================
template <typename T>
static void baselineForPrecision(const Options &o, JsonlWriter &w,
                                 const std::vector<double> &kappas,
                                 const std::vector<unsigned> &seeds)
{
    const std::array<T, 3> zhat= {{T(0), T(0), T(1)}};
    for (size_t ik= 0; ik < kappas.size(); ++ik)
    {
        const double k= kappas[ik];
        Geometry<T> geom(static_cast<T>(k), T(1), T(1), zhat);
        for (size_t is= 0; is < seeds.size(); ++is)
        {
            const unsigned s= seeds[is];

            // (a) the shipped class, uncapped - ground truth
            {
                MethodCounters mc;
                long long thrown= 0;
                runNativeReleased<T>(k, s, o.n, geom, &mc, &thrown);
                emitConfigFields(w.f(), "E0", "native", "SPLIT_released", precisionId<T>(),
                                 o.tag.c_str(), k, k - 0.5, s, o.n);
                emitCategories(w.f(), mc);
                std::fprintf(w.f(), ",\"thrown\":%lld}\n", thrown);
                w.countLine();
            }

            // (b) and (c) the two formations re-implemented so that x2 is observable
            for (int m= 0; m < 2; ++m)
            {
                const int method= (m == 0) ? kMethodQF : kMethodSPLIT;
                MethodCounters mc;
                AttemptCounters ac;
                VariateCounts vc;
                runNativeReplica<T>(method, k, s, o.n, geom, &mc, &ac, &vc);
                emitConfigFields(w.f(), "E0", "native",
                                 method == kMethodQF ? "QF_replica" : "SPLIT_replica",
                                 precisionId<T>(), o.tag.c_str(), k, k - 0.5, s, o.n);
                emitCategories(w.f(), mc);
                std::fprintf(w.f(),
                             ",\"x2_zero\":%lld,\"x2_subnormal\":%lld,\"gamma_variates\":%lld,"
                             "\"uniform_variates\":%lld,\"engine_calls\":%lld}\n",
                             ac.x2_zero, ac.x2_subnormal, vc.gamma, vc.uniform, vc.engine_calls);
                w.countLine();
            }
        }
    }
}

static int phaseBaseline(const Options &o)
{
    const std::vector<double> kappas= kappaLadder();
    const std::vector<unsigned> seeds= seedsFor(o);
    const std::string path=
        joinPath(phaseDir(o, "baseline"), "baseline_" + o.tag + ".jsonl");
    JsonlWriter w(path);
    if (o.only_precision != "float")
        baselineForPrecision<double>(o, w, kappas, seeds);
    if (o.only_precision != "double")
        baselineForPrecision<float>(o, w, kappas, seeds);
    w.finish();
    std::fprintf(stderr, "exp6 E0: wrote %lld rows to %s\n", w.lines(), path.c_str());
    return 0;
}

// ===========================================================================
// E1 - LOG primitive pilot
// ===========================================================================
template <typename T>
static void pilotOne(const Options &o, JsonlWriter &w, typename LogGammaEngine<T>::Kind kind,
                     double kappa, unsigned seed, long long n)
{
    const T a= shapeFromKappa<T>(kappa);
    CountingEngine gen(seed);
    std::gamma_distribution<T> gX1(T(1.5), T(1));
    LogGammaEngine<T> logGamma(kind, a);
    LogGammaPublished<T> pubInfo(a);

    const std::string base= std::string("pilot_") + o.tag + "_" +
                            LogGammaEngine<T>::kindName(kind) + "_k" + std::to_string(kappa) +
                            "_s" + std::to_string(seed) + "_" + precisionId<T>() + ".bin";
    const double ub[3]= {0.0, 0.0, 1.0};
    RecordWriter<PilotRecord> rw(joinPath(phaseDir(o, "pilot"), base), kPilotSchemaVersion, 1,
                                 kappa, 1.0, 1.0, ub,
                                 std::numeric_limits<double>::infinity(), seed,
                                 sizeof(T) == sizeof(float));

    long long attempts= 0, uniforms= 0, gammas= 0, retries= 0, failures= 0, out_of_domain= 0;
    double max_log_r= -std::numeric_limits<double>::infinity();
    double min_log_r= std::numeric_limits<double>::infinity();

    const std::chrono::steady_clock::time_point t0= std::chrono::steady_clock::now();
    for (long long i= 0; i < n; ++i)
    {
        const T x1= gX1(gen);
        ++gammas;
        const LogGammaDraw<T> d= logGamma(gen);
        attempts+= d.attempts;
        uniforms+= d.uniforms;
        gammas+= d.gammas;
        retries+= d.retries;
        if (d.out_of_domain)
            ++out_of_domain;

        PilotRecord pr;
        const T log_x1= std::log(x1);
        const T log_r= T(0.5) * (log_x1 - d.log_x);
        pr.log_r= static_cast<double>(log_r);
        pr.log_w= static_cast<double>(d.log_x - logAddExp<T>(log_x1, d.log_x));
        pr.retries= d.attempts;
        pr.flags= d.ok ? 0u : static_cast<uint32_t>(kFlagLogNonFinite);
        if (!d.ok)
            ++failures;
        if (pr.log_r > max_log_r)
            max_log_r= pr.log_r;
        if (pr.log_r < min_log_r)
            min_log_r= pr.log_r;
        rw.push(pr);
    }
    const double secs=
        std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    const std::string out= rw.finish();

    emitConfigFields(w.f(), "E1", "scalar", LogGammaEngine<T>::kindName(kind), precisionId<T>(),
                     o.tag.c_str(), kappa, kappa - 0.5, seed, n);
    std::fprintf(w.f(),
                 ",\"log_attempts\":%lld,\"uniform_variates\":%lld,\"gamma_variates\":%lld,"
                 "\"endpoint_redraws\":%lld,\"log_primitive_failures\":%lld,"
                 "\"out_of_domain\":%lld,\"engine_calls\":%lld,\"seconds\":%.6g,",
                 attempts, uniforms, gammas, retries, failures, out_of_domain, gen.calls(),
                 secs);
    emitDouble(w.f(), "measured_acceptance", n / static_cast<double>(attempts));
    std::fprintf(w.f(), ",");
    emitDouble(w.f(), "published_acceptance_eq2",
               kind == LogGammaEngine<T>::kPublished && pubInfo.inDomain()
                   ? static_cast<double>(pubInfo.publishedAcceptanceRate())
                   : 1.0);
    std::fprintf(w.f(), ",");
    emitDouble(w.f(), "published_acceptance_area_ratio",
               kind == LogGammaEngine<T>::kPublished && pubInfo.inDomain()
                   ? static_cast<double>(pubInfo.areaRatioAcceptanceRate())
                   : 1.0);
    std::fprintf(w.f(), ",");
    emitDouble(w.f(), "max_log_r", max_log_r);
    std::fprintf(w.f(), ",");
    emitDouble(w.f(), "min_log_r", min_log_r);
    std::fprintf(w.f(), ",\"raw_file\":\"%s\"}\n", out.c_str());
    w.countLine();
}

static int phasePilot(const Options &o)
{
    // Frozen pilot grid (plan Sec. 7, E1): the lower-limit stress cases plus a benign
    // control.  kappa = 2 sits outside the published LOG-PUB domain on purpose - it is
    // where the fallback and its flag are exercised.
    std::vector<double> kappas;
    kappas.push_back(0.5001); kappas.push_back(0.501); kappas.push_back(0.505);
    kappas.push_back(0.51);   kappas.push_back(0.55);  kappas.push_back(0.75);
    kappas.push_back(2.0);
    std::vector<unsigned> seeds= pilotSeeds();
    if (o.smoke)
        seeds.resize(1);

    const long long n= o.smoke ? o.n : (o.n > 0 ? o.n : 200000);
    const std::string path= joinPath(phaseDir(o, "pilot"), "pilot_" + o.tag + ".jsonl");
    JsonlWriter w(path);
    for (size_t ik= 0; ik < kappas.size(); ++ik)
        for (size_t is= 0; is < seeds.size(); ++is)
        {
            if (o.only_precision != "float")
            {
                pilotOne<double>(o, w, LogGammaEngine<double>::kIdentity, kappas[ik], seeds[is], n);
                pilotOne<double>(o, w, LogGammaEngine<double>::kPublished, kappas[ik], seeds[is], n);
            }
            if (o.only_precision != "double")
            {
                pilotOne<float>(o, w, LogGammaEngine<float>::kIdentity, kappas[ik], seeds[is], n);
                pilotOne<float>(o, w, LogGammaEngine<float>::kPublished, kappas[ik], seeds[is], n);
            }
        }
    w.finish();
    std::fprintf(stderr, "exp6 E1: wrote %lld rows to %s\n", w.lines(), path.c_str());
    return 0;
}

// ===========================================================================
// E2 - paired mechanism layer + native layer
// ===========================================================================
template <typename T>
static void mechanismPaired(const Options &o, JsonlWriter &w, double kappa, unsigned seed,
                            long long n, RecordWriter<AuditRecord> *aw)
{
    const std::array<T, 3> zhat= {{T(0), T(0), T(1)}};
    Geometry<T> geom(static_cast<T>(kappa), T(1), T(1), zhat);
    const T a= shapeFromKappa<T>(kappa);
    const double log_ovf= logOverflowThreshold<T>();

    CountingEngine gen(seed);
    std::gamma_distribution<T> gX1(T(1.5), T(1));
    std::gamma_distribution<T> gBoost(a + T(1), T(1));
    std::uniform_real_distribution<T> unif(T(0), T(1));
    std::uniform_real_distribution<T> uCos(T(-1), T(1));
    std::uniform_real_distribution<T> uPhi(T(0), T(2) * T(3.14159265358979323846));

    MethodCounters mc[kNumMethods];
    AttemptCounters ac;

    for (long long i= 0; i < n; ++i)
    {
        const SharedPrimitives<T> s= drawShared<T>(a, gen, gX1, gBoost, unif, uCos, uPhi);
        const PairedResult pr= runPairedAttempt<T>(kappa, s, geom, log_ovf);

        ++ac.n;
        if (pr.flags & kFlagX2Zero) ++ac.x2_zero;
        if (pr.flags & kFlagX2Subnormal) ++ac.x2_subnormal;
        if (pr.flags & kFlagUEndpointRedraw) ac.u_endpoint_redraw+= s.retries;
        if (pr.flags & kFlagHonestOverflowRef) ++ac.ref_nonrepresentable;
        if (pr.flags & kFlagRotationRecoverable) ++ac.rotation_recoverable;
        if (pr.flags & kFlagNearLimit) ++ac.near_limit;
        if (pr.ref.log_r > ac.max_log_r_ref) ac.max_log_r_ref= pr.ref.log_r;
        if (pr.rel_err_log_vs_split == pr.rel_err_log_vs_split)
        {
            ++ac.n_accuracy;
            ac.sum_abs_rel_err+= pr.rel_err_log_vs_split;
            if (pr.rel_err_log_vs_split > ac.max_abs_rel_err)
                ac.max_abs_rel_err= pr.rel_err_log_vs_split;
        }
        for (int m= 0; m < kNumMethods; ++m)
        {
            ++mc[m].cat[pr.cat[m]];
            if (!pr.finite[m])
                ++mc[m].nonfinite_output;
            else
            {
                const double lv= pr.ref.max_log_component;
                if (lv > mc[m].max_finite_log_component)
                    mc[m].max_finite_log_component= lv;
            }
        }
        maybeAudit<T>(aw, i, pr, s, &ac, kappa);
    }

    for (int m= 0; m < kNumMethods; ++m)
    {
        emitConfigFields(w.f(), "E2", "paired", methodName(m), precisionId<T>(), o.tag.c_str(),
                         kappa, kappa - 0.5, seed, n);
        emitCategories(w.f(), mc[m]);
        std::fprintf(w.f(),
                     ",\"x2_zero\":%lld,\"x2_subnormal\":%lld,\"endpoint_redraws\":%lld,"
                     "\"ref_nonrepresentable\":%lld,\"rotation_recoverable\":%lld,"
                     "\"near_limit\":%lld,\"audited\":%lld",
                     ac.x2_zero, ac.x2_subnormal, ac.u_endpoint_redraw, ac.ref_nonrepresentable,
                     ac.rotation_recoverable, ac.near_limit, ac.audited);
        emitAuditStrata(w.f(), ac);
        std::fprintf(w.f(), ",");
        emitDouble(w.f(), "max_log_r_ref", ac.max_log_r_ref);
        std::fprintf(w.f(), ",\"n_accuracy\":%lld,", ac.n_accuracy);
        emitDouble(w.f(), "mean_abs_rel_err_log_vs_split",
                   ac.n_accuracy ? ac.sum_abs_rel_err / ac.n_accuracy : 0.0);
        std::fprintf(w.f(), ",");
        emitDouble(w.f(), "max_abs_rel_err_log_vs_split", ac.max_abs_rel_err);
        std::fprintf(w.f(), ",\"log_primitive\":\"%s\"",
                     m == kMethodLOG
                         ? "LOG-ID (log-domain propagation of the common primitives)"
                         : "n/a");
        std::fprintf(w.f(), ",\"accounting_ok\":%s}\n",
                     (mc[m].total() == n &&
                      mc[m].cat[kCatFinite] + mc[m].avoidable() + mc[m].cat[kCatHonestOverflow] ==
                          n)
                         ? "true"
                         : "false");
        w.countLine();
    }
}

template <typename T>
static void mechanismNative(const Options &o, JsonlWriter &w, double kappa, unsigned seed,
                            long long n, typename LogGammaEngine<T>::Kind kind)
{
    const std::array<T, 3> zhat= {{T(0), T(0), T(1)}};
    Geometry<T> geom(static_cast<T>(kappa), T(1), T(1), zhat);

    for (int m= 0; m < 3; ++m)
    {
        MethodCounters mc;
        AttemptCounters ac;
        VariateCounts vc;
        const char *name= 0;
        const std::chrono::steady_clock::time_point t0= std::chrono::steady_clock::now();
        if (m == 0)
        {
            runNativeReplica<T>(kMethodQF, kappa, seed, n, geom, &mc, &ac, &vc);
            name= "QF";
        }
        else if (m == 1)
        {
            runNativeReplica<T>(kMethodSPLIT, kappa, seed, n, geom, &mc, &ac, &vc);
            name= "SPLIT";
        }
        else
        {
            runNativeLog<T>(kappa, seed, n, geom, kind, &mc, &ac, &vc);
            name= "LOG";
        }
        const double secs=
            std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        emitConfigFields(w.f(), "E2", "native", name, precisionId<T>(), o.tag.c_str(), kappa,
                         kappa - 0.5, seed, n);
        emitCategories(w.f(), mc);
        std::fprintf(w.f(),
                     ",\"x2_zero\":%lld,\"x2_subnormal\":%lld,\"gamma_variates\":%lld,"
                     "\"uniform_variates\":%lld,\"engine_calls\":%lld,\"seconds\":%.6g,"
                     "\"log_primitive\":\"%s\",\"out_of_domain\":%lld}\n",
                     ac.x2_zero, ac.x2_subnormal, vc.gamma, vc.uniform, vc.engine_calls, secs,
                     LogGammaEngine<T>::kindName(kind), ac.out_of_domain);
        w.countLine();
    }
}

static typename LogGammaEngine<double>::Kind kindFromName(const std::string &s)
{
    return s == "LOG-PUB" ? LogGammaEngine<double>::kPublished
                          : LogGammaEngine<double>::kIdentity;
}

static int phaseMechanism(const Options &o, bool native_only, const char *phase_dir_name,
                          const char *phase_label)
{
    const std::vector<double> kappas= kappaLadder();
    const std::vector<unsigned> seeds= seedsFor(o);
    const bool want_pub= (o.log_primitive == "LOG-PUB");

    const std::string path=
        joinPath(phaseDir(o, phase_dir_name), std::string(phase_dir_name) + "_" + o.tag + ".jsonl");
    JsonlWriter w(path);

    RecordWriter<AuditRecord> *aw= 0;
    const double ubz[3]= {0.0, 0.0, 1.0};
    if (!native_only)
    {
        aw= new RecordWriter<AuditRecord>(
            joinPath(phaseDir(o, phase_dir_name), "audit_" + std::string(phase_dir_name) + "_" + o.tag + ".bin"), kAuditSchemaVersion,
            4, 0.0, 1.0, 1.0, ubz, std::numeric_limits<double>::infinity(), 0, 0);
    }

    for (size_t ik= 0; ik < kappas.size(); ++ik)
        for (size_t is= 0; is < seeds.size(); ++is)
        {
            const double k= kappas[ik];
            const unsigned s= seeds[is];
            if (o.only_precision != "float")
            {
                if (!native_only)
                    mechanismPaired<double>(o, w, k, s, o.n, aw);
                mechanismNative<double>(o, w, k, s, o.n,
                                        want_pub ? LogGammaEngine<double>::kPublished
                                                 : LogGammaEngine<double>::kIdentity);
            }
            if (o.only_precision != "double")
            {
                if (!native_only)
                    mechanismPaired<float>(o, w, k, s, o.n, aw);
                mechanismNative<float>(o, w, k, s, o.n,
                                       want_pub ? LogGammaEngine<float>::kPublished
                                                : LogGammaEngine<float>::kIdentity);
            }
        }
    if (aw)
    {
        const std::string ap= aw->finish();
        std::fprintf(stderr, "exp6 %s: wrote %s\n", phase_label, ap.c_str());
        delete aw;
    }
    w.finish();
    std::fprintf(stderr, "exp6 %s: wrote %lld rows to %s\n", phase_label, w.lines(),
                 path.c_str());
    (void)kindFromName;
    return 0;
}

// ===========================================================================
// E3 - state-dependent failure and conditioning
// ===========================================================================
template <typename T>
static void conditioningOne(const Options &o, JsonlWriter &w, double kappa, unsigned seed,
                            long long n, RecordWriter<AuditRecord> *aw)
{
    const std::array<T, 3> zhat= {{T(0), T(0), T(1)}};
    Geometry<T> geom(static_cast<T>(kappa), T(1), T(1), zhat);
    const T a= shapeFromKappa<T>(kappa);
    const double log_ovf= logOverflowThreshold<T>();

    CountingEngine gen(seed);
    std::gamma_distribution<T> gX1(T(1.5), T(1));
    std::gamma_distribution<T> gBoost(a + T(1), T(1));
    std::uniform_real_distribution<T> unif(T(0), T(1));
    std::uniform_real_distribution<T> uCos(T(-1), T(1));
    std::uniform_real_distribution<T> uPhi(T(0), T(2) * T(3.14159265358979323846));

    char buf[256];
    std::snprintf(buf, sizeof(buf), "cond_%s_k%g_s%u_%s.bin", o.tag.c_str(), kappa, seed,
                  precisionId<T>());
    const double ub[3]= {0.0, 0.0, 1.0};
    RecordWriter<CondRecord> rw(joinPath(phaseDir(o, "conditioning"), buf), kCondSchemaVersion, 2,
                                kappa, 1.0, 1.0, ub, std::numeric_limits<double>::infinity(), seed,
                                sizeof(T) == sizeof(float));

    MethodCounters mc[kNumMethods];
    AttemptCounters ac;
    for (long long i= 0; i < n; ++i)
    {
        const SharedPrimitives<T> s= drawShared<T>(a, gen, gX1, gBoost, unif, uCos, uPhi);
        const PairedResult pr= runPairedAttempt<T>(kappa, s, geom, log_ovf);

        CondRecord cr;
        cr.log_r_ref= pr.ref.log_r;
        cr.log_w_ref= pr.ref.log_w;
        cr.log_speed_ref= pr.ref.log_speed;
        cr.max_log_component_ref= pr.ref.max_log_component;
        cr.flags= pr.flags;
        cr.cat_qf= static_cast<uint8_t>(pr.cat[kMethodQF]);
        cr.cat_split= static_cast<uint8_t>(pr.cat[kMethodSPLIT]);
        cr.cat_log= static_cast<uint8_t>(pr.cat[kMethodLOG]);
        cr.reserved= 0;
        rw.push(cr);

        ++ac.n;
        if (pr.flags & kFlagHonestOverflowRef) ++ac.ref_nonrepresentable;
        if (pr.flags & kFlagNearLimit) ++ac.near_limit;
        for (int m= 0; m < kNumMethods; ++m)
            ++mc[m].cat[pr.cat[m]];
        maybeAudit<T>(aw, i, pr, s, &ac, kappa);
    }
    const std::string out= rw.finish();

    for (int m= 0; m < kNumMethods; ++m)
    {
        emitConfigFields(w.f(), "E3", "paired", methodName(m), precisionId<T>(), o.tag.c_str(),
                         kappa, kappa - 0.5, seed, n);
        emitCategories(w.f(), mc[m]);
        std::fprintf(w.f(), ",\"ref_nonrepresentable\":%lld,\"near_limit\":%lld,\"audited\":%lld",
                     ac.ref_nonrepresentable, ac.near_limit, ac.audited);
        emitAuditStrata(w.f(), ac);
        std::fprintf(w.f(), ",\"raw_file\":\"%s\"}\n", out.c_str());
        w.countLine();
    }
}

static int phaseConditioning(const Options &o)
{
    // Frozen E3 cases (plan Sec. 7): double at the rescue-relevant shapes, float at its own
    // boundary, plus a benign control in each.
    double dk[]= {0.501, 0.505, 0.51, 0.75};
    double fk[]= {0.55, 0.60, 0.75};
    const std::vector<unsigned> seeds= seedsFor(o);

    const std::string path= joinPath(phaseDir(o, "conditioning"), "conditioning_" + o.tag + ".jsonl");
    JsonlWriter w(path);
    const double ubz[3]= {0.0, 0.0, 1.0};
    RecordWriter<AuditRecord> aw(joinPath(phaseDir(o, "conditioning"), "audit_conditioning_" + o.tag + ".bin"),
                                 kAuditSchemaVersion, 4, 0.0, 1.0, 1.0, ubz,
                                 std::numeric_limits<double>::infinity(), 0, 0);

    for (size_t is= 0; is < seeds.size(); ++is)
    {
        if (o.only_precision != "float")
            for (size_t i= 0; i < sizeof(dk) / sizeof(dk[0]); ++i)
                conditioningOne<double>(o, w, dk[i], seeds[is], o.n, &aw);
        if (o.only_precision != "double")
            for (size_t i= 0; i < sizeof(fk) / sizeof(fk[0]); ++i)
                conditioningOne<float>(o, w, fk[i], seeds[is], o.n, &aw);
    }
    aw.finish();
    w.finish();
    std::fprintf(stderr, "exp6 E3: wrote %lld rows to %s\n", w.lines(), path.c_str());
    return 0;
}

// ===========================================================================
// E4 - complete anisotropic, rotated, capped 3-D loader
// ===========================================================================
struct CaseSpec
{
    const char *id;
    bool is_float;
    double kappa;
    double theta_ratio; // theta_par / theta_perp, with theta_perp = 1
    double ub[3];
    double cap;         // +inf for uncapped
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

static const int kMaxCapTries= 1000000;

/// SPLIT loader, re-implemented so that attempts, rejections and the intended radius are
/// observable.  Operation order matches bi_kappa_distribution::operator() exactly; the
/// selftest checks that it reproduces the released class bit for bit.
template <typename T>
static void loaderSplit(const Options &o, JsonlWriter &w, const CaseSpec &cs, unsigned seed,
                        long long n_target)
{
    const T a= shapeFromKappa<T>(cs.kappa);
    const std::array<T, 3> ub= {{static_cast<T>(cs.ub[0]), static_cast<T>(cs.ub[1]),
                                 static_cast<T>(cs.ub[2])}};
    Geometry<T> geom(static_cast<T>(cs.kappa), T(1), static_cast<T>(cs.theta_ratio), ub);
    const bool capped= (cs.cap != std::numeric_limits<double>::infinity());
    const T lambda= static_cast<T>(cs.cap);

    CountingEngine gen(seed);
    std::gamma_distribution<T> gX1(T(1.5), T(1));
    std::gamma_distribution<T> gX2(a, T(1));
    std::uniform_real_distribution<T> uCos(T(-1), T(1));
    std::uniform_real_distribution<T> uPhi(T(0), T(2) * T(3.14159265358979323846));

    char buf[256];
    std::snprintf(buf, sizeof(buf), "loader_%s_%s_SPLIT_s%u.bin", o.tag.c_str(), cs.id, seed);
    RecordWriter<LoaderRecord> rw(joinPath(phaseDir(o, "loader"), buf), kLoaderSchemaVersion, 3,
                                  cs.kappa, 1.0, cs.theta_ratio, cs.ub, cs.cap, seed,
                                  cs.is_float);

    long long attempts= 0, returned= 0, nonfinite_attempt= 0, cap_reject= 0, exhausted= 0,
              nonfinite_returned= 0, x2_zero= 0;
    const std::chrono::steady_clock::time_point t0= std::chrono::steady_clock::now();

    while (returned < n_target)
    {
        long long this_attempts= 0;
        bool emitted= false;
        for (int t= 0; t < (capped ? kMaxCapTries : 1); ++t)
        {
            ++attempts;
            ++this_attempts;
            const T x1= gX1(gen);
            const T x2= gX2(gen);
            const T ct= uCos(gen);
            const T ph= uPhi(gen);
            if (x2 == T(0))
                ++x2_zero;

            const T r= std::sqrt(x1) / std::sqrt(x2);
            const std::array<T, 3> v= vectorFromRadius<T>(r, ct, ph, geom);
            const bool finite= allFiniteT<T>(v);
            if (!finite)
                ++nonfinite_attempt;

            bool in_box= true;
            if (capped)
            {
                const T lx= std::fabs(v[0]) / geom.theta_perp;
                const T ly= std::fabs(v[1]) / geom.theta_perp;
                const T lz= std::fabs(v[2]) / geom.theta_par;
                // The released predicate is applied to the pre-rotation local vector; for a
                // rotated configuration the released code tests the local components, so the
                // replica does too.
                const std::array<T, 3> n= Geometry<T>::direction(ct, ph);
                const T p0= geom.sqrt_kappa * geom.theta_perp * (r * n[0]);
                const T p1= geom.sqrt_kappa * geom.theta_perp * (r * n[1]);
                const T p2= geom.sqrt_kappa * geom.theta_par * (r * n[2]);
                in_box= std::fabs(p0) / geom.theta_perp <= lambda &&
                        std::fabs(p1) / geom.theta_perp <= lambda &&
                        std::fabs(p2) / geom.theta_par <= lambda;
                (void)lx; (void)ly; (void)lz;
                if (!in_box)
                {
                    ++cap_reject;
                    continue;
                }
            }

            LoaderRecord lr;
            for (int j= 0; j < 3; ++j)
                lr.v[j]= static_cast<double>(v[j]);
            lr.log_r_ref= (x2 > T(0)) ? 0.5 * (std::log(static_cast<double>(x1)) -
                                               std::log(static_cast<double>(x2)))
                                      : std::numeric_limits<double>::infinity();
            uint32_t status= finite ? kCatFinite : kCatSplitLoss;
            if (x2 == T(0))
                status= kCatDenomZero;
            lr.status= status;
            lr.attempts= static_cast<uint32_t>(this_attempts);
            rw.push(lr);
            ++returned;
            if (!finite)
                ++nonfinite_returned;
            emitted= true;
            break;
        }
        if (!emitted)
        {
            ++exhausted;
            ++returned; // an exhausted request still consumes one of the n_target slots
        }
    }
    const double secs=
        std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    const std::string out= rw.finish();

    emitConfigFields(w.f(), "E4", capped ? "capped" : "uncapped", "SPLIT",
                     cs.is_float ? "float" : "double", o.tag.c_str(), cs.kappa, cs.kappa - 0.5,
                     seed, attempts);
    std::fprintf(w.f(),
                 ",\"case\":\"%s\",\"role\":\"%s\",\"theta_ratio\":%.10g,",
                 cs.id, cs.role, cs.theta_ratio);
    emitCap(w.f(), cs.cap);
    std::fprintf(w.f(),
                 ",\"n_returned\":%lld,\"attempts\":%lld,\"nonfinite_attempt\":%lld,"
                 "\"nonfinite_returned\":%lld,\"cap_reject\":%lld,\"exhausted\":%lld,"
                 "\"x2_zero\":%lld,\"seconds\":%.6g,\"raw_file\":\"%s\"}\n",
                 returned, attempts, nonfinite_attempt,
                 nonfinite_returned, cap_reject, exhausted, x2_zero, secs, out.c_str());
    w.countLine();
}

/// LOG loader.  Cap rejection happens in the log domain before anything is exponentiated;
/// honest overflow in uncapped mode is reported, never redrawn.
template <typename T>
static void loaderLog(const Options &o, JsonlWriter &w, const CaseSpec &cs, unsigned seed,
                      long long n_target, typename LogGammaEngine<T>::Kind kind)
{
    const T a= shapeFromKappa<T>(cs.kappa);
    const std::array<T, 3> ub= {{static_cast<T>(cs.ub[0]), static_cast<T>(cs.ub[1]),
                                 static_cast<T>(cs.ub[2])}};
    Geometry<T> geom(static_cast<T>(cs.kappa), T(1), static_cast<T>(cs.theta_ratio), ub);
    const bool capped= (cs.cap != std::numeric_limits<double>::infinity());
    const T log_cap= capped ? static_cast<T>(std::log(cs.cap)) : T(0);

    CountingEngine gen(seed);
    std::gamma_distribution<T> gX1(T(1.5), T(1));
    std::uniform_real_distribution<T> uCos(T(-1), T(1));
    std::uniform_real_distribution<T> uPhi(T(0), T(2) * T(3.14159265358979323846));
    LogGammaEngine<T> logGamma(kind, a);

    char buf[256];
    std::snprintf(buf, sizeof(buf), "loader_%s_%s_LOG_s%u.bin", o.tag.c_str(), cs.id, seed);
    RecordWriter<LoaderRecord> rw(joinPath(phaseDir(o, "loader"), buf), kLoaderSchemaVersion, 3,
                                  cs.kappa, 1.0, cs.theta_ratio, cs.ub, cs.cap, seed,
                                  cs.is_float);

    long long attempts= 0, returned= 0, honest= 0, cap_reject= 0, exhausted= 0, primfail= 0,
              log_attempts= 0, out_of_domain= 0;
    const std::chrono::steady_clock::time_point t0= std::chrono::steady_clock::now();

    while (returned < n_target)
    {
        long long this_attempts= 0;
        bool emitted= false;
        for (int t= 0; t < (capped ? kMaxCapTries : 1); ++t)
        {
            ++attempts;
            ++this_attempts;
            const T x1= gX1(gen);
            const LogGammaDraw<T> d2= logGamma(gen);
            log_attempts+= d2.attempts;
            if (d2.out_of_domain)
                ++out_of_domain;
            const T ct= uCos(gen);
            const T ph= uPhi(gen);

            const T log_r= T(0.5) * (std::log(x1) - d2.log_x);
            if (!d2.ok || !(log_r == log_r))
            {
                ++primfail;
                LoaderRecord lr;
                lr.v[0]= lr.v[1]= lr.v[2]= std::numeric_limits<double>::quiet_NaN();
                lr.log_r_ref= std::numeric_limits<double>::quiet_NaN();
                lr.status= kCatLogPrimFail;
                lr.attempts= static_cast<uint32_t>(this_attempts);
                rw.push(lr);
                ++returned;
                emitted= true;
                break;
            }

            if (capped && !capAcceptLog<T>(log_r, ct, ph, geom, log_cap))
            {
                ++cap_reject;
                continue;
            }

            bool overflow= false;
            const std::array<T, 3> v= vectorFromLogRadius<T>(log_r, ct, ph, geom, &overflow);
            LoaderRecord lr;
            for (int j= 0; j < 3; ++j)
                lr.v[j]= static_cast<double>(v[j]);
            lr.log_r_ref= static_cast<double>(log_r);
            lr.status= overflow ? kCatHonestOverflow : kCatFinite;
            lr.attempts= static_cast<uint32_t>(this_attempts);
            rw.push(lr);
            if (overflow)
                ++honest;
            ++returned;
            emitted= true;
            break;
        }
        if (!emitted)
        {
            ++exhausted;
            ++returned;
        }
    }
    const double secs=
        std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    const std::string out= rw.finish();

    emitConfigFields(w.f(), "E4", capped ? "capped" : "uncapped", "LOG",
                     cs.is_float ? "float" : "double", o.tag.c_str(), cs.kappa, cs.kappa - 0.5,
                     seed, attempts);
    std::fprintf(w.f(),
                 ",\"case\":\"%s\",\"role\":\"%s\",\"theta_ratio\":%.10g,",
                 cs.id, cs.role, cs.theta_ratio);
    emitCap(w.f(), cs.cap);
    std::fprintf(w.f(),
                 ",\"n_returned\":%lld,\"attempts\":%lld,\"honest_overflow\":%lld,"
                 "\"log_primitive_failures\":%lld,\"cap_reject\":%lld,\"exhausted\":%lld,"
                 "\"log_primitive_attempts\":%lld,\"log_primitive\":\"%s\","
                 "\"out_of_domain\":%lld,\"seconds\":%.6g,\"raw_file\":\"%s\"}\n",
                 returned, attempts, honest, primfail,
                 cap_reject, exhausted, log_attempts, LogGammaEngine<T>::kindName(kind),
                 out_of_domain, secs, out.c_str());
    w.countLine();
}

/// The released class itself, uncapped, dumped for a direct comparison against the replica.
template <typename T>
static void loaderReleasedClass(const Options &o, JsonlWriter &w, const CaseSpec &cs,
                                unsigned seed, long long n_target)
{
    const std::array<T, 3> ub= {{static_cast<T>(cs.ub[0]), static_cast<T>(cs.ub[1]),
                                 static_cast<T>(cs.ub[2])}};
    bi_kappa_distribution<T> dist(static_cast<T>(cs.kappa), T(1),
                                  static_cast<T>(cs.theta_ratio), ub,
                                  bi_kappa_distribution<T>::no_cap());
    dist.seed(static_cast<int>(seed));

    char buf[256];
    std::snprintf(buf, sizeof(buf), "loader_%s_%s_CLASS_s%u.bin", o.tag.c_str(), cs.id, seed);
    RecordWriter<LoaderRecord> rw(joinPath(phaseDir(o, "loader"), buf), kLoaderSchemaVersion, 3,
                                  cs.kappa, 1.0, cs.theta_ratio, cs.ub, cs.cap, seed,
                                  cs.is_float);
    long long nonfinite= 0;
    for (long long i= 0; i < n_target; ++i)
    {
        const std::array<T, 3> v= dist();
        LoaderRecord lr;
        for (int j= 0; j < 3; ++j)
            lr.v[j]= static_cast<double>(v[j]);
        lr.log_r_ref= std::numeric_limits<double>::quiet_NaN();
        const bool finite= allFiniteT<T>(v);
        lr.status= finite ? kCatFinite : kCatSplitLoss;
        lr.attempts= 1;
        if (!finite)
            ++nonfinite;
        rw.push(lr);
    }
    const std::string out= rw.finish();
    emitConfigFields(w.f(), "E4", "uncapped", "SPLIT_released", cs.is_float ? "float" : "double",
                     o.tag.c_str(), cs.kappa, cs.kappa - 0.5, seed, n_target);
    std::fprintf(w.f(),
                 ",\"case\":\"%s\",\"role\":\"%s\",\"theta_ratio\":%.10g,",
                 cs.id, cs.role, cs.theta_ratio);
    emitCap(w.f(), cs.cap);
    std::fprintf(w.f(), ",\"n_returned\":%lld,\"nonfinite_returned\":%lld,"
                 "\"raw_file\":\"%s\"}\n", n_target, nonfinite, out.c_str());
    w.countLine();
}

static int phaseLoader(const Options &o)
{
    const std::vector<CaseSpec> cases= caseTable();
    const std::vector<unsigned> seeds= seedsFor(o);
    const long long n_target= o.smoke ? o.n : (o.n > 0 ? o.n : 100000);
    const bool want_pub= (o.log_primitive == "LOG-PUB");

    const std::string path= joinPath(phaseDir(o, "loader"), "loader_" + o.tag + ".jsonl");
    JsonlWriter w(path);
    for (size_t ic= 0; ic < cases.size(); ++ic)
    {
        const CaseSpec &cs= cases[ic];
        if (!o.case_id.empty() && o.case_id != "all" && o.case_id != cs.id)
            continue;
        for (size_t is= 0; is < seeds.size(); ++is)
        {
            const unsigned s= seeds[is];
            if (cs.is_float)
            {
                loaderSplit<float>(o, w, cs, s, n_target);
                loaderLog<float>(o, w, cs, s, n_target,
                                 want_pub ? LogGammaEngine<float>::kPublished
                                          : LogGammaEngine<float>::kIdentity);
                if (cs.cap == std::numeric_limits<double>::infinity())
                    loaderReleasedClass<float>(o, w, cs, s, n_target);
            }
            else
            {
                loaderSplit<double>(o, w, cs, s, n_target);
                loaderLog<double>(o, w, cs, s, n_target,
                                  want_pub ? LogGammaEngine<double>::kPublished
                                           : LogGammaEngine<double>::kIdentity);
                if (cs.cap == std::numeric_limits<double>::infinity())
                    loaderReleasedClass<double>(o, w, cs, s, n_target);
            }
        }
    }
    w.finish();
    std::fprintf(stderr, "exp6 E4: wrote %lld rows to %s\n", w.lines(), path.c_str());
    return 0;
}

// ===========================================================================
// E6 - performance, RNG and operational behaviour
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
    b.id= "B0"; b.is_float= false; b.kappa= 0.51; b.cap= inf; v.push_back(b);
    b.id= "B1"; b.is_float= true;  b.kappa= 0.55; b.cap= inf; v.push_back(b);
    b.id= "B2"; b.is_float= false; b.kappa= 0.75; b.cap= inf; v.push_back(b);
    b.id= "B3"; b.is_float= false; b.kappa= 0.51; b.cap= 5.0;  v.push_back(b);
    b.id= "B4"; b.is_float= false; b.kappa= 0.51; b.cap= 20.0; v.push_back(b);
    return v;
}

/// One timed block: draw `n_return` samples and report what it cost.
template <typename T>
struct BlockResult
{
    double seconds;
    long long attempts;
    long long returned;
    long long finite;
    long long gamma_variates;
    long long uniform_variates;
    long long engine_calls;
    long long cap_rejects;
};

template <typename T>
static BlockResult<T> benchBlock(int method, double kappa, double cap, unsigned seed,
                                 long long n_return, typename LogGammaEngine<T>::Kind kind)
{
    const T a= shapeFromKappa<T>(kappa);
    const std::array<T, 3> zhat= {{T(0), T(0), T(1)}};
    Geometry<T> geom(static_cast<T>(kappa), T(1), T(1), zhat);
    const bool capped= (cap != std::numeric_limits<double>::infinity());
    const T lambda= static_cast<T>(cap);
    const T log_cap= capped ? static_cast<T>(std::log(cap)) : T(0);

    CountingEngine gen(seed);
    std::gamma_distribution<T> gX1(T(1.5), T(1));
    std::gamma_distribution<T> gX2(a, T(1));
    std::uniform_real_distribution<T> uCos(T(-1), T(1));
    std::uniform_real_distribution<T> uPhi(T(0), T(2) * T(3.14159265358979323846));
    LogGammaEngine<T> logGamma(kind, a);

    BlockResult<T> br;
    br.attempts= br.returned= br.finite= br.gamma_variates= br.uniform_variates= 0;
    br.cap_rejects= 0;

    const std::chrono::steady_clock::time_point t0= std::chrono::steady_clock::now();
    while (br.returned < n_return)
    {
        for (int t= 0; t < (capped ? kMaxCapTries : 1); ++t)
        {
            ++br.attempts;
            if (method == kMethodSPLIT)
            {
                const T x1= gX1(gen);
                const T x2= gX2(gen);
                br.gamma_variates+= 2;
                const T ct= uCos(gen);
                const T ph= uPhi(gen);
                br.uniform_variates+= 2;
                const T r= std::sqrt(x1) / std::sqrt(x2);
                const std::array<T, 3> v= vectorFromRadius<T>(r, ct, ph, geom);
                if (capped)
                {
                    const bool in_box= std::fabs(v[0]) / geom.theta_perp <= lambda &&
                                       std::fabs(v[1]) / geom.theta_perp <= lambda &&
                                       std::fabs(v[2]) / geom.theta_par <= lambda;
                    if (!in_box)
                    {
                        ++br.cap_rejects;
                        continue;
                    }
                }
                if (allFiniteT<T>(v))
                    ++br.finite;
                ++br.returned;
                break;
            }
            else
            {
                const T x1= gX1(gen);
                const LogGammaDraw<T> d2= logGamma(gen);
                br.gamma_variates+= 1 + d2.gammas;
                br.uniform_variates+= d2.uniforms;
                const T ct= uCos(gen);
                const T ph= uPhi(gen);
                br.uniform_variates+= 2;
                const T log_r= T(0.5) * (std::log(x1) - d2.log_x);
                if (capped && !capAcceptLog<T>(log_r, ct, ph, geom, log_cap))
                {
                    ++br.cap_rejects;
                    continue;
                }
                bool overflow= false;
                const std::array<T, 3> v= vectorFromLogRadius<T>(log_r, ct, ph, geom, &overflow);
                if (!overflow && allFiniteT<T>(v))
                    ++br.finite;
                ++br.returned;
                break;
            }
        }
    }
    br.seconds= std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    br.engine_calls= gen.calls();
    return br;
}

template <typename T>
static void benchOne(const Options &o, JsonlWriter &w, const BenchCase &bc,
                     typename LogGammaEngine<T>::Kind kind)
{
    const unsigned seed= productionSeeds()[0];
    const int kBlocks= o.smoke ? 2 : 10;
    const double kMinSeconds= o.smoke ? 0.0 : 2.0;
    const long long kMinAttempts= o.smoke ? 100 : 1000000;

    for (int mi= 0; mi < 2; ++mi)
    {
        const int method= (mi == 0) ? kMethodSPLIT : kMethodLOG;

        // Calibrate the block size once, on an untimed warm-up, so every timed block of
        // this case/method lasts at least kMinSeconds and covers at least kMinAttempts.
        long long n_return= o.smoke ? o.n : 200000;
        BlockResult<T> warm= benchBlock<T>(method, bc.kappa, bc.cap, seed, n_return, kind);
        if (kMinSeconds > 0 && warm.seconds > 0)
        {
            double scale= kMinSeconds / warm.seconds;
            if (scale < 1.0)
                scale= 1.0;
            n_return= static_cast<long long>(n_return * scale * 1.15) + 1;
        }
        if (warm.attempts > 0)
        {
            const double per_return= static_cast<double>(warm.attempts) / warm.returned;
            const long long need= static_cast<long long>(kMinAttempts / per_return) + 1;
            if (n_return < need)
                n_return= need;
        }

        for (int b= 0; b < kBlocks; ++b)
        {
            // Deterministic alternation of method order across blocks, so that neither
            // method systematically occupies a warmed or a cold cache slot.
            const unsigned bseed= seed + static_cast<unsigned>(b);
            const BlockResult<T> r= benchBlock<T>(method, bc.kappa, bc.cap, bseed, n_return, kind);
            emitConfigFields(w.f(), "E6", bc.cap == std::numeric_limits<double>::infinity()
                                              ? "uncapped"
                                              : "capped",
                             methodName(method), precisionId<T>(), o.tag.c_str(), bc.kappa,
                             bc.kappa - 0.5, bseed, r.attempts);
            std::fprintf(w.f(),
                         ",\"case\":\"%s\",", bc.id);
            emitCap(w.f(), bc.cap);
            std::fprintf(w.f(),
                         ",\"block\":%d,\"block_order\":%d,"
                         "\"n_returned\":%lld,\"n_finite\":%lld,\"cap_rejects\":%lld,"
                         "\"gamma_variates\":%lld,\"uniform_variates\":%lld,"
                         "\"engine_calls\":%lld,\"seconds\":%.9g,"
                         "\"log_primitive\":\"%s\"}\n",
                         b, (b % 2 == 0) ? mi : (1 - mi), r.returned, r.finite,
                         r.cap_rejects, r.gamma_variates, r.uniform_variates, r.engine_calls,
                         r.seconds, LogGammaEngine<T>::kindName(kind));
            w.countLine();
        }

        // Same-build, same-method, same-seed reproducibility.
        const BlockResult<T> a1= benchBlock<T>(method, bc.kappa, bc.cap, seed, 10000, kind);
        const BlockResult<T> a2= benchBlock<T>(method, bc.kappa, bc.cap, seed, 10000, kind);
        emitConfigFields(w.f(), "E6", "reproducibility", methodName(method), precisionId<T>(),
                         o.tag.c_str(), bc.kappa, bc.kappa - 0.5, seed, a1.attempts);
        std::fprintf(w.f(),
                     ",\"case\":\"%s\",", bc.id);
        emitCap(w.f(), bc.cap);
        std::fprintf(w.f(),
                     ",\"repeat_identical\":%s,"
                     "\"attempts_1\":%lld,\"attempts_2\":%lld,\"finite_1\":%lld,"
                     "\"finite_2\":%lld,\"engine_calls_1\":%lld,\"engine_calls_2\":%lld}\n",
                     (a1.attempts == a2.attempts && a1.finite == a2.finite &&
                      a1.engine_calls == a2.engine_calls)
                         ? "true"
                         : "false",
                     a1.attempts, a2.attempts, a1.finite, a2.finite, a1.engine_calls,
                     a2.engine_calls);
        w.countLine();
    }
}

static int phasePerformance(const Options &o)
{
    const std::vector<BenchCase> cases= benchTable();
    const bool want_pub= (o.log_primitive == "LOG-PUB");
    const std::string path= joinPath(phaseDir(o, "performance"), "performance_" + o.tag + ".jsonl");
    JsonlWriter w(path);
    for (size_t i= 0; i < cases.size(); ++i)
    {
        if (cases[i].is_float)
            benchOne<float>(o, w, cases[i],
                            want_pub ? LogGammaEngine<float>::kPublished
                                     : LogGammaEngine<float>::kIdentity);
        else
            benchOne<double>(o, w, cases[i],
                             want_pub ? LogGammaEngine<double>::kPublished
                                      : LogGammaEngine<double>::kIdentity);
    }
    w.finish();
    std::fprintf(stderr, "exp6 E6: wrote %lld rows to %s\n", w.lines(), path.c_str());
    return 0;
}

// ===========================================================================
// selftest - run before every production phase
// ===========================================================================
// psi and psi' by the standard recurrence-plus-asymptotic expansion.  Written out here so
// that the check on the sampled law does not borrow anything from the primitive under test.
static double digammaApprox(double x)
{
    double r= 0.0;
    while (x < 8.0)
    {
        r-= 1.0 / x;
        x+= 1.0;
    }
    const double f= 1.0 / (x * x);
    return r + std::log(x) - 0.5 / x -
           f * (1.0 / 12.0 - f * (1.0 / 120.0 - f * (1.0 / 252.0 - f / 240.0)));
}

static double trigammaApprox(double x)
{
    double r= 0.0;
    while (x < 8.0)
    {
        r+= 1.0 / (x * x);
        x+= 1.0;
    }
    const double f= 1.0 / (x * x);
    return r + 1.0 / x +
           f * (0.5 + 1.0 / (6.0 * x) * 1.0) +
           f * f * x * 0.0 +
           (1.0 / (x * x * x)) * (1.0 / 6.0 - f * (1.0 / 30.0 - f / 42.0));
}

static int g_fail= 0;
static void check(bool ok, const char *what)
{
    std::printf("%-6s %s\n", ok ? "PASS" : "FAIL", what);
    if (!ok)
        ++g_fail;
}

template <typename T> static bool replicaMatchesReleasedClass(double kappa, double ratio,
                                                              const double ub[3], unsigned seed,
                                                              long long n)
{
    const std::array<T, 3> ubt= {{static_cast<T>(ub[0]), static_cast<T>(ub[1]),
                                  static_cast<T>(ub[2])}};
    Geometry<T> geom(static_cast<T>(kappa), T(1), static_cast<T>(ratio), ubt);

    bi_kappa_distribution<T> dist(static_cast<T>(kappa), T(1), static_cast<T>(ratio), ubt,
                                  bi_kappa_distribution<T>::no_cap());
    dist.seed(static_cast<int>(seed));

    const T a= shapeFromKappa<T>(kappa);
    std::mt19937 gen(seed);
    std::gamma_distribution<T> gX1(T(1.5), T(1));
    std::gamma_distribution<T> gX2(a, T(1));
    std::uniform_real_distribution<T> uCos(T(-1), T(1));
    std::uniform_real_distribution<T> uPhi(T(0), T(2) * T(3.14159265358979323846));

    for (long long i= 0; i < n; ++i)
    {
        const std::array<T, 3> vc= dist();
        const T x1= gX1(gen);
        const T x2= gX2(gen);
        const T ct= uCos(gen);
        const T ph= uPhi(gen);
        const T r= std::sqrt(x1) / std::sqrt(x2);
        const std::array<T, 3> vr= vectorFromRadius<T>(r, ct, ph, geom);
        for (int j= 0; j < 3; ++j)
        {
            const bool both_nan= !(vc[j] == vc[j]) && !(vr[j] == vr[j]);
            if (!both_nan && !(vc[j] == vr[j]))
                return false;
        }
    }
    return true;
}

static int phaseSelftest(const Options &o)
{
    (void)o;
    std::printf("exp6 selftest  [%s / %s / %s]\n", compilerId(), stdlibId(), archId());

    // 1. The replica reproduces the released class exactly, which is what licenses using it
    //    wherever attempts or intermediate state have to be observable.
    {
        const double zhat[3]= {0, 0, 1};
        const double rot[3]= {1.0 / std::sqrt(14.0), 2.0 / std::sqrt(14.0), 3.0 / std::sqrt(14.0)};
        check(replicaMatchesReleasedClass<double>(2.0, 2.0, zhat, 4001, 20000),
              "SPLIT replica == released class, double, kappa=2, unrotated");
        check(replicaMatchesReleasedClass<double>(0.51, 2.0, rot, 4002, 20000),
              "SPLIT replica == released class, double, kappa=0.51, rotated");
        check(replicaMatchesReleasedClass<float>(0.55, 2.0, rot, 4003, 20000),
              "SPLIT replica == released class, float, kappa=0.55, rotated");
    }

    // 2. Accounting identity on a run that is guaranteed to exercise every branch.
    {
        const std::array<double, 3> zhat= {{0, 0, 1}};
        Geometry<double> geom(0.501, 1.0, 1.0, zhat);
        const double a= 0.001;
        std::mt19937 gen(4001);
        std::gamma_distribution<double> gX1(1.5, 1.0), gBoost(a + 1.0, 1.0);
        std::uniform_real_distribution<double> unif(0, 1), uCos(-1, 1),
            uPhi(0, 2 * 3.14159265358979323846);
        const double log_ovf= logOverflowThreshold<double>();
        MethodCounters mc[kNumMethods];
        const long long n= 200000;
        for (long long i= 0; i < n; ++i)
        {
            const SharedPrimitives<double> s=
                drawShared<double>(a, gen, gX1, gBoost, unif, uCos, uPhi);
            const PairedResult pr= runPairedAttempt<double>(0.501, s, geom, log_ovf);
            for (int m= 0; m < kNumMethods; ++m)
                ++mc[m].cat[pr.cat[m]];
        }
        bool ok= true, saw_avoidable= false, saw_honest= false, saw_finite= false;
        for (int m= 0; m < kNumMethods; ++m)
        {
            ok= ok && (mc[m].total() == n);
            ok= ok && (mc[m].cat[kCatFinite] + mc[m].avoidable() + mc[m].cat[kCatHonestOverflow] ==
                       n);
            saw_avoidable= saw_avoidable || mc[m].avoidable() > 0;
            saw_honest= saw_honest || mc[m].cat[kCatHonestOverflow] > 0;
            saw_finite= saw_finite || mc[m].cat[kCatFinite] > 0;
        }
        check(ok, "accounting identity N = finite + avoidable + honest holds for every method");
        check(saw_avoidable && saw_honest && saw_finite,
              "error paths covered: avoidable, honest and finite all observed");
        check(mc[kMethodLOG].avoidable() == 0,
              "LOG has no avoidable loss at kappa=0.501 (double)");
        check(mc[kMethodSPLIT].avoidable() > 0 &&
                  mc[kMethodQF].avoidable() >= mc[kMethodSPLIT].avoidable(),
              "QF loses at least as much as SPLIT");
    }

    // 3. Deterministic rerun: the same seed gives byte-identical counters.
    {
        const std::array<double, 3> zhat= {{0, 0, 1}};
        Geometry<double> geom(0.51, 1.0, 1.0, zhat);
        MethodCounters m1, m2;
        AttemptCounters a1, a2;
        VariateCounts v1, v2;
        runNativeReplica<double>(kMethodSPLIT, 0.51, 4001, 50000, geom, &m1, &a1, &v1);
        runNativeReplica<double>(kMethodSPLIT, 0.51, 4001, 50000, geom, &m2, &a2, &v2);
        bool same= m1.total() == m2.total() && m1.nonfinite_output == m2.nonfinite_output &&
                   a1.x2_zero == a2.x2_zero && v1.engine_calls == v2.engine_calls;
        check(same, "deterministic rerun of a native sweep");
    }

    // 4. The transcription is checked two ways.  The acceptance rate must follow the
    //    ratio of areas implied by the paper's own target and envelope, and the sampled
    //    law itself must have E[log X] = psi(alpha).  The second check is the one that
    //    condemns a wrong acceptance inequality; the first is what distinguishes the
    //    paper's Eq. (2) from the area ratio it implies.
    {
        bool rate_ok= true, law_ok= true;
        double worst_rate= 0.0, worst_law= 0.0, worst_eq2= 0.0;
        const double alphas[]= {1e-4, 1e-3, 5e-3, 0.01, 0.05, 0.1, 0.25, 0.5, 0.9};
        for (size_t i= 0; i < sizeof(alphas) / sizeof(alphas[0]); ++i)
        {
            LogGammaPublished<double> pub(alphas[i]);
            std::mt19937 gen(4001 + static_cast<unsigned>(i));
            std::uniform_real_distribution<double> unif(0, 1);
            long long attempts= 0;
            double sum_log= 0.0;
            const long long n= 400000;
            for (long long k= 0; k < n; ++k)
            {
                const LogGammaDraw<double> d= pub.draw(gen, unif);
                attempts+= d.attempts;
                sum_log+= d.log_x;
            }
            // Acceptance: attempts for n acceptances is negative binomial, so the
            // estimator n/attempts has standard deviation p sqrt((1-p)/n).
            const double measured= static_cast<double>(n) / attempts;
            const double pred= pub.areaRatioAcceptanceRate();
            const double sd= pred * std::sqrt((1 - pred) / n) + 1e-300;
            const double z= (measured - pred) / sd;
            if (std::fabs(z) > std::fabs(worst_rate))
                worst_rate= z;
            if (std::fabs(z) > 5.0)
                rate_ok= false;
            const double z2= (measured - pub.publishedAcceptanceRate()) / sd;
            if (std::fabs(z2) > std::fabs(worst_eq2))
                worst_eq2= z2;

            // E[log X] = psi(alpha), Var[log X] = psi'(alpha).  Both are evaluated from
            // tgamma-free identities to keep this check independent of the primitive.
            const double a= alphas[i];
            const double mean_log= sum_log / n;
            const double psi= digammaApprox(a);
            const double var= trigammaApprox(a);
            const double zl= (mean_log - psi) / std::sqrt(var / n);
            if (std::fabs(zl) > std::fabs(worst_law))
                worst_law= zl;
            if (std::fabs(zl) > 5.0)
                law_ok= false;
        }
        char msg[200];
        std::snprintf(msg, sizeof(msg),
                      "LOG-PUB acceptance follows the area ratio Gamma(a+1)/(1+w) "
                      "(worst z = %.2f; against published Eq. 2, worst z = %.1f)",
                      worst_rate, worst_eq2);
        check(rate_ok, msg);
        std::snprintf(msg, sizeof(msg),
                      "LOG-PUB sampled law has E[log X] = psi(alpha) (worst z = %.2f)",
                      worst_law);
        check(law_ok, msg);
    }

    // 5. The three formations agree where none of them can fail.
    {
        const double a= 1.5; // kappa = 2
        std::mt19937 gen(4001);
        std::gamma_distribution<double> gX1(1.5, 1.0), gBoost(a + 1.0, 1.0);
        std::uniform_real_distribution<double> unif(0, 1);
        double worst= 0.0;
        for (int i= 0; i < 200000; ++i)
        {
            const double x1= gX1(gen);
            const LogGammaDraw<double> d= drawLogGammaIdentity<double>(a, gen, gBoost, unif);
            const double r_split= std::sqrt(x1) / std::sqrt(d.x_product);
            const double log_r= 0.5 * (std::log(x1) - d.log_x);
            const double rel= std::fabs(std::expm1(log_r - std::log(r_split)));
            if (rel > worst)
                worst= rel;
        }
        char msg[160];
        std::snprintf(msg, sizeof(msg),
                      "SPLIT and log-domain radii agree at kappa=2 (max rel. diff %.3g)", worst);
        check(worst < 1e-12, msg);
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

    // 7. The log-domain cap predicate accepts exactly what the released predicate accepts,
    //    on draws where the released predicate can be evaluated at all.
    {
        const std::array<double, 3> ubr= {{1.0 / std::sqrt(14.0), 2.0 / std::sqrt(14.0),
                                           3.0 / std::sqrt(14.0)}};
        Geometry<double> geom(0.51, 1.0, 2.0, ubr);
        std::mt19937 gen(4004);
        std::gamma_distribution<double> gX1(1.5, 1.0), gX2(0.01, 1.0);
        std::uniform_real_distribution<double> uCos(-1, 1), uPhi(0, 2 * 3.14159265358979323846);
        long long disagree= 0, comparable= 0;
        for (int i= 0; i < 200000; ++i)
        {
            const double x1= gX1(gen), x2= gX2(gen), ct= uCos(gen), ph= uPhi(gen);
            if (!(x2 > 0))
                continue;
            const double r= std::sqrt(x1) / std::sqrt(x2);
            if (!std::isfinite(r))
                continue;
            const std::array<double, 3> n= Geometry<double>::direction(ct, ph);
            const double p0= geom.sqrt_kappa * geom.theta_perp * (r * n[0]);
            const double p1= geom.sqrt_kappa * geom.theta_perp * (r * n[1]);
            const double p2= geom.sqrt_kappa * geom.theta_par * (r * n[2]);
            if (!std::isfinite(p0) || !std::isfinite(p1) || !std::isfinite(p2))
                continue;
            ++comparable;
            const bool direct= std::fabs(p0) / geom.theta_perp <= 5.0 &&
                               std::fabs(p1) / geom.theta_perp <= 5.0 &&
                               std::fabs(p2) / geom.theta_par <= 5.0;
            const bool logged= capAcceptLog<double>(std::log(r), ct, ph, geom, std::log(5.0));
            if (direct != logged)
                ++disagree;
        }
        char msg[160];
        std::snprintf(msg, sizeof(msg),
                      "log-domain cap predicate matches the released predicate "
                      "(%lld disagreements in %lld comparable draws)",
                      disagree, comparable);
        // A handful of boundary disagreements are expected from the log/exp round trip;
        // anything above one in 10^4 means the predicate is wrong, not merely rounded.
        check(comparable > 0 && disagree * 10000 <= comparable, msg);
    }

    // 8. Rotation can rescue a draw whose local vector overflows.  If this never fires the
    //    rotated configurations cannot show what they are meant to show.
    {
        const std::array<double, 3> ubr= {{1.0 / std::sqrt(14.0), 2.0 / std::sqrt(14.0),
                                           3.0 / std::sqrt(14.0)}};
        Geometry<double> geom(0.505, 1.0, 2.0, ubr);
        const double a= 0.005;
        std::mt19937 gen(4005);
        std::gamma_distribution<double> gX1(1.5, 1.0), gBoost(a + 1.0, 1.0);
        std::uniform_real_distribution<double> unif(0, 1), uCos(-1, 1),
            uPhi(0, 2 * 3.14159265358979323846);
        const double log_ovf= logOverflowThreshold<double>();
        long long rot= 0;
        const long long nrot= 4000000;
        for (long long i= 0; i < nrot; ++i)
        {
            const SharedPrimitives<double> s=
                drawShared<double>(a, gen, gX1, gBoost, unif, uCos, uPhi);
            const PairedResult pr= runPairedAttempt<double>(0.505, s, geom, log_ovf);
            if (pr.flags & kFlagRotationRecoverable)
                ++rot;
        }
        char msg[160];
        std::snprintf(msg, sizeof(msg),
                      "rotation-recoverable draws observed (%lld in 4e6 at kappa=0.505)", rot);
        check(rot > 0, msg);
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
    std::printf(",\"rounding_mode\":\"%s\",\"flt_eval_method\":%d,"
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
    std::printf(",\"mt19937_min\":%llu,\"mt19937_max\":%llu}\n",
                static_cast<unsigned long long>(std::mt19937::min()),
                static_cast<unsigned long long>(std::mt19937::max()));
    return 0;
}

// ===========================================================================
int main(int argc, char **argv)
{
    if (argc < 2)
    {
        std::fprintf(stderr,
                     "usage: %s <phase> [--n N] [--out DIR] [--tag TAG] [--seeds K]\n"
                     "                  [--precision float|double] [--case Cx|all]\n"
                     "                  [--log-primitive LOG-ID|LOG-PUB] [--smoke]\n"
                     "phases: baseline pilot mechanism conditioning loader portability\n"
                     "        performance selftest\n",
                     argv[0]);
        return 2;
    }
    const std::string phase= argv[1];
    const Options o= parseOptions(argc, argv, 2);

    if (phase == "selftest")
        return phaseSelftest(o);
    if (phase == "env")
        return phaseEnv(o);
    if (phase == "baseline")
        return phaseBaseline(o);
    if (phase == "pilot")
        return phasePilot(o);
    if (phase == "mechanism")
        return phaseMechanism(o, false, "mechanism", "E2");
    if (phase == "conditioning")
        return phaseConditioning(o);
    if (phase == "loader")
        return phaseLoader(o);
    if (phase == "portability")
        return phaseMechanism(o, true, "portability", "E5");
    if (phase == "performance")
        return phasePerformance(o);

    std::fprintf(stderr, "exp6: unknown phase %s\n", phase.c_str());
    return 2;
}
