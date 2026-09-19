// Experiment 6 - shared definitions for the low-kappa finite-precision study.
//
// Strict C++11.  Nothing in this header, or in anything it includes, may depend on a
// later standard: the point of the experiment is a header-only C++11 reference sampler,
// so the sampling code has to keep compiling as C++11.  The arbitrary-precision oracle
// lives in a separate translation unit (src/exp6_oracle.cpp, C++14) precisely so that it
// cannot leak a dependency in here, and so that it re-derives every classification from
// the recorded primitives rather than sharing code with the thing it audits.
//
// Three ideas carry the whole experiment and are worth stating once:
//
//  1. The radius is R = sqrt(X1)/sqrt(X2) with X1 ~ Ga(3/2,1) and X2 ~ Ga(a,1),
//     a = kappa - 1/2.  As a -> 0+ the denominator variate underflows and R is lost.
//     log R = (log X1 - log X2)/2 is finite whenever the mathematical draw exists, so it
//     is the arbiter that decides whether a lost draw was recoverable.
//
//  2. The final velocity is R times an order-unity vector.  Writing V = R g with
//     g = Q(ub) diag(sqrt(kappa) theta) n lets representability be decided from
//     log|V_j| = log R + log|g_j| before anything is exponentiated.  The released loader
//     instead forms R g_local and then rotates, so it can overflow on a draw whose final
//     components are perfectly representable.
//
//  3. "Failure" is three different things - avoidable arithmetic loss, honest output
//     overflow, and legal cap rejection - and they must never be pooled.

#ifndef EXP6_COMMON_H
#define EXP6_COMMON_H

#include "../../../cpp/bi_kappa_distribution.H"

#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <limits>
#include <random>
#include <string>
#include <vector>

namespace exp6
{

// ---------------------------------------------------------------------------
// Schema versions.  Bump when a record layout changes; analyze.py refuses to read
// a file whose header does not match what it was written to understand.
// ---------------------------------------------------------------------------
static const uint32_t kPilotSchemaVersion= 1;
static const uint32_t kCondSchemaVersion= 1;
static const uint32_t kLoaderSchemaVersion= 1;
static const uint32_t kAuditSchemaVersion= 1;

// ---------------------------------------------------------------------------
// Environment identification
// ---------------------------------------------------------------------------
inline const char *compilerId()
{
#if defined(__clang__)
    return "clang " __clang_version__;
#elif defined(__GNUC__)
    return "gcc " __VERSION__;
#else
    return "unknown";
#endif
}

inline const char *stdlibId()
{
#if defined(_LIBCPP_VERSION)
    return "libc++";
#elif defined(__GLIBCXX__)
    return "libstdc++";
#else
    return "unknown";
#endif
}

inline const char *archId()
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
template <> inline const char *precisionId<float>() { return "float"; }
template <> inline const char *precisionId<double>() { return "double"; }

// ---------------------------------------------------------------------------
// Exact floating-point limits, in the log domain.
//
// A real x rounds to infinity in type T exactly when |x| >= (2 - 2^-p) 2^emax, which is
// half an ulp above the largest finite value.  That threshold, not max(), is what decides
// honest overflow, and log of the two differs only in the 16th digit - which is precisely
// why draws near it are sent to the arbitrary-precision oracle instead of being
// classified here.
// ---------------------------------------------------------------------------
template <typename T> inline double logOverflowThreshold()
{
    const int p= std::numeric_limits<T>::digits;         // 24 (float), 53 (double)
    const int emax= std::numeric_limits<T>::max_exponent - 1; // 127, 1023
    const double kLn2= 0.69314718055994530942;
    return static_cast<double>(emax) * kLn2 + std::log(2.0 - std::ldexp(1.0, -p));
}

template <typename T> inline double logMaxFinite()
{
    return std::log(static_cast<double>(std::numeric_limits<T>::max()));
}

/// Draws whose intended log-component lands inside this band of the overflow threshold
/// are never classified locally; they are written to the audit stream for the oracle.
static const double kOracleBandLogUnits= 2.0;

// ---------------------------------------------------------------------------
// Stratified audit sampling.
//
// The independent oracle costs about 1.1 ms per record at 100 decimal digits.  At the
// frozen sample sizes the class "the released path failed" contains roughly 8e6 draws,
// which would be several CPU-hours of adjudication spent almost entirely on one trivially
// classified mode: the denominator materialized as exactly zero, thousands of log units
// away from any representability threshold.
//
// So the audit is stratified rather than exhaustive, with the strata and their rates fixed
// here in the committed source.  Every decision-relevant draw is still audited in full:
// anything within kOracleBandLogUnits of a type limit, where the classification could
// genuinely go either way.  The homogeneous classes are sampled deterministically, and
// both the audited and the total count of every stratum is reported, so the coverage is
// visible rather than implied.  A single disagreement in any stratum fails the gate.
// ---------------------------------------------------------------------------
enum AuditStratum
{
    kAuditNone= 0,
    kAuditNearLimit= 1,  ///< within kOracleBandLogUnits of the overflow threshold
    kAuditDisagree= 2,   ///< the formations reached different terminal categories
    kAuditSubnormal= 3,  ///< the denominator Gamma landed in the subnormal range
    kAuditBulkFail= 4,   ///< every formation failed and the draw is far from any threshold
    kAuditSample= 5,     ///< the unbiased deterministic sample of all attempts
    kNumAuditStrata= 6
};

inline const char *auditStratumName(int s)
{
    switch (s)
    {
    case kAuditNearLimit: return "near_limit";
    case kAuditDisagree: return "method_disagreement";
    case kAuditSubnormal: return "subnormal_denominator";
    case kAuditBulkFail: return "bulk_failure";
    case kAuditSample: return "uniform_sample";
    default: return "none";
    }
}

/// One in this many members of each stratum is audited.  1 means exhaustive.
inline long long auditStride(int stratum)
{
    switch (stratum)
    {
    case kAuditNearLimit: return 1;
    case kAuditDisagree: return 100;
    case kAuditSubnormal: return 10;
    case kAuditBulkFail: return 10000;
    case kAuditSample: return 10000;
    default: return 0;
    }
}

/// Retained name for the unbiased sample stride.
static const long long kOracleStride= 10000;

// ---------------------------------------------------------------------------
// Terminal categories.  Exactly one per method per attempt.
// ---------------------------------------------------------------------------
enum TerminalCategory
{
    kCatFinite= 0,        ///< a finite vector (or finite radius, in scalar phases) was produced
    kCatDenomZero= 1,     ///< denominator Gamma materialized as exactly zero
    kCatQuotientLoss= 2,  ///< the intermediate quotient X1/X2 overflowed; sqrt was representable
    kCatSplitLoss= 3,     ///< the split form still lost a representable draw
    kCatLogPrimFail= 4,   ///< the log-domain primitive produced a non-finite log
    kCatHonestOverflow= 5,///< at least one intended final component is not representable
    kCatCapReject= 6,     ///< intended sample lies outside the declared component box
    kCatExhausted= 7,     ///< capped mode hit the attempt limit
    kNumCategories= 8
};

inline const char *categoryName(int c)
{
    switch (c)
    {
    case kCatFinite: return "finite";
    case kCatDenomZero: return "denominator_zero";
    case kCatQuotientLoss: return "quotient_first_loss";
    case kCatSplitLoss: return "split_form_loss";
    case kCatLogPrimFail: return "log_primitive_failure";
    case kCatHonestOverflow: return "honest_overflow";
    case kCatCapReject: return "cap_reject";
    case kCatExhausted: return "cap_exhausted";
    default: return "unknown";
    }
}

/// Which categories count as avoidable loss in the accounting identity
/// N_attempts = N_finite + N_avoidable + N_honest (uncapped runs).
inline bool isAvoidable(int c)
{
    return c == kCatDenomZero || c == kCatQuotientLoss || c == kCatSplitLoss ||
           c == kCatLogPrimFail;
}

// ---------------------------------------------------------------------------
// Overlapping primitive flags.  These describe the draw, not a method's verdict, so
// several can be set at once.
// ---------------------------------------------------------------------------
enum PrimitiveFlag
{
    kFlagX2Zero= 1u << 0,
    kFlagX2Subnormal= 1u << 1,
    kFlagQfNonFinite= 1u << 2,
    kFlagSplitNonFinite= 1u << 3,
    kFlagLogNonFinite= 1u << 4,
    kFlagHonestOverflowRef= 1u << 5, ///< intended final vector not representable (reference)
    kFlagCapAccept= 1u << 6,
    kFlagAudited= 1u << 7,           ///< written to the oracle audit stream
    kFlagNearLimit= 1u << 8,         ///< within kOracleBandLogUnits of the overflow threshold
    kFlagUEndpointRedraw= 1u << 9,   ///< the open-interval uniform guard fired
    /// The *local* field-aligned vector overflowed while the rotated global components
    /// are all representable.  The released loader forms the local vector and then
    /// rotates, so it loses these draws; an orthogonal rotation can spread one
    /// over-large local component across three representable global ones.
    kFlagRotationRecoverable= 1u << 10,
    kFlagCapExhausted= 1u << 11
};

// ---------------------------------------------------------------------------
// Binary record layouts.  Every bulk file starts with this header.
// ---------------------------------------------------------------------------
struct FileHeader
{
    char magic[8];       // "EXP6REC\0"
    uint32_t schema;     // record schema version
    uint32_t record_kind;// 1 pilot, 2 conditioning, 3 loader, 4 audit
    uint32_t record_size;// sizeof(record)
    uint32_t reserved;
    double kappa;
    double theta_perp;
    double theta_par;
    double ub[3];
    double cap;          // +inf when uncapped
    uint32_t seed;
    uint32_t precision_is_float;
    int64_t n_records;   // filled in on close
};

struct PilotRecord
{
    double log_r;   // (log X1 - log X2) / 2, in the working precision, widened for storage
    double log_w;   // log X2 - logaddexp(log X1, log X2); W ~ Beta(a, 3/2)
    uint32_t retries;
    uint32_t flags;
};

struct CondRecord
{
    double log_r_ref;             // intended log R, recomputed in double from the primitives
    double log_w_ref;             // intended log W
    double log_speed_ref;         // intended log |V|
    double max_log_component_ref; // intended max_j log |V_j|
    uint32_t flags;
    uint8_t cat_qf;
    uint8_t cat_split;
    uint8_t cat_log;
    uint8_t reserved;
};

struct LoaderRecord
{
    double v[3];      // returned components in the global frame (may be non-finite)
    double log_r_ref; // intended log R for this attempt
    uint32_t status;  // terminal category in the low byte, flags above it
    uint32_t attempts;// attempts consumed to produce this record (capped mode)
};

/// One row handed to the arbitrary-precision oracle.  It carries the *primitives*, not
/// any derived quantity, so the oracle re-derives everything independently.
struct AuditRecord
{
    double x1;       // Ga(3/2,1) variate, exactly as produced in the working precision
    double y;        // Ga(a+1,1) variate of the boosting identity
    double u;        // U(0,1) variate of the boosting identity
    double cos_theta;
    double phi;
    double log_x2_working; // what the working-precision log path computed
    double kappa;          // carried per record: one audit stream spans many configurations
    uint32_t flags;
    uint32_t precision_is_float;
    uint8_t cat_qf;
    uint8_t cat_split;
    uint8_t cat_log;
    uint8_t reserved;
    uint32_t pad;
};

// ---------------------------------------------------------------------------
// Buffered binary writer with an atomic finish: the file is built under a .part name,
// its header is rewritten with the final record count, and only then is it renamed.
// A truncated job therefore never leaves a file that looks complete.
// ---------------------------------------------------------------------------
template <typename Record> class RecordWriter
{
  public:
    RecordWriter(const std::string &path, uint32_t schema, uint32_t kind, double kappa,
                 double theta_perp, double theta_par, const double ub[3], double cap,
                 uint32_t seed, bool is_float)
        : m_path(path), m_part(path + ".part"), m_f(0), m_n(0)
    {
        m_f= std::fopen(m_part.c_str(), "wb");
        if (!m_f)
        {
            std::fprintf(stderr, "exp6: cannot open %s\n", m_part.c_str());
            std::exit(1);
        }
        std::memset(&m_hdr, 0, sizeof(m_hdr));
        std::memcpy(m_hdr.magic, "EXP6REC", 8);
        m_hdr.schema= schema;
        m_hdr.record_kind= kind;
        m_hdr.record_size= static_cast<uint32_t>(sizeof(Record));
        m_hdr.kappa= kappa;
        m_hdr.theta_perp= theta_perp;
        m_hdr.theta_par= theta_par;
        m_hdr.ub[0]= ub[0];
        m_hdr.ub[1]= ub[1];
        m_hdr.ub[2]= ub[2];
        m_hdr.cap= cap;
        m_hdr.seed= seed;
        m_hdr.precision_is_float= is_float ? 1u : 0u;
        m_hdr.n_records= 0;
        std::fwrite(&m_hdr, sizeof(m_hdr), 1, m_f);
        m_buf.reserve(65536);
    }

    void push(const Record &r)
    {
        m_buf.push_back(r);
        ++m_n;
        if (m_buf.size() >= 65536)
            flush();
    }

    long long count() const { return m_n; }

    /// Returns the final path.  After this the object must not be used again.
    std::string finish()
    {
        flush();
        m_hdr.n_records= m_n;
        std::fseek(m_f, 0, SEEK_SET);
        std::fwrite(&m_hdr, sizeof(m_hdr), 1, m_f);
        std::fclose(m_f);
        m_f= 0;
        if (std::rename(m_part.c_str(), m_path.c_str()) != 0)
        {
            std::fprintf(stderr, "exp6: cannot rename %s -> %s\n", m_part.c_str(),
                         m_path.c_str());
            std::exit(1);
        }
        return m_path;
    }

  private:
    void flush()
    {
        if (!m_buf.empty())
        {
            std::fwrite(&m_buf[0], sizeof(Record), m_buf.size(), m_f);
            m_buf.clear();
        }
    }

    std::string m_path, m_part;
    FILE *m_f;
    long long m_n;
    FileHeader m_hdr;
    std::vector<Record> m_buf;
};

// ---------------------------------------------------------------------------
// The denominator shape, formed exactly as the released loader forms it.
//
// bi_kappa_distribution computes `kappa - kHalf` in RealType.  Forming kappa - 1/2 in
// double and narrowing afterwards gives a different float: at kappa = 0.55 the two differ
// in the last bits of the shape, which changes every subsequent Gamma variate.  A replica
// that gets this wrong diverges from the released class immediately and silently.
// ---------------------------------------------------------------------------
template <typename T> inline T shapeFromKappa(double kappa)
{
    return static_cast<T>(kappa) - T(0.5);
}

/// The same shape widened for the double-precision reference calculation.  The intended
/// mathematical draw is the one the working-precision shape defines, not the one an
/// infinitely precise kappa - 1/2 would define.
template <typename T> inline double shapeReference(double kappa)
{
    return static_cast<double>(shapeFromKappa<T>(kappa));
}

// ---------------------------------------------------------------------------
// Numerics used in the log domain
// ---------------------------------------------------------------------------
template <typename T> inline T logAddExp(T a, T b)
{
    if (a == b)
        return a + std::log(T(2));
    const T m= a > b ? a : b;
    const T d= a > b ? b - a : a - b;
    return m + std::log1p(std::exp(d));
}

inline double logAddExpD(double a, double b) { return logAddExp<double>(a, b); }

// ---------------------------------------------------------------------------
// Geometry.  g is the order-unity vector that multiplies R to give the final velocity,
// already rotated into the global frame.  Building it before exponentiating is the whole
// reason the log path can decide representability of the *returned* components rather
// than of the local ones.
// ---------------------------------------------------------------------------
template <typename T> struct Geometry
{
    T kappa;
    T sqrt_kappa;
    T theta_perp;
    T theta_par;
    std::array<T, 3> ub;
    bool identity_rotation;
    bi_kappa_distribution<T> rotator; // used only for rotate_from_fieldAligned_frame

    Geometry(T k, T tperp, T tpar, const std::array<T, 3> &b)
        : kappa(k), sqrt_kappa(std::sqrt(k)), theta_perp(tperp), theta_par(tpar), ub(b),
          identity_rotation(b[0] == T(0) && b[1] == T(0) && b[2] == T(1)),
          rotator(k, tperp, tpar, b, bi_kappa_distribution<T>::no_cap())
    {
    }

    /// Unit direction from the two released-loader uniforms, in the field-aligned frame.
    static std::array<T, 3> direction(T cos_theta, T phi)
    {
        const T sin_theta= std::sqrt(std::max(T(0), T(1) - cos_theta * cos_theta));
        std::array<T, 3> n;
        n[0]= sin_theta * std::cos(phi);
        n[1]= sin_theta * std::sin(phi);
        n[2]= cos_theta;
        return n;
    }

    /// g = Q(ub) diag(sqrt(kappa) theta_perp, sqrt(kappa) theta_perp, sqrt(kappa) theta_par) n
    std::array<T, 3> gVector(const std::array<T, 3> &n) const
    {
        std::array<T, 3> local;
        local[0]= sqrt_kappa * theta_perp * n[0];
        local[1]= sqrt_kappa * theta_perp * n[1];
        local[2]= sqrt_kappa * theta_par * n[2];
        if (identity_rotation)
            return local;
        return rotator.rotate_from_fieldAligned_frame(local, ub);
    }
};

// ---------------------------------------------------------------------------
// Small JSON emission helpers.  One object per line, machine-readable, no dependency.
// ---------------------------------------------------------------------------
inline void emitEnvFields(FILE *f, const char *precision)
{
    std::fprintf(f,
                 "\"compiler\":\"%s\",\"stdlib\":\"%s\",\"arch\":\"%s\",\"precision\":\"%s\","
                 "\"rng\":\"mt19937\"",
                 compilerId(), stdlibId(), archId(), precision);
}

/// A cap of +infinity means "uncapped".  JSON has no infinity literal, so it is emitted as
/// null: a bare `inf` is not parseable and quoting it would make the column a string in
/// some rows and a number in others.
inline void emitCap(FILE *f, double cap)
{
    if (cap == std::numeric_limits<double>::infinity())
        std::fprintf(f, "\"cap\":null");
    else
        std::fprintf(f, "\"cap\":%.10g", cap);
}

/// %.17g round-trips a double exactly; anything less silently loses the last digits of a
/// rate that the analysis then reports to more places than it actually has.
inline void emitDouble(FILE *f, const char *key, double v)
{
    if (v == std::numeric_limits<double>::infinity())
        std::fprintf(f, "\"%s\":\"inf\"", key);
    else if (v == -std::numeric_limits<double>::infinity())
        std::fprintf(f, "\"%s\":\"-inf\"", key);
    else if (v != v)
        std::fprintf(f, "\"%s\":\"nan\"", key);
    else
        std::fprintf(f, "\"%s\":%.17g", key, v);
}

// ---------------------------------------------------------------------------
// Frozen constants, spelled out in the committed source so that a run is determined by
// the command line recorded in README.md and nothing else.
// ---------------------------------------------------------------------------
inline std::vector<double> kappaLadder()
{
    std::vector<double> k;
    k.push_back(0.5001); k.push_back(0.501); k.push_back(0.505); k.push_back(0.51);
    k.push_back(0.55);   k.push_back(0.60);  k.push_back(0.75);  k.push_back(1.0);
    k.push_back(1.5);    k.push_back(2.0);   k.push_back(5.0);
    return k;
}

inline std::vector<unsigned> productionSeeds()
{
    std::vector<unsigned> s;
    s.push_back(4001); s.push_back(4002); s.push_back(4003); s.push_back(4004);
    s.push_back(4005);
    return s;
}

inline std::vector<unsigned> pilotSeeds()
{
    std::vector<unsigned> s;
    s.push_back(4001); s.push_back(4002); s.push_back(4003);
    return s;
}

} // namespace exp6

#endif // EXP6_COMMON_H
