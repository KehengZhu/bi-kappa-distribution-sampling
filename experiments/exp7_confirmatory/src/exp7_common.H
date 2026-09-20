// Experiment 7 - shared definitions for the confirmatory probe.
//
// Strict C++11.  Nothing in this header, or in anything it includes, may depend on a
// later standard: the released loader is a header-only C++11 reference sampler, so the
// sampling code has to keep compiling as C++11.  The arbitrary-precision oracle lives in a
// separate translation unit (src/exp7_oracle.cpp, C++14) precisely so that it cannot leak a
// dependency in here, and so that it re-derives every classification from the recorded
// primitives rather than sharing code with the thing it audits.
//
// Four ideas carry the whole experiment and are worth stating once:
//
//  1. The radius is R = sqrt(X1)/sqrt(X2) with X1 ~ Ga(3/2,1) and X2 ~ Ga(a,1),
//     a = kappa - 1/2.  As a -> 0+ the denominator variate underflows and R is lost.
//     log R = (log X1 - log X2)/2 is finite whenever the mathematical draw exists, so it
//     is the arbiter that decides whether a lost draw was recoverable.
//
//  2. The final velocity is R times an order-unity vector.  Writing V = R g with
//     g = Q(ub) diag(sqrt(kappa) theta) n lets each component be materialized on its own and
//     its representability decided from the component itself -- the number that will be
//     returned -- rather than from a logarithm compared against a rounded log(max()).
//     LEGACY instead forms R n, scales, and then rotates, so it can overflow on a draw whose
//     final components are perfectly representable.
//
//  3. "Failure" is three different things - avoidable arithmetic loss, honest output
//     overflow, and legal cap rejection - and they must never be pooled.  The accounting
//     identity N = finite + avoidable + honest is checked per seed and configuration.
//
//  4. CANDIDATE's stream is a function of the engine alone, so two builds against
//     different standard libraries must agree *byte for byte*.  That is tested with
//     digests over the returned bit patterns, not with interval overlap.

#ifndef EXP7_COMMON_H
#define EXP7_COMMON_H

#include "../../../cpp/bi_kappa_distribution.H"
#include "legacy/bi_kappa_distribution_v1.H"

#include "exp7_digest.H"
#include "exp7_protocol.H"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <limits>
#include <random>
#include <string>
#include <vector>

namespace exp7
{

// ---------------------------------------------------------------------------
// Schema versions.  Bump when a record layout changes; the analysis refuses to read a file
// whose header does not match what it was written to understand.
// ---------------------------------------------------------------------------
static const uint32_t kCondSchemaVersion= 1;
static const uint32_t kLoaderSchemaVersion= 1;
static const uint32_t kAuditSchemaVersion= 1;
static const uint32_t kTailSchemaVersion= 1;

/// Record kinds carried in the binary file header.
static const uint32_t kKindCond= 2;
static const uint32_t kKindLoader= 3;
static const uint32_t kKindAudit= 4;
static const uint32_t kKindTail= 5;

/// The manifest schema, so that a reader can tell an old manifest from a new one.
static const char *const kManifestSchema= "exp7-manifest-1";

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

/// The build's compiler flags, injected by the GNUmakefile.  Experiment 6 rebuilt its
/// manifest at analysis time and recorded whichever flags the analysis happened to think
/// the build used; 165 of its rows carried the wrong ones.  Here the string is compiled
/// into the binary that wrote the file.
#ifndef EXP7_CXXFLAGS
#define EXP7_CXXFLAGS "unrecorded"
#endif
inline const char *cxxFlagsId() { return EXP7_CXXFLAGS; }

template <typename T> const char *precisionId();
template <> inline const char *precisionId<float>() { return "float"; }
template <> inline const char *precisionId<double>() { return "double"; }

/// UTC timestamp, seconds resolution, ISO-8601.
inline std::string nowUtcIso()
{
    std::time_t t= std::time(0);
    std::tm g;
#if defined(_WIN32)
    gmtime_s(&g, &t);
#else
    gmtime_r(&t, &g);
#endif
    char buf[32];
    std::strftime(buf, sizeof(buf), "%Y-%m-%dT%H:%M:%SZ", &g);
    return std::string(buf);
}

// ---------------------------------------------------------------------------
// Exact floating-point limits, in the log domain.
//
// A real x rounds to infinity in type T exactly when |x| >= (2 - 2^-p) 2^emax, which is
// half an ulp above the largest finite value.  That threshold, not max(), is what decides
// whether the *intended* draw was representable, and log of the two differs only in the
// 16th digit - which is precisely why draws near it are sent to the arbitrary-precision
// oracle instead of being classified here.
//
// The released 2.0.0 sampler makes its own decision against log(max()) rather than against
// this threshold; that is the predicate a user actually gets, so the formations below
// reproduce it exactly and the reference keeps the exact threshold.  Draws between the two
// are inside the oracle band by construction.
// ---------------------------------------------------------------------------
template <typename T> inline double logOverflowThreshold()
{
    const int p= std::numeric_limits<T>::digits;              // 24 (float), 53 (double)
    const int emax= std::numeric_limits<T>::max_exponent - 1; // 127, 1023
    const double kLn2= 0.69314718055994530942;
    return static_cast<double>(emax) * kLn2 + std::log(2.0 - std::ldexp(1.0, -p));
}

/// The logarithm of the largest finite value, TAKEN IN T, which is the value the released
/// header's capped screen actually compares against.  Computing it in `double` and narrowing
/// gives a different number for `float` -- 88.72283905206835 against 88.72283935546875 -- and
/// reporting the one the code does not use is how a record of a finite-precision experiment
/// comes to disagree with the experiment.  The released loader no longer decides
/// representability from this quantity at all (see bikappa_detail::materializeComponents);
/// it remains a screening threshold and an environment record.
template <typename T> inline double logMaxFinite()
{
    return static_cast<double>(std::log(std::numeric_limits<T>::max()));
}

/// Draws whose intended log-component lands inside this band of a type limit are the
/// decision-relevant ones: they are written to the audit stream in full and never settled by
/// the working-precision reference alone.  The width is the protocol's, not this file's --
/// amendment 1.2.0 states audit coverage as a margin so that the claim rests on an
/// inequality rather than on a sampling rate.  Twenty natural-log units is about four
/// thousand times the worst reference-versus-oracle disagreement ever measured here.
inline double oracleBandLogUnits() { return protocol::kAuditMarginLogUnits; }

// ---------------------------------------------------------------------------
// Stratified audit sampling, by margin rather than by rate.
//
// The independent oracle costs of order a millisecond per record at 100 decimal digits, so
// which draws it sees is a real decision and amendment 1.2.0 takes it out of the
// implementation's hands.  An attempt is *decision-relevant* -- adjudicated in full -- when
// its intended log-component lies within `audit.margin_log_units` of a type limit, when the
// formations reach different terminal categories, when CANDIDATE returned a finite but
// inaccurate radius, when CANDIDATE suffered avoidable loss, or when a subnormal or zero
// denominator coincided with a representable target.  Everything else is *unambiguous*: the
// intended value lies so far outside the type's range that no arithmetic could have returned
// it.  Those are sampled, and the margin is asserted in working precision for every one of
// them, so the claim rests on an inequality rather than on a sampling rate.
//
// Experiment 6 sampled its bulk failures at 1 in 4957 while its plan said "every public-path
// failure"; with zero disagreements in 7039 of them the one-sided 95 % bound still admitted
// about 14 900 misclassified failures among 34.9 million.  Requiring every public-path
// failure instead is both infeasible and beside the point, because almost all of them are
// thousands of log units outside the range.  The margin is what closes the gap.  BOTH the
// audited and the total count of every stratum is reported, so the coverage is visible
// rather than implied, and a single disagreement in any stratum fails gate G2.
// ---------------------------------------------------------------------------
enum AuditStratum
{
    kAuditNone= 0,
    kAuditWithinMargin= 1,   ///< intended log-component within the margin of a type limit
    kAuditDisagree= 2,       ///< the formations reached different terminal categories
    kAuditFiniteButWrong= 3, ///< CANDIDATE returned a finite radius that is wrong
    kAuditAvoidableCandidate= 4, ///< CANDIDATE suffered avoidable loss; this must stay empty
    kAuditSubnormalRepresentable= 5, ///< subnormal or zero denominator, representable target
    kAuditUnambiguous= 6,    ///< a failure more than the margin beyond a type limit
    kAuditSample= 7,         ///< everything else, sampled without reference to its content
    kNumAuditStrata= 8
};

/// The stratum names are the protocol's own, so a reader can line the reported coverage up
/// against `audit.strata_rates` in config/protocol.json without a translation table.
inline const char *auditStratumName(int s)
{
    switch (s)
    {
    case kAuditWithinMargin: return protocol::kAuditStratumNames[0];
    case kAuditDisagree: return protocol::kAuditStratumNames[1];
    case kAuditFiniteButWrong: return protocol::kAuditStratumNames[2];
    case kAuditAvoidableCandidate: return protocol::kAuditStratumNames[3];
    case kAuditSubnormalRepresentable: return protocol::kAuditStratumNames[4];
    case kAuditUnambiguous: return protocol::kAuditStratumNames[5];
    case kAuditSample: return protocol::kAuditStratumNames[6];
    default: return "none";
    }
}

inline double auditStratumRate(int s)
{
    switch (s)
    {
    case kAuditWithinMargin: return protocol::kAuditStratumRates[0];
    case kAuditDisagree: return protocol::kAuditStratumRates[1];
    case kAuditFiniteButWrong: return protocol::kAuditStratumRates[2];
    case kAuditAvoidableCandidate: return protocol::kAuditStratumRates[3];
    case kAuditSubnormalRepresentable: return protocol::kAuditStratumRates[4];
    case kAuditUnambiguous: return protocol::kAuditStratumRates[5];
    case kAuditSample: return protocol::kAuditStratumRates[6];
    default: return 0.0;
    }
}

/// One in this many members of each stratum is audited; 1 means exhaustive.  The rates come
/// from config/protocol.json rather than from this file.
inline long long auditStride(int stratum)
{
    const double rate= auditStratumRate(stratum);
    if (!(rate > 0.0))
        return 0;
    if (rate >= 1.0)
        return 1;
    const double inv= 1.0 / rate;
    long long stride= static_cast<long long>(inv + 0.5);
    if (stride < 1)
        stride= 1;
    return stride;
}

// ---------------------------------------------------------------------------
// Terminal categories.  Exactly one per method per attempt.
// ---------------------------------------------------------------------------
enum TerminalCategory
{
    kCatFinite= 0,         ///< a finite vector was produced and its radius is accurate
    kCatDenomZero= 1,      ///< denominator Gamma materialized as exactly zero
    kCatQuotientLoss= 2,   ///< the intermediate quotient X1/X2 overflowed (QF diagnostic only)
    kCatLegacyLoss= 3,     ///< LEGACY's split form still lost a representable draw
    kCatLogPrimFail= 4,    ///< the log-domain primitive produced a non-finite log
    kCatHonestOverflow= 5, ///< intended vector not representable; a non-finite vector returned
    kCatCapReject= 6,      ///< intended sample lies outside the declared component box
    kCatExhausted= 7,      ///< capped mode hit the attempt limit
    /// Amendment 1.1.0.  The intended draw was representable and a finite vector came back,
    /// but its radius has lost more than half the significand of its type against the
    /// reference.  It is a loss, not a success: the value is a number, and the number is
    /// wrong.  Experiment 6 had no such category and scored these as successes.
    kCatFiniteButWrong= 8,
    /// The intended draw was NOT representable and the method nevertheless returned a finite
    /// vector.  Also a loss, but an *honest* one: no formation can rescue a draw that has no
    /// representation, so this is not evidence that the arithmetic threw anything away.
    kCatOverflowReturnedFinite= 9,
    kNumCategories= 10
};

inline const char *categoryName(int c)
{
    switch (c)
    {
    case kCatFinite: return "finite";
    case kCatDenomZero: return "denominator_zero";
    case kCatQuotientLoss: return "quotient_first_loss";
    case kCatLegacyLoss: return "legacy_form_loss";
    case kCatLogPrimFail: return "log_primitive_failure";
    case kCatHonestOverflow: return "honest_overflow";
    case kCatCapReject: return "cap_reject";
    case kCatExhausted: return "cap_exhausted";
    case kCatFiniteButWrong: return "finite_but_wrong";
    case kCatOverflowReturnedFinite: return "overflow_returned_finite";
    default: return "unknown";
    }
}

/// Which categories count as avoidable loss in the accounting identity
/// N_attempts = N_finite + N_avoidable + N_honest (uncapped runs).  A finite-but-wrong draw
/// belongs here because the intended draw *was* representable: a formation that carries the
/// radius in the log domain returns it accurately.
inline bool isAvoidable(int c)
{
    return c == kCatDenomZero || c == kCatQuotientLoss || c == kCatLegacyLoss ||
           c == kCatLogPrimFail || c == kCatFiniteButWrong;
}

/// The honest half of the same identity: the intended draw has no representation in the
/// type, whether the method admitted that or returned a number anyway.
inline bool isHonest(int c)
{
    return c == kCatHonestOverflow || c == kCatOverflowReturnedFinite;
}

/// The relative-radius error above which a finite return is FINITE_BUT_WRONG, per type.
template <typename T> inline double maxRelErrorFor();
template <> inline double maxRelErrorFor<double>() { return protocol::kMaxRelErrorDouble; }
template <> inline double maxRelErrorFor<float>() { return protocol::kMaxRelErrorFloat; }

// ---------------------------------------------------------------------------
// Overlapping primitive flags.  These describe the draw, not a method's verdict, so
// several can be set at once.
// ---------------------------------------------------------------------------
enum PrimitiveFlag
{
    kFlagX2Zero= 1u << 0,
    kFlagX2Subnormal= 1u << 1,
    kFlagQfNonFinite= 1u << 2,
    kFlagLegacyNonFinite= 1u << 3,
    kFlagCandidateNonFinite= 1u << 4,
    kFlagHonestOverflowRef= 1u << 5, ///< intended final vector not representable (reference)
    kFlagCapAccept= 1u << 6,
    kFlagAudited= 1u << 7,   ///< written to the oracle audit stream
    kFlagNearLimit= 1u << 8, ///< within the protocol margin of a type limit
    /// The *local* field-aligned vector overflowed while the rotated global components are
    /// all representable.  LEGACY forms the local vector and then rotates, so it loses
    /// these draws; an orthogonal rotation can spread one over-large local component across
    /// three representable global ones.
    kFlagRotationRecoverable= 1u << 9,
    kFlagCapExhausted= 1u << 10,
    kFlagLegacySuccess= 1u << 11,   ///< LEGACY returned a finite vector for this attempt
    kFlagCandidateSuccess= 1u << 12 ///< CANDIDATE returned a finite vector for this attempt
};

// ---------------------------------------------------------------------------
// Binary record layouts.  Every bulk file starts with this header.
// ---------------------------------------------------------------------------
struct FileHeader
{
    char magic[8];        // "EXP7REC\0"
    uint32_t schema;      // record schema version
    uint32_t record_kind; // 2 conditioning, 3 loader, 4 audit
    uint32_t record_size; // sizeof(record)
    uint32_t reserved;
    double kappa;
    double theta_perp;
    double theta_par;
    double ub[3];
    double cap; // +inf when uncapped
    uint32_t seed;
    uint32_t precision_is_float;
    int64_t n_records; // filled in on close
};

/// P3 conditioning.  One record per attempt.  log_w_ref is not stored: it is an exact
/// function of log_r_ref, log W = -logaddexp(0, 2 log R), and storing it would add eight
/// bytes per attempt to the largest phase in the run for no information.
struct CondRecord
{
    double log_r_ref;             // intended log R, recomputed in double from the primitives
    double log_speed_ref;         // intended log |V|
    double max_log_component_ref; // intended max_j log |V_j|
    uint32_t flags;               // PrimitiveFlag bits, incl. per-method success
    uint8_t cat_legacy;
    uint8_t cat_candidate;
    uint8_t cat_qf;
    uint8_t reserved;
};

/// P4 loader.  One record per *returned* sample.  `attempts` is how many attempts that
/// return consumed, which is 1 in uncapped mode and the cap's rejection run length
/// otherwise, so the attempt stream is recoverable from the record stream plus the
/// per-configuration rejection counters.
struct LoaderRecord
{
    double v[3];       // returned components in the global frame (may be non-finite)
    double log_r_ref;  // intended log R for this attempt
    uint32_t status;   // terminal category
    uint32_t attempts; // attempts consumed to produce this record
};

/// P1 upper tail.  One record per draw whose Z exceeded the *lowest* protocol threshold
/// z0 = -log(1e-2), which nests the other two, so F4's exceedance counts and its KS of the
/// excesses against Exp(1) are both exact for every q0 in the protocol from this file plus
/// the resolved count in the summary row.  Draws whose Z is not finite are not written --
/// they exceed every threshold by construction and are carried as a count instead, which
/// keeps the file at its expected 10^4 records rather than at 10^6.
struct TailRecord
{
    double z;
};

/// One row handed to the arbitrary-precision oracle.  It carries the *primitives*, not any
/// derived quantity, so the oracle re-derives everything independently.
struct AuditRecord
{
    double x1;             // Ga(3/2,1) variate, exactly as produced in the working precision
    double y;              // Ga(a+1,1) variate of the boosting identity
    double u;              // U(0,1) variate of the boosting identity
    // The paired layer's two declared direction primitives.  They are the layer's own
    // declaration, not a record of what either released sampler drew: 1.0.0 draws its
    // angles through std::uniform_real_distribution, and 2.0.0 draws no angle at all.
    // Both are exact functions of the engine, so this stream is identical across standard
    // libraries by construction rather than by result.
    double cos_theta;
    double phi;
    double log_x2_working; // what the working-precision log path computed
    double kappa;          // carried per record: one audit stream spans many configurations
    uint32_t flags;
    uint32_t precision_is_float;
    uint8_t cat_legacy;
    uint8_t cat_candidate;
    uint8_t cat_qf;
    uint8_t reserved;
    uint32_t pad;
};

// ---------------------------------------------------------------------------
// The run-time manifest.
//
// Experiment 6 rebuilt raw/manifest.csv from whatever was on disk when the analysis ran, so
// `probe_sha256` was the hash of the binary present at analysis time rather than of the one
// that produced the file.  Here the probe writes one row per output file, from inside the
// process that wrote it, after the atomic rename and after re-opening the file and checking
// its header, record count and size against what was intended.  Smoke output goes to its own
// manifest under raw/smoke/ and is never appended to the production one.
// ---------------------------------------------------------------------------
struct RunContext
{
    std::string manifest_path;
    std::string run_utc;
    std::string probe_path;
    std::string probe_sha256;
    std::string tag;       ///< environment tag: libcxx, libstdcxx, ...
    std::string execution; ///< native | translated (Rosetta, emulation)
    bool smoke;

    RunContext() : smoke(false) {}
};

inline std::string csvField(const std::string &s)
{
    bool needs= false;
    for (size_t i= 0; i < s.size(); ++i)
        if (s[i] == ',' || s[i] == '"' || s[i] == '\n')
            needs= true;
    if (!needs)
        return s;
    std::string out= "\"";
    for (size_t i= 0; i < s.size(); ++i)
    {
        if (s[i] == '"')
            out+= '"';
        out+= s[i];
    }
    out+= '"';
    return out;
}

/// Append one manifest row.  Fails the process rather than recording an unverified file.
inline void manifestAppend(const RunContext &rc, const char *phase, const char *layer,
                           const char *method, const char *precision, const char *case_id,
                           double kappa, long long seed, const std::string &path,
                           const char *format, long long n_records, long long record_size,
                           long long schema, long long kind, long long bytes,
                           const std::string &sha256)
{
    const bool fresh= (std::fopen(rc.manifest_path.c_str(), "rb") == 0);
    FILE *f= std::fopen(rc.manifest_path.c_str(), "a");
    if (!f)
    {
        std::fprintf(stderr, "exp7: cannot append to %s\n", rc.manifest_path.c_str());
        std::exit(1);
    }
    if (fresh)
        std::fprintf(f, "manifest_schema,run_utc,phase,layer,method,precision,case,kappa,seed,"
                        "file,format,n_records,record_size,schema,record_kind,bytes,sha256,"
                        "probe,probe_sha256,cxxflags,compiler,stdlib,arch,env_tag,execution,"
                        "protocol_sha256,smoke\n");
    std::fprintf(f, "%s,%s,%s,%s,%s,%s,%s,%.17g,%lld,%s,%s,%lld,%lld,%lld,%lld,%lld,%s,%s,%s,"
                    "%s,%s,%s,%s,%s,%s,%s,%s\n",
                 kManifestSchema, rc.run_utc.c_str(), phase, layer, method, precision, case_id,
                 kappa, seed, csvField(path).c_str(), format, n_records, record_size, schema,
                 kind, bytes, sha256.c_str(), csvField(rc.probe_path).c_str(),
                 rc.probe_sha256.c_str(), csvField(std::string(cxxFlagsId())).c_str(),
                 csvField(std::string(compilerId())).c_str(), stdlibId(), archId(),
                 rc.tag.c_str(), rc.execution.c_str(), protocol::kProtocolSha256,
                 rc.smoke ? "true" : "false");
    std::fclose(f);
}

/// Re-open a finished binary record file and check that it says what it was written to say.
/// Returns its size in bytes; exits on any mismatch.
inline long long verifyBinaryContract(const std::string &path, uint32_t schema, uint32_t kind,
                                      uint32_t record_size, long long expected_records)
{
    FILE *f= std::fopen(path.c_str(), "rb");
    if (!f)
    {
        std::fprintf(stderr, "exp7: cannot reopen %s for verification\n", path.c_str());
        std::exit(1);
    }
    FileHeader h;
    const bool read_ok= (std::fread(&h, sizeof(h), 1, f) == 1);
    std::fseek(f, 0, SEEK_END);
    const long long bytes= static_cast<long long>(std::ftell(f));
    std::fclose(f);
    const long long want= static_cast<long long>(sizeof(FileHeader)) +
                          static_cast<long long>(record_size) * expected_records;
    if (!read_ok || std::memcmp(h.magic, "EXP7REC", 7) != 0 || h.schema != schema ||
        h.record_kind != kind || h.record_size != record_size ||
        h.n_records != expected_records || bytes != want)
    {
        std::fprintf(stderr,
                     "exp7: %s fails its own contract "
                     "(schema %u/%u kind %u/%u size %u/%u records %lld/%lld bytes %lld/%lld)\n",
                     path.c_str(), h.schema, schema, h.record_kind, kind, h.record_size,
                     record_size, static_cast<long long>(h.n_records), expected_records, bytes,
                     want);
        std::exit(1);
    }
    return bytes;
}

/// The same for a line-oriented text file: the line count must be what the writer counted.
inline long long verifyJsonlContract(const std::string &path, long long expected_lines)
{
    FILE *f= std::fopen(path.c_str(), "rb");
    if (!f)
    {
        std::fprintf(stderr, "exp7: cannot reopen %s for verification\n", path.c_str());
        std::exit(1);
    }
    long long lines= 0, bytes= 0;
    std::vector<char> buf(1 << 16);
    for (;;)
    {
        const size_t got= std::fread(&buf[0], 1, buf.size(), f);
        for (size_t i= 0; i < got; ++i)
            if (buf[i] == '\n')
                ++lines;
        bytes+= static_cast<long long>(got);
        if (got < buf.size())
            break;
    }
    std::fclose(f);
    if (lines != expected_lines)
    {
        std::fprintf(stderr, "exp7: %s has %lld lines, expected %lld\n", path.c_str(), lines,
                     expected_lines);
        std::exit(1);
    }
    return bytes;
}

// ---------------------------------------------------------------------------
// Buffered binary writer with an atomic finish: the file is built under a .part name, its
// header is rewritten with the final record count, and only then is it renamed.  A
// truncated job therefore never leaves a file that looks complete.  This is Experiment 6's
// design, unchanged; what is new is that finish() also verifies and registers the file.
// ---------------------------------------------------------------------------
template <typename Record> class RecordWriter
{
  public:
    RecordWriter(const std::string &path, uint32_t schema, uint32_t kind, double kappa,
                 double theta_perp, double theta_par, const double ub[3], double cap,
                 uint32_t seed, bool is_float)
        : m_path(path), m_part(path + ".part"), m_f(0), m_n(0), m_schema(schema), m_kind(kind)
    {
        m_f= std::fopen(m_part.c_str(), "wb");
        if (!m_f)
        {
            std::fprintf(stderr, "exp7: cannot open %s\n", m_part.c_str());
            std::exit(1);
        }
        std::memset(&m_hdr, 0, sizeof(m_hdr));
        std::memcpy(m_hdr.magic, "EXP7REC", 8);
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

    /// Closes, renames, verifies the file against its own contract and appends the manifest
    /// row.  Returns the final path.  After this the object must not be used again.
    std::string finish(const RunContext &rc, const char *phase, const char *layer,
                       const char *method, const char *precision, const char *case_id,
                       double kappa, long long seed)
    {
        const std::string out= finishRaw();
        const long long bytes= verifyBinaryContract(out, m_schema, m_kind,
                                                    static_cast<uint32_t>(sizeof(Record)), m_n);
        manifestAppend(rc, phase, layer, method, precision, case_id, kappa, seed, out, "bin",
                       m_n, static_cast<long long>(sizeof(Record)), m_schema, m_kind, bytes,
                       Sha256::file(out));
        return out;
    }

    /// Close and rename without registering.  Used only by the audit streams, which span
    /// many configurations and are registered once by their owning phase.
    std::string finishRaw()
    {
        flush();
        m_hdr.n_records= m_n;
        std::fseek(m_f, 0, SEEK_SET);
        std::fwrite(&m_hdr, sizeof(m_hdr), 1, m_f);
        std::fclose(m_f);
        m_f= 0;
        if (std::rename(m_part.c_str(), m_path.c_str()) != 0)
        {
            std::fprintf(stderr, "exp7: cannot rename %s -> %s\n", m_part.c_str(),
                         m_path.c_str());
            std::exit(1);
        }
        return m_path;
    }

    uint32_t schema() const { return m_schema; }
    uint32_t kind() const { return m_kind; }

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
    uint32_t m_schema, m_kind;
    FileHeader m_hdr;
    std::vector<Record> m_buf;
};

// ---------------------------------------------------------------------------
// The denominator shape, formed exactly as the loaders form it.
//
// Both headers compute `kappa - 1/2` in RealType.  Forming kappa - 1/2 in double and
// narrowing afterwards gives a different float: at kappa = 0.55 the two differ in the last
// bits of the shape, which changes every subsequent Gamma variate.  A replica that gets
// this wrong diverges from the released class immediately and silently.
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

/// Pi exactly as release 1.0.0 spells it, so that the uniform on phi has the identical
/// bound.  Release 2.0.0 does not spell pi at all: its direction comes from a rejection
/// method on the unit disc and needs no angle.  This constant therefore serves the LEGACY
/// replica and the paired layer's declared direction primitive, and nothing else.
template <typename T> inline T piValue()
{
    return T(3.14159265358979323846264338327950288L);
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

/// log W for W = X2/(X1+X2) = expit(-2 log R).  Never forms R or R^2: R^2 overflows for
/// R > 1.3e154 while R is representable, which would discard exactly the heavy-tail draws
/// the radial diagnostic exists to test.
inline double logWFromLogR(double log_r)
{
    if (log_r == std::numeric_limits<double>::infinity())
        return -std::numeric_limits<double>::infinity();
    if (log_r == -std::numeric_limits<double>::infinity())
        return 0.0;
    return -logAddExpD(0.0, 2.0 * log_r);
}

// ---------------------------------------------------------------------------
// log of the regularized incomplete beta I_w(a, b).
//
// The radial diagnostic of record is Z = -log I_W(a, 3/2), which is standard exponential
// when W ~ Beta(a, 3/2) whatever the shape, and which stays finite at shapes where W itself
// rounds to zero.  Evaluating it needs log I_w from log w, so the evaluator is two-piece:
// a continued fraction where w is representable, and the leading small-w series below the
// frozen switch point, where w has underflowed but a log w has not.  The switch is at
// config/protocol.json's `log_q_switch`, carried in by the generated header.
//
// The analysis has its own evaluator, validated against an arbitrary-precision incomplete
// beta on the declared grid under gate G0.  This one is independent of it on purpose: the
// probe emits both the sufficient statistics of Z and the raw order statistics of log R
// from which the analysis can recompute Z with its own evaluator, so the two can be
// compared rather than assumed equal.
// ---------------------------------------------------------------------------
inline double logBetaFn(double a, double b)
{
    return std::lgamma(a) + std::lgamma(b) - std::lgamma(a + b);
}

/// Modified Lentz evaluation of the continued fraction of the incomplete beta, as in the
/// standard formulation; converges for x < (a+1)/(a+b+2).
inline double betaContinuedFraction(double x, double a, double b)
{
    const double kTiny= 1e-300;
    const double qab= a + b, qap= a + 1.0, qam= a - 1.0;
    double c= 1.0;
    double d= 1.0 - qab * x / qap;
    if (std::fabs(d) < kTiny)
        d= kTiny;
    d= 1.0 / d;
    double h= d;
    for (int m= 1; m <= 300; ++m)
    {
        const int m2= 2 * m;
        double aa= m * (b - m) * x / ((qam + m2) * (a + m2));
        d= 1.0 + aa * d;
        if (std::fabs(d) < kTiny)
            d= kTiny;
        c= 1.0 + aa / c;
        if (std::fabs(c) < kTiny)
            c= kTiny;
        d= 1.0 / d;
        h*= d * c;
        aa= -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
        d= 1.0 + aa * d;
        if (std::fabs(d) < kTiny)
            d= kTiny;
        c= 1.0 + aa / c;
        if (std::fabs(c) < kTiny)
            c= kTiny;
        d= 1.0 / d;
        const double del= d * c;
        h*= del;
        if (std::fabs(del - 1.0) < 3e-16)
            break;
    }
    return h;
}

/// log I_w(a, b) from log w.  Returns -inf for w = 0 and 0 for w = 1.
inline double logRegularizedIncompleteBeta(double log_w, double a, double b)
{
    if (!(log_w == log_w))
        return std::numeric_limits<double>::quiet_NaN();
    if (log_w == -std::numeric_limits<double>::infinity())
        return -std::numeric_limits<double>::infinity();
    if (log_w >= 0.0)
        return 0.0;

    if (log_w < protocol::kLogQSwitch)
    {
        // Small-w series: I_w = w^a / (a B(a,b)) * [1 + a(1-b)/(a+1) w + O(w^2)].  At the
        // switch point the correction is below 1e-13 relative, and it is kept only so that
        // the two pieces join smoothly rather than with a visible step.
        const double w= std::exp(log_w);
        return a * log_w - std::log(a) - logBetaFn(a, b) +
               std::log1p(a * (1.0 - b) / (a + 1.0) * w);
    }

    const double w= std::exp(log_w);
    const double lbt= a * log_w + b * std::log1p(-w) - logBetaFn(a, b);
    if (w < (a + 1.0) / (a + b + 2.0))
        return lbt - std::log(a) + std::log(betaContinuedFraction(w, a, b));
    const double q= std::exp(lbt - std::log(b)) * betaContinuedFraction(1.0 - w, b, a);
    if (q >= 1.0)
        return -std::numeric_limits<double>::infinity();
    return std::log1p(-q);
}

/// Z = -log I_W(a, 3/2), standard exponential under the target law.
inline double zFromLogR(double log_r, double a)
{
    const double log_w= logWFromLogR(log_r);
    const double li= logRegularizedIncompleteBeta(log_w, a, 1.5);
    if (li == -std::numeric_limits<double>::infinity())
        return std::numeric_limits<double>::infinity();
    return -li;
}

// ---------------------------------------------------------------------------
// Geometry.  g is the order-unity vector that multiplies R to give the final velocity,
// already rotated into the global frame.  Building it before exponentiating is the whole
// reason the log path can decide representability of the *returned* components rather than
// of the local ones.
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

    /// Loop invariants of the capped branch, hoisted exactly as the released class hoists
    /// them, so that the replica matches its cost profile and not only its arithmetic.
    /// `log_g_bound` is the rotation-invariant bound
    /// max_j |g_j| <= |g| <= sqrt(kappa) max(theta_perp, theta_par).
    T log_sqrt_kappa;
    T log_g_bound;
    T log_max_finite;

    Geometry(T k, T tperp, T tpar, const std::array<T, 3> &b)
        : kappa(k), sqrt_kappa(std::sqrt(k)), theta_perp(tperp), theta_par(tpar), ub(b),
          identity_rotation(b[0] == T(0) && b[1] == T(0) && b[2] == T(1)),
          rotator(k, tperp, tpar, b, bi_kappa_distribution<T>::no_cap()),
          log_sqrt_kappa(std::log(std::sqrt(k))),
          log_g_bound(std::log(std::sqrt(k)) + std::log(std::max(tperp, tpar))),
          log_max_finite(std::log(std::numeric_limits<T>::max()))
    {
    }

    /// Unit direction from an angle pair, in the field-aligned frame.  This is release
    /// 1.0.0's construction and the paired diagnostic layer's declared direction primitive;
    /// release 2.0.0 draws its direction by disc rejection instead and does not come
    /// through here.  It evaluates sin(phi) and cos(phi) of one argument, which a compiler
    /// may fuse into a combined routine whose sine differs by an ulp -- and it must keep
    /// doing so, because reproducing 1.0.0 exactly means reproducing that too.
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

inline void emitDoubleValue(FILE *f, double v)
{
    if (v == std::numeric_limits<double>::infinity())
        std::fprintf(f, "\"inf\"");
    else if (v == -std::numeric_limits<double>::infinity())
        std::fprintf(f, "\"-inf\"");
    else if (v != v)
        std::fprintf(f, "\"nan\"");
    else
        std::fprintf(f, "%.17g", v);
}

// ---------------------------------------------------------------------------
// Frozen constants, spelled out in the committed source so that a run is determined by the
// command line recorded in README.md and nothing else.  Every one of them is compared
// against config/protocol.json at start-up by checkAgainstProtocol(); the two statements of
// the matrix check each other instead of one being trusted.
// ---------------------------------------------------------------------------

/// PROTOCOL.md Sec. 4.  1.25 and 1.49 are new relative to Experiment 6 and are load-bearing:
/// they sample the region where the rejected primitive's cost diverges, and 1.49 and 1.5
/// straddle the shape a = 1.
inline std::vector<double> kappaLadder()
{
    std::vector<double> k;
    k.push_back(0.5001); k.push_back(0.501); k.push_back(0.505); k.push_back(0.51);
    k.push_back(0.55);   k.push_back(0.60);  k.push_back(0.75);  k.push_back(1.0);
    k.push_back(1.25);   k.push_back(1.49);  k.push_back(1.5);   k.push_back(2.0);
    k.push_back(5.0);
    return k;
}

/// PROTOCOL.md Sec. 3.  Declared, never derived.  Experiment 6 obtained five of its
/// performance seeds as `4001 + block`, which is why 4006-4010 appear in its manifest and in
/// no declaration; both blocks here are written out, and no seed anywhere in this
/// experiment is computed from another.
///
/// Protocol 2.0.0 governs the SECOND holdout and uses 8001-8010.  The block follows from the
/// rule in PROTOCOL.md Sec. 3 rather than from a choice: the highest seed declared anywhere
/// in this repository is 7505, rounded up to the next multiple of 1000 gives 8000, and the
/// next ten integers are taken.  7001-7010 belong to the first holdout, whose NO-GO result
/// is preserved in commit 45d3ef8; Sec. 8 forbids recomputing any result on them.
inline std::vector<unsigned> productionSeeds()
{
    std::vector<unsigned> s;
    s.push_back(8001u); s.push_back(8002u); s.push_back(8003u); s.push_back(8004u);
    s.push_back(8005u);
    return s;
}

inline std::vector<unsigned> performanceSeeds()
{
    std::vector<unsigned> s;
    s.push_back(8006u); s.push_back(8007u); s.push_back(8008u); s.push_back(8009u);
    s.push_back(8010u);
    return s;
}

/// The first holdout's block.  Declared so that a rerun on it is a compile-time-visible
/// mistake rather than an accident, and so that `spentSeeds()` can be asserted disjoint from
/// the block actually in use.
inline std::vector<unsigned> spentSeeds()
{
    std::vector<unsigned> s;
    for (unsigned v= 7001u; v <= 7010u; ++v)
        s.push_back(v);
    return s;
}

/// Fixture seeds for `make selftest`.  They write no data and no result depends on them,
/// but they are declared and kept disjoint from both blocks above for a reason that matters:
/// the selftest asserts that CANDIDATE has zero avoidable loss at kappa = 0.501 in double,
/// which is a gate G1 predicate.  Running that check on a production seed before the
/// holdout is read would be reading the holdout.
inline std::vector<unsigned> selftestSeeds()
{
    std::vector<unsigned> s;
    s.push_back(7501u); s.push_back(7502u); s.push_back(7503u); s.push_back(7504u);
    s.push_back(7505u);
    return s;
}

/// Compare the declarations above against the frozen protocol.  Called at start-up by every
/// phase; a mismatch is a hard stop, never a warning.
inline void checkAgainstProtocol()
{
    bool ok= true;
    const std::vector<double> k= kappaLadder();
    if (static_cast<int>(k.size()) != protocol::kNumKappa)
        ok= false;
    else
        for (int i= 0; i < protocol::kNumKappa; ++i)
            if (k[static_cast<size_t>(i)] != protocol::kKappaLadder[i])
                ok= false;

    const std::vector<unsigned> ps= productionSeeds();
    if (static_cast<int>(ps.size()) != protocol::kNumProductionSeeds)
        ok= false;
    else
        for (int i= 0; i < protocol::kNumProductionSeeds; ++i)
            if (ps[static_cast<size_t>(i)] != protocol::kProductionSeeds[i])
                ok= false;

    const std::vector<unsigned> qs= performanceSeeds();
    if (static_cast<int>(qs.size()) != protocol::kNumPerformanceSeeds)
        ok= false;
    else
        for (int i= 0; i < protocol::kNumPerformanceSeeds; ++i)
            if (qs[static_cast<size_t>(i)] != protocol::kPerformanceSeeds[i])
                ok= false;

    // The blocks must also be pairwise disjoint, which is the property Sec. 3 actually
    // relies on, and the fixture seeds must touch neither of them.
    const std::vector<unsigned> fs= selftestSeeds();
    for (size_t i= 0; i < ps.size(); ++i)
        for (size_t j= 0; j < qs.size(); ++j)
            if (ps[i] == qs[j])
                ok= false;
    for (size_t i= 0; i < fs.size(); ++i)
        for (size_t j= 0; j < ps.size(); ++j)
            if (fs[i] == ps[j])
                ok= false;
    for (size_t i= 0; i < fs.size(); ++i)
        for (size_t j= 0; j < qs.size(); ++j)
            if (fs[i] == qs[j])
                ok= false;

    if (!ok)
    {
        std::fprintf(stderr,
                     "exp7: the matrix declared in src/exp7_common.H disagrees with\n"
                     "      config/protocol.json (%s).  Refusing to run.\n",
                     protocol::kProtocolSha256);
        std::exit(3);
    }
}

// ---------------------------------------------------------------------------
// The fixed grid on which the empirical CDF of Z is reported.
//
// Z is standard exponential under the target law, so 25 natural-log units is a survival
// probability of 1.4e-11: the grid covers the whole range in which 10^6 draws can place
// anything, and the exceedance counts above the three protocol thresholds are reported
// exactly rather than read off the grid.
// ---------------------------------------------------------------------------
static const int kZGridPoints= 501;
static const double kZGridStep= 0.05;
inline double zGridPoint(int i) { return kZGridStep * static_cast<double>(i); }

} // namespace exp7

#endif // EXP7_COMMON_H
