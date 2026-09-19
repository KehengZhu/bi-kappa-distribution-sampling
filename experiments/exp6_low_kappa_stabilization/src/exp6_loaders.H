// Experiment 6 - the three radius formations and the complete loaders built on them.
//
// Strict C++11.
//
// Two layers live here and must never be conflated:
//
//   PAIRED DIAGNOSTIC LAYER (runPairedAttempt)
//       QF, SPLIT and the log-domain propagation all receive the *same* declared
//       primitives (x1, y, u, cos_theta, phi), so one mathematical attempt can be followed
//       through every formation and "was this draw recoverable?" is decided per draw
//       instead of by comparing histograms.  This is a mechanism experiment.  It is NOT a
//       draw-by-draw reconstruction of any standard library's Gamma generator: it declares
//       its own primitives and says so.
//
//   NATIVE LAYER (runNativeAttempt)
//       The released loader runs against its own std::gamma_distribution, and the log
//       loader runs against its own log primitive.  This is what a user actually receives.
//       It cannot separate avoidable from honest loss on a draw where the denominator
//       materialized as zero, because the intended value is then unknowable from the
//       output alone.  That separation is what the paired layer is for.

#ifndef EXP6_LOADERS_H
#define EXP6_LOADERS_H

#include "exp6_common.H"
#include "log_gamma_identity.H"
#include "log_gamma_published.H"

#include <cmath>
#include <cstdint>
#include <limits>
#include <random>

namespace exp6
{

enum Method
{
    kMethodQF= 0,
    kMethodSPLIT= 1,
    kMethodLOG= 2,
    kNumMethods= 3
};

inline const char *methodName(int m)
{
    switch (m)
    {
    case kMethodQF: return "QF";
    case kMethodSPLIT: return "SPLIT";
    case kMethodLOG: return "LOG";
    default: return "unknown";
    }
}

// ---------------------------------------------------------------------------
// Selectable log-domain primitive
// ---------------------------------------------------------------------------
template <typename T> class LogGammaEngine
{
  public:
    enum Kind
    {
        kIdentity= 0,  ///< LOG-ID: log X = log Y + log(U)/a
        kPublished= 1  ///< LOG-PUB: transcribed published small-shape log-scale generator
    };

    LogGammaEngine(Kind kind, T a)
        : m_kind(kind), m_a(a), m_gBoost(a + T(1), T(1)), m_unif(T(0), T(1)), m_pub(a)
    {
    }

    static const char *kindName(Kind k) { return k == kIdentity ? "LOG-ID" : "LOG-PUB"; }
    Kind kind() const { return m_kind; }

    template <typename Gen> LogGammaDraw<T> operator()(Gen &gen)
    {
        if (m_kind == kIdentity)
            return drawLogGammaIdentity<T>(m_a, gen, m_gBoost, m_unif);
        return m_pub.draw(gen, m_unif);
    }

  private:
    Kind m_kind;
    T m_a;
    std::gamma_distribution<T> m_gBoost;
    std::uniform_real_distribution<T> m_unif;
    LogGammaPublished<T> m_pub;
};

// ---------------------------------------------------------------------------
// The declared primitives of one paired attempt
// ---------------------------------------------------------------------------
template <typename T> struct SharedPrimitives
{
    T x1;
    T y;
    T u;
    T cos_theta;
    T phi;
    T x2;      ///< the product form y * u^(1/a), i.e. what a library would materialize
    T log_x1;  ///< working-precision log
    T log_x2;  ///< working-precision log via the identity; finite when the draw exists
    uint32_t retries;
};

template <typename T, typename Gen>
SharedPrimitives<T> drawShared(T a, Gen &gen, std::gamma_distribution<T> &gX1,
                               std::gamma_distribution<T> &gBoost,
                               std::uniform_real_distribution<T> &unif,
                               std::uniform_real_distribution<T> &uCos,
                               std::uniform_real_distribution<T> &uPhi)
{
    SharedPrimitives<T> s;
    s.x1= gX1(gen);
    const LogGammaDraw<T> d= drawLogGammaIdentity<T>(a, gen, gBoost, unif);
    s.y= d.y;
    s.u= d.u;
    s.retries= d.retries;
    s.x2= d.x_product;
    s.log_x2= d.log_x;
    s.log_x1= std::log(s.x1);
    s.cos_theta= uCos(gen);
    s.phi= uPhi(gen);
    return s;
}

// ---------------------------------------------------------------------------
// The intended sample, evaluated in double from the recorded primitives.
//
// This is the reference against which every working-precision formation is judged.  It is
// a *reference*, not an oracle: for a float run it is over thirty bits more accurate than
// the thing it judges, but for a double run it shares the type, so draws that land near a
// representability threshold are sent to the independent 100-digit oracle instead of being
// decided here.
// ---------------------------------------------------------------------------
struct Reference
{
    double log_x1;
    double log_x2;
    double log_r;
    double log_w;             ///< log W, W = X2/(X1+X2) ~ Beta(a, 3/2)
    double g[3];              ///< order-unity vector: V = R g, already rotated
    double log_component[3];  ///< log |V_j|
    double max_log_component;
    double log_speed;         ///< log |V|
    double log_local_max;     ///< log of the largest *local*-frame component
    bool representable;       ///< every returned component is representable in the run type
    bool local_representable; ///< the pre-rotation local vector is representable
    bool near_limit;          ///< within kOracleBandLogUnits of the overflow threshold
};

template <typename T>
Reference makeReference(double a, const SharedPrimitives<T> &s, const Geometry<T> &geom,
                        double log_ovf)
{
    Reference r;
    r.log_x1= std::log(static_cast<double>(s.x1));
    r.log_x2= logGammaIdentityReference(a, static_cast<double>(s.y), static_cast<double>(s.u));
    r.log_r= 0.5 * (r.log_x1 - r.log_x2);
    r.log_w= r.log_x2 - logAddExpD(r.log_x1, r.log_x2);

    // Geometry in double, from the same two uniforms.
    const double ct= static_cast<double>(s.cos_theta);
    const double st= std::sqrt(std::max(0.0, 1.0 - ct * ct));
    const double ph= static_cast<double>(s.phi);
    const double n[3]= {st * std::cos(ph), st * std::sin(ph), ct};
    const double sk= std::sqrt(static_cast<double>(geom.kappa));
    const double tp= static_cast<double>(geom.theta_perp);
    const double tq= static_cast<double>(geom.theta_par);
    const double local[3]= {sk * tp * n[0], sk * tp * n[1], sk * tq * n[2]};

    double loc_max= -std::numeric_limits<double>::infinity();
    for (int j= 0; j < 3; ++j)
    {
        const double lv= r.log_r + std::log(std::fabs(local[j]));
        if (lv > loc_max)
            loc_max= lv;
    }
    r.log_local_max= loc_max;

    if (geom.identity_rotation)
    {
        for (int j= 0; j < 3; ++j)
            r.g[j]= local[j];
    }
    else
    {
        // Same basis construction as the released header, evaluated in double.
        std::array<double, 3> loc3;
        loc3[0]= local[0];
        loc3[1]= local[1];
        loc3[2]= local[2];
        std::array<double, 3> ubd;
        ubd[0]= static_cast<double>(geom.ub[0]);
        ubd[1]= static_cast<double>(geom.ub[1]);
        ubd[2]= static_cast<double>(geom.ub[2]);
        bi_kappa_distribution<double> rotD(static_cast<double>(geom.kappa), tp, tq, ubd,
                                           bi_kappa_distribution<double>::no_cap());
        const std::array<double, 3> gg= rotD.rotate_from_fieldAligned_frame(loc3, ubd);
        for (int j= 0; j < 3; ++j)
            r.g[j]= gg[j];
    }

    double mx= -std::numeric_limits<double>::infinity();
    double sum_sq_log= 0.0;
    (void)sum_sq_log;
    for (int j= 0; j < 3; ++j)
    {
        r.log_component[j]= r.log_r + std::log(std::fabs(r.g[j]));
        if (r.log_component[j] > mx)
            mx= r.log_component[j];
    }
    r.max_log_component= mx;

    const double gnorm= std::sqrt(r.g[0] * r.g[0] + r.g[1] * r.g[1] + r.g[2] * r.g[2]);
    r.log_speed= r.log_r + std::log(gnorm);

    r.representable= (mx < log_ovf);
    r.local_representable= (loc_max < log_ovf);
    r.near_limit= std::fabs(mx - log_ovf) <= kOracleBandLogUnits ||
                  std::fabs(loc_max - log_ovf) <= kOracleBandLogUnits;
    return r;
}

// ---------------------------------------------------------------------------
// Released-style vector construction from a materialized radius.
//
// Deliberately identical in operation order to bi_kappa_distribution::operator(): the
// unit direction is multiplied by the radius first, the thermal scaling is applied second,
// and the rotation comes last.  That order is what makes a rotated configuration lose
// draws whose *returned* components are representable, so it must not be "improved" here.
// ---------------------------------------------------------------------------
template <typename T>
std::array<T, 3> vectorFromRadius(T r, T cos_theta, T phi, const Geometry<T> &geom)
{
    const T sin_theta= std::sqrt(std::max(T(0), T(1) - cos_theta * cos_theta));
    const T ux= r * sin_theta * std::cos(phi);
    const T uy= r * sin_theta * std::sin(phi);
    const T uz= r * cos_theta;

    std::array<T, 3> p;
    p[0]= geom.sqrt_kappa * geom.theta_perp * ux;
    p[1]= geom.sqrt_kappa * geom.theta_perp * uy;
    p[2]= geom.sqrt_kappa * geom.theta_par * uz;
    if (geom.identity_rotation)
        return p;
    return geom.rotator.rotate_from_fieldAligned_frame(p, geom.ub);
}

/// Log-domain vector construction: build the order-unity vector first, decide
/// representability from log|V_j| = log R + log|g_j|, and only then exponentiate.
template <typename T>
std::array<T, 3> vectorFromLogRadius(T log_r, T cos_theta, T phi, const Geometry<T> &geom,
                                     bool *any_overflow)
{
    const std::array<T, 3> n= Geometry<T>::direction(cos_theta, phi);
    const std::array<T, 3> g= geom.gVector(n);

    const T log_ovf= static_cast<T>(logOverflowThreshold<T>());
    const T log_max= static_cast<T>(logMaxFinite<T>());

    *any_overflow= false;
    std::array<T, 3> v;
    const bool radius_representable= (log_r <= log_max);
    const T R= radius_representable ? std::exp(log_r) : T(0);

    for (int j= 0; j < 3; ++j)
    {
        if (g[j] == T(0))
        {
            v[j]= T(0);
            continue;
        }
        const T lv= log_r + std::log(std::fabs(g[j]));
        if (lv >= log_ovf)
        {
            *any_overflow= true;
            v[j]= std::copysign(std::numeric_limits<T>::infinity(), g[j]);
            continue;
        }
        // Multiplying by the exponentiated radius keeps a full working-precision product
        // wherever the radius itself is representable; the per-component exponential is
        // used only where it is not.
        v[j]= radius_representable ? R * g[j]
                                   : std::copysign(std::exp(lv), g[j]);
    }
    return v;
}

/// Cap predicate in the log domain, in the field-aligned frame, matching the released
/// predicate |p_j| / theta_j <= lambda exactly:  sqrt(kappa) R |n_j| <= lambda.
template <typename T>
bool capAcceptLog(T log_r, T cos_theta, T phi, const Geometry<T> &geom, T log_cap)
{
    const std::array<T, 3> n= Geometry<T>::direction(cos_theta, phi);
    const T lsk= std::log(geom.sqrt_kappa);
    for (int j= 0; j < 3; ++j)
    {
        if (n[j] == T(0))
            continue;
        if (lsk + log_r + std::log(std::fabs(n[j])) > log_cap)
            return false;
    }
    return true;
}

// ---------------------------------------------------------------------------
// Per-method verdict on one attempt
// ---------------------------------------------------------------------------
inline bool allFinite(const double v[3])
{
    for (int j= 0; j < 3; ++j)
        if (!(v[j] == v[j]) || v[j] == std::numeric_limits<double>::infinity() ||
            v[j] == -std::numeric_limits<double>::infinity())
            return false;
    return true;
}

template <typename T> bool allFiniteT(const std::array<T, 3> &v)
{
    for (int j= 0; j < 3; ++j)
        if (!std::isfinite(v[j]))
            return false;
    return true;
}

/// Terminal category for one method on one attempt.  Exactly one category is returned;
/// the overlapping descriptive flags live separately.
inline int classifyTerminal(int method, bool method_finite, bool ref_representable,
                            uint32_t flags)
{
    if (!ref_representable)
        return method_finite ? kCatFinite : kCatHonestOverflow;
    if (method_finite)
        return kCatFinite;
    if (method == kMethodLOG)
        return kCatLogPrimFail;
    if (flags & kFlagX2Zero)
        return kCatDenomZero;
    return (method == kMethodQF) ? kCatQuotientLoss : kCatSplitLoss;
}

// ---------------------------------------------------------------------------
// Counters
// ---------------------------------------------------------------------------
struct MethodCounters
{
    long long cat[kNumCategories];
    long long nonfinite_output;
    double max_finite_log_component;

    MethodCounters()
        : nonfinite_output(0),
          max_finite_log_component(-std::numeric_limits<double>::infinity())
    {
        for (int i= 0; i < kNumCategories; ++i)
            cat[i]= 0;
    }

    long long total() const
    {
        long long t= 0;
        for (int i= 0; i < kNumCategories; ++i)
            t+= cat[i];
        return t;
    }
    long long avoidable() const
    {
        long long t= 0;
        for (int i= 0; i < kNumCategories; ++i)
            if (isAvoidable(i))
                t+= cat[i];
        return t;
    }
};

struct AttemptCounters
{
    long long n;
    long long x2_zero;
    long long x2_subnormal;
    long long u_endpoint_redraw;
    long long ref_nonrepresentable;
    long long rotation_recoverable;
    long long near_limit;
    long long audited;
    long long out_of_domain;
    double max_log_r_ref;
    // Accuracy of the log domain relative to the split form, on draws where both resolve.
    long long n_accuracy;
    double sum_abs_rel_err;
    double max_abs_rel_err;
    long long stratum_total[kNumAuditStrata];
    long long stratum_audited[kNumAuditStrata];

    AttemptCounters()
        : n(0), x2_zero(0), x2_subnormal(0), u_endpoint_redraw(0), ref_nonrepresentable(0),
          rotation_recoverable(0), near_limit(0), audited(0), out_of_domain(0),
          max_log_r_ref(-std::numeric_limits<double>::infinity()), n_accuracy(0),
          sum_abs_rel_err(0.0), max_abs_rel_err(0.0)
    {
        for (int i= 0; i < kNumAuditStrata; ++i)
        {
            stratum_total[i]= 0;
            stratum_audited[i]= 0;
        }
    }
};

} // namespace exp6

#endif // EXP6_LOADERS_H
