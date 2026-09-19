// Experiment 7 - the two radius formations under test, and the loaders built on them.
//
// Strict C++11.
//
// Two methods, and nothing else:
//
//   LEGACY     the shipped 1.0.0 loader.  R = sqrt(X1)/sqrt(X2) with X2 materialized by
//              std::gamma_distribution, the unit direction multiplied by R, the thermal
//              scaling applied, and the rotation last.
//   CANDIDATE  the released 2.0.0 loader.  log X2 = log Y + log(U)/a with Y ~ Ga(a+1,1)
//              (Ahrens-Dieter shape boosting), log R = (log X1 - log X2)/2, the order-unity
//              vector g = Q(ub) diag(sqrt(kappa) theta) n built first, and each component's
//              representability decided from log|V_j| = log R + log|g_j| before anything is
//              exponentiated.
//
// QF (`sqrt(X1/X2)`, the pre-1.0 formation whose intermediate quotient overflows even where
// its square root is representable) survives only as a third *diagnostic column of the
// paired layer*.  It is not a method, it is never run natively, and no claim rests on it; it
// is kept because it costs one square root per attempt and it is the only way to show that
// the 1.0.0 split form was itself already a repair.
//
// Two layers live here and must never be conflated:
//
//   PAIRED DIAGNOSTIC LAYER (runPairedAttempt)
//       LEGACY, CANDIDATE and QF all receive the *same* declared primitives
//       (x1, y, u, cos_theta, phi), so one mathematical attempt can be followed through
//       every formation and "was this draw recoverable?" is decided per draw instead of by
//       comparing histograms.  This is a mechanism experiment.  It is NOT a draw-by-draw
//       reconstruction of any standard library's Gamma generator: it declares its own
//       primitives and says so.  Those primitives are the candidate's own, which has the
//       useful side effect that the paired layer is itself identical across standard
//       libraries -- the mechanism decomposition does not depend on which library ran it.
//
//   NATIVE LAYER (legacyNativeAttempt / candidateNativeAttempt)
//       Each method runs its own real implementation end to end, against its own
//       primitives.  This is what a user receives.  It cannot separate avoidable from
//       honest loss on a draw whose denominator materialized as zero, because the intended
//       value is then unknowable from the output alone.  That separation is what the paired
//       layer is for.  `make selftest` checks each native replica against the released class
//       it stands for, bit for bit, which is what licenses using the replica wherever an
//       intermediate has to be observable.

#ifndef EXP7_LOADERS_H
#define EXP7_LOADERS_H

#include "exp7_common.H"

#include <cmath>
#include <cstdint>
#include <limits>
#include <random>

namespace exp7
{

enum Method
{
    kMethodLegacy= 0,
    kMethodCandidate= 1,
    kMethodQF= 2, ///< paired-layer diagnostic column only; never a native method
    kNumMethods= 3
};

inline const char *methodName(int m)
{
    switch (m)
    {
    case kMethodLegacy: return "LEGACY";
    case kMethodCandidate: return "CANDIDATE";
    case kMethodQF: return "QF";
    default: return "unknown";
    }
}

/// True for the two methods the protocol puts under test.  QF is excluded everywhere a
/// verdict is formed.
inline bool isTestedMethod(int m) { return m == kMethodLegacy || m == kMethodCandidate; }

// ---------------------------------------------------------------------------
// The declared primitives of one paired attempt.
//
// Drawn with the candidate's own primitives from bikappa_detail, so that the paired layer's
// verdicts are a property of the arithmetic rather than of the standard library that ran
// it.  x2 is the product form y * u^(1/a) -- what a library-backed Gamma(a) call would
// materialize, and what underflows.
// ---------------------------------------------------------------------------
template <typename T> struct SharedPrimitives
{
    T x1;
    T y;
    T u;
    T cos_theta;
    T phi;
    T x2;     ///< the product form y * u^(1/a)
    T log_x1; ///< working-precision log
    T log_x2; ///< working-precision log via the identity; finite when the draw exists
};

template <typename T, typename Gen>
SharedPrimitives<T> drawShared(T a, Gen &gen, T &ncache, bool &ncached)
{
    SharedPrimitives<T> s;
    s.x1= bikappa_detail::marsagliaTsangGamma<T, Gen>(T(1.5), gen, ncache, ncached);
    s.y= bikappa_detail::marsagliaTsangGamma<T, Gen>(a + T(1), gen, ncache, ncached);
    s.u= bikappa_detail::openCanonical<T, Gen>(gen);
    s.x2= s.y * std::pow(s.u, T(1) / a);
    s.log_x1= std::log(s.x1);
    s.log_x2= std::log(s.y) + std::log(s.u) / a;
    s.cos_theta= T(2) * bikappa_detail::openCanonical<T, Gen>(gen) - T(1);
    s.phi= T(2) * piValue<T>() * bikappa_detail::openCanonical<T, Gen>(gen);
    return s;
}

// ---------------------------------------------------------------------------
// The intended sample, evaluated in double from the recorded primitives.
//
// This is the reference against which every working-precision formation is judged.  It is a
// *reference*, not an oracle: for a float run it is over thirty bits more accurate than the
// thing it judges, and for a double run it shares the type but not the formation -- the
// reference never materializes X2, so it is unaffected by the underflow that degrades the
// legacy radius, and its own error is about eps |log R|, which is 1.6e-13 at the double
// overflow threshold against a 1.05e-8 accuracy threshold.  Draws that land near a
// representability threshold are nevertheless sent to the independent 100-digit oracle
// rather than being decided here.
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
    r.log_x2= std::log(static_cast<double>(s.y)) + std::log(static_cast<double>(s.u)) / a;
    r.log_r= 0.5 * (r.log_x1 - r.log_x2);
    r.log_w= logWFromLogR(r.log_r);

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
// LEGACY vector construction from a materialized radius.
//
// Deliberately identical in operation order to bikappa_v1::bi_kappa_distribution's
// operator(): the unit direction is multiplied by the radius first, the thermal scaling is
// applied second, and the rotation comes last.  That order is what makes a rotated
// configuration lose draws whose *returned* components are representable, so it must not be
// "improved" here.
// ---------------------------------------------------------------------------
template <typename T>
std::array<T, 3> legacyVectorFromRadius(T r, T cos_theta, T phi, const Geometry<T> &geom)
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

/// The LEGACY cap predicate, applied where LEGACY applies it: to the pre-rotation local
/// vector, in linear arithmetic, after the radius has been materialized.
template <typename T>
bool legacyCapAccept(T r, T cos_theta, T phi, const Geometry<T> &geom, T lambda)
{
    const std::array<T, 3> n= Geometry<T>::direction(cos_theta, phi);
    const T p0= geom.sqrt_kappa * geom.theta_perp * (r * n[0]);
    const T p1= geom.sqrt_kappa * geom.theta_perp * (r * n[1]);
    const T p2= geom.sqrt_kappa * geom.theta_par * (r * n[2]);
    return std::fabs(p0) / geom.theta_perp <= lambda &&
           std::fabs(p1) / geom.theta_perp <= lambda &&
           std::fabs(p2) / geom.theta_par <= lambda;
}

/// CANDIDATE vector construction, reproducing bi_kappa_distribution 2.0.0 step for step:
/// build the order-unity vector first, take the fast path when log R + log gBound is safely
/// below the type's limit, and otherwise decide each component from log|V_j| before
/// exponentiating.  The predicate is the released one -- log(max()) -- not the exact
/// overflow threshold, because that is the predicate a user actually gets.
template <typename T>
std::array<T, 3> candidateVectorFromLogRadius(T log_r, T cos_theta, T phi,
                                              const Geometry<T> &geom, bool *any_overflow)
{
    const std::array<T, 3> n= Geometry<T>::direction(cos_theta, phi);
    const std::array<T, 3> g= geom.gVector(n);

    const T log_max= static_cast<T>(std::log(std::numeric_limits<T>::max()));
    const T log_g_bound= std::log(geom.sqrt_kappa) +
                         std::log(std::max(geom.theta_perp, geom.theta_par));

    *any_overflow= false;
    const bool radius_representable= (log_r <= log_max);
    const T R= radius_representable ? std::exp(log_r) : T(0);

    std::array<T, 3> v;
    if (radius_representable && log_r + log_g_bound < log_max)
    {
        v[0]= R * g[0];
        v[1]= R * g[1];
        v[2]= R * g[2];
        return v;
    }
    for (int j= 0; j < 3; ++j)
    {
        if (g[j] == T(0))
        {
            v[j]= T(0);
            continue;
        }
        const T lv= log_r + std::log(std::fabs(g[j]));
        if (lv >= log_max)
        {
            *any_overflow= true;
            v[j]= (g[j] < T(0)) ? -std::numeric_limits<T>::infinity()
                                : std::numeric_limits<T>::infinity();
            continue;
        }
        if (radius_representable)
            v[j]= R * g[j];
        else
        {
            const T mag= std::exp(lv);
            v[j]= (g[j] < T(0)) ? -mag : mag;
        }
    }
    return v;
}

/// The CANDIDATE cap predicate, in the log domain, in the field-aligned frame.  Exactly the
/// released test: only the largest |n_j| has to be examined, so it costs one logarithm.
template <typename T>
bool candidateCapAcceptLog(T log_r, T cos_theta, T phi, const Geometry<T> &geom, T log_cap,
                           bool *unrepresentable)
{
    const std::array<T, 3> n= Geometry<T>::direction(cos_theta, phi);
    const T max_abs_n= std::max(std::fabs(n[0]), std::max(std::fabs(n[1]), std::fabs(n[2])));
    const T log_largest= log_r + std::log(geom.sqrt_kappa) + std::log(max_abs_n);
    if (log_largest > log_cap)
    {
        if (unrepresentable)
            *unrepresentable= (log_largest >= static_cast<T>(std::log(
                                  std::numeric_limits<T>::max())));
        return false;
    }
    if (unrepresentable)
        *unrepresentable= false;
    return true;
}

// ---------------------------------------------------------------------------
// Verdicts
// ---------------------------------------------------------------------------
template <typename T> bool allFiniteT(const std::array<T, 3> &v)
{
    for (int j= 0; j < 3; ++j)
        if (!std::isfinite(v[j]))
            return false;
    return true;
}

/// Relative error of a formation's radius against the reference, evaluated in the log domain
/// so that nothing overflows: |r/R - 1| = |expm1(log r - log R)|.  Returns NaN when the
/// comparison is not defined, and 1 when the two are more than a factor of e apart, which is
/// beyond any accuracy threshold and would otherwise overflow expm1.
inline double relativeRadiusError(double log_r_method, double log_r_ref)
{
    if (!(log_r_method == log_r_method) || !(log_r_ref == log_r_ref))
        return std::numeric_limits<double>::quiet_NaN();
    if (log_r_method == std::numeric_limits<double>::infinity() ||
        log_r_method == -std::numeric_limits<double>::infinity() ||
        log_r_ref == std::numeric_limits<double>::infinity() ||
        log_r_ref == -std::numeric_limits<double>::infinity())
        return std::numeric_limits<double>::quiet_NaN();
    const double d= log_r_method - log_r_ref;
    if (d > 1.0 || d < -1.0)
        return 1.0;
    return std::fabs(std::expm1(d));
}

/// Terminal category for one method on one attempt.  Exactly one category is returned; the
/// overlapping descriptive flags live separately.
///
/// `rel_err` is the relative radius error against the reference, NaN where it is not
/// defined, and `max_rel_err` the type's FINITE_BUT_WRONG threshold from the protocol.  A
/// finite return is a success only when it is also accurate: amendment 1.1.0 exists because
/// Experiment 6's classifier scored `method_finite ? finite : honest_overflow` and so
/// recorded a fabricated number as a success.
inline int classifyTerminal(int method, bool method_finite, bool ref_representable,
                            double rel_err, double max_rel_err, uint32_t flags)
{
    if (!ref_representable)
        return method_finite ? kCatOverflowReturnedFinite : kCatHonestOverflow;
    if (method_finite)
    {
        if (rel_err == rel_err && rel_err > max_rel_err)
            return kCatFiniteButWrong;
        return kCatFinite;
    }
    if (method == kMethodCandidate)
        return kCatLogPrimFail;
    if (flags & kFlagX2Zero)
        return kCatDenomZero;
    return (method == kMethodQF) ? kCatQuotientLoss : kCatLegacyLoss;
}

// ---------------------------------------------------------------------------
// One paired attempt: the same mathematical draw pushed through all three formations.
// ---------------------------------------------------------------------------
struct PairedResult
{
    Reference ref;
    uint32_t flags;
    int cat[kNumMethods];
    bool finite[kNumMethods];
    double rel_err[kNumMethods]; ///< |r_method / R_ref - 1|, NaN where not comparable
    double log_r_method[kNumMethods];
};

template <typename T>
PairedResult runPairedAttempt(double kappa, const SharedPrimitives<T> &s,
                              const Geometry<T> &geom, double log_ovf)
{
    PairedResult pr;
    pr.flags= 0u;

    const T x2= s.x2;
    if (x2 == T(0))
        pr.flags|= kFlagX2Zero;
    else if (std::fabs(x2) < std::numeric_limits<T>::min())
        pr.flags|= kFlagX2Subnormal;

    // The shape must be the one the run actually used: kappa - 1/2 formed in the working
    // precision, then widened.  Forming it in double instead shifts log X2 by roughly 0.03
    // at float kappa = 0.5001, which is enough to move draws across the representability
    // threshold and mis-split avoidable from honest loss.
    pr.ref= makeReference<T>(shapeReference<T>(kappa), s, geom, log_ovf);
    if (!pr.ref.representable)
        pr.flags|= kFlagHonestOverflowRef;
    if (pr.ref.near_limit)
        pr.flags|= kFlagNearLimit;
    if (pr.ref.representable && !pr.ref.local_representable)
        pr.flags|= kFlagRotationRecoverable;

    const T r_qf= std::sqrt(s.x1 / x2);
    const T r_legacy= std::sqrt(s.x1) / std::sqrt(x2);
    const T log_r_cand= T(0.5) * (s.log_x1 - s.log_x2);

    const std::array<T, 3> v_qf= legacyVectorFromRadius<T>(r_qf, s.cos_theta, s.phi, geom);
    const std::array<T, 3> v_legacy=
        legacyVectorFromRadius<T>(r_legacy, s.cos_theta, s.phi, geom);
    bool cand_overflow= false;
    const std::array<T, 3> v_cand= candidateVectorFromLogRadius<T>(log_r_cand, s.cos_theta,
                                                                   s.phi, geom, &cand_overflow);

    pr.finite[kMethodQF]= allFiniteT<T>(v_qf);
    pr.finite[kMethodLegacy]= allFiniteT<T>(v_legacy);
    pr.finite[kMethodCandidate]= !cand_overflow && allFiniteT<T>(v_cand);

    if (!std::isfinite(r_qf))
        pr.flags|= kFlagQfNonFinite;
    if (!std::isfinite(r_legacy))
        pr.flags|= kFlagLegacyNonFinite;
    if (!(log_r_cand == log_r_cand))
        pr.flags|= kFlagCandidateNonFinite;

    pr.log_r_method[kMethodQF]=
        (std::isfinite(r_qf) && r_qf > T(0)) ? std::log(static_cast<double>(r_qf))
                                             : std::numeric_limits<double>::quiet_NaN();
    pr.log_r_method[kMethodLegacy]=
        (std::isfinite(r_legacy) && r_legacy > T(0))
            ? std::log(static_cast<double>(r_legacy))
            : std::numeric_limits<double>::quiet_NaN();
    pr.log_r_method[kMethodCandidate]= static_cast<double>(log_r_cand);

    const double max_rel= maxRelErrorFor<T>();
    for (int m= 0; m < kNumMethods; ++m)
    {
        pr.rel_err[m]= relativeRadiusError(pr.log_r_method[m], pr.ref.log_r);
        pr.cat[m]= classifyTerminal(m, pr.finite[m], pr.ref.representable, pr.rel_err[m],
                                    max_rel, pr.flags);
    }
    if (pr.cat[kMethodLegacy] == kCatFinite)
        pr.flags|= kFlagLegacySuccess;
    if (pr.cat[kMethodCandidate] == kCatFinite)
        pr.flags|= kFlagCandidateSuccess;
    return pr;
}

// ---------------------------------------------------------------------------
// Native layer.
//
// Each of these reproduces one released class exactly, so that the intermediates the class
// hides -- the materialized denominator, the carried log radius, the variate counts -- are
// observable.  `make selftest` compares each against the class it stands for, component by
// component, over a run that reaches every branch.
// ---------------------------------------------------------------------------

/// What one native attempt produced, beyond the returned vector.
template <typename T> struct NativeAttempt
{
    std::array<T, 3> v;
    T x2;        ///< materialized denominator (LEGACY only; 0 for CANDIDATE, which never forms it)
    double log_r; ///< the radius the method carried or materialized, widened; +-inf/NaN possible
    bool finite;
    bool overflow; ///< the method itself detected an unrepresentable component
};

/// LEGACY, reproducing bikappa_v1::bi_kappa_distribution<T>::operator() including its
/// param-override Gamma call and its pre-rotation cap test.
template <typename T> class LegacyNative
{
  public:
    explicit LegacyNative(double kappa)
        : m_gX1(T(1.5), T(1)), m_gX2(T(1), T(1)), m_uCos(T(-1), T(1)),
          m_uPhi(T(0), T(2) * piValue<T>()),
          m_px2(shapeFromKappa<T>(kappa), T(1))
    {
    }

    template <typename Gen>
    NativeAttempt<T> attempt(Gen &gen, const Geometry<T> &geom)
    {
        NativeAttempt<T> a;
        const T x1= m_gX1(gen);
        const T x2= m_gX2(gen, m_px2);
        const T r= std::sqrt(x1) / std::sqrt(x2);
        const T ct= m_uCos(gen);
        const T ph= m_uPhi(gen);
        a.v= legacyVectorFromRadius<T>(r, ct, ph, geom);
        a.x2= x2;
        a.log_r= (std::isfinite(r) && r > T(0))
                     ? std::log(static_cast<double>(r))
                     : (r == std::numeric_limits<T>::infinity()
                            ? std::numeric_limits<double>::infinity()
                            : (r == T(0) ? -std::numeric_limits<double>::infinity()
                                         : std::numeric_limits<double>::quiet_NaN()));
        a.finite= allFiniteT<T>(a.v);
        a.overflow= !a.finite;
        m_lastCos= ct;
        m_lastPhi= ph;
        m_lastX1= x1;
        m_lastR= r;
        return a;
    }

    T lastCosTheta() const { return m_lastCos; }
    T lastPhi() const { return m_lastPhi; }
    T lastX1() const { return m_lastX1; }
    T lastRadius() const { return m_lastR; }

  private:
    std::gamma_distribution<T> m_gX1;
    std::gamma_distribution<T> m_gX2;
    std::uniform_real_distribution<T> m_uCos;
    std::uniform_real_distribution<T> m_uPhi;
    typename std::gamma_distribution<T>::param_type m_px2;
    T m_lastCos, m_lastPhi, m_lastX1, m_lastR;
};

/// CANDIDATE, reproducing bi_kappa_distribution<T>::operator() including the normal cache
/// the two Gamma calls share.  The cache is the reason this cannot be a free function: a
/// reseeded stream must start with `cached = false`, which is what 2.0.0's seed() now does
/// and 1.x did not.
template <typename T> class CandidateNative
{
  public:
    explicit CandidateNative(double kappa)
        : m_a(shapeFromKappa<T>(kappa)), m_ncache(T(0)), m_ncached(false)
    {
    }

    void reset()
    {
        m_ncache= T(0);
        m_ncached= false;
    }

    template <typename Gen>
    NativeAttempt<T> attempt(Gen &gen, const Geometry<T> &geom)
    {
        NativeAttempt<T> a;
        const T x1= bikappa_detail::marsagliaTsangGamma<T, Gen>(T(1.5), gen, m_ncache,
                                                                m_ncached);
        const T log_x2= bikappa_detail::logGamma<T, Gen>(m_a, gen, m_ncache, m_ncached);
        const T log_r= T(0.5) * (std::log(x1) - log_x2);
        const T ct= T(2) * bikappa_detail::openCanonical<T, Gen>(gen) - T(1);
        const T ph= T(2) * piValue<T>() * bikappa_detail::openCanonical<T, Gen>(gen);

        bool overflow= false;
        a.v= candidateVectorFromLogRadius<T>(log_r, ct, ph, geom, &overflow);
        a.x2= T(0);
        a.log_r= static_cast<double>(log_r);
        a.overflow= overflow;
        a.finite= !overflow && allFiniteT<T>(a.v);
        m_lastCos= ct;
        m_lastPhi= ph;
        m_lastX1= x1;
        m_lastLogX2= log_x2;
        return a;
    }

    T shape() const { return m_a; }
    T lastCosTheta() const { return m_lastCos; }
    T lastPhi() const { return m_lastPhi; }
    T lastX1() const { return m_lastX1; }
    T lastLogX2() const { return m_lastLogX2; }

  private:
    T m_a;
    T m_ncache;
    bool m_ncached;
    T m_lastCos, m_lastPhi, m_lastX1, m_lastLogX2;
};

// ---------------------------------------------------------------------------
// Counters
// ---------------------------------------------------------------------------
struct MethodCounters
{
    long long cat[kNumCategories];
    long long nonfinite_output;
    double max_finite_log_component;
    /// Accuracy of the returned radius against the reference, on draws where the comparison
    /// is defined.  Never pooled across precisions: every row that carries these also
    /// carries its `precision` field, as amendment 1.1.0 requires.
    long long n_rel_err;
    double sum_rel_err;
    double max_rel_err;

    MethodCounters()
        : nonfinite_output(0),
          max_finite_log_component(-std::numeric_limits<double>::infinity()), n_rel_err(0),
          sum_rel_err(0.0), max_rel_err(0.0)
    {
        for (int i= 0; i < kNumCategories; ++i)
            cat[i]= 0;
    }

    void observe(int category, bool finite, double rel_err, double ref_max_log_component)
    {
        ++cat[category];
        if (!finite)
            ++nonfinite_output;
        else if (ref_max_log_component > max_finite_log_component)
            max_finite_log_component= ref_max_log_component;
        if (rel_err == rel_err)
        {
            ++n_rel_err;
            sum_rel_err+= rel_err;
            if (rel_err > max_rel_err)
                max_rel_err= rel_err;
        }
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
    long long honest() const
    {
        long long t= 0;
        for (int i= 0; i < kNumCategories; ++i)
            if (isHonest(i))
                t+= cat[i];
        return t;
    }
    /// N = finite + avoidable + honest, with cap rejection and exhaustion excluded because
    /// they are not attempts of an uncapped run.
    bool accountingHolds(long long n) const
    {
        return total() == n && cat[kCatFinite] + avoidable() + honest() == n;
    }
};

struct AttemptCounters
{
    long long n;
    long long x2_zero;
    long long x2_subnormal;
    long long ref_nonrepresentable;
    long long rotation_recoverable;
    long long near_limit;
    long long audited;
    double max_log_r_ref;
    long long stratum_total[kNumAuditStrata];
    long long stratum_audited[kNumAuditStrata];

    AttemptCounters()
        : n(0), x2_zero(0), x2_subnormal(0), ref_nonrepresentable(0), rotation_recoverable(0),
          near_limit(0), audited(0), max_log_r_ref(-std::numeric_limits<double>::infinity())
    {
        for (int i= 0; i < kNumAuditStrata; ++i)
        {
            stratum_total[i]= 0;
            stratum_audited[i]= 0;
        }
    }
};

struct VariateCounts
{
    long long gamma;
    long long uniform;
    long long engine_calls;

    VariateCounts() : gamma(0), uniform(0), engine_calls(0) {}
};

} // namespace exp7

#endif // EXP7_LOADERS_H
