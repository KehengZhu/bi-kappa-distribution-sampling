// LOG-PUB: a published small-shape Gamma generator that returns its samples on the
// logarithmic scale, transcribed from the primary source.
//
// SOURCE
//   Liu C, Martin R, Syring N (2017). "Efficient simulation from a gamma distribution with
//   small shape parameter." Computational Statistics 32(4), 1767-1775,
//   doi:10.1007/s00180-016-0692-0.  Short Note.  Method given in Sect. 3, p. 1770;
//   acceptance rate in Eq. (2), p. 1771; limit result in Theorem 1, p. 1769.
//
// WHY THIS SOURCE RATHER THAN XI, TAN AND LIU (2013)
//   Both were read in full.  Xi, Tan and Liu (2013, J. Stat. Softw. 55(4),
//   doi:10.18637/jss.v055.i04) give a ratio-of-uniforms Gamma generator valid for every
//   alpha > 0 whose Algorithm 1 delivers log-scale output when alpha < 0.01, using
//   tabulated piecewise-linear bounds b_s(theta), b_w(theta) fitted at efficiency R = 0.9.
//   Liu, Martin and Syring (2017) is specific to the small-shape regime this experiment
//   lives in, returns the log scale by design rather than by a hard-coded 0.01 switch, and
//   has a closed-form acceptance rate that gives an exact, falsifiable check on the
//   transcription.  Xi et al. (2013) is recorded as the considered alternative.
//
// WHAT THE PAPER PRINTS, VERBATIM (Sect. 3, p. 1770)
//   For Y ~ Gamma(alpha,1) with alpha near zero, set Z = -alpha log Y.  Its density is
//
//       h_alpha(z) = c exp( -z - exp(-z/alpha) ),   z in (-inf, inf),   c^{-1} = Gamma(alpha+1),
//
//   and the optimal envelope is
//
//       eta_alpha(z) = c exp(-z)                     for z >= 0,
//                      c w lambda exp(lambda z)      for z <  0,
//
//   with  lambda = lambda(alpha) = alpha^{-1} - 1  and  w = w(alpha) = alpha / (e (1-alpha)).
//   The normalized envelope is the two-sided exponential mixture
//
//       (1/(1+w)) Exp(1) + (w/(1+w)) { -Exp(lambda) },
//
//   and the acceptance rate is  r(alpha) = {1 + w(alpha)}^{-1}   (Eq. 2).
//   On the output scale the paper states: "since numerical precision can be lost in
//   exponentiating log Y, we recommend returning the samples on the log-scale, which is
//   the default in our R function."
//
// WHAT THE PAPER DOES NOT PRINT
//   There is no algorithm box, no pseudocode and no code listing in the article; the
//   R implementation (rgamss) is referenced at an external URL and is not reproduced.
//   In particular the acceptance inequality is not printed.  It follows uniquely from the
//   two displayed formulae by the standard acceptance-rejection rule "accept with
//   probability h/eta" that the paper invokes (citing Lange 1999, Chap. 20.4).  Writing
//   m = -z/alpha and cancelling c:
//
//       z >= 0 :  h/eta = exp( -exp(m) )                      -> accept iff log U <= -exp(m)
//       z <  0 :  h/eta = e exp( m - exp(m) )                 -> accept iff log U <= 1 + m - exp(m)
//
//   using  w lambda = [alpha/(e(1-alpha))] [(1-alpha)/alpha] = 1/e  exactly.  Both branches
//   have supremum 1 (the z<0 branch at m = 0, the z>=0 branch as m -> -inf), so the
//   envelope is exactly tangent and the overall acceptance probability is exactly
//   {1+w}^{-1}, reproducing the paper's Eq. (2).  That identity is the transcription check:
//   the measured acceptance rate must match r(alpha), and the pilot tests it.
//
//   Because the inequality is derived rather than quoted, this file is labelled an
//   attributed transcription of a published method, never a quotation of one.
//
// DOMAIN
//   The construction needs w > 0 and lambda > 0, hence 0 < alpha < 1 strictly; at alpha = 1
//   it degenerates (lambda = 0, w undefined).  For alpha >= 1 - which in this experiment
//   means kappa >= 3/2, a numerically benign regime - the draw falls back to
//   std::gamma_distribution<T>(alpha,1) followed by a logarithm, and the draw is flagged
//   out_of_domain so that no result is ever attributed to the published method outside the
//   range its authors state.

#ifndef EXP6_LOG_GAMMA_PUBLISHED_H
#define EXP6_LOG_GAMMA_PUBLISHED_H

#include "log_gamma_identity.H"

#include <cmath>
#include <cstdint>
#include <limits>
#include <random>

namespace exp6
{

template <typename T> class LogGammaPublished
{
  public:
    explicit LogGammaPublished(T alpha)
        : m_alpha(alpha), m_in_domain(alpha > T(0) && alpha < T(1)),
          m_fallback(alpha > T(0) ? alpha : T(1), T(1))
    {
        if (m_in_domain)
        {
            const T e= T(2.71828182845904523536028747135266250L);
            m_lambda= T(1) / alpha - T(1);          // lambda = alpha^{-1} - 1
            m_w= alpha / (e * (T(1) - alpha));      // w = alpha / (e (1-alpha))
            m_ww= T(1) / (T(1) + m_w);              // acceptance rate r(alpha), Eq. (2)
        }
        else
        {
            m_lambda= T(0);
            m_w= T(0);
            m_ww= T(1);
        }
    }

    /// The acceptance rate exactly as published: r(alpha) = {1 + w(alpha)}^{-1}, Eq. (2).
    T publishedAcceptanceRate() const { return m_ww; }

    /// The acceptance rate implied by the paper's own two displays.  For an
    /// acceptance-rejection scheme with target h and envelope eta the acceptance
    /// probability is the ratio of areas, and with the paper's normalizations
    ///
    ///     integral h = 1,        integral eta = c (1 + w) = (1 + w) / Gamma(alpha+1),
    ///
    /// so the rate is Gamma(alpha+1) / {1 + w(alpha)}, not {1 + w(alpha)}^{-1}.
    /// The two agree to leading order as alpha -> 0, which is the regime the paper
    /// addresses - Gamma(alpha+1) = 1 - gamma_E alpha + O(alpha^2) - and the paper's own
    /// prose identifies 1/(1+w) as the probability of *proposing* from the Exp(1) branch,
    /// which is what Eq. (2) appears to reuse.  At alpha = 1/2 the difference is 11%.
    /// This implementation is unaffected: the sampled law is the same either way.  Both
    /// rates are recorded so the pilot can report which one the measurement follows.
    T areaRatioAcceptanceRate() const
    {
        if (!m_in_domain)
            return T(1);
        return static_cast<T>(std::tgamma(static_cast<double>(m_alpha) + 1.0)) * m_ww;
    }

    /// Retained for callers that just want "the" rate; returns the area ratio, which is
    /// what a measurement actually converges to.
    T acceptanceRate() const { return areaRatioAcceptanceRate(); }
    bool inDomain() const { return m_in_domain; }
    T lambda() const { return m_lambda; }
    T w() const { return m_w; }

    template <typename Gen>
    LogGammaDraw<T> draw(Gen &gen, std::uniform_real_distribution<T> &unif)
    {
        LogGammaDraw<T> d;
        d.retries= 0;
        d.attempts= 0;
        d.uniforms= 0;
        d.gammas= 0;
        d.out_of_domain= !m_in_domain;
        d.y= T(0);
        d.u= T(0);

        if (!m_in_domain)
        {
            // Outside the published domain.  Benign regime; a plain Gamma draw and its
            // logarithm, flagged so the analysis never attributes it to the method.
            const T x= m_fallback(gen);
            d.gammas= 1;
            d.attempts= 1;
            d.x_product= x;
            d.log_x= std::log(x);
            d.ok= isFiniteT(d.log_x);
            return d;
        }

        const T log_max= static_cast<T>(std::log(
            static_cast<double>(std::numeric_limits<T>::max())));

        for (;;)
        {
            ++d.attempts;

            // Inverse-transform draw from the normalized envelope mixture
            //     (1/(1+w)) Exp(1) + (w/(1+w)) { -Exp(lambda) }.
            T u1= unif(gen);
            ++d.uniforms;
            while (!(u1 > T(0)))
            {
                ++d.retries;
                ++d.uniforms;
                u1= unif(gen);
            }

            T z;
            if (u1 <= m_ww)
                z= -std::log(u1 / m_ww);                       // z >= 0, Exp(1) branch
            else
                z= std::log((u1 - m_ww) / (T(1) - m_ww)) / m_lambda; // z < 0, -Exp(lambda)

            // m = -z/alpha is log Y for the proposed draw, and is what is returned on
            // acceptance.  It is never exponentiated on the accept path.
            const T m= -z / m_alpha;

            T log_ratio;
            if (m > log_max)
            {
                log_ratio= -std::numeric_limits<T>::infinity(); // exp(m) overflows: reject
            }
            else
            {
                const T em= std::exp(m);
                log_ratio= (z < T(0)) ? (T(1) + m - em) : (-em);
            }

            T u2= unif(gen);
            ++d.uniforms;
            while (!(u2 > T(0)))
            {
                ++d.retries;
                ++d.uniforms;
                u2= unif(gen);
            }

            if (std::log(u2) <= log_ratio)
            {
                d.log_x= m;
                d.u= u1;
                // The natural-scale value the caller would have materialized.  Recorded
                // only so that the underflow it suffers can be counted; never used.
                d.x_product= (m > log_max) ? std::numeric_limits<T>::infinity() : std::exp(m);
                d.ok= isFiniteT(d.log_x);
                return d;
            }
        }
    }

  private:
    static bool isFiniteT(T v)
    {
        return (v == v) && v != std::numeric_limits<T>::infinity() &&
               v != -std::numeric_limits<T>::infinity();
    }

    T m_alpha;
    bool m_in_domain;
    T m_lambda;
    T m_w;
    T m_ww;
    std::gamma_distribution<T> m_fallback;
};

} // namespace exp6

#endif // EXP6_LOG_GAMMA_PUBLISHED_H
