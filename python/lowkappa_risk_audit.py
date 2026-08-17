"""Low-kappa numerical-risk audit for the bi-Kappa Gamma-ratio radial draw.

Question under test
-------------------
Does the bounded transform Y = T/(1+T) ~ Beta(3/2, kappa-1/2) remove the
zero-division failure mode that Zenitani et al. (2026) flag for
shape(X2) = kappa - 1/2 < 1, or does it only relocate it?

Answer: it relocates it AND makes it far worse. See main() output.

Radial construction under test (cpp/bi_kappa_distribution.H:239):
    X1 ~ Ga(3/2, 1),  X2 ~ Ga(kappa - 1/2, 1),  R = sqrt(X1 / X2)

Variants compared
    naive  : R = sqrt(X1 / X2)              <- what the C++ does today
    split  : R = sqrt(X1) / sqrt(X2)        <- avoids the intermediate T overflow
    logdom : R = exp(0.5*(log X1 - log X2)) <- log of the *already drawn* X2
    beta   : Y = X1/(X1+X2), T = Y/(1-Y), R = sqrt(T)  <- the proposed "fix"
    logga  : log-space gamma, log X2 never materialised as a float

Run:  python3 python/lowkappa_risk_audit.py
"""

import numpy as np
from scipy import stats

TINY = np.finfo(float).tiny  # smallest normal,   2.225e-308
EPS = np.finfo(float).eps  # 2.220e-16
LOG_MAX = np.log(np.finfo(float).max)  # 709.78


def draw(kappa, n, rng):
    x1 = rng.standard_gamma(1.5, n)
    x2 = rng.standard_gamma(kappa - 0.5, n)
    return x1, x2


def log_gamma_small_shape(b, n, rng):
    """log of a Ga(b,1) variate for b < 1, computed without underflow.

    numpy and libstdc++ both realise a small-shape gamma via the Marsaglia-Tsang
    boost  X = Ga(b+1) * U**(1/b).  Taking the log of that identity,
        log X = log Ga(b+1) + (1/b) * log U
    is exactly computable: log U is O(-1e2) and (1/b)*log U is O(-1e5), both far
    inside double range, even where X itself flushes to zero.
    """
    g = rng.standard_gamma(b + 1.0, n)
    u = rng.uniform(size=n)
    return np.log(g) + np.log(u) / b


def variants(x1, x2):
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        naive = np.sqrt(x1 / x2)
        split = np.sqrt(x1) / np.sqrt(x2)
        logdom = np.exp(0.5 * (np.log(x1) - np.log(x2)))
        y = x1 / (x1 + x2)
        beta = np.sqrt(y / (1.0 - y))
    return {"naive": naive, "split": split, "logdom": logdom, "beta": beta}


def bad(r):
    return int(np.count_nonzero(~np.isfinite(r)))


def main():
    rng = np.random.default_rng(20260815)
    n = 2_000_000
    kappas = [1.5, 1.0, 0.75, 0.6, 0.55, 0.51, 0.505, 0.501, 0.5001]

    print(f"N = {n:,} per kappa;  smallest normal = {TINY:.3e};  eps = {EPS:.3e}")

    # ------------------------------------------------------------------
    print("\n[1] Failure counts by formulation (non-finite radii)\n")
    print(
        f"{'kappa':>8} {'b=k-1/2':>9} {'X2==0':>10} {'naive':>10} {'split':>10} "
        f"{'logdom':>10} {'beta':>10} {'pred P(X2=0)':>13}"
    )
    print("-" * 94)
    for kappa in kappas:
        b = kappa - 0.5
        x1, x2 = draw(kappa, n, rng)
        v = variants(x1, x2)
        # X2 ~ U**(1/b) for small b  =>  P(X2 == 0) ~ 10**(-323.3 b)
        print(
            f"{kappa:>8.4f} {b:>9.4f} {int((x2 == 0).sum()):>10,} "
            f"{bad(v['naive']):>10,} {bad(v['split']):>10,} "
            f"{bad(v['logdom']):>10,} {bad(v['beta']):>10,} "
            f"{10.0 ** (-323.3 * b):>13.2e}"
        )

    # ------------------------------------------------------------------
    print("\n[2] Why 'beta' fails first: Y rounds to 1.0 long before X2 underflows")
    print(f"    Y = X1/(X1+X2) == 1.0 as soon as X2/X1 < eps/2 = {EPS/2:.3e},")
    print(f"    whereas X2 itself survives down to {5e-324:.1e}.")
    print("    The bounded variable therefore truncates the usable dynamic range")
    print(f"    of T from ~1e308 to ~1/eps = {1/EPS:.3e}: a loss of ~292 decades.")
    x1, x2 = draw(0.55, n, rng)
    y = x1 / (x1 + x2)
    print(f"    kappa=0.55: Y==1.0 exactly for {int((y == 1.0).sum()):,}/{n:,} draws, "
          f"while X2==0 for {int((x2 == 0).sum()):,}.")

    # ------------------------------------------------------------------
    print("\n[3] 'naive' vs 'split': intermediate overflow of T = X1/X2, kappa=0.51")
    x1, x2 = draw(0.51, n, rng)
    v = variants(x1, x2)
    lost = np.isfinite(v["split"]) & ~np.isfinite(v["naive"])
    print(f"    finite under split but Inf/NaN under naive: {int(lost.sum()):,}")
    if lost.any():
        r = v["split"][lost]
        print(f"    those radii span [{r.min():.3e}, {r.max():.3e}] - all representable")
        print("    T = R^2 overflows for every one of them; R does not.")
        print("    => sqrt(X1)/sqrt(X2) is a strict, free improvement over sqrt(X1/X2).")

    # ------------------------------------------------------------------
    print("\n[4] 'logdom' buys nothing once X2 has already been drawn")
    print("    log(0.0) = -inf  ->  exp(+inf) = inf.  Identical failure count to split.")
    print("    A log-domain route only helps if log X2 is produced by the generator.")

    # ------------------------------------------------------------------
    print("\n[5] Log-space gamma: how much of the failure is recoverable?\n")
    print("    Splitting failures into")
    print("      spurious - true R is representable, lost only to underflow of X2")
    print("      honest   - true R exceeds the double range; no reformulation helps\n")
    print(f"{'kappa':>8} {'log-space fail':>15} {'spurious':>10} {'honest':>10} "
          f"{'naive fail':>11}")
    print("-" * 60)
    for kappa in [0.51, 0.505, 0.501, 0.5001]:
        b = kappa - 0.5
        log_x1 = np.log(rng.standard_gamma(1.5, n))
        log_x2 = log_gamma_small_shape(b, n, rng)
        log_r = 0.5 * (log_x1 - log_x2)
        honest = int((log_r > LOG_MAX).sum())  # genuinely unrepresentable
        # what the float pipeline would have lost at the same parameters
        _, x2 = draw(kappa, n, rng)
        x1b = rng.standard_gamma(1.5, n)
        naive_fail = bad(variants(x1b, x2)["naive"])
        print(f"{kappa:>8.4f} {int((~np.isfinite(log_r)).sum()):>15,} "
              f"{max(naive_fail - honest, 0):>10,} {honest:>10,} {naive_fail:>11,}")

    # ------------------------------------------------------------------
    print("\n[6] Does the log-space gamma sample the right law? (KS on log X2)")
    for b in [0.5, 0.05, 0.01]:
        m = 200_000
        ref = rng.standard_gamma(b, m)
        ref = np.log(ref[ref > 0])
        alt = log_gamma_small_shape(b, m, rng)
        alt = alt[np.isfinite(alt)]
        lo = max(ref.min(), alt.min())
        ks = stats.ks_2samp(ref[ref >= lo], alt[alt >= lo])
        print(f"    b={b:<6} KS={ks.statistic:.5f}  p={ks.pvalue:.3f}  "
              f"(compared on the range numpy can represent)")


if __name__ == "__main__":
    main()
