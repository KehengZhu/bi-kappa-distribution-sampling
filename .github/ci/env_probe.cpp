// Experiment 7, P5 portability - the floating-point environment the run executes in.
//
// PROTOCOL.md 5.3 decides G4 on two claims: bitwise equality of the candidate's stream
// across standard libraries on one architecture, and equality of rates across
// architectures.  Neither claim can be read without knowing the arithmetic the numbers
// were produced by, so every environment records it rather than assuming it.  What is
// recorded here is measured, not declared: subnormal flushing (FTZ/DAZ on x86, the FPCR
// FZ bit on AArch64) turns a small denominator into an exact zero earlier than IEEE
// arithmetic would, and that is the mechanism the experiment is about.
//
// Strict C++11, standard library only.  Build with the same flags as the probe:
//
//     $CXX -std=c++11 -O2 -ffp-contract=off -I cpp -o env_probe .github/ci/env_probe.cpp
//
// Writes one JSON object to stdout.  It is merged into the environment record by
// .github/ci/make_env_record.py, which adds what only the shell can see (target triple,
// compiler version string, host identity, whether execution is native or translated).

#include <cfenv>
#include <cfloat>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <limits>
#include <random>

#include "bi_kappa_distribution.H"

#if defined(__x86_64__) || defined(__i386__)
#include <xmmintrin.h>
#endif

namespace
{

const char *compilerId()
{
#if defined(__clang__)
    return "clang " __clang_version__;
#elif defined(__GNUC__)
    return "gcc " __VERSION__;
#else
    return "unknown";
#endif
}

const char *compilerFamily()
{
#if defined(__clang__)
    return "clang";
#elif defined(__GNUC__)
    return "gcc";
#else
    return "unknown";
#endif
}

const char *stdlibId()
{
#if defined(_LIBCPP_VERSION)
    return "libc++";
#elif defined(__GLIBCXX__)
    return "libstdc++";
#else
    return "unknown";
#endif
}

/// The standard library's own version stamp: `_LIBCPP_VERSION` for libc++,
/// `__GLIBCXX__` (a release date) for libstdc++.  Reported as an integer so that two
/// records can be compared without parsing prose.
long stdlibVersion()
{
#if defined(_LIBCPP_VERSION)
    return static_cast<long>(_LIBCPP_VERSION);
#elif defined(__GLIBCXX__)
    return static_cast<long>(__GLIBCXX__);
#else
    return -1;
#endif
}

long stdlibRelease()
{
#if defined(_GLIBCXX_RELEASE)
    return static_cast<long>(_GLIBCXX_RELEASE);
#else
    return -1;
#endif
}

/// The architecture the *binary* was compiled for.  make_env_record.py compares it with
/// the architecture the host reports; a mismatch is translated execution, which
/// PROTOCOL.md 5.3 excludes from the F7 decision.
const char *archId()
{
#if defined(__aarch64__) || defined(_M_ARM64)
    return "arm64";
#elif defined(__x86_64__)
    return "x86_64";
#else
    return "unknown";
#endif
}

/// A real x rounds to infinity in type T exactly when |x| >= (2 - 2^-p) 2^emax, which is
/// half an ulp above the largest finite value.  That threshold, not max(), is what decides
/// honest overflow, so it is the number the environment record has to carry.
template <typename T> double logOverflowThreshold()
{
    const int p= std::numeric_limits<T>::digits;
    const int emax= std::numeric_limits<T>::max_exponent - 1;
    const double kLn2= 0.69314718055994530942;
    return static_cast<double>(emax) * kLn2 + std::log(2.0 - std::ldexp(1.0, -p));
}

template <typename T> double logMaxFinite()
{
    return std::log(static_cast<double>(std::numeric_limits<T>::max()));
}

/// %.17g round-trips a double exactly; anything less silently loses the last digits of a
/// threshold that the analysis then compares to more places than it actually has.
void emitDouble(const char *key, double v)
{
    std::printf("\"%s\":%.17g", key, v);
}

const char *roundingMode()
{
    switch (std::fegetround())
    {
    case FE_TONEAREST: return "FE_TONEAREST";
    case FE_DOWNWARD: return "FE_DOWNWARD";
    case FE_UPWARD: return "FE_UPWARD";
    case FE_TOWARDZERO: return "FE_TOWARDZERO";
    default: return "unknown";
    }
}

} // namespace

int main()
{
    // Measured, not declared: multiply the smallest subnormal by one and see whether it
    // survives.  This is the operative fact whatever the control register says.
    volatile double dsub= std::numeric_limits<double>::denorm_min();
    volatile float fsub= std::numeric_limits<float>::denorm_min();
    volatile double dprod= dsub * 1.0;
    volatile float fprod= fsub * 1.0f;
    const bool d_flushes= !(dprod == dsub) || dprod == 0.0;
    const bool f_flushes= !(fprod == fsub) || fprod == 0.0f;

    std::printf("{");
    std::printf("\"compiler\":\"%s\",\"compiler_family\":\"%s\",", compilerId(),
                compilerFamily());
    std::printf("\"stdlib\":\"%s\",\"stdlib_version\":%ld,\"stdlib_release\":%ld,", stdlibId(),
                stdlibVersion(), stdlibRelease());
    std::printf("\"arch_compiled_for\":\"%s\",\"cplusplus\":%ld,", archId(),
                static_cast<long>(__cplusplus));
    std::printf("\"bi_kappa_version\":\"%s\",", BI_KAPPA_VERSION_STRING);

    std::printf("\"rounding_mode\":\"%s\",\"flt_eval_method\":%d,", roundingMode(),
                static_cast<int>(FLT_EVAL_METHOD));
    std::printf("\"double_subnormals_flushed\":%s,\"float_subnormals_flushed\":%s,",
                d_flushes ? "true" : "false", f_flushes ? "true" : "false");

    // The control register itself where it is readable, so that "no flushing observed" can
    // be distinguished from "flushing is off".
#if defined(__x86_64__) || defined(__i386__)
    {
        const unsigned int csr= _mm_getcsr();
        std::printf("\"ftz_daz\":{\"source\":\"mxcsr\",\"raw\":\"0x%08x\",\"ftz\":%s,"
                    "\"daz\":%s},",
                    csr, (csr & 0x8000u) ? "true" : "false", (csr & 0x0040u) ? "true" : "false");
    }
#elif defined(__aarch64__)
    {
        unsigned long long fpcr= 0;
        __asm__ __volatile__("mrs %0, fpcr" : "=r"(fpcr));
        std::printf("\"ftz_daz\":{\"source\":\"fpcr\",\"raw\":\"0x%016llx\",\"fz\":%s,"
                    "\"fz16\":%s},",
                    fpcr, (fpcr & (1ULL << 24)) ? "true" : "false",
                    (fpcr & (1ULL << 19)) ? "true" : "false");
    }
#else
    std::printf("\"ftz_daz\":{\"source\":\"measured only\"},");
#endif

#if defined(FP_FAST_FMA)
    std::printf("\"fp_fast_fma\":true,");
#else
    std::printf("\"fp_fast_fma\":false,");
#endif

    std::printf("\"flt_radix\":%d,\"flt_mant_dig\":%d,\"dbl_mant_dig\":%d,"
                "\"ldbl_mant_dig\":%d,\"sizeof_long_double\":%zu,",
                static_cast<int>(FLT_RADIX), static_cast<int>(FLT_MANT_DIG),
                static_cast<int>(DBL_MANT_DIG), static_cast<int>(LDBL_MANT_DIG),
                sizeof(long double));

    std::printf("\"dbl_max_hex\":\"%a\",\"flt_max_hex\":\"%a\",",
                static_cast<double>(DBL_MAX), static_cast<double>(FLT_MAX));

    emitDouble("log_overflow_threshold_double", logOverflowThreshold<double>());
    std::printf(",");
    emitDouble("log_overflow_threshold_float", logOverflowThreshold<float>());
    std::printf(",");
    emitDouble("log_max_double", logMaxFinite<double>());
    std::printf(",");
    emitDouble("log_max_float", logMaxFinite<float>());
    std::printf(",");

    std::printf("\"rng\":\"mt19937\",\"mt19937_min\":%llu,\"mt19937_max\":%llu",
                static_cast<unsigned long long>(std::mt19937::min()),
                static_cast<unsigned long long>(std::mt19937::max()));
    std::printf("}\n");
    return 0;
}
