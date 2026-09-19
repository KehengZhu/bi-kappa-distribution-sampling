// Test-only entry point.
//
// main.cpp runs the suite and then writes about 25 MB of sample files, which makes
// it unusable as a continuous-integration step: the samples are examples, not
// assertions, and nothing reads them back.  This translation unit runs the suite
// and nothing else, so `make test` is a pure pass/fail with no side effects on the
// working directory.  Build it with `make test`; `make` still builds main.exe.

#include <cmath>
#include <iostream>
#include <random>

typedef double Real;

// Local headers
#include "bi_kappa_distribution.H"
#include "bi_maxwellian_distribution.H"
#include "general_position_generator.H"
#include "general_velocity_generator.H"
#include "field_aligned_velocity_generator.H"
#include "test_suite.H"

int main()
{
    return run_all_tests();
}
