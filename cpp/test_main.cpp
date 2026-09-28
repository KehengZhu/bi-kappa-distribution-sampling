// Entry point of the regression suite (test_suite.H).
//
// It runs the suite and nothing else: it writes no files, and its exit status is 0 if
// every check passes and 1 otherwise.  Build and run it with `make test`.  The demo
// program main.cpp, built by `make`, is separate: it runs no tests and writes five
// sample files of 50,000 lines each (about 6 MB).

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
