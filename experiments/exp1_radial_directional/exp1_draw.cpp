// Experiment 1 large-sample cell test: one seed of the uncapped sampler at one kappa,
// written raw.
//
// exp1_cells.py drives it one seed at a time and keeps only the cell counts; the raw
// draws are deleted as soon as they are counted.
//
// Usage: exp1_draw.exe kappa seed n_samples out.bin [theta_perp theta_par ub_x ub_y ub_z]
//        Default (theta_perp, theta_par) = (1, 2), field along z; uncapped, double
//        precision.  The optional arguments give the thermal speeds and the field
//        direction for the frame-invariance test (exp1_frame.py); ub is passed as given.

#include <array>
#include <cstdio>
#include <cstdlib>
#include <vector>

typedef double Real;

#include "../../cpp/bi_kappa_distribution.H"

int main(int argc, char **argv)
{
    if (argc != 5 && argc != 10)
    {
        std::fprintf(stderr,
                     "usage: %s kappa seed n_samples out.bin [theta_perp theta_par ub_x ub_y ub_z]\n",
                     argv[0]);
        return 2;
    }
    const Real kappa= std::atof(argv[1]);
    const int seed= std::atoi(argv[2]);
    const long n= std::atol(argv[3]);
    Real thetaPerp= 1.0, thetaPar= 2.0;
    std::array<Real, 3> ub= {0.0, 0.0, 1.0};
    if (argc == 10)
    {
        thetaPerp= std::atof(argv[5]);
        thetaPar= std::atof(argv[6]);
        ub= {std::atof(argv[7]), std::atof(argv[8]), std::atof(argv[9])};
    }

    bi_kappa_distribution<Real> dist;
    dist.define(kappa, thetaPerp, thetaPar, ub, bi_kappa_distribution<Real>::no_cap(), seed);

    std::vector<double> buffer(static_cast<size_t>(3 * n));
    for (long i= 0; i < n; ++i)
    {
        const bi_kappa_distribution<Real>::point_type v= dist();
        buffer[3 * i]= v[0];
        buffer[3 * i + 1]= v[1];
        buffer[3 * i + 2]= v[2];
    }

    std::FILE *f= std::fopen(argv[4], "wb");
    if (!f)
    {
        std::fprintf(stderr, "cannot write %s\n", argv[4]);
        return 1;
    }
    std::fwrite(buffer.data(), sizeof(double), buffer.size(), f);
    std::fclose(f);
    return 0;
}
