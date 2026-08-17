// Experiment 2 sample generator -- characterization of the component-wise
// velocity cap.
//
// Two families of runs are emitted from the *released* sampler:
//
//   mode = "uncapped"  max_normalized_velocity = no_cap().  The sampled law is
//                      the untruncated bi-Kappa distribution.  These runs are
//                      the reference target AND the instrument used to measure
//                      the rejected fraction: the analysis evaluates the box
//                      predicate on each uncapped draw, which gives P(accept)
//                      directly without touching cpp/bi_kappa_distribution.H.
//
//   mode = "capped"    max_normalized_velocity = lambda (finite).  The sampled
//                      law is the bi-Kappa distribution CONDITIONED on the
//                      component-wise box event
//                          |v_x|/theta_perp <= lambda AND
//                          |v_y|/theta_perp <= lambda AND
//                          |v_z|/theta_par  <= lambda.
//                      It is a different probability law, not an approximation
//                      of the same one.
//
// A capped run and the uncapped run at the same seed share the RNG stream: each
// loop iteration of operator() consumes x1, x2, cosTheta, phi in the same order
// regardless of the cap, so the capped output is *bitwise* the subsequence of
// the uncapped output that falls inside the box.  exp2_analyze.py verifies this,
// which turns "capped = uncapped conditioned on the box" from an assertion about
// the code into a measured fact about the shipped binary.
//
// Every run is fully specified by the manifest this writes.  Nothing in the
// analysis depends on state that is not recorded there.
//
// Build:  make
// Run:    ./exp2_sample.exe <output_dir>

#include <array>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

typedef double Real;

#include "../../cpp/bi_kappa_distribution.H"

namespace
{

struct Case
{
    std::string block; // "A" anisotropic (1,2) ladder, "C" isotropic (1,1) control
    Real kappa;
    Real theta_perp;
    Real theta_par;
};

// kappa ladder.  0.75 and 1.0 are the heavy-tailed cases where the untruncated
// second moment does not exist (that needs kappa > 3/2); 1.5 is the boundary
// case where it diverges; 2, 5, 10 are the manuscript's working values.
const Real kKappa[]= {0.75, 1.0, 1.5, 2.0, 5.0, 10.0};

// Isotropic control.  The cap is applied in *normalized* components, so
// |v_x|/theta_perp = sqrt(kappa)|u_x| etc. and the box event is a cube of
// half-side lambda/sqrt(kappa) in the isotropic u-coordinates -- independent of
// theta_perp and theta_par.  The control exists to verify that prediction
// rather than to assume it.
const Real kControlKappa[]= {0.75, 2.0};

// lambda ladder (max_normalized_velocity).  no_cap() is emitted separately as the
// reference mode.  20 is the library default and 100 is the value used in the
// manuscript's example, so both must stay visible in the sweep.
const Real kLambda[]= {3.0, 5.0, 10.0, 20.0, 50.0, 100.0};

const int kSeeds[]= {2001, 2002, 2003, 2004, 2005};
const int kNSamples= 100000;

// Field direction is fixed to +z throughout.  The cap is tested *before* the
// field-frame rotation (bi_kappa_distribution.H step 4 vs step 5), so the box
// event -- and therefore the acceptance probability and the entire capped law in
// the local frame -- does not depend on ub.  Experiment 1 already validates the
// rotation itself.
const std::array<Real, 3> kUb= {0.0, 0.0, 1.0};

std::string fmt(Real v)
{
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%.17g", static_cast<double>(v));
    return std::string(buf);
}

} // namespace

int main(int argc, char **argv)
{
    const std::string outDir= (argc > 1) ? argv[1] : "raw";

    std::vector<Case> cases;
    for (Real kappa : kKappa)
    {
        cases.push_back({"A", kappa, 1.0, 2.0});
    }
    for (Real kappa : kControlKappa)
    {
        cases.push_back({"C", kappa, 1.0, 1.0});
    }

    std::ofstream manifest(outDir + "/manifest.csv");
    if (!manifest)
    {
        std::cerr << "cannot write manifest in " << outDir << " (create the directory first)"
                  << std::endl;
        return 1;
    }
    manifest << "run_id,block,kappa,theta_perp,theta_par,ub_x,ub_y,ub_z,mode,lambda,seed,"
                "n_samples,n_nonfinite,file\n";

    int runId= 0;
    for (const Case &c : cases)
    {
        // Reference mode first, then the lambda ladder, so that a given
        // (case, seed) pair keeps its uncapped reference adjacent in the manifest.
        std::vector<Real> caps;
        caps.push_back(bi_kappa_distribution<Real>::no_cap());
        for (Real lam : kLambda)
        {
            caps.push_back(lam);
        }

        for (Real cap : caps)
        {
            const bool capped= std::isfinite(cap);
            for (int seed : kSeeds)
            {
                char nameBuf[256];
                std::snprintf(nameBuf, sizeof(nameBuf), "run_%04d.bin", runId);
                const std::string fileName= nameBuf;

                bi_kappa_distribution<Real> dist;
                dist.define(c.kappa, c.theta_perp, c.theta_par, kUb, cap, seed);

                if (dist.param().capped() != capped)
                {
                    std::cerr << "run " << runId << ": capped() disagrees with the requested mode"
                              << std::endl;
                    return 1;
                }

                std::vector<Real> buffer;
                buffer.reserve(static_cast<size_t>(kNSamples) * 3);

                long long nonFinite= 0;
                for (int i= 0; i < kNSamples; ++i)
                {
                    const bi_kappa_distribution<Real>::point_type v= dist();
                    if (!std::isfinite(v[0]) || !std::isfinite(v[1]) || !std::isfinite(v[2]))
                    {
                        ++nonFinite;
                    }
                    buffer.push_back(v[0]);
                    buffer.push_back(v[1]);
                    buffer.push_back(v[2]);
                }

                std::ofstream out(outDir + "/" + fileName, std::ios::binary);
                out.write(reinterpret_cast<const char *>(buffer.data()),
                          static_cast<std::streamsize>(buffer.size() * sizeof(Real)));
                out.close();

                manifest << runId << ',' << c.block << ',' << fmt(c.kappa) << ','
                         << fmt(c.theta_perp) << ',' << fmt(c.theta_par) << ',' << fmt(kUb[0])
                         << ',' << fmt(kUb[1]) << ',' << fmt(kUb[2]) << ','
                         << (capped ? "capped" : "uncapped") << ','
                         << (capped ? fmt(cap) : std::string("inf")) << ',' << seed << ','
                         << kNSamples << ',' << nonFinite << ',' << fileName << '\n';

                ++runId;
            }
        }
    }

    manifest.close();
    std::cout << "wrote " << runId << " runs (" << kNSamples << " samples each) to " << outDir
              << std::endl;
    return 0;
}
