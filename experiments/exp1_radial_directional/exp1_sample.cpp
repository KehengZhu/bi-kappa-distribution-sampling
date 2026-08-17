// Experiment 1 sample generator.
//
// Emits raw samples from the *released* uncapped bi_kappa_distribution so that
// the radial law, directional uniformity, radial-direction independence and
// frame invariance can be checked directly (reviewer comment R1.3).
//
// Every run is fully specified by the manifest this writes: kappa, theta_perp,
// theta_par, the magnetic-field direction, the sampling mode, the seed and N.
// Nothing about the analysis depends on state that is not recorded there.
//
// Build:  make -f GNUmakefile.exp1
// Run:    ./exp1_sample.exe <output_dir>

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
    std::string block;    // "A" (radial core) or "B" (anisotropy x frame)
    Real kappa;
    Real theta_perp;
    Real theta_par;
    std::array<Real, 3> ub;
    std::string ub_label;
};

// Block A -- the radial/directional core.  Covers the low-kappa range R1.3 asks
// about (1/2 < kappa <= 3/2) as well as the three kappa values already in the
// manuscript.  Anisotropic, field-aligned.
const Real kBlockAKappa[]= {0.51, 0.55, 0.75, 1.0, 1.25, 1.5, 2.0, 5.0, 10.0};

// Block B -- anisotropy and arbitrary field frame.  Three kappa values spanning
// the range, both isotropic and 1:2 scales, one axis-aligned and two
// non-axis-aligned field directions.
const Real kBlockBKappa[]= {0.55, 2.0, 10.0};

const int kSeeds[]= {1001, 1002, 1003, 1004, 1005};
const int kNSamples= 100000;

std::array<Real, 3> normalized(Real x, Real y, Real z)
{
    Real n= std::sqrt(x * x + y * y + z * z);
    return {x / n, y / n, z / n};
}

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

    // ---- Block A ----
    for (Real kappa : kBlockAKappa)
    {
        cases.push_back({"A", kappa, 1.0, 2.0, {0.0, 0.0, 1.0}, "z"});
    }

    // ---- Block B ----
    const std::array<Real, 3> bHat1= {0.0, 0.0, 1.0};
    const std::array<Real, 3> bHat2= normalized(1.0, 1.0, 1.0);
    const std::array<Real, 3> bHat3= normalized(0.3, -0.5, 0.8);
    for (Real kappa : kBlockBKappa)
    {
        for (int aniso= 0; aniso < 2; ++aniso)
        {
            const Real tPerp= 1.0;
            const Real tPar= (aniso == 0) ? 1.0 : 2.0;
            cases.push_back({"B", kappa, tPerp, tPar, bHat1, "z"});
            cases.push_back({"B", kappa, tPerp, tPar, bHat2, "diag111"});
            cases.push_back({"B", kappa, tPerp, tPar, bHat3, "oblique"});
        }
    }

    std::ofstream manifest(outDir + "/manifest.csv");
    if (!manifest)
    {
        std::cerr << "cannot write manifest in " << outDir
                  << " (create the directory first)" << std::endl;
        return 1;
    }
    manifest << "run_id,block,kappa,theta_perp,theta_par,ub_x,ub_y,ub_z,ub_label,"
                "mode,seed,n_samples,n_nonfinite,file\n";

    int runId= 0;
    for (const Case &c : cases)
    {
        for (int seed : kSeeds)
        {
            char nameBuf[256];
            std::snprintf(nameBuf, sizeof(nameBuf), "run_%04d.bin", runId);
            const std::string fileName= nameBuf;

            // Uncapped mode: no_cap() disables the component-wise box entirely,
            // so the sampled law is the full bi-Kappa distribution and the core
            // mapping runs exactly once per draw.
            bi_kappa_distribution<Real> dist;
            dist.define(c.kappa, c.theta_perp, c.theta_par, c.ub,
                        bi_kappa_distribution<Real>::no_cap(), seed);

            if (dist.param().capped())
            {
                std::cerr << "run " << runId << " is capped; expected uncapped" << std::endl;
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
                     << fmt(c.theta_perp) << ',' << fmt(c.theta_par) << ',' << fmt(c.ub[0])
                     << ',' << fmt(c.ub[1]) << ',' << fmt(c.ub[2]) << ',' << c.ub_label
                     << ",uncapped," << seed << ',' << kNSamples << ',' << nonFinite << ','
                     << fileName << '\n';

            ++runId;
        }
    }

    manifest.close();
    std::cout << "wrote " << runId << " runs (" << kNSamples << " samples each) to " << outDir
              << std::endl;
    return 0;
}
