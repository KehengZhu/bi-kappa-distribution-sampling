#include <vector>
#include <random>
#include <cmath>
#include <fstream>
#include <iostream>

// Chombo Headers
#include "REAL.H"
#include "RealVect.H"
#include "CH_Timer.H"
#include "MayDay.H"
#include "parstream.H"

// Local Headers
#include "bi_kappa_distribution.H"

std::vector<RealVect> sampleBiKappa(Real a_kappa,
                                    Real a_thetaPerp,
                                    Real a_thetaPara,
                                    int a_nSamples = 1)
{
    // 1. Validation
    if (a_kappa <= 0.5) {
        MayDay::Error("kappa must be > 0.5 for normalization");
    }
    if (a_thetaPerp <= 0.0 || a_thetaPara <= 0.0) {
        MayDay::Error("Temperature parameters must be positive");
    }

    // 2. Setup Random Number Generation
    std::random_device rd;
    std::mt19937 gen(rd());

    // Gamma distributions for the radius (R = sqrt(X1/X2))
    std::gamma_distribution<Real> distX1(1.5, 1.0);
    std::gamma_distribution<Real> distX2(a_kappa - 0.5, 1.0);

    // Uniform distributions for spherical directions
    std::uniform_real_distribution<Real> distCosTheta(-1.0, 1.0);
    std::uniform_real_distribution<Real> distPhi(0.0, 2.0 * M_PI);

    std::vector<RealVect> samples(a_nSamples);

    for (int i = 0; i < a_nSamples; ++i) {
        // Step 1: Sample radii
        Real X1 = distX1(gen);
        Real X2 = distX2(gen);
        Real T = X1 / X2;
        Real R = std::sqrt(T);

        // Step 2: Sample directions uniformly on S^2
        Real cosTheta = distCosTheta(gen);
        Real sinTheta = std::sqrt(std::max(0.0, 1.0 - cosTheta * cosTheta));
        Real phi = distPhi(gen);

        // Step 3: Compute unit vectors (u)
        Real ux = R * sinTheta * std::cos(phi);
        Real uy = R * sinTheta * std::sin(phi);
        Real uz = R * cosTheta;

        // Step 4: Map back to v-space and store in RealVect
        // RealVect components are accessed by index [0, 1, 2]
        Real sqrtKappa = std::sqrt(a_kappa);

        RealVect v;
        v[0] = sqrtKappa * a_thetaPerp * ux; // Vx
        v[1] = sqrtKappa * a_thetaPerp * uy; // Vy
        v[2] = sqrtKappa * a_thetaPara * uz; // Vz

        samples[i] = v;
    }

    return samples;
}

int main()
{
    // Set parameters
    Real kappa = 2.0;        // Shape parameter (must be > 0.5)
    Real theta_perp = 1.0;   // Perpendicular temperature
    Real theta_para = 2.0;   // Parallel temperature (anisotropic case)
    int n_samples = 100000;  // Number of samples

    // Example 1: Using the original sampleBiKappa function
    pout() << "=== Example 1: Using sampleBiKappa() ===" << std::endl;
    std::vector<RealVect> samples = sampleBiKappa(kappa, theta_perp, theta_para, n_samples);

    // Open output file
    std::ofstream outfile("samples.out");

    // Write samples to file (3 numbers per row)
    for (int i = 0; i < n_samples; ++i) {
        outfile << samples[i][0] << " "
                << samples[i][1] << " "
                << samples[i][2] << "\n";
    }

    outfile.close();
    pout() << "Successfully generated " << n_samples << " samples and wrote to samples.out" << std::endl;

    // Example 2: Using the template class BiKappaDistribution<>
    pout() << "\n=== Example 2: Using BiKappaDistribution<> ===" << std::endl;

    // Setup random number generator
    std::random_device rd;
    std::mt19937 gen(rd());

    // The class returns Point3D (std::array<Real, 3>)
    BiKappaDistribution<Real> distBikappa(kappa, theta_perp, theta_para);

    // Open second output file
    std::ofstream outfile2("samples_class.out");

    // Generate samples one at a time using the distribution class
    for (int i = 0; i < n_samples; ++i) {
        BiKappaDistribution<Real>::Point3D v = distBikappa(gen);
        outfile2 << v[0] << " " << v[1] << " " << v[2] << "\n";
    }

    outfile2.close();
    pout() << "Successfully generated " << n_samples << " samples using BiKappaDistribution<>" << std::endl;
    pout() << "Output written to samples_class.out" << std::endl;

    return 0;
}
