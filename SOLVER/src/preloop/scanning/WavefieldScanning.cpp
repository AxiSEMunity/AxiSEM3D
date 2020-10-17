//
//  WavefieldScanning.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 5/28/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  wavefield scanning

#include "WavefieldScanning.hpp"

#include "SE_Model.hpp"
#include "Domain.hpp"

#include "inparam.hpp"
#include "timer.hpp"
#include "io.hpp"
#include "bstring.hpp"

// setup
void WavefieldScanning::setup(double dt, double period, int numTotalSteps,
                              const SE_Model &sem, Domain &domain) {
    // create
    std::unique_ptr<WavefieldScanning> ws =
    std::make_unique<WavefieldScanning>();
    
    // inparam
    timer::gPreloopTimer.begin("Building from inparam");
    const InparamYAML &gm = inparam::gInparamNr;
    const std::string &rt = "wavefield_scanning";
    // file
    ws->mFileName = gm.get<std::string>(rt + ":output_file");
    // tolerance
    double tolFourierH1 =
    gm.getWithBounds(rt + ":threshold_Fourier_convergence", 1e-4, 1e-1);
    double relTolH1 =
    gm.getWithBounds(rt + ":relative_amplitude_skipped", 0., 1.);
    double absTolH1 =
    gm.getWithBounds(rt + ":advanced:absolute_amplitude_skipped", 1e-14, 1e-10);
    ws->mTolFourierH2 = tolFourierH1 * tolFourierH1;
    ws->mRelTolH2 = relTolH1 * relTolH1;
    ws->mAbsTolH2 = absTolH1 * absTolH1;
    ws->mMaxNumPeaks = gm.getWithBounds(rt + ":advanced:max_num_peaks", 1);
    // vertex
    bool vertexOnly = gm.get<bool>(rt + ":advanced:vertex_only");
    timer::gPreloopTimer.ended("Building from inparam");
    // time step
    int nStepsPerPeriod =
    gm.getWithBounds(rt + ":advanced:num_steps_per_mesh_period", 4);
    ws->mScanningInterval =
    std::max((int)round(period / nStepsPerPeriod / dt), 1);
    
    // prepare points for scanning
    timer::gPreloopTimer.begin("Initializing scanning on GLL points");
    sem.initScanningOnPoints(vertexOnly);
    timer::gPreloopTimer.ended("Initializing scanning on GLL points");
    
    // verbose
    if (io::gVerbose != io::VerboseLevel::None) {
        std::stringstream ss;
        using namespace bstring;
        ss << boxTitle("Wavefield Scanning");
        ss << boxEquals(0, 40, "output file for scanning result",
                        ws->mFileName);
        ss << boxEquals(0, 40, "threshold for Fourier series convergence",
                        tolFourierH1);
        ss << boxEquals(0, 40, "relative amplitude skipped for scanning",
                        relTolH1);
        ss << boxEquals(0, 40, "absolute amplitude skipped for scanning",
                        absTolH1);
        ss << boxEquals(0, 40, "maximum number of energy peaks",
                        ws->mMaxNumPeaks);
        ss << boxEquals(0, 40, "perform scanning only on vertex points",
                        vertexOnly);
        ss << boxEquals(0, 40, "# time steps scanned per mesh period",
                        nStepsPerPeriod);
        ss << boxEquals(0, 40, "# time steps between two scanning steps",
                        ws->mScanningInterval);
        ss << boxEquals(0, 40, "# time steps scanned in total",
                        numTotalSteps / ws->mScanningInterval +
                        (numTotalSteps % ws->mScanningInterval > 0));
        ss << boxBaseline() << "\n\n";
        io::cout << ss.str();
    }
    
    // release scanning to domain
    domain.setWavefieldScanning(ws);
}
