//
//  SolverFFTW.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/22/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  FFT solver based on FFTW3
//  v2.0 features:
//  * all logical sizes share the same static memory
//  * external data feeding and retrieving are enforced for readability
//  * template-based, available for both core and preloop

#include "SolverFFTW.hpp"

// create plans after adding ALL NRs
template <typename floatT, int HOWMANY>
void SolverFFTW<floatT, HOWMANY>::createPlans(double timeLimitForPlanning) {
    // clear if exists
    clearPlans();
    
    // empty
    if (mToBeNRs.size() == 0) {
        return;
    }
    
    // max size
    int maxNR = *std::max_element(mToBeNRs.begin(), mToBeNRs.end());
    int maxNC = maxNR / 2 + 1;
    
    // data allocation by fftw
    auto *rdata = InterfaceFFTW<floatT>::allocr(maxNR * HOWMANY);
    auto *cdata = InterfaceFFTW<floatT>::allocc(maxNC * HOWMANY);
    
    // hand over memory to Eigen::Map
    new (&mRMatXM) Eigen::Map<RMatXM>(rdata, maxNR, HOWMANY);
    new (&mCMatXM) Eigen::Map<CMatXM>(reinterpret_cast<cmplxT *>(cdata),
                                      maxNC, HOWMANY);
    
    // split time limit by NR
    double factor = 0.;
    for (const int &NR: mToBeNRs) {
        // use NR * log2(NR) for NR scaling
        factor += NR * std::log2(NR * 2.);
    }
    // factor * 2: R2C and C2R
    double timeUnit = timeLimitForPlanning / (factor * 2.);
    
    // create plans
    for (const int &NR: mToBeNRs) {
        // split time limit by NR
        InterfaceFFTW<floatT>::timelimit(timeUnit * NR * std::log2(NR * 2.));
        // r2c plan
        mPlansR2C.insert({NR,
            InterfaceFFTW<floatT>::planR2C(1, &NR, HOWMANY,
                                           rdata, NULL, 1, maxNR,
                                           cdata, NULL, 1, maxNC,
                                           FFTW_MEASURE)});
        // c2r plan
        mPlansC2R.insert({NR,
            InterfaceFFTW<floatT>::planC2R(1, &NR, HOWMANY,
                                           cdata, NULL, 1, maxNC,
                                           rdata, NULL, 1, maxNR,
                                           FFTW_MEASURE)});
    }
}

// clear plans
template <typename floatT, int HOWMANY>
void SolverFFTW<floatT, HOWMANY>::clearPlans() {
    // empty
    if (mPlansR2C.size() == 0) {
        return;
    }
    
    // destroy plans
    for (auto it = mPlansR2C.begin(); it != mPlansR2C.end(); ++it) {
        int NR = it->first;
        InterfaceFFTW<floatT>::destroy(mPlansR2C.at(NR));
        InterfaceFFTW<floatT>::destroy(mPlansC2R.at(NR));
    }
    mPlansR2C.clear();
    mPlansC2R.clear();
    
    // free memory
    InterfaceFFTW<floatT>::free(mRMatXM.data());
    InterfaceFFTW<floatT>::free(mCMatXM.data());
    
    // reset Eigen::Map
    new (&mRMatXM) Eigen::Map<RMatXM>(nullptr, 0, HOWMANY);
    new (&mCMatXM) Eigen::Map<CMatXM>(nullptr, 0, HOWMANY);
}

// time factor for planning
template <typename floatT, int HOWMANY>
double SolverFFTW<floatT, HOWMANY>::timeFactorForPlanning() const {
    // empty
    if (mToBeNRs.size() == 0) {
        return 0.;
    }
    
    // accumulate time factor
    double factor = 0.;
    for (const int &NR: mToBeNRs) {
        // use NR * log2(NR) for NR scaling
        factor += NR * std::log2(NR * 2.);
    }
    
    // use HOWMANY ^ 0.75 for HOWMANY scaling
    return factor * std::pow(HOWMANY * 1., 0.75);
}
