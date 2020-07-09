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
    
    // data allocation
    mRMatXM.resize(maxNR, HOWMANY);
    mCMatXM.resize(maxNC, HOWMANY);
    
    // split time limit by NR
    double factor = 0.;
    for (const int &NR: mToBeNRs) {
        // use NR * log2(NR) for NR scaling
        factor += NR * std::log2(NR * 1.);
    }
    // factor * 2: R2C and C2R
    double timeUnit = timeLimitForPlanning / (factor * 2.);
    
    // create plans
    for (const int &NR: mToBeNRs) {
        // split time limit by NR
        InterfaceFFTW<floatT>::timelimit(timeUnit * NR * std::log2(NR * 1.));
        // r2c plan
        mPlansR2C.insert({NR,
            InterfaceFFTW<floatT>::planR2C(1, &NR, HOWMANY,
                                           (mRMatXM.data()), NULL, 1, maxNR,
                                           reinterpret_cast<typename
                                           InterfaceFFTW<floatT>::fftw_cmplx *>
                                           (mCMatXM.data()), NULL, 1, maxNC,
                                           FFTW_PATIENT)});
        // c2r plan
        mPlansC2R.insert({NR,
            InterfaceFFTW<floatT>::planC2R(1, &NR, HOWMANY,
                                           reinterpret_cast<typename
                                           InterfaceFFTW<floatT>::fftw_cmplx *>
                                           (mCMatXM.data()), NULL, 1, maxNC,
                                           (mRMatXM.data()), NULL, 1, maxNR,
                                           FFTW_PATIENT)});
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
    
    // data
    mRMatXM.resize(0, HOWMANY);
    mCMatXM.resize(0, HOWMANY);
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
        factor += NR * std::log2(NR * 1.);
    }
    
    // use HOWMANY ^ 0.75 for HOWMANY scaling
    return factor * std::pow(HOWMANY * 1., 0.75);
}
