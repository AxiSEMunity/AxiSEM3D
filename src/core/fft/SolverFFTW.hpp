//
//  SolverFFTW.hpp
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

#ifndef SolverFFTW_hpp
#define SolverFFTW_hpp

extern "C" {
#include <fftw3.h>
}
#include <set>
#include <map>
#include <vector>
#include "mapPPvsN.hpp"

///////////////////// fftw3 interface /////////////////////
template <typename floatT>
struct InterfaceFFTW {
};

template <>
struct InterfaceFFTW<double> {
    typedef fftw_plan Plan;
    static constexpr auto planR2C = fftw_plan_many_dft_r2c;
    static constexpr auto planC2R = fftw_plan_many_dft_c2r;
    static constexpr auto destroy = fftw_destroy_plan;
    static constexpr auto execute = fftw_execute;
    static constexpr auto timelimit = fftw_set_timelimit;
    static constexpr auto allocr = fftw_alloc_real;
    static constexpr auto allocc = fftw_alloc_complex;
    static constexpr auto free = fftw_free;
};

template <>
struct InterfaceFFTW<float> {
    typedef fftwf_plan Plan;
    static constexpr auto planR2C = fftwf_plan_many_dft_r2c;
    static constexpr auto planC2R = fftwf_plan_many_dft_c2r;
    static constexpr auto destroy = fftwf_destroy_plan;
    static constexpr auto execute = fftwf_execute;
    static constexpr auto timelimit = fftwf_set_timelimit;
    static constexpr auto allocr = fftwf_alloc_real;
    static constexpr auto allocc = fftwf_alloc_complex;
    static constexpr auto free = fftwf_free;
};


///////////////////// template class for fftw3 /////////////////////
template <typename floatT, int HOWMANY>
class SolverFFTW {
    //////////////////// buffer type ////////////////////
    typedef std::complex<floatT> cmplxT;
    typedef Eigen::Matrix<floatT, Eigen::Dynamic, HOWMANY> RMatXM;
    typedef Eigen::Matrix<cmplxT, Eigen::Dynamic, HOWMANY> CMatXM;
    
    
    //////////////////// internal methods ////////////////////
public:
    // destructor
    ~SolverFFTW() {
        clearPlans();
    }
    
    // create plans after adding ALL NRs
    void createPlans(double timeLimitForPlanning);
    
    // clear plans
    void clearPlans();
    
    // time factor for planning
    double timeFactorForPlanning() const;
    
    // get planned NRs for verbose
    std::vector<int> getPlannedNRs() const {
        std::vector<int> plannedNRs;
        std::transform(mPlansR2C.begin(), mPlansR2C.end(),
                       std::back_inserter(plannedNRs),
                       [](auto pair) {return pair.first;});
        return plannedNRs;
    }
    
    
    //////////////////// external methods ////////////////////
public:
    // add an NR
    void addNR(int NR) {
        mToBeNRs.insert(NR);
    }
    
    // forward, real => complex, with flattened output
    void computeR2C(const RMatXM &r_in, CMatXM &c_out, int NR) {
        // feed input with the FFTW convention for normalization
        mRMatXM.topRows(NR) = r_in.topRows(NR) * (floatT)(1. / NR);
        
        // perform fft
        InterfaceFFTW<floatT>::execute(mPlansR2C.at(NR));
        
        // Nyquist
        int NC = NR / 2 + 1;
        if (NR % 2 == 0) {
            mCMatXM.row(NC - 1) *= (floatT).5;
        }
        
        // retrieve output
        c_out.topRows(NC) = mCMatXM.topRows(NC);
    }
    
    // forward, real => complex, with structured output
    template <typename vec_arD_MatPP_RM>
    void computeR2C(const RMatXM &r_in, vec_arD_MatPP_RM &c_out, int NR) {
        // feed input with the FFTW convention for normalisation
        mRMatXM.topRows(NR) = r_in.topRows(NR) * (floatT)(1. / NR);
        
        // perform fft
        InterfaceFFTW<floatT>::execute(mPlansR2C.at(NR));
        
        // Nyquist
        int NC = NR / 2 + 1;
        if (NR % 2 == 0) {
            mCMatXM.row(NC - 1) *= (floatT).5;
        }
        
        // retrieve output
        mapPPvsN::N2PP(mCMatXM, c_out, NC);
    }
    
    // backward, complex => real, with flattened input
    void computeC2R(const CMatXM &c_in, RMatXM &r_out, int NR) {
        // feed input
        int NC = NR / 2 + 1;
        mCMatXM.topRows(NC) = c_in.topRows(NC);
        
        // Nyquist
        if (NR % 2 == 0) {
            mCMatXM.row(NC - 1) *= (floatT)2.;
        }
        
        // perform fft
        InterfaceFFTW<floatT>::execute(mPlansC2R.at(NR));
        
        // retrieve output
        r_out.topRows(NR) = mRMatXM.topRows(NR);
    }
    
    // backward, complex => real, with structured input
    template <typename vec_arD_MatPP_RM>
    void computeC2R(const vec_arD_MatPP_RM &c_in, RMatXM &r_out, int NR) {
        // feed input
        int NC = NR / 2 + 1;
        mapPPvsN::PP2N(c_in, mCMatXM, NC);
        
        // Nyquist
        if (NR % 2 == 0) {
            mCMatXM.row(NC - 1) *= (floatT)2.;
        }
        
        // perform fft
        InterfaceFFTW<floatT>::execute(mPlansC2R.at(NR));
        
        // retrieve output
        r_out.topRows(NR) = mRMatXM.topRows(NR);
    }
    
    
    //////////////////// variables ////////////////////
private:
    // to-be NR list
    std::set<int> mToBeNRs;
    
    // plans
    std::map<int, typename InterfaceFFTW<floatT>::Plan> mPlansR2C;
    std::map<int, typename InterfaceFFTW<floatT>::Plan> mPlansC2R;
    
    // buffer data
    Eigen::Map<RMatXM> mRMatXM = Eigen::Map<RMatXM>(nullptr, 0, HOWMANY);
    Eigen::Map<CMatXM> mCMatXM = Eigen::Map<CMatXM>(nullptr, 0, HOWMANY);
};

#endif /* SolverFFTW_hpp */
