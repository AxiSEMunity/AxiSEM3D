//
//  Acoustic.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 2/25/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  acoustic material

#ifndef Acoustic_hpp
#define Acoustic_hpp

#include "eigen_element.hpp"
#include "FieldArithmetic.hpp"

class Acoustic {
public:
    // 1D constructor
    Acoustic(const eigen::DMatPP_RM &K):
    m1D(true), mK(K) {
        // nothing
    }
    
    // 3D constructor
    Acoustic(const eigen::DMatXN &K):
    m1D(false), mK(K) {
        // nothing
    }
    
    // copy constructor
    Acoustic(const Acoustic &other):
    m1D(other.m1D), mK(other.mK) {
        // nothing
    }
    
    // 1D operation
    bool is1D() const {
        return m1D;
    }
    
    // check compatibility
    void checkCompatibility(int nr, bool elemInFourier) const {
        // parameters
        mK.checkCompatibility(m1D, nr, elemInFourier, "mK", "Acoustic");
    }
    
    
    ///////////////////////// strain to stress /////////////////////////
    // strain => stress in Fourier space
    void strainToStress_FR(const eigen::vec_ar3_CMatPP_RM &strain,
                           eigen::vec_ar3_CMatPP_RM &stress,
                           int nu_1) const {
        for (int alpha = 0; alpha < nu_1; alpha++) {
            strainToStress<CaseFA::_1D_FR>(strain, stress, alpha, mK);
        }
    }
    
    // strain => stress in cardinal space
    void strainToStress_CD(const eigen::RMatXN3 &strain,
                           eigen::RMatXN3 &stress,
                           int nr) const {
        if (m1D) {
            strainToStress<CaseFA::_1D_CD>(strain, stress, nr, mK);
        } else {
            strainToStress<CaseFA::_3D_CD>(strain, stress, nr, mK);
        }
    }
    
    
    ///////////////////////// properties /////////////////////////
private:
    // 1D/3D flag
    const bool m1D;
    
    // K = 1 / rho
    const faN::PropertyN mK;
    
    
    /////////////////////////////////////////////////////////////////////
    /////////////////// template-based implementation ///////////////////
    /////////////////////////////////////////////////////////////////////
    
    // strain => stress
    template <CaseFA CASE, class FMat3> static void
    strainToStress(const FMat3 &strain, FMat3 &stress, int nx,
                   const faN::PropertyN &K) {
        typedef faN::FieldArithmeticN FA;
        FA::FxP<CASE>(nx, stress, 0, strain, 0, K);
        FA::FxP<CASE>(nx, stress, 1, strain, 1, K);
        FA::FxP<CASE>(nx, stress, 2, strain, 2, K);
    }
};

#endif /* Acoustic_hpp */
