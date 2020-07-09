//
//  TransverselyIsotropic.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 4/28/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  transversely isotropic

#ifndef TransverselyIsotropic_hpp
#define TransverselyIsotropic_hpp

// - N * 2
#ifdef _SAVE_MEMORY
// computed on the fly as static variables
#define xMinusN2 sMinusN2
#else
// precomputed and stored as member variables
#define xMinusN2 mMinusN2
#endif

#include "Elastic.hpp"

class TransverselyIsotropic: public Elastic {
public:
    // 1D constructor
    TransverselyIsotropic(std::unique_ptr<Attenuation> &attenuation,
                          const eigen::DMatPP_RM &A,
                          const eigen::DMatPP_RM &C,
                          const eigen::DMatPP_RM &F,
                          const eigen::DMatPP_RM &L,
                          const eigen::DMatPP_RM &N):
    Elastic(true, attenuation),
    mA(A), mC(C), mF(F), mL(L), mN(N)
#ifndef _SAVE_MEMORY
    , mMinusN2(eigen::DMatPP_RM(N * (-2.)))
#endif
    {
        // nothing
    }
    
    // 3D constructor
    TransverselyIsotropic(std::unique_ptr<Attenuation> &attenuation,
                          const eigen::DMatXN &A,
                          const eigen::DMatXN &C,
                          const eigen::DMatXN &F,
                          const eigen::DMatXN &L,
                          const eigen::DMatXN &N):
    Elastic(false, attenuation),
    mA(A), mC(C), mF(F), mL(L), mN(N)
#ifndef _SAVE_MEMORY
    , mMinusN2(eigen::DMatXN(N * (-2.)))
#endif
    {
        // nothing
    }
    
    // copy constructor
    TransverselyIsotropic(const TransverselyIsotropic &other):
    Elastic(other),
    mA(other.mA), mC(other.mC), mF(other.mF), mL(other.mL), mN(other.mN)
#ifndef _SAVE_MEMORY
    , mMinusN2(other.mMinusN2)
#endif
    {
        // nothing
    }
    
    // check compatibility
    void checkCompatibility(int nr, bool elemInFourier) const {
        // base
        Elastic::checkCompatibility(nr, elemInFourier);
        
        // parameters
        mA.checkCompatibility(m1D, nr, elemInFourier,
                              "mA", "TransverselyIsotropic");
        
        // workspace
#ifdef _SAVE_MEMORY
        sMinusN2.expandWorkspace(m1D, nr);
#endif
    }
    
    // clone for copy constructor
    std::unique_ptr<Elastic> clone() const {
        return std::make_unique<TransverselyIsotropic>(*this);
    }
    
    ///////////////////////// strain to stress /////////////////////////
    // RTZ coordinates
    bool inRTZ() const {
        return true;
    }
    
    // strain => stress in Fourier space
    void strainToStress_FR(const eigen::vec_ar6_CMatPP_RM &strain,
                           eigen::vec_ar6_CMatPP_RM &stress,
                           int nu_1) const {
        // elasticity
#ifdef _SAVE_MEMORY
        computeMinusN2();
#endif
        for (int alpha = 0; alpha < nu_1; alpha++) {
            strainToStress<CaseFA::_1D_FR>
            (strain, stress, alpha,
             mA, mC, mF, mL, mN, xMinusN2);
        }
        
        // attenuation
        applyAttenuation(strain, stress, nu_1);
    }
    
    // strain => stress in cardinal space
    void strainToStress_CD(const eigen::RMatXN6 &strain,
                           eigen::RMatXN6 &stress,
                           int nr) const {
        // elasticity
#ifdef _SAVE_MEMORY
        computeMinusN2();
#endif
        if (m1D) {
            strainToStress<CaseFA::_1D_CD>
            (strain, stress, nr,
             mA, mC, mF, mL, mN, xMinusN2);
        } else {
            strainToStress<CaseFA::_3D_CD>
            (strain, stress, nr,
             mA, mC, mF, mL, mN, xMinusN2);
        }
        
        // attenuation
        applyAttenuation(strain, stress, nr);
    }
    
    
    //////////////////////// properties ////////////////////////
private:
    // Cijkl
    const faN::PropertyN mA, mC, mF, mL, mN;
    
    // - N * 2
#ifdef _SAVE_MEMORY
    // computed on the fly as static variables
    inline static faN::PropertyN sMinusN2;
    void computeMinusN2() const {
        static const numerical::Real minusTwo = -2.;
        sMinusN2.set(m1D, mN, minusTwo);
    }
#else
    // precomputed and stored as member variables
    const faN::PropertyN mMinusN2;
#endif
    
    
    /////////////////////////////////////////////////////////////////////
    /////////////////// template-based implementation ///////////////////
    /////////////////////////////////////////////////////////////////////
    
    template <CaseFA CASE, class FMat6> static void
    strainToStress(const FMat6 &strain, FMat6 &stress, int nx,
                   const faN::PropertyN &A,
                   const faN::PropertyN &C,
                   const faN::PropertyN &F,
                   const faN::PropertyN &L,
                   const faN::PropertyN &N,
                   const faN::PropertyN &minusN2) {
        typedef faN::FieldArithmeticN FA;
        // use block 3, 4 as temp
        FA::F1_F2<CASE>(nx, stress, 3, strain, 0, 1);
        FA::R1xP1_F2xP2<CASE>(nx, stress, 4, 3, A, strain, 2, F);
        // diagonal
        FA::R1_FxP<CASE>(nx, stress, 0, 4, strain, 1, minusN2);
        FA::R1_FxP<CASE>(nx, stress, 1, 4, strain, 0, minusN2);
        FA::R1xP1_F2xP2<CASE>(nx, stress, 2, 3, F, strain, 2, C);
        // off-diagonal
        FA::FxP<CASE>(nx, stress, 3, strain, 3, L);
        FA::FxP<CASE>(nx, stress, 4, strain, 4, L);
        FA::FxP<CASE>(nx, stress, 5, strain, 5, N);
    }
};

#endif /* TransverselyIsotropic_hpp */
