//
//  Isotropic.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 4/28/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  isotropic

#ifndef Isotropic_hpp
#define Isotropic_hpp

// mu * 2
#ifdef _SAVE_MEMORY
// computed on the fly as static variables
#define xMu2 sMu2
#else
// precomputed and stored as member variables
#define xMu2 mMu2
#endif

#include "Elastic.hpp"

class Isotropic: public Elastic {
public:
    // 1D constructor
    Isotropic(std::unique_ptr<Attenuation> &attenuation,
              const eigen::DMatPP_RM &lambda,
              const eigen::DMatPP_RM &mu):
    Elastic(true, attenuation),
    mLambda(lambda), mMu(mu)
#ifndef _SAVE_MEMORY
    , mMu2(eigen::DMatPP_RM(mu * 2.))
#endif
    {
        // nothing
    }
    
    // 3D constructor
    Isotropic(std::unique_ptr<Attenuation> &attenuation,
              const eigen::DMatXN &lambda,
              const eigen::DMatXN &mu):
    Elastic(false, attenuation),
    mLambda(lambda), mMu(mu)
#ifndef _SAVE_MEMORY
    , mMu2(eigen::DMatXN(mu * 2.))
#endif
    {
        // nothing
    }
    
    // copy constructor
    Isotropic(const Isotropic &other):
    Elastic(other),
    mLambda(other.mLambda), mMu(other.mMu)
#ifndef _SAVE_MEMORY
    , mMu2(other.mMu2)
#endif
    {
        // nothing
    }
    
    // clone for copy constructor
    std::unique_ptr<Elastic> clone() const {
        return std::make_unique<Isotropic>(*this);
    }
    
    // check compatibility
    void checkCompatibility(int nr, bool elemInFourier) const {
        // base
        Elastic::checkCompatibility(nr, elemInFourier);
        
        // parameters
        mMu.checkCompatibility(m1D, nr, elemInFourier, "mMu", "Isotropic");
        
        // workspace
#ifdef _SAVE_MEMORY
        sMu2.expandWorkspace(m1D, nr);
#endif
    }
    
    
    ///////////////////////// strain to stress /////////////////////////
    // RTZ coordinates
    bool inRTZ() const {
        return false;
    }
    
    // strain => stress in Fourier space
    void strainToStress_FR(const eigen::vec_ar6_CMatPP_RM &strain,
                           eigen::vec_ar6_CMatPP_RM &stress,
                           int nu_1) const {
        // elasticity
#ifdef _SAVE_MEMORY
        computeMu2();
#endif
        for (int alpha = 0; alpha < nu_1; alpha++) {
            strainToStress<CaseFA::_1D_FR>
            (strain, stress, alpha, mLambda, mMu, xMu2);
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
        computeMu2();
#endif
        if (m1D) {
            strainToStress<CaseFA::_1D_CD>
            (strain, stress, nr, mLambda, mMu, xMu2);
        } else {
            strainToStress<CaseFA::_3D_CD>
            (strain, stress, nr, mLambda, mMu, xMu2);
        }
        
        // attenuation
        applyAttenuation(strain, stress, nr);
    }
    
    
    //////////////////////// properties ////////////////////////
private:
    // Cijkl
    const faN::PropertyN mLambda;
    const faN::PropertyN mMu;
    
    // mu * 2
#ifdef _SAVE_MEMORY
    // computed on the fly as static variables
    inline static faN::PropertyN sMu2;
    void computeMu2() const {
        static const numerical::Real two = 2.;
        sMu2.set(m1D, mMu, two);
    }
#else
    // precomputed and stored as member variables
    const faN::PropertyN mMu2;
#endif
    
    
    /////////////////////////////////////////////////////////////////////
    /////////////////// template-based implementation ///////////////////
    /////////////////////////////////////////////////////////////////////
    
    // make public for attenuation
public:
    template <CaseFA CASE, class FMat6> static void
    strainToStress(const FMat6 &strain, FMat6 &stress, int nx,
                   const faN::PropertyN &lambda,
                   const faN::PropertyN &mu,
                   const faN::PropertyN &mu2) {
        typedef faN::FieldArithmeticN FA;
        // use block 3 to store sii
        FA::F123xP<CASE>(nx, stress, 3, strain, 0, 1, 2, lambda);
        // diagonal
        FA::R1_FxP<CASE>(nx, stress, 0, 3, strain, 0, mu2);
        FA::R1_FxP<CASE>(nx, stress, 1, 3, strain, 1, mu2);
        FA::R1_FxP<CASE>(nx, stress, 2, 3, strain, 2, mu2);
        // off-diagonal
        FA::FxP<CASE>(nx, stress, 3, strain, 3, mu);
        FA::FxP<CASE>(nx, stress, 4, strain, 4, mu);
        FA::FxP<CASE>(nx, stress, 5, strain, 5, mu);
    }
};

#endif /* Isotropic_hpp */
