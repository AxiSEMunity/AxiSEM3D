//
//  AttenuationFull.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 8/13/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  attenuation on all GLL points

#ifndef AttenuationFull_hpp
#define AttenuationFull_hpp

#include "Attenuation.hpp"
#include "Isotropic.hpp"

class AttenuationFull: public Attenuation {
public:
    // 1D constructor
    AttenuationFull(const eigen::DMatPP_RM &dLambda,
                    const eigen::DMatPP_RM &dMu):
    Attenuation(true), mDLambda(dLambda), mDMu(dMu)
#ifndef _SAVE_MEMORY
    , mDMu2(eigen::DMatPP_RM(dMu * 2.))
#endif
    {
        // nothing
    }
    
    // 3D constructor
    AttenuationFull(const eigen::DMatXN &dLambda,
                    const eigen::DMatXN &dMu):
    Attenuation(false), mDLambda(dLambda), mDMu(dMu)
#ifndef _SAVE_MEMORY
    , mDMu2(eigen::DMatXN(dMu * 2.))
#endif
    {
        // nothing
    }
    
    // copy constructor
    AttenuationFull(const AttenuationFull &other):
    Attenuation(other), mDLambda(other.mDLambda), mDMu(other.mDMu)
#ifndef _SAVE_MEMORY
    , mDMu2(other.mDMu2)
#endif
    {
        // nothing
    }
    
    // clone for copy constructor
    std::unique_ptr<Attenuation> clone() const {
        return std::make_unique<AttenuationFull>(*this);
    }
    
    // check compatibility
    void checkCompatibility(int nr, bool elemInFourier, bool elastic1D) {
        // base
        Attenuation::checkCompatibility(nr, elemInFourier, elastic1D);
        
        // memory variables
        if (elemInFourier) {
            int nu_1 = nr / 2 + 1;
            mDStress_FR.resize(nu_1);
            for (int alpha = 0; alpha < nu_1; alpha++) {
                mDStress_FR[alpha].fill(eigen::CMatPP_RM::Zero());
            }
            mDStress_CD.resize(0, spectral::nPEM * 6);
        } else {
            mDStress_FR.clear();
            mDStress_FR.shrink_to_fit();
            mDStress_CD = eigen::RMatXN6::Zero(nr, spectral::nPEM * 6);
        }
        mMemVar_FR.assign(nsls(), mDStress_FR);
        mMemVar_CD.assign(nsls(), mDStress_CD);
        
        // parameters
        mDMu.checkCompatibility(m1D, nr, elemInFourier,
                                "mDMu", "AttenuationFull");
        
        // workspace
#ifdef _SAVE_MEMORY
        sDMu2.expandWorkspace(m1D, nr);
#endif
    }
    
    // reset to zero
    void resetToZero() {
        for (int alpha = 0; alpha < mDStress_FR.size(); alpha++) {
            for (int idim = 0; idim < 6; idim++) {
                mDStress_FR[alpha][idim].setZero();
            }
        }
        mDStress_CD.setZero();
        mMemVar_FR.assign(nsls(), mDStress_FR);
        mMemVar_CD.assign(nsls(), mDStress_CD);
    }
    
    
    //////////////////////// apply ////////////////////////
    // apply attenuation in Fourier space
    void apply(const eigen::vec_ar6_CMatPP_RM &strain,
               eigen::vec_ar6_CMatPP_RM &stress, int nu_1) {
#ifdef _SAVE_MEMORY
        computeDMu2();
#endif
        for (int alpha = 0; alpha < nu_1; alpha++) {
            apply<CaseFA::_1D_FR>(strain, stress, alpha,
                                  mDStress_FR, mMemVar_FR,
                                  mDLambda, mDMu, xDMu2);
        }
    }
    
    // apply attenuation in cardinal space
    void apply(const eigen::RMatXN6 &strain,
               eigen::RMatXN6 &stress, int nr) {
#ifdef _SAVE_MEMORY
        computeDMu2();
#endif
        if (m1D) {
            apply<CaseFA::_1D_CD>(strain, stress, nr,
                                  mDStress_CD, mMemVar_CD,
                                  mDLambda, mDMu, xDMu2);
        } else {
            apply<CaseFA::_3D_CD>(strain, stress, nr,
                                  mDStress_CD, mMemVar_CD,
                                  mDLambda, mDMu, xDMu2);
        }
    }
    
    
    ///////////////////////// properties /////////////////////////
private:
    // Cijkl
    const faN::PropertyN mDLambda;
    const faN::PropertyN mDMu;
    
    // dmu * 2
#ifdef _SAVE_MEMORY
    // computed on the fly as static variables
    inline static faN::PropertyN sDMu2;
    void computeDMu2() const {
        static const numerical::Real two = 2.;
        sDMu2.set(m1D, mDMu, two);
    }
#else
    // precomputed and stored as member variables
    const faN::PropertyN mDMu2;
#endif
    
    // memory variables in Fourier space
    eigen::vec_ar6_CMatPP_RM mDStress_FR;
    std::vector<eigen::vec_ar6_CMatPP_RM> mMemVar_FR;
    
    // memory variables in cardinal space
    eigen::RMatXN6 mDStress_CD;
    std::vector<eigen::RMatXN6> mMemVar_CD;
    
    
    /////////////////////////////////////////////////////////////////////
    /////////////////// template-based implementation ///////////////////
    /////////////////////////////////////////////////////////////////////
    
    /////////////////// stress -= mem ///////////////////
    // CASE == _1D_FR
    template <CaseFA CASE, class FMat6> static
    typename std::enable_if<CASE == CaseFA::_1D_FR, void>::type
    subtractMemFromStress(const std::vector<FMat6> &memVar, FMat6 &stress,
                          int nx) {
        for (int isls = 0; isls < nsls(); isls++) {
            for (int idim = 0; idim < 6; idim++) {
                stress[nx][idim] -= memVar[isls][nx][idim];
            }
        }
    }
    // CASE != _1D_FR
    template <CaseFA CASE, class FMat6> static
    typename std::enable_if<CASE != CaseFA::_1D_FR, void>::type
    subtractMemFromStress(const std::vector<FMat6> &memVar, FMat6 &stress,
                          int nx) {
        for (int isls = 0; isls < nsls(); isls++) {
            stress.topRows(nx) -= memVar[isls];
        }
    }
    
    
    /////////////////// apply ///////////////////
    // CASE == _1D_FR
    template <CaseFA CASE, class FMat6>
    static void apply(const FMat6 &strain, FMat6 &stress, int nx,
                      FMat6 &dStress, std::vector<FMat6> &memVar,
                      const faN::PropertyN &dLambda,
                      const faN::PropertyN &dMu,
                      const faN::PropertyN &dMu2) {
        // stress -= memory
        subtractMemFromStress<CASE>(memVar, stress, nx);
        
        // update memory variables, phase 1
        Attenuation::updateMemAlphaBeta<CASE>(dStress, memVar, nx);
        
        // update memory stress
        Isotropic::strainToStress<CASE>(strain, dStress, nx,
                                        dLambda, dMu, dMu2);
        
        // update memory variables, phase 2
        Attenuation::updateMemGamma<CASE>(dStress, memVar, nx);
    }
};

#endif /* AttenuationFull_hpp */
