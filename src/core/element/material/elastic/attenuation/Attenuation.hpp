//
//  Attenuation.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 2/25/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  attenuation

#ifndef Attenuation_hpp
#define Attenuation_hpp

// dmu * 2
#ifdef _SAVE_MEMORY
// computed on the fly as static variables
#define xDMu2 sDMu2
#else
// precomputed and stored as member variables
#define xDMu2 mDMu2
#endif

#include "eigen_element.hpp"
#include "FieldArithmetic.hpp"

namespace eigen {
    // alpha, beta, gamma
    typedef Eigen::Matrix<double, Eigen::Dynamic, 1> DColX;
    typedef Eigen::Matrix<numerical::Real, Eigen::Dynamic, 1> RColX;
}

class Attenuation {
public:
    // constructor
    Attenuation(bool is1D): m1D(is1D) {
        // nothing
    }
    
    // copy constructor
    Attenuation(const Attenuation &other): m1D(other.m1D) {
        // nothing
    }
    
    // destructor
    virtual ~Attenuation() = default;
    
    // clone for copy constructor
    virtual std::unique_ptr<Attenuation> clone() const = 0;
    
    // check compatibility
    virtual void
    checkCompatibility(int nr, bool elemInFourier, bool elastic1D) {
        // 1D/3D
        if (elastic1D != m1D) {
            throw std::runtime_error("AttenuationFull::checkCompatibility || "
                                     "Incompatible 1D/3D flags.");
        }
    }
    
    // reset to zero
    virtual void resetToZero() = 0;
    
    
    //////////////////////// apply ////////////////////////
    // apply attenuation in Fourier space
    virtual void apply(const eigen::vec_ar6_CMatPP_RM &strain,
                       eigen::vec_ar6_CMatPP_RM &stress, int nu_1) = 0;
    
    // apply attenuation in cardinal space
    virtual void apply(const eigen::RMatXN6 &strain,
                       eigen::RMatXN6 &stress, int nr) = 0;
    
    
    //////////////////////// data ////////////////////////
protected:
    // 1D/3D flag
    const bool m1D;
    
    
    /////////////////////////////////////////////////////////////////////
    /////////////////// template-based implementation ///////////////////
    /////////////////////////////////////////////////////////////////////
    
protected:
    // get number of SLS
    static int nsls() {
        return (int)sAlpha.rows();
    }
    
    ///////////////// update memory variable, phase 1 /////////////////
    // CASE == _1D_FR
    template <CaseFA CASE, class FMat6> static
    typename std::enable_if<CASE == CaseFA::_1D_FR, void>::type
    updateMemAlphaBeta(const FMat6 &dStress, std::vector<FMat6> &memVar,
                       int nx) {
        for (int isls = 0; isls < nsls(); isls++) {
            for (int idim = 0; idim < 6; idim++) {
                memVar[isls][nx][idim] =
                sAlpha[isls] * memVar[isls][nx][idim] +
                sBetta[isls] * dStress[nx][idim];
            }
        }
    }
    // CASE != _1D_FR
    template <CaseFA CASE, class FMat6> static
    typename std::enable_if<CASE != CaseFA::_1D_FR, void>::type
    updateMemAlphaBeta(const FMat6 &dStress, std::vector<FMat6> &memVar,
                       int nx) {
        for (int isls = 0; isls < nsls(); isls++) {
            memVar[isls] = (sAlpha[isls] * memVar[isls] +
                            sBetta[isls] * dStress);
        }
    }
    
    
    ///////////////// update memory variable, phase 2 /////////////////
    // CASE == _1D_FR
    template <CaseFA CASE, class FMat6> static
    typename std::enable_if<CASE == CaseFA::_1D_FR, void>::type
    updateMemGamma(const FMat6 &dStress, std::vector<FMat6> &memVar, int nx) {
        for (int isls = 0; isls < nsls(); isls++) {
            for (int idim = 0; idim < 6; idim++) {
                memVar[isls][nx][idim] += sGamma[isls] * dStress[nx][idim];
            }
        }
    }
    // CASE != _1D_FR
    template <CaseFA CASE, class FMat6> static
    typename std::enable_if<CASE != CaseFA::_1D_FR, void>::type
    updateMemGamma(const FMat6 &dStress, std::vector<FMat6> &memVar, int nx) {
        for (int isls = 0; isls < nsls(); isls++) {
            memVar[isls] += sGamma[isls] * dStress;
        }
    }
    
    
    ////////////////////////////////////////
    //////////////// static ////////////////
    ////////////////////////////////////////
    
public:
    // set alpha, beta, gamma
    static void setAlphaBetaGamma(const eigen::DColX &alpha,
                                  const eigen::DColX &betta,
                                  const eigen::DColX &gamma) {
        sAlpha = alpha.cast<numerical::Real>();
        sBetta = betta.cast<numerical::Real>();
        sGamma = gamma.cast<numerical::Real>();
    }
    
private:
    // coeffs
    inline static eigen::RColX sAlpha;
    inline static eigen::RColX sBetta;
    inline static eigen::RColX sGamma;
};

#endif /* Attenuation_hpp */
