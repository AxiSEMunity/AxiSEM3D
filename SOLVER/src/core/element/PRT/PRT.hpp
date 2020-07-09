//
//  PRT.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 2/25/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  particle relabelling transformation

#ifndef PRT_hpp
#define PRT_hpp

// J^-1 components scaled by |J|
#ifdef _SAVE_MEMORY
// computed on the fly as static variables
#define xX0_J sX0_J
#define xX1_J sX1_J
#define xX2_J sX2_J
#define xX3_J sX3_J
#else
// precomputed and stored as member variables
#define xX0_J mX0_J
#define xX1_J mX1_J
#define xX2_J mX2_J
#define xX3_J mX3_J
#endif

#include "eigen_element.hpp"
#include "FieldArithmetic.hpp"

class PRT {
public:
#ifdef _SAVE_MEMORY
    // 1D constructor
    PRT(const eigen::DMatPP_RM &X0, const eigen::DMatPP_RM &X1,
        const eigen::DMatPP_RM &X2, const eigen::DMatPP_RM &X3,
        const eigen::DMatPP_RM &XJ): m1D(true),
    mX0(X0), mX1(X1), mX2(X2), mX3(X3), mXJ(XJ) {
        // noting
    }
    
    // 3D constructor
    PRT(const eigen::DMatXN &X0, const eigen::DMatXN &X1,
        const eigen::DMatXN &X2, const eigen::DMatXN &X3,
        const eigen::DMatXN &XJ): m1D(false),
    mX0(X0), mX1(X1), mX2(X2), mX3(X3), mXJ(XJ) {
        // noting
    }
    
    // copy constructor
    PRT(const PRT &other): m1D(other.m1D),
    mX0(other.mX0), mX1(other.mX1), mX2(other.mX2), mX3(other.mX3),
    mXJ(other.mXJ) {
        // noting
    }
#else
    // 1D constructor
    PRT(const eigen::DMatPP_RM &X0, const eigen::DMatPP_RM &X1,
        const eigen::DMatPP_RM &X2, const eigen::DMatPP_RM &X3,
        const eigen::DMatPP_RM &XJ): m1D(true),
    mX0(X0), mX1(X1), mX2(X2), mX3(X3),
    mX0_J(eigen::DMatPP_RM(X0.cwiseProduct(XJ))),
    mX1_J(eigen::DMatPP_RM(X1.cwiseProduct(XJ))),
    mX2_J(eigen::DMatPP_RM(X2.cwiseProduct(XJ))),
    mX3_J(eigen::DMatPP_RM(X3.cwiseProduct(XJ))) {
        // noting
    }
    
    // 3D constructor
    PRT(const eigen::DMatXN &X0, const eigen::DMatXN &X1,
        const eigen::DMatXN &X2, const eigen::DMatXN &X3,
        const eigen::DMatXN &XJ): m1D(false),
    mX0(X0), mX1(X1), mX2(X2), mX3(X3),
    mX0_J(eigen::DMatXN(X0.cwiseProduct(XJ))),
    mX1_J(eigen::DMatXN(X1.cwiseProduct(XJ))),
    mX2_J(eigen::DMatXN(X2.cwiseProduct(XJ))),
    mX3_J(eigen::DMatXN(X3.cwiseProduct(XJ))) {
        // noting
    }
    
    // copy constructor
    PRT(const PRT &other): m1D(other.m1D),
    mX0(other.mX0), mX1(other.mX1), mX2(other.mX2), mX3(other.mX3),
    mX0_J(other.mX0_J),
    mX1_J(other.mX1_J),
    mX2_J(other.mX2_J),
    mX3_J(other.mX3_J) {
        // noting
    }
#endif
    
    // 1D operation
    bool is1D() const {
        return m1D;
    }
    
    // check compatibility
    void checkCompatibility(int nr, bool elemInFourier) const {
        // parameters
        mX0.checkCompatibility(m1D, nr, elemInFourier, "mX0", "PRT");
        
        // workspace
#ifdef _SAVE_MEMORY
        sX0_J.expandWorkspace(m1D, nr);
        sX1_J.expandWorkspace(m1D, nr);
        sX2_J.expandWorkspace(m1D, nr);
        sX3_J.expandWorkspace(m1D, nr);
#endif
    }
    
    
    //////////////////////// Fourier space ////////////////////////
    // fluid, sph -> und
    void sphericalToUndulated3_FR(const eigen::vec_ar3_CMatPP_RM &sph3,
                                  eigen::vec_ar3_CMatPP_RM &und3,
                                  int nu_1) const {
        for (int alpha = 0; alpha < nu_1; alpha++) {
            sphericalToUndulated3<CaseFA::_1D_FR>
            (sph3, und3, alpha, mX0, mX1, mX2, mX3);
        }
    }
    
    // solid, sph -> und
    void sphericalToUndulated6_FR(const eigen::vec_ar9_CMatPP_RM &sph9,
                                  eigen::vec_ar6_CMatPP_RM &und6,
                                  int nu_1) const {
        for (int alpha = 0; alpha < nu_1; alpha++) {
            sphericalToUndulated6<CaseFA::_1D_FR>
            (sph9, und6, alpha, mX0, mX1, mX2, mX3);
        }
    }
    
    // solid, sph -> und, for curl computation (9->9)
    void sphericalToUndulated9_FR(const eigen::vec_ar9_CMatPP_RM &sph9,
                                  eigen::vec_ar9_CMatPP_RM &und9,
                                  int nu_1) const {
        for (int alpha = 0; alpha < nu_1; alpha++) {
            sphericalToUndulated9<CaseFA::_1D_FR>
            (sph9, und9, alpha, mX0, mX1, mX2, mX3);
        }
    }
    
    // fluid, und -> sph
    void undulatedToSpherical3_FR(const eigen::vec_ar3_CMatPP_RM &und3,
                                  eigen::vec_ar3_CMatPP_RM &sph3,
                                  int nu_1) const {
#ifdef _SAVE_MEMORY
        computeScaledX();
#endif
        for (int alpha = 0; alpha < nu_1; alpha++) {
            undulatedToSpherical3<CaseFA::_1D_FR>
            (und3, sph3, alpha, xX0_J, xX1_J, xX2_J, xX3_J);
        }
    }
    
    // solid, und -> sph
    void undulatedToSpherical6_FR(const eigen::vec_ar6_CMatPP_RM &und6,
                                  eigen::vec_ar9_CMatPP_RM &sph9,
                                  int nu_1) const {
#ifdef _SAVE_MEMORY
        computeScaledX();
#endif
        for (int alpha = 0; alpha < nu_1; alpha++) {
            undulatedToSpherical6<CaseFA::_1D_FR>
            (und6, sph9, alpha, xX0_J, xX1_J, xX2_J, xX3_J);
        }
    }
    
    // solid, und -> sph, no integration, for moment tensor
    void
    undulatedToSpherical6_NoIntegration_FR(const eigen::vec_ar6_CMatPP_RM &und6,
                                           eigen::vec_ar9_CMatPP_RM &sph9,
                                           int nu_1) const {
        for (int alpha = 0; alpha < nu_1; alpha++) {
            undulatedToSpherical6<CaseFA::_1D_FR>
            (und6, sph9, alpha, mX0, mX1, mX2, mX3);
        }
    }
    
    
    //////////////////////// cardinal space ////////////////////////
    // fluid, sph -> und
    void sphericalToUndulated3_CD(const eigen::RMatXN3 &sph3,
                                  eigen::RMatXN3 &und3,
                                  int nr) const {
        if (m1D) {
            sphericalToUndulated3<CaseFA::_1D_CD>
            (sph3, und3, nr, mX0, mX1, mX2, mX3);
        } else {
            sphericalToUndulated3<CaseFA::_3D_CD>
            (sph3, und3, nr, mX0, mX1, mX2, mX3);
        }
    }
    
    // solid, sph -> und
    void sphericalToUndulated6_CD(const eigen::RMatXN9 &sph9,
                                  eigen::RMatXN6 &und6,
                                  int nr) const {
        if (m1D) {
            sphericalToUndulated6<CaseFA::_1D_CD>
            (sph9, und6, nr, mX0, mX1, mX2, mX3);
        } else {
            sphericalToUndulated6<CaseFA::_3D_CD>
            (sph9, und6, nr, mX0, mX1, mX2, mX3);
        }
    }
    
    // solid, sph -> und, for curl computation (9->9)
    void sphericalToUndulated9_CD(const eigen::RMatXN9 &sph9,
                                  eigen::RMatXN9 &und9,
                                  int nr) const {
        if (m1D) {
            sphericalToUndulated9<CaseFA::_1D_CD>
            (sph9, und9, nr, mX0, mX1, mX2, mX3);
        } else {
            sphericalToUndulated9<CaseFA::_3D_CD>
            (sph9, und9, nr, mX0, mX1, mX2, mX3);
        }
    }
    
    // fluid, und -> sph
    void undulatedToSpherical3_CD(const eigen::RMatXN3 &und3,
                                  eigen::RMatXN3 &sph3,
                                  int nr) const {
#ifdef _SAVE_MEMORY
        computeScaledX();
#endif
        if (m1D) {
            undulatedToSpherical3<CaseFA::_1D_CD>
            (und3, sph3, nr, xX0_J, xX1_J, xX2_J, xX3_J);
        } else {
            undulatedToSpherical3<CaseFA::_3D_CD>
            (und3, sph3, nr, xX0_J, xX1_J, xX2_J, xX3_J);
        }
    }
    
    // solid, und -> sph
    void undulatedToSpherical6_CD(const eigen::RMatXN6 &und6,
                                  eigen::RMatXN9 &sph9,
                                  int nr) const {
#ifdef _SAVE_MEMORY
        computeScaledX();
#endif
        if (m1D) {
            undulatedToSpherical6<CaseFA::_1D_CD>
            (und6, sph9, nr, xX0_J, xX1_J, xX2_J, xX3_J);
        } else {
            undulatedToSpherical6<CaseFA::_3D_CD>
            (und6, sph9, nr, xX0_J, xX1_J, xX2_J, xX3_J);
        }
    }
    
    // solid, und -> sph, no integration, for moment tensor
    void
    undulatedToSpherical6_NoIntegration_CD(const eigen::RMatXN6 &und6,
                                           eigen::RMatXN9 &sph9,
                                           int nr) const {
        if (m1D) {
            undulatedToSpherical6<CaseFA::_1D_CD>
            (und6, sph9, nr, mX0, mX1, mX2, mX3);
        } else {
            undulatedToSpherical6<CaseFA::_3D_CD>
            (und6, sph9, nr, mX0, mX1, mX2, mX3);
        }
    }
    
    
    ///////////////////////// properties /////////////////////////
private:
    // 1D/3D flag
    const bool m1D;
    
    // J^-1 components
    const faN::PropertyN mX0;
    const faN::PropertyN mX1;
    const faN::PropertyN mX2;
    const faN::PropertyN mX3;
    
#ifdef _SAVE_MEMORY
    // |J|
    const faN::PropertyN mXJ;
    
    // J^-1 components scaled by |J|
    inline static faN::PropertyN sX0_J;
    inline static faN::PropertyN sX1_J;
    inline static faN::PropertyN sX2_J;
    inline static faN::PropertyN sX3_J;
    
    // compute scaled J^-1 on the fly
    void computeScaledX() const {
        sX0_J.set(m1D, mX0, mXJ);
        sX1_J.set(m1D, mX1, mXJ);
        sX2_J.set(m1D, mX2, mXJ);
        sX3_J.set(m1D, mX3, mXJ);
    }
#else
    // J^-1 components scaled by |J|
    const faN::PropertyN mX0_J;
    const faN::PropertyN mX1_J;
    const faN::PropertyN mX2_J;
    const faN::PropertyN mX3_J;
#endif
    
    
    /////////////////////////////////////////////////////////////////////
    /////////////////// template-based implementation ///////////////////
    /////////////////////////////////////////////////////////////////////
    
private:
    // fluid, sph -> und
    template <CaseFA CASE, class FMat3>
    static void sphericalToUndulated3(const FMat3 &sph3, FMat3 &und3, int nx,
                                      const faN::PropertyN &X0,
                                      const faN::PropertyN &X1,
                                      const faN::PropertyN &X2,
                                      const faN::PropertyN &X3) {
        typedef faN::FieldArithmeticN FA;
        FA::F1xP1_F2xP2<CASE>(nx, und3, 0, sph3, 0, X0, 2, X1);
        FA::F1xP1_F2xP2<CASE>(nx, und3, 1, sph3, 1, X0, 2, X1);
        FA::FxP<CASE>        (nx, und3, 2, sph3, 2, X3);
    }
    
    // solid, sph -> und
    template <CaseFA CASE, class FMat9, class FMat6>
    static void sphericalToUndulated6(const FMat9 &sph9, FMat6 &und6, int nx,
                                      const faN::PropertyN &X0,
                                      const faN::PropertyN &X1,
                                      const faN::PropertyN &X2,
                                      const faN::PropertyN &X3) {
        typedef faN::FieldArithmeticN FA;
        FA::F1xP1_F2xP2<CASE>      (nx, und6, 0, sph9, 0, X0, 2, X1);
        FA::F1xP1_F2xP2<CASE>      (nx, und6, 1, sph9, 4, X0, 5, X2);
        FA::FxP<CASE>              (nx, und6, 2, sph9, 8, X3);
        FA::F1xP1_F2xP2_F3xP3<CASE>(nx, und6, 3, sph9, 7, X0, 8, X2, 5, X3);
        FA::F1xP1_F2xP2_F3xP3<CASE>(nx, und6, 4, sph9, 6, X0, 8, X1, 2, X3);
        FA::
        F1abxP1_F2xP2_F3xP3<CASE>  (nx, und6, 5, sph9, 3, 1, X0, 5, X1, 2, X2);
    }
    
    // solid, sph -> und, for curl computation (9->9)
    template <CaseFA CASE, class FMat9>
    static void sphericalToUndulated9(const FMat9 &sph9, FMat9 &und9, int nx,
                                      const faN::PropertyN &X0,
                                      const faN::PropertyN &X1,
                                      const faN::PropertyN &X2,
                                      const faN::PropertyN &X3) {
        typedef faN::FieldArithmeticN FA;
        FA::F1xP1_F2xP2<CASE>(nx, und9, 0, sph9, 0, X0, 2, X1);
        FA::F1xP1_F2xP2<CASE>(nx, und9, 1, sph9, 3, X0, 5, X1);
        FA::F1xP1_F2xP2<CASE>(nx, und9, 2, sph9, 6, X0, 8, X1);
        FA::F1xP1_F2xP2<CASE>(nx, und9, 3, sph9, 1, X0, 2, X2);
        FA::F1xP1_F2xP2<CASE>(nx, und9, 4, sph9, 4, X0, 5, X2);
        FA::F1xP1_F2xP2<CASE>(nx, und9, 5, sph9, 7, X0, 8, X2);
        FA::FxP<CASE>        (nx, und9, 6, sph9, 2, X3);
        FA::FxP<CASE>        (nx, und9, 7, sph9, 5, X3);
        FA::FxP<CASE>        (nx, und9, 8, sph9, 8, X3);
    }
    
    // fluid, und -> sph
    template <CaseFA CASE, class FMat3>
    static void undulatedToSpherical3(const FMat3 &und3, FMat3 &sph3, int nx,
                                      const faN::PropertyN &X0,
                                      const faN::PropertyN &X1,
                                      const faN::PropertyN &X2,
                                      const faN::PropertyN &X3) {
        typedef faN::FieldArithmeticN FA;
        FA::FxP<CASE>              (nx, sph3, 0, und3, 0, X0);
        FA::FxP<CASE>              (nx, sph3, 1, und3, 1, X0);
        FA::F1xP1_F2xP2_F3xP3<CASE>(nx, sph3, 2, und3, 0, X1, 1, X2, 2, X3);
    }
    
    // solid, und -> sph
    template <CaseFA CASE, class FMat6, class FMat9>
    static void undulatedToSpherical6(const FMat6 &und6, FMat9 &sph9, int nx,
                                      const faN::PropertyN &X0,
                                      const faN::PropertyN &X1,
                                      const faN::PropertyN &X2,
                                      const faN::PropertyN &X3) {
        typedef faN::FieldArithmeticN FA;
        FA::FxP<CASE>              (nx, sph9, 0, und6, 0, X0);
        FA::FxP<CASE>              (nx, sph9, 1, und6, 5, X0);
        FA::F1xP1_F2xP2_F3xP3<CASE>(nx, sph9, 2, und6, 0, X1, 4, X3, 5, X2);
        FA::F<CASE>                (nx, sph9, 3, sph9, 1);
        FA::FxP<CASE>              (nx, sph9, 4, und6, 1, X0);
        FA::F1xP1_F2xP2_F3xP3<CASE>(nx, sph9, 5, und6, 1, X2, 3, X3, 5, X1);
        FA::FxP<CASE>              (nx, sph9, 6, und6, 4, X0);
        FA::FxP<CASE>              (nx, sph9, 7, und6, 3, X0);
        FA::F1xP1_F2xP2_F3xP3<CASE>(nx, sph9, 8, und6, 2, X3, 3, X2, 4, X1);
    }
};

#endif /* PRT_hpp */
