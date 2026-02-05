//
//  GradientQuadrature.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 4/25/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  elemental gradient and quadrature
//  v2.0 features:
//  * template-based, available for both core and preloop
//  * element Jacobian will be considered here instead of in Cijkl

#ifndef GradientQuadrature_hpp
#define GradientQuadrature_hpp

#include "eigen_element.hpp"

// differential operators scaled by integral factor
#ifdef _SAVE_MEMORY
// computed on the fly as static variables
#define xDsDxii_IF sDsDxii_IF
#define xDsDeta_IF sDsDeta_IF
#define xDzDxii_IF sDzDxii_IF
#define xDzDeta_IF sDzDeta_IF
#define x1overS_IF s1overS_IF
#else
// precomputed and stored as member variables
#define xDsDxii_IF mDsDxii_IF
#define xDsDeta_IF mDsDeta_IF
#define xDzDxii_IF mDzDxii_IF
#define xDzDeta_IF mDzDeta_IF
#define x1overS_IF m1overS_IF
#endif

template <typename floatT>
class GradientQuadrature {
    /////////////////////////// typedef ///////////////////////////
    // "using" is not allowed inside a class scope
    static const int P = spectral::nPED;
    
    // complex
    typedef std::complex<floatT> cmplxT;
    inline static const cmplxT I = {0., 1.};
    
    // structured elemental field
    typedef Eigen::Matrix<floatT, P, P, Eigen::RowMajor> RMatPP_RM;
    typedef Eigen::Matrix<cmplxT, P, P, Eigen::RowMajor> CMatPP_RM;
    typedef std::vector<std::array<CMatPP_RM, 1>> vec_ar1_CMatPP_RM;
    typedef std::vector<std::array<CMatPP_RM, 3>> vec_ar3_CMatPP_RM;
    typedef std::vector<std::array<CMatPP_RM, 6>> vec_ar6_CMatPP_RM;
    typedef std::vector<std::array<CMatPP_RM, 9>> vec_ar9_CMatPP_RM;
    
    
    /////////////////////////// data ///////////////////////////
    // operators
    const RMatPP_RM mDsDxii;
    const RMatPP_RM mDsDeta;
    const RMatPP_RM mDzDxii;
    const RMatPP_RM mDzDeta;
    const RMatPP_RM m1overS;
    
    // for axis
    const bool mAxial;
    
#ifdef _SAVE_MEMORY
    // integral factor
    const RMatPP_RM mIntegralFactor;
    
    // operators scaled by integral factor
    inline static RMatPP_RM sDsDxii_IF;
    inline static RMatPP_RM sDsDeta_IF;
    inline static RMatPP_RM sDzDxii_IF;
    inline static RMatPP_RM sDzDeta_IF;
    inline static RMatPP_RM s1overS_IF;
    
    // compute scaled operators on the fly
    void computeScaledOperators() const {
        sDsDxii_IF = mDsDxii.cwiseProduct(mIntegralFactor);
        sDsDeta_IF = mDsDeta.cwiseProduct(mIntegralFactor);
        sDzDxii_IF = mDzDxii.cwiseProduct(mIntegralFactor);
        sDzDeta_IF = mDzDeta.cwiseProduct(mIntegralFactor);
        s1overS_IF = m1overS.cwiseProduct(mIntegralFactor);
    }
#else
    // operators scaled by integral factor
    const RMatPP_RM mDsDxii_IF;
    const RMatPP_RM mDsDeta_IF;
    const RMatPP_RM mDzDxii_IF;
    const RMatPP_RM mDzDeta_IF;
    const RMatPP_RM m1overS_IF;
#endif
    
    
    /////////////////////////// methods ///////////////////////////
public:
#ifdef _SAVE_MEMORY
    // constructor
    GradientQuadrature(const eigen::DMatPP_RM &dsdxii,
                       const eigen::DMatPP_RM &dsdeta,
                       const eigen::DMatPP_RM &dzdxii,
                       const eigen::DMatPP_RM &dzdeta,
                       const eigen::DMatPP_RM &_1ovrS, bool axial,
                       const eigen::DMatPP_RM &integralFactor):
    mDsDxii(dsdxii.cast<floatT>()), mDsDeta(dsdeta.cast<floatT>()),
    mDzDxii(dzdxii.cast<floatT>()), mDzDeta(dzdeta.cast<floatT>()),
    m1overS(_1ovrS.cast<floatT>()), mAxial(axial),
    mIntegralFactor(integralFactor.cast<floatT>()) {
        // nothing
    }
    
    // copy constructor
    GradientQuadrature(const GradientQuadrature &other):
    mDsDxii(other.mDsDxii), mDsDeta(other.mDsDeta),
    mDzDxii(other.mDzDxii), mDzDeta(other.mDzDeta),
    m1overS(other.m1overS), mAxial(other.mAxial),
    mIntegralFactor(other.mIntegralFactor) {
        // nothing
    }
#else
    // constructor
    GradientQuadrature(const eigen::DMatPP_RM &dsdxii,
                       const eigen::DMatPP_RM &dsdeta,
                       const eigen::DMatPP_RM &dzdxii,
                       const eigen::DMatPP_RM &dzdeta,
                       const eigen::DMatPP_RM &_1ovrS, bool axial,
                       const eigen::DMatPP_RM &integralFactor):
    mDsDxii(dsdxii.cast<floatT>()), mDsDeta(dsdeta.cast<floatT>()),
    mDzDxii(dzdxii.cast<floatT>()), mDzDeta(dzdeta.cast<floatT>()),
    m1overS(_1ovrS.cast<floatT>()), mAxial(axial),
    mDsDxii_IF((dsdxii.cwiseProduct(integralFactor)).cast<floatT>()),
    mDsDeta_IF((dsdeta.cwiseProduct(integralFactor)).cast<floatT>()),
    mDzDxii_IF((dzdxii.cwiseProduct(integralFactor)).cast<floatT>()),
    mDzDeta_IF((dzdeta.cwiseProduct(integralFactor)).cast<floatT>()),
    m1overS_IF((_1ovrS.cwiseProduct(integralFactor)).cast<floatT>()) {
        // nothing
    }
    
    // copy constructor
    GradientQuadrature(const GradientQuadrature &other):
    mDsDxii(other.mDsDxii), mDsDeta(other.mDsDeta),
    mDzDxii(other.mDzDxii), mDzDeta(other.mDzDeta),
    m1overS(other.m1overS), mAxial(other.mAxial),
    mDsDxii_IF(other.mDsDxii_IF),
    mDsDeta_IF(other.mDsDeta_IF),
    mDzDxii_IF(other.mDzDxii_IF),
    mDzDeta_IF(other.mDzDeta_IF),
    m1overS_IF(other.m1overS_IF) {
        // nothing
    }
#endif
    
    // check compatibility
    void checkCompatibility(int nr) const {
        // expand workspace if needed
        int nu_1 = nr / 2 + 1;
        if (sA.size() < nu_1) {
            sA.resize(nu_1);
        }
    }
    
    // grad 1 -> 3
    void computeGrad3(const vec_ar1_CMatPP_RM &p, vec_ar3_CMatPP_RM &p_i,
                      int nu_1) const {
        ////// higher-order terms //////
        computeGrad_dSdZ(nu_1, p, 0, p_i, 0, 2);
        
        ////// lower-order terms //////
        // i * alpha * p
        for (int alpha = 0; alpha < nu_1; alpha++) {
            cmplxT alpha_I = (floatT)alpha * I;
            sA[alpha][0] = alpha_I * p[alpha][0];
        }
        computeGrad_1overS_setTo(nu_1, sA, 0, p_i, 1);
    }
    
    // grad 3 -> 9
    // 00->0 01->1 02->2
    // 10->3 11->4 12->5
    // 20->6 21->7 22->8
    void computeGrad9(const vec_ar3_CMatPP_RM &ui, vec_ar9_CMatPP_RM &ui_j,
                      int nu_1) const {
        ////// higher-order terms //////
        computeGrad_dSdZ(nu_1, ui, 0, ui_j, 0, 2);
        computeGrad_dSdZ(nu_1, ui, 1, ui_j, 3, 5);
        computeGrad_dSdZ(nu_1, ui, 2, ui_j, 6, 8);
        
        ////// lower-order terms //////
        // i * alpha * u0 - u1 -> 1
        // i * alpha * u1 + u0 -> 4
        // i * alpha * u2 -> 7
        for (int alpha = 0; alpha < nu_1; alpha++) {
            cmplxT alpha_I = (floatT)alpha * I;
            sA[alpha][0] = alpha_I * ui[alpha][0] - ui[alpha][1];
            sA[alpha][1] = alpha_I * ui[alpha][1] + ui[alpha][0];
            sA[alpha][2] = alpha_I * ui[alpha][2];
        }
        computeGrad_1overS_setTo(nu_1, sA, 0, ui_j, 1);
        computeGrad_1overS_setTo(nu_1, sA, 1, ui_j, 4);
        computeGrad_1overS_setTo(nu_1, sA, 2, ui_j, 7);
    }
    
    // grad 3 -> 6
    // 00->0 01->5 02->4
    // 10->5 11->1 12->3
    // 20->4 21->3 22->2
    void computeGrad6(const vec_ar3_CMatPP_RM &ui, vec_ar6_CMatPP_RM &eij,
                      int nu_1) const {
        ////// higher-order terms //////
        computeGrad_dSdZ(nu_1, ui, 0, eij, 0, 4);
        computeGrad_dSdZ(nu_1, ui, 1, eij, 5, 3);
        computeGrad_dSdZ(nu_1, ui, 2, eij, 1, 2);
        // dim 1 used as a temporary block
        for (int alpha = 0; alpha < nu_1; alpha++) {
            eij[alpha][4] += eij[alpha][1];
        }
        
        ////// lower-order terms //////
        // i * alpha * u0 - u1 -> 5
        // i * alpha * u1 + u0 -> 1
        // i * alpha * u2 -> 3
        for (int alpha = 0; alpha < nu_1; alpha++) {
            cmplxT alpha_I = (floatT)alpha * I;
            sA[alpha][0] = alpha_I * ui[alpha][0] - ui[alpha][1];
            sA[alpha][1] = alpha_I * ui[alpha][1] + ui[alpha][0];
            sA[alpha][2] = alpha_I * ui[alpha][2];
        }
        computeGrad_1overS_addTo(nu_1, sA, 0, eij, 5);
        computeGrad_1overS_setTo(nu_1, sA, 1, eij, 1);
        computeGrad_1overS_addTo(nu_1, sA, 2, eij, 3);
    }
    
    // quad 3 -> 1
    void computeQuad3(const vec_ar3_CMatPP_RM &q_i, vec_ar1_CMatPP_RM &q,
                      int nu_1) const {
#ifdef _SAVE_MEMORY
        // compute scaled operators
        computeScaledOperators();
#endif
        ////// higher-order terms //////
        computeQuad_dSdZ(nu_1, q_i, 0, 2, q, 0,
                         xDsDxii_IF, xDsDeta_IF, xDzDxii_IF, xDzDeta_IF);
        
        ////// lower-order terms //////
        // i * alpha * q1
        for (int alpha = 0; alpha < nu_1; alpha++) {
            cmplxT alpha_I = (floatT)alpha * I;
            sA[alpha][0] = alpha_I * q_i[alpha][1];
        }
        computeQuad_1overS_addTo(nu_1, sA, 0, q, 0,
                                 x1overS_IF, xDzDxii_IF, xDzDeta_IF);
    }
    
    // quad 9 -> 3
    // 00->0 01->1 02->2
    // 10->3 11->4 12->5
    // 20->6 21->7 22->8
    void computeQuad9(const vec_ar9_CMatPP_RM &fi_j, vec_ar3_CMatPP_RM &fi,
                      int nu_1) const {
#ifdef _SAVE_MEMORY
        // compute scaled operators
        computeScaledOperators();
#endif
        ////// higher-order terms //////
        computeQuad_dSdZ(nu_1, fi_j, 0, 2, fi, 0,
                         xDsDxii_IF, xDsDeta_IF, xDzDxii_IF, xDzDeta_IF);
        computeQuad_dSdZ(nu_1, fi_j, 3, 5, fi, 1,
                         xDsDxii_IF, xDsDeta_IF, xDzDxii_IF, xDzDeta_IF);
        computeQuad_dSdZ(nu_1, fi_j, 6, 8, fi, 2,
                         xDsDxii_IF, xDsDeta_IF, xDzDxii_IF, xDzDeta_IF);
        
        ////// lower-order terms //////
        // i * alpha * df1 - df4 -> 0
        // i * alpha * df4 + df1 -> 1
        // i * alpha * df7 -> 2
        for (int alpha = 0; alpha < nu_1; alpha++) {
            cmplxT alpha_I = (floatT)alpha * I;
            sA[alpha][0] = alpha_I * fi_j[alpha][1] - fi_j[alpha][4];
            sA[alpha][1] = alpha_I * fi_j[alpha][4] + fi_j[alpha][1];
            sA[alpha][2] = alpha_I * fi_j[alpha][7];
        }
        computeQuad_1overS_addTo(nu_1, sA, 0, fi, 0,
                                 x1overS_IF, xDzDxii_IF, xDzDeta_IF);
        computeQuad_1overS_addTo(nu_1, sA, 1, fi, 1,
                                 x1overS_IF, xDzDxii_IF, xDzDeta_IF);
        computeQuad_1overS_addTo(nu_1, sA, 2, fi, 2,
                                 x1overS_IF, xDzDxii_IF, xDzDeta_IF);
    }
    
    // quad 6 -> 3
    // 00->0 01->5 02->4
    // 10->5 11->1 12->3
    // 20->4 21->3 22->2
    void computeQuad6(const vec_ar6_CMatPP_RM &sij, vec_ar3_CMatPP_RM &fi,
                      int nu_1) const {
#ifdef _SAVE_MEMORY
        // compute scaled operators
        computeScaledOperators();
#endif
        ////// higher-order terms //////
        computeQuad_dSdZ(nu_1, sij, 0, 4, fi, 0,
                         xDsDxii_IF, xDsDeta_IF, xDzDxii_IF, xDzDeta_IF);
        computeQuad_dSdZ(nu_1, sij, 5, 3, fi, 1,
                         xDsDxii_IF, xDsDeta_IF, xDzDxii_IF, xDzDeta_IF);
        computeQuad_dSdZ(nu_1, sij, 4, 2, fi, 2,
                         xDsDxii_IF, xDsDeta_IF, xDzDxii_IF, xDzDeta_IF);
        
        ////// lower-order terms //////
        // i * alpha * df5 - df1 -> 0
        // i * alpha * df1 + df5 -> 1
        // i * alpha * df3 -> 2
        for (int alpha = 0; alpha < nu_1; alpha++) {
            cmplxT alpha_I = (floatT)alpha * I;
            sA[alpha][0] = alpha_I * sij[alpha][5] - sij[alpha][1];
            sA[alpha][1] = alpha_I * sij[alpha][1] + sij[alpha][5];
            sA[alpha][2] = alpha_I * sij[alpha][3];
        }
        computeQuad_1overS_addTo(nu_1, sA, 0, fi, 0,
                                 x1overS_IF, xDzDxii_IF, xDzDeta_IF);
        computeQuad_1overS_addTo(nu_1, sA, 1, fi, 1,
                                 x1overS_IF, xDzDxii_IF, xDzDeta_IF);
        computeQuad_1overS_addTo(nu_1, sA, 2, fi, 2,
                                 x1overS_IF, xDzDxii_IF, xDzDeta_IF);
    }
    
    // quad 9 -> 3, no integration, for moment tensor
    // 00->0 01->1 02->2
    // 10->3 11->4 12->5
    // 20->6 21->7 22->8
    void computeQuad9_NoIntegration(const vec_ar9_CMatPP_RM &fi_j,
                                    vec_ar3_CMatPP_RM &fi,
                                    int nu_1) const {
        ////// higher-order terms //////
        computeQuad_dSdZ(nu_1, fi_j, 0, 2, fi, 0,
                         mDsDxii, mDsDeta, mDzDxii, mDzDeta);
        computeQuad_dSdZ(nu_1, fi_j, 3, 5, fi, 1,
                         mDsDxii, mDsDeta, mDzDxii, mDzDeta);
        computeQuad_dSdZ(nu_1, fi_j, 6, 8, fi, 2,
                         mDsDxii, mDsDeta, mDzDxii, mDzDeta);
        
        ////// lower-order terms //////
        // i * alpha * df1 - df4 -> 0
        // i * alpha * df4 + df1 -> 1
        // i * alpha * df7 -> 2
        for (int alpha = 0; alpha < nu_1; alpha++) {
            cmplxT alpha_I = (floatT)alpha * I;
            sA[alpha][0] = alpha_I * fi_j[alpha][1] - fi_j[alpha][4];
            sA[alpha][1] = alpha_I * fi_j[alpha][4] + fi_j[alpha][1];
            sA[alpha][2] = alpha_I * fi_j[alpha][7];
        }
        computeQuad_1overS_addTo(nu_1, sA, 0, fi, 0,
                                 m1overS, mDzDxii, mDzDeta);
        computeQuad_1overS_addTo(nu_1, sA, 1, fi, 1,
                                 m1overS, mDzDxii, mDzDeta);
        computeQuad_1overS_addTo(nu_1, sA, 2, fi, 2,
                                 m1overS, mDzDxii, mDzDeta);
    }
    
    // quad 6 -> 3, no integration, for moment tensor
    // 00->0 01->5 02->4
    // 10->5 11->1 12->3
    // 20->4 21->3 22->2
    void computeQuad6_NoIntegration(const vec_ar6_CMatPP_RM &sij,
                                    vec_ar3_CMatPP_RM &fi,
                                    int nu_1) const {
        ////// higher-order terms //////
        computeQuad_dSdZ(nu_1, sij, 0, 4, fi, 0,
                         mDsDxii, mDsDeta, mDzDxii, mDzDeta);
        computeQuad_dSdZ(nu_1, sij, 5, 3, fi, 1,
                         mDsDxii, mDsDeta, mDzDxii, mDzDeta);
        computeQuad_dSdZ(nu_1, sij, 4, 2, fi, 2,
                         mDsDxii, mDsDeta, mDzDxii, mDzDeta);
        
        ////// lower-order terms //////
        // i * alpha * df5 - df1 -> 0
        // i * alpha * df1 + df5 -> 1
        // i * alpha * df3 -> 2
        for (int alpha = 0; alpha < nu_1; alpha++) {
            cmplxT alpha_I = (floatT)alpha * I;
            sA[alpha][0] = alpha_I * sij[alpha][5] - sij[alpha][1];
            sA[alpha][1] = alpha_I * sij[alpha][1] + sij[alpha][5];
            sA[alpha][2] = alpha_I * sij[alpha][3];
        }
        computeQuad_1overS_addTo(nu_1, sA, 0, fi, 0,
                                 m1overS, mDzDxii, mDzDeta);
        computeQuad_1overS_addTo(nu_1, sA, 1, fi, 1,
                                 m1overS, mDzDxii, mDzDeta);
        computeQuad_1overS_addTo(nu_1, sA, 2, fi, 2,
                                 m1overS, mDzDxii, mDzDeta);
    }
    
    // get axial flag for cost signature
    bool axial() const {
        return mAxial;
    }
    
    
    ///////////////////////////////////////////////////////////////////////////
    ///////////////////////// template tool functions /////////////////////////
    ///////////////////////////////////////////////////////////////////////////
    
private:
    // higher-order terms, dF/dS and dF/dZ
    template <class VecArPP_F, class VecArPP_dF>
    void computeGrad_dSdZ(int nu_1, const VecArPP_F &F, int dimF,
                          VecArPP_dF &dF, int dim_dFdS, int dim_dFdZ) const {
        // GT_xii
        const RMatPP_RM &GT_xii = mAxial ? sGT_GLJ : sGT_GLL;
        // loop over alpha
        for (int alpha = 0; alpha < nu_1; alpha++) {
            sX.noalias() = GT_xii * F[alpha][dimF];
            sY.noalias() = F[alpha][dimF] * sG_GLL;
            dF[alpha][dim_dFdS] = (sX.cwiseProduct(mDzDeta) +
                                   sY.cwiseProduct(mDzDxii));
            dF[alpha][dim_dFdZ] = (sX.cwiseProduct(mDsDeta) +
                                   sY.cwiseProduct(mDsDxii));
        }
    }
    
    // lower-order terms, A/S
    template <class VecArPP_A, class VecArPP_AoS>
    void computeGrad_1overS_setTo(int nu_1, const VecArPP_A &A, int dimA,
                                  VecArPP_AoS &AoS, int dimAoS) const {
        // GT_xii
        const RMatPP_RM &GT_xii = mAxial ? sGT_GLJ : sGT_GLL;
        // loop over alpha
        for (int alpha = 0; alpha < nu_1; alpha++) {
            AoS[alpha][dimAoS] = A[alpha][dimA].cwiseProduct(m1overS);
            // axial
            if (mAxial) {
                sX.row(0).noalias() = GT_xii.row(0) * A[alpha][dimA];
                sY.row(0).noalias() = A[alpha][dimA].row(0) * sG_GLL;
                AoS[alpha][dimAoS].row(0) =
                (sX.row(0).cwiseProduct(mDzDeta.row(0)) +
                 sY.row(0).cwiseProduct(mDzDxii.row(0)));
            }
        }
    }
    
    // lower-order terms, A/S
    // same as above but add to the result
    template <class VecArPP_A, class VecArPP_AoS>
    void computeGrad_1overS_addTo(int nu_1, const VecArPP_A &A, int dimA,
                                  VecArPP_AoS &AoS, int dimAoS) const {
        // GT_xii
        const RMatPP_RM &GT_xii = mAxial ? sGT_GLJ : sGT_GLL;
        // loop over alpha
        for (int alpha = 0; alpha < nu_1; alpha++) {
            AoS[alpha][dimAoS] += A[alpha][dimA].cwiseProduct(m1overS);
            // axial
            if (mAxial) {
                sX.row(0).noalias() = GT_xii.row(0) * A[alpha][dimA];
                sY.row(0).noalias() = A[alpha][dimA].row(0) * sG_GLL;
                AoS[alpha][dimAoS].row(0) +=
                (sX.row(0).cwiseProduct(mDzDeta.row(0)) +
                 sY.row(0).cwiseProduct(mDzDxii.row(0)));
            }
        }
    }
    
    // higher-order terms, dF/dS and dF/dZ
    template <class VecArPP_dF, class VecArPP_F>
    void computeQuad_dSdZ(int nu_1, const VecArPP_dF &dF, int dim_dFdS,
                          int dim_dFdZ, VecArPP_F &F, int dimF,
                          const RMatPP_RM &dsdxii, const RMatPP_RM &dsdeta,
                          const RMatPP_RM &dzdxii, const RMatPP_RM &dzdeta)
    const {
        // G_xii
        const RMatPP_RM &G_xii = mAxial ? sG_GLJ : sG_GLL;
        // loop over alpha
        for (int alpha = 0; alpha < nu_1; alpha++) {
            sX = (dF[alpha][dim_dFdS].cwiseProduct(dzdeta) +
                  dF[alpha][dim_dFdZ].cwiseProduct(dsdeta));
            sY = (dF[alpha][dim_dFdS].cwiseProduct(dzdxii) +
                  dF[alpha][dim_dFdZ].cwiseProduct(dsdxii));
            F[alpha][dimF].noalias() = G_xii * sX + sY * sGT_GLL;
        }
    }
    
    // lower-order terms, A/S
    template <class VecArPP_A, class VecArPP_AoS>
    void computeQuad_1overS_addTo(int nu_1, const VecArPP_A &A, int dimA,
                                  VecArPP_AoS &AoS, int dimAoS,
                                  const RMatPP_RM &_1overS,
                                  const RMatPP_RM &dzdxii,
                                  const RMatPP_RM &dzdeta) const {
        // G_xii
        const RMatPP_RM &G_xii = mAxial ? sG_GLJ : sG_GLL;
        // loop over alpha
        for (int alpha = 0; alpha < nu_1; alpha++) {
            // A/S
            AoS[alpha][dimAoS] -= A[alpha][dimA].cwiseProduct(_1overS);
            if (mAxial) {
                sX.row(0) = A[alpha][dimA].row(0).cwiseProduct(dzdeta.row(0));
                sY.row(0) = A[alpha][dimA].row(0).cwiseProduct(dzdxii.row(0));
                AoS[alpha][dimAoS].noalias() -= G_xii.col(0) * sX.row(0);
                AoS[alpha][dimAoS].row(0).noalias() -= sY.row(0) * sGT_GLL;
            }
        }
    }
    
    
    ////////////////////////////////////////
    //////////////// static ////////////////
    ////////////////////////////////////////
    
    ///////////// G Mat /////////////
public:
    // set G Mat
    static void setGMat(const eigen::DMatPP_RM &G_GLL,
                        const eigen::DMatPP_RM &G_GLJ) {
        sG_GLL = G_GLL.cast<floatT>();
        sG_GLJ = G_GLJ.cast<floatT>();
        sGT_GLL = G_GLL.transpose().cast<floatT>();
        sGT_GLJ = G_GLJ.transpose().cast<floatT>();
    }
    
private:
    // G for GLL and GLJ
    inline static RMatPP_RM sG_GLL;
    inline static RMatPP_RM sG_GLJ;
    inline static RMatPP_RM sGT_GLL;
    inline static RMatPP_RM sGT_GLJ;
    
    
    ///////////// workspace /////////////
    inline static CMatPP_RM sX;
    inline static CMatPP_RM sY;
    inline static vec_ar3_CMatPP_RM sA;
};

#endif /* GradientQuadrature_hpp */
