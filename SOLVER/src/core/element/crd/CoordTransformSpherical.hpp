//
//  CoordTransformSpherical.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/10/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  coordinate transform between (s,phi,z) and (R,T,Z)
//  for vector and 2nd-order tensor fields

#ifndef CoordTransformSpherical_hpp
#define CoordTransformSpherical_hpp

#include "CoordTransform.hpp"

// sin(theta)
#ifdef _SAVE_MEMORY
// computed on the fly as static variables
#define xSin1t sSin1t
#define xCos1t sCos1t
#define xSin2t sSin2t
#define xCos2t sCos2t
#else
// precomputed and stored as member variables
#define xSin1t mSin1t
#define xCos1t mCos1t
#define xSin2t mSin2t
#define xCos2t mCos2t
#endif

class CoordTransformSpherical: public CoordTransform {
    inline static const numerical::Real half = .5;
    
public:
#ifdef _SAVE_MEMORY
    // constructor
    CoordTransformSpherical(const eigen::DMatPP_RM &theta):
    mTheta(theta.cast<numerical::Real>()) {
        // nothing
    }
#else
    // constructor
    CoordTransformSpherical(const eigen::DMatPP_RM &theta):
    mSin1t(theta.array().sin().cast<numerical::Real>()),
    mCos1t(theta.array().cos().cast<numerical::Real>()),
    mSin2t((2. * theta).array().sin().cast<numerical::Real>()),
    mCos2t((2. * theta).array().cos().cast<numerical::Real>()) {
        // nothing
    }
#endif
    
    // (s,phi,z) -> (R,T,Z)
    void transformSPZ_RTZ3(eigen::vec_ar3_CMatPP_RM &ui, int nu_1) const {
#ifdef _SAVE_MEMORY
        // compute sin(theta)
        computeSinOneTheta();
#endif
        for (int alpha = 0; alpha < nu_1; alpha++) {
            // copy u0
            eigen::CMatPP_RM &ui_alpha__0_ = sTemp0;
            ui_alpha__0_ = ui[alpha][0];
            
            // rtz0 = spz0 * cos(t) - spz2 * sin(t)
            ui[alpha][0] = (ui_alpha__0_.cwiseProduct(xCos1t) -
                            ui[alpha][2].cwiseProduct(xSin1t));
            
            // rtz2 = spz2 * cos(t) + spz0 * sin(t)
            ui[alpha][2] = (ui[alpha][2].cwiseProduct(xCos1t) +
                            ui_alpha__0_.cwiseProduct(xSin1t));
            
            // rtz1 = spz1
            // nothing
        }
    }
    
    // (R,T,Z) -> (s,phi,z)
    void transformRTZ_SPZ3(eigen::vec_ar3_CMatPP_RM &ui, int nu_1) const {
#ifdef _SAVE_MEMORY
        // compute sin(theta)
        computeSinOneTheta();
#endif
        for (int alpha = 0; alpha < nu_1; alpha++) {
            // copy u0
            eigen::CMatPP_RM &ui_alpha__0_ = sTemp0;
            ui_alpha__0_ = ui[alpha][0];
            
            // spz0 = rtz0 * cos(t) + rtz2 * sin(t)
            ui[alpha][0] = (ui_alpha__0_.cwiseProduct(xCos1t) +
                            ui[alpha][2].cwiseProduct(xSin1t));
            
            // spz2 = rtz2 * cos(t) - rtz0 * sin(t)
            ui[alpha][2] = (ui[alpha][2].cwiseProduct(xCos1t) -
                            ui_alpha__0_.cwiseProduct(xSin1t));
            
            // rtz1 = spz1
            // nothing
        }
    }
    
    // (s,phi,z) -> (R,T,Z) for nabla
    void transformSPZ_RTZ9(eigen::vec_ar9_CMatPP_RM &nij, int nu_1) const {
#ifdef _SAVE_MEMORY
        // compute sin(theta)
        computeSinOneTwoTheta();
#endif
        for (int alpha = 0; alpha < nu_1; alpha++) {
            // 2 * theta terms
            eigen::CMatPP_RM &s08 = sTemp0;
            eigen::CMatPP_RM &d08 = sTemp1;
            eigen::CMatPP_RM &s26 = sTemp2;
            eigen::CMatPP_RM &d26 = sTemp3;
            s08 = nij[alpha][0] + nij[alpha][8];
            d08 = nij[alpha][0] - nij[alpha][8];
            s26 = nij[alpha][2] + nij[alpha][6];
            d26 = nij[alpha][2] - nij[alpha][6];
            nij[alpha][0] = (s08 + d08.cwiseProduct(xCos2t) -
                             s26.cwiseProduct(xSin2t)) * half;
            nij[alpha][2] = (d26 + s26.cwiseProduct(xCos2t) +
                             d08.cwiseProduct(xSin2t)) * half;
            nij[alpha][6] = nij[alpha][2] - d26;
            nij[alpha][8] = s08 - nij[alpha][0];
            
            // 1 * theta terms
            eigen::CMatPP_RM &nij_alpha__1_ = sTemp0;
            nij_alpha__1_ = nij[alpha][1];
            nij[alpha][1] = (nij_alpha__1_.cwiseProduct(xCos1t) -
                             nij[alpha][7].cwiseProduct(xSin1t));
            nij[alpha][7] = (nij[alpha][7].cwiseProduct(xCos1t) +
                             nij_alpha__1_.cwiseProduct(xSin1t));
            eigen::CMatPP_RM &nij_alpha__3_ = sTemp0;
            nij_alpha__3_ = nij[alpha][3];
            nij[alpha][3] = (nij_alpha__3_.cwiseProduct(xCos1t) -
                             nij[alpha][5].cwiseProduct(xSin1t));
            nij[alpha][5] = (nij[alpha][5].cwiseProduct(xCos1t) +
                             nij_alpha__3_.cwiseProduct(xSin1t));
            
            // 0 * theta terms
            // nothing
        }
    }
    
    // (R,T,Z) -> (s,phi,z) for nabla
    void transformRTZ_SPZ9(eigen::vec_ar9_CMatPP_RM &nij, int nu_1) const {
#ifdef _SAVE_MEMORY
        // compute sin(theta)
        computeSinOneTwoTheta();
#endif
        for (int alpha = 0; alpha < nu_1; alpha++) {
            // 2 * theta terms
            eigen::CMatPP_RM &s08 = sTemp0;
            eigen::CMatPP_RM &d08 = sTemp1;
            eigen::CMatPP_RM &s26 = sTemp2;
            eigen::CMatPP_RM &d26 = sTemp3;
            s08 = nij[alpha][0] + nij[alpha][8];
            d08 = nij[alpha][0] - nij[alpha][8];
            s26 = nij[alpha][2] + nij[alpha][6];
            d26 = nij[alpha][2] - nij[alpha][6];
            nij[alpha][0] = (s08 + d08.cwiseProduct(xCos2t) +
                             s26.cwiseProduct(xSin2t)) * half;
            nij[alpha][2] = (d26 + s26.cwiseProduct(xCos2t) -
                             d08.cwiseProduct(xSin2t)) * half;
            nij[alpha][6] = nij[alpha][2] - d26;
            nij[alpha][8] = s08 - nij[alpha][0];
            
            // 1 * theta terms
            eigen::CMatPP_RM &nij_alpha__1_ = sTemp0;
            nij_alpha__1_ = nij[alpha][1];
            nij[alpha][1] = (nij_alpha__1_.cwiseProduct(xCos1t) +
                             nij[alpha][7].cwiseProduct(xSin1t));
            nij[alpha][7] = (nij[alpha][7].cwiseProduct(xCos1t) -
                             nij_alpha__1_.cwiseProduct(xSin1t));
            eigen::CMatPP_RM &nij_alpha__3_ = sTemp0;
            nij_alpha__3_ = nij[alpha][3];
            nij[alpha][3] = (nij_alpha__3_.cwiseProduct(xCos1t) +
                             nij[alpha][5].cwiseProduct(xSin1t));
            nij[alpha][5] = (nij[alpha][5].cwiseProduct(xCos1t) -
                             nij_alpha__3_.cwiseProduct(xSin1t));
            
            // 0 * theta terms
            // nothing
        }
    }
    
    // (s,phi,z) -> (R,T,Z) for Voigt
    void transformSPZ_RTZ6(eigen::vec_ar6_CMatPP_RM &eij, int nu_1) const {
#ifdef _SAVE_MEMORY
        // compute sin(theta)
        computeSinOneTwoTheta();
#endif
        for (int alpha = 0; alpha < nu_1; alpha++) {
            // 2 * theta terms
            eigen::CMatPP_RM &s02 = sTemp0;
            eigen::CMatPP_RM &d02 = sTemp1;
            s02 = eij[alpha][0] + eij[alpha][2];
            d02 = eij[alpha][0] - eij[alpha][2];
            eij[alpha][0] = (s02 + d02.cwiseProduct(xCos2t) -
                             eij[alpha][4].cwiseProduct(xSin2t)) * half;
            eij[alpha][2] = s02 - eij[alpha][0];
            eij[alpha][4] = (eij[alpha][4].cwiseProduct(xCos2t) +
                             d02.cwiseProduct(xSin2t));
            
            // 1 * theta terms
            eigen::CMatPP_RM &nij_alpha__3_ = sTemp0;
            nij_alpha__3_ = eij[alpha][3];
            eij[alpha][3] = (nij_alpha__3_.cwiseProduct(xCos1t) +
                             eij[alpha][5].cwiseProduct(xSin1t));
            eij[alpha][5] = (eij[alpha][5].cwiseProduct(xCos1t) -
                             nij_alpha__3_.cwiseProduct(xSin1t));
            
            // 0 * theta terms
            // nothing
        }
    }
    
    // (R,T,Z) -> (s,phi,z) for Voigt
    void transformRTZ_SPZ6(eigen::vec_ar6_CMatPP_RM &sij, int nu_1) const {
#ifdef _SAVE_MEMORY
        // compute sin(theta)
        computeSinOneTwoTheta();
#endif
        for (int alpha = 0; alpha < nu_1; alpha++) {
            // 2 * theta terms
            eigen::CMatPP_RM &s02 = sTemp0;
            eigen::CMatPP_RM &d02 = sTemp1;
            s02 = sij[alpha][0] + sij[alpha][2];
            d02 = sij[alpha][0] - sij[alpha][2];
            sij[alpha][0] = ((s02 + d02.cwiseProduct(xCos2t)) * half +
                             sij[alpha][4].cwiseProduct(xSin2t));
            sij[alpha][2] = s02 - sij[alpha][0];
            sij[alpha][4] = (sij[alpha][4].cwiseProduct(xCos2t) -
                             d02.cwiseProduct(xSin2t) * half);
            
            // 1 * theta terms
            eigen::CMatPP_RM &nij_alpha__3_ = sTemp0;
            nij_alpha__3_ = sij[alpha][3];
            sij[alpha][3] = (nij_alpha__3_.cwiseProduct(xCos1t) -
                             sij[alpha][5].cwiseProduct(xSin1t));
            sij[alpha][5] = (sij[alpha][5].cwiseProduct(xCos1t) +
                             nij_alpha__3_.cwiseProduct(xSin1t));
            
            // 0 * theta terms
            // nothing
        }
    }
    
private:
#ifdef _SAVE_MEMORY
    // theta
    const eigen::RMatPP_RM mTheta;
    
    // sin(theta)
    inline static eigen::RMatPP_RM sSin1t, sCos1t, sSin2t, sCos2t;
    
    // compute sin(theta) on the fly
    void computeSinOneTheta() const {
        sSin1t = mTheta.array().sin();
        sCos1t = mTheta.array().cos();
    }
    
    // compute sin(theta) and sin(2 theta) on the fly
    void computeSinOneTwoTheta() const {
        computeSinOneTheta();
        sSin2t = (2. * mTheta).array().sin();
        sCos2t = (2. * mTheta).array().cos();
    }
#else
    // sin(theta)
    const eigen::RMatPP_RM mSin1t, mCos1t, mSin2t, mCos2t;
#endif
    
    
    ////////////////////////////////////////
    //////////////// static ////////////////
    ////////////////////////////////////////
    
private:
    // workspace
    inline static eigen::CMatPP_RM sTemp0, sTemp1, sTemp2, sTemp3;
};

#endif /* CoordTransformSpherical_hpp */
