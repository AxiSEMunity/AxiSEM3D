//
//  ClaytonSolid1D.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 7/29/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  Clayton-Enquist ABC for solid points in 1D
//  theta: angle between surface normal and z-axis

#include "ClaytonSolid1D.hpp"
#include "SolidPoint.hpp"

// apply ABC
void ClaytonSolid1D::apply() const {
    // get fields
    const eigen::CMatX3 &veloc = mSolidPoint->getFields().mVeloc;
    eigen::CMatX3 &stiff = mSolidPoint->getFields().mStiff;
    
    // s, z
    stiff.col(0) -= (mRSA_CosT2_p_RPA_SinT2 * veloc.col(0) +
                     mRPA_m_RSA_x_CosT_SinT * veloc.col(2));
    stiff.col(2) -= (mRPA_m_RSA_x_CosT_SinT * veloc.col(0) +
                     mRSA_SinT2_p_RPA_CosT2 * veloc.col(2));
    
    // phi
    stiff.col(1) -= mRSA * veloc.col(1);
}
