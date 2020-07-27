//
//  Station.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 8/4/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  station

#include "Station.hpp"

// set element: inplane weights and Fourier exp
void Station::setElement(const eigen::DRowN &weights, int nu_1) {
    // indices of non-zero weights
    mNonZeroIndices.clear();
    double rEpsilon = (double)numerical::epsilon<numerical::Real>();
    for (int ipnt = 0; ipnt < spectral::nPEM; ipnt++) {
        if (std::abs(weights(ipnt)) >= rEpsilon) {
            mNonZeroIndices.push_back(ipnt);
        }
    }
    mNonZeroIndices.shrink_to_fit();
    
    // non-zero weights
    mNonZeroWeights = weights(mNonZeroIndices).cast<numerical::Real>();
    
    // 2 * exp(i * alpha * phi) for Fourier interpolation
#ifndef _SAVE_MEMORY
    // precompute
    m2ExpIAlphaPhi.resize(nu_1);
    eigen_tools::computeTwoExpIAlphaPhi(nu_1, mPhi, m2ExpIAlphaPhi);
#else
    // on-the-fly
    if (s2ExpIAlphaPhi.rows() < nu_1) {
        s2ExpIAlphaPhi.resize(nu_1);
    }
#endif
}

// rotate vector
void Station::rotate(eigen::RMatX3_RM &V, int nrow, int k, double angle) {
    // frame index
    int i = (k + 1 > 2) ? k - 2 : k + 1;
    int j = (k + 2 > 2) ? k - 1 : k + 2;
    // angle
    numerical::Real cosA = (numerical::Real)cos(angle);
    numerical::Real sinA = (numerical::Real)sin(angle);
    // rotate by k-axis
    sColBuffer.topRows(nrow) = V.block(0, i, nrow, 1);
    V.block(0, i, nrow, 1) = (sColBuffer.topRows(nrow) * cosA +
                              V.block(0, j, nrow, 1) * sinA);
    V.block(0, j, nrow, 1) = (-sColBuffer.topRows(nrow) * sinA +
                              V.block(0, j, nrow, 1) * cosA);
}

// rotate tensor
void Station::rotate(eigen::RMatX9_RM &T, int nrow, int k, double angle) {
    // frame index
    int i = (k + 1 > 2) ? k - 2 : k + 1;
    int j = (k + 2 > 2) ? k - 1 : k + 2;
    // angle
    numerical::Real cosA = (numerical::Real)cos(angle);
    numerical::Real sinA = (numerical::Real)sin(angle);
    // Smp = Qmn Tnp
    for (int n = 0; n < 3; n++) {
        int I = i * 3 + n;
        int J = j * 3 + n;
        sColBuffer.topRows(nrow) = T.block(0, I, nrow, 1);
        T.block(0, I, nrow, 1) = (sColBuffer.topRows(nrow) * cosA +
                                  T.block(0, J, nrow, 1) * sinA);
        T.block(0, J, nrow, 1) = (-sColBuffer.topRows(nrow) * sinA +
                                  T.block(0, J, nrow, 1) * cosA);
    }
    // Tmn = Smp Qnp
    for (int p = 0; p < 3; p++) {
        int I = p * 3 + i;
        int J = p * 3 + j;
        sColBuffer.topRows(nrow) = T.block(0, I, nrow, 1);
        T.block(0, I, nrow, 1) = (sColBuffer.topRows(nrow) * cosA +
                                  T.block(0, J, nrow, 1) * sinA);
        T.block(0, J, nrow, 1) = (-sColBuffer.topRows(nrow) * sinA +
                                  T.block(0, J, nrow, 1) * cosA);
    }
}

// rotate tensor in Voigt
// 0(0) 1(5) 2(4)
// 3(5) 4(1) 5(3)
// 6(4) 7(3) 8(2)
void Station::rotate(eigen::RMatX6_RM &S, int nrow, int k, double angle) {
    // voigt -> 3*3
    sTensor33.block(0, 0, nrow, 1) = S.block(0, 0, nrow, 1);
    sTensor33.block(0, 1, nrow, 1) = S.block(0, 5, nrow, 1);
    sTensor33.block(0, 2, nrow, 1) = S.block(0, 4, nrow, 1);
    sTensor33.block(0, 3, nrow, 1) = S.block(0, 5, nrow, 1);
    sTensor33.block(0, 4, nrow, 1) = S.block(0, 1, nrow, 1);
    sTensor33.block(0, 5, nrow, 1) = S.block(0, 3, nrow, 1);
    sTensor33.block(0, 6, nrow, 1) = S.block(0, 4, nrow, 1);
    sTensor33.block(0, 7, nrow, 1) = S.block(0, 3, nrow, 1);
    sTensor33.block(0, 8, nrow, 1) = S.block(0, 2, nrow, 1);
    // rotate 3*3
    rotate(sTensor33, nrow, k, angle);
    // 3*3 -> voigt
    S.block(0, 0, nrow, 1) = sTensor33.block(0, 0, nrow, 1);
    S.block(0, 1, nrow, 1) = sTensor33.block(0, 4, nrow, 1);
    S.block(0, 2, nrow, 1) = sTensor33.block(0, 8, nrow, 1);
    S.block(0, 3, nrow, 1) = sTensor33.block(0, 5, nrow, 1);
    S.block(0, 4, nrow, 1) = sTensor33.block(0, 2, nrow, 1);
    S.block(0, 5, nrow, 1) = sTensor33.block(0, 1, nrow, 1);
}
