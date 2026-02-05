//
//  Mass3D.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 1/25/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  3D mass

#include "Mass3D.hpp"
#include "fft.hpp"

// check compatibility
void Mass3D::checkCompatibility(int nr, bool solid) const {
    // check size
    if (mInvMass.rows() != nr) {
        throw std::runtime_error("Mass3D::checkCompatibility || "
                                 "Incompatible sizes.");
    }
    
    // workspace
    if (solid) {
        // expand workspace if needed
        if (sStiffR3.rows() < nr) {
            sStiffR3.resize(nr, 3);
        }
        // report request to FFT
        fft::gFFT_3.addNR(nr);
    } else {
        // expand workspace if needed
        if (sStiffR1.rows() < nr) {
            sStiffR1.resize(nr);
        }
        // report request to FFT
        fft::gFFT_1.addNR(nr);
    }
}

// compute accel in-place for fluid
void Mass3D::computeAccel(eigen::CColX &stiff1) const {
    // constants
    int nr = (int)mInvMass.rows();
    
    // FFT: Fourier => cardinal
    fft::gFFT_1.computeC2R(stiff1, sStiffR1, nr);
    
    // divide by mass in cardinal space
    sStiffR1.topRows(nr).array() *= mInvMass.array();
    
    // FFT: cardinal => Fourier
    fft::gFFT_1.computeR2C(sStiffR1, stiff1, nr);
}

// compute accel in-place for solid
void Mass3D::computeAccel(eigen::CMatX3 &stiff3) const {
    // constants
    int nr = (int)mInvMass.rows();
    
    // FFT: Fourier => cardinal
    fft::gFFT_3.computeC2R(stiff3, sStiffR3, nr);
    
    // divide by mass in cardinal space
    sStiffR3.topRows(nr).applyOnTheLeft(mInvMass.asDiagonal());
    
    // FFT: cardinal => Fourier
    fft::gFFT_3.computeR2C(sStiffR3, stiff3, nr);
}
