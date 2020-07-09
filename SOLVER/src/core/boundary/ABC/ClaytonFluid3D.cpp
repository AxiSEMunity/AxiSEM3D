//
//  ClaytonFluid3D.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 7/30/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  Clayton-Enquist ABC for fluid points in 3D

#include "ClaytonFluid3D.hpp"
#include "FluidPoint.hpp"
#include "fft.hpp"

// check compatibility
void ClaytonFluid3D::checkCompatibility() {
    // check size
    int nr = mFluidPoint->getNr();
    if (nr != mAreaOverRhoVp.rows()) {
        throw std::runtime_error("ClaytonFluid3D::checkCompatibility ||"
                                 "Incompatible sizes.");
    }
    
    // workspace
    if (sVecR.rows() < nr) {
        sVecR.resize(nr);
        sVecC.resize(nr / 2 + 1);
    }
    
    // report request to FFT
    fft::gFFT_1.addNR(nr);
}

// apply ABC
void ClaytonFluid3D::apply() const {
    // get fields
    const eigen::CColX &veloc = mFluidPoint->getFields().mVeloc;
    eigen::CColX &stiff = mFluidPoint->getFields().mStiff;
    
    // constants
    int nr = mFluidPoint->getNr();
    int nu_1 = nr / 2 + 1;
    
    // FFT: Fourier => cardinal
    fft::gFFT_1.computeC2R(veloc, sVecR, nr);
    
    // multiply by area / (rho * vp) in cardinal space
    sVecR.topRows(nr).array() *= mAreaOverRhoVp.array();
    
    // FFT: cardinal => Fourier
    fft::gFFT_1.computeR2C(sVecR, sVecC, nr);
    
    // subtract
    stiff -= sVecC.topRows(nu_1);
}
