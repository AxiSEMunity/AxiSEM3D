//
//  FluidPressure.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/5/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  pressure source on fluid element

#include "FluidPressure.hpp"
#include "FluidElement.hpp"

// constructor
FluidPressure::FluidPressure(std::unique_ptr<STF> &stf,
                             const std::shared_ptr<FluidElement> &element,
                             const eigen::CMatXN &pattern):
FluidSource(stf, element), mPattern(pattern) {
    // prepare
    element->preparePressureSource();
    
    // workspace
    if (sPattern.rows() < mPattern.rows()) {
        sPattern.resize(mPattern.rows(), spectral::nPEM);
    }
}

// apply source at a time step
void FluidPressure::apply(double time) const {
    int nu_1 = (int)mPattern.rows();
    sPattern.topRows(nu_1) = mPattern * mSTF->getValue(time);
    mElement->addPressureSource(sPattern, nu_1);
}
