//
//  SolidForce.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/5/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  force source on solid element

#include "SolidForce.hpp"
#include "SolidElement.hpp"

// constructor
SolidForce::SolidForce(std::unique_ptr<STF> &stf,
                       const std::shared_ptr<SolidElement> &element,
                       const eigen::CMatXN3 &pattern):
SolidSource(stf, element), mPattern(pattern) {
    // prepare
    element->prepareForceSource();
    
    // workspace
    if (sPattern.rows() < mPattern.rows()) {
        sPattern.resize(mPattern.rows(), spectral::nPEM * 3);
    }
}

// apply source at a time step
void SolidForce::apply(double time) const {
    int nu_1 = (int)mPattern.rows();
    sPattern.topRows(nu_1) = mPattern * mSTF->getValue(time);
    mElement->addForceSource(sPattern, nu_1);
}
