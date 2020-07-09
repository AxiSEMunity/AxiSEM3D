//
//  SolidPoint.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 1/31/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  solid GLL point

#include "SolidPoint.hpp"
#include "point_time.hpp"

// constructor
SolidPoint::SolidPoint(int nr, const eigen::DRow2 &crds, int meshTag,
                       std::unique_ptr<const Mass> &mass,
                       const TimeScheme &timeScheme):
Point(nr, crds, meshTag, mass) {
    // fields
    point_time::createFields(*this, timeScheme);
    // mass
    mMass->checkCompatibility(mNr, true);
}


/////////////////////////// time loop ///////////////////////////
// stiff to accel
void SolidPoint::computeStiffToAccel() {
    // Nyquist
    if (mNr % 2 == 0) {
        mFields.mStiff.bottomRows(1).imag().setZero();
    }
    
    // stiff to accel in-place
    mMass->computeAccel(mFields.mStiff);
}
