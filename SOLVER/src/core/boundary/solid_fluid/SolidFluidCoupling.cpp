//
//  SolidFluidCoupling.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 1/31/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  solid-fluid boundary condition

#include "SolidFluidCoupling.hpp"

// point
#include "SolidPoint.hpp"
#include "FluidPoint.hpp"

// constructor
SolidFluidCoupling::
SolidFluidCoupling(const std::shared_ptr<SolidPoint> &sp,
                   const std::shared_ptr<FluidPoint> &fp):
mSolidPoint(sp), mFluidPoint(fp) {
    if (mSolidPoint->getMeshTag() != mFluidPoint->getMeshTag()) {
        throw std::runtime_error("SolidFluidCoupling::SolidFluidCoupling || "
                                 "The coupled solid and fluid points have "
                                 "different mesh tags (positions).");
    }
}

// compute coupling
void SolidFluidCoupling::apply() const {
    // this order matters!
    coupleSolidToFluid(mSolidPoint->getFields().mDispl,
                       mFluidPoint->getFields().mStiff);
    coupleFluidToSolid(mFluidPoint->getFields().mStiff,
                       mSolidPoint->getFields().mStiff);
}


////////////////////////////// virtual //////////////////////////////
// check compatibility
void SolidFluidCoupling::checkCompatibility(int nr) const {
    if (mSolidPoint->getNr() != nr || mFluidPoint->getNr() != nr) {
        throw std::runtime_error("SolidFluidCoupling::checkCompatibility || "
                                 "Incompatible sizes.");
    }
}
