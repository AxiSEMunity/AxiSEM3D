//
//  FluidSurfaceBoundary.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 7/30/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  stress-free boundary condition on fluid

#include "FluidSurfaceBoundary.hpp"

// point
#include "FluidPoint.hpp"

// domain
#include "Messaging.hpp"

// apply stress-free boundary condition on fluid
void FluidSurfaceBoundary::apply() const {
    // pressure ≡ 0 or accel ≡ 0
    // so, veloc = disp = everything ≡ 0
    for (const std::shared_ptr<FluidPoint> &fp: mFluidPoints) {
        fp->getFields().mStiff.setZero();
    }
}

// count info
int FluidSurfaceBoundary::
countInfo(const Messaging &msg) const {
    int count = 0;
    for (const auto &point: mFluidPoints) {
        if (!msg.pointInSmallerRank(point)) {
            count++;
        }
    }
    return count;
}
