//
//  Sponge.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 5/29/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  sponge ABC

#ifndef Sponge_hpp
#define Sponge_hpp

#include "numerical.hpp"

template <class PointSF>
class Sponge {
public:
    // constructor
    Sponge(const std::shared_ptr<PointSF> &point, double gamma):
    mPoint(point), mGamma(gamma) {
        // nothing
    }
    
    // get point
    const std::shared_ptr<PointSF> &getPoint() const {
        return mPoint;
    }
    
    // apply ABC
    // must be called after point->computeStiffToAccel(),
    // so here "stiff" has been converted to acceleration
    void apply() const {
        auto &stiff = mPoint->getFields().mStiff;
        const auto &veloc = mPoint->getFields().mVeloc;
        const auto &displ = mPoint->getFields().mDispl;
        static const numerical::Real two = 2.;
        // update acceleration
        stiff -= (two * mGamma) * veloc + (mGamma * mGamma) * displ;
    }
    
private:
    // point
    const std::shared_ptr<PointSF> mPoint;
    
    // gamma
    const numerical::Real mGamma;
};

#endif /* Sponge_hpp */
