//
//  OceanLoad.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/30/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  ocean load on the suface

#ifndef OceanLoad_hpp
#define OceanLoad_hpp

#include "PhysicalProperty.hpp"

class OceanLoad {
public:
    // add sum(rho * depth)
    void addSumRhoDepth(const eigen::arP_DColX &sumRD) {
        mSumRhoDepth.addGLL(sumRD);
    }
    
    // get pointwise
    eigen::arP_DColX getPointwise() const {
        return mSumRhoDepth.getPointwise();
    }
    
    // bool
    operator bool() const {
        return mSumRhoDepth;
    }
    
private:
    // sum(rho * depth) over the water column
    PhysicalProperty<spectral::nPED> mSumRhoDepth;
};

#endif /* OceanLoad_hpp */
