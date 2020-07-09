//
//  Mass.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 1/24/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  mass for both solid and fluid

#ifndef Mass_hpp
#define Mass_hpp

#include "eigen_point.hpp"

class Mass {
public:
    // destructor
    virtual ~Mass() = default;
    
    // check compatibility
    virtual void checkCompatibility(int nr, bool solid) const {
        // nothing by default
    }
    
    // compute accel in-place for fluid
    virtual void computeAccel(eigen::CColX &stiff1) const = 0;
    
    // compute accel in-place for solid
    virtual void computeAccel(eigen::CMatX3 &stiff3) const = 0;
};

#endif /* Mass_hpp */
