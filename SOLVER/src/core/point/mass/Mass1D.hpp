//
//  Mass1D.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 1/25/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  1D mass

#ifndef Mass1D_hpp
#define Mass1D_hpp

#include "Mass.hpp"

class Mass1D: public Mass {
public:
    // constructor
    Mass1D(double mass):
    mInvMass((numerical::Real)(1. / mass)) {
        // nothing
    }
    
    // compute accel in-place for fluid
    void computeAccel(eigen::CColX &stiff1) const {
        stiff1 *= mInvMass;
    }
    
    // compute accel in-place for solid
    void computeAccel(eigen::CMatX3 &stiff3) const {
        stiff3 *= mInvMass;
    }
    
private:
    // inverse of mass
    const numerical::Real mInvMass;
};

#endif /* Mass1D_hpp */
