//
//  ClaytonFluid1D.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 7/29/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  Clayton-Enquist ABC for fluid points in 1D

#ifndef ClaytonFluid1D_hpp
#define ClaytonFluid1D_hpp

#include "ClaytonFluid.hpp"
#include "numerical.hpp"

class ClaytonFluid1D: public ClaytonFluid {
public:
    // constructor
    ClaytonFluid1D(const std::shared_ptr<FluidPoint> &fp,
                   double rhoVp, double area):
    ClaytonFluid(fp), mAreaOverRhoVp(area / rhoVp) {
        // nothing
    }
    
    // apply ABC
    void apply() const;
    
private:
    // area / (rho * vp)
    const numerical::Real mAreaOverRhoVp;
};

#endif /* ClaytonFluid1D_hpp */
