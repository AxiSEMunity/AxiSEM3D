//
//  FluidPressure.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/5/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  pressure source on fluid element

#ifndef FluidPressure_hpp
#define FluidPressure_hpp

#include "FluidSource.hpp"
#include "eigen_element.hpp"

class FluidPressure: public FluidSource {
public:
    // constructor
    FluidPressure(std::unique_ptr<STF> &stf,
                  const std::shared_ptr<FluidElement> &element,
                  const eigen::CMatXN &pattern);
    
    // apply source at a time step
    void apply(double time) const;
    
private:
    // source pattern
    const eigen::CMatXN mPattern;
    
    
    ////////////////////////////////////////
    //////////////// static ////////////////
    ////////////////////////////////////////
    
    // workspace
    inline static eigen::CMatXN sPattern = eigen::CMatXN(0, spectral::nPEM);
};

#endif /* FluidPressure_hpp */
