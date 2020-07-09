//
//  SolidMoment.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/5/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  moment source on solid element

#ifndef SolidMoment_hpp
#define SolidMoment_hpp

#include "SolidSource.hpp"
#include "eigen_element.hpp"

class SolidMoment: public SolidSource {
public:
    // constructor
    SolidMoment(std::unique_ptr<STF> &stf,
                const std::shared_ptr<SolidElement> &element,
                const eigen::CMatXN6 &pattern);
    
    // apply source at a time step
    void apply(double time) const;
    
private:
    // source pattern
    const eigen::CMatXN6 mPattern;
    
    
    ////////////////////////////////////////
    //////////////// static ////////////////
    ////////////////////////////////////////
    
    // workspace
    inline static eigen::CMatXN6 sPattern =
    eigen::CMatXN6(0, spectral::nPEM * 6);
};

#endif /* SolidMoment_hpp */
