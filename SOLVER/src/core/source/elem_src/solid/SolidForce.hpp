//
//  SolidForce.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/5/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  force source on solid element

#ifndef SolidForce_hpp
#define SolidForce_hpp

#include "SolidSource.hpp"
#include "eigen_element.hpp"

class SolidForce: public SolidSource {
public:
    // constructor
    SolidForce(std::unique_ptr<STF> &stf,
               const std::shared_ptr<SolidElement> &element,
               const eigen::CMatXN3 &pattern);
    
    // apply source at a time step
    void apply(double time) const;
    
private:
    // source pattern
    const eigen::CMatXN3 mPattern;
    
    
    ////////////////////////////////////////
    //////////////// static ////////////////
    ////////////////////////////////////////
    
    // workspace
    inline static eigen::CMatXN3 sPattern =
    eigen::CMatXN3(0, spectral::nPEM * 3);
};

#endif /* SolidForce_hpp */
