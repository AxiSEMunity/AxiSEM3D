//
//  ElementSource.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/4/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  source on element

#ifndef ElementSource_hpp
#define ElementSource_hpp

#include "STF.hpp"

class ElementSource {
public:
    // constructor
    ElementSource(std::unique_ptr<STF> &stf): mSTF(std::move(stf)) {
        // nothing
    }
    
    // destructor
    virtual ~ElementSource() = default;
    
    // apply source at a time step
    virtual void apply(double time) const = 0;
    
protected:
    // source-time function
    const std::unique_ptr<STF> mSTF;
};

#endif /* ElementSource_hpp */
