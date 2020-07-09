//
//  FluidSource.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/4/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  source on fluid element

#ifndef FluidSource_hpp
#define FluidSource_hpp

#include "ElementSource.hpp"
class FluidElement;

class FluidSource: public ElementSource {
public:
    // constructor
    FluidSource(std::unique_ptr<STF> &stf,
                const std::shared_ptr<const FluidElement> &element):
    ElementSource(stf), mElement(element) {
        // nothing
    }
    
    // destructor
    virtual ~FluidSource() = default;
    
protected:
    // element pointer
    const std::shared_ptr<const FluidElement> mElement;
};

#endif /* FluidSource_hpp */
