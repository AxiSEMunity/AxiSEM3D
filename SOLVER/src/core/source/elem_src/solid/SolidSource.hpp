//
//  SolidSource.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/4/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  source on solid element

#ifndef SolidSource_hpp
#define SolidSource_hpp

#include "ElementSource.hpp"
class SolidElement;

class SolidSource: public ElementSource {
public:
    // constructor
    SolidSource(std::unique_ptr<STF> &stf,
                const std::shared_ptr<const SolidElement> &element):
    ElementSource(stf), mElement(element) {
        // nothing
    }
    
    // destructor
    virtual ~SolidSource() = default;
    
protected:
    // element pointer
    const std::shared_ptr<const SolidElement> mElement;
};

#endif /* SolidSource_hpp */
