//
//  ClaytonSolid.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 7/29/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  Clayton-Enquist ABC for solid points

#ifndef ClaytonSolid_hpp
#define ClaytonSolid_hpp

// point
#include <memory>
class SolidPoint;

class ClaytonSolid {
public:
    // constructor
    ClaytonSolid(const std::shared_ptr<SolidPoint> &sp):
    mSolidPoint(sp) {
        // nothing
    }
    
    // get point
    const std::shared_ptr<SolidPoint> &getPoint() const {
        return mSolidPoint;
    }
    
    // destructor
    virtual ~ClaytonSolid() = default;
    
    // apply ABC
    virtual void apply() const = 0;
    
protected:
    // point
    const std::shared_ptr<SolidPoint> mSolidPoint;
};

#endif /* ClaytonSolid_hpp */
