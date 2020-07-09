//
//  Ellipticity.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 5/29/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  ellipticity

#ifndef Ellipticity_hpp
#define Ellipticity_hpp

#include "Geometric3D.hpp"

class Ellipticity: public Geometric3D {
public:
    // constructor
    Ellipticity(const std::string &modelName): Geometric3D(modelName) {
        // nothing
    }
    
private:
    // get undulation on points
    bool getUndulation(const eigen::DMatX3 &spz,
                       eigen::DColX &undulation) const;
    
    // verbose
    std::string verbose() const;
};

#endif /* Ellipticity_hpp */
