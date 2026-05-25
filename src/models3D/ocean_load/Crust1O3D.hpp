//
//  Crust1O3D.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 17/10/2024.
//  Copyright © 2024 Kuangdai Leng. All rights reserved.
//

#ifndef Crust1O3D_hpp
#define Crust1O3D_hpp

#include "OceanLoad3D.hpp"

class Crust1O3D: public OceanLoad3D {
public:
    // constructor
    Crust1O3D(const std::string &modelName, double mWaterDensity,
              bool includeIceAsWater, bool ellipticity,
              int gaussianOrder, double gaussianDev);
    
private:
    // get sum(rho * depth)
    bool getSumRhoDepth(const eigen::DMatX3 &spz,
                        const eigen::DMat24 &nodalSZ,
                        eigen::DColX &sumRhoDepth) const;
    
    // verbose
    std::string verbose() const;
    
private:
    // model constants
    const int sNLayer = 9;
    const int sNLat = 180;
    const int sNLon = 360;
    
    // water density
    const double mWaterDensity;
    
    // treat ice as water load
    const bool mIncludeIceAsWater;
    
    // cosider ellipticity
    const bool mEllipticity;
    
    // smoothening
    const int mGaussianOrder;
    const double mGaussianDev;
    
    // depth at grid points
    eigen::DMatXX mDepth;
    eigen::DColX mGridLat, mGridLon;
};


#endif /* Crust1O3D_hpp */
