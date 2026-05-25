//
//  Crust1G3D.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 10/14/24.
//  Copyright © 2024 Kuangdai Leng. All rights reserved.
//

//  3D geometric model based on Crust 1.0

#ifndef Crust1G3D_hpp
#define Crust1G3D_hpp

#include "Geometric3D.hpp"

class Crust1G3D: public Geometric3D {
public:
    // constructor
    Crust1G3D(const std::string &modelName,
              double rSurf, double rMoho, double rBase,
              bool includeSediment, bool includeIce, bool ellipticity,
              double surfaceFactor, double mohoFactor,
              int gaussianOrder, double gaussianDev);
    
private:
    // get undulation on an element
    bool getUndulation(const eigen::DMatX3 &spz,
                       const eigen::DMat24 &nodalSZ,
                       eigen::DColX &undulation) const;
    
    // get undulation on points
    bool getUndulation(const eigen::DMatX3 &spz,
                       eigen::DColX &undulation) const;
    
    // verbose
    std::string verbose() const;
    
    // get undulation on point
    double getUndulationPoint(double r, double lat, double lon,
                              double rElemCenter) const;
    
private:
    // model constants
    const int sNLayer = 9;
    const int sNLat = 180;
    const int sNLon = 360;
    
    // radii of reference sphere
    const double mRSurf;
    const double mRMoho;
    const double mRBase;
    
    // options
    // include sediment or not
    bool mIncludeSediment;
    // include ice or not
    const bool mIncludeIce;
    // cosider ellipticity
    const bool mEllipticity;
    
    // strengthening factor
    const double mSurfFactor;
    const double mMohoFactor;
    
    // smoothening
    const int mGaussianOrder;
    const double mGaussianDev;
    
    // deltaR at surface and moho
    eigen::DMatXX mDeltaRSurf;
    eigen::DMatXX mDeltaRMoho;
    eigen::DColX mGridLat, mGridLon;
};

#endif /* Crust1G3D_hpp */
