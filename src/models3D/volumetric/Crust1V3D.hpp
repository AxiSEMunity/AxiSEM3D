//
//  Crust1V3D.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 10/14/24.
//  Copyright © 2024 Kuangdai Leng. All rights reserved.
//

//  3D volumetric model based on Crust 1.0

#ifndef Crust1V3D_hpp
#define Crust1V3D_hpp

#include "Volumetric3D.hpp"

class Crust1V3D: public Volumetric3D {
public:
    // constructor
    Crust1V3D(const std::string &modelName, double rSurf, double rMoho,
              bool includeSediment, bool includeIce, bool ellipticity,
              const ExodusMesh &exodusMesh);
    
private:
    // using reference or undulated geometry
    bool usingUndulatedGeometry() const {
        return false;
    }
    
    // super-only: data stored only on super ranks
    bool isSuperOnly() const {
        return false;
    }
    
    // get property info
    void getPropertyInfo(std::vector<std::string> &propKeys,
                         std::vector<ReferenceKind> &refKinds) const;
    
    // get properties
    bool getProperties(const eigen::DMatX3 &spz,
                       const eigen::DMat24 &nodalSZ,
                       eigen::IMatXX &inScopes,
                       eigen::DMatXX &propValues) const;
    
    // verbose
    std::string verbose() const;
    
    // get properties at a point
    std::vector<double>
    getPropertiesPoint(double r, double lat, double lon) const;
    
private:
    // model constants
    const int sNLayer = 9;
    const int sNLat = 180;
    const int sNLon = 360;
    
    // radii of reference sphere
    const double mRSurf;
    const double mRMoho;
    
    // options
    // include sediment or not
    bool mIncludeSediment;
    // include ice or not
    const bool mIncludeIce;
    // cosider ellipticity
    const bool mEllipticity;
    
    // element boundaries in mesh
    std::vector<double> mElementBoundaries;
    
    // original data mapped onto reference sphere
    eigen::DMatXX mRl;
    eigen::DMatXX mVp;
    eigen::DMatXX mVs;
    eigen::DMatXX mRh;
    
    // thickness weighted data mapped onto reference sphere
    eigen::DColX mRlGLL;
    eigen::DMatXX mVpGLL;
    eigen::DMatXX mVsGLL;
    eigen::DMatXX mRhGLL;
    
    // lat and lon grid
    eigen::DColX mGridLat, mGridLon;
};

#endif /* Crust1V3D_hpp */
