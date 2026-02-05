//
//  StructuredGridV3D.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 4/15/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  3D volumetric models based on structured grid

#ifndef StructuredGridV3D_hpp
#define StructuredGridV3D_hpp

#include "Volumetric3D.hpp"
#include "sg_tools.hpp"

class StructuredGridV3D: public Volumetric3D {
public:
    // constructor
    StructuredGridV3D(const std::string &modelName, const std::string &fname,
                      const std::array<std::string, 3> &crdVarNames,
                      const std::array<int, 3> &shuffleData,
                      bool sourceCentered, bool xy, bool ellipticity,
                      bool useDepth, bool depthSolid, bool undulated,
                      double lengthUnit, double angleUnit, bool center,
                      const std::vector<std::tuple<std::string, std::string,
                      double, ReferenceKind>> &propertyInfo, bool superOnly);
    
private:
    // using reference or undulated geometry
    bool usingUndulatedGeometry() const {
        return mUndulatedGeometry;
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
    
    // super-only: data stored only on super ranks
    bool isSuperOnly() const {
        return mSuperOnly;
    }
    
private:
    // file
    const std::string mFileName;
    const std::array<std::string, 3> mCrdVarNames;
    
    // horizontal options
    const bool mSourceCentered;
    const bool mXY;
    const bool mEllipticity;
    bool mLon360 = false;
    
    // vertical options
    const bool mUseDepth;
    const bool mDepthSolid;
    const bool mUndulatedGeometry;
    
    // use element center for scope check
    const bool mElementCenter;
    
    // data info
    // not using std::map because one property may be set twice in a model
    std::vector<std::string> mPropertyKeys;
    std::vector<std::string> mPropertyVarNames;
    std::vector<double> mPropertyFactors;
    std::vector<ReferenceKind> mPropertyReferenceKinds;
    std::unique_ptr<eigen::IMat66> mIndexCIJ;
    
    // grid
    std::unique_ptr<StructuredGrid<3, double>> mGrid = nullptr;
    
    // super only
    const bool mSuperOnly;
};

#endif /* StructuredGridV3D_hpp */
