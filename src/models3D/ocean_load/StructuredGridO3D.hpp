//
//  StructuredGridO3D.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 4/15/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  3D ocean-load models based on structured grid

#ifndef StructuredGridO3D_hpp
#define StructuredGridO3D_hpp

#include "OceanLoad3D.hpp"
#include "sg_tools.hpp"

class StructuredGridO3D: public OceanLoad3D {
public:
    // constructor
    StructuredGridO3D(const std::string &modelName, const std::string &fname,
                      const std::array<std::string, 2> &crdVarNames,
                      const std::array<int, 2> &shuffleData,
                      bool sourceCentered, bool xy, bool ellipticity,
                      double lengthUnit, double angleUnit,
                      const std::string &dataVarName, double factor,
                      bool superOnly);
    
private:
    // get sum(rho * depth)
    bool getSumRhoDepth(const eigen::DMatX3 &spz,
                        const eigen::DMat24 &nodalSZ,
                        eigen::DColX &sumRhoDepth) const;
    
    // verbose
    std::string verbose() const;
    
    // super-only: data stored only on super ranks
    bool isSuperOnly() const {
        return mSuperOnly;
    }
    
private:
    // file
    const std::string mFileName;
    const std::array<std::string, 2> mCrdVarNames;
    
    // horizontal options
    const bool mSourceCentered;
    const bool mXY;
    const bool mEllipticity;
    bool mLon360 = false;
    
    // data
    const std::string mDataVarName;
    const double mFactor;
    
    // grid
    std::unique_ptr<StructuredGrid<2, double>> mGrid = nullptr;
    
    // super only
    const bool mSuperOnly;
};

#endif /* StructuredGridO3D_hpp */
