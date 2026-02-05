//
//  StructuredGridO3D.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 4/15/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  3D ocean-load models based on structured grid

#include "StructuredGridO3D.hpp"

// constructor
StructuredGridO3D::
StructuredGridO3D(const std::string &modelName, const std::string &fname,
                  const std::array<std::string, 2> &crdVarNames,
                  const std::array<int, 2> &shuffleData,
                  bool sourceCentered, bool xy, bool ellipticity,
                  double lengthUnit, double angleUnit,
                  const std::string &dataVarName, double factor,
                  bool superOnly):
OceanLoad3D(modelName), mFileName(fname), mCrdVarNames(crdVarNames),
mSourceCentered(sourceCentered), mXY(xy), mEllipticity(ellipticity),
mDataVarName(dataVarName), mFactor(factor), mSuperOnly(superOnly) {
    ////////////// init grid //////////////
    // info
    std::vector<std::pair<std::string, double>> dataInfo;
    dataInfo.push_back({mDataVarName, mFactor});
    
    // lambda
    auto initGrid = [this, &dataInfo, &shuffleData,
                     &lengthUnit, &angleUnit]() {
        // grid
        mGrid = std::make_unique<StructuredGrid<2, double>>
        (mFileName, mCrdVarNames, dataInfo, shuffleData);
        // coordinate units
        sg_tools::constructUnits(*mGrid, mSourceCentered, mXY, false,
                                 lengthUnit, angleUnit);
        // longitude range
        if (!mSourceCentered) {
            mLon360 = sg_tools::constructLon360(*mGrid, mModelName);
        }
    };
    
    // data
    if (mSuperOnly) {
        // constructor of StructuredGrid uses root + broadcast to read
        // right: mpi::super() after mpi::enterSuper()
        // wrong: mpi::root() after mpi::enterInfer()
        mpi::enterSuper();
        if (mpi::super()) {
            initGrid();
        }
        mpi::enterWorld();
    } else {
        initGrid();
    }
}

// get sum(rho * depth)
bool StructuredGridO3D::getSumRhoDepth(const eigen::DMatX3 &spz,
                                       const eigen::DMat24 &nodalSZ,
                                       eigen::DColX &sumRhoDepth) const {
    //////////////////////// coords ////////////////////////
    // check min/max
    const auto &gridCrds = mGrid->getGridCoords();
    if (!inplaneScope(nodalSZ, mSourceCentered && (!mXY),
                      gridCrds[0].front(), gridCrds[0].back(),
                      false, 0., 0., false, false)) {
        return false;
    }
    
    // compute grid coords
    const eigen::DMatX3 &crdGrid =
    coordsFromMeshToModel(spz, mSourceCentered, mXY, mEllipticity, mLon360,
                          false, false, mModelName);
    
    //////////////////////// values ////////////////////////
    // allocate and fill with zero
    int nCardinals = (int)spz.rows();
    sumRhoDepth = eigen::DColX::Zero(nCardinals);
    
    // point loop
    static const double err = std::numeric_limits<double>::lowest();
    bool oneInScope = false;
    for (int ipnt = 0; ipnt < nCardinals; ipnt++) {
        const eigen::DRow2 &horizontal = crdGrid.block(ipnt, 0, 1, 2);
        double val = mGrid->compute(horizontal, err);
        // check scope
        // ignore negative sum(rho * depth)
        if (val > 0.) {
            sumRhoDepth(ipnt) = val;
            oneInScope = true;
        }
    }
    return oneInScope;
}

// verbose
std::string StructuredGridO3D::verbose() const {
    if (!mpi::root()) {
        // grid uninitialized on infer ranks
        return "";
    }
    
    using namespace bstring;
    std::stringstream ss;
    // head
    ss << sg_tools::verboseHead(mModelName, "StructuredGridO3D", mFileName);
    
    // coords
    const auto &gcrds = mGrid->getGridCoords();
    ss << sg_tools::verboseCoords(mSourceCentered, mXY, false, false,
                                  {mCrdVarNames[0], mCrdVarNames[1]},
                                  {gcrds[0].front(), gcrds[1].front()},
                                  {gcrds[0].back(), gcrds[1].back()});
    
    // options
    if (!mSourceCentered) {
        ss << boxEquals(4, 22, "ellipticity correction", mEllipticity);
    }
    
    // data
    ss << boxSubTitle(2, "Data for sum(rho * depth)");
    ss << boxEquals(4, 19, "NetCDF variable", mDataVarName);
    const auto &minMax = mGrid->getDataRange();
    ss << boxEquals(4, 19, "data range", range(minMax(0, 0), minMax(0, 1)));
    ss << boxEquals(4, 19, "leader-only storage", mSuperOnly);
    return ss.str();
}
