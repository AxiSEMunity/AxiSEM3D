//
//  StructuredGridG3D.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 4/15/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  3D geometric models based on structured grid

#include "StructuredGridG3D.hpp"

// constructor
StructuredGridG3D::
StructuredGridG3D(const std::string &modelName, const std::string &fname,
                  const std::array<std::string, 2> &crdVarNames,
                  const std::array<int, 2> &shuffleData,
                  bool sourceCentered, bool xy, bool ellipticity,
                  bool useDepth, bool depthSolid,
                  double interface, double min, double max,
                  double lengthUnit, double angleUnit,
                  const std::string &dataVarName, double factor,
                  bool superOnly):
Geometric3D(modelName), mFileName(fname), mCrdVarNames(crdVarNames),
mSourceCentered(sourceCentered), mXY(xy), mEllipticity(ellipticity),
mUseDepth(useDepth), mDepthSolid(depthSolid),
mInterface(interface * lengthUnit),
mMin(min * lengthUnit), mMax(max * lengthUnit),
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

// get undulation on an element
bool StructuredGridG3D::getUndulation(const eigen::DMatX3 &spz,
                                      const eigen::DMat24 &nodalSZ,
                                      eigen::DColX &undulation) const {
    // check inplane scope
    const auto &gridCrds = mGrid->getGridCoords();
    if (!inplaneScope(nodalSZ, mSourceCentered && (!mXY),
                      gridCrds[0].front(), gridCrds[0].back(),
                      true, mMin, mMax, mUseDepth, mDepthSolid)) {
        return false;
    }
    return getUndulation(spz, undulation);
}

// get undulation on points
bool StructuredGridG3D::getUndulation(const eigen::DMatX3 &spz,
                                      eigen::DColX &undulation) const {
    // cannot use undulated geometry to locate source and receivers
    // for super-only storage
    if (mGrid == nullptr) {
        throw std::runtime_error
        ("StructuredGridG3D::getUndulation || "
         "Option store_grid_only_on_leaders for StructuredGridG3D || "
         "is incompatible with option undulated_geometry for || "
         "the vertical locations of sources and receivers.");
    }
    
    //////////////////////// coords ////////////////////////
    // compute grid coords
    const eigen::DMatX3 &crdGrid =
    coordsFromMeshToModel(spz, mSourceCentered, mXY, mEllipticity, mLon360,
                          mUseDepth, mDepthSolid, mModelName);
    
    //////////////////////// values ////////////////////////
    // allocate and fill with zero
    int nCardinals = (int)spz.rows();
    undulation = eigen::DColX::Zero(nCardinals);
    
    // point loop
    static const double err = std::numeric_limits<double>::lowest();
    bool oneInScope = false;
    for (int ipnt = 0; ipnt < nCardinals; ipnt++) {
        // check vertical scope
        double rOrD = crdGrid(ipnt, 2);
        if (rOrD < mMin || rOrD > mMax) {
            continue;
        }
        // horizontal interpolation
        const eigen::DRow2 &horizontal = crdGrid.block(ipnt, 0, 1, 2);
        double val = mGrid->compute(horizontal, err);
        // check horizontal scope
        if (val > err * .9) {
            // interpolation along vertical
            if (rOrD < mInterface) {
                undulation(ipnt) = val / (mInterface - mMin) * (rOrD - mMin);
            } else {
                undulation(ipnt) = val / (mMax - mInterface) * (mMax - rOrD);
            }
            oneInScope = true;
        }
    }
    return oneInScope;
}

// verbose
std::string StructuredGridG3D::verbose() const {
    if (!mpi::root()) {
        // grid uninitialized on infer ranks
        return "";
    }
    
    using namespace bstring;
    std::stringstream ss;
    // head
    ss << sg_tools::verboseHead(mModelName, "StructuredGridG3D", mFileName);
    
    // coords
    const auto &gcrds = mGrid->getGridCoords();
    ss << sg_tools::verboseCoords(mSourceCentered, mXY, true, mUseDepth,
                                  {mCrdVarNames[0], mCrdVarNames[1], "N/A"},
                                  {gcrds[0].front(), gcrds[1].front(), mMin},
                                  {gcrds[0].back(), gcrds[1].back(), mMax});
    // width
    int width = 19;
    if (!mSourceCentered) {
        width = 22;
    }
    if (mUseDepth) {
        width = 25;
    }
    // interface
    ss << boxEquals(4, width, mUseDepth ? "interface depth" :
                    "interface radius", mInterface);
    // options
    if (!mSourceCentered) {
        ss << boxEquals(4, width, "ellipticity correction", mEllipticity);
    }
    if (mUseDepth) {
        ss << boxEquals(4, width, "depth below solid surface", mDepthSolid);
    }
    
    // undulation
    ss << boxSubTitle(2, "Undulation data");
    ss << boxEquals(4, 19, "NetCDF variable", mDataVarName);
    const auto &minMax = mGrid->getDataRange();
    ss << boxEquals(4, 19, "data range", range(minMax(0, 0), minMax(0, 1)));
    ss << boxEquals(4, 19, "leader-only storage", mSuperOnly);
    return ss.str();
}
