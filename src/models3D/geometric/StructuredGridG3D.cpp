// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: MIT
// Copyright (C) 2019 - 2026 by the AxiSEM3D authors
//
// This file is part of the AxiSEM3D library. See the LICENSE file for details.
//
// -----------------------------------------------------------------------------

//  3D geometric models based on structured grid

#include "StructuredGridG3D.hpp"
#include <cmath>

namespace {
  double
  uniformSpacing(
      const std::vector<double>& coords, const std::string& axis, const std::string& modelName) {
    const double spacing = (coords.back() - coords.front()) / (coords.size() - 1);
    const double tolerance = std::max(1., std::abs(spacing)) * 1e-8;
    for (int index = 1; index < coords.size(); index++) {
      if (std::abs(coords[index] - coords[index - 1] - spacing) > tolerance) {
        throw std::runtime_error("StructuredGridG3D::StructuredGridG3D || "
                                 "Clamp and Gaussian smoothing requires a uniformly spaced "
                                 "horizontal grid. || Axis: " +
            axis + " || Model name: " + modelName);
      }
    }
    return spacing;
  }

} // namespace

// constructor
StructuredGridG3D::StructuredGridG3D(const std::string& modelName,
    const std::string& fname,
    const std::array<std::string, 2>& crdVarNames,
    const std::array<int, 2>& shuffleData,
    bool sourceCentered,
    bool xy,
    bool ellipticity,
    bool useDepth,
    bool depthSolid,
    double interface,
    double min,
    double max,
    double lengthUnit,
    double angleUnit,
    const std::string& dataVarName,
    double factor,
    bool superOnly,
    const ClampSmoothConfig& clampSmooth) :
    Geometric3D(modelName, clampSmooth), mFileName(fname), mCrdVarNames(crdVarNames),
    mSourceCentered(sourceCentered), mXY(xy), mEllipticity(ellipticity), mUseDepth(useDepth),
    mDepthSolid(depthSolid), mInterface(interface * lengthUnit), mMin(min * lengthUnit),
    mMax(max * lengthUnit), mDataVarName(dataVarName), mFactor(factor), mSuperOnly(superOnly) {
  ////////////// init grid //////////////
  // info
  std::vector<std::pair<std::string, double>> dataInfo;
  // With clamp enabled, materialise the factor before the nonlinear operation.
  dataInfo.push_back({mDataVarName, mClampSmooth.enabled ? 1. : mFactor});

  // lambda
  auto initGrid = [this, &dataInfo, &shuffleData, &lengthUnit, &angleUnit]() {
    // grid
    mGrid =
        std::make_unique<StructuredGrid<2, double>>(mFileName, mCrdVarNames, dataInfo, shuffleData);
    // coordinate units
    sg_tools::constructUnits(*mGrid, mSourceCentered, mXY, false, lengthUnit, angleUnit);
    // longitude range
    if (!mSourceCentered) {
      mLon360 = sg_tools::constructLon360(*mGrid, mModelName);
    }

    // clamp and smooth each unique data field on the native horizontal grid
    if (mClampSmooth.enabled) {
      const auto& coords = mGrid->getGridCoords();
      const std::array<double, 2> spacing = {uniformSpacing(coords[0], mCrdVarNames[0], mModelName),
          uniformSpacing(coords[1], mCrdVarNames[1], mModelName)};

      SmoothGridUnit gridUnit;
      if (geodesy::isCartesian()) {
        if (!mSourceCentered || !mXY) {
          throw std::runtime_error("StructuredGridG3D::StructuredGridG3D || "
                                   "Clamp and Gaussian smoothing in Cartesian geometry requires "
                                   "an XY grid in metres. || Model name: " +
              mModelName);
        }
        gridUnit = SmoothGridUnit::METER;
      } else if (!mSourceCentered) {
        // latitude and longitude are stored internally in degrees
        gridUnit = SmoothGridUnit::DEGREE;
      } else {
        // spherical distance and azimuth are stored internally in radians
        gridUnit = SmoothGridUnit::RADIAN;
      }

      auto& gridData = mGrid->getGridData();
      const int n0 = (int)gridData.dimension(1);
      const int n1 = (int)gridData.dimension(2);
      for (int idata = 0; idata < mGrid->numUniqueData(); idata++) {
        eigen::DMatXX physical(n0, n1);
        for (int i0 = 0; i0 < n0; i0++) {
          for (int i1 = 0; i1 < n1; i1++) {
            physical(i0, i1) = gridData(idata, i0, i1) * mFactor;
          }
        }
        applyClampSmooth(physical, spacing, gridUnit);
        for (int i0 = 0; i0 < n0; i0++) {
          for (int i1 = 0; i1 < n1; i1++) {
            gridData(idata, i0, i1) = physical(i0, i1);
          }
        }
      }
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
bool
StructuredGridG3D::getUndulation(
    const eigen::DMatX3& spz, const eigen::DMat24& nodalSZ, eigen::DColX& undulation) const {
  // check inplane scope
  const auto& gridCrds = mGrid->getGridCoords();
  if (!inplaneScope(nodalSZ,
          mSourceCentered && (!mXY),
          gridCrds[0].front(),
          gridCrds[0].back(),
          true,
          mMin,
          mMax,
          mUseDepth,
          mDepthSolid)) {
    return false;
  }
  return getUndulation(spz, undulation);
}

// get undulation on points
bool
StructuredGridG3D::getUndulation(const eigen::DMatX3& spz, eigen::DColX& undulation) const {
  // cannot use undulated geometry to locate source and receivers
  // for super-only storage
  if (mGrid == nullptr) {
    throw std::runtime_error("StructuredGridG3D::getUndulation || "
                             "Option store_grid_only_on_leaders for StructuredGridG3D || "
                             "is incompatible with option undulated_geometry for || "
                             "the vertical locations of sources and receivers.");
  }

  //////////////////////// coords ////////////////////////
  // compute grid coords
  const eigen::DMatX3& crdGrid = coordsFromMeshToModel(
      spz, mSourceCentered, mXY, mEllipticity, mLon360, mUseDepth, mDepthSolid, mModelName);

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
    const eigen::DRow2& horizontal = crdGrid.block(ipnt, 0, 1, 2);
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
std::string
StructuredGridG3D::verbose() const {
  if (!mpi::root()) {
    // grid uninitialized on infer ranks
    return "";
  }

  using namespace bstring;
  std::stringstream ss;
  // head
  ss << sg_tools::verboseHead(mModelName, "StructuredGridG3D", mFileName);

  // coords
  const auto& gcrds = mGrid->getGridCoords();
  ss << sg_tools::verboseCoords(mSourceCentered,
      mXY,
      true,
      mUseDepth,
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
  ss << boxEquals(4, width, mUseDepth ? "interface depth" : "interface radius", mInterface);
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
  const auto& minMax = mGrid->getDataRange();
  ss << boxEquals(4, 19, "data range", range(minMax(0, 0), minMax(0, 1)));
  ss << boxEquals(4, 19, "leader-only storage", mSuperOnly);
  ss << verboseClampSmooth(2, 19);
  return ss.str();
}
