// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: MIT
// Copyright (C) 2019 - 2026 by the AxiSEM3D authors
//
// This file is part of the AxiSEM3D library. See the LICENSE file for details.
//
// -----------------------------------------------------------------------------

//  3D ocean-load models based on structured grid

#ifndef StructuredGridO3D_hpp
#define StructuredGridO3D_hpp

#include "OceanLoad3D.hpp"
#include "sg_tools.hpp"

class StructuredGridO3D : public OceanLoad3D {
  public:
  // constructor
  StructuredGridO3D(const std::string& modelName,
      const std::string& fname,
      const std::array<std::string, 2>& crdVarNames,
      const std::array<int, 2>& shuffleData,
      bool sourceCentered,
      bool xy,
      bool ellipticity,
      double lengthUnit,
      double angleUnit,
      const std::string& dataVarName,
      double factor,
      bool superOnly);

  private:
  // get sum(rho * depth)
  bool
  getSumRhoDepth(const axisem3d::eigen::DMatX3& spz,
      const axisem3d::eigen::DMat24& nodalSZ,
      axisem3d::eigen::DColX& sumRhoDepth) const;

  // verbose
  std::string
  verbose() const;

  // super-only: data stored only on super ranks
  bool
  isSuperOnly() const {
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
