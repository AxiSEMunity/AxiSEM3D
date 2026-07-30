// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: MIT
// Copyright (C) 2019 - 2026 by the AxiSEM3D authors
//
// This file is part of the AxiSEM3D library. See the LICENSE file for details.
//
// -----------------------------------------------------------------------------

//  3D geometric models based on structured grid

#ifndef StructuredGridG3D_hpp
#define StructuredGridG3D_hpp

#include "Geometric3D.hpp"
#include "sg_tools.hpp"

class StructuredGridG3D : public Geometric3D {
  public:
  // constructor
  StructuredGridG3D(const std::string& modelName,
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
      const ClampSmoothConfig& clampSmooth);

  private:
  // get undulation on an element
  bool
  getUndulation(
      const eigen::DMatX3& spz, const eigen::DMat24& nodalSZ, eigen::DColX& undulation) const;

  // get undulation on points
  bool
  getUndulation(const eigen::DMatX3& spz, eigen::DColX& undulation) const;

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

  // vertical options
  const bool mUseDepth;
  const bool mDepthSolid;

  // undulation range
  const double mInterface, mMin, mMax;

  // data
  const std::string mDataVarName;
  const double mFactor;

  // grid
  std::unique_ptr<StructuredGrid<2, double>> mGrid = nullptr;

  // super only
  const bool mSuperOnly;
};

#endif /* StructuredGridG3D_hpp */
