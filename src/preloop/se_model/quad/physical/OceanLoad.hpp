// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: MIT
// Copyright (C) 2019 - 2026 by the AxiSEM3D authors
//
// This file is part of the AxiSEM3D library. See the LICENSE file for details.
//
// -----------------------------------------------------------------------------

//  ocean load on the suface

#ifndef OceanLoad_hpp
#define OceanLoad_hpp

#include "PhysicalProperty.hpp"

class OceanLoad {
  public:
  // add sum(rho * depth)
  void
  addSumRhoDepth(const axisem3d::eigen::arP_DColX& sumRD) {
    mSumRhoDepth.addGLL(sumRD);
  }

  // get pointwise
  axisem3d::eigen::arP_DColX
  getPointwise() const {
    return mSumRhoDepth.getPointwise();
  }

  // bool
  operator bool() const {
    return mSumRhoDepth;
  }

  private:
  // sum(rho * depth) over the water column
  PhysicalProperty<spectral::nPED> mSumRhoDepth;
};

#endif /* OceanLoad_hpp */
