// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: MIT
// Copyright (C) 2019 - 2026 by the AxiSEM3D authors
//
// This file is part of the AxiSEM3D library. See the LICENSE file for details.
//
// -----------------------------------------------------------------------------

//  mass for both solid and fluid

#ifndef Mass_hpp
#define Mass_hpp

#include "eigen_point.hpp"

class Mass {
  public:
  // destructor
  virtual ~Mass() = default;

  // check compatibility
  virtual void
  checkCompatibility(int nr, bool solid) const {
    // nothing by default
  }

  // compute accel in-place for fluid
  virtual void
  computeAccel(axisem3d::eigen::CColX& stiff1) const = 0;

  // compute accel in-place for solid
  virtual void
  computeAccel(axisem3d::eigen::CMatX3& stiff3) const = 0;
};

#endif /* Mass_hpp */
