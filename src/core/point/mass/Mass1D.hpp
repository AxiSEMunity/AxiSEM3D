// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: MIT
// Copyright (C) 2019 - 2026 by the AxiSEM3D authors
//
// This file is part of the AxiSEM3D library. See the LICENSE file for details.
//
// -----------------------------------------------------------------------------

//  1D mass

#ifndef Mass1D_hpp
#define Mass1D_hpp

#include "Mass.hpp"

class Mass1D : public Mass {
  public:
  // constructor
  Mass1D(double mass) : mInvMass((numerical::Real)(1. / mass)) {
    // nothing
  }

  // compute accel in-place for fluid
  void
  computeAccel(axisem3d::eigen::CColX& stiff1) const {
    stiff1 *= mInvMass;
  }

  // compute accel in-place for solid
  void
  computeAccel(axisem3d::eigen::CMatX3& stiff3) const {
    stiff3 *= mInvMass;
  }

  private:
  // inverse of mass
  const numerical::Real mInvMass;
};

#endif /* Mass1D_hpp */
