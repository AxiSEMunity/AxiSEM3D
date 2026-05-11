// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: MIT
// Copyright (C) 2019 - 2026 by the AxiSEM3D authors
//
// This file is part of the AxiSEM3D library. See the LICENSE file for details.
//
// -----------------------------------------------------------------------------

//  3D mass

#ifndef Mass3D_hpp
#define Mass3D_hpp

#include "Mass.hpp"

class Mass3D : public Mass {
  public:
  // constructor
  Mass3D(const axisem3d::eigen::DColX& mass) :
      mInvMass(mass.cwiseInverse().cast<numerical::Real>()) {
    // nothing
  }

  // check compatibility
  void
  checkCompatibility(int nr, bool solid) const;

  // compute accel in-place for fluid
  void
  computeAccel(axisem3d::eigen::CColX& stiff1) const;

  // compute accel in-place for solid
  void
  computeAccel(axisem3d::eigen::CMatX3& stiff3) const;

  private:
  // inverse of mass
  const axisem3d::eigen::RColX mInvMass;

  ////////////////////////////////////////
  //////////////// static ////////////////
  ////////////////////////////////////////

  // workspace
  inline static axisem3d::eigen::RColX sStiffR1 = axisem3d::eigen::RColX(0);
  inline static axisem3d::eigen::RMatX3 sStiffR3 = axisem3d::eigen::RMatX3(0, 3);
};

#endif /* Mass3D_hpp */
