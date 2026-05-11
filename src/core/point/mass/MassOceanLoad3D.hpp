// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: MIT
// Copyright (C) 2019 - 2026 by the AxiSEM3D authors
//
// This file is part of the AxiSEM3D library. See the LICENSE file for details.
//
// -----------------------------------------------------------------------------

//  3D mass with ocean load
//  used on solid-fluid boundary when the fluid is modelled as load
//  a = F.n n / (m + m0) + (F - F.n n) / m
//    = [m F.n n + (m + m0) (F - F.n n)] / [m (m + m0)]
//    = [(m + m0) F - m0 F.n n] / [m (m + m0)]
//    = F / m - m0 / [m (m + m0)] F.n n
//    = im F - F.k k, with k = sqrt(m0 / [m (m + m0)]) n
//  a -> acceleration
//  F -> force
//  n -> unit normal of surface
//  m -> mass of solid
//  m0 -> mass of water column above

#ifndef MassOceanLoad3D_hpp
#define MassOceanLoad3D_hpp

#include "Mass.hpp"

class MassOceanLoad3D : public Mass {
  public:
  // constructor
  MassOceanLoad3D(const axisem3d::eigen::DColX& mass,
      const axisem3d::eigen::DColX& massOcean,
      const axisem3d::eigen::DMatX3& unitNormal);

  // check compatibility
  void
  checkCompatibility(int nr, bool solid) const;

  // compute accel in-place for fluid
  void
  computeAccel(axisem3d::eigen::CColX& stiff1) const {
    throw std::runtime_error("MassOceanLoad3D::computeAccel || "
                             "Incompatible types: "
                             "ocean load on fluid point.");
  }

  // compute accel in-place for solid
  void
  computeAccel(axisem3d::eigen::CMatX3& stiff3) const;

  private:
  // im = 1 / m
  const axisem3d::eigen::RColX mIM;
  // k = sqrt(m0 / [m (m + m0)]) n
  const axisem3d::eigen::RMatX3 mK;

  ////////////////////////////////////////
  //////////////// static ////////////////
  ////////////////////////////////////////

  // workspace
  // F = FFT(stiff)
  inline static axisem3d::eigen::RMatX3 sF = axisem3d::eigen::RMatX3(0, 3);
  // a = im F - F.k k
  inline static axisem3d::eigen::RMatX3 sA = axisem3d::eigen::RMatX3(0, 3);
};

#endif /* MassOceanLoad3D_hpp */
