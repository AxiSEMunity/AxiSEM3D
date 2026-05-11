// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: MIT
// Copyright (C) 2019 - 2026 by the AxiSEM3D authors
//
// This file is part of the AxiSEM3D library. See the LICENSE file for details.
//
// -----------------------------------------------------------------------------

//  Clayton-Enquist ABC for solid points in 3D
//  f = v.n n * rpa + (v - v.n n) * rsa
//    = rsa v + v.k k, with k = sqrt(rpa - rsa) n
//  f -> stifness force
//  v -> velocity
//  n -> unit normal of surface
//  rsa -> rho * vs * area
//  rpa -> rho * vp * area

#ifndef ClaytonSolid3D_hpp
#define ClaytonSolid3D_hpp

#include "ClaytonSolid.hpp"
#include "eigen_point.hpp"

class ClaytonSolid3D : public ClaytonSolid {
  public:
  // constructor
  ClaytonSolid3D(const std::shared_ptr<SolidPoint>& sp,
      const axisem3d::eigen::DColX& rhoVp,
      const axisem3d::eigen::DColX& rhoVs,
      const axisem3d::eigen::DColX& area,
      const axisem3d::eigen::DMatX3& unitNormal) :
      ClaytonSolid(sp), mRSA(rhoVs.cwiseProduct(area).cast<numerical::Real>()),
      mK(((rhoVp - rhoVs).cwiseProduct(area).cwiseSqrt().asDiagonal() * unitNormal)
              .cast<numerical::Real>()) {
    // check compatibility
    checkCompatibility();
  }

  private:
  // check compatibility
  void
  checkCompatibility();

  public:
  // apply ABC
  void
  apply() const;

  private:
  // rsa = rho * vs * area
  const axisem3d::eigen::RColX mRSA;
  // k = sqrt(rpa - rsa) n
  const axisem3d::eigen::RMatX3 mK;

  ////////////////////////////////////////
  //////////////// static ////////////////
  ////////////////////////////////////////

  private:
  // workspace
  // V = FFT(velocity)
  inline static axisem3d::eigen::RMatX3 sVR = axisem3d::eigen::RMatX3(0, 3);
  // a = rsa V + V.k k
  inline static axisem3d::eigen::RMatX3 sAR = axisem3d::eigen::RMatX3(0, 3);
  inline static axisem3d::eigen::CMatX3 sAC = axisem3d::eigen::CMatX3(0, 3);
};

#endif /* ClaytonSolid3D_hpp */
