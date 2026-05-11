// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: MIT
// Copyright (C) 2019 - 2026 by the AxiSEM3D authors
//
// This file is part of the AxiSEM3D library. See the LICENSE file for details.
//
// -----------------------------------------------------------------------------

//  Clayton-Enquist ABC for fluid points in 3D

#ifndef ClaytonFluid3D_hpp
#define ClaytonFluid3D_hpp

#include "ClaytonFluid.hpp"
#include "eigen_point.hpp"

class ClaytonFluid3D : public ClaytonFluid {
  public:
  // constructor
  ClaytonFluid3D(const std::shared_ptr<FluidPoint>& fp,
      const axisem3d::eigen::DColX& rhoVp,
      const axisem3d::eigen::DColX& area) :
      ClaytonFluid(fp), mAreaOverRhoVp(area.cwiseQuotient(rhoVp).cast<numerical::Real>()) {
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
  // area / (rho * vp)
  const axisem3d::eigen::RColX mAreaOverRhoVp;

  ////////////////////////////////////////
  //////////////// static ////////////////
  ////////////////////////////////////////

  // workspace
  inline static axisem3d::eigen::RColX sVecR = axisem3d::eigen::RColX(0);
  inline static axisem3d::eigen::CColX sVecC = axisem3d::eigen::CColX(0);
};

#endif /* ClaytonFluid3D_hpp */
