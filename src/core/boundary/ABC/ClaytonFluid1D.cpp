// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: MIT
// Copyright (C) 2019 - 2026 by the AxiSEM3D authors
//
// This file is part of the AxiSEM3D library. See the LICENSE file for details.
//
// -----------------------------------------------------------------------------

//  Clayton-Enquist ABC for fluid points in 1D

#include "ClaytonFluid1D.hpp"
#include "FluidPoint.hpp"

// apply ABC
void
ClaytonFluid1D::apply() const {
  // get fields
  const axisem3d::eigen::CColX& veloc = mFluidPoint->getFields().mVeloc;
  axisem3d::eigen::CColX& stiff = mFluidPoint->getFields().mStiff;

  // apply
  stiff -= veloc * mAreaOverRhoVp;
}
