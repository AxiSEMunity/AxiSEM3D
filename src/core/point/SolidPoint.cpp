// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: MIT
// Copyright (C) 2019 - 2026 by the AxiSEM3D authors
//
// This file is part of the AxiSEM3D library. See the LICENSE file for details.
//
// -----------------------------------------------------------------------------

//  solid GLL point

#include "SolidPoint.hpp"
#include "point_time.hpp"

// constructor
SolidPoint::SolidPoint(int nr,
    const axisem3d::eigen::DRow2& crds,
    int meshTag,
    std::unique_ptr<const Mass>& mass,
    const TimeScheme& timeScheme) : Point(nr, crds, meshTag, mass) {
  // fields
  point_time::createFields(*this, timeScheme);
  // mass
  mMass->checkCompatibility(mNr, true);
}

/////////////////////////// time loop ///////////////////////////
// stiff to accel
void
SolidPoint::computeStiffToAccel() {
  // Nyquist
  if (mNr % 2 == 0) {
    mFields.mStiff.bottomRows(1).imag().setZero();
  }

  // stiff to accel in-place
  mMass->computeAccel(mFields.mStiff);
}
