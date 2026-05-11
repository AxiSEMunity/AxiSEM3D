// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: MIT
// Copyright (C) 2019 - 2026 by the AxiSEM3D authors
//
// This file is part of the AxiSEM3D library. See the LICENSE file for details.
//
// -----------------------------------------------------------------------------

//  fluid GLL point

#include "FluidPoint.hpp"
#include "point_time.hpp"

// constructor
FluidPoint::FluidPoint(int nr,
    const axisem3d::eigen::DRow2& crds,
    int meshTag,
    std::unique_ptr<const Mass>& mass,
    const TimeScheme& timeScheme) : Point(nr, crds, meshTag, mass) {
  // fields
  point_time::createFields(*this, timeScheme);
  // mass
  mMass->checkCompatibility(mNr, false);
}

/////////////////////////// time loop ///////////////////////////
// stiff to accel
void
FluidPoint::computeStiffToAccel() {
  // Nyquist
  if (mNr % 2 == 0) {
    mFields.mStiff.bottomRows(1).imag().setZero();
  }

  // store stiffness for delta output
  if (mFields.mDeltaStore.rows() > 0) {
    mFields.mDeltaStore = mFields.mStiff;
  }

  // stiff to accel in-place
  mMass->computeAccel(mFields.mStiff);

  // apply pressure source
  if (mFields.mPressureSource.rows() > 0) {
    // add pressure to new acceleration
    mFields.mStiff += mFields.mPressureSource;
    // zero pressure for the next time step
    mFields.mPressureSource.setZero();
  }

  // store acceleration for pressure output
  if (mFields.mPressureStore.rows() > 0) {
    mFields.mPressureStore = mFields.mStiff;
  }
}
