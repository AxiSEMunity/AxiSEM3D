// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: MIT
// Copyright (C) 2019 - 2026 by the AxiSEM3D authors
//
// This file is part of the AxiSEM3D library. See the LICENSE file for details.
//
// -----------------------------------------------------------------------------

//  pressure source on fluid element

#include "FluidPressure.hpp"
#include "FluidElement.hpp"

// constructor
FluidPressure::FluidPressure(std::unique_ptr<STF>& stf,
    const std::shared_ptr<FluidElement>& element,
    const axisem3d::eigen::CMatXN& pattern) : FluidSource(stf, element), mPattern(pattern) {
  // prepare
  element->preparePressureSource();

  // workspace
  if (sPattern.rows() < mPattern.rows()) {
    sPattern.resize(mPattern.rows(), spectral::nPEM);
  }
}

// apply source at a time step
void
FluidPressure::apply(double time) const {
  int nu_1 = (int)mPattern.rows();
  sPattern.topRows(nu_1) = mPattern * mSTF->getValue(time);
  mElement->addPressureSource(sPattern, nu_1);
}
