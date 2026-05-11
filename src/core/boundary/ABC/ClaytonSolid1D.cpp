// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: MIT
// Copyright (C) 2019 - 2026 by the AxiSEM3D authors
//
// This file is part of the AxiSEM3D library. See the LICENSE file for details.
//
// -----------------------------------------------------------------------------

//  Clayton-Enquist ABC for solid points in 1D
//  theta: angle between surface normal and z-axis

#include "ClaytonSolid1D.hpp"
#include "SolidPoint.hpp"

// apply ABC
void
ClaytonSolid1D::apply() const {
  // get fields
  const axisem3d::eigen::CMatX3& veloc = mSolidPoint->getFields().mVeloc;
  axisem3d::eigen::CMatX3& stiff = mSolidPoint->getFields().mStiff;

  // s, z
  stiff.col(0) -= (mRSA_CosT2_p_RPA_SinT2 * veloc.col(0) + mRPA_m_RSA_x_CosT_SinT * veloc.col(2));
  stiff.col(2) -= (mRPA_m_RSA_x_CosT_SinT * veloc.col(0) + mRSA_SinT2_p_RPA_CosT2 * veloc.col(2));

  // phi
  stiff.col(1) -= mRSA * veloc.col(1);
}
