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

#include "MassOceanLoad3D.hpp"
#include "fft.hpp"

// constructor
MassOceanLoad3D::MassOceanLoad3D(const axisem3d::eigen::DColX& mass,
    const axisem3d::eigen::DColX& massOcean,
    const axisem3d::eigen::DMatX3& unitNormal) :
    mIM(mass.cwiseInverse().cast<numerical::Real>()),
    mK((massOcean.cwiseQuotient(mass.cwiseProduct(mass + massOcean)).cwiseSqrt().asDiagonal() *
        unitNormal)
            .cast<numerical::Real>()) {
  // nothing
  // im = 1 / m
  // k = sqrt(m0 / [m (m + m0)]) n
}

// check compatibility
void
MassOceanLoad3D::checkCompatibility(int nr, bool solid) const {
  // must on solid
  if (!solid) {
    throw std::runtime_error("MassOceanLoad3D::checkCompatibility || "
                             "Incompatible types: "
                             "ocean load on fluid point.");
  }

  // check size
  if (mIM.rows() != nr) {
    throw std::runtime_error("MassOceanLoad3D::checkCompatibility || "
                             "Incompatible sizes.");
  }

  // expand workspace if needed
  if (sF.rows() < nr) {
    sF.resize(nr, 3);
    sA.resize(nr, 3);
  }

  // report request to FFT
  fft::gFFT_3.addNR(nr);
}

// compute accel in-place for solid
void
MassOceanLoad3D::computeAccel(axisem3d::eigen::CMatX3& stiff3) const {
  // constants
  int nr = (int)mIM.rows();

  // FFT: Fourier => cardinal
  fft::gFFT_3.computeC2R(stiff3, sF, nr);

  // a = im F
  sA.topRows(nr) = mIM.asDiagonal() * sF.topRows(nr);

  // a -= F.k k
  sA.topRows(nr) -= (sF.topRows(nr).cwiseProduct(mK).rowwise().sum().asDiagonal() * mK);

  // FFT: cardinal => Fourier
  fft::gFFT_3.computeR2C(sA, stiff3, nr);
}
