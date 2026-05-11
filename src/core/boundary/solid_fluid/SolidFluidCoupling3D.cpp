// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: MIT
// Copyright (C) 2019 - 2026 by the AxiSEM3D authors
//
// This file is part of the AxiSEM3D library. See the LICENSE file for details.
//
// -----------------------------------------------------------------------------

//  solid-fluid boundary condition in 3D

#include "SolidFluidCoupling3D.hpp"
#include "SolidPoint.hpp"
#include "fft.hpp"

// constructor
SolidFluidCoupling3D::SolidFluidCoupling3D(const std::shared_ptr<SolidPoint>& sp,
    const std::shared_ptr<FluidPoint>& fp,
    const axisem3d::eigen::DMatX3& n_unassmb,
    const axisem3d::eigen::DMatX3& n_assmb,
    const axisem3d::eigen::DColX& massFluid) :
    SolidFluidCoupling(sp, fp), mNormal_UnassembledMPI(n_unassmb.cast<numerical::Real>()),
    mNormal_AssembledMPI_InvMassFluid(
        (massFluid.cwiseInverse().asDiagonal() * n_assmb).cast<numerical::Real>()) {
  checkCompatibility(mSolidPoint->getNr());
}

// check compatibility
void
SolidFluidCoupling3D::checkCompatibility(int nr) const {
  // check size
  SolidFluidCoupling::checkCompatibility(nr);
  if (mNormal_UnassembledMPI.rows() != nr) {
    throw std::runtime_error("SolidFluidCoupling3D::checkCompatibility || "
                             "Incompatible sizes.");
  }

  // expand workspace if needed
  if (sSolidR.rows() < nr) {
    // real
    sSolidR.resize(nr, 3);
    sFluidR.resize(nr);
    // complex
    int nc = nr / 2 + 1;
    sSolidC.resize(nc, 3);
    sFluidC.resize(nc);
  }

  // report request to FFT
  fft::gFFT_1.addNR(nr);
  fft::gFFT_3.addNR(nr);
}

// solid => fluid
void
SolidFluidCoupling3D::coupleSolidToFluid(
    const axisem3d::eigen::CMatX3& solidDispl, axisem3d::eigen::CColX& fluidStiff) const {
  // constants
  int nr = (int)mNormal_UnassembledMPI.rows();
  int nu_1 = nr / 2 + 1;

  // FFT: Fourier => cardinal
  fft::gFFT_3.computeC2R(solidDispl, sSolidR, nr);

  // fluid => solid coupling term
  sFluidR.topRows(nr) = (mNormal_UnassembledMPI.cwiseProduct(sSolidR.topRows(nr)).rowwise().sum());

  // FFT: cardinal => Fourier
  fft::gFFT_1.computeR2C(sFluidR, sFluidC, nr);

  // add coupling term
  fluidStiff += sFluidC.topRows(nu_1);
}

// fluid => solid
void
SolidFluidCoupling3D::coupleFluidToSolid(
    const axisem3d::eigen::CColX& fluidStiff, axisem3d::eigen::CMatX3& solidStiff) const {
  // constants
  int nr = (int)mNormal_UnassembledMPI.rows();
  int nu_1 = nr / 2 + 1;

  // FFT: Fourier => cardinal
  fft::gFFT_1.computeC2R(fluidStiff, sFluidR, nr);

  // fluid => solid coupling term
  sSolidR.topRows(nr) = (sFluidR.topRows(nr).asDiagonal() * mNormal_AssembledMPI_InvMassFluid);

  // FFT: cardinal => Fourier
  fft::gFFT_3.computeR2C(sSolidR, sSolidC, nr);

  // add coupling term
  solidStiff -= sSolidC.topRows(nu_1);
}
