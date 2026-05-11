// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: MIT
// Copyright (C) 2019 - 2026 by the AxiSEM3D authors
//
// This file is part of the AxiSEM3D library. See the LICENSE file for details.
//
// -----------------------------------------------------------------------------

//  solid-fluid boundary condition in 3D

#ifndef SolidFluidCoupling3D_hpp
#define SolidFluidCoupling3D_hpp

#include "SolidFluidCoupling.hpp"

class SolidFluidCoupling3D : public SolidFluidCoupling {
  public:
  // constructor
  SolidFluidCoupling3D(const std::shared_ptr<SolidPoint>& sp,
      const std::shared_ptr<FluidPoint>& fp,
      const axisem3d::eigen::DMatX3& n_unassmb,
      const axisem3d::eigen::DMatX3& n_assmb,
      const axisem3d::eigen::DColX& massFluid);

  private:
  // check compatibility
  void
  checkCompatibility(int nr) const;

  public:
  // solid => fluid
  void
  coupleSolidToFluid(
      const axisem3d::eigen::CMatX3& solidDispl, axisem3d::eigen::CColX& fluidStiff) const;

  // fluid => solid
  void
  coupleFluidToSolid(
      const axisem3d::eigen::CColX& fluidStiff, axisem3d::eigen::CMatX3& solidStiff) const;

  private:
  // These two normal vectors enable isochronous MPI communication for solid
  // and fluid domains. Though it is bad practice to mix MPI and physics,
  // but this trick can lead to significant performance boost.
  const axisem3d::eigen::RMatX3 mNormal_UnassembledMPI;
  const axisem3d::eigen::RMatX3 mNormal_AssembledMPI_InvMassFluid;

  ////////////////////////////////////////
  //////////////// static ////////////////
  ////////////////////////////////////////

  // workspace
  inline static axisem3d::eigen::RMatX3 sSolidR = axisem3d::eigen::RMatX3(0, 3);
  inline static axisem3d::eigen::CMatX3 sSolidC = axisem3d::eigen::CMatX3(0, 3);
  inline static axisem3d::eigen::RColX sFluidR = axisem3d::eigen::RColX(0);
  inline static axisem3d::eigen::CColX sFluidC = axisem3d::eigen::CColX(0);
};

#endif /* SolidFluidCoupling3D_hpp */
