// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: MIT
// Copyright (C) 2019 - 2026 by the AxiSEM3D authors
//
// This file is part of the AxiSEM3D library. See the LICENSE file for details.
//
// -----------------------------------------------------------------------------

//  solid-fluid boundary condition

#ifndef SolidFluidCoupling_hpp
#define SolidFluidCoupling_hpp

// point
#include "eigen_point.hpp"
#include <memory>
class SolidPoint;
class FluidPoint;

class SolidFluidCoupling {
  public:
  // constructor
  SolidFluidCoupling(const std::shared_ptr<SolidPoint>& sp, const std::shared_ptr<FluidPoint>& fp);

  // destructor
  virtual ~SolidFluidCoupling() = default;

  // get solid point
  const std::shared_ptr<SolidPoint>&
  getSolidPoint() const {
    return mSolidPoint;
  }

  // apply coupling
  void
  apply() const;

  ////////////////////////////// virtual //////////////////////////////
  protected:
  // check compatibility
  virtual void
  checkCompatibility(int nr) const;

  public:
  // solid => fluid
  virtual void
  coupleSolidToFluid(
      const axisem3d::eigen::CMatX3& solidDispl, axisem3d::eigen::CColX& fluidStiff) const = 0;

  // fluid => solid
  virtual void
  coupleFluidToSolid(
      const axisem3d::eigen::CColX& fluidStiff, axisem3d::eigen::CMatX3& solidStiff) const = 0;

  protected:
  // coupled solid-fluid pair
  const std::shared_ptr<SolidPoint> mSolidPoint;
  const std::shared_ptr<FluidPoint> mFluidPoint;
};

#endif /* SolidFluidCoupling_hpp */
