// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: MIT
// Copyright (C) 2019 - 2026 by the AxiSEM3D authors
//
// This file is part of the AxiSEM3D library. See the LICENSE file for details.
//
// -----------------------------------------------------------------------------

//  vertical undulation
//  generator of PRT in core

#ifndef Undulation_hpp
#define Undulation_hpp

// data
#include "PhysicalProperty.hpp"

// gradient
#include "SolverFFTW.hpp"
class Quad;
namespace axisem3d::eigen {
  using numerical::ComplexD;
  using spectral::nPED;
  using spectral::nPEM;
  typedef Eigen::Matrix<ComplexD, nPED, nPED, Eigen::RowMajor> ZMatPP_RM;
  typedef std::vector<std::array<ZMatPP_RM, 1>> vec_ar1_ZMatPP_RM;
  typedef std::vector<std::array<ZMatPP_RM, 3>> vec_ar3_ZMatPP_RM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, nPEM * 3> DMatXN3;
} // namespace axisem3d::eigen

// release
#include "PRT.hpp"

class Undulation {
  public:
  // add undulation
  void
  addUndulation(const axisem3d::eigen::arN_DColX& und) {
    mDeltaZ.addGLL(und);
  }

  // get elemental
  axisem3d::eigen::DMatXN
  getElemental() const {
    return mDeltaZ.getElemental();
  }

  // get pointwise
  axisem3d::eigen::arN_DColX
  getPointwise() const {
    return mDeltaZ.getPointwise();
  }

  // finishing 3D properties
  void
  finishing3D() const;

  // finished 3D properties
  void
  finished3D(const Quad& myQuad);

  // get Jacobian for mass
  axisem3d::eigen::arN_DColX
  getMassJacobian(const axisem3d::eigen::DMat2N& sz) const;

  // create PRT
  std::unique_ptr<const PRT>
  createPRT(const axisem3d::eigen::DMat2N& sz) const;

  // compute 3D normal at a point
  axisem3d::eigen::DMatX3
  computeNormal3D(
      const axisem3d::eigen::DCol2& n1D, const axisem3d::eigen::DMat2N& sz, int ipnt) const;

  ///////////////////////// data /////////////////////////
  private:
  // delta Z
  PhysicalProperty<spectral::nPEM> mDeltaZ;

  // gradient of delta Z
  axisem3d::eigen::DMatXN3 mDeltaZ_RTZ;

  ///////////////////////// static /////////////////////////
  public:
  // finished 3D properties
  static void
  finished3D();

  private:
  // static fft variables
  inline static SolverFFTW<double, spectral::nPEM> sFFT_N1;
  inline static SolverFFTW<double, spectral::nPEM * 3> sFFT_N3;
  inline static axisem3d::eigen::vec_ar1_ZMatPP_RM sDeltaZ_Fourier;
  inline static axisem3d::eigen::vec_ar3_ZMatPP_RM sDeltaZ_SPZ_Fourier;
};

#endif /* Undulation_hpp */
