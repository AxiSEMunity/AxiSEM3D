// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: MIT
// Copyright (C) 2019 - 2026 by the AxiSEM3D authors
//
// This file is part of the AxiSEM3D library. See the LICENSE file for details.
//
// -----------------------------------------------------------------------------

//  eigen for station

#ifndef eigen_station_hpp
#define eigen_station_hpp

#include "eigen_tensor.hpp"
#include "eigen_element.hpp"

namespace axisem3d::eigen {
  // tensor
  typedef Eigen::Tensor<numerical::Real, 3, Eigen::RowMajor> RTensor3;
  typedef Eigen::Tensor<numerical::Real, 1, Eigen::RowMajor> RTensor1;
  typedef Eigen::array<Eigen::DenseIndex, 3> IArray3;
  typedef Eigen::array<Eigen::DenseIndex, 2> IArray2;
  typedef Eigen::array<Eigen::DenseIndex, 1> IArray1;

  // weights
  typedef Eigen::Matrix<numerical::Real, 1, Eigen::Dynamic> RRowX;
  typedef Eigen::Matrix<double, 1, nPEM> DRowN;

  // inplane interpolation
  typedef Eigen::Matrix<numerical::ComplexR, Eigen::Dynamic, 1> CMatX1;
  typedef Eigen::Matrix<numerical::ComplexR, Eigen::Dynamic, 3> CMatX3;
  typedef Eigen::Matrix<numerical::ComplexR, Eigen::Dynamic, 6> CMatX6;
  typedef Eigen::Matrix<numerical::ComplexR, Eigen::Dynamic, 9> CMatX9;

  // Fourier interpolation
  typedef Eigen::Matrix<numerical::Real, 1, 1> RRow1;
  typedef Eigen::Matrix<numerical::Real, 1, 3> RRow3;
  typedef Eigen::Matrix<numerical::Real, 1, 6> RRow6;
  typedef Eigen::Matrix<numerical::Real, 1, 9> RRow9;

  // buffers
  // alawys use double for time output
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> DColX;
  typedef Eigen::Matrix<numerical::Real, Eigen::Dynamic, 1> RColX;
  typedef Eigen::Matrix<numerical::ComplexR, Eigen::Dynamic, 1> CColX;
  typedef Eigen::Matrix<Real, Eigen::Dynamic, 1> RMatX1_RM;
  typedef Eigen::Matrix<Real, Eigen::Dynamic, 3, Eigen::RowMajor> RMatX3_RM;
  typedef Eigen::Matrix<Real, Eigen::Dynamic, 6, Eigen::RowMajor> RMatX6_RM;
  typedef Eigen::Matrix<Real, Eigen::Dynamic, 9, Eigen::RowMajor> RMatX9_RM;
} // namespace axisem3d::eigen

#endif /* eigen_station_hpp */
