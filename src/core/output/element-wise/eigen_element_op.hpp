// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: MIT
// Copyright (C) 2019 - 2026 by the AxiSEM3D authors
//
// This file is part of the AxiSEM3D library. See the LICENSE file for details.
//
// -----------------------------------------------------------------------------

//  eigen for element output

#ifndef eigen_element_op_hpp
#define eigen_element_op_hpp

#include "eigen_station.hpp"
#include "eigen_generic.hpp"

namespace axisem3d::eigen {
  using Eigen::Dynamic;
  using numerical::Real;
  using spectral::nPEM;

  // tensor
  typedef Eigen::Tensor<numerical::Real, 5, Eigen::RowMajor> RTensor5;
  typedef Eigen::Tensor<numerical::Real, 4, Eigen::RowMajor> RTensor4;
  typedef Eigen::array<Eigen::DenseIndex, 5> IArray5;
  typedef Eigen::array<Eigen::DenseIndex, 4> IArray4;

  // element-na info
  typedef Eigen::Matrix<int, Eigen::Dynamic, 4, Eigen::RowMajor> IMatX4_RM;
  typedef Eigen::Matrix<int, Eigen::Dynamic, 5, Eigen::RowMajor> IMatX5_RM;

  // azimuthal making real
  typedef Eigen::Matrix<Real, Dynamic, nPEM, Eigen::RowMajor> RMatXN_RM;
  typedef Eigen::Matrix<Real, Dynamic, nPEM * 3, Eigen::RowMajor> RMatXN3_RM;
  typedef Eigen::Matrix<Real, Dynamic, nPEM * 6, Eigen::RowMajor> RMatXN6_RM;
  typedef Eigen::Matrix<Real, Dynamic, nPEM * 9, Eigen::RowMajor> RMatXN9_RM;
} // namespace axisem3d::eigen

#endif /* eigen_element_op_hpp */
