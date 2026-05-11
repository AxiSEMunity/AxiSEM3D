// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: MIT
// Copyright (C) 2019 - 2026 by the AxiSEM3D authors
//
// This file is part of the AxiSEM3D library. See the LICENSE file for details.
//
// -----------------------------------------------------------------------------

//  eigen typedef for point

#ifndef eigen_point_hpp
#define eigen_point_hpp

#include "eigen.hpp"
#include "numerical.hpp"

namespace axisem3d::eigen {
  // coords
  typedef Eigen::Matrix<double, 1, 2> DRow2;

  // fields on a fluid point
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> DColX;
  typedef Eigen::Matrix<numerical::Real, Eigen::Dynamic, 1> RColX;
  typedef Eigen::Matrix<numerical::ComplexR, Eigen::Dynamic, 1> CColX;

  // fields on a solid point
  typedef Eigen::Matrix<double, Eigen::Dynamic, 3> DMatX3;
  typedef Eigen::Matrix<numerical::Real, Eigen::Dynamic, 3> RMatX3;
  typedef Eigen::Matrix<numerical::ComplexR, Eigen::Dynamic, 3> CMatX3;
} // namespace axisem3d::eigen

#endif /* eigen_point_hpp */
