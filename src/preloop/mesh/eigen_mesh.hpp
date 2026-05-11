// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: MIT
// Copyright (C) 2019 - 2026 by the AxiSEM3D authors
//
// This file is part of the AxiSEM3D library. See the LICENSE file for details.
//
// -----------------------------------------------------------------------------

//  eigen typedef for mesh

#ifndef eigen_mesh_hpp
#define eigen_mesh_hpp

#include "eigen_generic.hpp"
#include "spectral.hpp"

namespace axisem3d::eigen {
  using Eigen::Dynamic;
  using Eigen::RowMajor;

  ///////// nodal-level /////////
  // connectivity
  typedef Eigen::Matrix<int, Dynamic, 4, RowMajor> IMatX4_RM;
  // nodal coords
  typedef Eigen::Matrix<double, Dynamic, 2, RowMajor> DMatX2_RM;

  ///////// GLL-level /////////
  // GLL tag on elements
  typedef Eigen::Matrix<int, Dynamic, spectral::nPEM, RowMajor> IMatXN_RM;
} // namespace axisem3d::eigen

#endif /* eigen_mesh_hpp */
