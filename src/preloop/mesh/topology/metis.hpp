// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: MIT
// Copyright (C) 2019 - 2026 by the AxiSEM3D authors
//
// This file is part of the AxiSEM3D library. See the LICENSE file for details.
//
// -----------------------------------------------------------------------------

//  dual graph partioning by Metis
//  compatible with both 32-bit and 64-bit builds of Metis

#ifndef metis_hpp
#define metis_hpp

#include "eigen_mesh.hpp"

namespace axisem3d::metis {
  // form neighbourhood of connectivity
  void
  formNeighbourhood(const axisem3d::eigen::IMatX4_RM& connectivity,
      int ncommon,
      std::vector<axisem3d::eigen::IColX>& neighbours);

  // domain decomposition
  double
  decompose(const axisem3d::eigen::IMatX4_RM& connectivity,
      const axisem3d::eigen::IColX& solid_fluid,
      const axisem3d::eigen::DColX& weights,
      int npart,
      int rseed,
      axisem3d::eigen::IColX& elemRank);
} // namespace axisem3d::metis

#endif /* metis_hpp */
