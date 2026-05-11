// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: MIT
// Copyright (C) 2019 - 2026 by the AxiSEM3D authors
//
// This file is part of the AxiSEM3D library. See the LICENSE file for details.
//
// -----------------------------------------------------------------------------

//  coordinate transform between (s,phi,z) and (R,T,Z)
//  for vector and 2nd-order tensor fields
//  Nothing for Cartesian meshes.

#ifndef CoordTransformCartesian_hpp
#define CoordTransformCartesian_hpp

#include "CoordTransform.hpp"

class CoordTransformCartesian : public CoordTransform {
  public:
  // (s,phi,z) -> (R,T,Z)
  void
  transformSPZ_RTZ3(axisem3d::eigen::vec_ar3_CMatPP_RM& ui, int nu_1) const {
    // nothing
  }

  // (R,T,Z) -> (s,phi,z)
  void
  transformRTZ_SPZ3(axisem3d::eigen::vec_ar3_CMatPP_RM& ui, int nu_1) const {
    // nothing
  }

  // (s,phi,z) -> (R,T,Z) for nabla
  void
  transformSPZ_RTZ9(axisem3d::eigen::vec_ar9_CMatPP_RM& nij, int nu_1) const {
    // nothing
  }

  // (R,T,Z) -> (s,phi,z) for nabla
  void
  transformRTZ_SPZ9(axisem3d::eigen::vec_ar9_CMatPP_RM& nij, int nu_1) const {
    // nothing
  }

  // (s,phi,z) -> (R,T,Z) for Voigt
  void
  transformSPZ_RTZ6(axisem3d::eigen::vec_ar6_CMatPP_RM& eij, int nu_1) const {
    // nothing
  }

  // (R,T,Z) -> (s,phi,z) for Voigt
  void
  transformRTZ_SPZ6(axisem3d::eigen::vec_ar6_CMatPP_RM& sij, int nu_1) const {
    // nothing
  }
};

#endif /* CoordTransformCartesian_hpp */
