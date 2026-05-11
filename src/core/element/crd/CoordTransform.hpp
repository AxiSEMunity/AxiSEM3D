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

#ifndef CoordTransform_hpp
#define CoordTransform_hpp

#include "eigen_element.hpp"

class CoordTransform {
  public:
  // destructor
  virtual ~CoordTransform() = default;

  // (s,phi,z) -> (R,T,Z)
  virtual void
  transformSPZ_RTZ3(axisem3d::eigen::vec_ar3_CMatPP_RM& ui, int nu_1) const = 0;

  // (R,T,Z) -> (s,phi,z)
  virtual void
  transformRTZ_SPZ3(axisem3d::eigen::vec_ar3_CMatPP_RM& ui, int nu_1) const = 0;

  // (s,phi,z) -> (R,T,Z) for nabla
  virtual void
  transformSPZ_RTZ9(axisem3d::eigen::vec_ar9_CMatPP_RM& nij, int nu_1) const = 0;

  // (R,T,Z) -> (s,phi,z) for nabla
  virtual void
  transformRTZ_SPZ9(axisem3d::eigen::vec_ar9_CMatPP_RM& nij, int nu_1) const = 0;

  // (s,phi,z) -> (R,T,Z) for Voigt
  virtual void
  transformSPZ_RTZ6(axisem3d::eigen::vec_ar6_CMatPP_RM& eij, int nu_1) const = 0;

  // (R,T,Z) -> (s,phi,z) for Voigt
  virtual void
  transformRTZ_SPZ6(axisem3d::eigen::vec_ar6_CMatPP_RM& sij, int nu_1) const = 0;
};

#endif /* CoordTransform_hpp */
