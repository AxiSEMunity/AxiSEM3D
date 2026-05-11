// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: MIT
// Copyright (C) 2019 - 2026 by the AxiSEM3D authors
//
// This file is part of the AxiSEM3D library. See the LICENSE file for details.
//
// -----------------------------------------------------------------------------

//  linear mapping

#ifndef MappingLinear_hpp
#define MappingLinear_hpp

#include "Mapping.hpp"

class MappingLinear : public Mapping {
  public:
  // constructor
  MappingLinear(const axisem3d::eigen::DMat24& nodalSZ) : Mapping(nodalSZ) {
    // undefined curved outer
    mCurvedOuter = -1;
  }

  // forward mapping: (ξ,η) -> (s,z)
  axisem3d::eigen::DCol2
  mapping(const axisem3d::eigen::DCol2& xieta) const {
    // shape function
    axisem3d::eigen::DRow4 shp;
    double xip = 1. + xieta(0);
    double xim = 1. - xieta(0);
    double etp = 1. + xieta(1);
    double etm = 1. - xieta(1);
    shp(0) = xim * etm;
    shp(1) = xip * etm;
    shp(2) = xip * etp;
    shp(3) = xim * etp;
    return .25 * (mNodalSZ * shp.transpose());
  }

  // Jacobian: ∂(s,z) / ∂(ξ,η)
  axisem3d::eigen::DMat22
  jacobian(const axisem3d::eigen::DCol2& xieta) const {
    // derivative of shape function
    axisem3d::eigen::DMat24 dshp;
    double xip = 1. + xieta(0);
    double xim = 1. - xieta(0);
    double etp = 1. + xieta(1);
    double etm = 1. - xieta(1);
    dshp(0, 0) = -etm;
    dshp(0, 1) = etm;
    dshp(0, 2) = etp;
    dshp(0, 3) = -etp;
    dshp(1, 0) = -xim;
    dshp(1, 1) = -xip;
    dshp(1, 2) = xip;
    dshp(1, 3) = xim;
    return .25 * (mNodalSZ * dshp.transpose());
  }
};

#endif /* MappingLinear_hpp */
