// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: MIT
// Copyright (C) 2019 - 2026 by the AxiSEM3D authors
//
// This file is part of the AxiSEM3D library. See the LICENSE file for details.
//
// -----------------------------------------------------------------------------

//  ellipticity

#include "Ellipticity.hpp"
#include "sg_tools.hpp"

// get undulation on points
bool
Ellipticity::getUndulation(
    const axisem3d::eigen::DMatX3& spz, axisem3d::eigen::DColX& undulation) const {
  if (geodesy::isCartesian() || geodesy::getOuterFlattening() < numerical::dEpsilon) {
    return false;
  }
  // theta, phi, r
  const axisem3d::eigen::DMatX3& llr =
      coordsFromMeshToModel(spz, false, false, false, false, false, false, mModelName);
  const axisem3d::eigen::DMatX3& tpr = geodesy::llr2tpr(llr, false);
  // flattening
  typedef Eigen::Array<double, Eigen::Dynamic, 1> DColX_Array;
  const DColX_Array& r = tpr.col(2);
  const DColX_Array& t = tpr.col(0);
  const DColX_Array& f = geodesy::computeFlattening(r);
  // a and b
  const DColX_Array& b = (1. - f).pow(2. / 3.) * r;
  const DColX_Array& a = b / (1. - f);
  const DColX_Array& tmp = (a * t.cos()).square() + (b * t.sin()).square();
  undulation = a * b / tmp.sqrt().max(numerical::dEpsilon) - r;
  return true;
}

// verbose
std::string
Ellipticity::verbose() const {
  std::stringstream ss;
  ss << bstring::boxSubTitle(0, mModelName + " ", '~');
  ss << bstring::boxEquals(2, 21, "class name", "Ellipticity");
  ss << bstring::boxEquals(2, 21, "model scope", "global");
  ss << bstring::boxEquals(2, 21, "flattening on surface", geodesy::getOuterFlattening());
  return ss.str();
}
