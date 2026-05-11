// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: MIT
// Copyright (C) 2019 - 2026 by the AxiSEM3D authors
//
// This file is part of the AxiSEM3D library. See the LICENSE file for details.
//
// -----------------------------------------------------------------------------

//  constant Nr(s,z)

#include "NrFieldConstant.hpp"
#include "bstring.hpp"

// get nr by (s, z)
axisem3d::eigen::IColX
NrFieldConstant::getNrAtPoints(const axisem3d::eigen::DMatX2_RM& sz) const {
  return axisem3d::eigen::IColX::Constant(sz.rows(), mNr);
}

// verbose
std::string
NrFieldConstant::verbose() const {
  std::stringstream ss;
  ss << bstring::boxTitle("Nr(s,z)");
  ss << bstring::boxEquals(0, 5, "type", "CONSTANT");
  ss << bstring::boxEquals(0, 5, "value", mNr);
  ss << bstring::boxBaseline() << "\n\n";
  return ss.str();
}
