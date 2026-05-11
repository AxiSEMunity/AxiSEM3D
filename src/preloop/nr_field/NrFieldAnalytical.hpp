// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: MIT
// Copyright (C) 2019 - 2026 by the AxiSEM3D authors
//
// This file is part of the AxiSEM3D library. See the LICENSE file for details.
//
// -----------------------------------------------------------------------------

//  analytical Nr(s,z)

#ifndef NrFieldAnalytical_hpp
#define NrFieldAnalytical_hpp

#include "NrField.hpp"
#include <vector>

class NrFieldAnalytical : public NrField {
  public:
  // constructor
  NrFieldAnalytical();

  // get nr by (s, z)
  axisem3d::eigen::IColX
  getNrAtPoints(const axisem3d::eigen::DMatX2_RM& sz) const;

  // verbose
  std::string
  verbose() const;

  private:
  // code ID
  static std::string sCodeID;

  // TODO: add your data here
  // below are data for
  // sCodeID = "depth-dependent (AxiSEM3D default)"
  std::vector<double> mControlDepths;
  std::vector<double> mControlNrs;
};

#endif /* NrFieldAnalytical_hpp */
