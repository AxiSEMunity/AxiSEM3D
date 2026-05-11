// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: MIT
// Copyright (C) 2019 - 2026 by the AxiSEM3D authors
//
// This file is part of the AxiSEM3D library. See the LICENSE file for details.
//
// -----------------------------------------------------------------------------

//  ellipticity

#ifndef Ellipticity_hpp
#define Ellipticity_hpp

#include "Geometric3D.hpp"

class Ellipticity : public Geometric3D {
  public:
  // constructor
  Ellipticity(const std::string& modelName) : Geometric3D(modelName) {
    // nothing
  }

  private:
  // get undulation on points
  bool
  getUndulation(const axisem3d::eigen::DMatX3& spz, axisem3d::eigen::DColX& undulation) const;

  // verbose
  std::string
  verbose() const;
};

#endif /* Ellipticity_hpp */
