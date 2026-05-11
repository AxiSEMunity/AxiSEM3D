// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: MIT
// Copyright (C) 2019 - 2026 by the AxiSEM3D authors
//
// This file is part of the AxiSEM3D library. See the LICENSE file for details.
//
// -----------------------------------------------------------------------------

//  constant Nr(s,z)

#ifndef NrFieldConstant_hpp
#define NrFieldConstant_hpp

#include "NrField.hpp"

class NrFieldConstant : public NrField {
  public:
  // constructor
  NrFieldConstant(int nr) : mNr(nr) {
    // nothing
  }

  // get nr by (s, z)
  axisem3d::eigen::IColX
  getNrAtPoints(const axisem3d::eigen::DMatX2_RM& sz) const;

  // verbose
  std::string
  verbose() const;

  private:
  // the constant nr value
  const int mNr;
};

#endif /* NrFieldConstant_hpp */
