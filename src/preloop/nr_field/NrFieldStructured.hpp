// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: MIT
// Copyright (C) 2019 - 2026 by the AxiSEM3D authors
//
// This file is part of the AxiSEM3D library. See the LICENSE file for details.
//
// -----------------------------------------------------------------------------

//  structured Nr(s,z)

#ifndef NrFieldStructured_hpp
#define NrFieldStructured_hpp

#include "NrField.hpp"
#include "StructuredGrid.hpp"

class NrFieldStructured : public NrField {
  public:
  // constructor
  NrFieldStructured(const std::string& fname, int valOutOfRange);

  // get nr by (s, z)
  axisem3d::eigen::IColX
  getNrAtPoints(const axisem3d::eigen::DMatX2_RM& sz) const;

  // verbose
  std::string
  verbose() const;

  private:
  // file name
  const std::string mFilename;

  // factor
  const int mValueOutOfRange;

  // grid
  std::unique_ptr<const StructuredGrid<2, int>> mGrid = nullptr;
};

#endif /* NrFieldStructured_hpp */
