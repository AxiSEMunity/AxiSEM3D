// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: MIT
// Copyright (C) 2019 - 2026 by the AxiSEM3D authors
//
// This file is part of the AxiSEM3D library. See the LICENSE file for details.
//
// -----------------------------------------------------------------------------

//  base class of Nr(s,z)

#ifndef NrField_hpp
#define NrField_hpp

#include "eigen_generic.hpp"
#include <memory>

// need distance tolerance for NrFieldPointwise
class ExodusMesh;

namespace axisem3d::eigen {
  // coords
  typedef Eigen::Matrix<double, Eigen::Dynamic, 2, Eigen::RowMajor> DMatX2_RM;
} // namespace axisem3d::eigen

class NrField {
  public:
  // destructor
  virtual ~NrField() = default;

  // get nr by (s, z)
  virtual axisem3d::eigen::IColX
  getNrAtPoints(const axisem3d::eigen::DMatX2_RM& sz) const = 0;

  // verbose
  virtual std::string
  verbose() const = 0;

  ////////////////////////////// static //////////////////////////////
  public:
  // build from inparam
  static std::unique_ptr<const NrField>
  buildInparam(const ExodusMesh& exodusMesh);

  // is lucky number
  static bool
  isLuckyNumber(int n);

  // next lucky number
  static int
  nextLuckyNumber(int n);

  private:
  // lucky numbers
  static std::vector<int> sLuckyNumbers;
};

#endif /* NrField_hpp */
