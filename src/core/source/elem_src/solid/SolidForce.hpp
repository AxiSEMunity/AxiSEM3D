// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: MIT
// Copyright (C) 2019 - 2026 by the AxiSEM3D authors
//
// This file is part of the AxiSEM3D library. See the LICENSE file for details.
//
// -----------------------------------------------------------------------------

//  force source on solid element

#ifndef SolidForce_hpp
#define SolidForce_hpp

#include "SolidSource.hpp"
#include "eigen_element.hpp"

class SolidForce : public SolidSource {
  public:
  // constructor
  SolidForce(std::unique_ptr<STF>& stf,
      const std::shared_ptr<SolidElement>& element,
      const axisem3d::eigen::CMatXN3& pattern);

  // apply source at a time step
  void
  apply(double time) const;

  private:
  // source pattern
  const axisem3d::eigen::CMatXN3 mPattern;

  ////////////////////////////////////////
  //////////////// static ////////////////
  ////////////////////////////////////////

  // workspace
  inline static axisem3d::eigen::CMatXN3 sPattern = axisem3d::eigen::CMatXN3(0, spectral::nPEM * 3);
};

#endif /* SolidForce_hpp */
