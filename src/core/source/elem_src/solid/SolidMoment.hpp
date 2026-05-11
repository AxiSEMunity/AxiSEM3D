// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: MIT
// Copyright (C) 2019 - 2026 by the AxiSEM3D authors
//
// This file is part of the AxiSEM3D library. See the LICENSE file for details.
//
// -----------------------------------------------------------------------------

//  moment source on solid element

#ifndef SolidMoment_hpp
#define SolidMoment_hpp

#include "SolidSource.hpp"
#include "eigen_element.hpp"

class SolidMoment : public SolidSource {
  public:
  // constructor
  SolidMoment(std::unique_ptr<STF>& stf,
      const std::shared_ptr<SolidElement>& element,
      const axisem3d::eigen::CMatXN6& pattern);

  // apply source at a time step
  void
  apply(double time) const;

  private:
  // source pattern
  const axisem3d::eigen::CMatXN6 mPattern;

  ////////////////////////////////////////
  //////////////// static ////////////////
  ////////////////////////////////////////

  // workspace
  inline static axisem3d::eigen::CMatXN6 sPattern = axisem3d::eigen::CMatXN6(0, spectral::nPEM * 6);
};

#endif /* SolidMoment_hpp */
