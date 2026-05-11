// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: MIT
// Copyright (C) 2019 - 2026 by the AxiSEM3D authors
//
// This file is part of the AxiSEM3D library. See the LICENSE file for details.
//
// -----------------------------------------------------------------------------

//  elastic material

#ifndef Elastic_hpp
#define Elastic_hpp

#include "Attenuation.hpp"

class Elastic {
  public:
  // constructor
  Elastic(bool is1D, std::unique_ptr<Attenuation>& attenuation) :
      m1D(is1D), mAttenuation(attenuation.release()) {
    // nothing
  }

  // copy constructor
  Elastic(const Elastic& other) :
      m1D(other.m1D),
      mAttenuation((other.mAttenuation == nullptr) ? nullptr : other.mAttenuation->clone()) {
    // nothing
  }

  // destructor
  virtual ~Elastic() = default;

  // clone for copy constructor
  virtual std::unique_ptr<Elastic>
  clone() const = 0;

  // 1D operation
  bool
  is1D() const {
    return m1D;
  }

  // check compatibility
  virtual void
  checkCompatibility(int nr, bool elemInFourier) const {
    if (mAttenuation) {
      mAttenuation->checkCompatibility(nr, elemInFourier, m1D);
    }
  }

  // reset to zero
  void
  resetToZero() const {
    if (mAttenuation) {
      mAttenuation->resetToZero();
    }
  }

  ///////////////////////// strain to stress /////////////////////////
  // RTZ coordinates
  virtual bool
  inRTZ() const = 0;

  // strain => stress in Fourier space
  virtual void
  strainToStress_FR(const axisem3d::eigen::vec_ar6_CMatPP_RM& strain,
      axisem3d::eigen::vec_ar6_CMatPP_RM& stress,
      int nu_1) const = 0;

  // strain => stress in cardinal space
  virtual void
  strainToStress_CD(
      const axisem3d::eigen::RMatXN6& strain, axisem3d::eigen::RMatXN6& stress, int nr) const = 0;

  ///////////////////////// attenuation /////////////////////////
  protected:
  // attenuation in Fourier space
  void
  applyAttenuation(const axisem3d::eigen::vec_ar6_CMatPP_RM& strain,
      axisem3d::eigen::vec_ar6_CMatPP_RM& stress,
      int nu_1) const {
    if (mAttenuation) {
      mAttenuation->apply(strain, stress, nu_1);
    }
  }

  // attenuation in cardinal space
  void
  applyAttenuation(
      const axisem3d::eigen::RMatXN6& strain, axisem3d::eigen::RMatXN6& stress, int nr) const {
    if (mAttenuation) {
      mAttenuation->apply(strain, stress, nr);
    }
  }

  ////////////////////////// data //////////////////////////
  protected:
  // 1D/3D flag
  const bool m1D;

  private:
  // attenuation
  const std::unique_ptr<Attenuation> mAttenuation;
};

#endif /* Elastic_hpp */
