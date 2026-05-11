// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: MIT
// Copyright (C) 2019 - 2026 by the AxiSEM3D authors
//
// This file is part of the AxiSEM3D library. See the LICENSE file for details.
//
// -----------------------------------------------------------------------------

//  attenuation builder
//  generator of Attenuation in core

#ifndef AttBuilder_hpp
#define AttBuilder_hpp

#include "eigen_sem.hpp"
#include "Attenuation.hpp"
class ExodusMesh;

class AttBuilder {
  public:
  // constructor
  AttBuilder(const ExodusMesh& exodusMesh, bool cg4, double dt);

  // verbose
  static std::string
  verbose(const std::unique_ptr<const AttBuilder>& attBuilder);

  // create attenuation
  std::unique_ptr<Attenuation>
  createAttenuation(const axisem3d::eigen::DMatXN& QKp,
      const axisem3d::eigen::DMatXN& QMu,
      axisem3d::eigen::DMatXN& kp,
      axisem3d::eigen::DMatXN& mu,
      const axisem3d::eigen::DRow4& weightsCG4,
      bool elastic1D) const;

  private:
  // set alpha, beta, gamma
  // NOTE: with SPECFEM legency deprecated, alpha, beta and gamma
  //       all become element-independent
  void
  setAlphaBetaGamma() const;

  private:
  // from Exodus
  const double mFmin, mFmax, mFref;
  const int mNSLS;
  axisem3d::eigen::DColX mW, mY;

  // from inparam
  const bool mUseCG4;

  // dt
  const double mDt;
};

#endif /* AttBuilder_hpp */
