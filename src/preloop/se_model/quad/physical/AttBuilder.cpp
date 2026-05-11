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

#include "AttBuilder.hpp"
#include "ExodusMesh.hpp"

// verbose
#include "bstring.hpp"
#include "io.hpp"
#include "timer.hpp"

// attenuation
#include "AttenuationFull.hpp"
#include "AttenuationCG4.hpp"

// constructor
AttBuilder::AttBuilder(const ExodusMesh& exodusMesh, bool cg4, double dt) :
    mFmin(exodusMesh.getGlobalVariable("attenuation_f_min")),
    mFmax(exodusMesh.getGlobalVariable("attenuation_f_max")),
    mFref(exodusMesh.getGlobalVariable("attenuation_f_ref")),
    mNSLS(exodusMesh.getGlobalVariable("nr_lin_solids")), mUseCG4(cg4), mDt(dt) {
  timer::gPreloopTimer.begin("Creating attenuation builder");
  mW = mY = axisem3d::eigen::DColX(mNSLS);
  for (int isls = 0; isls < mNSLS; isls++) {
    mW[isls] = exodusMesh.getGlobalVariable("attenuation_w_" + bstring::toString(isls));
    mY[isls] = exodusMesh.getGlobalVariable("attenuation_y_" + bstring::toString(isls));
  }

  // check CG4
  if (mUseCG4 && spectral::nPol != 4) {
    throw std::runtime_error("AttBuilder::AttBuilder || "
                             "CG4 attenuation requires set(NPOL 4) "
                             "in CMakeLists.txt");
  }

  // frequency out of attenuation band, incorrect
  double freqMesh = 1. / exodusMesh.getGlobalVariable("minimum_period");
  if (freqMesh < mFmin || freqMesh > mFmax) {
    throw std::runtime_error("AttBuilder::AttBuilder || "
                             "Mesh frequency lies out of attenuation band. || "
                             "Change --attenuation.frequencies in mesher. || "
                             "Mesh frequency   = " +
        bstring::toString(freqMesh) +
        " || "
        "Attenuation band = " +
        bstring::range(mFmin, mFmax));
  }

  // frequency nearly out of attenuation band, inaccurate
  double margin = .1;
  double logMin = log10(mFmin);
  double logMax = log10(mFmax);
  double fminMargin = pow(10., logMin + (logMax - logMin) * margin);
  double fmaxMargin = pow(10., logMax - (logMax - logMin) * margin);
  if (freqMesh < fminMargin || freqMesh > fmaxMargin) {
    if (io::gVerboseWarnings) {
      io::cout << bstring::warning("AttBuilder::AttBuilder || "
                                   "Mesh frequency lies near the edge of attenuation band. || "
                                   "Consider changing --attenuation.frequencies in mesher. || "
                                   "Mesh frequency   = " +
          bstring::toString(freqMesh) +
          " || "
          "Attenuation band = " +
          bstring::range(mFmin, mFmax));
    }
  }
  timer::gPreloopTimer.ended("Creating attenuation builder");

  // set alpha, beta and gamma in Attenuation
  timer::gPreloopTimer.begin("Computing alpha, beta and gamma");
  setAlphaBetaGamma();
  timer::gPreloopTimer.ended("Computing alpha, beta and gamma");
}

// verbose
std::string
AttBuilder::verbose(const std::unique_ptr<const AttBuilder>& attBuilder) {
  using namespace bstring;
  std::stringstream ss;
  ss << bstring::boxTitle("Attenuation");
  if (attBuilder) {
    ss << boxEquals(0, 19, "frequency range", range(attBuilder->mFmin, attBuilder->mFmax));
    ss << boxEquals(0, 19, "reference frequency", attBuilder->mFref);
    ss << boxEquals(0, 19, "# std linear solids", attBuilder->mNSLS);
    ss << boxEquals(0, 19, "using CG4 mode", attBuilder->mUseCG4);
  } else {
    ss << "* Attenuation is turned off.\n";
  }
  ss << bstring::boxBaseline() << "\n\n";
  return ss.str();
}

// create attenuation
std::unique_ptr<Attenuation>
AttBuilder::createAttenuation(const axisem3d::eigen::DMatXN& QKp,
    const axisem3d::eigen::DMatXN& QMu,
    axisem3d::eigen::DMatXN& kp,
    axisem3d::eigen::DMatXN& mu,
    const axisem3d::eigen::DRow4& weightsCG4,
    bool elastic1D) const {
  // constants
  static double ysum = mY.sum();
  static double w0 = mFref * (2 * numerical::dPi);
  static double w1 = sqrt(mFmin * mFmax) * (2. * numerical::dPi);
  static double w10f = 2. * log(w1 / w0) / numerical::dPi;
  static double fact =
      (mY.array() * mW.array().square() / (mW.array().square() + w1 * w1)).sum() / ysum;
  const axisem3d::eigen::DMatXN& ones = axisem3d::eigen::DMatXN::Ones(QKp.rows(), QKp.cols());

  // kappa
  const axisem3d::eigen::DMatXN& kpNA = ones + w10f * QKp.cwiseInverse();
  axisem3d::eigen::DMatXN dkp = kpNA.cwiseQuotient(QKp / ysum + (1. - fact) * ones);
  const axisem3d::eigen::DMatXN& kpAT = kpNA + dkp * fact;
  dkp.array() *= kp.array();

  // mu
  const axisem3d::eigen::DMatXN& muNA = ones + w10f * QMu.cwiseInverse();
  axisem3d::eigen::DMatXN dmu = muNA.cwiseQuotient(QMu / ysum + (1. - fact) * ones);
  const axisem3d::eigen::DMatXN& muAT = muNA + dmu * fact;
  dmu.array() *= mu.array();

  // lambda
  const axisem3d::eigen::DMatXN& dlm = dkp - (2. / 3.) * dmu;

  // 1D or 3D
  if (mUseCG4) {
    // extract CG4 from full
    int nr = (int)dlm.rows();
    axisem3d::eigen::DMatX4 dlm4(nr, 4), dmu4(nr, 4);
    const std::vector<int> ipols = {1, 1, 3, 3};
    const std::vector<int> jpols = {1, 3, 1, 3};
    for (int ip4 = 0; ip4 < 4; ip4++) {
      int ipnt = ipols[ip4] * spectral::nPED + jpols[ip4];
      dlm4.col(ip4) = weightsCG4(ip4) * dlm.col(ipnt);
      dmu4.col(ip4) = weightsCG4(ip4) * dmu.col(ipnt);
    }

    // correct original
    // factor = NA + weightCG4 * (AT - NA), on CG4 points
    // factor = NA                        , on other points
    kp.array() *= kpNA.array();
    mu.array() *= muNA.array();
    const auto& ones = axisem3d::eigen::DColX::Ones(nr).array();
    for (int ip4 = 0; ip4 < 4; ip4++) {
      int ipnt = ipols[ip4] * spectral::nPED + jpols[ip4];
      kp.col(ipnt).array() *=
          ones + weightsCG4(ip4) * (kpAT.col(ipnt).array() / kpNA.col(ipnt).array() - ones);
      mu.col(ipnt).array() *=
          ones + weightsCG4(ip4) * (muAT.col(ipnt).array() / muNA.col(ipnt).array() - ones);
    }

    // return
    if (elastic1D) {
      return std::make_unique<AttenuationCG4>(
          Eigen::Map<const axisem3d::eigen::DMat22_RM>(dlm4.data()).eval(),
          Eigen::Map<const axisem3d::eigen::DMat22_RM>(dmu4.data()).eval());
    } else {
      return std::make_unique<AttenuationCG4>(dlm4, dmu4);
    }
  } else {
    // correct original
    kp.array() *= kpAT.array();
    mu.array() *= muAT.array();
    if (elastic1D) {
      return std::make_unique<AttenuationFull>(op1D_3D::toPP(dlm), op1D_3D::toPP(dmu));
    } else {
      return std::make_unique<AttenuationFull>(dlm, dmu);
    }
  }
}

// set alpha, beta, gamma
// NOTE: with SPECFEM legency deprecated, alpha, beta and gamma
//       all become element-independent
void
AttBuilder::setAlphaBetaGamma() const {
  // constants
  const axisem3d::eigen::DColX& yUn = mY / mY.sum();
  const axisem3d::eigen::DColX& wDt = mW * mDt;
  const axisem3d::eigen::DColX& nSLS1 = axisem3d::eigen::DColX::Ones(mNSLS);
  // alpha
  const axisem3d::eigen::DColX& alpha = (-wDt).array().exp();
  // beta
  const axisem3d::eigen::DColX& betta =
      yUn.cwiseProduct((nSLS1 - alpha).cwiseQuotient(wDt) - alpha);
  // gamma
  const axisem3d::eigen::DColX& gamma =
      yUn.cwiseProduct((alpha - nSLS1).cwiseQuotient(wDt) + nSLS1);
  // set statics in Attenuation
  Attenuation::setAlphaBetaGamma(alpha, betta, gamma);
}
