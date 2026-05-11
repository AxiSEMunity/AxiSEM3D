// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: MIT
// Copyright (C) 2019 - 2026 by the AxiSEM3D authors
//
// This file is part of the AxiSEM3D library. See the LICENSE file for details.
//
// -----------------------------------------------------------------------------

//  vertical undulation
//  generator of PRT in core

#include "Undulation.hpp"
#include "SolverFFTW.cpp"
#include "Quad.hpp"
#include "geodesy.hpp"

using spectral::nPEM;

// finishing 3D properties
void
Undulation::finishing3D() const {
  // no undulation
  if (!mDeltaZ) {
    return;
  }

  // initialize static
  int nr = (int)getElemental().rows();
  sFFT_N1.addNR(nr);
  sFFT_N3.addNR(nr);
  int nc = nr / 2 + 1;
  if (sDeltaZ_Fourier.size() < nc) {
    sDeltaZ_Fourier.resize(nc);
    sDeltaZ_SPZ_Fourier.resize(nc);
  }
}

// finished 3D properties
void
Undulation::finished3D(const Quad& myQuad) {
  // no undulation
  if (!mDeltaZ) {
    return;
  }

  ////////////////////// compute gradient //////////////////////
  // step 1: FFT
  const axisem3d::eigen::DMatXN& dz = getElemental();
  int nr = (int)dz.rows();
  sFFT_N1.computeR2C(dz, sDeltaZ_Fourier, nr);

  // step 2: gradient
  static axisem3d::eigen::DMat2N sz;
  static axisem3d::eigen::DMatPP_RM ifPP;
  const auto& grad = myQuad.createGradient<double>(sz, ifPP);
  grad->checkCompatibility(nr);
  grad->computeGrad3(sDeltaZ_Fourier, sDeltaZ_SPZ_Fourier, nr / 2 + 1);

  // step 3: IFFT
  mDeltaZ_RTZ.resize(nr, nPEM * 3);
  sFFT_N3.computeC2R(sDeltaZ_SPZ_Fourier, mDeltaZ_RTZ, nr);

  // step 4: rotate
  if (!geodesy::isCartesian()) {
    // theta
    const axisem3d::eigen::DMat2N& rt = geodesy::sz2rtheta(sz, false);
    const axisem3d::eigen::DRowN& cost = rt.array().row(1).cos();
    const axisem3d::eigen::DRowN& sint = rt.array().row(1).sin();
    // s, z
    const axisem3d::eigen::DMatXN& dZds = mDeltaZ_RTZ.block(0, nPEM * 0, nr, nPEM);
    const axisem3d::eigen::DMatXN& dZdz = mDeltaZ_RTZ.block(0, nPEM * 2, nr, nPEM);
    // s, z -> R, Z
    mDeltaZ_RTZ.block(0, nPEM * 0, nr, nPEM) =
        dZds.array() * cost.array().replicate(nr, 1) - dZdz.array() * sint.array().replicate(nr, 1);
    mDeltaZ_RTZ.block(0, nPEM * 2, nr, nPEM) =
        dZdz.array() * cost.array().replicate(nr, 1) + dZds.array() * sint.array().replicate(nr, 1);
  }
}

// get Jacobian for mass
axisem3d::eigen::arN_DColX
Undulation::getMassJacobian(const axisem3d::eigen::DMat2N& sz) const {
  // no undulation
  if (!mDeltaZ) {
    return axisem3d::eigen::arN_DColX();
  }

  // allocate dZdZ with the same size as dZ
  const axisem3d::eigen::arN_DColX& dZ = getPointwise();
  axisem3d::eigen::arN_DColX dZdZ = dZ;
  int maxNr = (int)mDeltaZ_RTZ.rows();
  for (int ipnt = 0; ipnt < nPEM; ipnt++) {
    int nr = (int)dZdZ[ipnt].rows();
    int col = nPEM * 2 + ipnt;
    // use linear interpolation to shrink dimension
    // because using fft may cause over-shooting
    PhysicalProperty<nPEM>::linInterpPhi(mDeltaZ_RTZ, dZdZ[ipnt], col, 0, maxNr, nr);
  }

  // compute Jacobian
  // based on eq.(6) in Leng et al., 2019
  axisem3d::eigen::arN_DColX J;
  for (int ipnt = 0; ipnt < nPEM; ipnt++) {
    const axisem3d::eigen::DColX& J3 = 1. + dZdZ[ipnt].array();
    axisem3d::eigen::DColX J0;
    if (geodesy::isCartesian()) {
      J0 = axisem3d::eigen::DColX::Ones(J3.rows());
    } else {
      double Z = sz.col(ipnt).norm();
      J0 = (Z > numerical::dEpsilon) ? (1. + dZ[ipnt].array() / Z) : J3;
    }
    J[ipnt] = J0.array().square() * J3.array();
    if (J[ipnt].minCoeff() < numerical::dEpsilon) {
      throw std::runtime_error("Undulation::getMassJacobian || "
                               "Negative Jacobian or distorted element occurred.");
    }
  }
  return J;
}

// create PRT
std::unique_ptr<const PRT>
Undulation::createPRT(const axisem3d::eigen::DMat2N& sz) const {
  // no undulation
  if (!mDeltaZ) {
    return nullptr;
  }

  // compute Jacobian
  // based on eq.(6) in Leng et al., 2019
  int nr = (int)mDeltaZ_RTZ.rows();
  const axisem3d::eigen::DMatXN& J1 = mDeltaZ_RTZ.block(0, nPEM * 0, nr, nPEM);
  const axisem3d::eigen::DMatXN& J2 = mDeltaZ_RTZ.block(0, nPEM * 1, nr, nPEM);
  const axisem3d::eigen::DMatXN& J3 = 1. + mDeltaZ_RTZ.block(0, nPEM * 2, nr, nPEM).array();
  axisem3d::eigen::DMatXN J0;
  if (geodesy::isCartesian()) {
    J0 = axisem3d::eigen::DMatXN::Ones(nr, nPEM);
  } else {
    const axisem3d::eigen::DMatXN& dZ = getElemental();
    const axisem3d::eigen::DMatXN& Z = sz.colwise().norm().replicate(nr, 1);
    J0 = (Z.array() > numerical::dEpsilon).select(1. + dZ.array() / Z.array(), J3);
  }

  // compute inverse Jacobian
  // based on eq.(13) in Leng et al., 2019
  const axisem3d::eigen::DMatXN& X0 = J0.cwiseInverse();
  const axisem3d::eigen::DMatXN& X1 = -J1.cwiseQuotient((J0.cwiseProduct(J3)));
  const axisem3d::eigen::DMatXN& X2 = -J2.cwiseQuotient((J0.cwiseProduct(J3)));
  const axisem3d::eigen::DMatXN& X3 = J3.cwiseInverse();
  const axisem3d::eigen::DMatXN& XJ = J0.array().square() * J3.array();
  if (XJ.minCoeff() < numerical::dEpsilon) {
    throw std::runtime_error("Undulation::createPRT || "
                             "Negative Jacobian or distorted element occurred.");
  }

  // create PRT
  if (nr == 1) {
    return std::make_unique<PRT>(op1D_3D::toPP(X0),
        op1D_3D::toPP(X1),
        op1D_3D::toPP(X2),
        op1D_3D::toPP(X3),
        op1D_3D::toPP(XJ));
  } else {
    return std::make_unique<PRT>(X0, X1, X2, X3, XJ);
  }
}

// compute 3D normal at a point
axisem3d::eigen::DMatX3
Undulation::computeNormal3D(
    const axisem3d::eigen::DCol2& n1D, const axisem3d::eigen::DMat2N& sz, int ipnt) const {
  // no undulation
  if (!mDeltaZ) {
    return (axisem3d::eigen::DMatX3(1, 3) << n1D(0), 0., n1D(1)).finished();
  }

  // dZ
  const axisem3d::eigen::arN_DColX& dZ = getPointwise();

  // extract gradient of dZ
  int nr = (int)dZ[ipnt].rows();
  int maxNr = (int)mDeltaZ_RTZ.rows();
  axisem3d::eigen::DColX dZdR(nr), dZdT(nr), dZdZ(nr);
  PhysicalProperty<nPEM>::linInterpPhi(mDeltaZ_RTZ, dZdR, nPEM * 0 + ipnt, 0, maxNr, nr);
  PhysicalProperty<nPEM>::linInterpPhi(mDeltaZ_RTZ, dZdT, nPEM * 1 + ipnt, 0, maxNr, nr);
  PhysicalProperty<nPEM>::linInterpPhi(mDeltaZ_RTZ, dZdZ, nPEM * 2 + ipnt, 0, maxNr, nr);

  // compute Jacobian
  // based on eq.(6) in Leng et al., 2019
  const axisem3d::eigen::DColX& J1 = dZdR;
  const axisem3d::eigen::DColX& J2 = dZdT;
  const axisem3d::eigen::DColX& J3 = 1. + dZdZ.array();
  axisem3d::eigen::DColX J0;
  if (geodesy::isCartesian()) {
    J0 = axisem3d::eigen::DColX::Ones(J3.rows());
  } else {
    double Z = sz.col(ipnt).norm();
    J0 = (Z > numerical::dEpsilon) ? (1. + dZ[ipnt].array() / Z) : J3;
  }

  // K = J^-T |J|
  const axisem3d::eigen::DColX& K0 = J3.cwiseProduct(J0);
  const axisem3d::eigen::DColX& K1 = -J1.cwiseProduct(J0);
  const axisem3d::eigen::DColX& K2 = -J2.cwiseProduct(J0);
  const axisem3d::eigen::DColX& K3 = J0.cwiseProduct(J0);

  // rotate n1D to RTZ
  double sint = 0., cost = 0.;
  axisem3d::eigen::DCol2 n1D_RZ = n1D;
  if (!geodesy::isCartesian()) {
    const axisem3d::eigen::DCol2& rt = geodesy::sz2rtheta(sz.col(ipnt).eval(), false);
    sint = sin(rt(1));
    cost = cos(rt(1));
    n1D_RZ(0) = n1D(0) * cost - n1D(1) * sint;
    n1D_RZ(1) = n1D(1) * cost + n1D(0) * sint;
  }

  // K * n
  // based on eq.(A4) in Leng et al., 2019
  axisem3d::eigen::DMatX3 n3D_RTZ(nr, 3);
  n3D_RTZ.col(0) = K1 * n1D_RZ(1) + K0 * n1D_RZ(0);
  n3D_RTZ.col(1) = K2 * n1D_RZ(1);
  n3D_RTZ.col(2) = K3 * n1D_RZ(1);

  // rotate back to SPZ
  if (!geodesy::isCartesian()) {
    axisem3d::eigen::DMatX3 n3D_SPZ(nr, 3);
    n3D_SPZ.col(0) = n3D_RTZ.col(0) * cost + n3D_RTZ.col(2) * sint;
    n3D_SPZ.col(2) = n3D_RTZ.col(2) * cost - n3D_RTZ.col(0) * sint;
    n3D_SPZ.col(1) = n3D_RTZ.col(1);
    return n3D_SPZ;
  } else {
    return n3D_RTZ;
  }
}

///////////////////////// static /////////////////////////
// finished 3D properties
void
Undulation::finished3D() {
  // fft
  sFFT_N1.createPlans(.1);
  sFFT_N3.createPlans(.1);
  // gradient
  GradientQuadrature<double>::setGMat(spectrals::gGMatrixGLL, spectrals::gGMatrixGLJ);
}
