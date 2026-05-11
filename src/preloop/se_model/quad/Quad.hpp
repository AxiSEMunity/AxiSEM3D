// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: MIT
// Copyright (C) 2019 - 2026 by the AxiSEM3D authors
//
// This file is part of the AxiSEM3D library. See the LICENSE file for details.
//
// -----------------------------------------------------------------------------

//  Quadrilateral
//  generator of Element in core

#ifndef Quad_hpp
#define Quad_hpp

// eigen
#include "eigen_sem.hpp"

// components
#include "Mapping.hpp"
#include "Material.hpp"
#include "Undulation.hpp"
#include "OceanLoad.hpp"

// external
class ExodusMesh;
class ABC;
class LocalMesh;
class GLLPoint;

// release
class Domain;
class AttBuilder;
class SolidElement;
class FluidElement;
#include "GradientQuadrature.hpp"

// measure
class Element;

class Quad {
  public:
  // constructor
  Quad(
      const ExodusMesh& exodusMesh, const LocalMesh& localMesh, int localTag, bool useLuckyNumbers);

  // finishing 3D properties
  void
  finishing3D() {
    mUndulation->finishing3D();
  }

  // finished 3D properties
  void
  finished3D() {
    mMaterial->finished3D();
    mUndulation->finished3D(*this);
  }

  // setup GLL
  void
  setupGLL(const ABC& abc, const LocalMesh& localMesh, std::vector<GLLPoint>& GLLPoints) const;

  // compute dt
  double
  computeDt(double courant, const ABC& abc) const;

  // get nodal sz
  const axisem3d::eigen::DMat24&
  getNodalSZ() const {
    return mMapping->getNodalSZ();
  }

  // release to domain
  void
  release(const LocalMesh& localMesh,
      const std::vector<GLLPoint>& GLLPoints,
      const std::unique_ptr<const AttBuilder>& attBuilder,
      Domain& domain);

  // inverse mapping: (s,z) -> (ξ,η)
  // return true if (s,z) is inside this element
  bool
  inverseMapping(const axisem3d::eigen::DCol2& sz,
      axisem3d::eigen::DCol2& xieta,
      double maxIter = 10,
      double tolerance = 1e-9) const {
    return mMapping->inverseMapping(sz, xieta, maxIter, tolerance);
  }

  private:
  //////////////////////// interal ////////////////////////
  // compute integral factor
  axisem3d::eigen::DRowN
  computeIntegralFactor(
      axisem3d::eigen::DMat2N& sz, std::array<axisem3d::eigen::DMat22, spectral::nPEM>& J) const;

  // get normal
  void
  computeNormal(int edge,
      const axisem3d::eigen::DMat2N& sz,
      const std::array<axisem3d::eigen::DMat22, spectral::nPEM>& J,
      std::vector<int>& ipnts,
      axisem3d::eigen::arP_DMatX3& normal) const;

  // weights for CG4 attenuation
  axisem3d::eigen::DRow4
  computeWeightsCG4(const axisem3d::eigen::DMatPP_RM& ifPP) const;

  public:
  //////////////////////// get ////////////////////////
  // get global tag
  inline int
  getGlobalTag() const {
    return mGlobalTag;
  }

  // get point nr
  inline const axisem3d::eigen::IRowN&
  getPointNr() const {
    return *mPointNr;
  }

  // get point sz
  axisem3d::eigen::DMat2N
  getPointSZ() const {
    axisem3d::eigen::DMat2N pointSZ;
    const axisem3d::eigen::DMat2N& xieta = spectrals::getXiEtaElement(axial());
    for (int ipnt = 0; ipnt < spectral::nPEM; ipnt++) {
      pointSZ.col(ipnt) = mMapping->mapping(xieta.col(ipnt));
    }
    return pointSZ;
  }

  // get undulation
  axisem3d::eigen::arN_DColX
  getUndulation() const {
    return mUndulation->getPointwise();
  }

  // get material
  std::unique_ptr<Material>&
  getMaterialPtr() {
    return mMaterial;
  }

  // get undulation
  std::unique_ptr<Undulation>&
  getUndulationPtr() {
    return mUndulation;
  }

  // get oceanload
  std::unique_ptr<OceanLoad>&
  getOceanLoadPtr() {
    return mOceanLoad;
  }

  // surface edge
  int
  getSurfaceEdge() const {
    return mEdgesOnBoundary.at("TOP");
  }

  // fluid
  bool
  fluid() const {
    return mFluid;
  }

  // axial
  bool
  axial() const {
    return mEdgesOnBoundary.at("LEFT") != -1;
  }

  // create gradient
  template <typename floatT>
  std::unique_ptr<const GradientQuadrature<floatT>>
  createGradient(axisem3d::eigen::DMat2N& sz, axisem3d::eigen::DMatPP_RM& ifPP) const {
    // integral factor
    static std::array<axisem3d::eigen::DMat22, spectral::nPEM> J;
    const axisem3d::eigen::DRowN& ifact = computeIntegralFactor(sz, J);
    // save integral factor for later use (CG4)
    ifPP = Eigen::Map<const axisem3d::eigen::DMatPP_RM>(ifact.data());

    // compute Jacobian and s^-1
    static axisem3d::eigen::DMatPP_RM dsdxii, dsdeta, dzdxii, dzdeta, inv_s;
    for (int ipol = 0; ipol < spectral::nPED; ipol++) {
      for (int jpol = 0; jpol < spectral::nPED; jpol++) {
        int ipnt = ipol * spectral::nPED + jpol;
        // Jacobian
        double detJ = J[ipnt].determinant();
        dsdxii(ipol, jpol) = J[ipnt](0, 0) / detJ;
        dsdeta(ipol, jpol) = -J[ipnt](0, 1) / detJ;
        dzdxii(ipol, jpol) = -J[ipnt](1, 0) / detJ;
        dzdeta(ipol, jpol) = J[ipnt](1, 1) / detJ;
        // s^-1, use zero when s=0
        double s = sz(0, ipnt);
        inv_s(ipol, jpol) = (axial() && ipol == 0) ? 0. : 1. / s;
      }
    }

    // return
    return std::make_unique<GradientQuadrature<floatT>>(
        dsdxii, dsdeta, dzdxii, dzdeta, inv_s, axial(), ifPP);
  }

  // get element
  std::shared_ptr<const Element>
  getElement() const;

  // get solid element
  const std::shared_ptr<SolidElement>&
  getSolidElement() const;

  // get fluid element
  const std::shared_ptr<FluidElement>&
  getFluidElement() const;

  //////////////////////////////////////////////////////
  //////////////////////// data ////////////////////////
  //////////////////////////////////////////////////////
  private:
  // tags
  const int mLocalTag;
  const int mGlobalTag;

  // solid-fluid
  const bool mFluid;

  // model boundary
  std::map<std::string, int> mEdgesOnBoundary;

  // nr field
  std::unique_ptr<axisem3d::eigen::IRowN> mPointNr = nullptr;

  ////////////// components //////////////
  // mapping
  std::unique_ptr<Mapping> mMapping = nullptr;
  // material
  std::unique_ptr<Material> mMaterial = nullptr;
  // undulation
  std::unique_ptr<Undulation> mUndulation = nullptr;
  // ocean load
  std::unique_ptr<OceanLoad> mOceanLoad = nullptr;

  // element pointers after release
  std::shared_ptr<SolidElement> mSolidElement = nullptr;
  std::shared_ptr<FluidElement> mFluidElement = nullptr;
};

#endif /* Quad_hpp */
