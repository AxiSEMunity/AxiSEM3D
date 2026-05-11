// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: MIT
// Copyright (C) 2019 - 2026 by the AxiSEM3D authors
//
// This file is part of the AxiSEM3D library. See the LICENSE file for details.
//
// -----------------------------------------------------------------------------

//  GLL point for preloop processing
//  generator of Point and boundary conditions in core

#ifndef GLLPoint_hpp
#define GLLPoint_hpp

#include "eigen_sem.hpp"
#include <map>
#include <vector>
#include <memory>

// release
class ABC;
class TimeScheme;
class Domain;
class SolidPoint;
class FluidPoint;

class GLLPoint {
  public:
  // setup by Quad
  void
  setup(int globalTag, int nr, const axisem3d::eigen::DCol2& coords, double distTol) {
    if (mGlobalTag == -1) {
      // set by the first element
      mGlobalTag = globalTag;
      mNr = nr;
      mCoords = coords;
    } else {
      // check subsequent elements
      if (mGlobalTag != globalTag || mNr != nr || (mCoords - coords).norm() > distTol) {
        throw std::runtime_error("GLLPoint::setup || "
                                 "Conflict in GLL-point setup.");
      }
    }
    // reference count
    mElementCount++;
  }

  // add mass
  void
  addMass(const axisem3d::eigen::DColX& mass, bool fluid) {
    axisem3d::eigen::DColX& myMass = fluid ? mMassFluid : mMassSolid;
    op1D_3D::addTo(mass, myMass);
  }

  // add solid-fluid
  void
  addNormalSF(const axisem3d::eigen::DMatX3& nSF) {
    op1D_3D::addTo(nSF, mNormalSFU);
    mNormalSFA = mNormalSFU;
  }

  // add Clayton ABC
  // NOTE: the normal is NOT assembled because velocity and density
  //       can be discontinuous; at a boundary point shared by M
  //       elements, M Clayton instances will be created
  void
  addClaytonABC(const std::string& key,
      bool fluid,
      const axisem3d::eigen::DMatX3& nABC,
      const axisem3d::eigen::DColX& rho,
      const axisem3d::eigen::DColX& vp,
      const axisem3d::eigen::DColX& vs) {
    axisem3d::eigen::DColX rhoVp, rhoVs;
    op1D_3D::times(rho, vp, rhoVp);
    op1D_3D::times(rho, vs, rhoVs);
    // init (key, vector<tuple>)
    mClaytonABC.insert({key,
        std::vector<std::tuple<bool,
            axisem3d::eigen::DMatX3,
            axisem3d::eigen::DColX,
            axisem3d::eigen::DColX>>()});
    mClaytonABC.at(key).push_back({fluid, nABC, rhoVp, rhoVs});
    // try reduce to 1D
    op1D_3D::tryReduceTo1D(std::get<1>(mClaytonABC.at(key).back()));
    op1D_3D::tryReduceTo1D(std::get<2>(mClaytonABC.at(key).back()));
    op1D_3D::tryReduceTo1D(std::get<3>(mClaytonABC.at(key).back()));
  }

  // add gamma
  void
  addGamma(const axisem3d::eigen::DColX& Gamma) {
    op1D_3D::addTo(Gamma, mGamma);
    mCountGammasAdded++;
  }

  // add ocean load
  void
  addOceanLoad(const axisem3d::eigen::DMatX3& nTop, const axisem3d::eigen::DColX& sumRhoDepth) {
    op1D_3D::addTo(nTop, mNormalTop);
    // this is done repeatedly at shared points
    mSumRhoDepth = sumRhoDepth;
  }

  // set axial
  void
  setAxial() {
    mAxial = true;
  }

  // set surface
  void
  setSurface() {
    mSurface = true;
  }

  // comm size
  int
  getCommSize() const {
    // must use this full size because the compact size
    // of recv buffer is unknown
    return 1 + 4 + mNr + mNr + mNr * 3 + mNr * 3 + (2 + mNr);
  }

  // feed comm
  void
  feedComm(axisem3d::eigen::DColX& buffer, int& row) const {
    // reference element count
    buffer(row, 0) = mElementCount;
    // size info
    // sizes must be sent as they can have different sizes on ranks
    buffer(row + 1, 0) = mMassFluid.size();
    buffer(row + 2, 0) = mMassSolid.size();
    buffer(row + 3, 0) = mNormalSFA.size();
    buffer(row + 4, 0) = mNormalTop.size();
    buffer(row + 5, 0) = mGamma.size();
    buffer(row + 6, 0) = mCountGammasAdded;
    row += 7;
    // mass
    buffer.block(row, 0, mMassFluid.size(), 1) = mMassFluid;
    row += mMassFluid.size();
    buffer.block(row, 0, mMassSolid.size(), 1) = mMassSolid;
    row += mMassSolid.size();
    // solid-fluid
    buffer.block(row, 0, mNormalSFA.size(), 1) =
        Eigen::Map<const axisem3d::eigen::DColX>(mNormalSFA.data(), mNormalSFA.size(), 1);
    row += mNormalSFA.size();
    // ocean load
    buffer.block(row, 0, mNormalTop.size(), 1) =
        Eigen::Map<const axisem3d::eigen::DColX>(mNormalTop.data(), mNormalTop.size(), 1);
    row += mNormalTop.size();
    // sponge boundary
    buffer.block(row, 0, mGamma.size(), 1) = mGamma;
    row += mGamma.size();
  }

  // extract comm
  void
  extractComm(const axisem3d::eigen::DColX& buffer, int& row) {
    // reference element count
    mElementCount += (int)round(buffer(row, 0));
    // size info
    int sizeMassFluid = (int)round(buffer(row + 1, 0));
    int sizeMassSolid = (int)round(buffer(row + 2, 0));
    int sizeNormalSFA = (int)round(buffer(row + 3, 0));
    int sizeNormalTop = (int)round(buffer(row + 4, 0));
    int sizeGamma = (int)round(buffer(row + 5, 0));
    int countGammasAdded = (int)round(buffer(row + 6, 0));
    row += 7;
    // mass
    op1D_3D::addTo(buffer.block(row, 0, sizeMassFluid, 1), mMassFluid);
    row += sizeMassFluid;
    op1D_3D::addTo(buffer.block(row, 0, sizeMassSolid, 1), mMassSolid);
    row += sizeMassSolid;
    // solid-fluid
    op1D_3D::addTo(Eigen::Map<const axisem3d::eigen::DMatX3>(
                       buffer.block(row, 0, sizeNormalSFA, 1).data(), sizeNormalSFA / 3, 3),
        mNormalSFA);
    row += sizeNormalSFA;
    // ocean load
    op1D_3D::addTo(Eigen::Map<const axisem3d::eigen::DMatX3>(
                       buffer.block(row, 0, sizeNormalTop, 1).data(), sizeNormalTop / 3, 3),
        mNormalTop);
    row += sizeNormalTop;
    // sponge boundary
    op1D_3D::addTo(buffer.block(row, 0, sizeGamma, 1), mGamma);
    mCountGammasAdded += countGammasAdded;
    row += sizeGamma;
  }

  // release to domain
  void
  release(const ABC& abc, const TimeScheme& timeScheme, Domain& domain);

  // get SolidPoint after release
  const std::shared_ptr<SolidPoint>&
  getSolidPoint() const {
    return mSolidPoint;
  }

  // get FluidPoint after release
  const std::shared_ptr<FluidPoint>&
  getFluidPoint() const {
    return mFluidPoint;
  }

  // get global tag
  int
  getGlobalTag() const {
    return mGlobalTag;
  }

  // get reference element count
  int
  getElementCount() const {
    return mElementCount;
  }

  ///////////////////////////// data /////////////////////////////
  private:
  // tag
  int mGlobalTag = -1;

  // nr
  int mNr = -1;

  // coords
  axisem3d::eigen::DCol2 mCoords = axisem3d::eigen::DCol2::Zero();

  // reference element count
  int mElementCount = 0;

  // mass
  axisem3d::eigen::DColX mMassFluid = axisem3d::eigen::DColX::Zero(0);
  axisem3d::eigen::DColX mMassSolid = axisem3d::eigen::DColX::Zero(0);

  // solid-fuild
  // unassembled
  axisem3d::eigen::DMatX3 mNormalSFU = axisem3d::eigen::DMatX3::Zero(0, 3);
  // assembled
  axisem3d::eigen::DMatX3 mNormalSFA = axisem3d::eigen::DMatX3::Zero(0, 3);

  // Clayton ABC {key, {fluid/solid, normal, rho * vp, rho * vs}}
  std::map<std::string,
      std::vector<std::
              tuple<bool, axisem3d::eigen::DMatX3, axisem3d::eigen::DColX, axisem3d::eigen::DColX>>>
      mClaytonABC;

  // Kosloff_Kosloff sponge boundary
  axisem3d::eigen::DColX mGamma = axisem3d::eigen::DColX::Zero(0);
  int mCountGammasAdded = 0;

  // ocean load
  axisem3d::eigen::DMatX3 mNormalTop = axisem3d::eigen::DMatX3::Zero(0, 3);
  axisem3d::eigen::DColX mSumRhoDepth = axisem3d::eigen::DColX::Zero(0);

  // axial
  bool mAxial = false;

  // surface
  bool mSurface = false;

  // point pointers after release
  std::shared_ptr<SolidPoint> mSolidPoint = nullptr;
  std::shared_ptr<FluidPoint> mFluidPoint = nullptr;
};

#endif /* GLLPoint_hpp */
