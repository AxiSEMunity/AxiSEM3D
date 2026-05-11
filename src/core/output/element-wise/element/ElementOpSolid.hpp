// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: MIT
// Copyright (C) 2019 - 2026 by the AxiSEM3D authors
//
// This file is part of the AxiSEM3D library. See the LICENSE file for details.
//
// -----------------------------------------------------------------------------

//  element output in solid

#ifndef ElementOpSolid_hpp
#define ElementOpSolid_hpp

#include "ElementOp.hpp"
#include "SolidElement.hpp"
#include "Point.hpp"

class ElementOpSolid : public ElementOp {
  public:
  // constructor
  ElementOpSolid(const std::vector<int>& ipnts) : ElementOp(ipnts) {
    // nothing
  }

  /////////////////////////// setup ///////////////////////////
  // set element
  void
  setElement(const std::shared_ptr<SolidElement>& element) {
    mElement = element;
  }

  // set in group
  void
  setInGroup(int dumpIntv, const channel::solid::ChannelOptions& chops, int nphis);

  /////////////////////////// get sizes ///////////////////////////
  // get na
  int
  getNa(int nphis) const {
    if (nphis == 0) {
      return mElement->getNr();
    } else {
      return nphis;
    }
  }

  // get nu_1
  int
  getNu_1() const {
    return mElement->getNu_1();
  }

  // get element tag
  int
  getElementTag() const {
    return mElement->getQuadTag();
  }

  // get coords
  axisem3d::eigen::DRowX
  getCoords() const {
    axisem3d::eigen::DRowX sz(mIPnts.size() * 2);
    for (int ip = 0; ip < mIPnts.size(); ip++) {
      sz.block(0, ip * 2, 1, 2) = mElement->getPoint(mIPnts[ip]).getCoords();
    }
    return sz;
  }

  /////////////////////////// record ///////////////////////////
  public:
  // record
  void
  record(int bufferLine,
      const channel::solid::ChannelOptions& chops,
      const axisem3d::eigen::CMatXX& expIAlphaPhi);

  private:
  // element
  std::shared_ptr<SolidElement> mElement = nullptr;

  // buffer
  axisem3d::eigen::RTensor4 mBufferU;
  axisem3d::eigen::RTensor4 mBufferG;
  axisem3d::eigen::RTensor4 mBufferE;
  axisem3d::eigen::RTensor4 mBufferR;
  axisem3d::eigen::RTensor4 mBufferS;

  /////////////////////////// process ///////////////////////////
  public:
  // process and report to group
  void
  processReport(int bufferLine,
      const channel::solid::ChannelOptions& chops,
      int elemIndexNaGrid,
      int naGridIndex,
      std::vector<axisem3d::eigen::RTensor5>& ioBuffers);

  ////////////////////////////////////////
  //////////////// static ////////////////
  ////////////////////////////////////////

  private:
  // expand workspace for record
  static void
  expandWorkspaceRecord(int nu_1, int na, const channel::solid::ChannelOptions& chops) {
    // use (na + 1) to handle Nyquist
    int na_1 = na + 1;

    if (chops.mNeedBufferU) {
      if (sCUXN3.rows() < nu_1) {
        sCUXN3.resize(nu_1, spectral::nPEM * 3);
      }
      if (sRUXN3.rows() < na_1) {
        sRUXN3.resize(na_1, spectral::nPEM * 3);
      }
    }
    if (chops.mNeedBufferG) {
      if (sCGXN9.rows() < nu_1) {
        sCGXN9.resize(nu_1, spectral::nPEM * 9);
      }
      if (sRGXN9.rows() < na_1) {
        sRGXN9.resize(na_1, spectral::nPEM * 9);
      }
    }
    if (chops.mNeedBufferE) {
      if (sCEXN6.rows() < nu_1) {
        sCEXN6.resize(nu_1, spectral::nPEM * 6);
      }
      if (sREXN6.rows() < na_1) {
        sREXN6.resize(na_1, spectral::nPEM * 6);
      }
    }
    if (chops.mNeedBufferR) {
      if (sCRXN3.rows() < nu_1) {
        sCRXN3.resize(nu_1, spectral::nPEM * 3);
      }
      if (sRRXN3.rows() < na_1) {
        sRRXN3.resize(na_1, spectral::nPEM * 3);
      }
    }
    if (chops.mNeedBufferS) {
      if (sCSXN6.rows() < nu_1) {
        sCSXN6.resize(nu_1, spectral::nPEM * 6);
      }
      if (sRSXN6.rows() < na_1) {
        sRSXN6.resize(na_1, spectral::nPEM * 6);
      }
    }
  }

  // workspace for record
  // get response from elememt
  inline static axisem3d::eigen::CMatXN3 sCUXN3 = axisem3d::eigen::CMatXN3(0, spectral::nPEM * 3);
  inline static axisem3d::eigen::CMatXN9 sCGXN9 = axisem3d::eigen::CMatXN9(0, spectral::nPEM * 9);
  inline static axisem3d::eigen::CMatXN6 sCEXN6 = axisem3d::eigen::CMatXN6(0, spectral::nPEM * 6);
  inline static axisem3d::eigen::CMatXN3 sCRXN3 = axisem3d::eigen::CMatXN3(0, spectral::nPEM * 3);
  inline static axisem3d::eigen::CMatXN6 sCSXN6 = axisem3d::eigen::CMatXN6(0, spectral::nPEM * 6);

  // making real
  inline static axisem3d::eigen::RMatXN3_RM sRUXN3 =
      axisem3d::eigen::RMatXN3_RM(0, spectral::nPEM * 3);
  inline static axisem3d::eigen::RMatXN9_RM sRGXN9 =
      axisem3d::eigen::RMatXN9_RM(0, spectral::nPEM * 9);
  inline static axisem3d::eigen::RMatXN6_RM sREXN6 =
      axisem3d::eigen::RMatXN6_RM(0, spectral::nPEM * 6);
  inline static axisem3d::eigen::RMatXN3_RM sRRXN3 =
      axisem3d::eigen::RMatXN3_RM(0, spectral::nPEM * 3);
  inline static axisem3d::eigen::RMatXN6_RM sRSXN6 =
      axisem3d::eigen::RMatXN6_RM(0, spectral::nPEM * 6);
};

#endif /* ElementOpSolid_hpp */
