// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: MIT
// Copyright (C) 2019 - 2026 by the AxiSEM3D authors
//
// This file is part of the AxiSEM3D library. See the LICENSE file for details.
//
// -----------------------------------------------------------------------------

//  station in solid

#ifndef StationSolid_hpp
#define StationSolid_hpp

#include "Station.hpp"
class SolidElement;

class StationSolid : public Station {
  public:
  // constructor
  StationSolid(const std::string& key, double phi, double theta, double backAzimuth) :
      Station(key, phi, theta, backAzimuth) {
    // nothing
  }

  /////////////////////////// setup ///////////////////////////
  // set element
  void
  setElement(const std::shared_ptr<SolidElement>& element, const axisem3d::eigen::DRowN& weights);

  // set in group
  void
  setInGroup(int dumpIntv, const channel::solid::ChannelOptions& chops);

  /////////////////////////// record ///////////////////////////
  public:
  // record
  void
  record(int bufferLine, const channel::solid::ChannelOptions& chops);

  private:
  // element
  std::shared_ptr<SolidElement> mElement = nullptr;

  // buffer
  axisem3d::eigen::RMatX3_RM mBufferU = axisem3d::eigen::RMatX3_RM(0, 3);
  axisem3d::eigen::RMatX9_RM mBufferG = axisem3d::eigen::RMatX9_RM(0, 9);
  axisem3d::eigen::RMatX6_RM mBufferE = axisem3d::eigen::RMatX6_RM(0, 6);
  axisem3d::eigen::RMatX3_RM mBufferR = axisem3d::eigen::RMatX3_RM(0, 3);
  axisem3d::eigen::RMatX6_RM mBufferS = axisem3d::eigen::RMatX6_RM(0, 6);

  /////////////////////////// process ///////////////////////////
  public:
  // process and report to group
  void
  processReport(int bufferLine,
      const channel::solid::ChannelOptions& chops,
      int stationIndex,
      axisem3d::eigen::RTensor3& bufferFields);

  private:
  // process 1: rotate
  void
  rotate(int bufferLine, const channel::solid::ChannelOptions& chops);

  ////////////////////////////////////////
  //////////////// static ////////////////
  ////////////////////////////////////////

  private:
  // expand workspace for record
  static void
  expandWorkspaceRecord(int nu_1, const channel::solid::ChannelOptions& chops) {
    // for element
    if (chops.mNeedBufferU && sUXN3.rows() < nu_1) {
      sUXN3.resize(nu_1, spectral::nPEM * 3);
      sUX3.resize(nu_1, 3);
    }

    if (chops.mNeedBufferG && sGXN9.rows() < nu_1) {
      sGXN9.resize(nu_1, spectral::nPEM * 9);
      sGX9.resize(nu_1, 9);
    }

    if (chops.mNeedBufferE && sEXN6.rows() < nu_1) {
      sEXN6.resize(nu_1, spectral::nPEM * 6);
      sEX6.resize(nu_1, 6);
    }

    if (chops.mNeedBufferR && sRXN3.rows() < nu_1) {
      sRXN3.resize(nu_1, spectral::nPEM * 3);
      sRX3.resize(nu_1, 3);
    }

    if (chops.mNeedBufferS && sSXN6.rows() < nu_1) {
      sSXN6.resize(nu_1, spectral::nPEM * 6);
      sSX6.resize(nu_1, 6);
    }
  }

  // workspace for record
  inline static axisem3d::eigen::RRow3 sU3;
  inline static axisem3d::eigen::RRow9 sG9;
  inline static axisem3d::eigen::RRow6 sE6;
  inline static axisem3d::eigen::RRow3 sR3;
  inline static axisem3d::eigen::RRow6 sS6;
  inline static axisem3d::eigen::CMatX3 sUX3 = axisem3d::eigen::CMatX3(0, 3);
  inline static axisem3d::eigen::CMatX9 sGX9 = axisem3d::eigen::CMatX9(0, 9);
  inline static axisem3d::eigen::CMatX6 sEX6 = axisem3d::eigen::CMatX6(0, 6);
  inline static axisem3d::eigen::CMatX3 sRX3 = axisem3d::eigen::CMatX3(0, 3);
  inline static axisem3d::eigen::CMatX6 sSX6 = axisem3d::eigen::CMatX6(0, 6);
  inline static axisem3d::eigen::CMatXN3 sUXN3 = axisem3d::eigen::CMatXN3(0, spectral::nPEM * 3);
  inline static axisem3d::eigen::CMatXN9 sGXN9 = axisem3d::eigen::CMatXN9(0, spectral::nPEM * 9);
  inline static axisem3d::eigen::CMatXN6 sEXN6 = axisem3d::eigen::CMatXN6(0, spectral::nPEM * 6);
  inline static axisem3d::eigen::CMatXN3 sRXN3 = axisem3d::eigen::CMatXN3(0, spectral::nPEM * 3);
  inline static axisem3d::eigen::CMatXN6 sSXN6 = axisem3d::eigen::CMatXN6(0, spectral::nPEM * 6);
};

#endif /* StationSolid_hpp */
