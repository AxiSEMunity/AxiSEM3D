// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: MIT
// Copyright (C) 2019 - 2026 by the AxiSEM3D authors
//
// This file is part of the AxiSEM3D library. See the LICENSE file for details.
//
// -----------------------------------------------------------------------------

//  parallel NetCDF IO for station output

#ifndef StationIO_ParNetCDF_hpp
#define StationIO_ParNetCDF_hpp

#include "StationIO.hpp"
#include "NetCDF_Writer.hpp"

class StationIO_ParNetCDF : public StationIO {
  public:
  // initialize
  void
  initialize(const std::string& groupName,
      int numRecordSteps,
      const std::vector<std::string>& channels,
      const std::vector<std::string>& stKeys);

  // finalize
  void
  finalize();

  // dump to file
  void
  dumpToFile(const axisem3d::eigen::DColX& bufferTime,
      const axisem3d::eigen::RTensor3& bufferFields,
      int bufferLine);

  private:
  //////////////////// const ////////////////////
  // file
  std::unique_ptr<NetCDF_Writer> mNcFile = nullptr;

  // variable id
  int mVarID_Time = -1;
  int mVarID_Data = -1;

  // current line in time dimension
  int mFileLineTime = 0;
};

#endif /* StationIO_ParNetCDF_hpp */
