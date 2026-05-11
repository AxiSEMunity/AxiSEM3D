// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: MIT
// Copyright (C) 2019 - 2026 by the AxiSEM3D authors
//
// This file is part of the AxiSEM3D library. See the LICENSE file for details.
//
// -----------------------------------------------------------------------------

//  parallel NetCDF IO for element output

#ifndef ElementIO_ParNetCDF_hpp
#define ElementIO_ParNetCDF_hpp

#include "ElementIO.hpp"
#include "NetCDF_Writer.hpp"

class ElementIO_ParNetCDF : public ElementIO {
  public:
  // initialize
  void
  initialize(const std::string& groupName,
      int numRecordSteps,
      const std::vector<std::string>& channels,
      int npnts,
      const std::vector<int>& naGrid,
      const axisem3d::eigen::IMatX4_RM& elemNaInfo,
      const axisem3d::eigen::DMatXX_RM& elemCoords);

  // finalize
  void
  finalize();

  // dump to file
  void
  dumpToFile(const axisem3d::eigen::DColX& bufferTime,
      const std::vector<axisem3d::eigen::RTensor5>& bufferFields,
      int bufferLine);

  private:
  //////////////////// const ////////////////////
  // file
  std::unique_ptr<NetCDF_Writer> mNcFile = nullptr;

  // variable id
  int mVarID_Time = -1;
  std::vector<int> mVarID_Data;

  // current line in time dimension
  int mFileLineTime = 0;

  // first element index on na-grid
  std::vector<int> mFirstElemNaGrid;
};

#endif /* ElementIO_ParNetCDF_hpp */
